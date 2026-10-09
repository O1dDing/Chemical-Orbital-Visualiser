#include "cov/compute_abi.h"

#include <hip/hip_runtime.h>

#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

#ifndef __HIP_PLATFORM_AMD__
#error "The COV HIP module must be built with the AMD HIP platform"
#endif

namespace {

void error_text(char* error, size_t capacity, const char* message) noexcept {
    if (error && capacity) std::snprintf(error, capacity, "%s", message);
}

void hip_check(hipError_t status, const char* action) {
    if (status != hipSuccess) {
        throw std::runtime_error(std::string(action) + ": " + hipGetErrorString(status));
    }
}

bool valid_request(const CovGridRequest& r) noexcept {
    if (!r.nx || !r.ny || !r.nz) return false;
    if (!std::isfinite(r.min_x) || !std::isfinite(r.min_y) || !std::isfinite(r.min_z) ||
        !std::isfinite(r.max_x) || !std::isfinite(r.max_y) || !std::isfinite(r.max_z) ||
        !std::isfinite(r.max_x - r.min_x) || !std::isfinite(r.max_y - r.min_y) ||
        !std::isfinite(r.max_z - r.min_z)) return false;
    const uint64_t max = std::numeric_limits<uint64_t>::max();
    if (r.nx > max / r.ny) return false;
    const uint64_t plane = uint64_t(r.nx) * r.ny;
    if (plane > max / r.nz) return false;
    const uint64_t total = plane * r.nz;
    return r.first_point <= total && r.point_count <= total - r.first_point;
}

bool valid_terms(const CovGaussianTerm* terms, uint32_t count) noexcept {
    if (count && !terms) return false;
    for (uint32_t i = 0; i < count; ++i) {
        const auto& t = terms[i];
        if (!std::isfinite(t.cx) || !std::isfinite(t.cy) || !std::isfinite(t.cz) ||
            !std::isfinite(t.exponent) || !std::isfinite(t.coefficient) || t.exponent < 0.0f)
            return false;
    }
    return true;
}

struct Context {
    int device = -1;
    hipStream_t stream = nullptr;
    CovGaussianTerm* terms = nullptr;
    float* values = nullptr;
    uint32_t term_count = 0;
    uint32_t term_capacity = 0;
    uint32_t value_capacity = 0;
    std::string name;

    ~Context() noexcept {
        if (device >= 0 && hipSetDevice(device) == hipSuccess) {
            if (stream) (void)hipStreamSynchronize(stream);
            if (values) (void)hipFree(values);
            if (terms) (void)hipFree(terms);
            if (stream) (void)hipStreamDestroy(stream);
        }
    }
};

__device__ float integer_power(float base, uint32_t power) {
    float result = 1.0f;
    while (power) {
        if (power & 1u) result *= base;
        base *= base;
        power >>= 1u;
    }
    return result;
}

__global__ void evaluate_kernel(CovGridRequest request, const CovGaussianTerm* terms,
                                uint32_t term_count, float* output) {
    const uint32_t local = blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= request.point_count) return;

    const uint64_t point = request.first_point + local;
    const uint64_t plane = uint64_t(request.nx) * request.ny;
    const uint32_t ix = uint32_t(point % request.nx);
    const uint32_t iy = uint32_t((point / request.nx) % request.ny);
    const uint32_t iz = uint32_t(point / plane);
    const float x = request.nx == 1 ? request.min_x :
        request.min_x + (request.max_x - request.min_x) * (float(ix) / float(request.nx - 1));
    const float y = request.ny == 1 ? request.min_y :
        request.min_y + (request.max_y - request.min_y) * (float(iy) / float(request.ny - 1));
    const float z = request.nz == 1 ? request.min_z :
        request.min_z + (request.max_z - request.min_z) * (float(iz) / float(request.nz - 1));

    float sum = 0.0f;
    for (uint32_t i = 0; i < term_count; ++i) {
        const CovGaussianTerm t = terms[i];
        if (t.coefficient == 0.0f) continue;
        const float dx = x - t.cx;
        const float dy = y - t.cy;
        const float dz = z - t.cz;
        const float r2 = dx * dx + dy * dy + dz * dz;
        const float radial = expf(-t.exponent * r2);
        if (radial == 0.0f) continue;
        sum += t.coefficient * integer_power(dx, t.ax) * integer_power(dy, t.ay) *
               integer_power(dz, t.az) * radial;
    }
    output[local] = sum;
}

void* create(int device_index, char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    try {
        if (device_index < -1) throw std::runtime_error("Invalid HIP device index");
        int count = 0;
        hip_check(hipGetDeviceCount(&count), "Enumerating AMD GPUs");
        if (count == 0) throw std::runtime_error("No AMD HIP GPU was found");
        const int selected = device_index == -1 ? 0 : device_index;
        if (selected >= count) throw std::runtime_error("Selected HIP device index is unavailable");
        hipDeviceProp_t properties{};
        hip_check(hipGetDeviceProperties(&properties, selected), "Reading HIP device properties");
        auto context = std::make_unique<Context>();
        context->device = selected;
        context->name = properties.name;
        hip_check(hipSetDevice(selected), "Selecting AMD GPU");
        hip_check(hipStreamCreateWithFlags(&context->stream, hipStreamNonBlocking),
                  "Creating HIP stream");
        return context.release();
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown HIP initialization error");
    }
    return nullptr;
}

void destroy(void* opaque) noexcept { delete static_cast<Context*>(opaque); }

int set_terms(void* opaque, const CovGaussianTerm* source, uint32_t count,
              char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    auto* context = static_cast<Context*>(opaque);
    if (!context || !valid_terms(source, count)) {
        error_text(error, capacity, "Invalid HIP Gaussian terms or context");
        return 0;
    }
    try {
        hip_check(hipSetDevice(context->device), "Selecting AMD GPU");
        if (count > context->term_capacity) {
            CovGaussianTerm* replacement = nullptr;
            hip_check(hipMalloc(reinterpret_cast<void**>(&replacement), size_t(count) * sizeof(*source)),
                      "Allocating HIP term buffer");
            if (context->terms) (void)hipFree(context->terms);
            context->terms = replacement;
            context->term_capacity = count;
        }
        context->term_count = 0;
        if (count) {
            hip_check(hipMemcpyAsync(context->terms, source, size_t(count) * sizeof(*source),
                                     hipMemcpyHostToDevice, context->stream), "Uploading HIP terms");
            hip_check(hipStreamSynchronize(context->stream), "Completing HIP term upload");
        }
        context->term_count = count;
        return 1;
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown HIP term upload error");
    }
    return 0;
}

int evaluate(void* opaque, const CovGridRequest* request, float* output,
             char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    auto* context = static_cast<Context*>(opaque);
    if (!context || !request || !valid_request(*request) ||
        (request->point_count && !output)) {
        error_text(error, capacity, "Invalid HIP grid request, output or context");
        return 0;
    }
    if (!request->point_count) return 1;
    try {
        hip_check(hipSetDevice(context->device), "Selecting AMD GPU");
        if (request->point_count > context->value_capacity) {
            float* replacement = nullptr;
            hip_check(hipMalloc(reinterpret_cast<void**>(&replacement),
                                size_t(request->point_count) * sizeof(float)),
                      "Allocating HIP grid buffer");
            if (context->values) (void)hipFree(context->values);
            context->values = replacement;
            context->value_capacity = request->point_count;
        }
        constexpr uint32_t block_size = 256;
        const uint32_t blocks = (request->point_count - 1) / block_size + 1;
        hipLaunchKernelGGL(evaluate_kernel, dim3(blocks), dim3(block_size), 0,
                           context->stream, *request, context->terms, context->term_count,
                           context->values);
        hip_check(hipGetLastError(), "Launching HIP grid kernel");
        hip_check(hipMemcpyAsync(output, context->values,
                                 size_t(request->point_count) * sizeof(float),
                                 hipMemcpyDeviceToHost, context->stream), "Reading HIP grid");
        hip_check(hipStreamSynchronize(context->stream), "Completing HIP grid evaluation");
        return 1;
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown HIP grid evaluation error");
    }
    return 0;
}

const char* device_name(void* opaque) noexcept {
    const auto* context = static_cast<const Context*>(opaque);
    return context ? context->name.c_str() : "";
}

const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "hip", create, destroy,
                           set_terms, evaluate, device_name};

} // namespace

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
