#include "cov/compute_abi.h"

#include <sycl/sycl.hpp>

#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void error_text(char* error, size_t capacity, const char* message) noexcept {
    if (error && capacity) std::snprintf(error, capacity, "%s", message);
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

std::vector<sycl::device> intel_level_zero_gpus() {
    std::vector<sycl::device> result;
    for (const auto& device : sycl::device::get_devices(sycl::info::device_type::gpu)) {
        if (device.get_backend() == sycl::backend::ext_oneapi_level_zero &&
            device.get_info<sycl::info::device::vendor_id>() == 0x8086u &&
            device.has(sycl::aspect::usm_device_allocations))
            result.push_back(device);
    }
    return result;
}

struct Context {
    sycl::queue queue;
    CovGaussianTerm* terms = nullptr;
    float* values = nullptr;
    uint32_t term_count = 0;
    uint32_t term_capacity = 0;
    uint32_t value_capacity = 0;
    std::string name;

    explicit Context(const sycl::device& device)
        : queue(device, sycl::property_list{sycl::property::queue::in_order{}}),
          name(device.get_info<sycl::info::device::name>()) {}

    ~Context() noexcept {
        try {
            queue.wait_and_throw();
        } catch (...) {
            // Destructors cannot report errors across the C ABI.
        }
        try {
            if (values) sycl::free(values, queue);
            if (terms) sycl::free(terms, queue);
        } catch (...) {
            // Device loss may also make deallocation fail.
        }
    }
};

void* create(int device_index, char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    try {
        if (device_index < -1) throw std::runtime_error("Invalid SYCL device index");
        const auto devices = intel_level_zero_gpus();
        if (devices.empty()) throw std::runtime_error("No Intel Level Zero GPU was found");
        const size_t selected = device_index == -1 ? 0 : size_t(device_index);
        if (selected >= devices.size())
            throw std::runtime_error("Selected Intel Level Zero GPU index is unavailable");
        return new Context(devices[selected]);
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown SYCL initialization error");
    }
    return nullptr;
}

void destroy(void* opaque) noexcept { delete static_cast<Context*>(opaque); }

int set_terms(void* opaque, const CovGaussianTerm* source, uint32_t count,
              char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    auto* context = static_cast<Context*>(opaque);
    if (!context || !valid_terms(source, count)) {
        error_text(error, capacity, "Invalid SYCL Gaussian terms or context");
        return 0;
    }
    try {
        if (count > context->term_capacity) {
            auto* replacement = sycl::malloc_device<CovGaussianTerm>(count, context->queue);
            if (!replacement) throw std::runtime_error("Allocating SYCL term buffer failed");
            if (context->terms) sycl::free(context->terms, context->queue);
            context->terms = replacement;
            context->term_capacity = count;
        }
        context->term_count = 0;
        if (count) {
            context->queue.memcpy(context->terms, source, size_t(count) * sizeof(*source));
            context->queue.wait_and_throw();
        }
        context->term_count = count;
        return 1;
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown SYCL term upload error");
    }
    return 0;
}

int evaluate(void* opaque, const CovGridRequest* request, float* output,
             char* error, size_t capacity) noexcept {
    error_text(error, capacity, "");
    auto* context = static_cast<Context*>(opaque);
    if (!context || !request || !valid_request(*request) ||
        (request->point_count && !output)) {
        error_text(error, capacity, "Invalid SYCL grid request, output or context");
        return 0;
    }
    if (!request->point_count) return 1;
    try {
        if (request->point_count > context->value_capacity) {
            auto* replacement = sycl::malloc_device<float>(request->point_count, context->queue);
            if (!replacement) throw std::runtime_error("Allocating SYCL grid buffer failed");
            if (context->values) sycl::free(context->values, context->queue);
            context->values = replacement;
            context->value_capacity = request->point_count;
        }

        const CovGridRequest grid = *request;
        const CovGaussianTerm* terms = context->terms;
        const uint32_t term_count = context->term_count;
        float* values = context->values;
        context->queue.parallel_for(sycl::range<1>(size_t(grid.point_count)), [=](sycl::id<1> id) {
            const uint32_t local = uint32_t(id[0]);
            const uint64_t point = grid.first_point + local;
            const uint64_t plane = uint64_t(grid.nx) * grid.ny;
            const uint32_t ix = uint32_t(point % grid.nx);
            const uint32_t iy = uint32_t((point / grid.nx) % grid.ny);
            const uint32_t iz = uint32_t(point / plane);
            const float x = grid.nx == 1 ? grid.min_x :
                grid.min_x + (grid.max_x - grid.min_x) * (float(ix) / float(grid.nx - 1));
            const float y = grid.ny == 1 ? grid.min_y :
                grid.min_y + (grid.max_y - grid.min_y) * (float(iy) / float(grid.ny - 1));
            const float z = grid.nz == 1 ? grid.min_z :
                grid.min_z + (grid.max_z - grid.min_z) * (float(iz) / float(grid.nz - 1));

            float sum = 0.0f;
            for (uint32_t i = 0; i < term_count; ++i) {
                const CovGaussianTerm t = terms[i];
                const float dx = x - t.cx;
                const float dy = y - t.cy;
                const float dz = z - t.cz;
                const float r2 = dx * dx + dy * dy + dz * dz;
                auto power = [](float base, uint32_t exponent) {
                    float result = 1.0f;
                    while (exponent) {
                        if (exponent & 1u) result *= base;
                        base *= base;
                        exponent >>= 1u;
                    }
                    return result;
                };
                sum += t.coefficient * power(dx, t.ax) * power(dy, t.ay) *
                       power(dz, t.az) * sycl::exp(-t.exponent * r2);
            }
            values[local] = sum;
        });
        context->queue.memcpy(output, context->values, size_t(request->point_count) * sizeof(float));
        context->queue.wait_and_throw();
        return 1;
    } catch (const std::exception& ex) {
        error_text(error, capacity, ex.what());
    } catch (...) {
        error_text(error, capacity, "Unknown SYCL grid evaluation error");
    }
    return 0;
}

const char* device_name(void* opaque) noexcept {
    const auto* context = static_cast<const Context*>(opaque);
    return context ? context->name.c_str() : "";
}

const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "sycl", create, destroy,
                           set_terms, evaluate, device_name};

} // namespace

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
