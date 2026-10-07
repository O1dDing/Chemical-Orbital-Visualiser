#include "cov/compute_abi.h"
#include "shader.hpp"
#include <webgpu/wgpu.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <thread>
#include <vector>

namespace {
constexpr uint32_t kNativeVersion = 0x1d000101u; // 29.0.1.1
constexpr uint32_t kWorkgroup = 64;
constexpr uint32_t kMaxChunkPoints = 1u << 20; // bound staging memory to 4 MiB
constexpr auto kTimeout = std::chrono::seconds(30);

void error_text(char* dst, size_t cap, const std::string& value) noexcept {
    if (!dst || !cap) return;
    const size_t n = std::min(cap - 1, value.size());
    std::memcpy(dst, value.data(), n);
    dst[n] = 0;
}
std::string str(WGPUStringView value) {
    if (!value.data) return {};
    return {value.data, value.length == WGPU_STRLEN ? std::strlen(value.data) : value.length};
}
WGPUStringView view(const char* s) { return {s, WGPU_STRLEN}; }
bool finite(float x) { return std::isfinite(x); }

struct AsyncState {
    bool done = false;
    int status = 0;
    std::string message;
    WGPUDevice device = nullptr;
    ~AsyncState() { if (device) wgpuDeviceRelease(device); }
};

// The heap-held shared reference stays valid even if a timed-out operation later calls back.
// A permanently stalled vendor callback retains only this small state object.
using CallbackRef = std::shared_ptr<AsyncState>*;
CallbackRef callback_ref(const std::shared_ptr<AsyncState>& state) { return new std::shared_ptr<AsyncState>(state); }
void device_ready(WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message,
                  void* userdata, void*) noexcept {
    std::unique_ptr<std::shared_ptr<AsyncState>> holder(static_cast<CallbackRef>(userdata));
    auto& state = **holder;
    state.status = static_cast<int>(status);
    state.device = device;
    try { state.message = str(message); } catch (...) { state.status = 0; }
    state.done = true;
}
void map_ready(WGPUMapAsyncStatus status, WGPUStringView message, void* userdata, void*) noexcept {
    std::unique_ptr<std::shared_ptr<AsyncState>> holder(static_cast<CallbackRef>(userdata));
    auto& state = **holder;
    state.status = static_cast<int>(status);
    try { state.message = str(message); } catch (...) { state.status = 0; }
    state.done = true;
}
void scope_ready(WGPUPopErrorScopeStatus status, WGPUErrorType type, WGPUStringView message,
                 void* userdata, void*) noexcept {
    std::unique_ptr<std::shared_ptr<AsyncState>> holder(static_cast<CallbackRef>(userdata));
    auto& state = **holder;
    state.status = status == WGPUPopErrorScopeStatus_Success && type == WGPUErrorType_NoError ? 1 : 0;
    try { state.message = str(message); } catch (...) { state.status = 0; }
    state.done = true;
}
bool await(WGPUInstance instance, WGPUDevice device, const std::shared_ptr<AsyncState>& state) {
    const auto deadline = std::chrono::steady_clock::now() + kTimeout;
    while (!state->done && std::chrono::steady_clock::now() < deadline) {
        if (device) wgpuDevicePoll(device, WGPU_FALSE, nullptr);
        wgpuInstanceProcessEvents(instance);
        if (!state->done) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return state->done;
}

struct Params {
    float min_xyz[4], max_xyz[4];
    uint32_t dims[4], base_xyz_count[4];
};
static_assert(sizeof(Params) == 64);

struct Context {
    WGPUInstance instance = nullptr;
    WGPUAdapter adapter = nullptr;
    WGPUDevice device = nullptr;
    WGPUQueue queue = nullptr;
    WGPUBindGroupLayout layout = nullptr;
    WGPUPipelineLayout pipeline_layout = nullptr;
    WGPUComputePipeline pipeline = nullptr;
    WGPUBuffer terms_buffer = nullptr;
    uint32_t term_count = 0;
    WGPULimits limits = WGPU_LIMITS_INIT;
    std::string name;

    ~Context() {
        if (terms_buffer) { wgpuBufferDestroy(terms_buffer); wgpuBufferRelease(terms_buffer); }
        if (pipeline) wgpuComputePipelineRelease(pipeline);
        if (pipeline_layout) wgpuPipelineLayoutRelease(pipeline_layout);
        if (layout) wgpuBindGroupLayoutRelease(layout);
        if (queue) wgpuQueueRelease(queue);
        if (device) { wgpuDeviceDestroy(device); wgpuDeviceRelease(device); }
        if (adapter) wgpuAdapterRelease(adapter);
        if (instance) wgpuInstanceRelease(instance);
    }
};

bool pop_scope(Context& ctx, std::string& reason) {
    auto state = std::make_shared<AsyncState>();
    WGPUPopErrorScopeCallbackInfo info = WGPU_POP_ERROR_SCOPE_CALLBACK_INFO_INIT;
    info.mode = WGPUCallbackMode_AllowProcessEvents;
    info.callback = scope_ready;
    info.userdata1 = callback_ref(state);
    wgpuDevicePopErrorScope(ctx.device, info);
    if (!await(ctx.instance, ctx.device, state)) { reason = "WebGPU validation timed out"; return false; }
    if (state->status != 1) { reason = state->message.empty() ? "WebGPU validation failed" : state->message; return false; }
    return true;
}

bool setup_pipeline(Context& ctx, std::string& reason) {
    wgpuDevicePushErrorScope(ctx.device, WGPUErrorFilter_Validation);
    WGPUShaderSourceWGSL source = WGPU_SHADER_SOURCE_WGSL_INIT;
    source.code = view(cov_webgpu_shader);
    WGPUShaderModuleDescriptor shader_desc = WGPU_SHADER_MODULE_DESCRIPTOR_INIT;
    shader_desc.nextInChain = &source.chain;
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(ctx.device, &shader_desc);
    if (!shader) { pop_scope(ctx, reason); if (reason.empty()) reason = "WebGPU shader creation failed"; return false; }
    WGPUBindGroupLayoutEntry entries[3] = { WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT,
                                            WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT,
                                            WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT };
    const WGPUBufferBindingType types[3] = {WGPUBufferBindingType_ReadOnlyStorage,
                                            WGPUBufferBindingType_Uniform,
                                            WGPUBufferBindingType_Storage};
    for (uint32_t i = 0; i < 3; ++i) {
        entries[i].binding = i;
        entries[i].visibility = WGPUShaderStage_Compute;
        entries[i].buffer.type = types[i];
    }
    WGPUBindGroupLayoutDescriptor layout_desc = WGPU_BIND_GROUP_LAYOUT_DESCRIPTOR_INIT;
    layout_desc.entryCount = 3;
    layout_desc.entries = entries;
    ctx.layout = wgpuDeviceCreateBindGroupLayout(ctx.device, &layout_desc);
    WGPUPipelineLayoutDescriptor pipeline_layout_desc = WGPU_PIPELINE_LAYOUT_DESCRIPTOR_INIT;
    pipeline_layout_desc.bindGroupLayoutCount = 1;
    pipeline_layout_desc.bindGroupLayouts = &ctx.layout;
    if (ctx.layout) ctx.pipeline_layout = wgpuDeviceCreatePipelineLayout(ctx.device, &pipeline_layout_desc);
    WGPUComputePipelineDescriptor pipeline_desc = WGPU_COMPUTE_PIPELINE_DESCRIPTOR_INIT;
    pipeline_desc.layout = ctx.pipeline_layout;
    pipeline_desc.compute.module = shader;
    pipeline_desc.compute.entryPoint = view("main");
    if (ctx.pipeline_layout) ctx.pipeline = wgpuDeviceCreateComputePipeline(ctx.device, &pipeline_desc);
    wgpuShaderModuleRelease(shader);
    const bool valid = pop_scope(ctx, reason);
    if (!valid || !ctx.pipeline) {
        if (reason.empty()) reason = "WebGPU compute pipeline creation failed";
        return false;
    }
    return true;
}

WGPUBuffer buffer(Context& ctx, uint64_t bytes, WGPUBufferUsage usage) {
    WGPUBufferDescriptor desc = WGPU_BUFFER_DESCRIPTOR_INIT;
    desc.size = bytes;
    desc.usage = usage;
    return wgpuDeviceCreateBuffer(ctx.device, &desc);
}

void* create_impl(int index, char* error, size_t capacity) {
    if (index < -1) { error_text(error, capacity, "Invalid WebGPU device index"); return nullptr; }
    if (wgpuGetVersion() != kNativeVersion) {
        error_text(error, capacity, "wgpu-native runtime version differs from 29.0.1.1"); return nullptr;
    }
    auto ctx = std::make_unique<Context>();
    WGPUInstanceDescriptor instance_desc = WGPU_INSTANCE_DESCRIPTOR_INIT;
    ctx->instance = wgpuCreateInstance(&instance_desc);
    if (!ctx->instance) { error_text(error, capacity, "Could not create WebGPU instance"); return nullptr; }
    const size_t count = wgpuInstanceEnumerateAdapters(ctx->instance, nullptr, nullptr);
    std::vector<WGPUAdapter> adapters(count);
    if (count) wgpuInstanceEnumerateAdapters(ctx->instance, nullptr, adapters.data());
    int gpu_index = 0;
    for (auto adapter : adapters) {
        WGPUAdapterInfo info = WGPU_ADAPTER_INFO_INIT;
        if (wgpuAdapterGetInfo(adapter, &info) != WGPUStatus_Success) { wgpuAdapterRelease(adapter); continue; }
        const bool gpu = info.adapterType == WGPUAdapterType_DiscreteGPU ||
                         info.adapterType == WGPUAdapterType_IntegratedGPU;
        if (gpu && ((index == -1 && !ctx->adapter) || index == gpu_index)) {
            ctx->adapter = adapter;
            ctx->name = str(info.device);
        } else {
            wgpuAdapterRelease(adapter);
        }
        if (gpu) ++gpu_index;
        wgpuAdapterInfoFreeMembers(info);
    }
    if (!ctx->adapter) { error_text(error, capacity, "No suitable WebGPU GPU adapter found"); return nullptr; }
    auto state = std::make_shared<AsyncState>();
    WGPURequestDeviceCallbackInfo callback = WGPU_REQUEST_DEVICE_CALLBACK_INFO_INIT;
    callback.mode = WGPUCallbackMode_AllowProcessEvents;
    callback.callback = device_ready;
    callback.userdata1 = callback_ref(state);
    WGPUDeviceDescriptor device_desc = WGPU_DEVICE_DESCRIPTOR_INIT;
    wgpuAdapterRequestDevice(ctx->adapter, &device_desc, callback);
    if (!await(ctx->instance, nullptr, state)) { error_text(error, capacity, "WebGPU device request timed out"); return nullptr; }
    if (state->status != WGPURequestDeviceStatus_Success || !state->device) {
        error_text(error, capacity, state->message.empty() ? "WebGPU device request failed" : state->message);
        return nullptr;
    }
    ctx->device = state->device;
    state->device = nullptr;
    ctx->queue = wgpuDeviceGetQueue(ctx->device);
    if (!ctx->queue || wgpuDeviceGetLimits(ctx->device, &ctx->limits) != WGPUStatus_Success) {
        error_text(error, capacity, "Could not query WebGPU device limits"); return nullptr;
    }
    std::string reason;
    if (!setup_pipeline(*ctx, reason)) { error_text(error, capacity, reason); return nullptr; }
    error_text(error, capacity, "");
    return ctx.release();
}

int set_terms_impl(void* opaque, const CovGaussianTerm* terms, uint32_t count,
                   char* error, size_t capacity) {
    auto* ctx = static_cast<Context*>(opaque);
    if (!ctx || (count && !terms)) { error_text(error, capacity, "Invalid WebGPU term input"); return 0; }
    for (uint32_t i = 0; i < count; ++i) {
        const auto& t = terms[i];
        if (!finite(t.cx) || !finite(t.cy) || !finite(t.cz) || !finite(t.exponent) ||
            !finite(t.coefficient) || t.exponent <= 0 || t.ax > 4 || t.ay > 4 || t.az > 4) {
            error_text(error, capacity, "Invalid Gaussian term"); return 0;
        }
    }
    const uint64_t bytes = uint64_t(count) * sizeof(CovGaussianTerm);
    if (bytes > ctx->limits.maxStorageBufferBindingSize || bytes > ctx->limits.maxBufferSize) {
        error_text(error, capacity, "Gaussian terms exceed WebGPU storage buffer limit"); return 0;
    }
    // WGSL runtime arrays cannot be bound with an empty zero-byte buffer.
    wgpuDevicePushErrorScope(ctx->device, WGPUErrorFilter_Validation);
    WGPUBuffer next = buffer(*ctx, std::max<uint64_t>(bytes, sizeof(CovGaussianTerm)),
                             WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst);
    if (next && bytes) wgpuQueueWriteBuffer(ctx->queue, next, 0, terms, static_cast<size_t>(bytes));
    std::string validation_reason;
    if (!pop_scope(*ctx, validation_reason)) {
        if (next) { wgpuBufferDestroy(next); wgpuBufferRelease(next); }
        error_text(error, capacity, validation_reason); return 0;
    }
    if (!next) { error_text(error, capacity, "Could not allocate WebGPU term buffer"); return 0; }
    if (ctx->terms_buffer) { wgpuBufferDestroy(ctx->terms_buffer); wgpuBufferRelease(ctx->terms_buffer); }
    ctx->terms_buffer = next;
    ctx->term_count = count;
    error_text(error, capacity, "");
    return 1;
}

int evaluate_impl(void* opaque, const CovGridRequest* req, float* output,
                  char* error, size_t capacity) {
    auto* ctx = static_cast<Context*>(opaque);
    if (!ctx || !req || (req->point_count && !output) || !ctx->terms_buffer) {
        error_text(error, capacity, "Invalid WebGPU grid request or terms not set"); return 0;
    }
    const float coords[] = {req->min_x, req->min_y, req->min_z, req->max_x, req->max_y, req->max_z};
    for (float x : coords) if (!finite(x)) { error_text(error, capacity, "Non-finite grid coordinate"); return 0; }
    if (req->min_x > req->max_x || req->min_y > req->max_y || req->min_z > req->max_z ||
        !finite(req->max_x - req->min_x) || !finite(req->max_y - req->min_y) ||
        !finite(req->max_z - req->min_z)) {
        error_text(error, capacity, "Invalid grid bounds"); return 0;
    }
    if (!req->nx || !req->ny || !req->nz || req->nx > 0x7fffffffu ||
        req->ny > 0x7fffffffu || req->nz > 0x7fffffffu) {
        error_text(error, capacity, "Invalid WebGPU grid dimensions"); return 0;
    }
    const uint64_t nx = req->nx, ny = req->ny, nz = req->nz;
    if (nx > UINT64_MAX / ny || nx * ny > UINT64_MAX / nz) {
        error_text(error, capacity, "Grid dimensions overflow"); return 0;
    }
    const uint64_t total = nx * ny * nz;
    if (req->first_point > total || req->point_count > total - req->first_point) {
        error_text(error, capacity, "Grid chunk is outside the volume"); return 0;
    }
    if (!req->point_count) { error_text(error, capacity, ""); return 1; }
    const uint64_t chunk_limit = std::min<uint64_t>({kMaxChunkPoints,
        ctx->limits.maxStorageBufferBindingSize / sizeof(float),
        ctx->limits.maxBufferSize / sizeof(float),
        uint64_t(ctx->limits.maxComputeWorkgroupsPerDimension) * kWorkgroup});
    if (!chunk_limit || ctx->limits.maxUniformBufferBindingSize < sizeof(Params)) {
        error_text(error, capacity, "WebGPU device buffer or dispatch limit is too small"); return 0;
    }
    for (uint64_t done = 0; done < req->point_count;) {
        const uint32_t n = static_cast<uint32_t>(std::min<uint64_t>(chunk_limit, req->point_count - done));
        const uint64_t first = req->first_point + done;
        Params params{};
        params.min_xyz[0] = req->min_x; params.min_xyz[1] = req->min_y; params.min_xyz[2] = req->min_z;
        params.max_xyz[0] = req->max_x; params.max_xyz[1] = req->max_y; params.max_xyz[2] = req->max_z;
        params.dims[0] = req->nx; params.dims[1] = req->ny;
        params.dims[2] = req->nz; params.dims[3] = ctx->term_count;
        params.base_xyz_count[0] = static_cast<uint32_t>(first % nx);
        params.base_xyz_count[1] = static_cast<uint32_t>((first / nx) % ny);
        params.base_xyz_count[2] = static_cast<uint32_t>(first / (nx * ny));
        params.base_xyz_count[3] = n;
        const uint64_t bytes = uint64_t(n) * sizeof(float);
        WGPUBuffer param = buffer(*ctx, sizeof(Params), WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst);
        WGPUBuffer result = buffer(*ctx, bytes, WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc);
        WGPUBuffer readback = buffer(*ctx, bytes, WGPUBufferUsage_CopyDst | WGPUBufferUsage_MapRead);
        WGPUBindGroup group = nullptr;
        WGPUCommandEncoder encoder = nullptr;
        WGPUCommandBuffer command = nullptr;
        bool mapped = false;
        std::string reason;
        wgpuDevicePushErrorScope(ctx->device, WGPUErrorFilter_Validation);
        if (!param || !result || !readback) reason = "Could not allocate WebGPU grid buffers";
        if (reason.empty()) {
            wgpuQueueWriteBuffer(ctx->queue, param, 0, &params, sizeof(params));
            WGPUBindGroupEntry entries[3] = { WGPU_BIND_GROUP_ENTRY_INIT,
                                              WGPU_BIND_GROUP_ENTRY_INIT,
                                              WGPU_BIND_GROUP_ENTRY_INIT };
            entries[0].binding = 0; entries[0].buffer = ctx->terms_buffer;
            entries[0].size = std::max<uint64_t>(uint64_t(ctx->term_count) * sizeof(CovGaussianTerm), sizeof(CovGaussianTerm));
            entries[1].binding = 1; entries[1].buffer = param; entries[1].size = sizeof(Params);
            entries[2].binding = 2; entries[2].buffer = result; entries[2].size = bytes;
            WGPUBindGroupDescriptor desc = WGPU_BIND_GROUP_DESCRIPTOR_INIT;
            desc.layout = ctx->layout; desc.entryCount = 3; desc.entries = entries;
            group = wgpuDeviceCreateBindGroup(ctx->device, &desc);
            if (!group) reason = "Could not bind WebGPU grid buffers";
        }
        if (reason.empty()) {
            encoder = wgpuDeviceCreateCommandEncoder(ctx->device, nullptr);
            if (!encoder) reason = "Could not create WebGPU command encoder";
        }
        if (reason.empty()) {
            WGPUComputePassEncoder pass = wgpuCommandEncoderBeginComputePass(encoder, nullptr);
            if (!pass) reason = "Could not begin WebGPU compute pass";
            else {
                wgpuComputePassEncoderSetPipeline(pass, ctx->pipeline);
                wgpuComputePassEncoderSetBindGroup(pass, 0, group, 0, nullptr);
                wgpuComputePassEncoderDispatchWorkgroups(pass, (n + kWorkgroup - 1) / kWorkgroup, 1, 1);
                wgpuComputePassEncoderEnd(pass);
                wgpuComputePassEncoderRelease(pass);
                wgpuCommandEncoderCopyBufferToBuffer(encoder, result, 0, readback, 0, bytes);
                command = wgpuCommandEncoderFinish(encoder, nullptr);
                if (!command) reason = "Could not finish WebGPU command buffer";
            }
        }
        if (reason.empty()) {
            wgpuQueueSubmit(ctx->queue, 1, &command);
            auto state = std::make_shared<AsyncState>();
            WGPUBufferMapCallbackInfo info = WGPU_BUFFER_MAP_CALLBACK_INFO_INIT;
            info.mode = WGPUCallbackMode_AllowProcessEvents;
            info.callback = map_ready;
            info.userdata1 = callback_ref(state);
            wgpuBufferMapAsync(readback, WGPUMapMode_Read, 0, static_cast<size_t>(bytes), info);
            if (!await(ctx->instance, ctx->device, state)) reason = "WebGPU readback timed out";
            else if (state->status != WGPUMapAsyncStatus_Success)
                reason = state->message.empty() ? "WebGPU readback failed" : state->message;
            else {
                mapped = true;
                const void* src = wgpuBufferGetConstMappedRange(readback, 0, static_cast<size_t>(bytes));
                if (!src) reason = "WebGPU mapped readback is unavailable";
                else std::memcpy(output + done, src, static_cast<size_t>(bytes));
            }
        }
        if (mapped) wgpuBufferUnmap(readback);
        if (command) wgpuCommandBufferRelease(command);
        if (encoder) wgpuCommandEncoderRelease(encoder);
        if (group) wgpuBindGroupRelease(group);
        if (readback) { wgpuBufferDestroy(readback); wgpuBufferRelease(readback); }
        if (result) { wgpuBufferDestroy(result); wgpuBufferRelease(result); }
        if (param) { wgpuBufferDestroy(param); wgpuBufferRelease(param); }
        std::string validation_reason;
        if (!pop_scope(*ctx, validation_reason) && reason.empty()) reason = validation_reason;
        if (!reason.empty()) { error_text(error, capacity, reason); return 0; }
        done += n;
    }
    error_text(error, capacity, "");
    return 1;
}

void* create(int index, char* error, size_t capacity) noexcept {
    try { return create_impl(index, error, capacity); }
    catch (const std::exception& ex) { error_text(error, capacity, ex.what()); return nullptr; }
    catch (...) { error_text(error, capacity, "WebGPU initialization failed"); return nullptr; }
}
void destroy(void* context) noexcept { try { delete static_cast<Context*>(context); } catch (...) {} }
int set_terms(void* c, const CovGaussianTerm* t, uint32_t n, char* e, size_t cap) noexcept {
    try { return set_terms_impl(c, t, n, e, cap); }
    catch (const std::exception& ex) { error_text(e, cap, ex.what()); return 0; }
    catch (...) { error_text(e, cap, "WebGPU term upload failed"); return 0; }
}
int evaluate(void* c, const CovGridRequest* r, float* o, char* e, size_t cap) noexcept {
    try { return evaluate_impl(c, r, o, e, cap); }
    catch (const std::exception& ex) { error_text(e, cap, ex.what()); return 0; }
    catch (...) { error_text(e, cap, "WebGPU evaluation failed"); return 0; }
}
const char* device_name(void* c) noexcept {
    return c ? static_cast<Context*>(c)->name.c_str() : "";
}
const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "webgpu", create, destroy, set_terms, evaluate, device_name};
} // namespace

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
