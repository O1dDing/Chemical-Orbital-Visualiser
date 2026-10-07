#include "cov/compute_abi.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <mach-o/getsect.h>
#include <mach-o/ldsyms.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>

namespace {
constexpr uint32_t kChunkPoints = 1u << 20;

void error_text(char* dst, size_t cap, const std::string& message) noexcept {
    if (!dst || !cap) return;
    const size_t n = std::min(cap - 1, message.size());
    std::memcpy(dst, message.data(), n);
    dst[n] = 0;
}
bool finite(float x) { return std::isfinite(x); }
std::string ns_error(NSError* error, const char* fallback) {
    return error ? std::string([[error localizedDescription] UTF8String]) : fallback;
}

struct Context {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> pipeline;
    id<MTLBuffer> terms;
    uint32_t term_count = 0;
    std::string name;
    MTLResourceOptions storage = MTLResourceStorageModeShared;
};

void* create_impl(int index, char* error, size_t capacity) {
    if (index < -1) { error_text(error, capacity, "Invalid Metal device index"); return nullptr; }
    NSArray<id<MTLDevice>>* devices = MTLCopyAllDevices();
    if (!devices || !devices.count || (index >= 0 && static_cast<NSUInteger>(index) >= devices.count)) {
        error_text(error, capacity, "No suitable Metal GPU found"); return nullptr;
    }
    auto ctx = std::make_unique<Context>();
    ctx->device = devices[index == -1 ? 0 : index];
    ctx->name = [[ctx->device name] UTF8String];
    ctx->storage = [ctx->device respondsToSelector:@selector(hasUnifiedMemory)] &&
                   [ctx->device hasUnifiedMemory] ? MTLResourceStorageModeShared : MTLResourceStorageModeManaged;
    ctx->queue = [ctx->device newCommandQueue];
    if (!ctx->queue) { error_text(error, capacity, "Could not create Metal command queue"); return nullptr; }
    unsigned long size = 0;
    const uint8_t* bytes = getsectiondata(&_mh_dylib_header, "__TEXT", "__cov_metal", &size);
    if (!bytes || !size) { error_text(error, capacity, "Embedded Metal shader is missing"); return nullptr; }
    dispatch_data_t data = dispatch_data_create(bytes, size, dispatch_get_global_queue(QOS_CLASS_DEFAULT, 0),
                                                 ^{});
    NSError* ns_err = nil;
    id<MTLLibrary> library = [ctx->device newLibraryWithData:data error:&ns_err];
    if (!library) { error_text(error, capacity, ns_error(ns_err, "Could not load Metal shader")); return nullptr; }
    id<MTLFunction> function = [library newFunctionWithName:@"cov_orbital"];
    if (!function) { error_text(error, capacity, "Metal shader entry point is missing"); return nullptr; }
    ctx->pipeline = [ctx->device newComputePipelineStateWithFunction:function error:&ns_err];
    if (!ctx->pipeline) {
        error_text(error, capacity, ns_error(ns_err, "Could not create Metal compute pipeline")); return nullptr;
    }
    error_text(error, capacity, "");
    return ctx.release();
}

int set_terms_impl(void* opaque, const CovGaussianTerm* terms, uint32_t count,
                   char* error, size_t capacity) {
    auto* ctx = static_cast<Context*>(opaque);
    if (!ctx || (count && !terms)) { error_text(error, capacity, "Invalid Metal term input"); return 0; }
    for (uint32_t i = 0; i < count; ++i) {
        const auto& t = terms[i];
        if (!finite(t.cx) || !finite(t.cy) || !finite(t.cz) || !finite(t.exponent) ||
            !finite(t.coefficient) || t.exponent <= 0 || t.ax > 4 || t.ay > 4 || t.az > 4) {
            error_text(error, capacity, "Invalid Gaussian term"); return 0;
        }
    }
    const uint64_t bytes = uint64_t(count) * sizeof(CovGaussianTerm);
    if (bytes > [ctx->device maxBufferLength]) {
        error_text(error, capacity, "Gaussian terms exceed Metal buffer limit"); return 0;
    }
    id<MTLBuffer> next = [ctx->device newBufferWithLength:std::max<uint64_t>(bytes, sizeof(CovGaussianTerm))
                                                 options:ctx->storage];
    if (!next) { error_text(error, capacity, "Could not allocate Metal term buffer"); return 0; }
    if (bytes) std::memcpy([next contents], terms, static_cast<size_t>(bytes));
    if (ctx->storage == MTLResourceStorageModeManaged) [next didModifyRange:NSMakeRange(0, next.length)];
    ctx->terms = next;
    ctx->term_count = count;
    error_text(error, capacity, "");
    return 1;
}

int evaluate_impl(void* opaque, const CovGridRequest* req, float* output,
                  char* error, size_t capacity) {
    auto* ctx = static_cast<Context*>(opaque);
    if (!ctx || !req || (req->point_count && !output) || !ctx->terms) {
        error_text(error, capacity, "Invalid Metal grid request or terms not set"); return 0;
    }
    const float coords[] = {req->min_x, req->min_y, req->min_z, req->max_x, req->max_y, req->max_z};
    for (float x : coords) if (!finite(x)) { error_text(error, capacity, "Non-finite grid coordinate"); return 0; }
    if (req->min_x > req->max_x || req->min_y > req->max_y || req->min_z > req->max_z ||
        !finite(req->max_x - req->min_x) || !finite(req->max_y - req->min_y) ||
        !finite(req->max_z - req->min_z)) {
        error_text(error, capacity, "Invalid grid bounds"); return 0;
    }
    if (!req->nx || !req->ny || !req->nz) { error_text(error, capacity, "Invalid grid dimensions"); return 0; }
    const uint64_t nx = req->nx, ny = req->ny, nz = req->nz;
    if (nx > UINT64_MAX / ny || nx * ny > UINT64_MAX / nz) {
        error_text(error, capacity, "Grid dimensions overflow"); return 0;
    }
    const uint64_t total = nx * ny * nz;
    if (req->first_point > total || req->point_count > total - req->first_point) {
        error_text(error, capacity, "Grid chunk is outside the volume"); return 0;
    }
    for (uint64_t done = 0; done < req->point_count;) {
        const uint32_t count = static_cast<uint32_t>(std::min<uint64_t>(kChunkPoints, req->point_count - done));
        const size_t bytes = size_t(count) * sizeof(float);
        if (bytes > [ctx->device maxBufferLength]) {
            error_text(error, capacity, "Grid chunk exceeds Metal buffer limit"); return 0;
        }
        id<MTLBuffer> result = [ctx->device newBufferWithLength:bytes options:ctx->storage];
        if (!result) { error_text(error, capacity, "Could not allocate Metal output buffer"); return 0; }
        CovGridRequest chunk = *req;
        chunk.first_point += done;
        chunk.point_count = count;
        id<MTLCommandBuffer> command = [ctx->queue commandBuffer];
        if (!command) { error_text(error, capacity, "Could not create Metal command buffer"); return 0; }
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (!encoder) { error_text(error, capacity, "Could not create Metal compute encoder"); return 0; }
        [encoder setComputePipelineState:ctx->pipeline];
        [encoder setBuffer:ctx->terms offset:0 atIndex:0];
        [encoder setBytes:&ctx->term_count length:sizeof(ctx->term_count) atIndex:1];
        [encoder setBytes:&chunk length:sizeof(chunk) atIndex:2];
        [encoder setBuffer:result offset:0 atIndex:3];
        const NSUInteger width = std::min<NSUInteger>(ctx->pipeline.threadExecutionWidth,
                                                        ctx->pipeline.maxTotalThreadsPerThreadgroup);
        [encoder dispatchThreadgroups:MTLSizeMake((count + width - 1) / width, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
        [encoder endEncoding];
        if (ctx->storage == MTLResourceStorageModeManaged) {
            id<MTLBlitCommandEncoder> blit = [command blitCommandEncoder];
            if (!blit) { error_text(error, capacity, "Could not synchronize Metal output"); return 0; }
            [blit synchronizeResource:result];
            [blit endEncoding];
        }
        [command commit];
        [command waitUntilCompleted];
        if (command.status != MTLCommandBufferStatusCompleted) {
            error_text(error, capacity, ns_error(command.error, "Metal compute failed")); return 0;
        }
        std::memcpy(output + done, [result contents], bytes);
        done += count;
    }
    error_text(error, capacity, "");
    return 1;
}

void* create(int i, char* e, size_t n) noexcept {
    @autoreleasepool {
        try { return create_impl(i, e, n); }
        catch (const std::exception& ex) { error_text(e, n, ex.what()); return nullptr; }
        catch (...) { error_text(e, n, "Metal initialization failed"); return nullptr; }
    }
}
void destroy(void* c) noexcept { @autoreleasepool { try { delete static_cast<Context*>(c); } catch (...) {} } }
int set_terms(void* c, const CovGaussianTerm* t, uint32_t n, char* e, size_t cap) noexcept {
    @autoreleasepool {
        try { return set_terms_impl(c, t, n, e, cap); }
        catch (const std::exception& ex) { error_text(e, cap, ex.what()); return 0; }
        catch (...) { error_text(e, cap, "Metal term upload failed"); return 0; }
    }
}
int evaluate(void* c, const CovGridRequest* r, float* o, char* e, size_t cap) noexcept {
    @autoreleasepool {
        try { return evaluate_impl(c, r, o, e, cap); }
        catch (const std::exception& ex) { error_text(e, cap, ex.what()); return 0; }
        catch (...) { error_text(e, cap, "Metal evaluation failed"); return 0; }
    }
}
const char* device_name(void* c) noexcept {
    return c ? static_cast<Context*>(c)->name.c_str() : "";
}
const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "metal", create, destroy, set_terms, evaluate, device_name};
} // namespace

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
