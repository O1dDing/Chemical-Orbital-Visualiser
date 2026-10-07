#pragma once

/* Private, versioned boundary for separately compiled GPU runtimes. */
#include <stddef.h>
#include <stdint.h>

#define COV_COMPUTE_ABI_VERSION 1u
#ifdef _WIN32
#define COV_COMPUTE_EXPORT __declspec(dllexport)
#else
#define COV_COMPUTE_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* One normalized Cartesian Gaussian monomial, including the MO coefficient.
 * Exactly 32 bytes; scalar fields keep the same layout in C++, MSL and WGSL. */
typedef struct CovGaussianTerm {
    float cx, cy, cz, exponent;
    float coefficient;
    uint32_t ax, ay, az;
} CovGaussianTerm;

typedef struct CovGridRequest {
    float min_x, min_y, min_z;
    float max_x, max_y, max_z;
    uint32_t nx, ny, nz;
    uint64_t first_point;
    uint32_t point_count;
} CovGridRequest;

/* The owner keeps term data valid until set_terms returns. evaluate writes
 * exactly point_count float values, in x-fastest order, before returning.
 * Device index -1 chooses the first suitable GPU. No exception crosses ABI.
 * Error buffers are NUL-terminated whenever their capacity is nonzero. */
typedef struct CovComputeApi {
    uint32_t abi_version;
    const char* backend_name;
    void* (*create)(int device_index, char* error, size_t error_capacity);
    void (*destroy)(void* context);
    int (*set_terms)(void* context, const CovGaussianTerm* terms, uint32_t count,
                     char* error, size_t error_capacity);
    int (*evaluate)(void* context, const CovGridRequest* request, float* output,
                    char* error, size_t error_capacity);
    const char* (*device_name)(void* context);
} CovComputeApi;

typedef const CovComputeApi* (*CovGetComputeApi)(uint32_t version);

#ifdef __cplusplus
}
static_assert(sizeof(CovGaussianTerm) == 32);
#endif
