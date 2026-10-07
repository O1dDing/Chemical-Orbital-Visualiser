#include "cov/compute_abi.h"

#include <algorithm>
#include <cstring>

namespace {
int fake_context = 1;

void message(char* error, size_t capacity, const char* text) noexcept {
    if (!error || !capacity) return;
    const size_t count = std::min(capacity - 1, std::strlen(text));
    std::memcpy(error, text, count);
    error[count] = 0;
}
void* create(int device_index, char* error, size_t capacity) noexcept {
    if (device_index != -1 && device_index != 0) {
        message(error, capacity, "fake HIP device index unavailable");
        return nullptr;
    }
    message(error, capacity, "");
    return &fake_context;
}
void destroy(void*) noexcept {}
int set_terms(void*, const CovGaussianTerm*, uint32_t, char* error, size_t capacity) noexcept {
    message(error, capacity, "");
    return 1;
}
int evaluate(void*, const CovGridRequest*, float*, char* error, size_t capacity) noexcept {
    message(error, capacity, "fake HIP evaluate failure");
    return 0;
}
const char* device_name(void*) noexcept { return "fake HIP GPU"; }
const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "hip", create, destroy,
                           set_terms, evaluate, device_name};
}

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
