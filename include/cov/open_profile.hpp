#pragma once
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace cov {
// Optional wall-clock checkpoints. They do not alter parsing or UI state.
class OpenProfile {
    bool enabled_=std::getenv("COV_PROFILE_OPEN")!=nullptr;
    std::chrono::steady_clock::time_point previous_=std::chrono::steady_clock::now();
public:
    void stage(const char* label) {
        if(!enabled_)return;
        const auto now=std::chrono::steady_clock::now();
        std::fprintf(stderr,"COV open detail: %s seconds=%.6f\n",label,
            std::chrono::duration<double>(now-previous_).count());
        std::fflush(stderr);previous_=now;
    }
};
}
