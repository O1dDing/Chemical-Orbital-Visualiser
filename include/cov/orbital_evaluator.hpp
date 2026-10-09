#pragma once

#include "cov/model.hpp"
#include <memory>
#include <string>

namespace cov {

struct ComputeOptions {
    std::string backend = "auto";
    int device_index = -1;
};

// Only the UI thread touches OpenGL. Portable compute runs in a cancellable
// worker and publishes a complete grid before the texture is replaced.
class OrbitalEvaluator {
public:
    OrbitalEvaluator(const Wavefunction&, ComputeOptions options = {});
    ~OrbitalEvaluator();
    OrbitalEvaluator(const OrbitalEvaluator&) = delete;
    OrbitalEvaluator& operator=(const OrbitalEvaluator&) = delete;
    void attach_gl_texture(unsigned int texture);
    void detach_gl_texture();
    void evaluate(std::size_t mo, const GridBox&, int nx, int ny, int nz);
    void begin_evaluate(std::size_t mo, const GridBox&, int nx, int ny, int nz);
    bool poll();
    void cancel();
    bool busy() const noexcept;
    bool ready() const noexcept;
    const char* device_name() const noexcept;
    double last_kernel_ms() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cov
