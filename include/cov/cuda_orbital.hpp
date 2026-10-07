#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace cov {

class CudaOrbitalEvaluator {
public:
    explicit CudaOrbitalEvaluator(const Wavefunction& wavefunction, int device_index = -1);
    ~CudaOrbitalEvaluator();

    CudaOrbitalEvaluator(const CudaOrbitalEvaluator&) = delete;
    CudaOrbitalEvaluator& operator=(const CudaOrbitalEvaluator&) = delete;
    CudaOrbitalEvaluator(CudaOrbitalEvaluator&&) noexcept;
    CudaOrbitalEvaluator& operator=(CudaOrbitalEvaluator&&) noexcept;

    void attach_gl_texture(unsigned int texture);
    void detach_gl_texture();

    void evaluate(std::size_t mo_index,
                  const GridBox& box,
                  int nx,
                  int ny,
                  int nz);

    [[nodiscard]] const char* device_name() const noexcept;
    [[nodiscard]] double last_kernel_ms() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cov
