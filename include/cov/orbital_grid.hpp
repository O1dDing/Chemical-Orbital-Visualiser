#pragma once

#include "cov/compute_abi.h"
#include "cov/model.hpp"
#include <span>
#include <stop_token>
#include <vector>

namespace cov {

// This is a rendering representation. Wavefunction scientific data stay double.
std::vector<CovGaussianTerm> pack_orbital_terms(const Wavefunction& wavefunction,
                                               std::size_t mo_index);
std::size_t checked_grid_size(int nx, int ny, int nz);
CovGridRequest make_grid_request(const GridBox& box, int nx, int ny, int nz,
                                std::size_t first, std::size_t count);
void evaluate_cpu_grid(std::span<const CovGaussianTerm> terms,
                       const CovGridRequest& request, std::span<float> output,
                       std::stop_token stop = {});

}
