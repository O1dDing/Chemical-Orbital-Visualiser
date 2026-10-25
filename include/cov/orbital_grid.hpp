/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#pragma once

#include "cov/compute_abi.h"
#include "cov/model.hpp"
#include "cov/threading.hpp"
#include <span>
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
                       cov::stop_token stop = {});

}
