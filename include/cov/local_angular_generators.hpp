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
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace cov {
struct LocalAngularGenerators {
    bool available = false;
    int shell_degree = -1;
    int angular_degree = -1;
    bool pure_shell = false;
    std::size_t basis_count = 0;
    std::size_t columns = 0;
    // Row-major AO coefficients, one generator for each real angular function.
    // The order is m=0, real m=1, imaginary m=1, real m=2, ... with Condon--Shortley phase.
    std::vector<double> coefficients;
    double rotation_orthogonality_error = 0;
    double polynomial_reconstruction_error = 0;
    std::string detail;
};

// Functions are defined in the supplied local frame. Cartesian degree L also
// contains r^(L-l) Y_lm for l=L-2,L-4,...; these directions remain explicit.
// Radial contractions and overlap normalization are handled by the caller's
// common AO metric, not by independent coefficient-squared weights.
LocalAngularGenerators local_angular_generators(
    int shell_degree, bool pure_shell, int angular_degree,
    const std::array<double,9>& rotation_reference_to_input);
}
