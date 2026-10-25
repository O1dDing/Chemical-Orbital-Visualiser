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

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cov {

// phi_internal[i] = basis_scale * phi_Gaussian[source_index].  Gaussian and
// COV use the same primitive normalization; the scale here is a phase (+/-1).
// Keep coefficient_scale explicit so a future normalization change cannot
// accidentally apply a covariant matrix transform to contravariant data.
struct GaussianAoTransformEntry {
    std::size_t source_index = 0;
    double basis_scale = 1.0;
    double coefficient_scale = 1.0;
};

inline std::vector<GaussianAoTransformEntry> gaussian_ao_transform(
    const int shell_type, const std::size_t source_offset = 0) {
    const int l = shell_type < 0 ? -shell_type : shell_type;
    if (l > 4) throw std::runtime_error("Gaussian AO transform above g is unavailable");
    const bool pure = shell_type <= -2;
    const std::size_t count = shell_type == -1 ? 4u : pure
        ? static_cast<std::size_t>(2*l+1)
        : static_cast<std::size_t>((l+1)*(l+2)/2);
    static constexpr std::size_t cartesian_g[15] = {
        14,4,0,13,12,8,3,5,1,11,9,2,10,7,6};
    std::vector<GaussianAoTransformEntry> result;
    result.reserve(count);
    for (std::size_t i=0; i<count; ++i) {
        // Real solid harmonics use m=0,+1,-1,+2,-2,... . COV's
        // polynomials include the Condon--Shortley (-1)^|m| phase;
        // Gaussian FCHK pure shells omit that phase.
        const std::size_t m = (i+1)/2;
        const double phase = pure && m%2 != 0 ? -1.0 : 1.0;
        const std::size_t source = !pure && l==4 ? cartesian_g[i] : i;
        result.push_back({source_offset+source, phase, 1.0/phase});
    }
    return result;
}

} // namespace cov
