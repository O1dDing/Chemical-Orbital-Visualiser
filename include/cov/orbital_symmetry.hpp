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

#include "cov/model.hpp"

#include <cstddef>
#include <string>

namespace cov {

struct OrbitalSymmetryOptions {
    double degeneracy_tolerance_hartree = 1.0e-5;
    double character_tolerance = 0.18;
    double minimum_subspace_retention = 0.985;
};

struct OrbitalSymmetryResult {
    std::string point_group;
    std::size_t groups_examined = 0;
    std::size_t groups_labelled = 0;
    std::size_t orbitals_labelled = 0;
    double worst_subspace_retention = 1.0;

    [[nodiscard]] bool applied() const noexcept {
        return orbitals_labelled != 0;
    }
};

// Derive conventional MO irreps from COV's own validated molecular symmetry
// operations. The implementation transforms the actual AO basis, evaluates
// each degenerate MO subspace with the AO-overlap metric, and only assigns a
// label when the resulting characters match a supported point-group irrep.
// Producer-reported MO labels are never overwritten.
OrbitalSymmetryResult derive_orbital_symmetry(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options = {});

} // namespace cov
