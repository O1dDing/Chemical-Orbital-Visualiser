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

namespace cov {

struct OrbitalChemistryOptions {
    double degeneracy_tolerance_hartree = 1.0e-5;
    double pair_distance_ceiling_bohr = 6.0;
    double pair_mayer_floor = 0.02;
    double pair_atom_weight_floor = 0.08;
    double overlap_bonding_threshold = 0.015;
    double nonbonding_overlap_threshold = 0.004;
    double determined_fraction = 0.82;
    double maximum_undetermined_for_percentages = 0.55;
};

// Builds a molecule-independent, element/ECP-aware minimal chemical-valence
// reference inside the working AO basis, projects every canonical MO onto that
// reference with the AO-overlap metric, selects the canonical valence manifold,
// and derives per-interaction chemistry that FCHK can support.
//
// This is COV's own FCHK-derived reference analysis. It is not called IAO,
// NAO, NBO, Wiberg or producer data.
void derive_orbital_chemistry(
    Wavefunction& wavefunction,
    const OrbitalChemistryOptions& options = {});

[[nodiscard]] const char* chemistry_status_name(ChemistryStatus value) noexcept;
[[nodiscard]] const char* orbital_angular_family_name(
    OrbitalAngularFamily value) noexcept;
[[nodiscard]] const char* orbital_bonding_role_name(
    OrbitalBondingRole value) noexcept;

} // namespace cov
