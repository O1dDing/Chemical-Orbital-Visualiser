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

#include <map>

namespace cov {

struct BondAnalysisOptions {
    // Retain weak positive/negative pairwise indices for diagnostics. Rendering
    // applies its own stricter connectivity threshold.
    double record_abs_mayer_floor = 0.02;

    // Multicentre candidates are participation diagnostics only. Atomic weights
    // are absolute Mulliken-like MO gross populations normalized over atoms.
    double multicentre_atom_participation = 0.12;
    double multicentre_min_covered_weight = 0.72;

    // Strict 3-centre active-subspace classification. These thresholds are
    // intentionally conservative: no assignment is emitted unless a seeded
    // three-atom set has a unique three-orbital active subspace, at least two
    // electronically supported pair connections, and a clean 2/0/0 or 2/2/0
    // restricted-occupation pattern.
    double assignment_triple_coverage = 0.88;
    double assignment_three_atom_floor = 0.08;
    double assignment_pair_atom_floor = 0.18;
    double assignment_pair_third_ceiling = 0.08;
    double assignment_min_mayer_support = 0.06;
    double assignment_occupation_tolerance = 0.15;
    double assignment_min_uniqueness_margin = 0.05;
    double assignment_min_energy_hartree = -3.0;
    double assignment_max_energy_hartree = 3.0;
};

// Requires wavefunction.ao_overlap and total density. Populates pairwise Mayer
// bond orders, occupied-MO multicentre participation candidates, and only when
// the stricter active-subspace evidence is unique, provenance-aware 3c2e/3c4e
// assignments. No molecule-specific templates are used.
void derive_bond_and_multicentre_analysis(
    Wavefunction& wavefunction,
    const BondAnalysisOptions& options = {});

} // namespace cov
