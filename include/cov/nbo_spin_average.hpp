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
#include "cov/nbo_salc.hpp"

namespace cov {
struct NboRoCommonEnergyOrbital {
    std::size_t canonical_index=0, source_orbital_index=0; // both zero based
    Spin spin=Spin::Alpha;
    double source_energy_hartree=0, occupation=0;
    std::optional<double> alpha_energy_hartree, beta_energy_hartree, common_energy_hartree;
    bool available=false;
    std::string status="unavailable";
};
struct NboRoCommonEnergyModel {
    bool available=false;
    // Association proves numerical compatibility, not the original SCF job's
    // operator identity. No current manifest verifier certifies the latter.
    bool source_step_identity_verified=false;
    std::string status="unavailable", detail, dataset_id, canonical_fingerprint;
    std::string operator_semantics="associated_archive_spin_average_expectation";
    std::string provenance_status;
    NboRestrictedOpenShellEvidence restricted_open_shell;
    std::vector<NboSalcEnergyEvidence> spin_evidence;
    std::vector<NboRoCommonEnergyOrbital> orbitals;
    double canonical_linkage_error=0, metric_error=0, projection_residual=0;
    double hermiticity_error=0, density_error=0, maximum_offdiagonal_hartree=0;
};
// Evaluates <source MO|(F_alpha+F_beta)/2|source MO>. These are
// expectations, not re-canonicalized eigenvalues. Inputs are never modified.
NboRoCommonEnergyModel build_nbo_ro_common_energy(
    const Wavefunction&, const NboIntegration&, const NboSalcModel& raw);
std::string serialize_nbo_ro_common_energy_json(const NboRoCommonEnergyModel&);
} // namespace cov
