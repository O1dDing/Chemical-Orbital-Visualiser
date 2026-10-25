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
#include "cov/interaction_graph.hpp"
#include <utility>

namespace cov {
// A fixed set of atom pairs describes one overlap-analysis scope. It is never
// selected from a particular MO's atom weights and is not a bond-energy model.
struct FixedBondingScope {
    std::string id="whole_skeleton",source,status="unavailable";
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    std::size_t atom_count=0,invalid_edges=0;
};
enum class OrbitalGroupBondingStatus {
    Positive,Negative,Mixed,Unresolved,NotApplicable,Unavailable,IncompleteGroup
};
struct OrbitalGroupBondingOptions {
    // Zero means the caller supplied the complete member list. A known group
    // dimension should be passed explicitly to reject an incomplete list.
    std::size_t expected_dimension=0;
    double metric_rank_relative_tolerance=1e-10;
    double absolute_error=1e-7,relative_error=1e-5;
};
struct OrbitalGroupBondingResult {
    OrbitalGroupBondingStatus status=OrbitalGroupBondingStatus::Unavailable;
    std::string scope_id,scope_source,detail;
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    std::vector<std::size_t> source_members;
    std::size_t dimension=0,metric_rank=0;
    // Row-major matrices in a derived orthonormal group basis. Original source
    // coefficients, labels, energies, occupations and ordering remain intact.
    std::vector<double> source_metric,normalization,operator_matrix;
    std::vector<double> occupation_matrix,eigenvalues;
    bool occupation_available=false;
    double metric_error=0,hermiticity_error=0,eigensolver_residual=0;
    double trace=0,electron_weighted_trace=0,electron_count=0,error_bound=0;
    double minimum_eigenvalue=0,maximum_eigenvalue=0;
    double source_energy_span_hartree=0;
};
[[nodiscard]] FixedBondingScope make_fixed_bonding_scope(const Wavefunction&);
[[nodiscard]] FixedBondingScope make_fixed_bonding_scope(
    const Wavefunction&,const InteractionGraph&);
[[nodiscard]] FixedBondingScope make_fixed_bonding_scope(
    const Wavefunction&,const std::vector<std::pair<std::size_t,std::size_t>>&,
    const std::string& scope_id);
// Source members must all belong to one physical spin block. RO spatial copies
// are analysed in their original alpha/beta blocks, never stacked as duplicate
// columns. A supplied occupation matrix is transformed with the group metric.
[[nodiscard]] OrbitalGroupBondingResult analyse_orbital_group_bonding(
    const Wavefunction&,const FixedBondingScope&,
    const std::vector<std::size_t>& complete_members,
    const OrbitalGroupBondingOptions& options={},
    const std::vector<double>& occupation_matrix={});
[[nodiscard]] const char* orbital_group_bonding_status_name(
    OrbitalGroupBondingStatus) noexcept;
} // namespace cov
