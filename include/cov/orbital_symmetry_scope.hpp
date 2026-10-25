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
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace cov {
struct LocalIrrepAssignment;
struct PiPartnerAssessment;
struct LocalAngularProjection;
class LocalAngularProjectionWorkspace;
struct LigandFieldEnvironment;

enum class OrbitalSymmetryOrigin {
    Unavailable,
    Producer,
    MolecularOperations,
    LocalMetricProjection,
    LocalDimensionCandidate,
    PiPartnerCandidate,
    SpinCounterpartCandidate,
    MixedMembers,
    LocalMetricDecomposition,
};

// A label and its scope form one value. A local/template/candidate label is
// never a replacement for the full-orbital label stored on MolecularOrbital.
// Candidate sources are immutable snapshots; target MO identities remain here.
struct OrbitalSymmetryExplanation {
    std::string label;
    OrbitalSymmetryOrigin origin = OrbitalSymmetryOrigin::Unavailable;
    std::string point_group;
    std::string point_group_basis = "unresolved";
    std::vector<std::size_t> orbital_indices;
    std::vector<std::size_t> atom_indices;
    std::array<double,9> rotation_reference_to_input{};
    bool axes_available = false;
    std::string geometry_id;
    double geometry_angular_rms = std::numeric_limits<double>::quiet_NaN();
    double geometry_shape_measure = std::numeric_limits<double>::quiet_NaN();
    std::shared_ptr<const LocalIrrepAssignment> local_assignment;
    std::shared_ptr<const LocalAngularProjection> local_decomposition;
    std::string assessment_detail;
    std::shared_ptr<const PiPartnerAssessment> pi_partner_evidence;
    std::shared_ptr<const DerivedOrbitalSymmetryAssignment> molecular_assignment;
    std::shared_ptr<const OrbitalSymmetryExplanation> candidate_source;
    double candidate_score = std::numeric_limits<double>::quiet_NaN();
    std::string source_path;
    std::size_t source_line_begin = 0;
    std::size_t source_line_end = 0;
    std::size_t source_job_segment = 0;
    // Context reported by the producer is kept without assuming that a printed
    // MO label belongs to either the detected group or its Abelian subgroup.
    std::string producer_detected_group;
    std::string producer_abelian_group;
};

[[nodiscard]] const char* orbital_symmetry_origin_name(OrbitalSymmetryOrigin) noexcept;
[[nodiscard]] bool orbital_symmetry_is_local(const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] bool orbital_symmetry_is_candidate(const OrbitalSymmetryExplanation&) noexcept;
// Numerical absence is distinct from a computed mixed, zero or invalid span.
[[nodiscard]] bool orbital_symmetry_missing_local_input(const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] OrbitalSymmetryExplanation molecular_orbital_symmetry(
    const Wavefunction&, std::size_t orbital_index);
[[nodiscard]] OrbitalSymmetryExplanation local_orbital_symmetry(
    const LocalIrrepAssignment&, const LigandFieldEnvironment&,
    std::span<const std::size_t> orbital_indices);
[[nodiscard]] OrbitalSymmetryExplanation evaluate_local_orbital_symmetry(
    const LocalAngularProjectionWorkspace&,const LigandFieldEnvironment&,
    std::span<const std::size_t> orbital_indices,int admitted_shell_family);
[[nodiscard]] OrbitalSymmetryExplanation candidate_orbital_symmetry(
    const OrbitalSymmetryExplanation& source, OrbitalSymmetryOrigin,
    std::span<const std::size_t> target_indices, double score);
[[nodiscard]] bool compatible_symmetry_scopes(
    const OrbitalSymmetryExplanation&, const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] OrbitalSymmetryExplanation aggregate_molecular_symmetry(
    const Wavefunction&, std::span<const std::size_t> orbital_indices);
// Human-facing compact text always marks local and candidate explanations.
// More detailed values use the same explanation object in UI and exports.
[[nodiscard]] std::string orbital_symmetry_compact_text(const OrbitalSymmetryExplanation&);
[[nodiscard]] std::string orbital_symmetry_json(const OrbitalSymmetryExplanation&);
}
