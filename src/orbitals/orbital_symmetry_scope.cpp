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

#include "cov/orbital_symmetry_scope.hpp"
#include "cov/ligand_field.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cov {
const char* orbital_symmetry_origin_name(OrbitalSymmetryOrigin origin) noexcept {
    switch (origin) {
        case OrbitalSymmetryOrigin::Producer: return "producer";
        case OrbitalSymmetryOrigin::MolecularOperations: return "molecular-operations";
        case OrbitalSymmetryOrigin::LocalMetricProjection: return "local-metric-projection";
        case OrbitalSymmetryOrigin::LocalDimensionCandidate: return "local-dimension-candidate";
        case OrbitalSymmetryOrigin::PiPartnerCandidate: return "pi-partner-candidate";
        case OrbitalSymmetryOrigin::SpinCounterpartCandidate: return "spin-counterpart-candidate";
        case OrbitalSymmetryOrigin::MixedMembers: return "mixed-members";
        case OrbitalSymmetryOrigin::LocalMetricDecomposition: return "local-metric-decomposition";
        default: return "unavailable";
    }
}
bool orbital_symmetry_is_candidate(const OrbitalSymmetryExplanation& value) noexcept {
    return value.origin == OrbitalSymmetryOrigin::LocalDimensionCandidate ||
           value.origin == OrbitalSymmetryOrigin::PiPartnerCandidate ||
           value.origin == OrbitalSymmetryOrigin::SpinCounterpartCandidate;
}
bool orbital_symmetry_is_local(const OrbitalSymmetryExplanation& value) noexcept {
    return value.origin == OrbitalSymmetryOrigin::LocalMetricProjection ||
           value.origin == OrbitalSymmetryOrigin::LocalMetricDecomposition ||
           value.origin == OrbitalSymmetryOrigin::LocalDimensionCandidate ||
           value.origin == OrbitalSymmetryOrigin::PiPartnerCandidate ||
           (value.origin == OrbitalSymmetryOrigin::SpinCounterpartCandidate &&
            value.candidate_source && orbital_symmetry_is_local(*value.candidate_source));
}
OrbitalSymmetryExplanation molecular_orbital_symmetry(const Wavefunction& wf, std::size_t index) {
    OrbitalSymmetryExplanation result;
    result.orbital_indices = {index};
    if (index >= wf.orbitals.size()) return result;
    const auto& mo = wf.orbitals[index];
    result.label = mo.symmetry;
    if (result.label.empty()) return result;
    if (mo.symmetry_provenance == DataProvenance::Producer) {
        result.origin = OrbitalSymmetryOrigin::Producer;
        result.source_path = wf.enrichment_source;
        for (const auto& record : wf.orbital_symmetry_source_records) {
            const auto found = std::find(record.orbital_indices.begin(), record.orbital_indices.end(), index);
            if (found == record.orbital_indices.end()) continue;
            const auto offset = static_cast<std::size_t>(found - record.orbital_indices.begin());
            if (offset >= record.labels.size() || record.labels[offset] != mo.symmetry) continue;
            result.source_path = record.source_path;
            result.source_line_begin = record.line_begin;
            result.source_line_end = record.line_end;
            result.source_job_segment = record.job_segment;
            result.producer_detected_group = record.detected_group_context;
            result.producer_abelian_group = record.abelian_group_context;
        }
    } else if (mo.symmetry_provenance == DataProvenance::Derived) {
        result.origin = OrbitalSymmetryOrigin::MolecularOperations;
        for (const auto& assignment : wf.derived_orbital_symmetry_assignments) {
            if (assignment.label != mo.symmetry ||
                std::find(assignment.orbital_indices.begin(), assignment.orbital_indices.end(), index) ==
                    assignment.orbital_indices.end()) continue;
            result.point_group = assignment.point_group;
            result.point_group_basis = "derived-molecular-operations";
            result.molecular_assignment = std::make_shared<const DerivedOrbitalSymmetryAssignment>(assignment);
            break;
        }
    }
    return result;
}
OrbitalSymmetryExplanation local_orbital_symmetry(const LocalIrrepAssignment& assignment,
    const LigandFieldEnvironment& environment, std::span<const std::size_t> indices) {
    OrbitalSymmetryExplanation result;
    result.label = assignment.label;
    result.origin = assignment.source == LocalIrrepSource::MetricAngularProjection
        ? OrbitalSymmetryOrigin::LocalMetricProjection : OrbitalSymmetryOrigin::LocalDimensionCandidate;
    result.point_group = assignment.point_group;
    result.point_group_basis = "local-coordination-template";
    result.orbital_indices.assign(indices.begin(), indices.end());
    result.atom_indices = {environment.metal_atom};
    result.atom_indices.insert(result.atom_indices.end(), environment.ligand_atoms.begin(), environment.ligand_atoms.end());
    result.rotation_reference_to_input = environment.rotation_reference_to_input;
    result.axes_available = true;
    result.geometry_id = environment.geometry_machine_id();
    result.geometry_angular_rms = environment.angular_rms;
    result.geometry_shape_measure = environment.shape_measure;
    result.local_assignment = std::make_shared<const LocalIrrepAssignment>(assignment);
    result.local_decomposition = assignment.projection;
    return result;
}
bool orbital_symmetry_missing_local_input(const OrbitalSymmetryExplanation& value) noexcept {
    return value.local_decomposition &&
        value.local_decomposition->status == MetricSubspaceStatus::MissingInput;
}
OrbitalSymmetryExplanation evaluate_local_orbital_symmetry(
    const LocalAngularProjectionWorkspace& workspace,const LigandFieldEnvironment& environment,
    std::span<const std::size_t> indices,int family) {
    const auto assessed=assess_local_metal_irrep(workspace,indices,environment.local_point_group());
    OrbitalSymmetryExplanation result;
    if (assessed.assignment) result=local_orbital_symmetry(*assessed.assignment,environment,indices);
    else {
        // Only absent numerical input permits the catalogue-only candidate.
        // Existing mixed/zero/invalid evidence must not be overridden by count.
        std::optional<LocalIrrepAssignment> candidate;
        if (assessed.projection->status==MetricSubspaceStatus::MissingInput && family>=0 && family<=2)
            candidate=classify_local_irrep_by_dimension(environment.local_point_group(),
                static_cast<MetalAOShell>(family),indices.size());
        if (candidate) result=local_orbital_symmetry(*candidate,environment,indices);
        else {
            LocalIrrepAssignment scope;
            scope.point_group=environment.local_point_group();
            result=local_orbital_symmetry(scope,environment,indices);
            result.origin=OrbitalSymmetryOrigin::LocalMetricDecomposition;
            result.local_assignment.reset();
            result.label=assessed.numerically_resolved?"mixed":"UND";
        }
    }
    result.local_decomposition=assessed.projection;
    result.assessment_detail=assessed.detail;
    return result;
}
OrbitalSymmetryExplanation candidate_orbital_symmetry(const OrbitalSymmetryExplanation& source,
    OrbitalSymmetryOrigin origin, std::span<const std::size_t> indices, double score) {
    if (origin != OrbitalSymmetryOrigin::PiPartnerCandidate &&
        origin != OrbitalSymmetryOrigin::SpinCounterpartCandidate)
        throw std::invalid_argument("A propagated label must retain its candidate origin");
    OrbitalSymmetryExplanation result = source;
    result.origin = origin;
    result.orbital_indices.assign(indices.begin(), indices.end());
    result.local_assignment.reset();
    result.local_decomposition.reset();
    result.assessment_detail.clear();
    result.pi_partner_evidence.reset();
    result.molecular_assignment.reset();
    result.source_path.clear();
    result.source_line_begin = result.source_line_end = result.source_job_segment = 0;
    result.candidate_source = std::make_shared<const OrbitalSymmetryExplanation>(source);
    result.candidate_score = score;
    return result;
}
bool compatible_symmetry_scopes(const OrbitalSymmetryExplanation& a,
                               const OrbitalSymmetryExplanation& b) noexcept {
    if (a.label.empty() || b.label.empty() || a.origin == OrbitalSymmetryOrigin::Unavailable ||
        b.origin == OrbitalSymmetryOrigin::Unavailable || a.origin == OrbitalSymmetryOrigin::MixedMembers ||
        b.origin == OrbitalSymmetryOrigin::MixedMembers) return false;
    if (orbital_symmetry_is_local(a) != orbital_symmetry_is_local(b)) return false;
    if (a.point_group != b.point_group || a.point_group_basis != b.point_group_basis ||
        a.axes_available != b.axes_available || a.atom_indices != b.atom_indices) return false;
    if (a.axes_available && a.rotation_reference_to_input != b.rotation_reference_to_input) return false;
    if (orbital_symmetry_is_local(a)) return !a.point_group.empty();
    if (a.origin != b.origin) return false;
    if (a.origin == OrbitalSymmetryOrigin::Producer) {
        return a.source_path == b.source_path && a.source_line_begin == b.source_line_begin &&
               a.source_line_end == b.source_line_end && a.source_job_segment == b.source_job_segment &&
               a.producer_detected_group == b.producer_detected_group && a.producer_abelian_group == b.producer_abelian_group;
    }
    if (a.origin == OrbitalSymmetryOrigin::MolecularOperations) {
        if (!a.molecular_assignment || !b.molecular_assignment) return false;
        const auto& x = *a.molecular_assignment; const auto& y = *b.molecular_assignment;
        return x.centre_bohr == y.centre_bohr && x.axes_available == y.axes_available &&
               x.principal_axis == y.principal_axis && x.secondary_axis == y.secondary_axis &&
               x.axis_convention == y.axis_convention;
    }
    return false;
}
OrbitalSymmetryExplanation aggregate_molecular_symmetry(const Wavefunction& wf,
                                                       std::span<const std::size_t> indices) {
    if (indices.empty()) return {};
    auto result = molecular_orbital_symmetry(wf, indices.front());
    bool uniform = true;
    for (const auto index : indices.subspan(1)) {
        const auto other = molecular_orbital_symmetry(wf, index);
        if (result.label != other.label || (!result.label.empty() && !compatible_symmetry_scopes(result, other)))
            uniform = false;
    }
    if (!uniform) {
        result = {};
        result.label = "mixed";
        result.origin = OrbitalSymmetryOrigin::MixedMembers;
    }
    result.orbital_indices.assign(indices.begin(), indices.end());
    return result;
}
std::string orbital_symmetry_compact_text(const OrbitalSymmetryExplanation& value) {
    if (value.label.empty()) return "N/A";
    const std::string group = value.point_group.empty() ? std::string{} : ":" + value.point_group;
    if (orbital_symmetry_is_candidate(value)) return value.label + " [candidate" + group + "]";
    if (orbital_symmetry_is_local(value)) return value.label + " [local" + group + "]";
    return value.label;
}
} // namespace cov
