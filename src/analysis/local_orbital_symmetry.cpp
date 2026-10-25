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

#include "cov/local_orbital_symmetry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace cov {
namespace {
using ComponentWeights=std::array<double,5>;
constexpr double kWeightEpsilon=1.0e-14;
constexpr double kMinimumCopyConfidence=0.55;
std::optional<std::size_t> basis_component(
    const MetalAOShell shell,
    std::string_view token) noexcept {
    if (shell == MetalAOShell::S) return token == "s" ? std::optional<std::size_t>{0}
                                                       : std::nullopt;
    if (shell == MetalAOShell::P) {
        if (token == "p_x") return 0;
        if (token == "p_y") return 1;
        if (token == "p_z") return 2;
        return std::nullopt;
    }
    if (token == "d_z2") return 0;
    if (token == "d_xz") return 1;
    if (token == "d_yz") return 2;
    if (token == "d_x2-y2") return 3;
    if (token == "d_xy") return 4;
    return std::nullopt;
}

double score_basis_functions(const ComponentWeights& weights,
                             const MetalAOShell shell,
                             std::string_view functions) noexcept {
    double score = 0.0;
    while (!functions.empty()) {
        const std::size_t separator = functions.find(',');
        std::string_view token = functions.substr(0, separator);
        while (!token.empty() && token.front() == ' ') token.remove_prefix(1);
        while (!token.empty() && token.back() == ' ') token.remove_suffix(1);
        if (const auto component = basis_component(shell, token)) {
            score += weights[*component];
        }
        if (separator == std::string_view::npos) break;
        functions.remove_prefix(separator + 1);
    }
    return score;
}

struct ShellEvidence {
    MetalAOShell shell = MetalAOShell::S;
    ComponentWeights components{};
    double weight = 0.0;
    double fraction = 0.0;
};

double evidence_floor(const MetalAOShell shell) noexcept {
    return shell == MetalAOShell::D ? 0.02 : 0.08;
}

std::optional<LocalIrrepAssignment> classify_copy(
    const ShellEvidence& evidence,
    std::string_view point_group) {
    const auto decomposition = decompose_metal_ao_shell(point_group, evidence.shell);
    if (!decomposition || decomposition->empty()) return std::nullopt;

    struct ScoredCopy {
        const IrrepCopy* copy = nullptr;
        double score = 0.0;
    };
    struct ScoredLabel {
        std::string_view label;
        double score = 0.0;
        std::vector<ScoredCopy> copies;
    };
    std::vector<ScoredLabel> labels;
    double candidate_sum = 0.0;
    for (const auto& copy : *decomposition) {
        const double score = score_basis_functions(
            evidence.components, evidence.shell, copy.basis_functions);
        candidate_sum += score;
        auto label=std::find_if(labels.begin(),labels.end(),[&](const auto& item) {
            return item.label==copy.label;
        });
        if (label==labels.end()) {
            labels.push_back(ScoredLabel{copy.label,score,{{&copy,score}}});
        } else {
            label->score+=score;
            label->copies.push_back({&copy,score});
        }
    }
    if (labels.empty() || candidate_sum <= kWeightEpsilon) return std::nullopt;
    const auto best_label=std::max_element(
        labels.begin(),labels.end(),[](const auto& left,const auto& right) {
            return left.score<right.score;
        });
    const double confidence = best_label->score / candidate_sum;
    if (!std::isfinite(confidence) || confidence < kMinimumCopyConfidence) {
        return std::nullopt;
    }

    const auto best_copy=std::max_element(
        best_label->copies.begin(),best_label->copies.end(),
        [](const auto& left,const auto& right) {return left.score<right.score;});
    const double copy_confidence=best_label->score>kWeightEpsilon
        ?best_copy->score/best_label->score:0.0;
    const bool copy_determined=std::isfinite(copy_confidence) &&
        copy_confidence>=kMinimumCopyConfidence;
    return LocalIrrepAssignment{
        std::string(best_label->label),evidence.shell,
        copy_determined?best_copy->copy->copy_index:std::uint8_t{0},
        confidence,copy_determined
            ?std::string(best_copy->copy->basis_functions):std::string{}};
}

bool resolved_projection(const LocalAngularProjection& projection) {
    if(projection.status!=MetricSubspaceStatus::Available || projection.represented_spin_orbital_rank==0 ||
       !(projection.centre_projection_trace>kWeightEpsilon) || !std::isfinite(projection.centre_mean_fraction)) return false;
    for(const auto& spin:projection.spins) {
        const auto& centre=spin.centre;
        const double tolerance=512.0*std::numeric_limits<double>::epsilon()*
            static_cast<double>(std::max<std::size_t>(1,centre.basis_dimension))*
            static_cast<double>(std::max<std::size_t>(1,centre.orbitals.numerical_rank));
        if(centre.orbitals.numerical_rank!=spin.orbital_indices.size() ||
           !std::isfinite(spin.angular_partition_residual) || std::abs(spin.angular_partition_residual)>tolerance ||
           centre.metric_eigen_residual>tolerance || centre.reference.orthonormality_residual>tolerance ||
           centre.orbitals.orthonormality_residual>tolerance) return false;
        for(const auto& component:spin.components) {
            const auto& metric=component.metric;
            if(metric.status!=MetricSubspaceStatus::Available || !std::isfinite(metric.subspace_overlap_trace) ||
               metric.reference.orthonormality_residual>tolerance || metric.orbitals.orthonormality_residual>tolerance ||
               metric.reference.gram_inverse_residual>tolerance || metric.orbitals.gram_inverse_residual>tolerance) return false;
        }
    }
    return true;
}
} // namespace

std::optional<LocalIrrepAssignment> classify_local_irrep_by_dimension(
    const std::string_view point_group,
    const MetalAOShell shell,
    const std::size_t dimension) {
    if (dimension==0u) return std::nullopt;
    const auto decomposition=decompose_metal_ao_shell(point_group,shell);
    if (!decomposition) return std::nullopt;
    const IrrepCopy* match=nullptr;
    std::string_view label;
    bool repeated_copy=false;
    for (const auto& copy:*decomposition) {
        if (copy.dimension!=dimension) continue;
        if (match==nullptr) {
            match=&copy;
            label=copy.label;
        } else if (copy.label!=label) {
            return std::nullopt;
        } else {
            repeated_copy=true;
        }
    }
    if (match==nullptr) return std::nullopt;
    LocalIrrepAssignment assignment{
        std::string(match->label),shell,
        repeated_copy?std::uint8_t{0}:match->copy_index,1.0,
        repeated_copy?std::string{}:std::string(match->basis_functions)};
    assignment.source=LocalIrrepSource::DimensionCandidate;
    assignment.point_group=std::string(point_group);
    assignment.confidence=std::numeric_limits<double>::quiet_NaN();
    return assignment;
}


static std::optional<LocalIrrepAssignment> classify_projected_local_irrep(
    std::shared_ptr<const LocalAngularProjection> projection,
    std::string_view point_group) {
    if(!find_point_group(point_group))return std::nullopt;
    if(!resolved_projection(*projection))return std::nullopt;
    std::array<ShellEvidence,3> evidence{{
        {MetalAOShell::S,{},0,0},{MetalAOShell::P,{},0,0},{MetalAOShell::D,{},0,0}}};
    const auto& traces=projection->component_projection_traces;
    evidence[0].components[0]=traces[0];
    evidence[1].components={traces[2],traces[3],traces[1],0,0};
    for(std::size_t i=0;i<5;++i)evidence[2].components[i]=traces[4+i];
    for(auto& item:evidence) {
        for(double trace:item.components)item.weight+=trace;
        // The denominator includes all represented l=0..4 content, including
        // Cartesian lower-l directions. Nothing is renormalized away as s/p/d.
        item.fraction=item.weight/projection->centre_projection_trace;
    }
    std::stable_sort(evidence.begin(),evidence.end(),[](const auto& a,const auto& b){return a.fraction>b.fraction;});
    for(const auto& item:evidence) {
        if(item.fraction+kWeightEpsilon<evidence_floor(item.shell))continue;
        auto assignment=classify_copy(item,point_group);
        if(!assignment)continue;
        assignment->source=LocalIrrepSource::MetricAngularProjection;
        assignment->point_group=std::string(point_group);
        assignment->centre_mean_fraction=projection->centre_mean_fraction;
        assignment->angular_fraction_within_centre=item.fraction;
        assignment->labelled_fraction_of_target=assignment->confidence*item.weight/
            static_cast<double>(projection->represented_spin_orbital_rank);
        assignment->projection=std::move(projection);
        return assignment;
    }
    return std::nullopt;
}

LocalIrrepAssessment assess_local_metal_irrep(
    const LocalAngularProjectionWorkspace& workspace,std::span<const std::size_t> orbital_indices,
    std::string_view point_group) {
    LocalIrrepAssessment result;
    result.projection=std::make_shared<const LocalAngularProjection>(workspace.project(orbital_indices));
    result.numerically_resolved=resolved_projection(*result.projection);
    result.assignment=classify_projected_local_irrep(result.projection,point_group);
    if (result.assignment) result.detail="Local angular decomposition supports the displayed conditional label";
    else if (!find_point_group(point_group)) result.detail="Local decomposition retained; point-group catalogue mapping is unavailable";
    else if (!result.numerically_resolved) result.detail="Local decomposition retained with its zero, rank or numerical diagnostic; no irrep assigned";
    else result.detail="Local decomposition is mixed or outside the readable s/p/d label model; no single irrep assigned";
    return result;
}

std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const LocalAngularProjectionWorkspace& workspace,std::span<const std::size_t> orbital_indices,
    std::string_view point_group) {
    return assess_local_metal_irrep(workspace,orbital_indices,point_group).assignment;
}

std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const Wavefunction& wavefunction,std::span<const std::size_t> orbital_indices,
    std::size_t atom,std::string_view point_group,
    const std::array<double,9>& rotation_reference_to_input) {
    if(!find_point_group(point_group))return std::nullopt;
    const LocalAngularProjectionWorkspace workspace(wavefunction,atom,rotation_reference_to_input);
    return classify_local_metal_irrep(workspace,orbital_indices,point_group);
}
} // namespace cov
