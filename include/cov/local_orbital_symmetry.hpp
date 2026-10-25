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
#include "cov/local_angular_projection.hpp"
#include "cov/point_group_catalog.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace cov {

// Local s/p/d explanation, independent of all producer/whole-molecule labels.
// copy_index is retained because lower-symmetry groups can contain multiple
// independent copies carrying the same textual irrep label.  Zero means that
// the irrep label is determined but equivalent copies are internally mixed.
enum class LocalIrrepSource { MetricAngularProjection, DimensionCandidate };

struct LocalIrrepAssignment {
    std::string label;
    MetalAOShell shell = MetalAOShell::S;
    std::uint8_t copy_index = 1;
    double confidence = 0.0;
    std::string basis_functions;
    LocalIrrepSource source = LocalIrrepSource::DimensionCandidate;
    std::string point_group;
    double centre_mean_fraction = std::numeric_limits<double>::quiet_NaN();
    double angular_fraction_within_centre = std::numeric_limits<double>::quiet_NaN();
    double labelled_fraction_of_target = std::numeric_limits<double>::quiet_NaN();
    // Full target identities, alpha/beta blocks, centre, axes, reference sources,
    // ranks and residuals belong to this local explanation only.
    std::shared_ptr<const LocalAngularProjection> projection;
};

struct LocalIrrepAssessment {
    std::shared_ptr<const LocalAngularProjection> projection;
    std::optional<LocalIrrepAssignment> assignment;
    bool numerically_resolved = false;
    std::string detail;
};

// Retains the computed decomposition even when no single readable irrep wins.
// A failed/mixed numerical assessment is not replaced by dimensional certainty.
[[nodiscard]] LocalIrrepAssessment assess_local_metal_irrep(
    const LocalAngularProjectionWorkspace&,std::span<const std::size_t>,std::string_view point_group);

// A dimension-based candidate, without a numerical confidence. It succeeds only
// when all irrep copies with the supplied dimension carry one textual label;
// genuinely different same-dimensional irreps remain unassigned.
[[nodiscard]] std::optional<LocalIrrepAssignment>
classify_local_irrep_by_dimension(
    std::string_view point_group,
    MetalAOShell shell,
    std::size_t dimension);

// rotation_reference_to_input is a row-major orthogonal matrix mapping a
// vector in the point-group reference frame into the input molecular frame.
// The classifier projects actual MO subspaces through S onto joint radial
// angular references. Its confidence is conditional on the stated local shell,
// not evidence that the full MO transforms as this irrep. Repeated equal labels
// aggregate; unresolved copy identity is retained as zero.
[[nodiscard]] std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const Wavefunction& wavefunction,
    std::span<const std::size_t> orbital_indices,
    std::size_t metal_atom,
    std::string_view point_group,
    const std::array<double, 9>& rotation_reference_to_input);

// Reuse a workspace prepared outside the group loop.
[[nodiscard]] std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const LocalAngularProjectionWorkspace& workspace,
    std::span<const std::size_t> orbital_indices,
    std::string_view point_group);

} // namespace cov
