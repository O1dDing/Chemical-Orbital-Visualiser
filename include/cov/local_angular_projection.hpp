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
#include "cov/metric_subspace.hpp"
#include "cov/model.hpp"
#include <array>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace cov {
struct LocalAngularProjectionComponent {
    int angular_degree = 0;
    std::size_t component = 0; // Real spherical order, including its phase convention.
    std::vector<std::size_t> source_shell_indices;
    MetricSubspaceOverlap metric;
    double fraction_of_local_projection = std::numeric_limits<double>::quiet_NaN();
};
struct LocalAngularSpinProjection {
    Spin spin = Spin::Alpha;
    std::vector<std::size_t> orbital_indices;
    MetricSubspaceOverlap centre;
    std::vector<LocalAngularProjectionComponent> components;
    double angular_partition_residual = std::numeric_limits<double>::quiet_NaN();
};
struct LocalAngularProjection {
    MetricSubspaceStatus status = MetricSubspaceStatus::MissingInput;
    std::size_t atom_index = 0;
    std::array<double,9> rotation_reference_to_input{};
    std::vector<std::size_t> source_shell_indices;
    std::vector<LocalAngularSpinProjection> spins;
    std::size_t represented_spin_orbital_rank = 0;
    double centre_projection_trace = std::numeric_limits<double>::quiet_NaN();
    double centre_mean_fraction = std::numeric_limits<double>::quiet_NaN();
    std::array<double,25> component_projection_traces = [] {
        std::array<double,25> values;
        values.fill(std::numeric_limits<double>::quiet_NaN());
        return values;
    }(); // l*l + component; unavailable values are not zero.
    double angular_partition_residual = std::numeric_limits<double>::quiet_NaN();
    std::string detail;
};

// Owns an immutable snapshot of the metric, AO definitions and MO coefficients.
// Prepare once outside group loops. Alpha and beta are separate orthogonal spin
// blocks; their traces and represented ranks are added after spatial projection.
class LocalAngularProjectionWorkspace {
public:
    LocalAngularProjectionWorkspace(const Wavefunction&,std::size_t atom_index,
        const std::array<double,9>& rotation_reference_to_input);
    LocalAngularProjection project(std::span<const std::size_t> orbital_indices) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
}
