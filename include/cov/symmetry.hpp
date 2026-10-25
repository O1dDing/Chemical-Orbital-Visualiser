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

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cov {

enum class SymmetryOperationKind : std::uint8_t {
    Identity = 0,
    ProperRotation,
    ImproperRotation,
    Reflection,
    Inversion,
};

struct SymmetryOperation {
    SymmetryOperationKind kind = SymmetryOperationKind::Identity;
    int order = 1;
    int power = 0;
    std::array<double, 3> axis_or_normal{0.0, 0.0, 1.0};
    std::array<double, 9> matrix{
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
    std::vector<std::size_t> atom_permutation;
    double max_mapping_error_bohr = 0.0;
};

struct MolecularSymmetry {
    std::string point_group;
    std::array<double, 3> centre_bohr{0.0, 0.0, 0.0};
    double molecular_radius_bohr = 0.0;
    double tolerance_bohr = 0.0;
    double max_mapping_error_bohr = 0.0;
    bool linear = false;
    bool has_inversion = false;
    std::vector<SymmetryOperation> operations;

    [[nodiscard]] bool available() const noexcept {
        return !point_group.empty();
    }
};

struct SymmetryOptions {
    // Symmetry is validated by atom permutation, never by chemistry-specific
    // pattern matching. The effective tolerance is the larger of the absolute
    // and radius-scaled tolerances.
    double absolute_tolerance_bohr = 2.0e-3;
    double relative_tolerance = 2.0e-4;
    int maximum_rotation_order = 12;
    // Retained for source compatibility. Verified molecular isomorphisms are
    // enumerated completely and are not truncated by this former search limit.
    std::size_t maximum_candidate_axes = 768;
};

// Independent geometry-symmetry engine. Candidate axes/planes are generated
// from the molecular geometry, every operation is validated by a same-element
// atom permutation, and the point-group family is classified from the accepted
// operation set. No external symmetry library or molecule-specific templates
// are used.
[[nodiscard]] MolecularSymmetry analyse_molecular_symmetry(
    const Wavefunction& wavefunction,
    const SymmetryOptions& options = {});

// Populate the Wavefunction point-group fields only when producer-reported
// symmetry is absent. Producer metadata always outranks geometry-derived data.
void derive_point_group_from_geometry(
    Wavefunction& wavefunction,
    const SymmetryOptions& options = {});

} // namespace cov
