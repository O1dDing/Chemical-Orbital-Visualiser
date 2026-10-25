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

#include "cov/symmetry.hpp"

#include <optional>
#include <string_view>

namespace cov {

struct PointGroupIrrepRow {
    std::string label;
    std::size_t dimension = 0;
    // <chi,chi> over all operations. Real conjugate-pair chemical E rows in
    // cyclic groups and T have norm 2, not 1; their dimension is the real 2D.
    double character_norm = 1.0;
    std::vector<double> characters;
};

struct PointGroupIrrepOptions {
    std::optional<std::array<double, 3>> principal_axis;
    std::optional<std::array<double, 3>> reference_axis;
    double tolerance = 2.0e-5;
};

struct PointGroupIrrepTable {
    std::string point_group;
    bool valid = false;
    std::string reason;
    std::string axis_detail;
    std::array<double, 3> principal_axis{0, 0, 1};
    std::array<double, 3> reference_axis{1, 0, 0};
    std::size_t identity_operation = 0;
    // Every row follows the caller's operation ordering. Matrix data, rather
    // than potentially stale kind/order/axis metadata, determines characters.
    std::vector<PointGroupIrrepRow> rows;
};

struct PointGroupCharacterDecomposition {
    bool valid = false;
    std::string reason;
    // Same order as table.rows; number of complete real chemical irrep copies.
    std::vector<std::size_t> multiplicities;
    std::size_t dimension = 0;
    double maximum_reconstruction_error = 0;
};

// Zero means unsupported, continuous/linear, atomic, or invalid group name.
// Parameterized axial families have no catalogue-dependent n limit.
[[nodiscard]] std::size_t finite_point_group_order(std::string_view group);
[[nodiscard]] PointGroupIrrepTable finite_point_group_irreps(
    const MolecularSymmetry& symmetry,
    const PointGroupIrrepOptions& options = {});
[[nodiscard]] PointGroupCharacterDecomposition decompose_point_group_characters(
    const PointGroupIrrepTable& table,
    const std::vector<double>& characters,
    double tolerance = 2.0e-4);

} // namespace cov
