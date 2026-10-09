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
