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
#include <optional>
#include <string_view>
#include <vector>

namespace cov {

// Stable semantic identifiers.  Display names and point groups live in the
// catalogue, so persistence and tests never depend on translated UI text.
enum class GeometryId {
    Unknown = 0,
    Linear2,
    Angular2,
    TrigonalPlanar3,
    TrigonalPyramidal3,
    TShaped3,
    Tetrahedral4,
    SquarePlanar4,
    Seesaw4,
    TrigonalBipyramidal5,
    SquarePyramidal5,
    Octahedral6,
    TrigonalPrismatic6,
    PentagonalBipyramidal7,
    CappedOctahedral7,
    CappedTrigonalPrismatic7,
    SquareAntiprismatic8,
    TriangularDodecahedral8,
    BicappedTrigonalPrismatic8,
    CappedSquareAntiprismatic9,
    TricappedTrigonalPrismatic9,
    BicappedSquareAntiprismatic10,
    Sphenocorona10,
    Tetradecahedral10,
    // Appended to preserve the numeric values of the original public IDs.
    PentagonalPrismatic10,
    PentagonalAntiprismatic10,
    TrigonalPyramidal4,
};

struct CoordinationDirection {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

enum class ReferenceGeometryStatus {
    Available,
    PendingAuthoritativeCoordinates,
};

struct CoordinationGeometryDescriptor {
    GeometryId id = GeometryId::Unknown;
    std::string_view machine_id;
    std::string_view name;
    std::string_view point_group;
    std::size_t coordination_number = 0;
    ReferenceGeometryStatus reference_status =
        ReferenceGeometryStatus::PendingAuthoritativeCoordinates;
    std::string_view reference_note;
    std::vector<CoordinationDirection> reference_directions;

    [[nodiscard]] bool matchable() const noexcept {
        return reference_status == ReferenceGeometryStatus::Available &&
               reference_directions.size() == coordination_number;
    }
};

[[nodiscard]] const std::vector<CoordinationGeometryDescriptor>&
coordination_geometry_catalog();

[[nodiscard]] const CoordinationGeometryDescriptor*
coordination_geometry_descriptor(GeometryId id) noexcept;

[[nodiscard]] const CoordinationGeometryDescriptor*
coordination_geometry_descriptor(std::string_view machine_id) noexcept;

[[nodiscard]] std::vector<const CoordinationGeometryDescriptor*>
coordination_geometries_for_cn(std::size_t coordination_number);

struct GeometryCandidate {
    GeometryId id = GeometryId::Unknown;
    // RMS angular displacement after the optimal assignment/alignment, radians.
    double angular_rms = 0.0;
    // COV's direction-only continuous score: 100 * angular_rms^2.  This is
    // deliberately not labelled as the full-coordinate SHAPE CShM value.
    double shape_measure = 0.0;
    // Coefficient of variation of the input metal-ligand radii.
    double radial_cv = 0.0;
    // Row-major proper rotation satisfying v_input = R * v_reference.
    std::array<double, 9> rotation_reference_to_input{
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
    // assignment[input_index] is the matched reference-vertex index.
    std::vector<std::size_t> assignment;
};

struct GeometryMatchOptions {
    double max_angular_rms = 0.40;
    double ambiguity_margin = 0.035;
    std::size_t max_iterations = 12;
};

struct GeometryMatch {
    std::size_t coordination_number = 0;
    std::optional<GeometryCandidate> best;
    std::optional<GeometryCandidate> runner_up;
    std::vector<GeometryId> unavailable_candidates;
    bool accepted = false;
    bool ambiguous = false;
};

// Input vector magnitudes are used only for radial_cv; angular matching uses
// their unit directions.  Zero-length or non-finite inputs produce no match.
[[nodiscard]] GeometryMatch match_coordination_geometry(
    const std::vector<CoordinationDirection>& input_vectors,
    const GeometryMatchOptions& options = {});

struct CoordinationContact {
    std::size_t atom_index = 0;
    double distance = 0.0;
    double mayer_bond_order = 0.0;
    CoordinationDirection direction;
};

struct CoordinationShell {
    std::size_t centre_atom = 0;
    std::vector<CoordinationContact> contacts;
    std::size_t radial_candidate_count = 0;
    bool used_low_mayer_retry = false;
};

struct CoordinationShellOptions {
    double minimum_mayer_bond_order = 0.02;
    double absolute_electronic_floor = 0.035;
    double relative_electronic_floor = 0.08;
    double radial_shadow_cosine = 0.94;
    GeometryMatchOptions match;
};

// Extracts the nearest contact on each ligand ray.  Unlike the legacy Td/Oh
// path, this function has no six-contact truncation; CN 7--10 remains intact.
[[nodiscard]] CoordinationShell extract_coordination_shell(
    const Wavefunction& wavefunction,
    std::size_t centre_atom,
    const CoordinationShellOptions& options = {});

[[nodiscard]] GeometryMatch analyse_coordination_shell(
    const CoordinationShell& shell,
    const GeometryMatchOptions& options = {});

} // namespace cov
