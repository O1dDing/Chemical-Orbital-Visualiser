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

#include "cov/coordination_geometry.hpp"
#include "cov/model.hpp"
#include "cov/orbital_view.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace cov {

// Admission floor for a readable local d-shell explanation. Dimensions alone
// yield a candidate; only the S-metric projection provides angular evidence.
inline constexpr double kLocalDIrrepWeightFloor=0.02;

enum class LigandFieldGeometry {
    Unknown,
    Tetrahedral,
    Octahedral,
    General,
};

// A local coordination environment is deliberately separate from the full
// molecular point group.  Substituent orientations can lower the latter while
// leaving an unambiguous first-shell ligand field around a transition metal.
struct LigandFieldEnvironment {
    LigandFieldGeometry geometry = LigandFieldGeometry::Unknown;
    GeometryId geometry_id = GeometryId::Unknown;
    std::size_t metal_atom = 0;
    std::vector<std::size_t> ligand_atoms;
    double confidence = 0.0;
    double angular_rms = 0.0;
    double shape_measure = 0.0;
    double radial_cv = 0.0;
    std::array<double, 9> rotation_reference_to_input{
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
    bool equivalent_ligand_elements = false;
    bool ambiguous = false;

    [[nodiscard]] bool available() const noexcept {
        return geometry_id != GeometryId::Unknown && !ambiguous &&
               !ligand_atoms.empty();
    }

    [[nodiscard]] std::size_t coordination_number() const noexcept {
        return ligand_atoms.size();
    }

    [[nodiscard]] std::string geometry_machine_id() const;
    [[nodiscard]] std::string geometry_name() const;
    [[nodiscard]] std::string local_point_group() const;
};

// Uses density-derived Mayer connectivity to isolate the first coordination
// shell, removes collinear through-ligand contacts, and compares CN 2--10 with
// the shared ideal-geometry catalogue.  No molecule name, atom ordering, or
// canonical-MO number is consulted.
[[nodiscard]] LigandFieldEnvironment analyse_ligand_field_environment(
    const Wavefunction& wavefunction);

// Add local explanations to symmetry_view, with explicit target membership,
// geometry/axes and evidence. MolecularOrbital::symmetry and metadata.symmetry
// retain the full-orbital label; dimension/partner/spin candidates stay marked.
void apply_local_ligand_field_symmetry(
    const Wavefunction& wavefunction,
    std::vector<OrbitalMetadata>& metadata);

} // namespace cov
