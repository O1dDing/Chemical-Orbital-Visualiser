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
#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/orbital_view.hpp"
#include "cov/ui.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace cov::ui {
struct NboUIState;

struct OrbitalUIDiagramCache {
    const Wavefunction* wavefunction = nullptr;
    const MolecularOrbital* orbital_data = nullptr;
    std::size_t atom_count = 0;
    std::size_t orbital_count = 0;
    std::optional<MODiagramOptions> options;
    std::optional<MODiagramData> data;
    std::shared_ptr<const MODiagramViewSnapshot> snapshot;
};

struct OrbitalUIBrowserCache {
    const Wavefunction* wavefunction = nullptr;
    const MolecularOrbital* orbital_data = nullptr;
    std::size_t atom_count = 0;
    std::size_t orbital_count = 0;
    std::optional<DegeneracySettings> degeneracy;
    std::optional<OrbitalFilterSettings> filter;
    std::optional<FrontierOrbitals> frontier;
    std::vector<OrbitalMetadata> metadata;
    std::string routed_identity;
};

struct OrbitalUIState {
    std::optional<ActiveOrbitalView> active_view; // current scene identity, refreshed before drawing
    const NboSelectionView* inspection = nullptr; // main-owned; refreshed with active_view
    NboUIState* nbo_ui = nullptr; // owned by the application; same lifetime as this UI state
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    EnergyAxisMode energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    bool use_ro_common_energy=true;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    bool grouped_labels = true;
    // Chemical subsets are an explicit choice, independent of compact layout.
    MODiagramMode diagram_mode = MODiagramMode::ValenceCentral;
    bool hide_ligand_centred_intermediates = true;
    bool show_diagram_details = false;
    // Only an explicit open/reopen action requests focus. Background controls
    // never cover this floating window or repeatedly steal keyboard focus.
    bool focus_diagram_details = false;
    // Last drawn details rectangle (x, y, width, height), in ImGui coordinates.
    std::optional<std::array<float, 4>> diagram_details_bounds;
    int diagram_neighbourhood = 12;
    std::array<char, 96> search{};
    OrbitalUIBrowserCache browser_cache{};
    OrbitalUIDiagramCache diagram_cache{};
};

struct OrbitalUIActions {
    std::optional<std::size_t> select_orbital;
    bool export_diagram = false;
    bool export_analysis = false;
    std::shared_ptr<const MODiagramViewSnapshot> drawn_diagram;
};

void draw_orbital_browser(const Wavefunction& wavefunction,
                          std::size_t selected_index,
                          OrbitalUIState& state,
                          Language language,
                          float ui_scale,
                          OrbitalUIActions& actions);

void draw_energy_diagram(const Wavefunction& wavefunction,
                         std::size_t selected_index,
                         OrbitalUIState& state,
                         Language language,
                         float ui_scale,
                         OrbitalUIActions& actions);

} // namespace cov::ui
