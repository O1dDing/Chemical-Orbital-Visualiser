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
#include "cov/molecule_style.hpp"
#include "cov/molecular_overlay.hpp"

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace cov {

struct OrbitCamera {
    float yaw = 0.65f;
    float pitch = 0.35f;
    float distance = 2.2f;
    float fov_degrees = 42.0f; // Field of view across the shorter viewport dimension.
};

enum class OrbitalMaterial : std::uint8_t {
    Standard = 0,
    Glass,
};

enum class OrbitalSurfaceMode : std::uint8_t {
    Solid = 0,
    Wire,
    SolidWire,
};

class VolumeRenderer {
public:
    VolumeRenderer();
    ~VolumeRenderer();

    VolumeRenderer(const VolumeRenderer&) = delete;
    VolumeRenderer& operator=(const VolumeRenderer&) = delete;

    void resize_volume(int nx, int ny, int nz);
    // Call after replacing or mutating the wavefunction stored at an existing
    // address.  Interaction analysis is intentionally cached because the
    // layered graph contains an O(N^2) weak-contact pass and must not run once
    // per rendered frame.
    void invalidate_geometry_cache() noexcept;

    [[nodiscard]] unsigned int volume_texture() const noexcept { return texture_; }
    [[nodiscard]] int nx() const noexcept { return nx_; }
    [[nodiscard]] int ny() const noexcept { return ny_; }
    [[nodiscard]] int nz() const noexcept { return nz_; }

    void render_volume(int framebuffer_width,
                       int framebuffer_height,
                       float isovalue,
                       const OrbitCamera& camera,
                       float opacity = 1.0f,
                       OrbitalMaterial material = OrbitalMaterial::Standard,
                       OrbitalSurfaceMode surface_mode = OrbitalSurfaceMode::Solid,
                       int phase_palette = 0);

    void render_geometry(const Wavefunction& wavefunction,
                         const GridBox& box,
                         int framebuffer_width,
                         int framebuffer_height,
                         const OrbitCamera& camera,
                         const MoleculeRenderSettings& settings = {},
                         const MoleculeOverlay* overlay = nullptr,
                         const InteractionGraph* interactions = nullptr);

    [[nodiscard]] const std::vector<GeometryTarget>& geometry_targets() const noexcept { return geometry_targets_; }
    [[nodiscard]] std::optional<GeometryTarget> pick_geometry(float x, float y) const;

private:
    unsigned int texture_ = 0;
    unsigned int program_ = 0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    const Wavefunction* geometry_cache_wavefunction_ = nullptr;
    const InteractionGraph* geometry_graph_source_ = nullptr;
    std::vector<BondVisual> geometry_bonds_;
    std::map<std::pair<std::size_t,std::size_t>,std::size_t>
        geometry_bond_indices_;
    InteractionGraph geometry_interactions_;
    std::vector<GeometryTarget> geometry_targets_;
    int geometry_width_=1, geometry_height_=1;
};

} // namespace cov
