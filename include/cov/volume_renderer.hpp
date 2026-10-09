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
