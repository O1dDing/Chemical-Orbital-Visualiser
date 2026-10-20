#pragma once

#include "cov/interaction_graph.hpp"
#include "cov/model.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace cov {

enum class MoleculeStyle {
    MediumBallAndStick = 0, // UI: Enhanced Ball-and-Stick
    StickDelocalisation,
};

struct MoleculeRenderSettings {
    MoleculeStyle style = MoleculeStyle::MediumBallAndStick;
    // 1.0 is the human-facing reference size. The renderer's base sphere and
    // cylinder radii are calibrated so 1.0 reproduces the proportions that
    // were manually accepted during pre.2 validation.
    float atom_scale = 1.00f;
    float bond_scale = 1.00f;
    float molecule_opacity = 1.00f;
    float orbital_opacity = 0.92f;
    bool show_hydrogens = true;
    bool show_numbers=false, number_atoms=true, number_fragments=true;
    bool number_ignore_h=true;
    // Strong coordination connectivity uses the same solid structural line.
    // Its chemical semantics remain independent of the visual style.
    bool show_coordination_contacts = true;
    // Keep the default ball-and-stick skeleton clean. The independent
    // multicentre layer remains available for explicit inspection.
    bool show_multicentre_support = false;
    // Polyhedral cage support is a global skeletal/topological inference, not
    // a 3c electron assignment, and therefore has its own visual layer.
    bool show_polyhedral_cage_support = true;
    // Hydrogen bonds, generic non-covalent contacts and ion-pair contacts are
    // an optional context layer.  Ambiguous contacts remain hidden even when
    // this layer is enabled.
    bool show_weak_interactions = false;
};

enum class InteractionVisualStyle : std::uint8_t {
    Hidden = 0,
    OrdinaryBond,
    CoordinationDash,
    MulticentreDash,
    PolyhedralCageDash,
    HydrogenBondDots,
    NoncovalentDots,
    IonicDash,
};

// A single, testable semantic-to-visual contract used by the renderer.  In
// particular, weak or ambiguous contacts can never fall through to the
// ordinary-bond style.
[[nodiscard]] constexpr InteractionVisualStyle interaction_visual_style(
    const InteractionKind kind,
    const MoleculeRenderSettings& settings) noexcept {
    switch (kind) {
        case InteractionKind::CovalentConnectivity:
            return InteractionVisualStyle::OrdinaryBond;
        case InteractionKind::CoordinationContact:
            return settings.show_coordination_contacts
                       ? InteractionVisualStyle::OrdinaryBond
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::MulticentreSupport:
            return settings.show_multicentre_support
                       ? InteractionVisualStyle::MulticentreDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::PolyhedralCageSupport:
            return settings.show_polyhedral_cage_support
                       ? InteractionVisualStyle::PolyhedralCageDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::HydrogenBond:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::HydrogenBondDots
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::NoncovalentContact:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::NoncovalentDots
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::IonicContact:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::IonicDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::AmbiguousContact:
            return InteractionVisualStyle::Hidden;
    }
    return InteractionVisualStyle::Hidden;
}

struct BondVisual {
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double distance_bohr = 0.0;
    bool delocalised = false;
    double bond_order = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
    std::string connectivity_method, connectivity_source_path;
};

// Structural radii through Cm follow Cordero et al. (Dalton Trans. 2008,
// DOI:10.1039/B801115J); the remaining super-heavy values use Pyykko--Atsumi
// single-bond radii (Chem. Eur. J. 2009, DOI:10.1002/chem.200800987). Keeping
// this complete Z=1..118 table inline gives analysis and rendering exactly the
// same first-neighbour envelope without a cross-library link dependency.
[[nodiscard]] inline double covalent_radius_angstrom(
    const int atomic_number) noexcept {
    static constexpr std::array<double,119> radii{
        0.85,
        0.31,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,0.58,
        1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,
        2.03,1.76,1.70,1.60,1.53,1.39,1.39,1.32,1.26,1.24,
        1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,
        2.20,1.95,1.90,1.75,1.64,1.54,1.47,1.46,1.42,1.39,
        1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,
        2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
        1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,
        1.51,1.44,1.41,1.36,1.36,1.32,1.45,1.46,1.48,1.40,
        1.50,1.50,
        2.60,2.21,2.15,2.06,2.00,1.96,1.90,1.87,1.80,1.69,
        1.68,1.68,1.65,1.67,1.73,1.76,1.61,1.57,1.49,1.43,
        1.41,1.34,1.29,1.28,1.21,1.22,1.36,1.43,1.62,1.75,
        1.65,1.57,
    };
    if (atomic_number>0 &&
        static_cast<std::size_t>(atomic_number)<radii.size()) {
        return radii[static_cast<std::size_t>(atomic_number)];
    }
    return radii[0];
}

// Wavefunction-derived Mayer evidence is used when available, inside a
// conservative structural-neighbour envelope. This keeps non-local electronic
// coupling from being drawn as a bond. Covalent-radius geometry is the fallback
// only when electronic evidence is unavailable. Pairwise Mayer evidence is not
// an electron-counting or multicentre claim.
[[nodiscard]] std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction,
    const std::vector<InteractionBondEvidence>& nbo_bonds = {});

} // namespace cov
