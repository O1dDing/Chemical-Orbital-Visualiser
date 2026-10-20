#pragma once

#include "cov/model.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace cov {

// Pairwise interaction semantics are deliberately separate from BondVisual.
// Only the first four kinds are strong molecular connectivity. The remaining
// kinds are contacts and must never be rendered as ordinary covalent bonds or
// used to merge molecular fragments.
enum class InteractionKind : std::uint8_t {
    CovalentConnectivity = 0,
    CoordinationContact,
    MulticentreSupport,
    PolyhedralCageSupport,
    HydrogenBond,
    NoncovalentContact,
    IonicContact,
    AmbiguousContact,
};

enum class InteractionStrength : std::uint8_t {
    StrongConnectivity = 0,
    WeakContact,
};

struct InteractionEdge {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    InteractionKind kind = InteractionKind::AmbiguousContact;
    InteractionStrength strength = InteractionStrength::WeakContact;

    double distance_bohr = 0.0;
    // R_AB/(r_cov,A+r_cov,B). This dimensionless value is useful for
    // diagnostics only; it is not by itself a bond-order claim.
    double covalent_radius_ratio = 0.0;
    double mayer_order = 0.0;
    double confidence = 0.0;
    DataProvenance electronic_provenance = DataProvenance::Unavailable;
    // Distinct validated NBO observable; never substituted for Mayer order.
    std::optional<double> wiberg_index;
    std::string wiberg_source_path;
    std::string connectivity_method, connectivity_source_path;
};

// A multicentre assignment is a hyperedge. Pairwise MulticentreSupport edges
// merely keep its spatially adjacent atoms in one strong fragment; they do not
// reduce the electronic assignment to a collection of two-centre bonds.
struct MulticentreSupportGroup {
    MulticentreKind kind = MulticentreKind::Unclassified;
    std::vector<std::uint32_t> atoms;
    std::vector<std::uint32_t> orbitals;
    double electron_count = 0.0;
    std::string source_subspace_id;
    double source_subspace_electron_count = 0.0;
    double source_subspace_fraction = 1.0;
    double confidence = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

enum class PolyhedralCageKind : std::uint8_t {
    Icosahedral = 0,
};

struct PolyhedralCage {
    PolyhedralCageKind kind = PolyhedralCageKind::Icosahedral;
    std::vector<std::uint32_t> atoms;
    std::vector<std::array<std::uint32_t, 2>> support_edges;
    double confidence = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

struct MolecularFragment {
    std::uint32_t index = 0;
    std::vector<std::uint32_t> atoms;
    std::vector<std::uint32_t> strong_edge_indices;
    // Populated only when compatible atom-resolved charge evidence is passed
    // to build_interaction_graph().
    double evidenced_charge = 0.0;
    DataProvenance charge_provenance = DataProvenance::Unavailable;
};

struct FragmentAnalysis {
    static constexpr std::uint32_t no_fragment =
        std::numeric_limits<std::uint32_t>::max();

    std::vector<MolecularFragment> fragments;
    std::vector<std::uint32_t> atom_to_fragment;
    std::vector<std::uint32_t> weak_interfragment_edge_indices;
};

// Charges are optional analysis evidence, not inferred from element names or
// from the total molecular charge. This keeps isolated weakly coordinating
// anions and contact ion pairs representation-independent.
struct InteractionBondEvidence {
    std::size_t atom_a=0, atom_b=0;
    // Only associated, independently validated total WBI records belong here.
    // A missing pair is absent, never an implicit zero or a spin sum.
    double wiberg_index=0;
    std::string source_path;
};
struct InteractionEvidence {
    std::vector<double> atomic_partial_charges;
    DataProvenance atomic_charge_provenance = DataProvenance::Unavailable;
    std::vector<InteractionBondEvidence> bonds;

    [[nodiscard]] bool has_atomic_charges(std::size_t atom_count) const noexcept;
};

struct InteractionGraphOptions {
    // Low-order metal--donor connections may be chemically structural even
    // below the renderer's ordinary Mayer threshold.
    double coordination_mayer_floor = 0.02;
    double coordination_distance_factor = 1.65;
    double coordination_shield_margin_angstrom = 0.25;

    double multicentre_distance_factor = 1.65;

    // Weak contacts are first collected geometrically and then restricted to
    // the closest surface shell for each pair of strong fragments.
    double weak_contact_distance_factor = 2.75;
    double weak_contact_max_angstrom = 4.0;
    double weak_contact_shell_angstrom = 0.35;
    double weak_record_abs_mayer_floor = 0.01;

    double hydrogen_bond_max_h_acceptor_angstrom = 2.65;
    double hydrogen_bond_min_angle_degrees = 120.0;

    double ionic_atomic_charge_floor = 0.25;
    double ionic_fragment_charge_floor = 0.50;
    // Strong electronic evidence outranks a polar atom-pair charge pattern.
    double ionic_strong_bond_override_mayer = 0.20;
};

struct InteractionGraph {
    std::size_t atom_count = 0;
    // Distinct semantic layers may intentionally share one atom pair. For
    // example, a cage edge can simultaneously support a local multicentre
    // hyperedge and the global polyhedral skeleton. Consumers must filter by
    // InteractionKind rather than assuming atom-pair uniqueness.
    std::vector<InteractionEdge> edges;
    std::vector<MulticentreSupportGroup> multicentre_groups;
    std::vector<PolyhedralCage> polyhedral_cages;
    FragmentAnalysis fragment_analysis;
};

// Geometry-derived global cage topology is separate from a local CN centre.
// The current detector recognises an icosahedral 12-vertex nearest-neighbour
// framework after terminal substituents are removed by a graph k-core. It does
// not label the support edges as ordinary two-centre bonds.
[[nodiscard]] std::vector<PolyhedralCage> analyse_polyhedral_cages(
    const Wavefunction& wavefunction);

[[nodiscard]] bool interaction_merges_fragments(InteractionKind kind) noexcept;
[[nodiscard]] bool interaction_is_ordinary_bond(InteractionKind kind) noexcept;
[[nodiscard]] const char* interaction_kind_name(InteractionKind kind) noexcept;

// Builds a molecule-name-independent layered graph from the existing Mayer,
// structural-distance and multicentre evidence. Strong edges define fragments;
// hydrogen-bond, ionic, noncovalent and ambiguous contacts never do.
[[nodiscard]] InteractionGraph build_interaction_graph(
    const Wavefunction& wavefunction,
    const InteractionEvidence& evidence = {},
    const InteractionGraphOptions& options = {});

} // namespace cov
