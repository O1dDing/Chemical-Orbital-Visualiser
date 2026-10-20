#pragma once
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cov {
// Rendering input only: all indices refer to the immutable canonical molecule.
// Classification and numerical evidence are supplied by the analysis layer.
enum class OverlayBondStyle { Covalent, Coordination, Delocalised, Multicentre, Unresolved };
enum class AtomScalarMode { Element, NaturalCharge, SpinPopulation };
enum class NboBondDisplayMode { DefaultSkeleton, LewisStructure };
struct MoleculeOverlayBond {
    std::size_t atom_a=0, atom_b=0, evidence_index=0;
    int multiplicity=1; // display decision; never rounded from a continuous index
    std::optional<unsigned> lewis_bond_count; // source evidence, independent of display
    std::vector<std::size_t> evidence_indices; // all associated spin/source records
    bool density_equivalence_checked=false, lewis_equivalence_conflict=false;
    OverlayBondStyle style=OverlayBondStyle::Covalent;
    std::optional<double> continuous_index;
    bool selected=false;
};
struct MoleculeOverlayGroup {
    std::size_t evidence_index=0;
    std::vector<std::size_t> atoms;
    bool selected=false;
    // Several spin/channel records may support the same spatial hyperedge.
    std::vector<std::size_t> evidence_indices;
};
struct MoleculeOverlayRelation {
    std::size_t evidence_index=0;
    std::vector<std::size_t> donor_atoms, acceptor_atoms;
    bool selected=false;
};
struct MoleculeOverlay {
    std::vector<MoleculeOverlayBond> bonds;
    std::vector<MoleculeOverlayGroup> multicentre;
    std::vector<MoleculeOverlayRelation> relations;
    std::vector<std::optional<double>> atom_values;
    std::vector<std::size_t> selected_atoms;
    AtomScalarMode colour_mode=AtomScalarMode::Element;
    double scalar_range=1;
    bool show_bond_indices=false;
    NboBondDisplayMode bond_mode=NboBondDisplayMode::DefaultSkeleton;
};
enum class GeometryTargetKind { Atom, Bond, Relation, Multicentre };
struct GeometryTarget {
    GeometryTargetKind kind=GeometryTargetKind::Atom;
    std::size_t index=0;
    float x=0, y=0, end_x=0, end_y=0, radius=0, depth=0;
    // Coordinates are normalized to the scene viewport, top-left origin;
    // radius is in framebuffer pixels, allowing aspect-correct picking.
    bool segment=false;
};
}
