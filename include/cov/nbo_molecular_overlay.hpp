#pragma once
#include "cov/molecular_overlay.hpp"
#include "cov/nbo_integration.hpp"
#include "cov/interaction_graph.hpp"
#include "cov/chemistry_route.hpp"
namespace cov {
struct NboElectronicSymmetryScope;
MoleculeOverlay make_nbo_molecule_overlay(const NboIntegration&,const InteractionGraph&,
    std::size_t atom_count,const std::vector<std::size_t>& selected_atoms,
    std::optional<std::size_t> selected_structure,AtomScalarMode,bool show_indices,bool show_e2,
    const RoutedAnalysis* routed=nullptr,const Wavefunction* canonical=nullptr,
    NboBondDisplayMode bond_mode=NboBondDisplayMode::DefaultSkeleton,
    const NboElectronicSymmetryScope* electronic_scope=nullptr);
std::string serialize_molecule_overlay_scalars_json(const MoleculeOverlay&,
                                                     const RoutedAnalysis* routed=nullptr);
}
