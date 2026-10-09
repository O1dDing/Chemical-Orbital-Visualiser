#include "cov/molecule_style.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

namespace {

void add_all_pair_mayer_records(cov::Wavefunction& wf, const double order) {
    wf.bond_order_provenance = cov::DataProvenance::Derived;
    for (std::size_t i = 0; i < wf.atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < wf.atoms.size(); ++j) {
            wf.bond_orders.push_back({
                static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j),
                order, cov::DataProvenance::Derived});
        }
    }
}

bool check_carbon_graph(const cov::Wavefunction& wf,
                        const std::vector<cov::BondVisual>& bonds,
                        const std::size_t expected_edges,
                        const std::vector<int>& expected_degrees) {
    if (bonds.size() != expected_edges || expected_degrees.size() != wf.atoms.size()) {
        return false;
    }
    std::vector<int> degrees(wf.atoms.size(), 0);
    for (const auto& bond : bonds) {
        if (bond.atom_a >= degrees.size() || bond.atom_b >= degrees.size()) return false;
        ++degrees[bond.atom_a];
        ++degrees[bond.atom_b];
    }
    return degrees == expected_degrees;
}

bool check_real_topology(const std::filesystem::path& path,
                         const std::size_t expected_bonds,
                         const std::size_t expected_carbon_bonds) {
    const auto wf = cov::parse_wavefunction(path);
    const auto bonds = cov::analyse_bonds(wf);
    std::size_t carbon_bonds = 0;
    for (const auto& bond : bonds) {
        const auto& a = wf.atoms[bond.atom_a];
        const auto& b = wf.atoms[bond.atom_b];
        if (a.atomic_number == 6 && b.atomic_number == 6) ++carbon_bonds;
        const double radii_bohr =
            (cov::covalent_radius_angstrom(a.atomic_number) +
             cov::covalent_radius_angstrom(b.atomic_number)) * cov::kAngstromToBohr;
        if (bond.distance_bohr > 1.50 * radii_bohr + 1.0e-10) return false;
    }
    return bonds.size() == expected_bonds && carbon_bonds == expected_carbon_bonds;
}

bool has_central_degree(const cov::Wavefunction& wf,
                        const std::size_t expected_degree) {
    const auto bonds = cov::analyse_bonds(wf);
    std::size_t degree = 0;
    for (const auto& bond : bonds) {
        if (bond.atom_a == 0u || bond.atom_b == 0u) ++degree;
    }
    return bonds.size() == expected_degree && degree == expected_degree;
}

} // namespace

int main(const int argc, char** argv) {
    const cov::MoleculeRenderSettings defaults;
    if (std::abs(defaults.atom_scale - 1.0f) > 1.0e-6f ||
        std::abs(defaults.bond_scale - 1.0f) > 1.0e-6f ||
        defaults.molecule_opacity < 0.95f ||
        !defaults.show_coordination_contacts ||
        defaults.show_multicentre_support ||
        !defaults.show_polyhedral_cage_support ||
        defaults.show_weak_interactions) {
        std::cerr << "1.0 molecular-size reference defaults regressed\n";
        return 1;
    }

    if (cov::interaction_visual_style(cov::InteractionKind::CovalentConnectivity,
                                      defaults) !=
            cov::InteractionVisualStyle::OrdinaryBond ||
        cov::interaction_visual_style(cov::InteractionKind::CoordinationContact,
                                      defaults) !=
            cov::InteractionVisualStyle::OrdinaryBond ||
        cov::interaction_visual_style(cov::InteractionKind::MulticentreSupport,
                                       defaults) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::PolyhedralCageSupport,
                                       defaults) !=
            cov::InteractionVisualStyle::PolyhedralCageDash ||
        cov::interaction_visual_style(cov::InteractionKind::HydrogenBond,
                                      defaults) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::NoncovalentContact,
                                      defaults) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::IonicContact,
                                      defaults) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::AmbiguousContact,
                                      defaults) !=
            cov::InteractionVisualStyle::Hidden) {
        std::cerr << "default interaction visual mapping regressed\n";
        return 19;
    }

    cov::MoleculeRenderSettings optional_layers = defaults;
    optional_layers.show_coordination_contacts = false;
    optional_layers.show_multicentre_support = true;
    optional_layers.show_polyhedral_cage_support = false;
    optional_layers.show_weak_interactions = true;
    if (cov::interaction_visual_style(cov::InteractionKind::CoordinationContact,
                                      optional_layers) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::MulticentreSupport,
                                       optional_layers) !=
            cov::InteractionVisualStyle::MulticentreDash ||
        cov::interaction_visual_style(cov::InteractionKind::PolyhedralCageSupport,
                                       optional_layers) !=
            cov::InteractionVisualStyle::Hidden ||
        cov::interaction_visual_style(cov::InteractionKind::HydrogenBond,
                                      optional_layers) !=
            cov::InteractionVisualStyle::HydrogenBondDots ||
        cov::interaction_visual_style(cov::InteractionKind::NoncovalentContact,
                                      optional_layers) !=
            cov::InteractionVisualStyle::NoncovalentDots ||
        cov::interaction_visual_style(cov::InteractionKind::IonicContact,
                                      optional_layers) !=
            cov::InteractionVisualStyle::IonicDash ||
        cov::interaction_visual_style(cov::InteractionKind::AmbiguousContact,
                                      optional_layers) !=
            cov::InteractionVisualStyle::Hidden) {
        std::cerr << "optional interaction visual mapping regressed\n";
        return 20;
    }

    cov::Wavefunction wf;
    constexpr double pi = 3.14159265358979323846;
    constexpr double side_angstrom = 1.40;
    const double radius_angstrom = side_angstrom / (2.0 * std::sin(pi / 5.0));

    for (int i = 0; i < 5; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 5.0;
        cov::Atom atom;
        atom.symbol = "C";
        atom.atomic_number = 6;
        atom.x = radius_angstrom * std::cos(angle) * cov::kAngstromToBohr;
        atom.y = radius_angstrom * std::sin(angle) * cov::kAngstromToBohr;
        atom.z = 0.0;
        wf.atoms.push_back(atom);
    }

    const auto bonds = cov::analyse_bonds(wf);
    if (bonds.size() != 5) {
        std::cerr << "expected five ring bonds, got " << bonds.size() << '\n';
        return 2;
    }
    for (const auto& bond : bonds) {
        if (!bond.delocalised) {
            std::cerr << "compact equal five-member ring was not marked delocalised\n";
            return 3;
        }
    }

    // The exact five-membered-ring failure mode requires electronic evidence:
    // all-pair delocalised Mayer coupling places each short diagonal just
    // inside the broad 1.50 radius envelope. The two-hop lune filter must keep
    // only the five structural perimeter edges.
    cov::Wavefunction electronic_five_ring = wf;
    add_all_pair_mayer_records(electronic_five_ring, 0.06);
    const auto electronic_five_bonds = cov::analyse_bonds(electronic_five_ring);
    if (!check_carbon_graph(electronic_five_ring, electronic_five_bonds, 5u,
                            std::vector<int>(5u, 2))) {
        std::cerr << "five-membered electronic coupling created ring chords\n";
        return 23;
    }

    // A weak but electronically coherent, distorted three-membered ring must
    // not lose its longest edge merely because that edge has a two-hop route.
    // These sides are 1.2/1.2/1.8 A and all three Mayer indices are 0.10.
    cov::Wavefunction distorted_c3;
    const double distorted_height=std::sqrt(1.20*1.20-0.90*0.90);
    distorted_c3.atoms = {
        {"C", 6, -0.90*cov::kAngstromToBohr, 0.0, 0.0},
        {"C", 6,  0.90*cov::kAngstromToBohr, 0.0, 0.0},
        {"C", 6, 0.0, distorted_height*cov::kAngstromToBohr, 0.0},
    };
    add_all_pair_mayer_records(distorted_c3,0.10);
    const auto distorted_c3_bonds=cov::analyse_bonds(distorted_c3);
    if (!check_carbon_graph(distorted_c3,distorted_c3_bonds,3u,{2,2,2})) {
        std::cerr << "weak distorted three-membered ring lost an edge\n";
        return 30;
    }

    // A linear XY2 candidate can carry the same weak all-pair electronic
    // coupling, but its terminal--terminal pair is a non-local two-hop chord.
    // Only the two X--Y legs may survive.
    cov::Wavefunction linear_xy2;
    linear_xy2.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"Si", 14, -1.60*cov::kAngstromToBohr, 0.0, 0.0},
        {"Si", 14,  1.60*cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_all_pair_mayer_records(linear_xy2,0.10);
    const auto linear_xy2_bonds=cov::analyse_bonds(linear_xy2);
    if (linear_xy2_bonds.size()!=2u ||
        std::any_of(linear_xy2_bonds.begin(),linear_xy2_bonds.end(),
                    [](const cov::BondVisual& bond) {
                        return bond.atom_a!=0u && bond.atom_b!=0u;
                    })) {
        std::cerr << "linear XY2 retained a terminal--terminal chord\n";
        return 31;
    }

    // A square four-membered ring has two especially short diagonals. They
    // are non-local coupling, not structural bonds, even when every pair has a
    // small positive Mayer value.
    cov::Wavefunction electronic_square;
    electronic_square.atoms = {
        {"C", 6, -0.72 * cov::kAngstromToBohr,
         -0.72 * cov::kAngstromToBohr, 0.0},
        {"C", 6,  0.72 * cov::kAngstromToBohr,
         -0.72 * cov::kAngstromToBohr, 0.0},
        {"C", 6,  0.72 * cov::kAngstromToBohr,
          0.72 * cov::kAngstromToBohr, 0.0},
        {"C", 6, -0.72 * cov::kAngstromToBohr,
          0.72 * cov::kAngstromToBohr, 0.0},
    };
    add_all_pair_mayer_records(electronic_square, 0.06);
    const auto square_bonds = cov::analyse_bonds(electronic_square);
    if (!check_carbon_graph(electronic_square, square_bonds, 4u,
                            std::vector<int>(4u, 2))) {
        std::cerr << "four-membered electronic coupling created diagonals\n";
        return 24;
    }

    cov::Wavefunction saturated = wf;
    // Stretch the ring into a single-bond-like regime. The conservative
    // detector must stop claiming delocalisation.
    for (std::size_t i = 0; i < saturated.atoms.size(); ++i) {
        saturated.atoms[i].x *= 1.10;
        saturated.atoms[i].y *= 1.10;
    }
    const auto saturated_bonds = cov::analyse_bonds(saturated);
    bool any_delocalised = false;
    for (const auto& bond : saturated_bonds) any_delocalised = any_delocalised || bond.delocalised;
    if (any_delocalised) {
        std::cerr << "long single-bond-like ring was incorrectly marked delocalised\n";
        return 4;
    }

    // Once pairwise electronic evidence is available it must outrank the old
    // 1.25*r_cov geometry rule. This synthetic C--C separation lies outside the
    // geometry fallback cutoff but inside the deliberately generous electronic
    // sanity distance; a positive Mayer index should therefore retain the bond.
    cov::Wavefunction electronic;
    electronic.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 2.25 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    electronic.bond_order_provenance = cov::DataProvenance::Derived;
    electronic.bond_orders.push_back({0, 1, 0.18, cov::DataProvenance::Derived});
    const auto electronic_bonds = cov::analyse_bonds(electronic);
    if (electronic_bonds.size() != 1u ||
        electronic_bonds[0].provenance != cov::DataProvenance::Derived ||
        std::abs(electronic_bonds[0].bond_order - 0.18) > 1.0e-12) {
        std::cerr << "Mayer connectivity did not outrank geometry fallback\n";
        return 5;
    }

    // Signed Mayer decompositions may assign a negative contribution to a
    // short, directly bonded pair. Its magnitude remains structural evidence,
    // while the same contribution at a remote distance must remain invisible.
    cov::Wavefunction signed_direct;
    signed_direct.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"F", 9, 1.34 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    signed_direct.bond_order_provenance = cov::DataProvenance::Derived;
    signed_direct.bond_orders.push_back(
        {0, 1, -0.15, cov::DataProvenance::Derived});
    const auto signed_direct_bonds=cov::analyse_bonds(signed_direct);
    if (signed_direct_bonds.size()!=1u ||
        std::abs(signed_direct_bonds.front().bond_order-0.15)>1.0e-12) {
        std::cerr << "short signed Mayer contact was not retained\n";
        return 27;
    }

    cov::Wavefunction signed_crowded=signed_direct;
    signed_crowded.atoms[1].x=1.70*cov::kAngstromToBohr;
    if (!cov::analyse_bonds(signed_crowded).empty()) {
        std::cerr << "negative Mayer coupling escaped the strict geometry gate\n";
        return 29;
    }

    cov::Wavefunction signed_remote=signed_direct;
    signed_remote.atoms[1].x=3.50*cov::kAngstromToBohr;
    if (!cov::analyse_bonds(signed_remote).empty()) {
        std::cerr << "remote signed Mayer coupling became a structural bond\n";
        return 28;
    }

    // Capability is per pair. A global provenance flag with no actual pair
    // value describes missing evidence, rather than a measured zero.
    cov::Wavefunction electronic_none;
    electronic_none.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 1.40 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    electronic_none.bond_order_provenance = cov::DataProvenance::Derived;
    const auto missing_pair_bonds=cov::analyse_bonds(electronic_none);
    if (missing_pair_bonds.size()!=1u ||
        missing_pair_bonds.front().provenance!=cov::DataProvenance::Unavailable) {
        std::cerr << "missing pair was converted into an electronic zero\n";
        return 32;
    }

    // An actual zero must veto geometry even at a conventional bond length.
    // Keep this negative assertion distinct from the missing-data case.
    electronic_none.bond_orders.push_back(
        {0, 1, 0.0, cov::DataProvenance::Derived});
    if (!cov::analyse_bonds(electronic_none).empty()) {
        std::cerr << "geometry fallback overrode an actual electronic zero\n";
        return 6;
    }

    // A partial matrix cannot suppress an unrelated missing neighbour. Once
    // that neighbour receives an explicit zero, the same geometry is vetoed.
    cov::Wavefunction partial_electronic;
    partial_electronic.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 1.40*cov::kAngstromToBohr, 0.0, 0.0},
        {"C", 6, 2.80*cov::kAngstromToBohr, 0.0, 0.0},
    };
    partial_electronic.bond_order_provenance=cov::DataProvenance::Derived;
    partial_electronic.bond_orders.push_back(
        {0, 1, 1.0, cov::DataProvenance::Derived});
    const auto partial_bonds=cov::analyse_bonds(partial_electronic);
    if (partial_bonds.size()!=2u ||
        std::count_if(partial_bonds.begin(),partial_bonds.end(),
            [](const cov::BondVisual& bond) {
                return bond.atom_a==1u && bond.atom_b==2u &&
                    bond.provenance==cov::DataProvenance::Unavailable;
            })!=1) {
        std::cerr << "partial matrix suppressed a missing neighbouring pair\n";
        return 33;
    }
    partial_electronic.bond_orders.push_back(
        {1, 2, 0.0, cov::DataProvenance::Derived});
    const auto partial_zero_bonds=cov::analyse_bonds(partial_electronic);
    if (partial_zero_bonds.size()!=1u ||
        partial_zero_bonds.front().atom_a!=0u ||
        partial_zero_bonds.front().atom_b!=1u) {
        std::cerr << "partial explicit zero lost its geometry veto\n";
        return 34;
    }

    // Delocalised Mayer coupling is not structural adjacency. Even when every
    // carbon pair carries a small positive index, benzene must remain a C6
    // perimeter rather than a complete graph with centre-crossing chords.
    cov::Wavefunction benzene;
    constexpr double benzene_radius_angstrom = side_angstrom;
    for (int i = 0; i < 6; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 6.0;
        benzene.atoms.push_back({
            "C", 6,
            benzene_radius_angstrom * std::cos(angle) * cov::kAngstromToBohr,
            benzene_radius_angstrom * std::sin(angle) * cov::kAngstromToBohr,
            0.0});
    }
    add_all_pair_mayer_records(benzene, 0.06);
    const auto benzene_bonds = cov::analyse_bonds(benzene);
    if (!check_carbon_graph(benzene, benzene_bonds, 6u,
                            std::vector<int>(6u, 2))) {
        std::cerr << "benzene electronic coupling created non-neighbour chords\n";
        return 7;
    }
    for (const auto& bond : benzene_bonds) {
        if (!bond.delocalised) {
            std::cerr << "benzene perimeter was not marked delocalised\n";
            return 8;
        }
    }
    const auto benzene_interactions = cov::build_interaction_graph(benzene);
    if (benzene_interactions.edges.size() != 6u ||
        !std::all_of(benzene_interactions.edges.begin(),
                     benzene_interactions.edges.end(),
                     [&](const cov::InteractionEdge& edge) {
                         return cov::interaction_visual_style(edge.kind, defaults) ==
                                cov::InteractionVisualStyle::OrdinaryBond;
                     })) {
        std::cerr << "benzene perimeter changed in the layered render graph\n";
        return 21;
    }

    // Two regular fused six-membered rings: the carbon skeleton has exactly
    // eleven edges. Only the two shared-edge atoms have degree three.
    cov::Wavefunction naphthalene;
    const double root3 = std::sqrt(3.0);
    const std::vector<std::pair<double, double>> naphthalene_xy = {
        {0.0, 0.5 * side_angstrom},
        {0.0, -0.5 * side_angstrom},
        {-0.5 * root3 * side_angstrom, side_angstrom},
        {-root3 * side_angstrom, 0.5 * side_angstrom},
        {-root3 * side_angstrom, -0.5 * side_angstrom},
        {-0.5 * root3 * side_angstrom, -side_angstrom},
        {0.5 * root3 * side_angstrom, side_angstrom},
        {root3 * side_angstrom, 0.5 * side_angstrom},
        {root3 * side_angstrom, -0.5 * side_angstrom},
        {0.5 * root3 * side_angstrom, -side_angstrom},
    };
    for (const auto [x_angstrom, y_angstrom] : naphthalene_xy) {
        naphthalene.atoms.push_back({
            "C", 6, x_angstrom * cov::kAngstromToBohr,
            y_angstrom * cov::kAngstromToBohr, 0.0});
    }
    add_all_pair_mayer_records(naphthalene, 0.06);
    const auto naphthalene_bonds = cov::analyse_bonds(naphthalene);
    const std::vector<int> naphthalene_degrees = {3, 3, 2, 2, 2, 2, 2, 2, 2, 2};
    if (!check_carbon_graph(naphthalene, naphthalene_bonds, 11u,
                            naphthalene_degrees)) {
        std::cerr << "naphthalene electronic coupling created non-neighbour chords\n";
        return 9;
    }
    for (const auto& bond : naphthalene_bonds) {
        if (!bond.delocalised) {
            std::cerr << "naphthalene fused-ring perimeter was not marked delocalised\n";
            return 10;
        }
    }
    const auto naphthalene_interactions = cov::build_interaction_graph(naphthalene);
    if (naphthalene_interactions.edges.size() != 11u ||
        !std::all_of(naphthalene_interactions.edges.begin(),
                     naphthalene_interactions.edges.end(),
                     [&](const cov::InteractionEdge& edge) {
                         return cov::interaction_visual_style(edge.kind, defaults) ==
                                cov::InteractionVisualStyle::OrdinaryBond;
                     })) {
        std::cerr << "naphthalene fused perimeter changed in the layered render graph\n";
        return 22;
    }

    // All real elements supported by the parser must have a defined radius.
    // In particular, 4d/5d centres must never fall back to the old 0.85 A
    // sentinel: doing so removed long but electronically supported
    // coordination bonds.
    if (std::abs(cov::covalent_radius_angstrom(40) - 1.75) > 1.0e-12 ||
        std::abs(cov::covalent_radius_angstrom(72) - 1.75) > 1.0e-12) {
        std::cerr << "Zr/Hf covalent-radius coverage regressed\n";
        return 11;
    }
    for (int z = 1; z <= 118; ++z) {
        if (cov::covalent_radius_angstrom(z) <= 0.0) {
            std::cerr << "missing covalent radius for atomic number " << z << '\n';
            return 12;
        }
    }

    // Synthetic pentagonal-bipyramidal ZrF7: two axial 2.10 A contacts and
    // five equatorial 2.15 A contacts, each backed by a strong Mayer index.
    cov::Wavefunction zrf7;
    zrf7.atoms.push_back({"Zr", 40, 0.0, 0.0, 0.0});
    zrf7.atoms.push_back({"F", 9, 0.0, 0.0, 2.10 * cov::kAngstromToBohr});
    zrf7.atoms.push_back({"F", 9, 0.0, 0.0, -2.10 * cov::kAngstromToBohr});
    for (int i = 0; i < 5; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 5.0;
        zrf7.atoms.push_back({
            "F", 9,
            2.15 * std::cos(angle) * cov::kAngstromToBohr,
            2.15 * std::sin(angle) * cov::kAngstromToBohr,
            0.0});
    }
    zrf7.bond_order_provenance = cov::DataProvenance::Derived;
    for (std::uint32_t ligand = 1; ligand <= 7; ++ligand) {
        zrf7.bond_orders.push_back({0, ligand, 0.80, cov::DataProvenance::Derived});
    }
    if (!has_central_degree(zrf7, 7u)) {
        std::cerr << "electronically supported ZrF7 first shell lost bonds\n";
        return 13;
    }

    // Synthetic D5d pentagonal-antiprismatic HfO10. Every donor lies 2.50 A
    // from Hf and has the same moderate Mayer evidence as an aqua complex.
    cov::Wavefunction hfo10;
    hfo10.atoms.push_back({"Hf", 72, 0.0, 0.0, 0.0});
    constexpr double ring_z_angstrom = 1.0;
    const double ring_radius_angstrom =
        std::sqrt(2.50 * 2.50 - ring_z_angstrom * ring_z_angstrom);
    for (int ring = 0; ring < 2; ++ring) {
        const double z_angstrom = ring == 0 ? ring_z_angstrom : -ring_z_angstrom;
        const double offset = ring == 0 ? 0.0 : pi / 5.0;
        for (int i = 0; i < 5; ++i) {
            const double angle = 2.0 * pi * static_cast<double>(i) / 5.0 + offset;
            hfo10.atoms.push_back({
                "O", 8,
                ring_radius_angstrom * std::cos(angle) * cov::kAngstromToBohr,
                ring_radius_angstrom * std::sin(angle) * cov::kAngstromToBohr,
                z_angstrom * cov::kAngstromToBohr});
        }
    }
    hfo10.bond_order_provenance = cov::DataProvenance::Derived;
    for (std::uint32_t ligand = 1; ligand <= 10; ++ligand) {
        hfo10.bond_orders.push_back({0, ligand, 0.37, cov::DataProvenance::Derived});
    }
    if (!has_central_degree(hfo10, 10u)) {
        std::cerr << "electronically supported HfO10 first shell lost bonds\n";
        return 14;
    }

    // An equilibrium one-electron H2+ bond is longer than the ordinary H--H
    // covalent-radius envelope but has direct Mayer evidence. The adaptive
    // envelope must keep it as one molecular fragment without relaxing weak
    // remote couplings in general.
    cov::Wavefunction h2_plus;
    h2_plus.atoms = {
        {"H", 1, -0.53 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1,  0.53 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    h2_plus.bond_order_provenance = cov::DataProvenance::Derived;
    h2_plus.bond_orders.push_back({0, 1, 0.50, cov::DataProvenance::Derived});
    if (cov::analyse_bonds(h2_plus).size() != 1u) {
        std::cerr << "electronically supported one-electron H2+ bond was lost\n";
        return 15;
    }

    // A complete heavy-element table must not make genuinely remote Mayer
    // coupling drawable. The existing adjacency envelope remains authoritative.
    cov::Wavefunction remote_hf_o;
    remote_hf_o.atoms = {
        {"Hf", 72, 0.0, 0.0, 0.0},
        {"O", 8, 5.0 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    remote_hf_o.bond_order_provenance = cov::DataProvenance::Derived;
    remote_hf_o.bond_orders.push_back({0, 1, 0.06, cov::DataProvenance::Derived});
    if (!cov::analyse_bonds(remote_hf_o).empty()) {
        std::cerr << "remote heavy-element Mayer coupling became a structural bond\n";
        return 16;
    }

    if (argc == 3 || argc == 5) {
        if (!check_real_topology(argv[1], 12u, 6u)) {
            std::cerr << "real benzene topology regression\n";
            return 17;
        }
        if (!check_real_topology(argv[2], 19u, 11u)) {
            std::cerr << "real naphthalene topology regression\n";
            return 18;
        }
    }
    if (argc == 5) {
        if (!check_real_topology(argv[3], 10u, 5u)) {
            std::cerr << "real cyclopentadienyl topology regression\n";
            return 25;
        }
        if (!check_real_topology(argv[4], 19u, 11u)) {
            std::cerr << "real azulene topology regression\n";
            return 26;
        }
    }

    std::cout << "molecule_style_smoke ok\n";
    return 0;
}
