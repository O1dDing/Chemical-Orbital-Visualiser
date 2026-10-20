#include "cov/wavefunction_io.hpp"
#include "cov/molecule_style.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/pi_topology_evidence.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// Read-only numerical instrumentation. Production parsing and analysis are
// called unchanged; the independent Python reference consumes this stdout.
template<class Range> void numeric_array(const Range& values) {
    std::cout << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) std::cout << ',';
        first = false;
        if (std::isfinite(static_cast<double>(value))) std::cout << value;
        else std::cout << "null";
    }
    std::cout << ']';
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3 || (argc == 3 && std::string(argv[2]) != "--interactions")) return 1;
    const bool detailed_interactions = argc == 3;
    try {
        cov::WavefunctionParseOptions options;
        options.max_atoms = 1000;
        const auto wavefunction = cov::parse_wavefunction(argv[1], options);
        std::cout << std::setprecision(17) << "{\"nbasis\":" << wavefunction.basis_count
                  << ",\"overlap\":";
        numeric_array(wavefunction.ao_overlap);
        std::cout << ",\"numerical_diagnostics\":";
        cov::write_numerical_diagnostics_json(std::cout,wavefunction);
        std::cout << ",\"density_evidence\":";
        cov::write_density_evidence_json(std::cout,wavefunction);
        std::cout << ",\"pi_topology_assignments\":";
        cov::write_pi_topology_assignments_json(std::cout,wavefunction);
        std::cout << ",\"point_group_source_records\":[";
        for (std::size_t i=0;i<wavefunction.point_group_source_records.size();++i) {
            if (i) std::cout << ',';
            const auto& record=wavefunction.point_group_source_records[i];
            std::cout << "{\"field\":";
            cov::numerical_json::string(std::cout,record.field);
            std::cout << ",\"raw_record\":";
            cov::numerical_json::string(std::cout,record.raw_record);
            std::cout << ",\"value\":";
            cov::numerical_json::string(std::cout,record.value);
            std::cout << ",\"field_text\":";
            cov::numerical_json::string(std::cout,record.field_text);
            std::cout << ",\"status\":";
            cov::numerical_json::string(std::cout,record.status);
            std::cout << ",\"line_number\":" << record.line_number
                      << ",\"job_segment\":" << record.job_segment
                      << ",\"valid\":" << (record.valid?"true":"false") << '}';
        }
        std::cout << ']';
        std::cout << ",\"derived_symmetry_assignments\":[";
        for (std::size_t i=0;i<wavefunction.derived_orbital_symmetry_assignments.size();++i) {
            if (i) std::cout << ',';
            const auto& record=wavefunction.derived_orbital_symmetry_assignments[i];
            std::cout << "{\"point_group\":";
            cov::numerical_json::string(std::cout,record.point_group);
            std::cout << ",\"label\":";
            cov::numerical_json::string(std::cout,record.label);
            std::cout << ",\"orbital_indices\":";numeric_array(record.orbital_indices);
            std::cout << ",\"subspace_retention\":";
            cov::numerical_json::number(std::cout,record.subspace_retention);
            std::cout << ",\"maximum_character_error\":";
            cov::numerical_json::number(std::cout,record.maximum_character_error);
            std::cout << ",\"centre_bohr\":";numeric_array(record.centre_bohr);
            std::cout << ",\"axis_frame\":";
            if (record.axes_available) {
                std::cout << "{\"principal\":";numeric_array(record.principal_axis);
                std::cout << ",\"secondary\":";numeric_array(record.secondary_axis);
                std::cout << ",\"convention\":";
                cov::numerical_json::string(std::cout,record.axis_convention);
                std::cout << '}';
            } else std::cout << "null";
            std::cout << '}';
        }
        std::cout << ']';
        std::cout << ",\"producer_overlap\":";
        numeric_array(wavefunction.producer_ao_overlap);
        std::cout << ",\"ao_transform\":[";
        for (std::size_t i=0;i<wavefunction.gaussian_ao_transform.size();++i) {
            if (i) std::cout << ',';
            const auto& entry=wavefunction.gaussian_ao_transform[i];
            numeric_array(std::array<double,4>{static_cast<double>(i),static_cast<double>(entry.source_index),
                          entry.basis_scale,entry.coefficient_scale});
        }
        std::cout << ']';
        std::cout << ",\"bond_orders\":[";
        bool first = true;
        for (const auto& bond : wavefunction.bond_orders) {
            if (!first) std::cout << ',';
            first = false;
            numeric_array(std::array<double,3>{static_cast<double>(bond.atom_a),
                          static_cast<double>(bond.atom_b), bond.mayer_order});
        }
        std::cout << "],\"orbitals\":[";
        first = true;
        for (const auto& orbital : wavefunction.orbitals) {
            if (!first) std::cout << ',';
            first=false;
            std::cout << "{\"energy_hartree\":" << orbital.energy_hartree
                      << ",\"occupation\":" << orbital.occupation
                      << ",\"spin\":" << static_cast<int>(orbital.spin)
                      << ",\"coefficients\":";
            numeric_array(orbital.coefficients);
            std::cout << ",\"gaussian_source_coefficients\":";
            numeric_array(orbital.gaussian_source_coefficients);
            std::cout << '}';
        }
        std::cout << "],\"chemistry\":[";
        first = true;
        for (const auto& orbital : wavefunction.orbitals) {
            if (!first) std::cout << ',';
            first = false;
            const auto& c = orbital.chemistry;
            std::cout << "{\"available\":" << (c.available ? "true" : "false")
                      << ",\"spin\":" << static_cast<int>(orbital.spin)
                      << ",\"symmetry\":" << std::quoted(orbital.symmetry)
                      << ",\"symmetry_provenance\":" << static_cast<int>(orbital.symmetry_provenance)
                      << ",\"bonding_status\":" << static_cast<int>(c.bonding.status)
                      << ",\"channel_status\":" << static_cast<int>(c.channel.status)
                      << ",\"bonding\":";
            numeric_array(std::array<double,4>{c.bonding.bonding,c.bonding.antibonding,
                          c.bonding.nonbonding,c.bonding.undetermined});
            std::cout << ",\"channel\":";
            numeric_array(std::array<double,5>{c.channel.sigma,c.channel.pi,
                          c.channel.delta,c.channel.phi,c.channel.undetermined});
            std::cout << ",\"manifold\":";
            numeric_array(std::array<double,4>{c.deep_core_weight,c.semicore_weight,
                          c.valence_weight,c.unresolved_weight});
            if (detailed_interactions) {
                std::cout << ",\"interactions\":[";
                bool first_pair = true;
                for (const auto& pair : c.interactions) {
                    if (!first_pair) std::cout << ',';
                    first_pair = false;
                    std::cout << "{\"atoms\":[" << pair.atom_a << ',' << pair.atom_b
                              << "],\"overlap\":" << pair.overlap_character
                              << ",\"channel\":";
                    numeric_array(std::array<double,5>{pair.channel.sigma,pair.channel.pi,
                                  pair.channel.delta,pair.channel.phi,pair.channel.undetermined});
                    std::cout << ",\"bonding\":";
                    numeric_array(std::array<double,4>{pair.bonding.bonding,pair.bonding.antibonding,
                                  pair.bonding.nonbonding,pair.bonding.undetermined});
                    std::cout << '}';
                }
                std::cout << ']';
            }
            std::cout << '}';
        }
        std::cout << "],\"atoms\":[";
        first = true;
        for (const auto& atom : wavefunction.atoms) {
            if (!first) std::cout << ',';
            first = false;
            numeric_array(std::array<double,4>{static_cast<double>(atom.atomic_number),
                          atom.x,atom.y,atom.z});
        }
        // Use exactly the renderer's default graph/evidence/settings contract.
        const auto graph = cov::build_interaction_graph(wavefunction);
        const auto bonds = cov::analyse_bonds(wavefunction);
        const cov::MoleculeRenderSettings settings;
        std::cout << "],\"render_edges\":[";
        first = true;
        for (const auto& edge : graph.edges) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"atoms\":[" << edge.atom_a << ',' << edge.atom_b
                      << "],\"kind\":" << std::quoted(cov::interaction_kind_name(edge.kind))
                      << ",\"default_style\":" << static_cast<int>(cov::interaction_visual_style(edge.kind,settings))
                      << ",\"distance_bohr\":" << edge.distance_bohr
                      << ",\"radius_ratio\":" << edge.covalent_radius_ratio
                      << ",\"mayer\":" << edge.mayer_order
                      << ",\"confidence\":" << edge.confidence << '}';
        }
        std::cout << "],\"bond_visuals\":[";
        first = true;
        for (const auto& bond : bonds) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"atoms\":[" << bond.atom_a << ',' << bond.atom_b
                      << "],\"delocalised\":" << (bond.delocalised ? "true" : "false")
                      << ",\"order\":" << bond.bond_order << '}';
        }
        // Read-only compact builder, using the live UI's selection sentinel.
        cov::MODiagramOptions diagram_options;
        diagram_options.mode = cov::preferred_compact_mo_diagram_mode(wavefunction,true);
        diagram_options.hide_ligand_centred_intermediates = true;
        diagram_options.selected_index = wavefunction.orbitals.size();
        diagram_options.nonlinear_minimum_gap_weight = 0.055;
        const auto diagram = cov::build_mo_diagram_data(wavefunction,diagram_options);
        std::cout << "],\"compact_mode\":" << std::quoted(cov::mo_diagram_mode_name(diagram.mode))
                  << ",\"point_group\":" << std::quoted(diagram.electronic_state.point_group_detected)
                  << ",\"compact_rows\":[";
        first = true;
        for (const auto& row : diagram.levels) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"members\":";
            numeric_array(row.member_indices);
            std::cout << ",\"counterparts\":";
            numeric_array(row.member_spin_counterparts);
            std::cout << ",\"symmetry\":" << std::quoted(row.metadata.symmetry_view.label)
                      << ",\"source_symmetry\":" << std::quoted(row.metadata.symmetry)
                      << ",\"symmetry_explanation\":" << cov::orbital_symmetry_json(row.metadata.symmetry_view)
                      << ",\"energy\":" << row.layout_energy_hartree
                      << ",\"occupation\":" << row.total_occupation
                      << ",\"degeneracy\":" << row.metadata.degeneracy_size
                      << ",\"family\":" << std::quoted(row.annotation.family)
                      << ",\"bonding\":" << std::quoted(cov::bonding_class_name(row.annotation.bonding_class))
                      << ",\"pi_atoms\":";
            numeric_array(row.annotation.delocalised_pi.atom_indices);
            std::cout << ",\"pi_orbitals\":";
            numeric_array(row.annotation.delocalised_pi.orbital_indices);
            std::cout << '}';
        }
        std::cout << "],\"pi_partner_candidates\":" << cov::pi_partner_candidates_json(diagram.pi_partner_candidates)
                  << ",\"pi_interactions\":" << cov::orbital_energy_gap_array_json(diagram.pi_interactions)
                  << ",\"crystal_field_gaps\":" << cov::orbital_energy_gap_array_json(diagram.crystal_field_gaps) << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
