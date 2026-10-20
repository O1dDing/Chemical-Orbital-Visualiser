#include "cov/mo_diagram.hpp"
#include "cov/ligand_field.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/molecule_style.hpp"
#include "cov/pi_topology.hpp"
#include "cov/pi_pair_evidence.hpp"
#include "cov/point_group_catalog.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

cov::Wavefunction synthetic_oh_environment() {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(7u);
    wavefunction.atoms[0].symbol="Cr";
    wavefunction.atoms[0].atomic_number=24;
    const std::array<std::array<double,3>,6> ligand_positions{{
        {{4.0,0.0,0.0}},{{-4.0,0.0,0.0}},
        {{0.0,4.0,0.0}},{{0.0,-4.0,0.0}},
        {{0.0,0.0,4.0}},{{0.0,0.0,-4.0}},
    }};
    for (std::size_t i=0;i<ligand_positions.size();++i) {
        auto& atom=wavefunction.atoms[i+1u];
        atom.symbol="F";
        atom.atomic_number=9;
        atom.x=ligand_positions[i][0];
        atom.y=ligand_positions[i][1];
        atom.z=ligand_positions[i][2];
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+1u);
        bond.mayer_order=0.20;
        wavefunction.bond_orders.push_back(bond);
    }
    return wavefunction;
}

cov::Wavefunction synthetic_td_environment() {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(5u);
    wavefunction.atoms[0].symbol="Cr";
    wavefunction.atoms[0].atomic_number=24;
    const std::array<std::array<double,3>,4> ligand_positions{{
        {{2.3,2.3,2.3}},{{2.3,-2.3,-2.3}},
        {{-2.3,2.3,-2.3}},{{-2.3,-2.3,2.3}},
    }};
    for (std::size_t i=0;i<ligand_positions.size();++i) {
        auto& atom=wavefunction.atoms[i+1u];
        atom.symbol="F";
        atom.atomic_number=9;
        atom.x=ligand_positions[i][0];
        atom.y=ligand_positions[i][1];
        atom.z=ligand_positions[i][2];
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+1u);
        bond.mayer_order=0.20;
        wavefunction.bond_orders.push_back(bond);
    }
    return wavefunction;
}

bool validate_radial_first_shell_retry() {
    auto wavefunction=synthetic_oh_environment();
    for (auto& bond:wavefunction.bond_orders) bond.mayer_order=0.025;
    const std::size_t near_count=wavefunction.atoms.size();
    wavefunction.atoms.resize(13u);
    for (std::size_t i=0;i<6u;++i) {
        const auto& near=wavefunction.atoms[i+1u];
        auto& remote=wavefunction.atoms[i+7u];
        remote.symbol="O";
        remote.atomic_number=8;
        remote.x=2.0*near.x;
        remote.y=2.0*near.y;
        remote.z=2.0*near.z;
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+7u);
        bond.mayer_order=0.08;
        wavefunction.bond_orders.push_back(bond);
    }
    const auto environment=cov::analyse_ligand_field_environment(wavefunction);
    if (near_count!=7u || environment.local_point_group()!="Oh" ||
        environment.ligand_atoms.size()!=6u ||
        std::any_of(environment.ligand_atoms.begin(),
                    environment.ligand_atoms.end(),[&](const auto atom) {
            return atom==0u || atom>=near_count ||
                   wavefunction.atoms[atom].atomic_number!=9;
        })) {
        std::cerr<<"radial-first low-Mayer Oh retry selected a remote atom\n";
        return false;
    }
    return true;
}

bool validate_local_symmetry_missing_markers() {
    auto make_block=[](const int angular_momentum,const double weight) {
        auto wavefunction=synthetic_oh_environment();
        for (std::size_t i=0;i<3u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.10+0.001*static_cast<double>(i);
            orbital.occupation=1.0f;
            orbital.chemistry.available=true;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=angular_momentum;
            contribution.weight=weight;
            orbital.chemistry.ao_contributions.push_back(contribution);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
        return wavefunction;
    };
    const std::array<std::string,3> missing{{"","?","N/A"}};
    auto metadata_for=[&](const auto& wavefunction) {
        std::vector<cov::OrbitalMetadata> metadata(wavefunction.orbitals.size());
        for (std::size_t i=0;i<metadata.size();++i) {
            metadata[i].energy_hartree=wavefunction.orbitals[i].energy_hartree;
            metadata[i].degeneracy_size=1u;
            metadata[i].symmetry=missing[i%missing.size()];
        }
        return metadata;
    };

    auto above=make_block(2,0.021);
    auto above_metadata=metadata_for(above);
    cov::apply_local_ligand_field_symmetry(above,above_metadata);
    if (!std::all_of(above_metadata.begin(),above_metadata.end(),[](const auto& item) {
            return item.symmetry_view.label=="T2g" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::LocalDimensionCandidate;
        })) {
        std::cerr<<"missing local input did not retain a scoped dimension candidate\n";
        return false;
    }
    for (std::size_t i=0;i<above_metadata.size();++i) {
        if (above_metadata[i].symmetry!=missing[i]) {
            std::cerr<<"local dimension candidate overwrote the original label\n";
            return false;
        }
    }

    auto below=make_block(2,0.019);
    auto below_metadata=metadata_for(below);
    cov::apply_local_ligand_field_symmetry(below,below_metadata);
    if (std::any_of(below_metadata.begin(),below_metadata.end(),[](const auto& item) {
            return item.symmetry_view.label=="T2g";
        })) {
        std::cerr<<"sub-threshold metal-d noise was grouped\n";
        return false;
    }

    auto low_p=make_block(1,0.021);
    auto low_p_metadata=metadata_for(low_p);
    cov::apply_local_ligand_field_symmetry(low_p,low_p_metadata);
    if (std::any_of(low_p_metadata.begin(),low_p_metadata.end(),[](const auto& item) {
            return item.symmetry_view.label=="T1u";
        })) {
        std::cerr<<"low metal-p noise used the metal-d grouping floor\n";
        return false;
    }
    return true;
}

bool validate_ligand_only_noise_preserves_producer_irrep() {
    auto wavefunction=synthetic_oh_environment();
    wavefunction.point_group_detected="Oh";
    wavefunction.basis_count=5u;
    cov::Shell metal_d;
    metal_d.atom_index=0u;
    metal_d.basis_offset=0u;
    metal_d.angular_momentum=2u;
    metal_d.pure=1u;
    wavefunction.shells.push_back(metal_d);

    for (std::size_t member=0;member<3u;++member) {
        cov::MolecularOrbital orbital;
        orbital.energy_hartree=-0.08;
        orbital.occupation=1.0f;
        orbital.symmetry="T1g";
        orbital.symmetry_provenance=cov::DataProvenance::Producer;
        orbital.coefficients.assign(wavefunction.basis_count,0.0f);
        // A deliberately misleading, numerically tiny central-metal dz2
        // coefficient would be classified as Eg if it were normalised on its
        // own.  The chemically meaningful population is far below the d-floor.
        orbital.coefficients[0]=1.0e-5f;
        orbital.chemistry.available=true;
        orbital.chemistry.valence_manifold=true;
        orbital.chemistry.valence_weight=1.0;
        orbital.chemistry.unresolved_weight=0.0;

        cov::OrbitalAOContribution metal_noise;
        metal_noise.atom_index=0u;
        metal_noise.angular_momentum=2;
        metal_noise.weight=1.0e-12;
        orbital.chemistry.ao_contributions.push_back(metal_noise);

        cov::OrbitalAOContribution ligand;
        ligand.atom_index=static_cast<std::uint32_t>(member+1u);
        ligand.angular_momentum=1;
        ligand.weight=0.85;
        orbital.chemistry.ao_contributions.push_back(ligand);
        wavefunction.orbitals.push_back(std::move(orbital));
    }

    const auto environment=cov::analyse_ligand_field_environment(wavefunction);
    const std::array<std::size_t,3> members{{0u,1u,2u}};
    const auto noise_only=cov::classify_local_metal_irrep(
        wavefunction,members,environment.metal_atom,
        environment.local_point_group(),environment.rotation_reference_to_input);
    if (noise_only) {
        std::cerr<<"missing metric and repeated tiny coefficients produced an irrep proof\n";
        return false;
    }

    cov::MODiagramOptions options;
    options.selected_index=0u;
    const auto metadata=cov::build_orbital_metadata(
        wavefunction,options.selected_index,options.degeneracy,options.filter);
    if (metadata.size()!=3u ||
        !std::all_of(metadata.begin(),metadata.end(),[](const auto& item) {
            return item.symmetry=="T1g" && item.degeneracy_size==3u &&
                   item.symmetry_view.label=="T1g" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::Producer &&
                   item.symmetry_view.point_group.empty();
        })) {
        std::cerr<<"ligand-only metal noise overwrote producer T1g metadata\n";
        return false;
    }

    const auto diagram=cov::build_mo_diagram_data(wavefunction,options);
    const auto level=std::find_if(
        diagram.levels.begin(),diagram.levels.end(),[](const auto& item) {
            return item.member_indices.size()==3u &&
                   item.member_indices.front()==0u;
        });
    const auto* point_group=cov::find_point_group(diagram.ligand_field_point_group);
    if (level==diagram.levels.end() || level->metadata.symmetry!="T1g" ||
        point_group==nullptr) {
        std::cerr<<"ligand-only noise produced an illegal local-irrep row dimension\n";
        return false;
    }
    const auto irrep=std::find_if(
        point_group->irreps.begin(),point_group->irreps.end(),
        [&](const auto& item) {return item.label==level->metadata.symmetry;});
    if (irrep==point_group->irreps.end() ||
        irrep->dimension!=level->member_indices.size()) {
        std::cerr<<"ligand-only noise produced an illegal local-irrep row dimension\n";
        return false;
    }
    return true;
}

bool validate_resolved_five_d_runs() {
    const auto make_case=[](cov::Wavefunction wavefunction,
                            const std::array<double,5>& energies) {
        std::vector<cov::OrbitalMetadata> metadata(energies.size());
        for (std::size_t i=0;i<energies.size();++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energies[i];
            orbital.occupation=1.0f;
            orbital.chemistry.available=true;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=2;
            contribution.weight=0.30;
            orbital.chemistry.ao_contributions.push_back(contribution);
            wavefunction.orbitals.push_back(std::move(orbital));
            metadata[i].energy_hartree=energies[i];
            metadata[i].degeneracy_size=1u;
            metadata[i].symmetry=i%2u==0u?"N/A":"?";
        }
        cov::apply_local_ligand_field_symmetry(wavefunction,metadata);
        return metadata;
    };
    const auto td=make_case(
        synthetic_td_environment(),{{-0.1000,-0.0998,-0.0968,-0.0966,-0.0964}});
    const bool td_ok=std::all_of(td.begin(),td.begin()+2u,[](const auto& item) {
            return item.symmetry_view.label=="E" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::LocalDimensionCandidate;
        }) && std::all_of(td.begin()+2u,td.end(),[](const auto& item) {
            return item.symmetry_view.label=="T2" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::LocalDimensionCandidate;
        });
    const auto oh=make_case(
        synthetic_oh_environment(),{{-0.1000,-0.0998,-0.0996,-0.0966,-0.0964}});
    const bool oh_ok=std::all_of(oh.begin(),oh.begin()+3u,[](const auto& item) {
            return item.symmetry_view.label=="T2g" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::LocalDimensionCandidate;
        }) && std::all_of(oh.begin()+3u,oh.end(),[](const auto& item) {
            return item.symmetry_view.label=="Eg" &&
                   item.symmetry_view.origin==cov::OrbitalSymmetryOrigin::LocalDimensionCandidate;
        });
    if (!td_ok || !oh_ok) {
        std::cerr<<"resolved weak-field five-d run was greedily mispartitioned\n";
        return false;
    }
    for (std::size_t i=0;i<td.size();++i) {
        const std::string raw=i%2u==0u?"N/A":"?";
        if (td[i].symmetry!=raw || oh[i].symmetry!=raw) {
            std::cerr<<"five-d dimension candidate overwrote a source label\n";
            return false;
        }
    }
    return true;
}

bool validate_minimal_compact_ligand_field_framework() {
    const auto add_group=[](cov::Wavefunction& wavefunction,
                            const std::string& symmetry,
                            const std::size_t degeneracy,
                            const double energy,
                            const int metal_angular_momentum,
                            const double metal_weight,
                            const double ligand_p_weight,
                            const bool sigma) {
        const std::size_t ligand_count=wavefunction.atoms.size()-1u;
        for (std::size_t member=0;member<degeneracy;++member) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energy;
            orbital.occupation=2.0f;
            orbital.symmetry=symmetry;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.valence_weight=1.0;
            orbital.chemistry.unresolved_weight=0.0;
            if (metal_weight>0.0) {
                cov::OrbitalAOContribution metal;
                metal.atom_index=0u;
                metal.angular_momentum=metal_angular_momentum;
                metal.weight=metal_weight;
                orbital.chemistry.ao_contributions.push_back(metal);
            }
            cov::OrbitalAOContribution ligand;
            ligand.atom_index=static_cast<std::uint32_t>(
                1u+member%ligand_count);
            ligand.angular_momentum=1;
            ligand.weight=ligand_p_weight;
            orbital.chemistry.ao_contributions.push_back(ligand);

            cov::OrbitalPairInteraction interaction;
            interaction.atom_a=0u;
            interaction.atom_b=ligand.atom_index;
            interaction.overlap_character=0.06;
            interaction.channel.sigma=sigma?1.0:0.0;
            interaction.channel.pi=sigma?0.0:1.0;
            interaction.channel.undetermined=0.0;
            orbital.chemistry.interactions.push_back(interaction);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
    };
    const auto diagram_symmetries=[](const cov::Wavefunction& wavefunction,
                                     const bool compact) {
        cov::MODiagramOptions options;
        options.hide_ligand_centred_intermediates=compact;
        options.selected_index=wavefunction.orbitals.size();
        const auto data=cov::build_mo_diagram_data(wavefunction,options);
        std::set<std::string> result;
        for (const auto& level:data.levels) {
            result.insert(level.metadata.symmetry);
        }
        return std::pair{data.ligand_field_point_group,std::move(result)};
    };

    auto oh=synthetic_oh_environment();
    oh.point_group_detected="Oh";
    add_group(oh,"Eg",2u,-0.20,2,0.70,0.20,true);
    add_group(oh,"T2g",3u,-0.10,2,0.70,0.20,false);
    add_group(oh,"T1g",3u,-0.05,2,0.00,0.80,false);
    add_group(oh,"T1u",3u,0.05,1,0.25,0.30,true);
    const auto [oh_expanded_group,oh_expanded]=diagram_symmetries(oh,false);
    const auto [oh_group,oh_symmetries]=diagram_symmetries(oh,true);
    if (oh_expanded_group!="Oh" || oh_expanded.count("T1g")==0u ||
        oh_group!="Oh" || oh_symmetries.count("Eg")==0u ||
        oh_symmetries.count("T2g")==0u ||
        oh_symmetries.count("T1u")==0u ||
        oh_symmetries.count("T1g")!=0u) {
        std::cerr<<"minimal compact Oh framework pruning failed\n";
        return false;
    }

    auto td=synthetic_td_environment();
    td.point_group_detected="Td";
    add_group(td,"E",2u,-0.20,2,0.70,0.20,false);
    add_group(td,"T2",3u,-0.10,2,0.70,0.20,true);
    add_group(td,"T1",3u,-0.05,2,0.00,0.80,false);
    const auto [td_expanded_group,td_expanded]=diagram_symmetries(td,false);
    const auto [td_group,td_symmetries]=diagram_symmetries(td,true);
    if (td_expanded_group!="Td" || td_expanded.count("T1")==0u ||
        td_group!="Td" || td_symmetries.count("E")==0u ||
        td_symmetries.count("T2")==0u ||
        td_symmetries.count("T1")!=0u) {
        std::cerr<<"minimal compact Td framework pruning failed\n";
        return false;
    }
    return true;
}

cov::Wavefunction synthetic_oh_pi_pair(const int donor_z,
                                       const int external_z,
                                       const double external_mayer,
                                       const double lower_metal_d,
                                       const double lower_ligand_p,
                                       const double upper_metal_d,
                                       const double upper_ligand_p) {
    auto wavefunction=synthetic_oh_environment();
    wavefunction.point_group_detected="Oh";
    for (std::size_t i=1u;i<wavefunction.atoms.size();++i) {
        wavefunction.atoms[i].atomic_number=donor_z;
        wavefunction.atoms[i].symbol=donor_z==6?"C":(donor_z==7?"N":"X");
    }
    if (external_z>0) {
        wavefunction.atoms.resize(13u);
        for (std::size_t i=0;i<6u;++i) {
            const auto& direct=wavefunction.atoms[i+1u];
            auto& external=wavefunction.atoms[i+7u];
            external.atomic_number=external_z;
            external.symbol=external_z==6?"C":(external_z==7?"N":"X");
            external.x=1.7*direct.x;
            external.y=1.7*direct.y;
            external.z=1.7*direct.z;
            cov::BondOrderRecord bond;
            bond.atom_a=static_cast<std::uint32_t>(i+1u);
            bond.atom_b=static_cast<std::uint32_t>(i+7u);
            bond.mayer_order=external_mayer;
            wavefunction.bond_orders.push_back(bond);
        }
    }
    const auto add_group=[&](const double energy,
                             const float occupation,
                             const double metal_d,
                             const double ligand_p,
                             const double overlap) {
        for (std::size_t member=0;member<3u;++member) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energy;
            orbital.occupation=occupation;
            orbital.symmetry="T2g";
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.confidence=1.0;
            cov::OrbitalAOContribution metal;
            metal.atom_index=0u;
            metal.angular_momentum=2;
            metal.weight=metal_d;
            orbital.chemistry.ao_contributions.push_back(metal);
            cov::OrbitalAOContribution ligand;
            ligand.atom_index=static_cast<std::uint32_t>(member+1u);
            ligand.angular_momentum=1;
            ligand.weight=ligand_p;
            orbital.chemistry.ao_contributions.push_back(ligand);
            cov::OrbitalPairInteraction interaction;
            interaction.atom_a=0u;
            interaction.atom_b=static_cast<std::uint32_t>(member+1u);
            interaction.overlap_character=overlap;
            interaction.channel.pi=1.0;
            interaction.channel.undetermined=0.0;
            interaction.channel.dominant=cov::OrbitalAngularFamily::Pi;
            interaction.channel.status=cov::ChemistryStatus::Determined;
            orbital.chemistry.interactions.push_back(interaction);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
    };
    add_group(-0.20,2.0f,lower_metal_d,lower_ligand_p,1.0);
    add_group(0.10,0.0f,upper_metal_d,upper_ligand_p,-1.0);
    return wavefunction;
}

bool validate_pi_topology_and_two_sided_composition() {
    cov::MODiagramOptions options;
    options.selected_index=0u;
    options.max_levels=0u;
    const auto unverified_candidate=[](const cov::MODiagramData& data) {
        return data.pi_interactions.empty() &&
            std::any_of(data.pi_partner_candidates.begin(),data.pi_partner_candidates.end(),
                [](const auto& candidate) {
                    return candidate.input_valid && !candidate.accepted && !candidate.weak &&
                        candidate.direction==cov::PiPairDirection::Unresolved &&
                        candidate.channel.channel_id.empty() &&
                        !candidate.channel.frozen_operator.available &&
                        candidate.detail=="composition-only-candidate-needs-verified-local-channel-and-calibrated-sensitivity";
                });
    };
    const auto one_sided=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(6,7,1.50,0.70,0.10,0.10,0.12),options);
    if (!one_sided.pi_interactions.empty()) {
        std::cerr<<"one-sided pi composition produced a strong pair\n";
        return false;
    }
    const auto alkyl=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(6,6,0.80,0.10,0.70,0.70,0.10),options);
    if (!unverified_candidate(alkyl)) {
        std::cerr<<"carbon composition certified a direction without a verified local operator\n";
        return false;
    }
    const auto terminal_n=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(7,0,0.0,0.10,0.70,0.70,0.10),options);
    if (!unverified_candidate(terminal_n)) {
        std::cerr<<"terminal nitrogen candidate lacked an explicit unverified-operator decision\n";
        return false;
    }
    // Controlled P-centred call-path evidence. Three external single bonds
    // trigger the catalogue's sigma-only prior. Neither orientation of the
    // supplied composition contrast establishes a physical direction without
    // an independently verified local operator.
    for (const bool donor:{true,false}) {
        auto phosphorus=synthetic_oh_pi_pair(15,0,0.0,
            donor?0.10:0.70,donor?0.70:0.10,
            donor?0.70:0.10,donor?0.10:0.70);
        for (std::size_t centre=1u;centre<=6u;++centre) {
            for (std::size_t substituent=0u;substituent<3u;++substituent) {
                const auto atom_index=phosphorus.atoms.size();
                auto atom=phosphorus.atoms[centre];
                atom.symbol="H";atom.atomic_number=1;
                atom.x+=1.0+0.1*static_cast<double>(substituent);
                phosphorus.atoms.push_back(atom);
                cov::BondOrderRecord bond;
                bond.atom_a=static_cast<std::uint32_t>(centre);
                bond.atom_b=static_cast<std::uint32_t>(atom_index);
                bond.mayer_order=0.80;
                phosphorus.bond_orders.push_back(bond);
            }
        }
        const auto data=cov::build_mo_diagram_data(phosphorus,options);
        if (!unverified_candidate(data) ||
            std::none_of(data.pi_partner_candidates.begin(),data.pi_partner_candidates.end(),
                [](const auto& candidate){return candidate.prior==cov::LigandPiPrior::SigmaOnly;})) {
            std::cerr<<"P-centred composition or catalogue prior invented a directed counterpart\n";
            return false;
        }
        for (auto& orbital:phosphorus.orbitals) {
            for (auto& interaction:orbital.chemistry.interactions) {
                interaction.channel.pi=0.0;interaction.channel.sigma=1.0;
                interaction.channel.dominant=cov::OrbitalAngularFamily::Sigma;
            }
        }
        if (!cov::build_mo_diagram_data(phosphorus,options).pi_interactions.empty()) {
            std::cerr<<"P-centred composition without a pi channel produced a pi partner\n";
            return false;
        }
    }
    return true;
}

bool validate_equal_count_unrestricted_collapse() {
    auto unrestricted=synthetic_oh_environment();
    unrestricted.basis_count=3u;
    unrestricted.alpha_electrons=3u;
    unrestricted.beta_electrons=3u;
    unrestricted.ao_overlap={1.0,0.0,0.0,
                             0.0,1.0,0.0,
                             0.0,0.0,1.0};
    for (const auto spin:{cov::Spin::Alpha,cov::Spin::Beta}) {
        for (std::size_t i=0;i<3u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.10+0.001*static_cast<double>(i);
            orbital.occupation=1.0f;
            orbital.spin=spin;
            orbital.coefficients.assign(3u,0.0f);
            orbital.coefficients[i]=1.0f;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.confidence=1.0;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=2;
            contribution.weight=0.60;
            orbital.chemistry.ao_contributions.push_back(contribution);
            unrestricted.orbitals.push_back(std::move(orbital));
        }
    }
    cov::MODiagramOptions options;
    options.selected_index=1u;
    options.max_levels=0u;
    const auto unrestricted_data=cov::build_mo_diagram_data(
        unrestricted,options);
    if (!unrestricted_data.spin_counterparts_collapsed ||
        unrestricted_data.spin_counterparts_partial ||
        unrestricted_data.spin_counterpart_pair_count==0u ||
        unrestricted_data.spin_counterpart_unmatched_visible!=0u) {
        std::cerr<<"equal-count unrestricted spin sets were not collapsed\n";
        return false;
    }

    auto restricted=unrestricted;
    restricted.orbitals.resize(3u);
    for (auto& orbital:restricted.orbitals) orbital.occupation=2.0f;
    const auto restricted_data=cov::build_mo_diagram_data(restricted,options);
    if (restricted_data.spin_counterpart_pair_count!=0u ||
        restricted_data.spin_counterparts_collapsed ||
        restricted_data.spin_counterparts_partial) {
        std::cerr<<"alpha-only restricted orbitals were treated as UHF\n";
        return false;
    }
    return true;
}

cov::Wavefunction synthetic_spin_matching_wavefunction(
    const std::size_t basis_count,const std::string& point_group) {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(1u);
    wavefunction.atoms[0].symbol="Cr";
    wavefunction.atoms[0].atomic_number=24;
    wavefunction.point_group_detected=point_group;
    wavefunction.basis_count=basis_count;
    wavefunction.ao_overlap.assign(basis_count*basis_count,0.0);
    for (std::size_t i=0;i<basis_count;++i) {
        wavefunction.ao_overlap[i*basis_count+i]=1.0;
    }
    return wavefunction;
}

void add_spin_matching_orbital(cov::Wavefunction& wavefunction,
                               const cov::Spin spin,
                               const double energy,
                               const std::string& symmetry,
                               const std::vector<double>& coefficients,
                               const int metal_angular_momentum) {
    cov::MolecularOrbital orbital;
    orbital.energy_hartree=energy;
    orbital.occupation=1.0f;
    orbital.spin=spin;
    orbital.symmetry=symmetry;
    // This matching fixture supplies a declared representation and common
    // axes. A bare producer point-group string does not define the label's
    // scope. Molecular symmetry detection itself is tested separately.
    orbital.symmetry_provenance=cov::DataProvenance::Derived;
    cov::DerivedOrbitalSymmetryAssignment assignment;
    assignment.point_group=wavefunction.point_group_detected;
    assignment.label=symmetry;
    assignment.orbital_indices={wavefunction.orbitals.size()};
    assignment.axes_available=true;
    assignment.principal_axis={0.0,0.0,1.0};
    assignment.secondary_axis={1.0,0.0,0.0};
    assignment.axis_convention="declared synthetic representation";
    wavefunction.derived_orbital_symmetry_assignments.push_back(assignment);
    orbital.coefficients=coefficients;
    orbital.chemistry.available=true;
    orbital.chemistry.valence_manifold=true;
    orbital.chemistry.valence_weight=1.0;
    orbital.chemistry.confidence=1.0;
    cov::OrbitalAOContribution metal;
    metal.atom_index=0u;
    metal.angular_momentum=metal_angular_momentum;
    metal.weight=0.60;
    orbital.chemistry.ao_contributions.push_back(metal);
    wavefunction.orbitals.push_back(std::move(orbital));
    if (spin==cov::Spin::Beta) ++wavefunction.beta_electrons;
    else ++wavefunction.alpha_electrons;
}

bool validate_general_spin_counterpart_matching() {
    cov::MODiagramOptions options;
    options.selected_index=0u;
    options.max_levels=0u;

    // A valid irrep identifies the same spatial block even when spin
    // polarisation changes which metal angular family has the largest
    // Mulliken population.
    auto family=synthetic_spin_matching_wavefunction(1u,"C3v");
    add_spin_matching_orbital(
        family,cov::Spin::Alpha,-0.20,"A1",{1.0f},1);
    add_spin_matching_orbital(
        family,cov::Spin::Beta,-0.19,"A1",{1.0f},2);
    const auto family_data=cov::build_mo_diagram_data(family,options);
    if (!family_data.spin_counterparts_collapsed ||
        family_data.spin_counterpart_pair_count!=1u ||
        family_data.spin_counterpart_unmatched_visible!=0u) {
        std::cerr<<"same-irrep spin pair was rejected by metal-family drift\n";
        return false;
    }
    auto unresolved=family;
    unresolved.derived_orbital_symmetry_assignments.clear();
    for (auto& orbital:unresolved.orbitals)
        orbital.symmetry_provenance=cov::DataProvenance::Unavailable;
    const auto unresolved_data=cov::build_mo_diagram_data(unresolved,options);
    if (unresolved_data.spin_counterparts_collapsed) {
        std::cerr<<"unscoped matching text bypassed the metal-family guard\n";
        return false;
    }

    // The largest single edge is a trap: choosing it leaves the second alpha
    // group unmatched, while the two crossed edges form a complete matching.
    auto global=synthetic_spin_matching_wavefunction(2u,"C3v");
    add_spin_matching_orbital(
        global,cov::Spin::Alpha,-0.30,"A1",{1.0f,0.0f},2);
    add_spin_matching_orbital(
        global,cov::Spin::Alpha,-0.20,"A1",{0.0f,1.0f},2);
    add_spin_matching_orbital(
        global,cov::Spin::Beta,-0.29,"A1",{0.9f,0.8f},2);
    add_spin_matching_orbital(
        global,cov::Spin::Beta,-0.19,"A1",{0.8f,0.1f},2);
    const auto global_data=cov::build_mo_diagram_data(global,options);
    if (!global_data.spin_counterparts_collapsed ||
        global_data.spin_counterpart_pair_count!=2u ||
        global_data.spin_counterpart_unmatched_visible!=0u) {
        std::cerr<<"spin matching did not maximise counterpart cardinality\n";
        return false;
    }

    // D4d E2 may be printed as two non-adjacent alpha singlets while beta is
    // one formal two-dimensional block.  Their joint S-metric subspace is
    // exact and must become one spatial row in either partition direction.
    auto split=synthetic_spin_matching_wavefunction(3u,"D4d");
    add_spin_matching_orbital(
        split,cov::Spin::Alpha,-0.30,"E2",{1.0f,0.0f,0.0f},2);
    add_spin_matching_orbital(
        split,cov::Spin::Alpha,-0.20,"A1",{0.0f,0.0f,1.0f},0);
    add_spin_matching_orbital(
        split,cov::Spin::Alpha,-0.10,"E2",{0.0f,1.0f,0.0f},2);
    add_spin_matching_orbital(
        split,cov::Spin::Beta,-0.30,"E2",{1.0f,0.0f,0.0f},2);
    add_spin_matching_orbital(
        split,cov::Spin::Beta,-0.30,"E2",{0.0f,1.0f,0.0f},2);
    add_spin_matching_orbital(
        split,cov::Spin::Beta,-0.20,"A1",{0.0f,0.0f,1.0f},0);
    const auto split_data=cov::build_mo_diagram_data(split,options);
    const auto e2=std::find_if(
        split_data.levels.begin(),split_data.levels.end(),[](const auto& level) {
            return level.metadata.symmetry=="E2" &&
                   level.member_indices.size()==2u;
        });
    if (!split_data.spin_counterparts_collapsed ||
        split_data.spin_counterpart_pair_count!=2u ||
        split_data.spin_counterpart_unmatched_visible!=0u ||
        e2==split_data.levels.end() ||
        e2->member_spin_counterparts.size()!=2u) {
        std::cerr<<"split D4d E2 spin partition was not jointly collapsed\n";
        return false;
    }
    return true;
}

bool same_compact_structure(const cov::MODiagramData& left,
                            const cov::MODiagramData& right) {
    if (left.ligand_field_point_group!=right.ligand_field_point_group ||
        left.levels.size()!=right.levels.size() ||
        left.pi_interactions.size()!=right.pi_interactions.size()) {
        return false;
    }
    for (std::size_t i=0;i<left.levels.size();++i) {
        const auto& a=left.levels[i];
        const auto& b=right.levels[i];
        if (a.member_indices!=b.member_indices ||
            a.member_spin_counterparts!=b.member_spin_counterparts ||
            a.metadata.symmetry!=b.metadata.symmetry ||
            a.approximate_nonbonding!=b.approximate_nonbonding ||
            std::abs(a.layout_energy_hartree-b.layout_energy_hartree)>1.0e-12 ||
            std::abs(a.energy_spread_hartree-b.energy_spread_hartree)>1.0e-12 ||
            std::abs(cov::energy_display_coordinate(
                a.layout_energy_hartree,left.energy_transform)-
                cov::energy_display_coordinate(
                    b.layout_energy_hartree,right.energy_transform))>1.0e-12) {
            return false;
        }
    }
    for (std::size_t i=0;i<left.pi_interactions.size();++i) {
        const auto& a=left.pi_interactions[i];
        const auto& b=right.pi_interactions[i];
        if (a.lower_level!=b.lower_level || a.upper_level!=b.upper_level ||
            a.lower_orbitals!=b.lower_orbitals ||
            a.upper_orbitals!=b.upper_orbitals ||
            a.symmetry!=b.symmetry || a.kind!=b.kind ||
            std::abs(a.splitting_hartree-b.splitting_hartree)>1.0e-12 ||
            a.lower_visible!=b.lower_visible ||
            a.upper_visible!=b.upper_visible ||
            a.retained_level!=b.retained_level) {
            return false;
        }
    }
    return true;
}

std::size_t compact_hidden_intermediate_count(
    const cov::DiagramSelectionPlan& selection) {
    constexpr std::string_view marker="intermediate groups hidden=";
    const auto marker_position=selection.summary.find(marker);
    if (marker_position==std::string::npos) return 0u;
    std::size_t value=0u;
    bool found_digit=false;
    for (std::size_t i=marker_position+marker.size();
         i<selection.summary.size();++i) {
        const char character=selection.summary[i];
        if (character<'0' || character>'9') break;
        found_digit=true;
        value=10u*value+static_cast<std::size_t>(character-'0');
    }
    return found_digit?value:0u;
}

std::size_t coordination_number_encoded_in_filename(
    const std::string_view filename) {
    const auto marker=filename.find("_CN");
    if (marker==std::string_view::npos) return 0u;
    std::size_t value=0u;
    bool found_digit=false;
    for (std::size_t i=marker+3u;i<filename.size();++i) {
        const char character=filename[i];
        if (character<'0' || character>'9') break;
        found_digit=true;
        value=10u*value+static_cast<std::size_t>(character-'0');
    }
    return found_digit?value:0u;
}

bool inspect_real_fchk(const std::filesystem::path& path) {
    const auto wf=cov::parse_wavefunction(path);
    cov::MODiagramOptions options;
    options.selected_index=0u;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].occupation>1.0e-4f) options.selected_index=i;
    }
    options.max_levels=0u;
    const auto data=cov::build_mo_diagram_data(wf,options);
    // This path reads FCHK only and deliberately supplies no independently
    // validated local Fock channels. Directed-pair acceptance belongs to the
    // saved-package routed tests; geometry, sigma, spin and weak checks below remain.
    if(std::any_of(data.pi_interactions.begin(),data.pi_interactions.end(),[](const auto& pair){
        return pair.kind==cov::PiInteractionKind::Donor||pair.kind==cov::PiInteractionKind::Acceptor;
    })) {
        std::cerr<<path.filename().string()<<": FCHK-only data asserted a strong directed pi pair\n";
        return false;
    }
    const auto browser_metadata=cov::build_orbital_metadata(
        wf,options.selected_index,options.degeneracy,options.filter);
    if (wf.delocalised_pi_assignments.empty()) {
        std::vector<std::pair<std::uint32_t,std::uint32_t>> topology_edges;
        for (const auto& bond:cov::analyse_bonds(wf)) {
            topology_edges.emplace_back(
                static_cast<std::uint32_t>(bond.atom_a),
                static_cast<std::uint32_t>(bond.atom_b));
        }
        const auto networks=cov::infer_oriented_pi_networks(
            wf,topology_edges);
        std::cout<<"  diagnostic-pi-networks="<<networks.size();
        for (const auto& network:networks) {
            std::cout<<" [atoms=";
            for (const auto atom:network.atoms) std::cout<<atom+1u<<',';
            std::cout<<" edges="<<network.edges.size()
                     <<" cyclic="<<network.cyclic
                     <<" coherence="<<network.orientation_coherence<<']';
        }
        std::cout<<'\n';
    }
    std::cout<<"DIAGNOSTIC "<<path.filename().string()
             <<" point-group="<<wf.point_group_detected
             <<" alpha="<<wf.alpha_electrons
             <<" beta="<<wf.beta_electrons
             <<" orbitals="<<wf.orbitals.size()
             <<" rows="<<data.levels.size()
             <<" selected="<<data.selection.included_indices.size()
             <<" local-field="<<data.ligand_field_point_group
             <<" geometry="<<data.ligand_field_geometry_id
             <<" CN="<<data.ligand_field_coordination_number
             <<" first-shell="<<data.ligand_field_ligand_atoms.size()
             <<" spin-collapsed="<<data.spin_counterparts_collapsed
             <<" spin-partial="<<data.spin_counterparts_partial
             <<" spin-pairs="<<data.spin_counterpart_pair_count
             <<" spin-unmatched="<<data.spin_counterpart_unmatched_visible
             <<" pi-pairs="<<data.pi_interactions.size()<<'\n';
    if (!wf.orbitals.empty() && !wf.orbitals.front().chemistry.available) {
        std::cout<<"  diagnostic-chemistry-unavailable: "
                 <<wf.orbitals.front().chemistry.note<<'\n';
    }
    std::cout<<"  diagnostic-first-shell metal="<<data.ligand_field_metal_atom;
    for (const auto atom:data.ligand_field_ligand_atoms) {
        if (atom<wf.atoms.size()) {
            std::cout<<' '<<atom<<':'<<wf.atoms[atom].symbol;
        }
    }
    std::cout<<'\n';
    for (const auto& bond:wf.bond_orders) {
        if (bond.atom_a>=wf.atoms.size() || bond.atom_b>=wf.atoms.size()) continue;
        const int za=wf.atoms[bond.atom_a].atomic_number;
        const int zb=wf.atoms[bond.atom_b].atomic_number;
        const auto tm=[](const int z) {
            return (z>=21 && z<=30) || (z>=39 && z<=48) ||
                   (z>=72 && z<=80) || (z>=104 && z<=112);
        };
        if (!tm(za) && !tm(zb)) continue;
        std::cout<<"  diagnostic-bond "<<bond.atom_a<<'-'<<bond.atom_b
                 <<' '<<wf.atoms[bond.atom_a].symbol<<'-'
                 <<wf.atoms[bond.atom_b].symbol
                 <<" Mayer="<<bond.mayer_order<<'\n';
    }
    for (std::size_t i=0;i<data.levels.size();++i) {
        const auto& level=data.levels[i];
        std::cout<<"  diagnostic-row "<<i
                 <<" MO"<<level.metadata.raw_mo_number
                 <<" spin="<<(level.member_indices.empty()
                        ?-1:static_cast<int>(
                            wf.orbitals[level.member_indices.front()].spin))
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" occ="<<level.total_occupation
                 <<" s="<<level.metal_s_weight
                 <<" p="<<level.metal_p_weight
                 <<" d="<<level.metal_d_weight
                 <<" direct-Lp="<<level.direct_ligand_p_weight
                 <<" Lp="<<level.ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" pi="<<level.pi_fraction
                 <<" ML="<<level.metal_ligand_overlap
                 <<" raw="<<level.raw_data_fallback<<'\n';
    }
    if (data.levels.empty() || data.selection.included_indices.empty()) {
        std::cerr<<path.filename().string()<<": empty ligand-field diagram\n";
        return false;
    }
    for (std::size_t level_index=0;level_index<data.levels.size();++level_index) {
        const auto& level=data.levels[level_index];
        const double allowed_spread=data.ligand_field_point_group.empty()
            ?options.degeneracy.tolerance_hartree:0.01501;
        if (level.member_indices.empty() ||
            level.member_indices.size()!=level.metadata.degeneracy_size ||
            level.energy_spread_hartree>allowed_spread) {
            std::cerr<<path.filename().string()<<": broken degenerate row "
                     <<level_index<<" sym="<<level.metadata.symmetry
                     <<" members="<<level.member_indices.size()
                     <<" declared="<<level.metadata.degeneracy_size
                     <<" spread="<<level.energy_spread_hartree
                     <<" allowed="<<allowed_spread<<'\n';
            return false;
        }
    }
    const std::size_t row_budget=options.max_levels>0u
        ?std::max<std::size_t>(4u,options.max_levels)
        :std::clamp<std::size_t>(2u*options.neighbourhood+1u,10u,48u);
    if (data.levels.size()>data.selection.included_indices.size()) {
        std::cerr<<path.filename().string()
                 <<": valence reduction produced an invalid row set\n";
        return false;
    }
    // Producer-derived ligand-field manifolds can be complete without using
    // the generic raw-MO supplement.  Raw fallback is therefore an optional
    // recovery path, not a correctness requirement for every fixture.
    auto compact_options=options;
    compact_options.hide_ligand_centred_intermediates=true;
    compact_options.mode=cov::preferred_compact_mo_diagram_mode(wf,true);
    const auto compact=cov::build_mo_diagram_data(wf,compact_options);
    const std::size_t actual_overflow=compact.levels.size()>row_budget
        ?compact.levels.size()-row_budget:0u;
    if (compact.levels.size()>compact.selection.included_indices.size() ||
        compact.selection.protected_overflow_count!=actual_overflow) {
        std::cerr<<path.filename().string()
                 <<": compact row budget/protected-overflow contract failed\n";
        return false;
    }
    if ((compact_options.mode==cov::MODiagramMode::DelocalisedPiFamilyOnly ||
         compact_options.mode==cov::MODiagramMode::MulticentreActiveSpaceOnly) &&
        compact.levels.empty()) {
        std::cerr<<path.filename().string()
                 <<": preferred active-space mode produced no rows\n";
        return false;
    }
    const std::string filename=path.filename().string();
    auto pi_contract=[&](const std::size_t atoms,
                         const std::size_t channels,
                         const std::size_t orbitals,
                         const cov::DelocalisedPiTopology topology,
                         const double electrons) {
        return std::any_of(
            wf.delocalised_pi_assignments.begin(),
            wf.delocalised_pi_assignments.end(),
            [&](const auto& assignment) {
                return assignment.atoms.size()==atoms &&
                       assignment.orientation_channels.size()==channels &&
                       assignment.orbitals.size()==orbitals &&
                       assignment.topology==topology &&
                       std::abs(assignment.electron_count-electrons)<0.25;
            });
    };
    auto require_pi_mode=[&](const bool expected,const char* contract) {
        const bool actual=compact.mode==
            cov::MODiagramMode::DelocalisedPiFamilyOnly;
        if (actual==expected) return true;
        std::cerr<<filename<<": "<<contract<<" compact-mode contract failed\n";
        return false;
    };
    for (const auto& assignment:wf.delocalised_pi_assignments) {
        const auto& graph=assignment.topology_graph;
        if(!graph.channel_association_search_complete || graph.atom_count!=wf.atoms.size()) {
            std::cerr<<filename<<": pi family structural graph scope is incomplete\n";
            return false;
        }
        if(assignment.topology==cov::DelocalisedPiTopology::Spiro && graph.channel_ring_witnesses.empty()) {
            std::cerr<<filename<<": spiro label has no structural ring paths\n";
            return false;
        }
        for(const auto& evidence:graph.channel_ring_witnesses) {
            if(!cov::verify_cycle_union(static_cast<std::uint32_t>(wf.atoms.size()),graph.edges,evidence.ring_union)) {
                std::cerr<<filename<<": invalid structural ring paths\n";
                return false;
            }
            for(const auto channel:evidence.channel_indices)if(channel>=assignment.orientation_channels.size()) {
                std::cerr<<filename<<": structural ring references an absent pi channel\n";
                return false;
            }
        }
    }
    if (filename.find("057_hexatriene")!=std::string::npos &&
        (!require_pi_mode(true,"conjugated path") ||
         !pi_contract(6u,1u,6u,cov::DelocalisedPiTopology::Path,6.0))) {
        std::cerr<<filename<<": complete path pi family missing\n";
        return false;
    }
    if ((filename.find("086_carbonate")!=std::string::npos ||
         filename.find("087_nitrate")!=std::string::npos ||
         filename.find("117_guanidinium")!=std::string::npos) &&
        (!require_pi_mode(true,"branched resonance") ||
         !pi_contract(4u,1u,4u,
             cov::DelocalisedPiTopology::BranchedResonance,6.0))) {
        std::cerr<<filename<<": occupied branched-resonance family missing\n";
        return false;
    }
    if ((filename.find("029_bf3")!=std::string::npos ||
         filename.find("039_xef4")!=std::string::npos) &&
        !require_pi_mode(false,"donor/full ligand-p star")) {
        return false;
    }
    if (filename.find("069_spiropentadiene")!=std::string::npos &&
        (!require_pi_mode(true,"orthogonal spiro") ||
         !pi_contract(4u,2u,4u,
             cov::DelocalisedPiTopology::Spiro,4.0))) {
        std::cerr<<filename<<": orthogonal spiro direct-sum family missing\n";
        return false;
    }
    if (filename.find("092_ferrocene")!=std::string::npos &&
        (!require_pi_mode(true,"haptic-metal") ||
         !pi_contract(10u,2u,10u,
             cov::DelocalisedPiTopology::HapticMetal,12.0))) {
        std::cerr<<filename<<": joint eta5-Cp active subspace missing\n";
        return false;
    }
    if (filename.find("133_tetrakis_pentafluorophenyl_borate")!=
            std::string::npos &&
        (!require_pi_mode(true,"multi-ring direct sum") ||
         !pi_contract(24u,4u,24u,
             cov::DelocalisedPiTopology::SymmetryDirectSum,24.0))) {
        std::cerr<<filename<<": four-ring pi direct sum missing\n";
        return false;
    }
    if (filename.find("001_so2")!=std::string::npos &&
        (!require_pi_mode(true,"SO2 valence-p") ||
         !pi_contract(3u,1u,3u,cov::DelocalisedPiTopology::Path,4.0))) {
        std::cerr<<filename<<": S-core-orthogonal SO2 family missing\n";
        return false;
    }
    if (filename.find("006_thiourea")!=std::string::npos &&
        (!require_pi_mode(true,"thiourea valence-p") ||
         !pi_contract(4u,1u,4u,
             cov::DelocalisedPiTopology::BranchedResonance,6.0))) {
        std::cerr<<filename<<": S-core-orthogonal thiourea family missing\n";
        return false;
    }
    if ((filename.find("009_urea_planar_start")!=std::string::npos ||
         filename.find("010_urea_nonplanar_start")!=std::string::npos) &&
        (!require_pi_mode(true,"optimised urea") ||
         !pi_contract(4u,1u,4u,
             cov::DelocalisedPiTopology::BranchedResonance,6.0))) {
        std::cerr<<filename<<": optimised urea pi family missing\n";
        return false;
    }
    if ((filename.find("090_feco5")!=std::string::npos ||
         filename.find("091_crco6")!=std::string::npos) &&
        (!wf.delocalised_pi_assignments.empty() ||
         compact.mode!=cov::MODiagramMode::ValenceCentral)) {
        std::cerr<<filename<<": separate CO ligands were bundled as one pi family\n";
        return false;
    }
    const std::size_t encoded_coordination_number=
        coordination_number_encoded_in_filename(filename);
    if (encoded_coordination_number>0u &&
        (data.ligand_field_geometry_id.empty() ||
         data.ligand_field_point_group.empty() ||
         data.ligand_field_coordination_number!=encoded_coordination_number ||
         data.ligand_field_ligand_atoms.size()!=encoded_coordination_number)) {
        std::cerr<<filename
                 <<": encoded coordination shell was not resolved (expected CN="
                 <<encoded_coordination_number<<", geometry="
                 <<data.ligand_field_geometry_id<<", point-group="
                 <<data.ligand_field_point_group<<", CN="
                 <<data.ligand_field_coordination_number<<", shell="
                 <<data.ligand_field_ligand_atoms.size()<<")\n";
        return false;
    }
    const cov::CoordinationGeometryDescriptor* encoded_geometry=nullptr;
    for (const auto& descriptor:cov::coordination_geometry_catalog()) {
        const std::string token="_"+std::string(descriptor.machine_id)+"_";
        if (filename.find(token)!=std::string::npos) {
            encoded_geometry=&descriptor;
            break;
        }
    }
    if (encoded_geometry!=nullptr &&
        (data.ligand_field_geometry_id!=encoded_geometry->machine_id ||
         data.ligand_field_point_group!=encoded_geometry->point_group ||
         data.ligand_field_coordination_number!=
             encoded_geometry->coordination_number)) {
        std::cerr<<filename<<": encoded geometry "
                 <<encoded_geometry->machine_id<<'/'
                 <<encoded_geometry->point_group<<" was resolved as "
                 <<data.ligand_field_geometry_id<<'/'
                 <<data.ligand_field_point_group<<" (CN="
                 <<data.ligand_field_coordination_number<<")\n";
        return false;
    }
    // The diagram and the 3-D molecule consume the same wavefunction but have
    // historically used independent adjacency filters. Every donor accepted
    // into a resolved first coordination shell must also remain connected to
    // its metal in the rendered molecular graph. This generic invariant
    // catches missing element radii and future threshold drift for every CN.
    if (!data.ligand_field_ligand_atoms.empty()) {
        const auto rendered_bonds=cov::analyse_bonds(wf);
        for (const auto ligand_atom:data.ligand_field_ligand_atoms) {
            const bool connected=std::any_of(
                rendered_bonds.begin(),rendered_bonds.end(),
                [&](const auto& bond) {
                    return (bond.atom_a==data.ligand_field_metal_atom &&
                            bond.atom_b==ligand_atom) ||
                           (bond.atom_b==data.ligand_field_metal_atom &&
                            bond.atom_a==ligand_atom);
                });
            if (!connected) {
                std::cerr<<filename<<": first-shell atom "<<ligand_atom
                         <<" is absent from the rendered metal-ligand graph\n";
                return false;
            }
        }
    }
    for (std::size_t i=0;i<compact.levels.size();++i) {
        const auto& level=compact.levels[i];
        std::cout<<"  compact-row "<<i
                 <<" MO"<<level.metadata.raw_mo_number
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" coord="<<cov::energy_display_coordinate(
                        level.layout_energy_hartree,compact.energy_transform)
                 <<" occ="<<level.total_occupation
                 <<" d="<<level.metal_d_weight
                 <<" direct-Lp="<<level.direct_ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" ML="<<level.metal_ligand_overlap<<'\n';
    }
    std::cout<<"COMPACT "<<path.filename().string()
             <<" mode="<<cov::mo_diagram_mode_name(compact.mode)
             <<" geometry="<<compact.ligand_field_geometry_id
             <<" CN="<<compact.ligand_field_coordination_number
             <<" point-group="<<compact.ligand_field_point_group
             <<" symmetries=";
    for (const auto& level:compact.levels) {
        std::cout<<' '<<level.metadata.symmetry
                 <<'('<<level.member_indices.size()<<')';
    }
    std::cout<<'\n';
    const std::size_t hidden_intermediates=
        compact_hidden_intermediate_count(compact.selection);
    if (compact.levels.size()>data.levels.size() ||
        (hidden_intermediates>0u &&
         compact.levels.size()>=data.levels.size()) ||
        compact.pi_interactions.size()!=data.pi_interactions.size()) {
        std::cerr<<path.filename().string()
                 <<": intermediate-orbital toggle did not reduce safely ("
                 <<data.levels.size()<<" -> "<<compact.levels.size()
                 <<" rows; "<<data.pi_interactions.size()<<" -> "
                 <<compact.pi_interactions.size()<<" pairs; hidden="
                 <<hidden_intermediates<<")\n";
        return false;
    }
    const bool selection_invariance_fixture=
        filename.find("TiF6")!=std::string::npos ||
        filename.find("ZnCl4")!=std::string::npos ||
        filename.find("19_FeF6")!=std::string::npos;
    if (selection_invariance_fixture) {
        auto canonical_options=compact_options;
        canonical_options.selected_index=wf.orbitals.size();
        const auto canonical=cov::build_mo_diagram_data(
            wf,canonical_options);
        if (!same_compact_structure(compact,canonical)) {
            std::cerr<<filename
                     <<": compact row set changed with the initial selection\n";
            return false;
        }
        const bool octahedral=canonical.ligand_field_point_group=="Oh";
        const std::string hidden_symmetry=octahedral?"T1g":"T1";
        const std::set<std::string> probe_symmetries=octahedral
            ?std::set<std::string>{"T1g","T1u","T2g"}
            :std::set<std::string>{"T1","T2","E"};
        std::vector<std::size_t> probes;
        bool source_contains_hidden=false;
        for (std::size_t i=0;i<data.metadata.size();++i) {
            const auto& symmetry=data.metadata[i].symmetry;
            if (symmetry==hidden_symmetry) source_contains_hidden=true;
            if (probe_symmetries.count(symmetry)) probes.push_back(i);
        }
        std::sort(probes.begin(),probes.end());
        probes.erase(std::unique(probes.begin(),probes.end()),probes.end());
        for (const auto selected:probes) {
            auto selected_options=canonical_options;
            selected_options.selected_index=selected;
            const auto selected_data=cov::build_mo_diagram_data(
                wf,selected_options);
            if (!same_compact_structure(canonical,selected_data)) {
                std::cerr<<filename
                         <<": compact row set changed after selecting MO"
                         <<selected+1u<<'\n';
                return false;
            }
        }
        const auto compact_has=[&](const std::string& symmetry) {
            return std::any_of(
                canonical.levels.begin(),canonical.levels.end(),
                [&](const auto& level) {
                    return level.metadata.symmetry==symmetry;
                });
        };
        const bool required_framework=octahedral
            ?compact_has("T1u") && compact_has("T2g")
            :compact_has("T2");
        if (!source_contains_hidden || compact_has(hidden_symmetry) ||
            !required_framework) {
            std::cerr<<filename
                     <<": compact T1/T1g exclusion damaged the metal framework\n";
            return false;
        }
    }
    const auto valid_symmetry=[](const std::string& symmetry) {
        return !symmetry.empty() && symmetry!="?" &&
               symmetry!="N/A" && symmetry!="n/a";
    };
    const auto level_metrics_are_bounded=[](const auto& levels) {
        return std::all_of(levels.begin(),levels.end(),[](const auto& level) {
            return std::isfinite(level.metal_ligand_overlap) &&
                   std::abs(level.metal_ligand_overlap)<=1.000001;
        });
    };
    if (!level_metrics_are_bounded(data.levels) ||
        !level_metrics_are_bounded(compact.levels)) {
        std::cerr<<path.filename().string()
                 <<": first-shell metal-ligand metric is not finite/bounded\n";
        return false;
    }
    if (!data.ligand_field_point_group.empty()) {
        const bool bad_shell=data.ligand_field_metal_atom>=wf.atoms.size() ||
            data.ligand_field_ligand_atoms.empty() ||
            std::any_of(data.ligand_field_ligand_atoms.begin(),
                        data.ligand_field_ligand_atoms.end(),[&](const auto atom) {
                return atom>=wf.atoms.size() ||
                       atom==data.ligand_field_metal_atom ||
                       std::count(data.ligand_field_ligand_atoms.begin(),
                                  data.ligand_field_ligand_atoms.end(),atom)!=1;
            });
        const bool bad_compact_label=std::any_of(
            compact.levels.begin(),compact.levels.end(),[&](const auto& level) {
                if (!valid_symmetry(level.metadata.symmetry)) return true;
                return std::any_of(
                    level.member_indices.begin(),level.member_indices.end(),
                    [&](const auto member) {
                        return member>=wf.orbitals.size() ||
                               member>=data.metadata.size() ||
                               member>=browser_metadata.size() ||
                               data.metadata[member].symmetry!=
                                   level.metadata.symmetry ||
                               browser_metadata[member].symmetry!=
                                   level.metadata.symmetry ||
                               std::count(level.member_indices.begin(),
                                          level.member_indices.end(),member)!=1;
                    });
            });
        if (bad_shell || bad_compact_label) {
            for (const auto& level:compact.levels) {
                for (const auto member:level.member_indices) {
                    if (member>=data.metadata.size() ||
                        member>=browser_metadata.size() ||
                        data.metadata[member].symmetry!=level.metadata.symmetry ||
                        browser_metadata[member].symmetry!=level.metadata.symmetry) {
                        std::cerr<<"  member "<<member
                                 <<" level="<<level.metadata.symmetry
                                 <<" diagram="<<(member<data.metadata.size()
                                     ?data.metadata[member].symmetry:"<out>")
                                 <<" browser="<<(member<browser_metadata.size()
                                     ?browser_metadata[member].symmetry:"<out>")
                                 <<'\n';
                    }
                }
            }
            std::cerr<<path.filename().string()
                     <<": local-field shell/compact member metadata regression\n";
            return false;
        }
    }
    for (const auto& interaction:data.pi_interactions) {
        const auto validate_members=[&](const auto& members) {
            return std::all_of(members.begin(),members.end(),[&](const auto member) {
                if (member>=data.metadata.size()) return false;
                return valid_symmetry(data.metadata[member].symmetry);
            });
        };
        if (!validate_members(interaction.lower_orbitals) ||
            !validate_members(interaction.upper_orbitals)) {
            std::cerr<<path.filename().string()
                     <<": split-orbital member symmetry fell back to N/A\n";
            return false;
        }
    }

    const auto expected_pair=[&](const cov::PiInteractionKind kind,
                                 const std::string& symmetry) {
        return std::find_if(
            data.pi_interactions.begin(),data.pi_interactions.end(),
            [&](const auto& pair) {
                return pair.kind==kind && pair.symmetry==symmetry;
            });
    };
    const auto all_members_have=[&](const auto& members,
                                     const std::string& symmetry) {
        return std::all_of(members.begin(),members.end(),[&](const auto member) {
            return member<data.metadata.size() &&
                   data.metadata[member].symmetry==symmetry;
        });
    };
    const auto compact_group=[&](const std::string& symmetry,
                                 const std::size_t degeneracy) {
        return std::find_if(
            compact.levels.begin(),compact.levels.end(),[&](const auto& level) {
                return level.metadata.symmetry==symmetry &&
                       level.member_indices.size()==degeneracy;
            });
    };
    const auto expected_shell=[&](const std::string& point_group,
                                  const std::size_t ligand_count,
                                  const int metal_z,
                                  const int ligand_z) {
        return data.ligand_field_point_group==point_group &&
            data.ligand_field_metal_atom<wf.atoms.size() &&
            wf.atoms[data.ligand_field_metal_atom].atomic_number==metal_z &&
            data.ligand_field_ligand_atoms.size()==ligand_count &&
            std::all_of(data.ligand_field_ligand_atoms.begin(),
                        data.ligand_field_ligand_atoms.end(),[&](const auto atom) {
                return atom<wf.atoms.size() &&
                       wf.atoms[atom].atomic_number==ligand_z;
            });
    };
    const auto fully_collapsed_open_shell=[&]() {
        return data.spin_counterparts_collapsed &&
               !data.spin_counterparts_partial &&
               data.spin_counterpart_pair_count>0u &&
               data.spin_counterpart_unmatched_visible==0u;
    };
    if (filename.find("TiF6")!=std::string::npos) {
        if (!expected_shell("Oh",6u,22,9)) {
            std::cerr<<filename<<": FCHK-only donor geometry regression\n";
            return false;
        }
    } else if (filename.find("ZnCl4")!=std::string::npos) {
        // FCHK composition and a close E/T2 gap do not provide the frozen
        // local operator needed to certify a negligible pi relation.
        if (!expected_shell("Td",4u,30,17) ||
            !data.pi_interactions.empty() ||
            std::any_of(data.pi_partner_candidates.begin(),data.pi_partner_candidates.end(),
                [](const auto& candidate){return candidate.accepted||candidate.weak;})) {
            std::cerr<<filename
                     <<": FCHK-only close E/T2 gap invented a calibrated weak pi relation\n";
            return false;
        }
    } else if (filename.find("CrCO6")!=std::string::npos) {
        const bool sigma_bonding=std::any_of(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.sigma_fraction>=0.55 &&
                       level.metal_s_weight+level.metal_p_weight>=0.03 &&
                       level.metal_ligand_overlap>0.0;
            });
        const bool sigma_antibonding=std::any_of(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.sigma_fraction>=0.55 &&
                       level.metal_s_weight+level.metal_p_weight>=0.03 &&
                       level.metal_ligand_overlap<0.0;
            });
        if (!expected_shell("Oh",6u,24,6) ||
            compact.levels.size()>10u || !sigma_bonding || !sigma_antibonding) {
            std::cerr<<filename<<": FCHK-only carbonyl geometry/sigma regression\n";
            return false;
        }
    } else if (filename.find("CrNH3_6")!=std::string::npos) {
        const auto find_compact=[&](const std::string& symmetry,
                                    const std::size_t degeneracy,
                                    const bool antibonding) {
            return std::find_if(
                compact.levels.begin(),compact.levels.end(),
                [&](const auto& level) {
                    return level.metadata.symmetry==symmetry &&
                           level.member_indices.size()==degeneracy &&
                           (level.metal_ligand_overlap<0.0)==antibonding;
                });
        };
        const auto t2g=std::find_if(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.metadata.symmetry=="T2g" &&
                       level.member_indices.size()==3u &&
                       level.total_occupation>2.9 &&
                       level.total_occupation<3.1;
            });
        const bool t2g_arrows=t2g!=compact.levels.end() &&
            t2g->member_electrons.size()==3u &&
            std::all_of(t2g->member_electrons.begin(),
                        t2g->member_electrons.end(),[](const auto& electrons) {
                return electrons.alpha==1 && electrons.beta==0;
            }) && t2g->member_spin_counterparts.size()==3u &&
            std::all_of(t2g->member_spin_counterparts.begin(),
                        t2g->member_spin_counterparts.end(),[&](const auto index) {
                return index<wf.orbitals.size() &&
                       wf.orbitals[index].spin==cov::Spin::Beta &&
                       index<data.metadata.size() &&
                       data.metadata[index].symmetry=="T2g";
            });
        const auto labelled_count=[&](const std::string& symmetry) {
            return static_cast<std::size_t>(std::count_if(
                data.metadata.begin(),data.metadata.end(),
                [&](const auto& item) {return item.symmetry==symmetry;}));
        };
        const auto browser_labelled_count=[&](const std::string& symmetry) {
            return static_cast<std::size_t>(std::count_if(
                browser_metadata.begin(),browser_metadata.end(),
                [&](const auto& item) {return item.symmetry==symmetry;}));
        };
        if (wf.point_group_detected=="Oh" ||
            !expected_shell("Oh",6u,24,7) ||
            data.ligand_field_confidence<0.90 ||
            !data.spin_counterparts_collapsed ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            find_compact("A1g",1u,false)==compact.levels.end() ||
            find_compact("A1g",1u,true)==compact.levels.end() ||
            find_compact("T1u",3u,false)==compact.levels.end() ||
            find_compact("T1u",3u,true)==compact.levels.end() ||
            find_compact("Eg",2u,false)==compact.levels.end() ||
            find_compact("Eg",2u,true)==compact.levels.end() ||
            !t2g_arrows || labelled_count("T2g")<6u ||
            labelled_count("Eg")<8u ||
            browser_labelled_count("T2g")<6u ||
            browser_labelled_count("Eg")<8u) {
            std::cerr<<filename
                     <<": local-Oh open-shell browser/diagram metadata regression\n";
            return false;
        }
    } else if (filename.find("19_FeF6")!=std::string::npos) {
        if (!expected_shell("Oh",6u,26,9) ||
            !fully_collapsed_open_shell() ||
            !data.pi_interactions.empty() || compact.levels.size()>10u ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": weak-field open-shell Oh reduction regression\n";
            return false;
        }
    } else if (filename.find("20_CoCl4")!=std::string::npos) {
        if (!expected_shell("Td",4u,27,17) ||
            !fully_collapsed_open_shell() ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("E",2u)==compact.levels.end() ||
            compact_group("T2",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": weak-field open-shell Td reduction regression\n";
            return false;
        }
    } else if (filename.find("21_CoNH3_6")!=std::string::npos) {
        if (wf.point_group_detected=="Oh" ||
            !expected_shell("Oh",6u,27,7) ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("A1g",1u)==compact.levels.end() ||
            compact_group("T1u",3u)==compact.levels.end() ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": sigma-only local-Oh framework regression\n";
            return false;
        }
    } else if (filename.find("22_ZnNH3_4")!=std::string::npos) {
        if (wf.point_group_detected=="Td" ||
            !expected_shell("Td",4u,30,7) ||
            !data.pi_interactions.empty() || compact.levels.size()>9u ||
            compact_group("A1",1u)==compact.levels.end() ||
            compact_group("E",2u)==compact.levels.end() ||
            compact_group("T2",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": sigma-only local-Td framework regression\n";
            return false;
        }
    } else if (filename.find("23_NiCO4")!=std::string::npos) {
        if (!expected_shell("Td",4u,28,6) ||
            compact.levels.size()>8u) {
            std::cerr<<filename
                     <<": strong-field Td pi-acceptor regression\n";
            return false;
        }
    } else if (filename.find("24_CrCN6")!=std::string::npos) {
        if (!expected_shell("Oh",6u,24,6) ||
            !fully_collapsed_open_shell() ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end() ||
            compact.levels.size()>11u) {
            std::cerr<<filename
                     <<": strong-field open-shell Oh pi-acceptor regression\n";
            return false;
        }
    } else if (filename.find("25_AgNH3_2")!=std::string::npos) {
        if (!expected_shell("Dinfh",2u,47,7) ||
            data.ligand_field_geometry_id!="L-2" ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("Sigma_g+",1u)==compact.levels.end() ||
            compact_group("Sigma_u+",1u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": linear CN2 sigma-framework regression\n";
            return false;
        }
    } else if (filename.find("26_FeCl3")!=std::string::npos) {
        if (!expected_shell("D3h",3u,26,17) ||
            data.ligand_field_geometry_id!="TP-3" ||
            !fully_collapsed_open_shell() || compact.levels.size()>8u ||
            compact_group("A1'",1u)==compact.levels.end() ||
            compact_group("E'",2u)==compact.levels.end() ||
            compact_group("E''",2u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": trigonal-planar CN3 open-shell regression\n";
            return false;
        }
    } else if (filename.find("27_FeCO5")!=std::string::npos) {
        if (!expected_shell("D3h",5u,26,6) ||
            data.ligand_field_geometry_id!="TBPY-5" ||
            compact.levels.size()>10u) {
            std::cerr<<filename
                     <<": trigonal-bipyramidal CN5 pi-acceptor regression\n";
            return false;
        }
    } else if (filename.find("28_ZrF7")!=std::string::npos) {
        if (!expected_shell("D5h",7u,40,9) ||
            data.ligand_field_geometry_id!="PBPY-7" ||
            compact.levels.size()>11u) {
            std::cerr<<filename
                     <<": pentagonal-bipyramidal CN7 pi-donor regression\n";
            return false;
        }
    } else if (filename.find("29_MoCN8")!=std::string::npos) {
        if (!expected_shell("D4d",8u,42,6) ||
            data.ligand_field_geometry_id!="SAPR-8" ||
            !fully_collapsed_open_shell() ||
            compact.levels.size()>11u) {
            std::cerr<<filename
                     <<": square-antiprismatic CN8 open-shell regression\n";
            return false;
        }
    } else if (filename.find("30_ReH9")!=std::string::npos) {
        if (!expected_shell("D3h",9u,75,1) ||
            data.ligand_field_geometry_id!="TCTPR-9" ||
            !data.pi_interactions.empty() || compact.levels.size()>11u ||
            compact_group("A1'",1u)==compact.levels.end() ||
            compact_group("E'",2u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": tricapped-trigonal-prismatic CN9 regression\n";
            return false;
        }
    } else if (filename.find("31_HfH2O10")!=std::string::npos) {
        if (!expected_shell("D5d",10u,72,8) ||
            data.ligand_field_geometry_id!="PAPR-10" ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("A1g",1u)==compact.levels.end() ||
            compact_group("E1u",2u)==compact.levels.end() ||
            compact_group("E2g",2u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": pentagonal-antiprismatic CN10 regression\n";
            return false;
        }
    } else if (filename.find("32_CuCl4")!=std::string::npos) {
        if (!expected_shell("C3v",4u,29,17) ||
            data.ligand_field_geometry_id!="vTBPY-4" ||
            !fully_collapsed_open_shell() ||
            !data.pi_interactions.empty() || compact.levels.size()>6u ||
            compact_group("A1",1u)==compact.levels.end() ||
            compact_group("E",2u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": vacant-trigonal-bipyramidal CN4 regression\n";
            return false;
        }
    }

    std::cout<<"REAL "<<path.filename().string()
             <<" point-group="<<wf.point_group_detected
             <<" rows="<<data.levels.size()
             <<" compact-rows="<<compact.levels.size()
             <<" raw-MOs="<<data.selection.included_indices.size()
             <<" pi-pairs="<<data.pi_interactions.size()<<'\n';
    for (std::size_t i=0;i<data.levels.size();++i) {
        const auto& level=data.levels[i];
        std::cout<<"  row "<<i<<" MO"<<level.metadata.raw_mo_number
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" d="<<level.metal_d_weight
                 <<" Lp="<<level.ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" pi="<<level.pi_fraction
                 <<" ML="<<level.metal_ligand_overlap
                 <<" raw="<<level.raw_data_fallback
                 <<" approxNB="<<level.approximate_nonbonding<<'\n';
    }
    for (const auto& pair:data.pi_interactions) {
        std::cout<<"  PI "<<cov::pi_interaction_kind_name(pair.kind)
                 <<" sym="<<pair.symmetry
                 <<" split="<<pair.splitting_hartree
                 <<" confidence="<<pair.confidence
                 <<" lower-visible="<<pair.lower_visible
                 <<" upper-visible="<<pair.upper_visible<<'\n';
    }
    return true;
}

} // namespace

int main(int argc,char** argv) {
    std::vector<std::string> failed_real_fchk;
    for (int i=1;i<argc;++i) {
        const std::filesystem::path path=argv[i];
        const auto started=std::chrono::steady_clock::now();
        bool passed=false;
        std::string exception_message;
        try {
            passed=inspect_real_fchk(path);
        } catch (const std::exception& exception) {
            exception_message=exception.what();
        } catch (...) {
            exception_message="unknown exception";
        }
        const auto elapsed=std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now()-started);
        std::cout<<"CASE "<<path.filename().string()<<' '
                 <<(passed?"PASS":"FAIL")
                 <<" elapsed_ms="<<elapsed.count()<<'\n';
        if (!passed) {
            failed_real_fchk.push_back(path.string());
            if (!exception_message.empty()) {
                std::cerr<<path.filename().string()<<": "
                         <<exception_message<<'\n';
            }
        }
    }
    if (!failed_real_fchk.empty()) {
        std::cerr<<"REAL FCHK FAILURES "<<failed_real_fchk.size()<<'/'
                 <<(argc-1)<<'\n';
        for (const auto& path:failed_real_fchk) {
            std::cerr<<"  "<<path<<'\n';
        }
    }

    if (!validate_radial_first_shell_retry()) return 14;
    if (!validate_local_symmetry_missing_markers()) return 12;
    if (!validate_ligand_only_noise_preserves_producer_irrep()) return 18;
    if (!validate_resolved_five_d_runs()) return 15;
    if (!validate_minimal_compact_ligand_field_framework()) return 17;
    if (!validate_pi_topology_and_two_sided_composition()) return 16;
    if (!validate_equal_count_unrestricted_collapse()) return 13;
    if (!validate_general_spin_counterpart_matching()) return 19;

    // A modern main-group row has no central-metal/ligand scope.  Zero in the
    // metal-ligand accumulator is therefore unavailable evidence, not a
    // 0%-confidence nonbonding classification.
    {
        cov::Wavefunction main_group;
        main_group.atoms.push_back({"C",6,0.0,0.0,0.0});
        cov::MolecularOrbital valence;
        valence.energy_hartree=-0.45;
        valence.occupation=2.0f;
        valence.chemistry.available=true;
        valence.chemistry.valence_manifold=true;
        valence.chemistry.valence_weight=1.0;
        main_group.orbitals.push_back(std::move(valence));
        cov::MODiagramOptions main_group_options;
        const auto main_group_data=cov::build_mo_diagram_data(
            main_group,main_group_options);
        if (main_group_data.levels.size()!=1u ||
            main_group_data.levels.front().annotation.bonding_class!=
                cov::BondingClass::Unclassified ||
            main_group_data.levels.front().annotation.bonding_confidence!=0.0) {
            std::cerr<<"unavailable metal-ligand evidence became nonbonding\n";
            return 21;
        }
    }

    // One canonical MO may contribute to both a geometry-qualified/generic
    // multicentre subspace and a wider delocalised-pi family.  Their atom and
    // electron counts are independent evidence and must never overwrite one
    // another (the real acrolein fixture exercises this overlap).
    {
        cov::MolecularOrbital overlap;
        auto& chemistry=overlap.chemistry;
        chemistry.available=true;
        chemistry.confidence=0.88;
        chemistry.channel.dominant=cov::OrbitalAngularFamily::Pi;
        chemistry.channel.status=cov::ChemistryStatus::Determined;
        chemistry.multicentre_label="3c2e";
        chemistry.multicentre_assignment_available=true;
        chemistry.multicentre_participating_atoms=3u;
        chemistry.multicentre_participating_electrons=2.0;
        chemistry.multicentre_participating_atom_indices={0u,1u,3u};
        chemistry.multicentre_confidence=0.61;
        chemistry.delocalised_family_id="pi-overlap";
        chemistry.delocalised_family_orbitals={1u,3u,5u,6u};
        chemistry.delocalised_participating_atoms=4u;
        chemistry.delocalised_participating_electrons=4.0;
        chemistry.delocalised_participating_atom_indices={0u,1u,2u,3u};
        chemistry.delocalised_pi_confidence=0.93;
        chemistry.delocalised_pi_weight=0.8;
        const auto annotation=cov::annotate_orbital(overlap);
        if (!annotation.multicentre.available ||
            annotation.multicentre.centres!=3u ||
            std::abs(annotation.multicentre.electrons-2.0)>1.0e-8 ||
            annotation.multicentre.atom_indices!=
                std::vector<std::size_t>{0u,1u,3u} ||
            std::abs(annotation.multicentre.confidence-0.61)>1.0e-8 ||
            !annotation.delocalised_pi.available ||
            annotation.delocalised_pi.participating_atoms!=4u ||
            std::abs(annotation.delocalised_pi.participating_electrons-4.0)>
                1.0e-8 ||
            annotation.delocalised_pi.atom_indices!=
                std::vector<std::size_t>{0u,1u,2u,3u} ||
            std::abs(annotation.delocalised_pi.confidence-0.93)>1.0e-8) {
            std::cerr<<"multicentre/pi overlap metadata was conflated\n";
            return 22;
        }
    }

    cov::Wavefunction wf;
    wf.atoms.resize(4);

    cov::MolecularOrbital core;
    core.energy_hartree = -3.0;
    core.occupation = 2.0f;
    core.symmetry = "core";
    wf.orbitals.push_back(core);

    cov::MolecularOrbital occupied;
    occupied.energy_hartree = -0.5;
    occupied.occupation = 2.0f;
    occupied.symmetry = "sigma_g bonding";
    wf.orbitals.push_back(occupied);

    cov::MolecularOrbital virtual_a;
    virtual_a.energy_hartree = 0.2;
    virtual_a.occupation = 0.0f;
    virtual_a.symmetry = "pi_u";
    wf.orbitals.push_back(virtual_a);

    cov::MolecularOrbital virtual_b = virtual_a;
    virtual_b.energy_hartree = 0.200000004;
    wf.orbitals.push_back(virtual_b);

    cov::MolecularOrbital producer_multicentre;
    producer_multicentre.energy_hartree = 0.45;
    producer_multicentre.occupation = 0.0f;
    producer_multicentre.symmetry = "sigma 3c2e antibonding";
    wf.orbitals.push_back(producer_multicentre);

    cov::MolecularOrbital very_high;
    very_high.energy_hartree = 7.0;
    very_high.occupation = 0.0f;
    very_high.symmetry = "?";
    wf.orbitals.push_back(very_high);

    cov::MODiagramOptions options;
    options.selected_index = 1;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 760;
    options.height = 720;
    options.filter.core_energy_cutoff_hartree = -1.5;
    options.filter.virtual_window_hartree = 1.5;
    options.degeneracy.tolerance_hartree = 1.0e-5;
    // Keep all three low-lying virtuals so the explicit producer 3c2e fixture
    // participates in the human and machine export checks.
    options.max_levels = 7;

    const auto data = cov::build_mo_diagram_data(wf, options);
    if (!data.frontier.homo || !data.frontier.lumo) return 1;
    if (data.mode != cov::MODiagramMode::ValenceCentral || data.levels.size() != 4) {
        std::cerr << "valence-central selection failed\n";
        return 2;
    }
    if (data.levels[1].metadata.display_label != "3-a" ||
        data.levels[2].metadata.display_label != "3-b") {
        std::cerr << "degenerate labels failed\n";
        return 3;
    }
    if (data.levels[0].electrons.alpha != 1 || data.levels[0].electrons.beta != 1) {
        std::cerr << "electron population failed\n";
        return 4;
    }
    // Preserve producer text as parsed metadata, but this fixture supplies no
    // coefficients/operator that could certify the displayed group's role.
    const auto parsed_occupied=cov::annotate_orbital(occupied);
    if (parsed_occupied.bonding_class!=cov::BondingClass::Bonding ||
        parsed_occupied.bonding_source!=cov::AnnotationSource::ParsedLabel ||
        data.levels[0].annotation.family != "sigma" ||
        data.levels[0].annotation.bonding_class != cov::BondingClass::Unclassified ||
        data.levels[0].annotation.bonding_source != cov::AnnotationSource::Unavailable ||
        data.levels[1].annotation.family != "pi") {
        std::cerr << "parsed family metadata or missing-operator group-role separation failed\n";
        return 5;
    }
    if (!data.levels[3].annotation.multicentre.available ||
        data.levels[3].annotation.multicentre.label != "3c2e" ||
        data.levels[3].annotation.multicentre.source != cov::AnnotationSource::ParsedLabel) {
        std::cerr << "explicit 3c2e producer annotation failed\n";
        return 6;
    }

    // Compact active-space modes contain exactly the derived family members;
    // changing the MO inspected in 3-D must not change the compact row set.
    {
        cov::Wavefunction active;
        for (int i=0;i<4;++i) {
            cov::Atom atom;
            atom.symbol="C";
            atom.atomic_number=6;
            atom.x=2.5*static_cast<double>(i);
            active.atoms.push_back(atom);
        }
        for (std::size_t i=0;i<6u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.6+0.18*static_cast<double>(i);
            orbital.occupation=i<3u?2.0f:0.0f;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.valence_weight=0.9;
            active.orbitals.push_back(std::move(orbital));
        }
        for (std::uint32_t i=0;i<3u;++i) {
            cov::BondOrderRecord bond;
            bond.atom_a=i;
            bond.atom_b=i+1u;
            bond.mayer_order=1.0;
            active.bond_orders.push_back(bond);
        }
        active.bond_order_provenance=cov::DataProvenance::Derived;
        cov::DelocalisedPiAssignment pi;
        pi.atoms={0u,1u,2u,3u};
        pi.orbitals={1u,3u,5u};
        pi.electron_count=2.0;
        pi.topology=cov::DelocalisedPiTopology::Path;
        pi.provenance=cov::DataProvenance::Derived;
        active.delocalised_pi_assignments.push_back(pi);
        cov::MulticentreAssignment generic_fallback;
        generic_fallback.kind=cov::MulticentreKind::ThreeCentreTwoElectron;
        generic_fallback.atoms={0u,1u,2u};
        generic_fallback.orbitals={1u,4u,5u};
        generic_fallback.electron_count=2.0;
        generic_fallback.provenance=cov::DataProvenance::Derived;
        generic_fallback.rationale=
            "unique three-orbital active subspace; generic fallback";
        active.multicentre_assignments.push_back(generic_fallback);
        if (cov::preferred_compact_mo_diagram_mode(active,true)!=
            cov::MODiagramMode::DelocalisedPiFamilyOnly) {
            std::cerr<<"pi compact mode was not selected\n";
            return 30;
        }
        cov::MODiagramOptions compact;
        compact.mode=cov::MODiagramMode::DelocalisedPiFamilyOnly;
        compact.hide_ligand_centred_intermediates=true;
        compact.selected_index=0u;
        const auto first=cov::build_mo_diagram_data(active,compact);
        compact.selected_index=4u;
        const auto second=cov::build_mo_diagram_data(active,compact);
        const std::vector<std::size_t> expected_pi{1u,3u,5u};
        if (first.mode!=cov::MODiagramMode::DelocalisedPiFamilyOnly ||
            first.selection.included_indices!=expected_pi ||
            second.selection.included_indices!=expected_pi ||
            first.levels.size()!=3u || second.levels.size()!=3u) {
            std::cerr<<"pi compact family membership or selection invariance failed\n";
            return 31;
        }

        auto mixed_degeneracy=active;
        mixed_degeneracy.point_group_detected="C3V";
        mixed_degeneracy.orbitals[1].energy_hartree=-0.24;
        mixed_degeneracy.orbitals[2].energy_hartree=-0.24;
        mixed_degeneracy.orbitals[1].symmetry="E";
        mixed_degeneracy.orbitals[2].symmetry="E";
        mixed_degeneracy.orbitals[1].symmetry_provenance=
            cov::DataProvenance::Producer;
        mixed_degeneracy.orbitals[2].symmetry_provenance=
            cov::DataProvenance::Producer;
        const auto exact=cov::build_mo_diagram_data(
            mixed_degeneracy,compact);
        std::vector<std::size_t> exact_rows;
        for (const auto& level:exact.levels) {
            exact_rows.insert(exact_rows.end(),level.member_indices.begin(),
                              level.member_indices.end());
        }
        std::sort(exact_rows.begin(),exact_rows.end());
        if (exact_rows!=expected_pi ||
            std::any_of(exact.levels.begin(),exact.levels.end(),
                [](const auto& level) {
                    return std::find(level.member_indices.begin(),
                                     level.member_indices.end(),2u)!=
                           level.member_indices.end();
                })) {
            std::cerr<<"near-degenerate outsider leaked into compact family rows\n";
            return 38;
        }
        const auto exact_json=std::filesystem::temp_directory_path()/
            "cov_exact_active_membership.json";
        std::string exact_error;
        if (!cov::write_mo_diagram_json(
                exact,compact,exact_json,&exact_error)) {
            std::cerr<<"exact active-space JSON export failed: "
                     <<exact_error<<'\n';
            return 39;
        }
        std::ifstream exact_json_file(exact_json,std::ios::binary);
        std::ostringstream exact_json_buffer;
        exact_json_buffer<<exact_json_file.rdbuf();
        const auto exact_json_text=exact_json_buffer.str();
        std::error_code exact_ec;
        std::filesystem::remove(exact_json,exact_ec);
        if (exact_json_text.find("\"member_orbitals\":[1,2]")!=
                std::string::npos ||
            exact_json_text.find("\"member_orbitals\":[1]")==
                std::string::npos) {
            std::cerr<<"compact JSON diagram rows leaked a family outsider\n";
            return 40;
        }

        // Exact-family rebuilding happens after unrestricted alpha/beta
        // spatial-row collapse.  The compact row must retain the matched beta
        // occupation and local symmetry as well as the two arrow glyphs.
        auto active_spin=synthetic_spin_matching_wavefunction(2u,"C2v");
        active_spin.atoms.resize(3u);
        for (std::size_t atom=0u;atom<active_spin.atoms.size();++atom) {
            active_spin.atoms[atom].symbol="C";
            active_spin.atoms[atom].atomic_number=6;
            active_spin.atoms[atom].x=2.5*static_cast<double>(atom);
        }
        add_spin_matching_orbital(
            active_spin,cov::Spin::Alpha,-0.20,"A1",{1.0f,0.0f},1);
        add_spin_matching_orbital(
            active_spin,cov::Spin::Alpha, 0.10,"B1",{0.0f,1.0f},1);
        add_spin_matching_orbital(
            active_spin,cov::Spin::Beta,-0.19,"A1",{1.0f,0.0f},1);
        add_spin_matching_orbital(
            active_spin,cov::Spin::Beta, 0.11,"B1",{0.0f,1.0f},1);
        active_spin.orbitals[1].occupation=0.0f;
        active_spin.orbitals[3].occupation=0.0f;
        active_spin.alpha_electrons=1u;
        active_spin.beta_electrons=1u;
        for (std::uint32_t atom=0u;atom<2u;++atom) {
            cov::BondOrderRecord bond;
            bond.atom_a=atom;
            bond.atom_b=atom+1u;
            bond.mayer_order=1.0;
            active_spin.bond_orders.push_back(bond);
        }
        active_spin.bond_order_provenance=cov::DataProvenance::Derived;
        cov::DelocalisedPiAssignment spin_pi;
        spin_pi.family_id="udft-pi-family";
        spin_pi.atoms={0u,1u,2u};
        spin_pi.orbitals={0u,1u,2u,3u};
        spin_pi.electron_count=2.0;
        spin_pi.topology=cov::DelocalisedPiTopology::Path;
        spin_pi.provenance=cov::DataProvenance::Derived;
        active_spin.delocalised_pi_assignments.push_back(spin_pi);
        const auto spin_compact=cov::build_mo_diagram_data(
            active_spin,compact);
        const auto occupied_spatial=std::find_if(
            spin_compact.levels.begin(),spin_compact.levels.end(),
            [](const auto& level) {
                return !level.member_indices.empty() &&
                       level.member_indices.front()==0u;
            });
        if (!spin_compact.spin_counterparts_collapsed ||
            occupied_spatial==spin_compact.levels.end() ||
            occupied_spatial->member_spin_counterparts.size()!=1u ||
            occupied_spatial->member_spin_counterparts.front()!=2u ||
            occupied_spatial->electrons.alpha!=1 ||
            occupied_spatial->electrons.beta!=1 ||
            std::abs(occupied_spatial->total_occupation-2.0)>1.0e-6 ||
            std::abs(occupied_spatial->metadata.occupation-2.0f)>1.0e-6f ||
            occupied_spatial->metadata.symmetry!="A1") {
            std::cerr<<"UDFT compact row lost beta occupation or symmetry\n";
            return 41;
        }

        active.delocalised_pi_assignments.clear();
        active.multicentre_assignments.clear();
        active.multicentre_assignments.push_back(generic_fallback);
        if (cov::preferred_compact_mo_diagram_mode(active,true)!=
            cov::MODiagramMode::ValenceCentral) {
            std::cerr<<"generic three-centre candidate incorrectly selected "
                       "multicentre-only mode\n";
            return 32;
        }
        active.multicentre_assignments.clear();
        cov::MulticentreAssignment multicentre;
        multicentre.kind=cov::MulticentreKind::ThreeCentreFourElectron;
        multicentre.atoms={0u,1u,2u};
        multicentre.orbitals={0u,2u,4u};
        multicentre.electron_count=4.0;
        multicentre.provenance=cov::DataProvenance::Derived;
        multicentre.source_subspace_id="synthetic-3c4e";
        multicentre.geometry_qualified_framework=true;
        multicentre.rationale=
            "geometry-qualified three-centre framework; synthetic fixture";
        active.multicentre_assignments.push_back(multicentre);
        if (cov::preferred_compact_mo_diagram_mode(active,true)!=
            cov::MODiagramMode::MulticentreActiveSpaceOnly) {
            std::cerr<<"multicentre compact mode was not selected\n";
            return 33;
        }
        compact.mode=cov::MODiagramMode::MulticentreActiveSpaceOnly;
        const auto mc=cov::build_mo_diagram_data(active,compact);
        const std::vector<std::size_t> expected_mc{0u,2u,4u};
        if (mc.selection.included_indices!=expected_mc || mc.levels.size()!=3u) {
            std::cerr<<"multicentre compact membership failed\n";
            return 34;
        }

        // A non-cyclic, star-shaped ligand p family is not an automatic pi
        // focus.  This protects BF3/XeF4-like main-group geometries.
        active.multicentre_assignments.clear();
        active.bond_orders.clear();
        active.atoms[0].x=0.0; active.atoms[0].y=0.0;
        active.atoms[1].x=2.5; active.atoms[1].y=0.0;
        active.atoms[2].x=-1.25; active.atoms[2].y=2.165063509461;
        active.atoms[3].x=-1.25; active.atoms[3].y=-2.165063509461;
        for (std::uint32_t leaf=1u;leaf<4u;++leaf) {
            cov::BondOrderRecord bond;
            bond.atom_a=0u;
            bond.atom_b=leaf;
            bond.mayer_order=1.0;
            active.bond_orders.push_back(bond);
        }
        active.delocalised_pi_assignments.push_back(pi);
        active.delocalised_pi_assignments.back().topology=
            cov::DelocalisedPiTopology::BranchedResonance;
        active.delocalised_pi_assignments.back().
            branch_centre_projected_occupation=0.78;
        if (cov::preferred_compact_mo_diagram_mode(active,true)!=
            cov::MODiagramMode::ValenceCentral) {
            std::cerr<<"star-shaped p family incorrectly selected pi-only mode\n";
            return 35;
        }
        active.delocalised_pi_assignments.back().
            branch_centre_projected_occupation=1.10;
        if (cov::preferred_compact_mo_diagram_mode(active,true)!=
            cov::MODiagramMode::DelocalisedPiFamilyOnly) {
            std::cerr<<"electronically occupied resonance star was not focused\n";
            return 36;
        }

        // Haptic pi chemistry takes precedence over a geometrically available
        // metal first shell; otherwise a valid eta-ring assignment can be
        // hidden merely because the same contacts also fit a CN template.
        auto haptic_lf=synthetic_oh_environment();
        haptic_lf.bond_order_provenance=cov::DataProvenance::Derived;
        for (std::size_t i=0u;i<6u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.5+0.2*static_cast<double>(i);
            orbital.occupation=i<3u?2.0f:0.0f;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.valence_weight=0.9;
            haptic_lf.orbitals.push_back(std::move(orbital));
        }
        cov::DelocalisedPiAssignment haptic;
        haptic.atoms={1u,2u,3u,4u,5u,6u};
        haptic.orbitals={0u,2u,4u};
        haptic.electron_count=4.0;
        haptic.topology=cov::DelocalisedPiTopology::HapticMetal;
        haptic.provenance=cov::DataProvenance::Derived;
        haptic_lf.delocalised_pi_assignments.push_back(std::move(haptic));
        if (!cov::analyse_ligand_field_environment(haptic_lf).available() ||
            cov::preferred_compact_mo_diagram_mode(haptic_lf,true)!=
                cov::MODiagramMode::DelocalisedPiFamilyOnly) {
            std::cerr<<"haptic pi focus lost to ligand-field precedence\n";
            return 37;
        }
    }

    const auto build_id=std::hash<std::string>{}(
        std::filesystem::absolute(std::filesystem::current_path()).string());
    const auto temp_stem=std::string("cov_mo_diagram_smoke_")+
                         std::to_string(build_id);
    const auto temp=std::filesystem::temp_directory_path()/
                    (temp_stem+".fchk");
    const auto normalised_temp=std::filesystem::temp_directory_path()/temp_stem;
    const auto result = cov::export_mo_diagram_bundle(wf, options, temp);
    if (!result.svg || !result.png || !result.json || !result.csv) {
        std::cerr << "diagram export failed: " << result.error << '\n';
        return 7;
    }

    const auto png_path = std::filesystem::path(normalised_temp.string() + ".mo.png");
    const auto svg_path = std::filesystem::path(normalised_temp.string() + ".mo.svg");
    const auto json_path = std::filesystem::path(normalised_temp.string() + ".mo.json");
    const auto csv_path = std::filesystem::path(normalised_temp.string() + ".mo.csv");
    if (result.png_path!=png_path || result.svg_path!=svg_path ||
        result.json_path!=json_path || result.csv_path!=csv_path ||
        result.png_path.string().find(".fchk.mo.")!=std::string::npos) {
        std::cerr << "normalised export bundle paths failed\n";
        return 20;
    }
    std::ifstream png(png_path, std::ios::binary);
    std::array<std::uint8_t, 8> signature{};
    png.read(reinterpret_cast<char*>(signature.data()), static_cast<std::streamsize>(signature.size()));
    const std::array<std::uint8_t, 8> expected{137,80,78,71,13,10,26,10};
    if (signature != expected) return 8;

    std::ifstream svg_file(svg_path, std::ios::binary);
    std::ostringstream svg_buffer;
    svg_buffer << svg_file.rdbuf();
    const std::string svg = svg_buffer.str();
    if (svg.find("Valence MO diagram") == std::string::npos ||
        svg.find("MO numbering is intentionally omitted") != std::string::npos ||
        svg.find("data-orbital-index=") == std::string::npos ||
        svg.find(">3-a<") != std::string::npos) {
        std::cerr << "human export policy failed\n";
        return 9;
    }

    std::ifstream json_file(json_path, std::ios::binary);
    std::ostringstream json_buffer;
    json_buffer << json_file.rdbuf();
    const std::string json = json_buffer.str();
    if (json.find("\"mode\": \"valence-central\"") == std::string::npos ||
        json.find("\"label\": \"3-a\"") == std::string::npos ||
        json.find("\"multicentre_label\": \"3c2e\"") == std::string::npos ||
        json.find("\"strict_salc_claimed\": false") == std::string::npos) {
        std::cerr << "machine-readable markers missing\n";
        return 10;
    }

    // A bundle is all-or-nothing.  A bad destination discovered during
    // preflight must not replace any already-existing sibling export.
    {
        std::ofstream sentinel(svg_path,std::ios::binary|std::ios::trunc);
        sentinel<<"preserve-existing-svg";
    }
    std::error_code ec;
    std::filesystem::remove(png_path,ec);
    ec.clear();
    std::filesystem::create_directory(png_path,ec);
    const auto rejected=cov::export_mo_diagram_bundle(wf,options,temp);
    std::ifstream sentinel(svg_path,std::ios::binary);
    std::ostringstream sentinel_buffer;
    sentinel_buffer<<sentinel.rdbuf();
    if (rejected.svg || rejected.png || rejected.json || rejected.csv ||
        rejected.error.empty() ||
        sentinel_buffer.str()!="preserve-existing-svg") {
        std::cerr<<"transactional export preflight modified an existing sibling\n";
        return 21;
    }
    std::filesystem::remove(png_path,ec);

    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(csv_path, ec);

    if (!failed_real_fchk.empty()) {
        return 11;
    }

    std::cout << "mo_diagram_smoke ok\n";
    return 0;
}
