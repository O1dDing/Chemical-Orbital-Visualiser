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

#include "cov/mo_diagram.hpp"
#include "cov/orbital_symmetry_scope.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>

int main(int argc,char** argv) {
    using namespace cov;
    if(argc!=2)return 2;
    const auto check=[](bool ok,const char* name) {
        std::cout<<"{\"check\":\""<<name<<"\",\"passed\":"<<(ok?"true":"false")<<"}\n";
    };
    PiPartnerComponents a{-.2,.9,.05,.1,.001},b{-.195,.85,.07,.1,-.001};
    auto weak=assess_weak_crystal_field(a,b);
    check(weak.input_valid&&weak.accepted,"weak-d-components-accepted-without-pi-channel");
    check(std::abs(weak.splitting_hartree-.005)<1e-15,"gap-from-actual-energies");
    auto reversed=assess_weak_crystal_field(b,a);
    check(reversed.accepted&&reversed.splitting_hartree==weak.splitting_hartree,"energy-endpoint-order-invariant");
    auto no_d=a;no_d.metal_d=.59;
    check(!assess_weak_crystal_field(no_d,b).accepted,"insufficient-metal-d-rejected");
    auto strong=a;strong.metal_ligand_overlap=.026;
    check(!assess_weak_crystal_field(strong,b).accepted,"strong-local-mixing-not-weak");
    auto far=b;far.energy_hartree=.1;
    check(!assess_weak_crystal_field(a,far).accepted,"large-gap-not-weak");
    check(!assess_weak_crystal_field(a,b,.004,.025).accepted,"cf-gap-threshold-applied");
    check(!assess_weak_crystal_field(a,b,.020,.0005).accepted,"cf-mixing-threshold-applied");
    check(assess_weak_crystal_field(a,a,0,.025).accepted,"zero-threshold-exact-zero-gap");
    check(!assess_weak_crystal_field(a,b,0,.025).accepted,"zero-threshold-nonzero-gap-rejected");
    check(!assess_weak_crystal_field(a,b,-1,.025).input_valid,"negative-threshold-invalid");
    auto invalid=a;invalid.energy_hartree=std::numeric_limits<double>::quiet_NaN();
    check(!assess_weak_crystal_field(invalid,b).input_valid,"nan-energy-invalid");
    invalid=a;invalid.metal_d=std::numeric_limits<double>::infinity();
    check(!assess_weak_crystal_field(invalid,b).input_valid,"infinite-weight-invalid");
    invalid=a;invalid.metal_ligand_overlap=std::numeric_limits<double>::quiet_NaN();
    check(!assess_weak_crystal_field(invalid,b).input_valid,"nan-overlap-invalid");
    PiPartnerComponents donor{-.4,.1,.8,.8,.05},acceptor{-.15,.8,.1,.8,-.05};
    auto pi=assess_pi_partner(donor,acceptor,LigandPiPrior::SigmaOnly);
    check(pi.accepted&&pi.direction==PiPairDirection::Donor&&pi.prior_relation=="contradicted",
          "actual-pi-partner-survives-sigma-only-prior");
    check(!assess_pi_partner(a,b,LigandPiPrior::Donor).accepted,"cf-only-components-are-not-pi-partners");
    MODiagramData data;
    MODiagramOptions options;options.selected_index=0;options.width=1200;options.height=900;
    const double energies[]={-.2,-.195,-.4,-.15};
    const char* labels[]={"E","T2","T2","T2"};
    for(std::size_t i=0;i<4;++i) {
        OrbitalMetadata metadata;metadata.orbital_index=i;metadata.raw_mo_number=i+1;metadata.energy_hartree=energies[i];
        metadata.symmetry=labels[i];
        metadata.symmetry_view.label=labels[i];
        metadata.symmetry_view.point_group="Td";
        metadata.symmetry_view.point_group_basis="local-coordination-template";
        metadata.symmetry_view.origin=OrbitalSymmetryOrigin::LocalDimensionCandidate;
        metadata.symmetry_view.orbital_indices={i};
        metadata.symmetry_view.atom_indices={0,1,2,3,4};
        metadata.symmetry_view.axes_available=true;
        metadata.symmetry_view.rotation_reference_to_input={1,0,0,0,1,0,0,0,1};
        data.metadata.push_back(metadata);
        data.annotations.emplace_back();
        MODiagramLevel level;level.metadata=metadata;level.layout_energy_hartree=energies[i];
        level.member_indices={i};data.levels.push_back(level);
    }
    data.energy_transform=build_energy_transform({-.4,-.2,-.195,-.15},EnergyAxisMode::Linear,.055);
    options.energy_axis_mode=EnergyAxisMode::Linear;
    auto make_gap=[&](std::size_t lower,std::size_t upper,OrbitalEnergyGapKind kind) {
        OrbitalEnergyGapDescriptor gap;gap.gap_kind=kind;gap.lower_level=lower;gap.upper_level=upper;
        gap.retained_level=lower;gap.lower_orbitals={lower};gap.upper_orbitals={upper};
        gap.lower_symmetry_scope=data.metadata[lower].symmetry_view;
        gap.upper_symmetry_scope=data.metadata[upper].symmetry_view;
        gap.lower_energy_hartree=energies[lower];gap.upper_energy_hartree=energies[upper];
        gap.splitting_hartree=energies[upper]-energies[lower];
        gap.lower_energy_spread_hartree=gap.upper_energy_spread_hartree=0;
        gap.weak_split_threshold_hartree=.020;gap.weak_overlap_threshold=.025;
        return gap;
    };
    auto cf_gap=make_gap(0,1,OrbitalEnergyGapKind::CrystalField);
    cf_gap.symmetry="E/T2";cf_gap.kind=PiInteractionKind::WeakNearNonbonding;
    cf_gap.confidence=weak.support_score;
    cf_gap.crystal_field_evidence=std::make_shared<const WeakCrystalFieldAssessment>(weak);
    data.crystal_field_gaps.push_back(cf_gap);
    auto pi_gap=make_gap(2,3,OrbitalEnergyGapKind::PiPartner);
    pi_gap.symmetry="T2";pi_gap.kind=PiInteractionKind::Donor;pi_gap.confidence=pi.support_score;
    pi_gap.orbital_evidence=std::make_shared<const PiPartnerAssessment>(pi);
    data.pi_interactions.push_back(pi_gap);
    check(data.pi_interactions.size()==1&&data.crystal_field_gaps.size()==1,"separate-canonical-counts");
    auto all=orbital_energy_gaps(data);
    check(all.size()==2&&all[0]==&data.pi_interactions[0]&&all[1]==&data.crystal_field_gaps[0],
          "render-union-borrows-current-records");
    const auto snapshot=make_mo_diagram_view_snapshot(data,options,0,"energy-gap-control");
    const auto bundle=export_mo_diagram_bundle(snapshot,std::filesystem::path(argv[1])/"gap-control");
    check(bundle.svg&&bundle.png&&bundle.json&&bundle.csv&&bundle.error.empty(),"actual-shared-export-bundle");
    auto special=pi;special.detail="line\n\"quoted\"\tcontrol";special.prior_relation="a\rb";
    std::cout<<"{\"record\":\"escaped-pi\",\"value\":"<<pi_partner_assessment_json(special)<<"}\n";
    std::cout<<"{\"record\":\"cf\",\"value\":"<<orbital_energy_gap_json(cf_gap)<<"}\n";
    std::cout<<"{\"record\":\"pi\",\"value\":"<<orbital_energy_gap_json(pi_gap)<<"}\n";
    std::cout<<"{\"record\":\"invalid-cf\",\"value\":"
             <<weak_crystal_field_assessment_json(assess_weak_crystal_field(invalid,b))<<"}\n";
    const EnergyUnit units[]={EnergyUnit::Hartree,EnergyUnit::ElectronVolt,EnergyUnit::JoulePerMol,
        EnergyUnit::KilojoulePerMol,EnergyUnit::CaloriePerMol,EnergyUnit::KilocaloriePerMol};
    const char* labels_by_unit[]={"Ha","eV","J/mol","kJ/mol","cal/mol","kcal/mol"};
    for(std::size_t i=0;i<6;++i) {
        const auto value=orbital_energy_gap_json(pi_gap,units[i]);
        const std::string needle=std::string("\"display_unit\":\"")+labels_by_unit[i]+'"';
        check(value.find(needle)!=std::string::npos,(std::string("gap-unit-label-")+labels_by_unit[i]).c_str());
        std::cout<<"{\"record\":\"unit-"<<i<<"\",\"value\":"<<value<<"}\n";
        auto unit_options=options;unit_options.energy_unit=units[i];
        const auto unit_snapshot=make_mo_diagram_view_snapshot(data,unit_options,0,"energy-gap-unit-control");
        const auto unit_bundle=export_mo_diagram_bundle(unit_snapshot,
            std::filesystem::path(argv[1])/("gap-unit-"+std::to_string(i)));
        check(unit_bundle.svg&&unit_bundle.png&&unit_bundle.json&&unit_bundle.csv&&unit_bundle.error.empty(),
            (std::string("gap-unit-export-")+labels_by_unit[i]).c_str());
    }
    return 0;
}
