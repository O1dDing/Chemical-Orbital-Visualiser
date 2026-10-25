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

#include "cov/orbital_symmetry_scope.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/ligand_field.hpp"
#include "cov/gaussian_log.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok,const char* name) {
    std::cout<<"{\"check\":\""<<name<<"\",\"passed\":"<<(ok?"true":"false")<<"}\n";
}
}
int main() {
    using namespace cov;
    Wavefunction wf;
    wf.atoms.resize(3);
    wf.orbitals.resize(3);
    wf.orbitals[0].symmetry="A1"; wf.orbitals[1].symmetry="B2"; wf.orbitals[2].symmetry="A1";
    for(auto& mo:wf.orbitals)mo.symmetry_provenance=DataProvenance::Producer;
    wf.orbitals[2].spin=Spin::Beta;
    wf.point_group_detected="Oh"; wf.point_group_used="D2h";
    OrbitalSymmetrySourceRecord record;
    record.source_path="producer\n\"source\".log"; record.line_begin=12;record.line_end=18;record.job_segment=4;
    record.detected_group_context="Oh";record.abelian_group_context="D2h";
    record.orbital_indices={0,1,2};record.labels={"A1","B2","A1"};
    wf.orbital_symmetry_source_records={record};
    auto original=molecular_orbital_symmetry(wf,0);
    require(original.label=="A1" && original.origin==OrbitalSymmetryOrigin::Producer,"producer-label-retained");
    require(original.point_group.empty() && original.point_group_basis=="unresolved","printed-group-binding-not-invented");
    require(original.producer_detected_group=="Oh" && original.producer_abelian_group=="D2h","distinct-producer-group-context");
    require(original.source_path==record.source_path && original.source_line_begin==12 &&
            original.source_line_end==18 && original.source_job_segment==4,"exact-source-block");
    const std::array<std::size_t,2> pair{0,1};
    auto mixed=aggregate_molecular_symmetry(wf,pair);
    require(mixed.origin==OrbitalSymmetryOrigin::MixedMembers && mixed.orbital_indices==std::vector<std::size_t>{0,1},
            "mixed-global-members-not-one-representative-label");
    LigandFieldEnvironment environment;
    environment.geometry_id=GeometryId::Octahedral6; environment.metal_atom=0;
    environment.ligand_atoms={1,2};environment.angular_rms=.12;environment.shape_measure=.04;
    auto dimension=classify_local_irrep_by_dimension("Oh",MetalAOShell::D,2);
    if(!dimension)throw std::runtime_error("catalogue control unavailable");
    auto local=local_orbital_symmetry(*dimension,environment,pair);
    require(local.origin==OrbitalSymmetryOrigin::LocalDimensionCandidate && orbital_symmetry_is_candidate(local),"dimension-marked-candidate");
    require(local.local_assignment && std::isnan(local.local_assignment->confidence) &&
            !local.local_assignment->projection,"dimension-no-false-numerical-proof");
    require(wf.orbitals[0].symmetry=="A1" && original.label=="A1","local-candidate-does-not-change-source");
    require(local.point_group=="Oh" && local.point_group_basis=="local-coordination-template" &&
            local.geometry_angular_rms==.12 && local.geometry_shape_measure==.04,"local-template-and-deviation");
    require(local.orbital_indices==std::vector<std::size_t>{0,1} && local.atom_indices==std::vector<std::size_t>{0,1,2} &&
            local.axes_available,"local-members-centre-neighbours-and-axes");
    require(!compatible_symmetry_scopes(original,local),"matching-text-does-not-join-global-and-local");
    auto rotated=local;rotated.rotation_reference_to_input={0,-1,0,1,0,0,0,0,1};
    require(!compatible_symmetry_scopes(local,rotated),"different-local-axes-incompatible");
    auto other_atom=local;other_atom.atom_indices={1,0,2};
    require(!compatible_symmetry_scopes(local,other_atom),"different-central-atom-incompatible");
    auto other_group=local;other_group.point_group="D4h";
    require(!compatible_symmetry_scopes(local,other_group),"different-local-group-incompatible");
    require(compatible_symmetry_scopes(local,local),"same-local-scope-compatible");
    const std::array<std::size_t,1> target{2};
    auto spin=candidate_orbital_symmetry(local,OrbitalSymmetryOrigin::SpinCounterpartCandidate,target,.93);
    require(spin.origin==OrbitalSymmetryOrigin::SpinCounterpartCandidate && spin.orbital_indices==std::vector<std::size_t>{2},
            "spin-candidate-has-target-identity");
    require(!spin.local_assignment && spin.candidate_source && spin.candidate_source->orbital_indices==std::vector<std::size_t>{0,1},
            "source-projection-not-claimed-for-other-spin");
    local.label="changed-after-copy";
    require(spin.candidate_source->label=="Eg","candidate-source-immutable-snapshot");
    require(orbital_symmetry_compact_text(spin).find("candidate")!=std::string::npos &&
            orbital_symmetry_compact_text(spin).find("Oh")!=std::string::npos,"compact-candidate-and-group-visible");
    bool rejected=false;try { (void)candidate_orbital_symmetry(original,OrbitalSymmetryOrigin::Producer,target,1); }
    catch(const std::invalid_argument&) {rejected=true;}
    require(rejected,"propagation-cannot-promote-to-producer");
    auto pi=candidate_orbital_symmetry(*spin.candidate_source,OrbitalSymmetryOrigin::PiPartnerCandidate,target,2.7);
    require(pi.origin==OrbitalSymmetryOrigin::PiPartnerCandidate && !pi.local_assignment &&
            pi.candidate_score==2.7,"pi-ranking-score-is-not-purity");
    DerivedOrbitalSymmetryAssignment derived;
    derived.point_group="C2v";derived.label="B2";derived.orbital_indices={1};derived.axes_available=true;
    derived.principal_axis={0,0,1};derived.secondary_axis={1,0,0};derived.subspace_retention=.999;
    wf.orbitals[1].symmetry_provenance=DataProvenance::Derived;wf.derived_orbital_symmetry_assignments={derived};
    auto full=molecular_orbital_symmetry(wf,1);
    require(full.origin==OrbitalSymmetryOrigin::MolecularOperations && full.point_group=="C2v" &&
            full.molecular_assignment && full.molecular_assignment->subspace_retention==.999,"derived-label-bound-to-actual-operations");
    Wavefunction mixed_wf;mixed_wf.atoms.resize(3);mixed_wf.basis_count=5;
    mixed_wf.shells.push_back({0,0,0,0,2,1});
    mixed_wf.ao_overlap.assign(25,0);for(std::size_t i=0;i<5;++i)mixed_wf.ao_overlap[i*5+i]=1;
    mixed_wf.orbitals.resize(2);
    mixed_wf.orbitals[0].coefficients={1,0,0,0,0};
    mixed_wf.orbitals[1].coefficients={0,1,0,0,0};
    LocalAngularProjectionWorkspace mixed_workspace(mixed_wf,0,environment.rotation_reference_to_input);
    const auto decomposition=evaluate_local_orbital_symmetry(mixed_workspace,environment,pair,2);
    require(decomposition.origin==OrbitalSymmetryOrigin::LocalMetricDecomposition && decomposition.label=="mixed",
            "computed-Eg-T2g-mixture-retained");
    require(decomposition.local_decomposition && decomposition.local_decomposition->represented_spin_orbital_rank==2 &&
            std::abs(decomposition.local_decomposition->centre_projection_trace-2)<2e-11,"mixed-decomposition-keeps-numerical-proof");
    require(dimension->label=="Eg" && !decomposition.local_assignment && !orbital_symmetry_is_candidate(decomposition),
            "computed-mixture-not-overridden-by-dimension");
    auto no_metric=mixed_wf;no_metric.ao_overlap.clear();
    LocalAngularProjectionWorkspace no_metric_workspace(no_metric,0,environment.rotation_reference_to_input);
    const auto candidate=evaluate_local_orbital_symmetry(no_metric_workspace,environment,pair,2);
    require(candidate.origin==OrbitalSymmetryOrigin::LocalDimensionCandidate && candidate.local_decomposition &&
            candidate.local_decomposition->status==MetricSubspaceStatus::MissingInput,"missing-input-candidate-keeps-missing-status");
    auto duplicate=mixed_wf;duplicate.orbitals[1].coefficients=duplicate.orbitals[0].coefficients;
    LocalAngularProjectionWorkspace duplicate_workspace(duplicate,0,environment.rotation_reference_to_input);
    const auto unresolved=evaluate_local_orbital_symmetry(duplicate_workspace,environment,pair,2);
    require(unresolved.origin==OrbitalSymmetryOrigin::LocalMetricDecomposition && !unresolved.local_assignment &&
            unresolved.local_decomposition->represented_spin_orbital_rank==1,"rank-deficiency-not-overridden-by-dimension");
    std::cout<<"{\"record\":\"producer\",\"value\":"<<orbital_symmetry_json(original)<<"}\n";
    std::cout<<"{\"record\":\"spin-candidate\",\"value\":"<<orbital_symmetry_json(spin)<<"}\n";
    std::cout<<"{\"record\":\"molecular-derived\",\"value\":"<<orbital_symmetry_json(full)<<"}\n";
    std::cout<<"{\"record\":\"mixed-local\",\"value\":"<<orbital_symmetry_json(decomposition)<<"}\n";
}
