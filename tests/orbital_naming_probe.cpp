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

#include "cov/nbo_aomo_labels.hpp"
#include "cov/wavefunction_io.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
void stats(const std::vector<cov::ui::NboAomoName>& names){std::size_t unknown=0,unordered=0,ordered=0,question=0,decomposed=0,unanalysed=0;for(const auto& n:names){if(n.decomposition_verified)++decomposed;if(!n.verified&&!n.decomposition_verified)++unanalysed;if(!n.verified||n.irrep.empty())++unknown;else if(n.ordinal)++ordered;else ++unordered;if(n.verified&&!n.irrep.empty()&&n.label.starts_with('?'))++question;}std::cout<<"{\"total\":"<<names.size()<<",\"decomposition_verified\":"<<decomposed<<",\"unanalysed\":"<<unanalysed<<",\"irrep_unknown\":"<<unknown<<",\"ordinal_unknown\":"<<unordered<<",\"proven_ordinals\":"<<ordered<<",\"question_plus_known_irrep\":"<<question<<",\"unknown_irrep_indices1\":[";bool comma=false;for(std::size_t i=0;i<names.size();++i)if(!names[i].verified||names[i].irrep.empty()){if(comma)std::cout<<',';comma=true;std::cout<<i+1;}std::cout<<"],\"unordered\":[";comma=false;for(std::size_t i=0;i<names.size();++i)if(names[i].verified&&!names[i].irrep.empty()&&!names[i].ordinal){if(comma)std::cout<<',';comma=true;std::cout<<"{\"index1\":"<<i+1<<",\"label\":"<<std::quoted(names[i].label)<<",\"irrep\":"<<std::quoted(names[i].irrep)<<",\"multiplicity\":"<<names[i].representation_multiplicity<<'}';}std::cout<<"],\"rows\":[";for(std::size_t i=0;i<names.size();++i){const auto& n=names[i];if(i)std::cout<<',';std::cout<<"{\"index1\":"<<i+1<<",\"label\":"<<std::quoted(n.label)<<",\"irrep\":"<<std::quoted(n.irrep)<<",\"status\":"<<std::quoted(n.status)<<",\"ordinal\":"<<n.ordinal<<",\"multiplicity\":"<<n.representation_multiplicity<<",\"detail\":"<<std::quoted(n.detail)<<",\"containing_members0\":[";for(std::size_t j=0;j<n.containing_members.size();++j){if(j)std::cout<<',';std::cout<<n.containing_members[j];}std::cout<<"],\"content\":[";for(std::size_t j=0;j<n.containing_irreps.size();++j){if(j)std::cout<<',';const auto& c=n.containing_irreps[j];std::cout<<"{\"irrep\":"<<std::quoted(c.irrep)<<",\"dimension\":"<<c.dimension<<",\"copies\":"<<c.multiplicity<<'}';}std::cout<<"],\"decomposition_verified\":"<<(n.decomposition_verified?"true":"false")<<",\"decomposition_status\":"<<std::quoted(n.decomposition_status)<<",\"components\":[";for(std::size_t j=0;j<n.components.size();++j){if(j)std::cout<<',';const auto& c=n.components[j];std::cout<<"{\"irrep\":"<<std::quoted(c.irrep)<<",\"dimension\":"<<c.dimension<<",\"weight\":"<<std::setprecision(17)<<c.weight<<'}';}std::cout<<"]}";}std::cout<<"]}";}
int main(int argc,char** argv){try{if(argc!=3)throw std::runtime_error("canonical.fchk analysis_directory required");const auto w=cov::parse_wavefunction(argv[1]);const auto before=w;const auto found=cov::discover_nbo_inputs({argv[2]});if(found.candidates.size()!=1)throw std::runtime_error("exactly one attachment required");const auto data=cov::read_nbo_integration(w,found.candidates.front());const auto source=cov::build_nbo_salc_model(w,data);const auto names=cov::ui::build_nbo_aomo_names(w,data,&source);const bool two=std::any_of(source.orbitals.begin(),source.orbitals.end(),[](const auto& o){return o.spin!=cov::NboSpin::Total;});std::cout<<"{\"group\":"<<std::quoted(source.point_group)<<",\"source\":{\"canonical\":";stats(names.canonical);std::cout<<",\"salc\":";stats(names.salc);std::cout<<"},\"teaching\":";if(two){const auto teaching=cov::build_nbo_spin_averaged_model(w,data,source);const auto simple=cov::ui::build_nbo_aomo_names(w,data,&teaching);std::cout<<"{\"spin_averaged\":"<<(teaching.spin_averaged?"true":"false")<<",\"canonical\":";stats(simple.canonical);std::cout<<",\"salc\":";stats(simple.salc);std::cout<<'}';}else std::cout<<"null";bool unchanged=w.ao_overlap==before.ao_overlap;for(std::size_t i=0;i<w.orbitals.size();++i){const auto& a=w.orbitals[i];const auto& b=before.orbitals[i];unchanged=unchanged&&a.coefficients==b.coefficients&&a.gaussian_source_coefficients==b.gaussian_source_coefficients&&a.energy_hartree==b.energy_hartree&&a.occupation==b.occupation&&a.spin==b.spin&&a.source_orbital_index==b.source_orbital_index&&a.symmetry==b.symmetry&&a.symmetry_provenance==b.symmetry_provenance;}std::cout<<",\"canonical_immutable\":"<<(unchanged?"true":"false")<<'}';return unchanged?0:2;}catch(const std::exception& e){std::cerr<<e.what();return 1;}}
