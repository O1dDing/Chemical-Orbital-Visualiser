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

#include "cov/local_orbital_symmetry.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>

void number(double x) {if(std::isfinite(x))std::cout<<x;else std::cout<<"null";}
void output(const char* name,const std::optional<cov::LocalIrrepAssignment>& a) {
    std::cout<<"{\"name\":\""<<name<<"\",\"assigned\":"<<(a?"true":"false");
    if(a) {
        std::cout<<",\"label\":\""<<a->label<<"\",\"shell\":"<<static_cast<int>(a->shell)<<",\"confidence\":";number(a->confidence);
#ifndef COV_OLD_LOCAL_BASELINE
        std::cout<<",\"source\":"<<static_cast<int>(a->source)<<",\"centre_fraction\":";number(a->centre_mean_fraction);
        std::cout<<",\"angular_local_fraction\":";number(a->angular_fraction_within_centre);
        std::cout<<",\"labelled_target_fraction\":";number(a->labelled_fraction_of_target);
        std::cout<<",\"projection_present\":"<<(a->projection?"true":"false");
        if(a->projection)std::cout<<",\"represented_spin_rank\":"<<a->projection->represented_spin_orbital_rank;
#endif
    }
    std::cout<<"}\n";
}
cov::Wavefunction p_basis(bool radial=false) {
    cov::Wavefunction wf;wf.atoms.resize(1);wf.basis_count=radial?6:3;
    wf.shells.push_back({0,0,0,0,1,0});if(radial)wf.shells.push_back({0,0,0,3,1,0});
    wf.ao_overlap.assign(wf.basis_count*wf.basis_count,0);
    for(std::size_t i=0;i<wf.basis_count;++i)wf.ao_overlap[i*wf.basis_count+i]=1;
    if(radial)for(std::size_t i=0;i<3;++i)wf.ao_overlap[i*6+i+3]=wf.ao_overlap[(i+3)*6+i]=.99;
    return wf;
}
void orbital(cov::Wavefunction& wf,std::vector<double> coefficients,cov::Spin spin=cov::Spin::Alpha) {
    cov::MolecularOrbital mo;mo.coefficients=std::move(coefficients);mo.spin=spin;wf.orbitals.push_back(std::move(mo));
}
int main() {
    std::cout<<std::setprecision(17);
    const std::array<double,9> I{1,0,0,0,1,0,0,0,1};
    const std::array<std::size_t,1> one{0};const std::array<std::size_t,2> two{0,1};
    auto radial=p_basis(true);orbital(radial,{.5,-1,0,0,1,0});
    output("radial-cancellation",cov::classify_local_metal_irrep(radial,one,0,"C2v",I));
    auto missing=radial;missing.ao_overlap.clear();
    output("missing-metric",cov::classify_local_metal_irrep(missing,one,0,"C2v",I));
    auto spin=p_basis();orbital(spin,{1,0,0});orbital(spin,{1,0,0},cov::Spin::Beta);
    output("alpha-beta-separate-rank",cov::classify_local_metal_irrep(spin,two,0,"C2v",I));
    spin.orbitals[1].spin=cov::Spin::Alpha;
    output("same-spin-duplicate-not-complete-MO-set",cov::classify_local_metal_irrep(spin,two,0,"C2v",I));
    cov::Wavefunction cart;cart.atoms.resize(1);cart.basis_count=6;cart.shells.push_back({0,0,0,0,2,0});
    cart.ao_overlap.assign(36,0);for(std::size_t i=0;i<6;++i)cart.ao_overlap[i*6+i]=1;
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)if(i!=j)cart.ao_overlap[i*6+j]=1.0/3.0;
    orbital(cart,{1,1,1,0,0,0});
    output("Cartesian-d-scalar-is-local-l0",cov::classify_local_metal_irrep(cart,one,0,"D4h",I));
    auto p=p_basis();orbital(p,{1,0,0});
    output("unknown-parent-group",cov::classify_local_metal_irrep(p,one,0,"missing",I));
    output("dimension-is-candidate-only",cov::classify_local_irrep_by_dimension("Oh",cov::MetalAOShell::D,2));
    output("same-dimension-different-labels-ambiguous",cov::classify_local_irrep_by_dimension("D4h",cov::MetalAOShell::D,1));
    cov::Wavefunction weak;weak.atoms.resize(1);weak.basis_count=8;
    weak.shells.push_back({0,0,0,0,1,0});weak.shells.push_back({0,0,0,3,2,1});
    weak.ao_overlap.assign(64,0);for(std::size_t i=0;i<8;++i)weak.ao_overlap[i*8+i]=1;
    orbital(weak,{std::sqrt(.5),std::sqrt(.5),0,std::sqrt(.03),0,0,0,0});
    output("weak-local-explanation-retains-absolute-coverage",cov::classify_local_metal_irrep(weak,one,0,"C2v",I));
    return 0;
}
