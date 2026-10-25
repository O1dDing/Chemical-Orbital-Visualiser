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

#include "cov/local_angular_projection.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

void number(double value) { if(std::isfinite(value))std::cout<<value;else std::cout<<"null"; }
int main() {
    std::cout<<std::setprecision(17);
    int id=0,atoms=0,n=0,ns=0,nmo=0,selected=0,atom=0,mutate=0;
    while(std::cin>>id>>atoms>>n>>ns>>nmo>>selected>>atom>>mutate) {
        cov::Wavefunction wf;wf.atoms.resize(atoms);wf.basis_count=n;
        std::array<double,9> rotation{};for(auto& v:rotation)if(!(std::cin>>v))return 2;
        for(int i=0;i<ns;++i) {
            int owner=0,offset=0,L=0,pure=0;std::cin>>owner>>offset>>L>>pure;
            cov::Shell shell; shell.atom_index=owner;shell.basis_offset=offset;
            shell.angular_momentum=static_cast<std::uint8_t>(L);shell.pure=static_cast<std::uint8_t>(pure);wf.shells.push_back(shell);
        }
        wf.ao_overlap.resize(static_cast<std::size_t>(n)*n);for(auto& v:wf.ao_overlap)std::cin>>v;
        for(int j=0;j<nmo;++j) {
            int spin=0;std::cin>>spin;cov::MolecularOrbital mo;mo.spin=static_cast<cov::Spin>(spin);
            mo.coefficients.resize(n);for(auto& v:mo.coefficients)std::cin>>v;wf.orbitals.push_back(std::move(mo));
        }
        std::vector<std::size_t> indices(selected);for(auto& v:indices)std::cin>>v;
        cov::LocalAngularProjectionWorkspace workspace(wf,atom,rotation);
        if(mutate) {wf.ao_overlap.clear();wf.shells.clear();wf.orbitals.clear();wf.atoms.clear();}
        const auto r=workspace.project(indices);
        std::cout<<"{\"id\":"<<id<<",\"status\":"<<static_cast<int>(r.status)<<",\"rank\":"<<r.represented_spin_orbital_rank
            <<",\"trace\":";number(r.centre_projection_trace);std::cout<<",\"fraction\":";number(r.centre_mean_fraction);
        std::cout<<",\"partition_residual\":";number(r.angular_partition_residual);std::cout<<",\"components\":[";
        bool first=true;for(double v:r.component_projection_traces){if(!first)std::cout<<',';first=false;number(v);}
        std::cout<<"],\"spins\":[";first=true;
        for(const auto& b:r.spins) {
            if(!first)std::cout<<',';first=false;
            std::cout<<"{\"spin\":"<<static_cast<int>(b.spin)<<",\"rank\":"<<b.centre.orbitals.numerical_rank
                <<",\"reference_rank\":"<<b.centre.reference.numerical_rank<<",\"trace\":";number(b.centre.subspace_overlap_trace);
            std::cout<<",\"components\":[";
            bool start=true;for(const auto& c:b.components) {
                if(!start)std::cout<<',';start=false;
                std::cout<<"{\"l\":"<<c.angular_degree<<",\"component\":"<<c.component<<",\"reference_rank\":"<<c.metric.reference.numerical_rank
                    <<",\"trace\":";number(c.metric.subspace_overlap_trace);std::cout<<",\"fraction_local\":";number(c.fraction_of_local_projection);
                std::cout<<",\"reference_residual\":";number(c.metric.reference.orthonormality_residual);
                std::cout<<",\"orbital_residual\":";number(c.metric.orbitals.orthonormality_residual);std::cout<<'}';
            }
            std::cout<<"]}";
        }
        std::cout<<"]}\n";
    }
    return std::cin.eof()?0:2;
}
