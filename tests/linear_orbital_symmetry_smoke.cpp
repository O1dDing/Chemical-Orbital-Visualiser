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

#include "cov/orbital_symmetry.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void identity_overlap(cov::Wavefunction& wf) {
    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;
}

void add_shell(cov::Wavefunction& wf,
               const std::size_t atom,
               const std::uint8_t l,
               const bool pure) {
    const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back({1.0f,1.0f});
    cov::Shell shell;
    shell.atom_index=static_cast<std::uint32_t>(atom);
    shell.primitive_offset=primitive_offset;
    shell.primitive_count=1;
    shell.basis_offset=wf.basis_count;
    shell.angular_momentum=l;
    shell.pure=pure?1u:0u;
    wf.shells.push_back(shell);
    wf.basis_count+=cov::shell_basis_count(shell);
}

void add_mo(cov::Wavefunction& wf,
            const double energy,
            const std::vector<std::pair<std::size_t,double>>& entries) {
    cov::MolecularOrbital mo;
    mo.energy_hartree=energy;
    mo.spin=cov::Spin::Alpha;
    mo.coefficients.assign(wf.basis_count,0.0f);
    for (const auto& [index,value]:entries) {
        mo.coefficients.at(index)=static_cast<float>(value);
    }
    wf.orbitals.push_back(std::move(mo));
}

bool check_label(const cov::Wavefunction& wf,
                 const std::size_t begin,
                 const std::size_t end,
                 const std::string& label) {
    for (std::size_t i=begin;i<end;++i) {
        if (wf.orbitals[i].symmetry!=label ||
            wf.orbitals[i].symmetry_provenance!=cov::DataProvenance::Derived) {
            std::cerr << "expected " << label << " at " << i
                      << ", got '" << wf.orbitals[i].symmetry << "'\n";
            return false;
        }
    }
    return true;
}

cov::Wavefunction make_linear() {
    cov::Wavefunction wf;
    wf.atoms={
        {"F",9,0.0,0.0,-2.0},
        {"F",9,0.0,0.0, 2.0},
    };

    // Two Cartesian p shells followed by two pure d shells.
    add_shell(wf,0,1,false); // basis 0..2 = px,py,pz
    add_shell(wf,1,1,false); // basis 3..5
    add_shell(wf,0,2,true);  // basis 6..10 = m0,+1,-1,+2,-2
    add_shell(wf,1,2,true);  // basis 11..15

    constexpr double s=0.70710678118654752440;

    // pz(A)-pz(B): inversion-even, sigma_v-even => Sigma_g+
    add_mo(wf,-1.0,{{2,s},{5,-s}});

    // px(A)+px(B), py(A)+py(B): inversion-odd, |m|=1 => Pi_u
    add_mo(wf,-0.4,{{0,s},{3,s}});
    add_mo(wf,-0.4,{{1,s},{4,s}});

    // Pure d m=+/-2 symmetric across centres: inversion-even => Delta_g
    add_mo(wf,0.1,{{9,s},{14,s}});
    add_mo(wf,0.1,{{10,s},{15,s}});

    identity_overlap(wf);
    return wf;
}

} // namespace

int main() {
    auto wf=make_linear();
    const auto result=cov::derive_orbital_symmetry(wf);
    if (result.point_group!="Dinfh" || result.orbitals_labelled!=5u ||
        !check_label(wf,0,1,"Σg+") ||
        !check_label(wf,1,3,"Πu") ||
        !check_label(wf,3,5,"Δg")) {
        std::cerr << "Dinfh Sigma/Pi/Delta regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "linear orbital symmetry smoke test passed\n";
    return EXIT_SUCCESS;
}
