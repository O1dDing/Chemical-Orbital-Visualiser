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

#include "cov/overlap.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(const bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
void close(const double a,const double b,const double tolerance=2.0e-12) {
    if (!std::isfinite(a) || std::abs(a-b)>tolerance)
        throw std::runtime_error("Analytic AO integral differs from independently known value");
}
cov::Wavefunction one_shell(const int l,const bool pure) {
    cov::Wavefunction wf;
    wf.atoms.push_back({"H",1,0,0,0});
    wf.primitives.push_back({0.7312345678912345,1.0});
    cov::Shell shell;
    shell.primitive_count=1;shell.angular_momentum=static_cast<std::uint8_t>(l);shell.pure=pure?1:0;
    wf.shells.push_back(shell);wf.basis_count=cov::shell_basis_count(shell);
    return wf;
}
cov::Wavefunction h2() {
    auto wf=one_shell(0,false);
    wf.atoms[0].z=-0.7;
    wf.atoms.push_back({"H",1,0,0,0.7});
    wf.primitives[0].exponent=1;
    wf.primitives.push_back({1,1});
    auto second=wf.shells[0];second.atom_index=1;second.basis_offset=1;second.primitive_offset=1;
    wf.shells.push_back(second);wf.basis_count=2;
    const double s=std::exp(-0.98);
    cov::MolecularOrbital bonding;
    bonding.coefficients={1/std::sqrt(2*(1+s)),1/std::sqrt(2*(1+s))};
    wf.orbitals.push_back(bonding);
    return wf;
}
}

int main() {
    try {
        // Primitive norms are exactly one; real spherical components form an
        // orthonormal set. Cartesian monomials are not mutually orthogonal.
        for (int l=0;l<=4;++l) for (const bool pure:{false,true}) {
            const auto wf=one_shell(l,pure);
            const auto s=cov::derive_ao_overlap_from_basis(wf);
            require(s.available(),"Supported s--g shell integral unavailable");
            const auto n=wf.basis_count;
            for (std::size_t i=0;i<n;++i) {
                close(s.matrix[i*n+i],1.0);
                for (std::size_t j=0;j<n;++j) {
                    close(s.matrix[i*n+j],s.matrix[j*n+i]);
                    if (pure && i!=j) close(s.matrix[i*n+j],0.0);
                }
            }
        }
        const auto d=cov::derive_ao_overlap_from_basis(one_shell(2,false));
        close(d.matrix[1],1.0/3.0); // <x^2|y^2> for normalized monomials
        const auto g=cov::derive_ao_overlap_from_basis(one_shell(4,false));
        close(g.matrix[1],3.0/35.0); // <x^4|y^4>

        // Independent, explicit d-shell pure-to-Cartesian polynomial map.
        auto mixed=one_shell(2,false);
        auto pure=mixed.shells[0];pure.pure=1;pure.basis_offset=6;
        mixed.shells.push_back(pure);mixed.basis_count=11;
        const auto both=cov::derive_ao_overlap_from_basis(mixed);
        require(both.available(),"Mixed shell overlap unavailable");
        const double t[5][6]={{-.5,-.5,1,0,0,0},{0,0,0,0,-1,0},
            {0,0,0,0,0,-1},{std::sqrt(3.0)/2,-std::sqrt(3.0)/2,0,0,0,0},{0,0,0,1,0,0}};
        for (std::size_t i=0;i<5;++i) for (std::size_t j=0;j<6;++j) {
            double expected=0;
            for (std::size_t k=0;k<6;++k) expected+=t[i][k]*d.matrix[k*6+j];
            close(both.matrix[(i+6)*11+j],expected);
        }

        auto wf=h2();
        cov::establish_ao_metric(wf);
        close(wf.ao_overlap[1],std::exp(-0.98));
        require(cov::orbital_metric_usable(wf),"Rectangular 2 AO / 1 MO block rejected");
        require(wf.ao_metric_diagnostics.numerical_rank==2,"Full AO rank lost with truncated MO input");
        close(wf.ao_overlap_orthonormality_error,0.0);
        require(std::isnan(wf.ao_metric_diagnostics.orbital_blocks[0].maximum_off_diagonal_error),
                "A single MO has no off-diagonal Gram entries to measure");
        const auto unshifted=wf.ao_overlap;
        for (auto& atom:wf.atoms) {atom.x+=1024;atom.y-=64;atom.z+=16;}
        close(cov::derive_ao_overlap_from_basis(wf).matrix[1],unshifted[1]);
        wf.orbitals.clear();
        cov::establish_ao_metric(wf);
        require(wf.ao_overlap.size()==4 && wf.ao_metric_diagnostics.status==cov::NumericalStatus::Available,
                "AO integrals incorrectly require MO coefficients");
        require(std::isnan(wf.ao_overlap_orthonormality_error),"Missing MO block became a zero residual");

        for (const double reported:{std::exp(-.98),.2,1.2,std::numeric_limits<double>::quiet_NaN()}) {
            wf=h2();wf.producer_ao_overlap={1,reported,reported,1};
            const auto original=wf.orbitals[0].coefficients;
            cov::establish_ao_metric(wf);
            close(wf.ao_overlap[1],std::exp(-.98));
            require(wf.orbitals[0].coefficients==original,"Input coefficients were altered to fit S");
            require(cov::orbital_metric_usable(wf),"A bad optional producer S suppressed valid independent basis/MO data");
            if (reported==std::exp(-.98)) {
                require(wf.producer_ao_overlap_basis_status==cov::NumericalStatus::Available,"Consistent producer S rejected");
                close(wf.producer_ao_overlap_basis_error,0);
            } else {
                require(wf.producer_ao_overlap_basis_status==cov::NumericalStatus::InvalidInput,"Contradictory producer S accepted");
                if (!std::isfinite(reported))
                    require(std::isnan(wf.producer_ao_overlap_basis_error),"Nonfinite comparison became a zero residual");
            }
        }

        wf=h2();
        auto beta=wf.orbitals[0];beta.spin=cov::Spin::Beta;
        wf.orbitals.push_back(beta);
        cov::establish_ao_metric(wf);
        require(cov::orbital_metric_usable(wf),"Separate rectangular alpha and beta blocks rejected");
        wf.orbitals.back().coefficients[0]*=1.2;
        cov::establish_ao_metric(wf);
        require(!cov::orbital_metric_usable(wf) && wf.ao_overlap_orthonormality_error>0.01,
                "Invalid beta normalization was not detected");
        wf.orbitals.back().coefficients[0]=std::numeric_limits<double>::quiet_NaN();
        cov::establish_ao_metric(wf);
        require(wf.ao_metric_diagnostics.orbital_blocks[1].status==cov::NumericalStatus::InvalidInput,
                "Nonfinite MO coefficient was not detected");

        wf=h2();wf.primitives[0].exponent=0;
        require(cov::derive_ao_overlap_from_basis(wf).status==cov::NumericalStatus::InvalidInput,
                "Nonpositive Gaussian exponent accepted");
        wf=h2();wf.shells[1].basis_offset=0;
        require(cov::derive_ao_overlap_from_basis(wf).status==cov::NumericalStatus::InvalidInput,
                "Overlapping AO ranges accepted");
        wf=h2();
        require(cov::inspect_ao_metric(wf,{1,1.2,1.2,1}).status==cov::NumericalStatus::InvalidInput,
                "Indefinite metric accepted");
        require(cov::inspect_ao_metric(wf,{1,.2,.3,1}).status==cov::NumericalStatus::InvalidInput,
                "Nonsymmetric metric accepted");
        const auto singular=cov::inspect_ao_metric(wf,{1,1,1,1});
        require(singular.status==cov::NumericalStatus::Available && singular.numerical_rank==1,
                "Semidefinite metric rank not reported");
        require(std::isnan(cov::AoMetricDiagnostics{}.eigen_residual),"Unmeasured spectral residual became zero");
        std::cout << "Independent basis integrals, rectangular spin metrics and invalid-input diagnostics passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';return 1;
    }
}
