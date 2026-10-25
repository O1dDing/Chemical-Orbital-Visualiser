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
#include "cov/fchk_parser.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t large_basis=384u;

cov::Wavefunction make_dense_wavefunction(const bool singular) {
    cov::Wavefunction wavefunction;
    wavefunction.basis_count=large_basis;
    wavefunction.orbitals.reserve(large_basis);
    for (std::size_t orbital=0u;orbital<large_basis;++orbital) {
        cov::MolecularOrbital mo;
        mo.spin=cov::Spin::Alpha;
        mo.coefficients.assign(large_basis,0.0f);
        mo.coefficients[orbital]=static_cast<float>(
            1.0+2.0e-4*static_cast<double>(orbital%17u));
        mo.coefficients[(orbital+1u)%large_basis]=0.04f;
        mo.coefficients[(orbital+7u)%large_basis]=-0.012f;
        wavefunction.orbitals.push_back(std::move(mo));
    }
    if (singular) {
        wavefunction.orbitals[1u].coefficients=
            wavefunction.orbitals[0u].coefficients;
    }
    return wavefunction;
}

double independent_complete_metric_error(
    const cov::Wavefunction& wavefunction,
    const std::vector<double>& overlap) {
    const std::size_t n=wavefunction.basis_count;
    std::vector<double> coefficients(n*n,0.0);
    for (std::size_t orbital=0u;orbital<n;++orbital) {
        for (std::size_t basis=0u;basis<n;++basis) {
            coefficients[basis*n+orbital]=static_cast<double>(
                wavefunction.orbitals[orbital].coefficients[basis]);
        }
    }
    std::vector<double> sc(n*n,0.0);
    for (std::size_t row=0u;row<n;++row) {
        for (std::size_t k=0u;k<n;++k) {
            const double scale=overlap[row*n+k];
            for (std::size_t column=0u;column<n;++column) {
                sc[row*n+column]+=scale*coefficients[k*n+column];
            }
        }
    }
    double maximum_error=0.0;
    for (std::size_t row=0u;row<n;++row) {
        for (std::size_t column=0u;column<n;++column) {
            double value=0.0;
            for (std::size_t k=0u;k<n;++k) {
                value+=coefficients[k*n+row]*sc[k*n+column];
            }
            const double target=row==column?1.0:0.0;
            maximum_error=std::max(
                maximum_error,std::abs(value-target));
        }
    }
    return maximum_error;
}

} // namespace

int main(const int argc,char** argv) {
    if (argc>1) {
        cov::FchkParseOptions options;
        options.max_atoms=256u;
        options.keep_density=false;
        options.reconstruct_density_if_missing=false;
        const auto wavefunction=cov::parse_fchk(argv[1],options);
        const auto start=std::chrono::steady_clock::now();
        const auto overlap=cov::derive_ao_overlap_from_mos(wavefunction);
        const double seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
        if (!overlap.available()) {
            std::cerr << "real FCHK overlap recovery failed: n="
                      << wavefunction.basis_count << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "real FCHK overlap recovery passed: n="
                  << wavefunction.basis_count << " seconds=" << seconds
                  << " max_error=" << overlap.max_orthonormality_error << '\n';
        return EXIT_SUCCESS;
    }

    const auto wavefunction=make_dense_wavefunction(false);
    const auto start=std::chrono::steady_clock::now();
    const auto overlap=cov::derive_ao_overlap_from_mos(wavefunction);
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    if (!overlap.available() ||
        overlap.matrix.size()!=large_basis*large_basis) {
        std::cerr << "large reversible overlap recovery failed\n";
        return EXIT_FAILURE;
    }
    double symmetry_error=0.0;
    for (std::size_t i=0u;i<large_basis;++i) {
        for (std::size_t j=0u;j<large_basis;++j) {
            symmetry_error=std::max(
                symmetry_error,
                std::abs(overlap.matrix[i*large_basis+j]-
                         overlap.matrix[j*large_basis+i]));
        }
    }
    const double independent_error=independent_complete_metric_error(
        wavefunction,overlap.matrix);
    if (symmetry_error>1.0e-10 || independent_error>5.0e-6 ||
        std::abs(independent_error-overlap.max_orthonormality_error)>
            1.0e-9) {
        std::cerr << "large complete C^T*S*C validation regression: symmetry="
                  << symmetry_error << " independent=" << independent_error
                  << " reported=" << overlap.max_orthonormality_error << '\n';
        return EXIT_FAILURE;
    }

    const auto singular=make_dense_wavefunction(true);
    const auto singular_result=cov::derive_ao_overlap_from_mos(singular);
    if (singular_result.available()) {
        std::cerr << "large singular coefficient matrix was accepted\n";
        return EXIT_FAILURE;
    }

    std::cout << "large overlap smoke test passed: n=" << large_basis
              << " seconds=" << seconds
              << " max_error=" << overlap.max_orthonormality_error << '\n';
    return EXIT_SUCCESS;
}
