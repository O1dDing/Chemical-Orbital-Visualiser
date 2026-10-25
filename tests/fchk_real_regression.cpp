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

#include "cov/wavefunction_io.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

namespace {

bool close(double a, double b, double eps = 1.0e-6) {
    return std::abs(a - b) <= eps;
}

} // namespace

int main(int argc, char** argv) {
    // The real external fixture is downloaded by CI from a pinned cclib commit.
    // A normal offline/local ctest run should remain usable without network.
    if (argc < 2 || !std::filesystem::exists(argv[1])) {
        std::cout << "real Gaussian FCHK fixture unavailable; skipping\n";
        return 77;
    }

    cov::WavefunctionParseOptions options;
    options.max_atoms = 100;
    options.require_orbitals = true;
    options.keep_density = true;
    options.reconstruct_density_if_missing = true;
    options.auto_enrich_gaussian_log = false;

    // Exercise the same FCHK-first dispatch/post-processing path used by the
    // viewer, not only the low-level parser.
    const cov::Wavefunction wf = cov::parse_wavefunction(argv[1], options);

    // Pinned fixture:
    // cclib/cclib@21daa960123d28aa21eeaaacc2a1dea39e136829
    // data/FChk/basicGaussian16/dvb_sp.fchk
    // Gaussian 16, restricted B3LYP/STO-3G, 20 atoms, 60 basis functions,
    // 35 alpha + 35 beta electrons.
    if (wf.source != cov::WavefunctionSource::Fchk ||
        wf.atoms.size() != 20u ||
        wf.basis_count != 60u ||
        wf.orbitals.size() != 60u ||
        wf.alpha_electrons != 35u ||
        wf.beta_electrons != 35u) {
        std::cerr << "real Gaussian16 FCHK headline metadata regression\n";
        return EXIT_FAILURE;
    }

    for (std::size_t i = 0; i < wf.orbitals.size(); ++i) {
        const auto& mo = wf.orbitals[i];
        if (mo.spin != cov::Spin::Alpha || mo.coefficients.size() != 60u) {
            std::cerr << "real Gaussian16 restricted MO block regression at " << i << '\n';
            return EXIT_FAILURE;
        }
        const double expected_occ = i < 35u ? 2.0 : 0.0;
        if (!close(mo.occupation, expected_occ)) {
            std::cerr << "real Gaussian16 occupation reconstruction regression at " << i << '\n';
            return EXIT_FAILURE;
        }
        if (!std::isfinite(mo.energy_hartree)) {
            std::cerr << "real Gaussian16 non-finite orbital energy at " << i << '\n';
            return EXIT_FAILURE;
        }
    }

    const std::size_t packed = 60u * 61u / 2u;
    if (wf.total_density_packed.size() != packed ||
        wf.total_density_provenance == cov::DataProvenance::Unavailable) {
        std::cerr << "real Gaussian16 density retention/reconstruction regression\n";
        return EXIT_FAILURE;
    }

    if (wf.ao_overlap.size() != 60u * 60u ||
        wf.ao_overlap_provenance == cov::DataProvenance::Unavailable) {
        std::cerr << "real Gaussian16 AO overlap import/reconstruction regression\n";
        return EXIT_FAILURE;
    }

    std::size_t shell_basis_total = 0;
    for (const auto& shell : wf.shells) shell_basis_total += cov::shell_basis_count(shell);
    if (shell_basis_total != wf.basis_count) {
        std::cerr << "real Gaussian16 shell/basis dimension regression\n";
        return EXIT_FAILURE;
    }

    if (wf.point_group_detected.empty() ||
        wf.point_group_provenance == cov::DataProvenance::Unavailable) {
        std::cerr << "real Gaussian16 geometry point-group regression\n";
        return EXIT_FAILURE;
    }

    if (wf.bond_order_provenance != cov::DataProvenance::Derived) {
        std::cerr << "real Gaussian16 electronic bond-analysis regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "real Gaussian16 FCHK full-pipeline regression passed\n";
    return EXIT_SUCCESS;
}
