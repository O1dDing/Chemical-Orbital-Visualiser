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

#include "cov/molden_parser.hpp"
#include "cov/wavefunction_io.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

namespace {

bool test_no_sym_mo_boundaries() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_no_sym_regression.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) {
            std::cerr << "unable to create temporary no-Sym Molden fixture\n";
            return false;
        }
        out << "[Molden Format]\n"
               "[Atoms] Angs\n"
               "H 1 1 0.0 0.0 0.0\n"
               "H 2 1 0.0 0.0 0.74\n"
               "[GTO]\n"
               "1 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "2 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "[MO]\n"
               "Ene= -0.500000\n"
               "Spin= Alpha\n"
               "Occup= 2.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n"
               "Ene= 0.200000\n"
               "Spin= Alpha\n"
               "Occup= 0.0\n"
               "1 0.70710678\n"
               "2 -0.70710678\n";
    }

    try {
        const auto wf = cov::parse_molden(path);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (wf.orbitals.size() != 2 || wf.basis_count != 2) {
            std::cerr << "no-Sym MO boundary regression: expected 2 MOs / 2 basis, got "
                      << wf.orbitals.size() << " / " << wf.basis_count << '\n';
            return false;
        }
        if (wf.orbitals[0].occupation != 2.0f || wf.orbitals[1].occupation != 0.0f) {
            std::cerr << "no-Sym MO metadata regression\n";
            return false;
        }
        if (wf.alpha_electrons != 1u || wf.beta_electrons != 1u ||
            wf.electron_counts_provenance != cov::DataProvenance::Derived) {
            std::cerr << "restricted Molden electron-count provenance regression\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::cerr << "no-Sym MO parsing failed: " << e.what() << '\n';
        return false;
    }
}

bool test_unrestricted_electron_count_derivation() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_unrestricted_count_regression.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) return false;
        out << "[Molden Format]\n"
               "[Atoms] Angs\n"
               "H 1 1 0.0 0.0 0.0\n"
               "H 2 1 0.0 0.0 0.74\n"
               "[GTO]\n"
               "1 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "2 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "[MO]\n"
               "Ene= -0.500000\n"
               "Spin= Alpha\n"
               "Occup= 1.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n"
               "Ene= -0.480000\n"
               "Spin= Beta\n"
               "Occup= 1.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n";
    }
    try {
        const auto wf = cov::parse_wavefunction(path);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (wf.alpha_electrons != 1u || wf.beta_electrons != 1u ||
            wf.electron_counts_provenance != cov::DataProvenance::Derived) {
            std::cerr << "unrestricted Molden electron-count provenance regression\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::cerr << "unrestricted Molden count derivation failed: " << e.what() << '\n';
        return false;
    }
}

bool test_unrestricted_double_occupation_stays_unavailable() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_unrestricted_double_occupation.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) return false;
        out << "[Molden Format]\n"
               "[Atoms] Angs\n"
               "H 1 1 0.0 0.0 0.0\n"
               "H 2 1 0.0 0.0 0.74\n"
               "[GTO]\n"
               "1 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "2 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "[MO]\n"
               "Ene= -0.500000\n"
               "Spin= Alpha\n"
               "Occup= 2.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n"
               "Ene= -0.480000\n"
               "Spin= Beta\n"
               "Occup= 0.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n";
    }
    try {
        const auto wf=cov::parse_wavefunction(path);
        std::error_code ec;
        std::filesystem::remove(path,ec);
        if (wf.electron_counts_provenance!=cov::DataProvenance::Unavailable) {
            std::cerr<<"invalid explicit-spin Occup=2 produced electron counts\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path,ec);
        std::cerr<<"explicit-spin occupation validity test failed: "<<e.what()<<'\n';
        return false;
    }
}

bool test_missing_occupation_stays_unavailable() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_missing_occupation_regression.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) return false;
        out << "[Molden Format]\n"
               "[Atoms] Angs\n"
               "H 1 1 0.0 0.0 0.0\n"
               "H 2 1 0.0 0.0 0.74\n"
               "[GTO]\n"
               "1 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "2 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "[MO]\n"
               "Ene= -0.500000\n"
               "Spin= Alpha\n"
               "Occup= 2.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n"
               "Ene= 0.200000\n"
               "Spin= Alpha\n"
               "1 0.70710678\n"
               "2 -0.70710678\n";
    }
    try {
        const auto wf = cov::parse_molden(path);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (wf.orbitals.size()!=2u ||
            wf.orbitals[0].occupation_provenance!=cov::DataProvenance::Producer ||
            wf.orbitals[1].occupation_provenance!=cov::DataProvenance::Unavailable ||
            wf.electron_counts_provenance!=cov::DataProvenance::Unavailable) {
            std::cerr << "missing Molden Occup was promoted to derived electron counts\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::cerr << "missing Molden occupation test failed: " << e.what() << '\n';
        return false;
    }
}

bool test_omitted_7f_marker_inference() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_omitted_7f_regression.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) {
            std::cerr << "unable to create temporary omitted-7F Molden fixture\n";
            return false;
        }
        out << "[Molden Format]\n"
               "[Atoms] AU\n"
               "C 1 6 0.0 0.0 0.0\n"
               "[GTO]\n"
               "1 0\n"
               "d 1 1.0\n"
               "1.0 1.0\n"
               "f 1 1.0\n"
               "1.0 1.0\n"
               "[5D]\n"
               "[MO]\n"
               "Ene= -0.100000\n"
               "Spin= Alpha\n"
               "Occup= 2.0\n";
        for (int i = 1; i <= 12; ++i) {
            out << i << " " << (i == 1 ? 1.0 : 0.0) << "\n";
        }
    }

    try {
        const auto wf = cov::parse_molden(path);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (wf.basis_count != 12u || !wf.pure_d || !wf.pure_f) {
            std::cerr << "omitted-7F regression: expected 5D+7F = 12 basis, got basis="
                      << wf.basis_count << " pure_d=" << wf.pure_d
                      << " pure_f=" << wf.pure_f << '\n';
            return false;
        }
        if (wf.orbitals.size() != 1u || wf.orbitals[0].coefficients.size() != 12u) {
            std::cerr << "omitted-7F regression: MO dimension mismatch\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::cerr << "omitted-7F inference failed: " << e.what() << '\n';
        return false;
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: cov_parser_smoke <file.molden> [--inspect]\n";
        return EXIT_FAILURE;
    }

    try {
        const auto wf = cov::parse_molden(argv[1]);
        std::cout << "atoms=" << wf.atoms.size()
                  << " shells=" << wf.shells.size()
                  << " basis=" << wf.basis_count
                  << " orbitals=" << wf.orbitals.size() << '\n';
        if (argc == 3 && std::string_view(argv[2]) == "--inspect") {
            return EXIT_SUCCESS;
        }
        if (wf.atoms.size() != 2) {
            std::cerr << "expected 2 atoms, got " << wf.atoms.size() << '\n';
            return EXIT_FAILURE;
        }
        if (wf.basis_count != 2) {
            std::cerr << "expected 2 basis functions, got " << wf.basis_count << '\n';
            return EXIT_FAILURE;
        }
        if (wf.orbitals.size() != 2) {
            std::cerr << "expected 2 orbitals, got " << wf.orbitals.size() << '\n';
            return EXIT_FAILURE;
        }
        if (wf.orbitals[0].coefficients.size() != wf.basis_count) {
            std::cerr << "MO coefficient consistency check failed\n";
            return EXIT_FAILURE;
        }
        if (!test_no_sym_mo_boundaries()) {
            return EXIT_FAILURE;
        }
        if (!test_omitted_7f_marker_inference()) {
            return EXIT_FAILURE;
        }
        if (!test_unrestricted_electron_count_derivation()) {
            return EXIT_FAILURE;
        }
        if (!test_missing_occupation_stays_unavailable()) {
            return EXIT_FAILURE;
        }
        if (!test_unrestricted_double_occupation_stays_unavailable()) {
            return EXIT_FAILURE;
        }
        std::cout << "parser smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
