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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void field_prefix(std::ofstream& out,const std::string& label,const char type) {
    out << std::left << std::setw(40) << label << type;
}

void scalar_i(std::ofstream& out,const std::string& label,const long long value) {
    field_prefix(out,label,'I');
    out << std::right << std::setw(16) << value << '\n';
}

void array_i(std::ofstream& out,const std::string& label,const std::vector<long long>& values) {
    field_prefix(out,label,'I');
    out << "   N=" << std::right << std::setw(12) << values.size() << '\n';
    for (std::size_t i=0;i<values.size();++i) {
        out << std::setw(12) << values[i];
        if ((i+1u)%6u==0u || i+1u==values.size()) out << '\n';
    }
}

void array_r(std::ofstream& out,const std::string& label,const std::vector<double>& values) {
    field_prefix(out,label,'R');
    out << "   N=" << std::right << std::setw(12) << values.size() << '\n';
    out << std::uppercase << std::scientific << std::setprecision(8);
    for (std::size_t i=0;i<values.size();++i) {
        out << std::setw(16) << values[i];
        if ((i+1u)%5u==0u || i+1u==values.size()) out << '\n';
    }
}

std::filesystem::path make_fchk() {
    const auto path=std::filesystem::temp_directory_path()/"cov_equivalent_h2.fchk";
    std::ofstream out(path,std::ios::binary);
    out << "Equivalent H2\nSP        RHF STO-3G\n";
    scalar_i(out,"Number of atoms",2);
    scalar_i(out,"Number of alpha electrons",1);
    scalar_i(out,"Number of beta electrons",1);
    scalar_i(out,"Number of basis functions",2);
    scalar_i(out,"Number of independent functions",2);
    array_i(out,"Atomic numbers",{1,1});
    array_r(out,"Current cartesian coordinates",{0.0,0.0,-0.7,0.0,0.0,0.7});
    array_i(out,"Shell types",{0,0});
    array_i(out,"Number of primitives per shell",{1,1});
    array_i(out,"Shell to atom map",{1,2});
    array_r(out,"Primitive exponents",{1.0,1.0});
    array_r(out,"Contraction coefficients",{1.0,1.0});
    array_r(out,"Alpha Orbital Energies",{-0.5,0.2});
    array_r(out,"Alpha MO coefficients",{
        0.7071067811865475,0.7071067811865475,
        0.7071067811865475,-0.7071067811865475
    });
    return path;
}

std::filesystem::path make_molden() {
    const auto path=std::filesystem::temp_directory_path()/"cov_equivalent_h2.molden";
    std::ofstream out(path,std::ios::binary);
    out << "[Molden Format]\n"
           "[Atoms] AU\n"
           "H 1 1 0.0 0.0 -0.7\n"
           "H 2 1 0.0 0.0  0.7\n"
           "[GTO]\n"
           "1 0\n"
           "s 1 1.0\n"
           " 1.0 1.0\n\n"
           "2 0\n"
           "s 1 1.0\n"
           " 1.0 1.0\n\n"
           "[MO]\n"
           "Ene= -0.5\n"
           "Spin= Alpha\n"
           "Occup= 2.0\n"
           "1 0.7071067811865475\n"
           "2 0.7071067811865475\n"
           "Ene= 0.2\n"
           "Spin= Alpha\n"
           "Occup= 0.0\n"
           "1 0.7071067811865475\n"
           "2 -0.7071067811865475\n";
    return path;
}

bool close(double a,double b,double tol=2.0e-6) {
    return std::abs(a-b)<=tol;
}

bool equivalent(const cov::Wavefunction& a,const cov::Wavefunction& b) {
    if (a.atoms.size()!=b.atoms.size() || a.shells.size()!=b.shells.size() ||
        a.basis_count!=b.basis_count || a.orbitals.size()!=b.orbitals.size()) return false;
    for (std::size_t i=0;i<a.atoms.size();++i) {
        if (a.atoms[i].atomic_number!=b.atoms[i].atomic_number ||
            !close(a.atoms[i].x,b.atoms[i].x) || !close(a.atoms[i].y,b.atoms[i].y) ||
            !close(a.atoms[i].z,b.atoms[i].z)) return false;
    }
    for (std::size_t i=0;i<a.shells.size();++i) {
        if (a.shells[i].atom_index!=b.shells[i].atom_index ||
            a.shells[i].angular_momentum!=b.shells[i].angular_momentum ||
            a.shells[i].pure!=b.shells[i].pure ||
            a.shells[i].primitive_count!=b.shells[i].primitive_count) return false;
    }
    for (std::size_t i=0;i<a.orbitals.size();++i) {
        const auto& x=a.orbitals[i];
        const auto& y=b.orbitals[i];
        if (!close(x.energy_hartree,y.energy_hartree) || !close(x.occupation,y.occupation) ||
            x.spin!=y.spin || x.coefficients.size()!=y.coefficients.size()) return false;
        for (std::size_t j=0;j<x.coefficients.size();++j) {
            if (!close(x.coefficients[j],y.coefficients[j])) return false;
        }
    }
    return true;
}

} // namespace

int main() {
    const auto fchk=make_fchk();
    const auto molden=make_molden();
    std::error_code ec;
    try {
        const auto a=cov::parse_wavefunction(fchk);
        const auto b=cov::parse_wavefunction(molden);
        const bool ok=a.source==cov::WavefunctionSource::Fchk &&
                      b.source==cov::WavefunctionSource::Molden && equivalent(a,b);
        std::filesystem::remove(fchk,ec);
        std::filesystem::remove(molden,ec);
        if (!ok) {
            std::cerr << "FCHK/Molden semantic equivalence regression\n";
            return EXIT_FAILURE;
        }
        std::cout << "FCHK/Molden semantic equivalence smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::filesystem::remove(fchk,ec);
        std::filesystem::remove(molden,ec);
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
