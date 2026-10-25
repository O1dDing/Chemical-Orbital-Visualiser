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

#include "cov/fchk_parser.hpp"
#include "cov/fchk_overlap.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {
void scalar(std::ostream& out,const std::string& name,int value) {
    out<<std::left<<std::setw(40)<<name<<"I"<<std::right<<std::setw(16)<<value<<'\n';
}
template<class T>
void array(std::ostream& out,const std::string& name,char kind,const std::vector<T>& values) {
    out<<std::left<<std::setw(40)<<name<<kind<<"   N="<<std::right<<std::setw(12)<<values.size()<<'\n';
    const std::size_t per_line=kind=='I'?6u:5u;
    for (std::size_t i=0;i<values.size();++i) {
        out<<std::scientific<<std::setprecision(8)<<std::setw(16)<<values[i];
        if ((i+1)%per_line==0 || i+1==values.size()) out<<'\n';
    }
}
std::size_t packed(std::size_t i,std::size_t j) {
    if (i<j) std::swap(i,j);
    return i*(i+1)/2+j;
}
bool close(double a,double b) {return std::abs(a-b)<1e-10;}
}

int main() {
    // These are transport controls, not normalized SCF reference molecules.
    // Independent expected phase lists cover every pure d/f/g component.
    const std::vector<int> shell_types{0,1,-1,2,-2,3,-3,4,-4};
    const std::vector<std::vector<double>> phases{
        {1},{1,1,1},{1,1,1,1},{1,1,1,1,1,1},{1,-1,-1,1,1},
        {1,1,1,1,1,1,1,1,1,1},{1,-1,-1,1,1,-1,-1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{1,-1,-1,1,1,-1,-1,1,1}};
    const std::vector<std::size_t> cartesian_g{14,4,0,13,12,8,3,5,1,11,9,2,10,7,6};
    const auto directory=std::filesystem::temp_directory_path()/"cov-fchk-ao-contract";
    std::filesystem::create_directories(directory);
    try {
        for (std::size_t test=0;test<shell_types.size();++test) {
            const auto n=phases[test].size();
            const auto path=directory/(std::to_string(test)+".fch");
            std::vector<double> alpha(n),beta(n),density(n*(n+1)/2),overlap(n*(n+1)/2);
            for (std::size_t i=0;i<n;++i) {
                alpha[i]=static_cast<double>(i+1)+0.123456789;
                beta[i]=-static_cast<double>(i+1)-0.314159265;
                for (std::size_t j=0;j<=i;++j) {
                    density[packed(i,j)]=static_cast<double>((i+1)*(j+1));
                    overlap[packed(i,j)]=i==j?1.0:0.01*static_cast<double>(i+j+1);
                }
            }
            // Match the decimal values written by the FCHK fixture formatter.
            alpha[0]=0.123456789;
            beta[0]=-0.314159265;
            std::ofstream out(path,std::ios::binary);
            out<<"Gaussian AO representation transport control\nSP        UHF synthetic\n";
            scalar(out,"Number of atoms",1);scalar(out,"Number of basis functions",static_cast<int>(n));
            scalar(out,"Number of independent functions",1);
            scalar(out,"Number of alpha electrons",1);scalar(out,"Number of beta electrons",1);
            array(out,"Atomic numbers",'I',std::vector<int>{2});
            array(out,"Current cartesian coordinates",'R',std::vector<double>{0,0,0});
            array(out,"Shell types",'I',std::vector<int>{shell_types[test]});
            array(out,"Number of primitives per shell",'I',std::vector<int>{1});
            array(out,"Shell to atom map",'I',std::vector<int>{1});
            array(out,"Primitive exponents",'R',std::vector<double>{0.123456789});
            array(out,"Contraction coefficients",'R',std::vector<double>{0.314159265});
            if (shell_types[test]==-1) array(out,"P(S=P) Contraction coefficients",'R',std::vector<double>{0.271828182});
            array(out,"Alpha Orbital Energies",'R',std::vector<double>{-0.5});
            array(out,"Alpha MO coefficients",'R',alpha);
            array(out,"Beta Orbital Energies",'R',std::vector<double>{-0.4});
            array(out,"Beta MO coefficients",'R',beta);
            array(out,"Total SCF Density",'R',density);array(out,"Spin SCF Density",'R',density);
            array(out,"Overlap Matrix",'R',overlap);out.close();
            auto wf=cov::parse_fchk(path);
            if (!cov::enrich_fchk_overlap_from_file(wf,path)) throw std::runtime_error("missing overlap transport");
            if (wf.primitives[0].exponent!=0.123456789 || wf.primitives[0].coefficient!=0.314159265 ||
                wf.orbitals[0].gaussian_source_coefficients[0]!=0.123456789 ||
                wf.orbitals[1].gaussian_source_coefficients[0]!=-0.314159265) {
                throw std::runtime_error("input precision lost before analysis");
            }
            if (wf.gaussian_ao_transform.size()!=n) throw std::runtime_error("missing shared AO convention");
            for (std::size_t i=0;i<n;++i) {
                const auto si=shell_types[test]==4?cartesian_g[i]:i;
                for (std::size_t spin=0;spin<2;++spin) {
                    const auto& mo=wf.orbitals[spin];
                    if (!close(mo.coefficients[i],phases[test][i]*mo.gaussian_source_coefficients[si]))
                        throw std::runtime_error("alpha/beta coefficient phase/order mismatch");
                }
                for (std::size_t j=0;j<=i;++j) {
                    const auto sj=shell_types[test]==4?cartesian_g[j]:j;
                    const double phase=phases[test][i]*phases[test][j];
                    if (!close(wf.total_density_packed[packed(i,j)],phase*density[packed(si,sj)]) ||
                        !close(wf.spin_density_packed[packed(i,j)],phase*density[packed(si,sj)]) ||
                        !close(wf.ao_overlap[i*n+j],phase*overlap[packed(si,sj)]))
                        throw std::runtime_error("density/metric AO covariance mismatch");
                }
            }
            std::filesystem::remove(path);
        }
        std::filesystem::remove(directory);
        std::cout<<"FCHK AO transport and retained precision: 9 shell representations passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
}
