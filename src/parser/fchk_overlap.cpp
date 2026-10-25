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

#include "cov/fchk_overlap.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cov {
namespace {

std::string trim(std::string value) {
    auto not_space=[](unsigned char c){ return !std::isspace(c); };
    value.erase(value.begin(),std::find_if(value.begin(),value.end(),not_space));
    value.erase(std::find_if(value.rbegin(),value.rend(),not_space).base(),value.end());
    return value;
}

double parse_real(std::string token) {
    for (char& c:token) if (c=='D'||c=='d') c='E';
    std::size_t used=0;
    const double value=std::stod(token,&used);
    if (used!=token.size()) throw std::runtime_error("Invalid FCHK overlap value: "+token);
    return value;
}

std::size_t packed_index(std::size_t i,std::size_t j) {
    if (j>i) std::swap(i,j);
    return i*(i+1u)/2u+j;
}

std::vector<GaussianAoTransformEntry> internal_to_fchk_basis_map(const Wavefunction& wf) {
    if (wf.gaussian_ao_transform.size()==wf.basis_count) return wf.gaussian_ao_transform;
    std::vector<GaussianAoTransformEntry> map;
    map.reserve(wf.basis_count);
    std::size_t source_offset=0;
    for (const auto& shell:wf.shells) {
        const std::size_t count=shell_basis_count(shell);
        const int l=static_cast<int>(shell.angular_momentum);
        const auto local=gaussian_ao_transform(shell.pure ? -l : l, source_offset);
        map.insert(map.end(),local.begin(),local.end());
        source_offset+=count;
    }
    if (map.size()!=wf.basis_count || source_offset!=wf.basis_count) return {};
    return map;
}

} // namespace

bool enrich_fchk_overlap_from_file(Wavefunction& wavefunction,
                                   const std::filesystem::path& fchk_path) {
    if (wavefunction.basis_count==0) return false;
    std::ifstream input(fchk_path,std::ios::binary);
    if (!input) return false;

    const std::size_t n=wavefunction.basis_count;
    const std::size_t expected=n*(n+1u)/2u;
    std::vector<double> packed;
    std::string line;
    while (std::getline(input,line)) {
        if (line.size()<41u) continue;
        const std::string label=trim(line.substr(0,40));
        if (label!="Overlap Matrix") continue;
        const std::string tail=line.substr(40);
        if (tail.find('R')==std::string::npos) return false;
        const auto npos=tail.find("N=");
        if (npos==std::string::npos) return false;
        std::istringstream count_stream(tail.substr(npos+2u));
        std::size_t count=0;
        count_stream>>count;
        if (!count_stream || count!=expected) {
            throw std::runtime_error("FCHK Overlap Matrix dimension does not match AO basis");
        }
        packed.reserve(count);
        while (packed.size()<count && std::getline(input,line)) {
            std::istringstream values(line);
            std::string token;
            while (values>>token) {
                packed.push_back(parse_real(token));
                if (packed.size()==count) break;
            }
        }
        if (packed.size()!=count) {
            throw std::runtime_error("Unexpected end of FCHK Overlap Matrix");
        }
        break;
    }
    if (packed.empty()) return false;

    const auto basis_map=internal_to_fchk_basis_map(wavefunction);
    if (basis_map.size()!=n) {
        throw std::runtime_error("FCHK overlap AO-order transform is unavailable for this basis");
    }

    wavefunction.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<=i;++j) {
            const auto& left=basis_map[i];
            const auto& right=basis_map[j];
            const double value=left.basis_scale*right.basis_scale*
                packed[packed_index(left.source_index,right.source_index)];
            wavefunction.ao_overlap[i*n+j]=value;
            wavefunction.ao_overlap[j*n+i]=value;
        }
    }
    wavefunction.ao_overlap_provenance=DataProvenance::Producer;
    wavefunction.producer_ao_overlap=wavefunction.ao_overlap;
    wavefunction.ao_overlap_orthonormality_error=std::numeric_limits<double>::quiet_NaN();
    return true;
}

} // namespace cov
