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

#include "cov/fchk_nuclear_charge.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace cov {
namespace {

std::string trim(std::string value) {
    const auto not_space=[](const unsigned char c) {
        return !std::isspace(c);
    };
    value.erase(value.begin(),std::find_if(value.begin(),value.end(),not_space));
    value.erase(std::find_if(value.rbegin(),value.rend(),not_space).base(),value.end());
    return value;
}

double parse_real(std::string token) {
    for (char& c:token) if (c=='D' || c=='d') c='E';
    return std::stod(token);
}

} // namespace

bool enrich_fchk_nuclear_charges_from_file(
    Wavefunction& wavefunction,
    const std::filesystem::path& path) {
    for (auto& atom:wavefunction.atoms) {
        if (!(atom.nuclear_charge>0.0)) {
            atom.nuclear_charge=static_cast<double>(atom.atomic_number);
        }
    }

    std::ifstream input(path,std::ios::binary);
    if (!input) return false;
    std::string line;
    while (std::getline(input,line)) {
        if (line.size()<41u) continue;
        const std::string label=trim(line.substr(0,40));
        if (label!="Nuclear charges") continue;
        const std::string tail=line.substr(40);
        const auto marker=tail.find("N=");
        if (marker==std::string::npos) return false;
        std::istringstream count_stream(tail.substr(marker+2u));
        std::size_t count=0;
        count_stream>>count;
        if (!count_stream || count!=wavefunction.atoms.size()) return false;

        std::vector<double> values;
        values.reserve(count);
        while (values.size()<count && std::getline(input,line)) {
            std::istringstream row(line);
            std::string token;
            while (row>>token && values.size()<count) {
                values.push_back(parse_real(token));
            }
        }
        if (values.size()!=count) return false;
        for (std::size_t i=0;i<count;++i) {
            if (values[i]>0.0) wavefunction.atoms[i].nuclear_charge=values[i];
        }
        return true;
    }
    return false;
}

} // namespace cov
