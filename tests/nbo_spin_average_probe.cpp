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

#include "cov/nbo_spin_average.hpp"
#include "cov/wavefunction_io.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv){try{
    if(argc!=4&&argc!=5)throw std::runtime_error("Usage: cov_nbo_spin_average_probe canonical.fchk analysis_directory output.json [common-only]");
    const bool common_only=argc==5&&std::string(argv[4])=="common-only";
    if(argc==5&&!common_only)throw std::runtime_error("Unknown probe mode");
    const auto start=std::chrono::steady_clock::now();
    const auto w=cov::parse_wavefunction(argv[1]);const auto original=w;
    const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1)throw std::runtime_error("Expected exactly one analysis candidate");
    const auto data=cov::read_nbo_integration(w,found.candidates.front());
    const auto raw=cov::build_nbo_salc_model(w,data);const auto before=cov::serialize_nbo_salc_json(raw);
    const auto common=cov::build_nbo_ro_common_energy(w,data,raw);
    if(common_only){const std::filesystem::path output(argv[3]);if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
        std::ofstream out(output);out<<cov::serialize_nbo_ro_common_energy_json(common);if(!out)throw std::runtime_error("Failed writing common energy output");
        std::cout<<common.status<<" count="<<common.orbitals.size()<<'\n';return 0;}
    const auto combined=cov::build_nbo_spin_averaged_model(w,data,raw);
    bool preserved=w.ao_overlap==original.ao_overlap&&w.orbitals.size()==original.orbitals.size();
    for(std::size_t i=0;i<w.orbitals.size();++i){const auto& a=w.orbitals[i];const auto& b=original.orbitals[i];
        preserved=preserved&&a.coefficients==b.coefficients&&a.gaussian_source_coefficients==b.gaussian_source_coefficients&&
            a.energy_hartree==b.energy_hartree&&a.occupation==b.occupation&&a.spin==b.spin&&a.source_orbital_index==b.source_orbital_index;}
    std::size_t selection_failures=0,with_energy=0,with_occupation=0;double norm_error=0;
    for(std::size_t i=0;i<combined.orbitals.size();++i){const auto& side=combined.orbitals[i];
        with_energy+=side.energy_hartree.has_value();with_occupation+=side.occupation.has_value();
        const auto selection=cov::nbo_salc_selection(combined,i);const auto view=cov::make_nbo_selection_view(data,w,selection);
        if(!view.available)++selection_failures;else norm_error=std::max(norm_error,std::abs(view.metric_norm2[0]-1));}
    const std::filesystem::path output(argv[3]);if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output);out<<std::boolalpha<<std::setprecision(17)<<"{\"canonical_preserved\":"<<preserved
        <<",\"source_model_preserved\":"<<(before==cov::serialize_nbo_salc_json(raw))<<",\"source_count\":"<<raw.orbitals.size()
        <<",\"display_count\":"<<combined.orbitals.size()<<",\"selection_failures\":"<<selection_failures<<",\"norm_error\":"<<norm_error
        <<",\"with_energy\":"<<with_energy<<",\"with_occupation\":"<<with_occupation
        <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"common_energy\":"<<cov::serialize_nbo_ro_common_energy_json(common)<<",\"model\":"<<cov::serialize_nbo_salc_json(combined)<<'}';
    if(!out)throw std::runtime_error("Failed writing probe output");
    std::cout<<"source="<<raw.orbitals.size()<<" display="<<combined.orbitals.size()<<" merged="<<combined.merged_spatial_count
        <<" separate="<<combined.separate_spin_count<<" energy="<<with_energy<<" selections_failed="<<selection_failures<<'\n';
    return preserved&&selection_failures==0?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
