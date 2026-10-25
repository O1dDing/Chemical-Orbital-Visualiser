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

#include "cov/chemistry_route.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/wavefunction_io.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv){try{
    if(argc!=5&&argc!=6)throw std::runtime_error("Usage: route_relation_projection_probe fchk package output-prefix source-index [expanded]");
    const auto start=std::chrono::steady_clock::now();
    const auto canonical=cov::parse_wavefunction(argv[1]);
    const auto identity=cov::nbo_canonical_fingerprint(canonical);
    const auto found=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
    if(found.candidates.size()!=1)throw std::runtime_error("Exactly one associated source analysis required");
    auto integration=cov::read_nbo_integration(canonical,found.candidates.front());
    cov::annotate_nbo_bond_channels(integration,canonical);
    const auto route=cov::route_chemistry(canonical,&integration);
    const auto index=static_cast<std::size_t>(std::stoull(argv[4]));
    const auto save=[&](const char* suffix,const std::string& contents){
        std::ofstream out(std::string(argv[3])+suffix);out<<contents;
        if(!out)throw std::runtime_error("Could not save probe output");return contents.size();};
    const auto normalized=save(".normalized.json",cov::serialize_routed_analysis_json(route));
    const auto scoped=save(".scoped.json",cov::serialize_routed_analysis_json(route,false,false,std::vector<std::size_t>{index}));
    std::size_t expanded=0;
    if(argc==6&&std::string(argv[5])=="expanded")
        expanded=save(".expanded.json",cov::serialize_routed_analysis_json(route,false,true));
    std::size_t stored_tuples=0,relation_count=0,measurements=0;
    for(const auto& result:route.mo_relations)if(result.value)stored_tuples+=result.value->size();
    for(const auto count:route.mo_relation_counts)relation_count+=count;
    for(const auto& column:route.local_source_projections)for(const auto& weight:column.canonical_weights)if(weight)++measurements;
    const auto selected=cov::routed_mo_relations(route,index);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream stats(std::string(argv[3])+".stats.json");
    stats<<"{\"normalized_bytes\":"<<normalized<<",\"scoped_bytes\":"<<scoped<<",\"expanded_bytes\":"<<expanded
         <<",\"stored_expanded_tuples\":"<<stored_tuples<<",\"recoverable_tuples\":"<<relation_count
         <<",\"unique_projection_measurements\":"<<measurements<<",\"selected_relation_count\":"
         <<(selected.available()?selected.value->size():0)<<",\"seconds\":"<<seconds
         <<",\"canonical_immutable\":"<<(identity==cov::nbo_canonical_fingerprint(canonical)?"true":"false")<<'}';
    if(!stats)throw std::runtime_error("Could not save probe stats");
    std::cout<<"normalized="<<normalized<<" scoped="<<scoped<<" stored-tuples="<<stored_tuples<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
