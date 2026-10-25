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

#include "cov/nbo_aomo_labels.hpp"
#include "cov/wavefunction_io.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
}
int main(int argc,char** argv){try{
    require(argc==4,"Usage: cov_nbo_salc_names_probe canonical.fchk analysis_directory output.json");
    const auto w=cov::parse_wavefunction(argv[1]);const auto identity=cov::nbo_canonical_fingerprint(w);
    const auto inputs=cov::discover_nbo_inputs({argv[2]});require(inputs.candidates.size()==1,"Expected one analysis input");
    const auto data=cov::read_nbo_integration(w,inputs.candidates.front());
    const auto raw=cov::build_nbo_salc_model(w,data);const auto averaged=cov::build_nbo_spin_averaged_model(w,data,raw);
    std::vector<std::string> failures;
    std::ofstream out(argv[3]);out<<"{\"models\":[";bool first_model=true;
    for(const auto* model:{&raw,&averaged}){
        if(!first_model&&!averaged.spin_averaged)continue;
        const auto names=cov::ui::build_nbo_aomo_names(w,data,model);
        std::vector<std::size_t> all_mo(w.orbitals.size()),all_side(model->orbitals.size());std::iota(all_mo.begin(),all_mo.end(),0);std::iota(all_side.begin(),all_side.end(),0);
        for(const auto* collection:{&names.canonical,&names.salc})for(const auto& name:*collection)if(name.approximate_dominant_label)
            require(name.label.starts_with("\xE2\x89\x88"),"Approximate name metadata lost its main label marker");
        const auto full=cov::ui::nbo_aomo_names_for_view(w,names,all_mo,all_side,model,"all copies");
        std::vector<std::size_t> partial;std::size_t certified=0;
        for(const auto& sub:model->subspaces)if(sub.symmetry_verified&&sub.multiplicity==1&&sub.irrep_dimension>1){
            ++certified;const auto first=sub.orbital_indices.front();const auto& ref=names.salc[first];
            if(!ref.verified||ref.partner_block_id.empty()||ref.partner_block_size!=sub.dimension){failures.push_back("Mathematical copy lost strict name/partners: "+sub.id+"; "+cov::ui::serialize_orbital_name_json(ref));continue;}
            for(auto i:sub.orbital_indices)require(names.salc[i].partner_block_id==ref.partner_block_id,"Members of one derived irreducible copy received different partner identities");
            partial.push_back(first);
        }
        // Keep all 1D entries and one member of each multidimensional copy so
        // the occurrence set stays fixed while visible member counts change.
        std::set<std::string> selected;for(auto i:partial)selected.insert(names.salc[i].partner_block_id);
        for(auto i:all_side)if(!selected.contains(names.salc[i].partner_block_id))partial.push_back(i);
        const auto filtered=cov::ui::nbo_aomo_names_for_view(w,names,all_mo,partial,model,"one partner visible");
        for(auto i:partial)if(selected.contains(names.salc[i].partner_block_id)){
            require(filtered.salc[i].ordinal==full.salc[i].ordinal,"Hiding partners changed the full-copy mean display order");
            require(filtered.salc[i].visible_partner_count==1,"Partial visible partner count lost");
            require(filtered.salc[i].complete_set_ordinal==names.salc[i].complete_set_ordinal,"View overwrote complete-set ordinal");
        }
        if(!first_model)out<<',';first_model=false;
        out<<"{\"spin_averaged\":"<<(model->spin_averaged?"true":"false")<<",\"geometry_group\":"<<std::quoted(model->point_group)<<",\"naming_group\":"<<std::quoted(model->used_group)<<",\"certified_multidimensional_copies\":"<<certified;
        const auto rows=[&](const char* key,const auto& entries){out<<",\""<<key<<"\":[";for(std::size_t i=0;i<entries.size();++i){if(i)out<<',';const auto& name=entries[i];out<<"{\"index0\":"<<i<<",\"label\":"<<std::quoted(name.label)<<",\"strict\":"<<(name.verified?"true":"false")<<",\"approximate\":"<<(name.approximate_dominant_label?"true":"false")<<",\"dominant_weight\":"<<std::setprecision(17)<<name.dominant_weight<<",\"irrep\":"<<std::quoted(name.irrep)<<",\"ordinal\":"<<name.ordinal<<",\"complete_ordinal\":"<<name.complete_set_ordinal<<",\"partner_id\":"<<std::quoted(name.partner_block_id)<<",\"partner_size\":"<<name.partner_block_size<<'}';}out<<']';};
        rows("canonical",full.canonical);rows("salc",full.salc);out<<'}';
        std::cout<<(model->spin_averaged?"spatial":"source-spin")<<" actual strict copies="<<certified<<"; partial-view numbering preserved\n";
    }
    require(cov::nbo_canonical_fingerprint(w)==identity,"Source canonical changed");out<<"],\"source_unchanged\":true}";for(const auto& failure:failures)std::cerr<<failure<<'\n';return failures.empty()?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
