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
#include <iostream>
#include <stdexcept>

namespace { void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);} }
int main(){try{
    cov::RoutedAnalysis route;route.local_relations_normalized=true;
    const cov::NboOrbitalRef a{cov::NboOrbitalKind::NBO,cov::NboSpin::Total,0};
    const cov::NboOrbitalRef b{cov::NboOrbitalKind::NBO,cov::NboSpin::Total,1};
    const cov::NboOrbitalRef c{cov::NboOrbitalKind::NLMO,cov::NboSpin::Total,0};
    cov::RoutedLocalRelation first;first.id="ordered-first";first.orbitals={a,b};
    first.kind="donor_acceptor";first.printed_value=1.25;first.units="kcal/mol";
    first.source.path="preserved/source.log";first.source.line_begin=123;
    auto second=first;second.id="ordered-second";second.orbitals={a,c,b};
    route.local_relation_registry={first,second};
    route.local_source_projections={{a,{.2,0.,std::nullopt}},
        {b,{1e-8,0.,1e-8}},{c,{std::nullopt,.3,std::nullopt}}};
    route.mo_relation_counts={2,1,0};route.mo_relations.resize(3);
    for(std::size_t i=0;i<2;++i){route.mo_relations[i].status=cov::RoutedStatus::Available;
        route.mo_relations[i].provider=cov::RoutedProvider::Nbo;
        route.mo_relations[i].value=std::vector<cov::RoutedRelationProjection>{};}
    route.mo_relations[2].status=cov::RoutedStatus::NotReportedAboveThreshold;
    const auto zero=cov::routed_mo_relations(route,0);
    require(zero.available()&&zero.value->size()==2,"complete local-record order restored");
    require((*zero.value)[0].relation_id==first.id&&(*zero.value)[1].relation_id==second.id,
        "registry order preserved");
    require((*zero.value)[0].canonical_projection_weights==std::vector<std::optional<double>>{.2,1e-8},
        "below-cutoff endpoint retained once a relation is supported");
    require((*zero.value)[1].canonical_projection_weights==std::vector<std::optional<double>>{.2,std::nullopt,1e-8},
        "missing endpoint remains null");
    const auto one=cov::routed_mo_relations(route,1);
    require(one.available()&&one.value->size()==1&&
        one.value->front().canonical_projection_weights==std::vector<std::optional<double>>{0.,.3,0.},
        "measured zero is distinct from missing projection");
    require(!cov::routed_mo_relations(route,2).available(),"weak-only evidence cannot invent relation support");
    auto legacy=route;legacy.local_relations_normalized=false;
    legacy.mo_relations[0]=zero;legacy.mo_relations[1]=one;
    require(cov::serialize_routed_analysis_json(route,false,true)==cov::serialize_routed_analysis_json(legacy),
        "explicit v1 expansion preserves every old serialized field");
    const auto compact=cov::serialize_routed_analysis_json(route);
    require(compact.find("cov.chemistry.route.v2")!=std::string::npos&&
        compact.find("shared-source-projections-v1")!=std::string::npos&&
        compact.find("canonical_projection_weights")==std::string::npos,
        "default v2 shares measured source weights instead of repeating relation tuples");
    const auto scoped=cov::serialize_routed_analysis_json(route,false,false,std::vector<std::size_t>{1});
    require(scoped.find("ordered-first")==std::string::npos&&scoped.find("ordered-second")!=std::string::npos&&
        scoped.find("\"canonical_indices\":[1]")!=std::string::npos,
        "current-source scope retains only supported local records and global source indices");
    const auto empty=cov::serialize_routed_analysis_json(route,false,false,std::vector<std::size_t>{});
    require(empty.find("\"local_relation_registry\":[]")!=std::string::npos&&
        empty.find("\"mo_relations\":[]")!=std::string::npos,"empty explicit display scope is not full analysis");
    bool invalid_scope=false;
    try{cov::serialize_routed_analysis_json(route,false,false,std::vector<std::size_t>{3});}
    catch(const std::invalid_argument&){invalid_scope=true;}
    require(invalid_scope,"out-of-source scope must not silently select a different MO");
    auto wrong=route;wrong.mo_relation_counts[0]=1;
    require(cov::routed_mo_relations(wrong,0).status==cov::RoutedStatus::Rejected,"incorrect count rejected");
    wrong=route;wrong.local_source_projections.push_back(wrong.local_source_projections.front());
    require(cov::routed_mo_relations(wrong,0).status==cov::RoutedStatus::Rejected,"ambiguous source column rejected");
    wrong=route;wrong.local_source_projections.pop_back();
    require(cov::routed_mo_relations(wrong,0).status==cov::RoutedStatus::Rejected,"missing registry source column rejected");
    std::cout<<"Normalized local relation projections preserve exact expanded semantics\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
