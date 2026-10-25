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

#include "cov/pi_topology_evidence.hpp"
#include <iostream>
#include <stdexcept>

int main() {
    static_assert(static_cast<int>(cov::DelocalisedPiTopology::Spiro)==4);
    static_assert(static_cast<int>(cov::DelocalisedPiTopology::SymmetryDirectSum)==6);
    static_assert(static_cast<int>(cov::DelocalisedPiTopology::MultiChannel)==7);
    cov::Wavefunction wf;
    cov::DelocalisedPiAssignment pi;
    pi.family_id="graph-evidence-only\"quoted";pi.atoms={1,2,3,4};pi.orbitals={0,1,2,3};pi.electron_count=4;
    pi.topology=cov::DelocalisedPiTopology::Spiro;pi.provenance=cov::DataProvenance::Derived;
    pi.orientation_channels={{{1,2},{0,0,1},1,false},{{3,4},{0,1,0},1,false}};
    auto& graph=pi.topology_graph;graph.atom_count=5;graph.source=cov::PiTopologyGraphSource::MayerDistanceModel;
    graph.bond_order_provenance=cov::DataProvenance::Derived;
    graph.minimum_mayer_order=.05;graph.maximum_covalent_radius_factor=1.45;
    graph.edges={{0,1},{1,2},{2,0},{0,3},{3,4},{4,0}};
    graph.examined_hubs={0};graph.unscoped_cycle_union_hubs={0};graph.channel_association_search_complete=true;
    const auto search=cov::find_channel_cycle_union_at_hub(5,graph.edges,0,{{{1,2},{3,4}}});
    if(!search.witness)throw std::logic_error("Expected two ring paths");
    graph.channel_ring_witnesses.push_back({{{0,1}},*search.witness});
    wf.delocalised_pi_assignments.push_back(pi);
    std::cout<<"{\"labels\":[";
    for(int i=0;i<8;++i) {if(i)std::cout<<',';cov::numerical_json::string(std::cout,cov::pi_topology_name(static_cast<cov::DelocalisedPiTopology>(i)));}
    std::cout<<"],\"unavailable\":";cov::write_pi_topology_graph_json(std::cout,{});
    std::cout<<",\"graph\":";cov::write_pi_topology_graph_json(std::cout,graph);
    std::cout<<",\"assignments\":";cov::write_pi_topology_assignments_json(std::cout,wf);
    graph.source=cov::PiTopologyGraphSource::CovalentDistanceModel;
    graph.bond_order_provenance=cov::DataProvenance::Unavailable;
    graph.minimum_mayer_order=std::numeric_limits<double>::quiet_NaN();graph.maximum_covalent_radius_factor=1.22;
    std::cout<<",\"geometry_model\":";cov::write_pi_topology_graph_json(std::cout,graph);
    std::cout<<"}\n";
}
