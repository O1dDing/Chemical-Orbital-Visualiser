#include "cov/orbital_group_bonding.hpp"
#include "cov/nbo_integration.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/wavefunction_io.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
int main(int argc,char** argv){try{
    if(argc!=4 && argc!=5)throw std::runtime_error("canonical.fchk package output.json [inventory.tsv]");
    const auto w=cov::parse_wavefunction(argv[1]);const auto candidates=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
    if(candidates.candidates.size()!=1)throw std::runtime_error("Expected complete immutable source");
    const auto data=cov::read_nbo_integration(w,candidates.candidates.front());const auto route=cov::route_chemistry(w,&data);
    if(!route.interaction_graph.available())throw std::runtime_error("Validated graph unavailable");
    const auto scope=cov::make_fixed_bonding_scope(w,*route.interaction_graph.value);
    std::ofstream out(argv[3]);out.precision(17);out<<"{\"scope_id\":\""<<scope.id<<"\",\"scope_source\":\""<<scope.source<<"\",\"edge_count\":"<<scope.edges.size()<<",\"groups\":[";
    const bool pg006=std::filesystem::path(argv[2]).filename()=="PG-006";
    if(pg006){
        // The historical 48-pair diagnostic included every |Mayer|>=.02
        // pair, including non-skeletal pairs. Its independent reference must
        // not be relabelled as the current 36-edge molecular skeleton.
        std::vector<std::pair<std::size_t,std::size_t>> legacy_pairs;
        for(const auto& b:w.bond_orders)if(b.provenance!=cov::DataProvenance::Unavailable && std::isfinite(b.mayer_order) && std::abs(b.mayer_order)>=.02)
            legacy_pairs.emplace_back(b.atom_a,b.atom_b);
        const auto legacy_scope=cov::make_fixed_bonding_scope(w,legacy_pairs,"legacy_reference_all_mayer_abs_ge_0p02");
        if(legacy_scope.edges.size()!=48)throw std::runtime_error("PG006 original all-Mayer diagnostic must have 48 pairs");
        if(scope.edges.size()!=36)throw std::runtime_error("PG006 validated structural skeleton must have 36 edges");
        const std::array<std::array<std::size_t,2>,3> members{{{{49,50}},{{66,67}},{{69,70}}}};
        const std::array<std::array<double,2>,3> expected{{{{.310022102,.310023014}},{{.254704695,.254704746}},{{-.636717983,-.636706940}}}};
        for(std::size_t i=0;i<3;++i){cov::OrbitalGroupBondingOptions options;options.expected_dimension=2;
            const auto r=cov::analyse_orbital_group_bonding(w,legacy_scope,{members[i][0],members[i][1]},options);
            const auto skeletal=cov::analyse_orbital_group_bonding(w,scope,{members[i][0],members[i][1]},options);
            if(r.eigenvalues.size()!=2)throw std::runtime_error("PG006 complete group unavailable: "+r.detail);
            for(unsigned k=0;k<2;++k)if(std::abs(r.eigenvalues[k]-expected[i][k])>2e-5)throw std::runtime_error("PG006 independent reference spectrum disagrees");
            if(skeletal.eigenvalues.size()!=2)throw std::runtime_error("PG006 current skeletal group unavailable");
            if(i)out<<',';out<<"{\"members\":["<<members[i][0]+1<<','<<members[i][1]+1<<"],\"reference_scope_id\":\""<<legacy_scope.id
                <<"\",\"reference_pair_count\":"<<legacy_scope.edges.size()<<",\"status\":\""<<cov::orbital_group_bonding_status_name(r.status)<<"\",\"spectrum\":["<<r.eigenvalues[0]<<','<<r.eigenvalues[1]<<"],\"trace\":"<<r.trace<<",\"error_bound\":"<<r.error_bound<<",\"metric_error\":"<<r.metric_error<<",\"eigensolver_residual\":"<<r.eigensolver_residual
                <<",\"current_skeleton_scope_id\":\""<<skeletal.scope_id<<"\",\"current_skeleton_edge_count\":"<<scope.edges.size()<<",\"current_skeleton_status\":\""<<cov::orbital_group_bonding_status_name(skeletal.status)
                <<"\",\"current_skeleton_spectrum\":["<<skeletal.eigenvalues[0]<<','<<skeletal.eigenvalues[1]<<"]}";
        }
    }
    out<<"],\"single_member_statuses\":[";std::map<std::string,std::size_t> counts;
    if(argc==4)for(std::size_t i=0;i<w.orbitals.size();++i){const auto r=cov::analyse_orbital_group_bonding(w,scope,{i});const std::string status=cov::orbital_group_bonding_status_name(r.status);++counts[status];if(i)out<<',';out<<"{\"member\":"<<i+1<<",\"status\":\""<<status<<"\",\"spectrum\":";
        if(r.eigenvalues.empty())out<<"null";else out<<r.eigenvalues.front();out<<'}';}
    out<<"],\"inventory\":[";bool row_comma=false;
    if(argc==5){std::ifstream input(argv[4]);if(!input)throw std::runtime_error("Cannot read indexed inventory");std::string line;
        while(std::getline(input,line)){if(line.empty())continue;const auto separator=line.find('\t');if(separator==std::string::npos)throw std::runtime_error("Invalid indexed inventory row");
            const auto id=line.substr(0,separator);if(id.find_first_of("\"\\")!=std::string::npos)throw std::runtime_error("Unsafe record ID");
            std::istringstream columns(line.substr(separator+1));std::vector<std::size_t> members;std::size_t member;while(columns>>member)members.push_back(member);
            cov::OrbitalGroupBondingOptions options;options.expected_dimension=members.size();const auto r=cov::analyse_orbital_group_bonding(w,scope,members,options);
            const std::string status=cov::orbital_group_bonding_status_name(r.status);++counts[status];if(row_comma)out<<',';row_comma=true;
            out<<"{\"record_id\":\""<<id<<"\",\"members0\":[";for(std::size_t i=0;i<members.size();++i){if(i)out<<',';out<<members[i];}
            out<<"],\"status\":\""<<status<<"\",\"detail\":\""<<r.detail<<"\",\"spectrum\":[";for(std::size_t i=0;i<r.eigenvalues.size();++i){if(i)out<<',';out<<r.eigenvalues[i];}
            out<<"],\"trace\":"<<r.trace<<",\"electron_weighted_trace\":"<<r.electron_weighted_trace<<",\"electron_count\":"<<r.electron_count<<",\"error_bound\":"<<r.error_bound
               <<",\"metric_error\":"<<r.metric_error<<",\"metric_rank\":"<<r.metric_rank<<",\"energy_span_hartree\":"<<r.source_energy_span_hartree<<'}';
        }
    }
    out<<"],\"status_counts\":{";bool comma=false;for(const auto& [key,count]:counts){if(comma)out<<',';comma=true;out<<'"'<<key<<"\":"<<count;}out<<"}}";
    if(!out)throw std::runtime_error("Cannot write result");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
