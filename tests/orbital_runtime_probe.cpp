#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_tracking.hpp"
#include "cov/wavefunction_io.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main(int argc,char** argv){try{
    if(argc<3)throw std::runtime_error("names input.fchk [index1 ...], or tracking from.fchk to.fchk");
    const auto started=std::chrono::steady_clock::now();
    const auto first=cov::parse_wavefunction(argv[2]);
    if(std::string(argv[1])=="names"){
        const auto names=cov::ui::canonical_mo_names(first);std::size_t ordered=0,unordered=0,mixed=0,unanalysed=0;
        for(const auto& name:names->canonical){if(name.verified){if(name.ordinal)++ordered;else ++unordered;}
            else if(name.decomposition_verified)++mixed;else ++unanalysed;}
        std::cout<<"{\"total\":"<<first.orbitals.size()<<",\"ordered\":"<<ordered<<",\"unordered\":"<<unordered
            <<",\"mixed_decomposed\":"<<mixed<<",\"unanalysed\":"<<unanalysed<<",\"rows\":[";
        for(int a=3;a<argc;++a){const auto i=std::stoul(argv[a])-1;if(i>=names->canonical.size())throw std::runtime_error("source index out of range");
            if(a>3)std::cout<<',';std::cout<<"{\"index1\":"<<i+1<<",\"energy_hartree\":"<<std::setprecision(17)<<first.orbitals[i].energy_hartree
                <<",\"name\":"<<cov::ui::serialize_orbital_name_json(names->canonical[i])<<'}';}
        std::cout<<"]}";
    }else if(std::string(argv[1])=="tracking"){
        if(argc!=4)throw std::runtime_error("tracking requires two source files");
        const auto second=cov::parse_wavefunction(argv[3]);
        const auto before=std::chrono::steady_clock::now();const auto result=cov::track_orbital_subspaces(first,second);
        const double matching=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
        const auto members=[](const auto& groups){std::size_t count=0;for(const auto& g:groups)count+=g.size();return count;};
        std::size_t matched_from=0,matched_to=0;for(const auto& match:result.matches){matched_from+=match.from_members.size();matched_to+=match.to_members.size();}
        const bool complete=matched_from+members(result.unmatched_from)+members(result.unresolved_from)==first.orbitals.size()&&
            matched_to+members(result.unmatched_to)+members(result.unresolved_to)==second.orbitals.size();
        std::cout<<std::setprecision(17)<<"{\"from_mos\":"<<first.orbitals.size()<<",\"to_mos\":"<<second.orbitals.size()
            <<",\"matching_seconds\":"<<matching<<",\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
            <<",\"budget_exhausted\":"<<(result.tracking_budget_exhausted?"true":"false")<<",\"budget_reason\":"<<int(result.budget_exhaustion_reason)
            <<",\"budget_stage\":"<<int(result.budget_exhausted_stage)<<",\"work_units\":"<<result.tracking_work_units
            <<",\"matches\":"<<result.matches.size()<<",\"unmatched_from_members\":"<<members(result.unmatched_from)
            <<",\"unmatched_to_members\":"<<members(result.unmatched_to)<<",\"unresolved_from_members\":"<<members(result.unresolved_from)
            <<",\"unresolved_to_members\":"<<members(result.unresolved_to)<<",\"all_members_retained\":"<<(complete?"true":"false")<<'}';
        if(!complete)return 2;
    }else throw std::runtime_error("unknown mode");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
