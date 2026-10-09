#include "cov/nbo_molecular_overlay.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/wavefunction_io.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
int main(int argc,char** argv){try{
    if(argc!=4)throw std::runtime_error("canonical.fchk package output.json");
    const auto w=cov::parse_wavefunction(argv[1]);
    const auto candidates=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
    if(candidates.candidates.size()!=1)throw std::runtime_error("Expected one immutable source");
    const auto data=cov::read_nbo_integration(w,candidates.candidates[0]);
    const auto route=cov::route_chemistry(w,&data);
    if(!route.interaction_graph.available())throw std::runtime_error("Validated graph unavailable");
    const auto scope=cov::analyse_nbo_electronic_symmetry_scope(w);
    const auto skeleton=cov::make_nbo_molecule_overlay(data,*route.interaction_graph.value,
        w.atoms.size(),{},std::nullopt,cov::AtomScalarMode::Element,false,false,&route,&w,
        cov::NboBondDisplayMode::DefaultSkeleton,&scope);
    const auto lewis=cov::make_nbo_molecule_overlay(data,*route.interaction_graph.value,
        w.atoms.size(),{},std::nullopt,cov::AtomScalarMode::Element,false,false,&route,&w,
        cov::NboBondDisplayMode::LewisStructure,&scope);
    const auto verify_same=[&](const cov::MoleculeOverlay& candidate,bool reversed_atoms){
        std::map<std::pair<std::size_t,std::size_t>,const cov::MoleculeOverlayBond*> observed;
        for(const auto& b:candidate.bonds){const auto a=reversed_atoms?w.atoms.size()-1-b.atom_a:b.atom_a;
            const auto c=reversed_atoms?w.atoms.size()-1-b.atom_b:b.atom_b;
            if(!observed.emplace(std::minmax(a,c),&b).second)throw std::runtime_error("Permutation duplicated pair");}
        if(observed.size()!=skeleton.bonds.size())throw std::runtime_error("Permutation changed skeleton size");
        for(const auto& b:skeleton.bonds){const auto found=observed.find(std::minmax(b.atom_a,b.atom_b));
            if(found==observed.end())throw std::runtime_error("Permutation lost physical atom pair");
            const auto& other=*found->second;
            if(other.multiplicity!=b.multiplicity || other.style!=b.style || other.continuous_index!=b.continuous_index ||
               other.evidence_indices.size()!=b.evidence_indices.size() || other.lewis_equivalence_conflict!=b.lewis_equivalence_conflict)
                throw std::runtime_error("Permutation changed physical skeleton evidence");}
        if(candidate.multicentre.size()!=skeleton.multicentre.size())throw std::runtime_error("Permutation changed multicentre group count");
    };
    auto reordered_data=data;std::reverse(reordered_data.structure.begin(),reordered_data.structure.end());
    verify_same(cov::make_nbo_molecule_overlay(reordered_data,*route.interaction_graph.value,w.atoms.size(),{},
        std::nullopt,cov::AtomScalarMode::Element,false,false,&route,&w,cov::NboBondDisplayMode::DefaultSkeleton,&scope),false);
    auto permuted_w=w;std::reverse(permuted_w.atoms.begin(),permuted_w.atoms.end());
    for(auto& s:permuted_w.shells)s.atom_index=w.atoms.size()-1-s.atom_index;
    auto permuted_data=data;for(auto& e:permuted_data.structure)for(auto& a:e.atoms)a=w.atoms.size()-1-a;
    auto permuted_graph=*route.interaction_graph.value;for(auto& e:permuted_graph.edges){e.atom_a=w.atoms.size()-1-e.atom_a;e.atom_b=w.atoms.size()-1-e.atom_b;}
    auto permuted_scope=scope;for(auto& op:permuted_scope.operations){auto old=op.atom_permutation;
        for(std::size_t i=0;i<old.size();++i)op.atom_permutation[i]=old.size()-1-old[old.size()-1-i];}
    verify_same(cov::make_nbo_molecule_overlay(permuted_data,permuted_graph,w.atoms.size(),{},
        std::nullopt,cov::AtomScalarMode::Element,false,false,&route,&permuted_w,cov::NboBondDisplayMode::DefaultSkeleton,&permuted_scope),true);
    std::ofstream out(argv[3]);out.precision(17);
    out<<"{\"atom_permutation_invariant\":true,\"source_order_invariant\":true,\"density_checked\":"<<(scope.density_checked?"true":"false")
       <<",\"density_scope_verified\":"<<(scope.naming_scope_verified?"true":"false")
       <<",\"scope_status\":\""<<scope.status<<"\",\"geometry_group\":\""<<scope.geometry_group
       <<"\",\"naming_group\":\""<<scope.naming_group<<"\",\"operation_count\":"<<scope.operations.size()
       <<",\"total_density_residual_max\":";
    if(scope.total_density_residuals.empty())out<<"null";
    else out<<*std::max_element(scope.total_density_residuals.begin(),scope.total_density_residuals.end());
    out<<",\"spin_density_residual_max\":";
    if(scope.spin_density_residuals.empty())out<<"null";
    else out<<*std::max_element(scope.spin_density_residuals.begin(),scope.spin_density_residuals.end());
    out<<",\"skeleton\":[";
    std::set<std::pair<std::size_t,std::size_t>> pairs;bool comma=false;
    for(const auto& b:skeleton.bonds){
        if(!pairs.insert(std::minmax(b.atom_a,b.atom_b)).second)throw std::runtime_error("Duplicate skeleton pair");
        if(comma)out<<',';comma=true;
        out<<"{\"atoms\":["<<b.atom_a<<','<<b.atom_b<<"],\"atomic_numbers\":["
           <<w.atoms[b.atom_a].atomic_number<<','<<w.atoms[b.atom_b].atomic_number
           <<"],\"multiplicity\":"<<b.multiplicity<<",\"lewis_count\":";
        if(b.lewis_bond_count)out<<*b.lewis_bond_count;else out<<"null";
        out<<",\"wiberg\":";if(b.continuous_index)out<<*b.continuous_index;else out<<"null";
        out<<",\"style\":"<<int(b.style)<<",\"equivalence_conflict\":"
           <<(b.lewis_equivalence_conflict?"true":"false")<<",\"source_records\":"<<b.evidence_indices.size()<<'}';
    }
    out<<"],\"lewis\":[";comma=false;
    for(const auto& b:lewis.bonds){if(comma)out<<',';comma=true;
        out<<"{\"atoms\":["<<b.atom_a<<','<<b.atom_b<<"],\"multiplicity\":"<<b.multiplicity<<'}';}
    out<<"],\"multicentre_groups\":"<<skeleton.multicentre.size()
       <<",\"e2_relations_default\":"<<skeleton.relations.size()<<'}';
    if(!out)throw std::runtime_error("Cannot write observation");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
