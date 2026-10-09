#include "cov/nbo_molecular_overlay.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/molecule_style.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace {
void require(bool condition,const char* text){if(!condition)throw std::runtime_error(text);}
cov::NboStructureEvidence bond(std::size_t atom,unsigned count,double wbi){
    cov::NboStructureEvidence e;e.kind="bond";e.spin=cov::NboSpin::Total;
    e.atoms={0,atom};e.lewis_bond_count=count;e.wiberg=wbi;return e;
}
void source_bd(cov::NboIntegration& data,std::size_t evidence){
    auto& e=data.structure[evidence];
    for(unsigned i=0;i<e.lewis_bond_count.value_or(0);++i){
        cov::NboOrbital raw;raw.spin=e.spin;raw.id=data.dataset.orbitals.size()+1;raw.kind="BD";
        data.dataset.orbitals.push_back(raw);
        cov::NboOrbitalDescriptor d;d.ref={cov::NboOrbitalKind::NBO,e.spin,raw.id-1};
        d.id="source-"+std::to_string(raw.id);d.occupation=1.99;d.metric_norm2=1;d.coefficients={1};
        data.orbitals.push_back(d);e.orbitals.push_back(d.ref);
    }
}
}
int main(){try{
    cov::Wavefunction w;w.atoms={{"C",6,0,0,0,6},{"O",8,2.4,0,0,8},{"O",8,-2.4,0,0,8}};
    cov::InteractionGraph g;
    for(unsigned atom:{1u,2u}){cov::InteractionEdge e;e.atom_a=0;e.atom_b=atom;
        e.kind=cov::InteractionKind::CovalentConnectivity;e.strength=cov::InteractionStrength::StrongConnectivity;e.wiberg_index=1.4;g.edges.push_back(e);}
    cov::NboIntegration data;data.structure={bond(1,2,1.4),bond(2,1,1.4)};
    source_bd(data,0);source_bd(data,1);
    cov::NboElectronicSymmetryScope scope;scope.density_checked=true;scope.naming_scope_verified=true;
    cov::SymmetryOperation identity;identity.atom_permutation={0,1,2};
    cov::SymmetryOperation swap;swap.atom_permutation={0,2,1};scope.operations={identity,swap};
    const auto make=[&](cov::NboBondDisplayMode mode,const cov::NboElectronicSymmetryScope* s){
        return cov::make_nbo_molecule_overlay(data,g,3,{},std::nullopt,cov::AtomScalarMode::Element,false,false,nullptr,&w,mode,s);};
    auto x=make(cov::NboBondDisplayMode::DefaultSkeleton,&scope);
    require(x.bonds.size()==2,"deduplicated skeleton lost connections");
    require(std::all_of(x.bonds.begin(),x.bonds.end(),[](const auto& b){return b.multiplicity==1 && b.lewis_equivalence_conflict;}),"density-equivalent Lewis alternatives broke default equivalence");
    x=make(cov::NboBondDisplayMode::LewisStructure,&scope);
    require(x.bonds[0].multiplicity==2 && x.bonds[1].multiplicity==1,"explicit Lewis view lost original integers");
    std::reverse(data.structure.begin(),data.structure.end());
    x=make(cov::NboBondDisplayMode::DefaultSkeleton,&scope);
    require(std::all_of(x.bonds.begin(),x.bonds.end(),[](const auto& b){return b.multiplicity==1 && b.lewis_equivalence_conflict;}),"source order changed the equivalent default skeleton");
    std::reverse(data.structure.begin(),data.structure.end());
    scope.operations={identity}; // contract fixture: physical inequivalence proved by caller
    x=make(cov::NboBondDisplayMode::DefaultSkeleton,&scope);
    require(x.bonds[0].multiplicity==2 && x.bonds[1].multiplicity==1,"stable local multiple bond flattened despite verified evidence");
    scope.naming_scope_verified=false;
    x=make(cov::NboBondDisplayMode::DefaultSkeleton,&scope);
    require(x.bonds[0].multiplicity==1,"ambiguous density scope certified a multiple bond");
    auto duplicate=data.structure[0];duplicate.spin=cov::NboSpin::Alpha;duplicate.lewis_bond_count=1;
    duplicate.orbitals.clear();data.structure.push_back(duplicate);
    x=make(cov::NboBondDisplayMode::DefaultSkeleton,nullptr);
    require(x.bonds.size()==2 && x.bonds[0].evidence_indices.size()==2,"spin source records duplicated skeleton or were discarded");
    require(x.bonds[0].lewis_bond_count==2,"spin source overwrote total Lewis evidence");
    const auto select_alpha=[&](std::size_t selected,cov::NboBondDisplayMode mode){
        return cov::make_nbo_molecule_overlay(data,g,3,{},selected,cov::AtomScalarMode::Element,false,false,nullptr,&w,mode,nullptr);};
    x=select_alpha(2,cov::NboBondDisplayMode::LewisStructure);
    require(x.bonds[0].selected && x.bonds[0].evidence_index==2 && x.bonds[0].multiplicity==1 &&
        x.bonds[0].lewis_bond_count==1 && x.bonds[0].evidence_indices.size()==2,
        "explicit Lewis selected spin source inherited another source's integer");
    std::reverse(data.structure.begin(),data.structure.end());
    x=select_alpha(0,cov::NboBondDisplayMode::LewisStructure);
    require(x.bonds[0].evidence_index==0 && x.bonds[0].multiplicity==1 && x.bonds[0].lewis_bond_count==1,
        "later total record overwrote selected explicit spin Lewis integer");
    x=select_alpha(0,cov::NboBondDisplayMode::DefaultSkeleton);
    require(x.bonds[0].multiplicity==1 && x.bonds[0].lewis_bond_count==2,
        "explicit source selection leaked into the default total skeleton");
    std::reverse(data.structure.begin(),data.structure.end());
    data.structure[0].wiberg=0;g.edges[0].wiberg_index=0;g.edges[0].strength=cov::InteractionStrength::WeakContact;
    x=make(cov::NboBondDisplayMode::DefaultSkeleton,nullptr);
    require(x.bonds.size()==1,"Lewis integer revived a pair rejected by validated connectivity");
    w.bond_order_provenance=cov::DataProvenance::Derived;
    w.bond_orders={{0,1,0.0,cov::DataProvenance::Derived}};
    const auto geometry=cov::analyse_bonds(w);
    require(std::none_of(geometry.begin(),geometry.end(),[](const auto& b){return b.atom_b==1;}),"real zero became geometry fallback");
    require(std::any_of(geometry.begin(),geometry.end(),[](const auto& b){return b.atom_a==0 && b.atom_b==2;}),"missing pair became an implicit zero");
    std::cout<<"nbo_skeleton_policy_smoke ok\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
