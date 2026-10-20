#include "cov/chemistry_route.hpp"
#include "cov/nbo_molecular_overlay.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}
cov::Wavefunction canonical(){
    cov::Wavefunction wf;
    wf.atoms={{"O",8,0,0,0,8},{"O",8,0,0,2.4,8}};
    return wf;
}
cov::NboStructureEvidence value(std::string kind,cov::NboSpin spin,
                                std::size_t atom,double number){
    cov::NboStructureEvidence row;
    row.kind=std::move(kind);row.spin=spin;row.atoms={atom};row.value=number;
    return row;
}
void verify(const cov::Wavefunction& wf,const cov::NboIntegration& data){
    const auto routed=cov::route_chemistry(wf,&data);
    require(routed.total_atomic_charge.available() &&
            routed.total_atomic_charge.provider==cov::RoutedProvider::Nbo,
            "complete independently validated total NPA charges must route");
    require(routed.atomic_spin.available(),"complete NPA spin must route");
    for(auto atom:std::array<std::size_t,2>{0,1}){
        require(std::abs((*routed.total_atomic_charge.value)[atom])<1e-12,
                "alpha/beta charge rows must never overwrite total NPA charge");
        require(std::abs((*routed.atomic_spin.value)[atom]-1.0)<1e-12,
                "spin must retain its independent validated value");
    }
    const auto charge=cov::make_nbo_molecule_overlay(data,*routed.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::NaturalCharge,false,false,&routed);
    const auto spin=cov::make_nbo_molecule_overlay(data,*routed.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::SpinPopulation,false,false,&routed);
    for(auto atom:std::array<std::size_t,2>{0,1}){
        require(charge.atom_values[atom] && std::abs(*charge.atom_values[atom])<1e-12,
                "overlay must show routed total charge including actual zero");
        require(spin.atom_values[atom] && std::abs(*spin.atom_values[atom]-1)<1e-12,
                "overlay must show routed spin independently");
    }
}
void connectivity_overlay_regression(){
    cov::NboIntegration data;
    cov::InteractionGraph graph;
    const auto edge=[&](std::size_t a,std::size_t b,cov::InteractionKind kind){
        cov::InteractionEdge e;e.atom_a=static_cast<std::uint32_t>(a);e.atom_b=static_cast<std::uint32_t>(b);
        e.kind=kind;e.strength=cov::InteractionStrength::StrongConnectivity;graph.edges.push_back(e);
    };
    edge(0,1,cov::InteractionKind::CovalentConnectivity);
    edge(0,1,cov::InteractionKind::MulticentreSupport);
    edge(1,2,cov::InteractionKind::MulticentreSupport);
    cov::NboStructureEvidence bond;bond.kind="bond";bond.atoms={0,1};bond.wiberg=.9;
    data.structure.push_back(bond); // No total BD: preserve existing covalent skeleton.
    bond.atoms={1,2};bond.lewis_bond_count=1;data.structure.push_back(bond);
    bond.atoms={0,3};bond.lewis_bond_count.reset();data.structure.push_back(bond); // remote WBI only
    const auto add=[&](const char* kind,cov::NboSpin spin,const std::vector<std::size_t>& atoms){
        cov::NboOrbital row;row.kind=kind;row.spin=spin;row.id=data.dataset.orbitals.size()+1;
        data.dataset.orbitals.push_back(row);
        cov::NboStructureEvidence e;e.kind=std::string(kind).starts_with("3C")?"multicentre":"bond";
        e.atoms=atoms;e.spin=spin;e.orbitals={{{cov::NboOrbitalKind::NBO},spin,row.id-1}};
        data.structure.push_back(e);
    };
    add("BD",cov::NboSpin::Alpha,{2,3}); // validated spin BD is still bond evidence
    add("BD*",cov::NboSpin::Beta,{0,3}); // antibond alone is not a bond
    add("3C",cov::NboSpin::Alpha,{0,1,2});
    add("3C",cov::NboSpin::Beta,{2,0,1});
    add("3Cn",cov::NboSpin::Total,{0,1,2});
    add("3C*",cov::NboSpin::Total,{0,1,2});
    const auto check=[&](){
        const auto overlay=cov::make_nbo_molecule_overlay(data,graph,4,{},std::nullopt,
            cov::AtomScalarMode::Element,false,false,nullptr,nullptr,
            cov::NboBondDisplayMode::LewisStructure);
        require(overlay.bonds.size()==3,"spin BD or original skeleton lost; remote WBI/antibond added");
        for(const auto& b:overlay.bonds)require(b.style==cov::OverlayBondStyle::Covalent,
            "multicentre support or absent total Lewis count downgraded a covalent bond");
        require(overlay.multicentre.size()==1 && overlay.multicentre[0].evidence_indices.size()==2,
            "3Cn/3C* created bonds or spin records duplicated a spatial hyperedge");
    };
    const auto skeleton=cov::make_nbo_molecule_overlay(data,graph,4,{},std::nullopt,
        cov::AtomScalarMode::Element,false,false);
    require(skeleton.bonds.size()==1 && skeleton.bonds.front().multiplicity==1,
        "default skeleton accepted isolated Lewis/multicentre pairs as ordinary bonds");
    check();std::reverse(graph.edges.begin(),graph.edges.end());
    std::reverse(data.structure.begin(),data.structure.end());check();
}
void validated_bond_route_regression(){
    auto wf=canonical();
    wf.bond_order_provenance=cov::DataProvenance::Derived;
    wf.bond_orders={{0,1,0.002,cov::DataProvenance::Derived}};
    cov::NboIntegration data;
    data.canonical_fingerprint=cov::nbo_canonical_fingerprint(wf);
    data.dataset.association.compatible=true;
    data.capabilities={{"wiberg",cov::NboCapabilityState::Available,"complete",{}}};
    cov::NboStructureEvidence bond;
    bond.kind="bond";bond.atoms={0,1};bond.wiberg=1.1;bond.source.path="validated.nbo";
    data.structure={bond};
    const auto has_strong=[](const auto& route){
        return std::any_of(route.interaction_graph.value->edges.begin(),
            route.interaction_graph.value->edges.end(),[](const auto& edge){
                return edge.strength==cov::InteractionStrength::StrongConnectivity;
            });
    };
    const auto actual=cov::route_chemistry(wf,&data);
    require(has_strong(actual),"validated WBI must participate in connectivity, not just overlay text");
    const auto& edge=actual.interaction_graph.value->edges.front();
    require(edge.wiberg_index && std::abs(*edge.wiberg_index-1.1)<1e-12 &&
            std::abs(edge.mayer_order-.002)<1e-12 &&
            edge.connectivity_method.find("Wiberg")!=std::string::npos,
            "NBO bond route must preserve distinct physical indices and provenance");
    data.capabilities[0].state=cov::NboCapabilityState::Rejected;
    require(!has_strong(cov::route_chemistry(wf,&data)),"rejected WBI enabled false connectivity");
    data.capabilities[0].state=cov::NboCapabilityState::Available;
    data.canonical_fingerprint="wrong-source";
    require(!has_strong(cov::route_chemistry(wf,&data)),"unassociated WBI enabled connectivity");
    data.canonical_fingerprint=cov::nbo_canonical_fingerprint(wf);
    data.structure[0].wiberg=0.0;
    require(!has_strong(cov::route_chemistry(wf,&data)),"actual zero WBI replaced by geometry");
    data.structure[0].wiberg=1.1;
    bond.wiberg=.1;data.structure.push_back(bond);
    require(!has_strong(cov::route_chemistry(wf,&data)),"conflicting duplicate WBI not rejected");
    std::reverse(data.structure.begin(),data.structure.end());
    require(!has_strong(cov::route_chemistry(wf,&data)),"WBI conflict depends on input order");
    require(wf.bond_orders[0].mayer_order==.002,"source Mayer modified by routing");
}
} // namespace

int main(){
    connectivity_overlay_regression();
    validated_bond_route_regression();
    const auto wf=canonical();
    cov::NboIntegration data;
    data.canonical_fingerprint=cov::nbo_canonical_fingerprint(wf);
    data.dataset.association.compatible=true;
    data.capabilities={{"charges",cov::NboCapabilityState::Available,"complete",{}},
                       {"spin",cov::NboCapabilityState::Available,"complete",{}}};
    for(std::size_t atom=0;atom<2;++atom){
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Total,atom,0));
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Alpha,atom,-0.5));
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Beta,atom,0.5));
        data.structure.push_back(value("atomic_spin",cov::NboSpin::Total,atom,1));
    }
    verify(wf,data);
    std::reverse(data.structure.begin(),data.structure.end());verify(wf,data);
    std::rotate(data.structure.begin(),data.structure.begin()+3,data.structure.end());verify(wf,data);
    data.capabilities[1].state=cov::NboCapabilityState::Rejected;
    auto rejected=cov::route_chemistry(wf,&data);
    require(rejected.total_atomic_charge.available(),"rejected spin must not reject total charge");
    require(!rejected.atomic_spin.available() &&
            rejected.atomic_spin.status==cov::RoutedStatus::Rejected,
            "rejected spin must remain rejected");
    const auto overlay=cov::make_nbo_molecule_overlay(data,*rejected.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::SpinPopulation,false,false,&rejected);
    require(!overlay.atom_values[0]&&!overlay.atom_values[1],
            "rejected spin must not render a pseudo-value");
    data.structure.erase(std::find_if(data.structure.begin(),data.structure.end(),[](const auto& row){
        return row.kind=="atomic_charge"&&row.spin==cov::NboSpin::Total&&row.atoms==std::vector<std::size_t>{0};
    }));
    const auto missing=cov::route_chemistry(wf,&data);
    require(!missing.total_atomic_charge.available() &&
            !missing.total_atomic_charge.value,
            "missing one total record is incomplete, never zero");
}
