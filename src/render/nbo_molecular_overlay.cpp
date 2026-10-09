#include "cov/nbo_molecular_overlay.hpp"
#include "cov/nbo_salc.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>
namespace cov {
MoleculeOverlay make_nbo_molecule_overlay(const NboIntegration& data,const InteractionGraph& graph,
    std::size_t atom_count,const std::vector<std::size_t>& selected_atoms,
    std::optional<std::size_t> selected_structure,AtomScalarMode mode,bool show_indices,bool show_e2,
    const RoutedAnalysis* routed,const Wavefunction* canonical,NboBondDisplayMode bond_mode,
    const NboElectronicSymmetryScope* electronic_scope) {
    MoleculeOverlay out;out.colour_mode=mode;out.selected_atoms=selected_atoms;
    out.bond_mode=bond_mode;
    out.show_bond_indices=show_indices;out.atom_values.resize(atom_count);out.scalar_range=0;
    // Scalar values come from the same complete, capability-checked result as
    // the UI. In particular alpha/beta NPA charges cannot overwrite total NPA.
    if(routed) {
        const auto& values=mode==AtomScalarMode::NaturalCharge?
            routed->total_atomic_charge:routed->atomic_spin;
        if(mode!=AtomScalarMode::Element && values.available() &&
           values.value->size()==atom_count)
            for(std::size_t atom=0;atom<atom_count;++atom) {
                out.atom_values[atom]=(*values.value)[atom];
                out.scalar_range=std::max(out.scalar_range,std::abs(*out.atom_values[atom]));
            }
    }
    struct Connectivity { bool covalent=false,coordination=false; };
    std::map<std::pair<std::size_t,std::size_t>,Connectivity> connectivity;
    for(const auto& e:graph.edges)if(e.strength==InteractionStrength::StrongConnectivity) {
        auto& pair=connectivity[std::minmax(e.atom_a,e.atom_b)];
        pair.covalent|=e.kind==InteractionKind::CovalentConnectivity;
        pair.coordination|=e.kind==InteractionKind::CoordinationContact;
    }
    std::map<std::pair<NboSpin,std::size_t>,const std::string*> raw_types;
    for(const auto& row:data.dataset.orbitals)if(row.id)raw_types[{row.spin,row.id-1}]=&row.kind;
    const auto has_type=[&](const NboStructureEvidence& e,const char* kind) {
        return std::any_of(e.orbitals.begin(),e.orbitals.end(),[&](const auto& ref) {
            if(ref.kind!=NboOrbitalKind::NBO)return false;
            const auto found=raw_types.find({ref.spin,ref.index});
            return found!=raw_types.end() && *found->second==kind;
        });
    };
    std::map<std::vector<std::size_t>,std::size_t> hyperedges;
    using Pair=std::pair<std::size_t,std::size_t>;
    std::map<Pair,std::size_t> bonds;
    std::map<Pair,std::map<NboSpin,unsigned>> lewis_counts;
    std::map<Pair,bool> local_lewis_verified, conflicting_records;
    for(std::size_t i=0;i<data.structure.size();++i){
        const auto& e=data.structure[i];const bool selected=selected_structure && *selected_structure==i;
        if(!routed &&
           ((e.kind=="atomic_charge" && e.spin==NboSpin::Total &&
             mode==AtomScalarMode::NaturalCharge && nbo_capability(data,"charges") &&
             nbo_capability(data,"charges")->available()) ||
            (e.kind=="atomic_spin" && mode==AtomScalarMode::SpinPopulation &&
             nbo_capability(data,"spin") && nbo_capability(data,"spin")->available()))){
            if(e.value && e.atoms.size()==1 && e.atoms[0]<atom_count){
                out.atom_values[e.atoms[0]]=e.value;
                out.scalar_range=std::max(out.scalar_range,std::abs(*e.value));
            }
        }else if(e.kind=="bond" && e.atoms.size()==2){
            if(e.atoms[0]>=atom_count || e.atoms[1]>=atom_count || e.atoms[0]==e.atoms[1])continue;
            const auto key=std::minmax(e.atoms[0],e.atoms[1]);const auto found=connectivity.find(key);
            // A continuous index or membership in a multicentre group does not
            // create a two-centre bond. Conversely, absent total BD counts do
            // not invalidate existing connectivity or validated spin BD records.
            const bool lewis=(e.lewis_bond_count && *e.lewis_bond_count>0) || has_type(e,"BD");
            const bool covalent=found!=connectivity.end() && found->second.covalent;
            const bool coordination=found!=connectivity.end() && found->second.coordination;
            // Default connectivity comes from the validated graph. An isolated
            // Lewis BD cannot restore a pair rejected by a real zero index or
            // create a remote ordinary bond. Its independent explicit view is
            // still available.
            if(!covalent && !coordination && !(bond_mode==NboBondDisplayMode::LewisStructure && lewis))continue;
            MoleculeOverlayBond b;b.atom_a=e.atoms[0];b.atom_b=e.atoms[1];b.evidence_index=i;
            b.continuous_index=e.wiberg;b.selected=selected;
            b.lewis_bond_count=e.lewis_bond_count;b.evidence_indices={i};
            if(!covalent && coordination)b.style=OverlayBondStyle::Coordination;
            if(found!=connectivity.end())for(const auto& edge:graph.edges)
                if(std::minmax(std::size_t(edge.atom_a),std::size_t(edge.atom_b))==key && edge.wiberg_index) {
                    b.continuous_index=edge.wiberg_index;break;
                }
            const auto [it,inserted]=bonds.emplace(key,out.bonds.size());
            if(inserted)out.bonds.push_back(b);
            else {
                auto& original=out.bonds[it->second];original.evidence_indices.push_back(i);
                original.selected|=selected;
                if(selected)original.evidence_index=i;
                if(!original.continuous_index && b.continuous_index)original.continuous_index=b.continuous_index;
                if(original.continuous_index && b.continuous_index &&
                   std::abs(*original.continuous_index-*b.continuous_index)>2e-4)
                    conflicting_records[key]=true;
                if(e.spin==NboSpin::Total && e.lewis_bond_count)original.lewis_bond_count=e.lewis_bond_count;
                else if(!original.lewis_bond_count)original.lewis_bond_count=e.lewis_bond_count;
            }
            if(e.lewis_bond_count) {
                auto& counts=lewis_counts[key];const auto count=counts.find(e.spin);
                if(count!=counts.end() && count->second!=*e.lewis_bond_count)conflicting_records[key]=true;
                counts[e.spin]=*e.lewis_bond_count;
                // A source integer alone is not a stable local multiple-bond
                // certificate. Require the corresponding associated occupied
                // BD orbitals and a validated metric, retaining their original
                // spin identity. No integer is derived from WBI or Mayer.
                unsigned occupied_bd=0;
                for(const auto& ref:e.orbitals) {
                    const auto raw=raw_types.find({ref.spin,ref.index});
                    if(raw==raw_types.end() || *raw->second!="BD")continue;
                    const auto* orbital=nbo_orbital(data,ref);
                    const double full=ref.spin==NboSpin::Total?2.0:1.0;
                    if(orbital && orbital->occupation && *orbital->occupation>=0.95*full &&
                       !orbital->coefficients.empty() && std::abs(orbital->metric_norm2-1.0)<2e-4)
                        ++occupied_bd;
                }
                local_lewis_verified[key]|=occupied_bd==*e.lewis_bond_count && occupied_bd>1;
            }
        }else if(e.kind=="multicentre" && e.atoms.size()>2){
            // 3Cn and 3C* remain inspectable orbitals, not extra chemical bonds.
            if(!has_type(e,"3C"))continue;
            auto atoms=e.atoms;std::sort(atoms.begin(),atoms.end());
            atoms.erase(std::unique(atoms.begin(),atoms.end()),atoms.end());
            if(atoms.size()<3 || atoms.back()>=atom_count)continue;
            const auto [it,inserted]=hyperedges.emplace(atoms,out.multicentre.size());
            if(inserted)out.multicentre.push_back({i,atoms,selected,{i}});
            else {auto& group=out.multicentre[it->second];group.selected|=selected;group.evidence_indices.push_back(i);}
        }else if(e.kind=="donor_acceptor" && (show_e2 || selected) && e.orbitals.size()>=2){
            const auto* donor=nbo_orbital(data,e.orbitals[0]);const auto* acceptor=nbo_orbital(data,e.orbitals[1]);
            if(donor && acceptor && !donor->atoms.empty() && !acceptor->atoms.empty())
                out.relations.push_back({i,donor->atoms,acceptor->atoms,selected});
        }
    }
    // A physical density scope is cached at dataset attachment by the caller.
    // Its operations already passed complete geometry mappings and total/spin
    // density invariance; WBI agreement provides a separate observable check.
    // Without this evidence default connectivity remains a single solid line.
    const bool density_scope=canonical && canonical->atoms.size()==atom_count &&
        electronic_scope && electronic_scope->naming_scope_verified &&
        !electronic_scope->operations.empty();
    for(const auto& [key,index]:bonds) {
        auto& b=out.bonds[index];
        const auto& counts=lewis_counts[key];
        if(!counts.empty()) {
            const unsigned first=counts.begin()->second;
            for(const auto& [spin,count]:counts)if(count!=first)conflicting_records[key]=true;
        }
        b.density_equivalence_checked=density_scope;
        b.lewis_equivalence_conflict=conflicting_records[key];
        if(density_scope && b.continuous_index)for(const auto& op:electronic_scope->operations) {
            if(op.atom_permutation.size()!=atom_count)continue;
            const auto a=op.atom_permutation[key.first],c=op.atom_permutation[key.second];
            if(a>=atom_count || c>=atom_count)continue;
            const auto equivalent=bonds.find(std::minmax(a,c));
            if(equivalent==bonds.end())continue;
            const auto& other=out.bonds[equivalent->second];
            if(other.continuous_index &&
               std::abs(*other.continuous_index-*b.continuous_index)<=2e-4 &&
               other.lewis_bond_count!=b.lewis_bond_count)
                b.lewis_equivalence_conflict=true;
        }
        if(bond_mode==NboBondDisplayMode::LewisStructure && selected_structure &&
            *selected_structure<data.structure.size()) {
            const auto& source=data.structure[*selected_structure];
            // Explicit inspection keeps the selected original spin record's
            // literal integer together with its source identity. A later total
            // record cannot replace it. Default total connectivity is unchanged.
            if(source.kind=="bond" && source.atoms.size()==2 && source.lewis_bond_count &&
                std::pair<std::size_t,std::size_t>(std::minmax(source.atoms[0],source.atoms[1]))==key)
                b.lewis_bond_count=source.lewis_bond_count;
        }
        if(bond_mode==NboBondDisplayMode::LewisStructure)
            b.multiplicity=static_cast<int>(b.lewis_bond_count.value_or(1));
        else if(density_scope && b.continuous_index && local_lewis_verified[key] &&
                !b.lewis_equivalence_conflict && b.style==OverlayBondStyle::Covalent)
            b.multiplicity=static_cast<int>(b.lewis_bond_count.value_or(1));
    }
    out.scalar_range=std::max(out.scalar_range,1e-12);return out;
}
std::string serialize_molecule_overlay_scalars_json(const MoleculeOverlay& overlay,
                                                     const RoutedAnalysis* routed) {
    std::ostringstream out;out<<std::setprecision(17)
        <<"{\"schema\":\"cov.overlay.scalars.v1\",\"mode\":\"";
    const RoutedResult<std::vector<double>>* result=nullptr;
    switch(overlay.colour_mode) {
        case AtomScalarMode::Element:out<<"element";break;
        case AtomScalarMode::NaturalCharge:out<<"natural_charge";
            if(routed)result=&routed->total_atomic_charge;break;
        case AtomScalarMode::SpinPopulation:out<<"spin_population";
            if(routed)result=&routed->atomic_spin;break;
    }
    out<<"\",\"status\":\""<<(result?routed_status_name(result->status):"not_analysed")
       <<"\",\"provider\":\""<<(result?routed_provider_name(result->provider):"none")
       <<"\",\"values\":[";
    for(std::size_t i=0;i<overlay.atom_values.size();++i) {
        if(i)out<<',';
        if(overlay.atom_values[i])out<<*overlay.atom_values[i];else out<<"null";
    }
    out<<"]}";return out.str();
}
}
