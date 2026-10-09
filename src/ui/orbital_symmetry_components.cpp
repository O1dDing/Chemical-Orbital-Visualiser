#include "cov/orbital_symmetry_components.hpp"
#include <algorithm>
#include <cmath>

namespace cov::ui {
namespace {
NboSpin source_spin(const Wavefunction& w,std::size_t index) {
    return w.orbitals[index].spin==Spin::Beta?NboSpin::Beta:
        w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin?NboSpin::Alpha:NboSpin::Total;
}
}
NboIntegration canonical_component_dataset(const Wavefunction& w) {
    NboIntegration out;
    out.canonical_fingerprint=nbo_canonical_fingerprint(w);
    out.id="canonical-components:"+out.canonical_fingerprint;
    for(std::size_t i=0;i<w.orbitals.size();++i) {
        const auto& mo=w.orbitals[i];
        if(mo.coefficients.size()!=w.basis_count)continue;
        NboOrbitalDescriptor d;
        d.ref={NboOrbitalKind::Canonical,source_spin(w,i),i};
        d.id="canonical:"+std::string(nbo_spin_name(d.ref.spin))+":"+std::to_string(i);
        d.label=canonical_mo_source_label(w,i);
        d.coefficients=mo.coefficients;
        d.energy_hartree=mo.energy_hartree;
        if(mo.occupation_provenance!=DataProvenance::Unavailable)d.occupation=mo.occupation;
        d.energy_semantics="unchanged source canonical eigenvalue";
        out.orbitals.push_back(std::move(d));
    }
    return out;
}
bool is_symmetry_component(const NboOrbitalSelection& s) {
    return s.semantic_kind=="canonical_symmetry_component" || s.semantic_kind=="salc_symmetry_component";
}
std::optional<NboOrbitalSelection> symmetry_component_selection(
    const Wavefunction& w,const NboIntegration* data,const NboSalcModel* model,
    const NboAomoName& name,std::size_t component,
    std::optional<std::size_t> canonical_index,std::optional<std::size_t> salc_index) {
    if(!name.decomposition_verified || component>=name.components.size() ||
       canonical_index.has_value()==salc_index.has_value())return {};
    const auto& part=name.components[component];
    if(!(part.weight>0) || !std::isfinite(part.weight) || part.irrep.empty() ||
       part.source_coefficients.size()!=name.component_source_members.size())return {};
    NboOrbitalSelection selection;
    selection.dataset_id=data?data->id:"canonical-components:"+nbo_canonical_fingerprint(w);
    selection.mode=NboSelectionMode::Combination;
    selection.normalize=false;
    selection.group_id="symmetry:"+name.point_group+":"+part.irrep;
    NboAomoName display;display.verified=true;display.irrep=part.irrep;display.point_group=name.point_group;
    const auto irrep=orbital_irrep_display_label(display);
    if(canonical_index) {
        if(*canonical_index>=w.orbitals.size())return {};
        selection.semantic_kind="canonical_symmetry_component";
        selection.target_canonical_index=*canonical_index;
        selection.source_id="symmetry:canonical:"+std::to_string(*canonical_index)+":"+part.irrep;
        selection.label=canonical_mo_source_label(w,*canonical_index)+" · "+irrep;
    } else {
        if(!data || !model || *salc_index>=model->orbitals.size() || model->dataset_id!=data->id)return {};
        selection.semantic_kind="salc_symmetry_component";
        selection.source_id=model->orbitals[*salc_index].id;
        selection.label="SALC "+std::to_string(*salc_index+1)+" · "+irrep;
    }
    const auto append=[&](const NboOrbitalRef& ref,double coefficient) {
        auto found=std::find_if(selection.terms.begin(),selection.terms.end(),[&](const auto& t){return t.orbital==ref;});
        if(found==selection.terms.end())selection.terms.push_back({ref,coefficient});
        else found->coefficient+=coefficient;
    };
    for(std::size_t i=0;i<part.source_coefficients.size();++i) {
        const double coefficient=part.source_coefficients[i];
        if(!std::isfinite(coefficient))return {};
        if(coefficient==0)continue;
        const auto member=name.component_source_members[i];
        if(name.component_source_kind=="canonical") {
            if(member>=w.orbitals.size())return {};
            NboOrbitalRef ref{NboOrbitalKind::Canonical,source_spin(w,member),member};
            if(data && !nbo_orbital(*data,ref))return {};
            append(ref,coefficient);
        } else if(name.component_source_kind=="salc") {
            if(!data || !model || model->dataset_id!=data->id || member>=model->orbitals.size())return {};
            const auto source=nbo_salc_selection(*model,member);
            for(const auto& term:source.terms) {
                if(!nbo_orbital(*data,term.orbital))return {};
                append(term.orbital,coefficient*term.coefficient);
            }
        } else return {};
    }
    selection.terms.erase(std::remove_if(selection.terms.begin(),selection.terms.end(),
        [](const auto& term){return term.coefficient==0;}),selection.terms.end());
    if(selection.terms.empty())return {};
    return selection;
}
} // namespace cov::ui
