#pragma once
#include "cov/orbital_ui.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_symmetry_components.hpp"
#include "cov/validation.hpp"
#include <imgui.h>
#include <cmath>
#include <string>

namespace cov::ui {
// A component's target MO is context, not the identity of the displayed field.
inline bool uses_inspection_details(const ActiveOrbitalView& view) {
    if(view.kind==ActiveOrbitalKind::Canonical)return false;
    if(view.selection){
        const auto& s=*view.selection;
        if(s.mode==NboSelectionMode::Orbital && s.terms.size()==1 &&
           s.terms.front().orbital.kind==NboOrbitalKind::Canonical)return false;
    }
    return true;
}
inline bool uses_inspection_details(const OrbitalUIState& state) {
    return state.active_view && uses_inspection_details(*state.active_view);
}
inline std::string inspection_copy_metadata(const ActiveOrbitalView& view) {
    return serialize_active_orbital_view_json(view);
}
inline const char* inspection_text(Language language,const char* en,const char* zh,
                                   const char* ja,const char* fr) {
    switch(language){case Language::ChineseSimplified:return zh;
        case Language::Japanese:return ja;case Language::French:return fr;default:return en;}
}
inline std::string inspection_glyph_seed() {
    return "对称性组成显示成分返回正则轨道返回源 "
        "対称性の組成成分を表示正準軌道に戻る元に戻る Composition de symétrie Afficher Retour MO Retour source "
        "当前检视源轨道带系数分量部分和组合叠加系数成分范数平方归一化目标 "
        "現在の表示元の軌道重み付き成分部分和組合せ重ね表示係数成分ノルムの二乗正規化対象 "
        "Inspection actuelle Orbitales sources Composante pondérée Somme partielle Combinaison "
        "Superposition Coefficient Norme au carré Normalisation Cible";
}
inline void draw_symmetry_composition(const Wavefunction& canonical,const NboAomoName* name,
    const OrbitalUIState& state,Language language,std::optional<std::size_t> canonical_index,
    std::optional<std::size_t> salc_index={},bool interactive=true) {
    if(!name || !name->decomposition_verified || name->verified)return;
    ImGui::TextUnformatted(inspection_text(language,"Symmetry composition","对称性组成",
        "対称性の組成","Composition de symétrie"));
    if(!name->point_group.empty())ImGui::SameLine(),ImGui::TextUnformatted(point_group_display(name->point_group).c_str());
    validation::field("details.symmetry.composition",serialize_orbital_name_json(*name));
    for(std::size_t i=0;i<name->components.size();++i) {
        const auto& part=name->components[i];
        // Only suppress numerical dust in the visible list; every coefficient
        // and weight, including these entries, remains in structured exports.
        if(part.weight<1e-10)continue;
        NboAomoName display;display.verified=true;display.irrep=part.irrep;display.point_group=name->point_group;
        ImGui::PushID(static_cast<int>(i));
        ImGui::Text("%s: %.8g%%",orbital_irrep_display_label(display).c_str(),100*part.weight);
        if(interactive && state.nbo_ui) {
            const auto selected=symmetry_component_selection(canonical,state.nbo_ui->integration,
                state.nbo_ui->aomo.salc_model.get(),*name,i,canonical_index,salc_index);
            if(selected) {
                ImGui::SameLine();
                if(ImGui::SmallButton(inspection_text(language,"Show component","显示成分",
                    "成分を表示","Afficher")))state.nbo_ui->aomo.pending_selection=*selected;
                validation::item("details.symmetry.component."+std::to_string(i));
            }
        }
        ImGui::PopID();
    }
}
inline void draw_inspection_details(const Wavefunction& canonical,OrbitalUIState& state,
                                    Language language) {
    if(!uses_inspection_details(state))return;
    const auto& view=*state.active_view;
    const auto* s=view.selection?&*view.selection:nullptr;
    const auto* integration=state.nbo_ui?state.nbo_ui->integration:nullptr;
    const auto* model=state.nbo_ui?state.nbo_ui->aomo.salc_model.get():nullptr;
    const NboSalcOrbital* side=nullptr;
    const NboAomoName* name=nullptr;
    std::optional<std::size_t> side_index;
    if(s && model && s->dataset_id==model->dataset_id && !s->source_id.empty()){
        for(std::size_t i=0;i<model->orbitals.size();++i)if(model->orbitals[i].id==s->source_id){
            side=&model->orbitals[i];
            side_index=i;
            if(state.nbo_ui->aomo.names && i<state.nbo_ui->aomo.names->salc.size())
                name=&state.nbo_ui->aomo.names->salc[i];
            break;
        }
    }
    const bool single=s && s->mode==NboSelectionMode::Orbital && s->terms.size()==1;
    const auto* descriptor=single && integration?nbo_orbital(*integration,s->terms[0].orbital):nullptr;
    const bool base_side=s && (s->semantic_kind=="salc" ||
        s->semantic_kind=="spin_averaged_spatial_orbital");
    ImGui::TextWrapped("%s",view.label.c_str());
    validation::field("details.active.semantic",view.semantic_kind);
    validation::field("details.active.selection",inspection_copy_metadata(view));
    if(s){
        const char* mode=nullptr;
        switch(s->mode){
        case NboSelectionMode::WeightedComponent:mode=inspection_text(language,"Weighted component","带系数分量","重み付き成分","Composante pondérée");break;
        case NboSelectionMode::PartialSum:mode=inspection_text(language,"Partial sum","部分和","部分和","Somme partielle");break;
        case NboSelectionMode::Combination:mode=inspection_text(language,"Combination","组合","組合せ","Combinaison");break;
        case NboSelectionMode::Overlay:mode=inspection_text(language,"Overlay","叠加","重ね表示","Superposition");break;
        default:break;
        }
        if(mode)ImGui::TextUnformatted(mode);
        if(s->target_canonical_index && *s->target_canonical_index<canonical.orbitals.size()){
            const auto target=*s->target_canonical_index;
            const auto* target_name=state.nbo_ui && state.nbo_ui->aomo.names &&
                target<state.nbo_ui->aomo.names->canonical.size()?&state.nbo_ui->aomo.names->canonical[target]:nullptr;
            ImGui::TextWrapped("%s: %s",inspection_text(language,"Target","目标","対象","Cible"),
                canonical_mo_display_label(canonical,target,target_name).c_str());
        }
    }
    // Display only physical values attached to this object, never a target MO's values.
    std::optional<double> energy,occupation;
    if(base_side && side){energy=side->energy_hartree;occupation=side->occupation;}
    else if(single && descriptor){energy=descriptor->energy_hartree;occupation=descriptor->occupation;}
    if(base_side && s->spatial_spin){energy=s->spatial_spin->energy_hartree;occupation=s->spatial_spin->occupation;}
    if(energy && std::isfinite(*energy))ImGui::Text("%s: %s",tr(Text::Energy,language),
        format_energy(*energy,state.energy_unit,6).c_str());
    if(occupation && std::isfinite(*occupation))ImGui::Text("%s: %.8g",tr(Text::Occupation,language),*occupation);
    if(base_side && name && name->verified && !name->irrep.empty()) {
        ImGui::Text("%s: %s",tr(Text::Symmetry,language),orbital_irrep_display_label(*name).c_str());
        if(!name->point_group.empty())ImGui::TextUnformatted(point_group_display(name->point_group).c_str());
    }
    if(base_side && side_index)draw_symmetry_composition(canonical,name,state,language,{},side_index);
    if(s && is_symmetry_component(*s) && s->target_canonical_index && state.nbo_ui) {
        if(ImGui::Button(inspection_text(language,"Return to source MO","返回正则轨道",
            "正準軌道に戻る","Retour MO")))
            state.nbo_ui->focus.pending_canonical_selection=*s->target_canonical_index;
        validation::item("details.symmetry.return-canonical");
    }
    if(s && s->semantic_kind=="salc_symmetry_component" && side_index && model && state.nbo_ui) {
        if(ImGui::Button(inspection_text(language,"Return to source SALC","返回源 SALC",
            "元の SALC に戻る","Retour SALC source")))
            state.nbo_ui->aomo.pending_selection=nbo_salc_selection(*model,*side_index);
        validation::item("details.symmetry.return-salc");
    }
    if(state.inspection && state.inspection->selection.dataset_id==view.source_id){
        for(std::size_t i=0;i<state.inspection->metric_norm2.size();++i)
            ImGui::Text("%s %zu: %.8g",inspection_text(language,"Component norm squared","成分范数平方",
                "成分ノルムの二乗","Norme au carré"),i+1,state.inspection->metric_norm2[i]);
    }
    if(ImGui::Button((std::string(tr(Text::CopyMetadata,language))+"###details.inspection.copy").c_str())){
        const auto copied=inspection_copy_metadata(view);ImGui::SetClipboardText(copied.c_str());
        const char* clipboard=ImGui::GetClipboardText();
        validation::record("details.copy","{\"text\":"+validation::quote(copied)+
            ",\"clipboard_matches\":"+(clipboard && copied==clipboard?"true":"false")+"}");
    }
    validation::item("details.inspection.copy");
    const bool sources_open=s && ImGui::TreeNode(inspection_text(language,"Source orbitals","源轨道","元の軌道","Orbitales sources"));
    if(s)validation::item("details.inspection.sources-tree");
    if(sources_open){
        ImGuiListClipper clipper;clipper.Begin(static_cast<int>(s->terms.size()),ImGui::GetTextLineHeightWithSpacing());
        while(clipper.Step())for(int i=clipper.DisplayStart;i<clipper.DisplayEnd;++i){
            const auto& t=s->terms[static_cast<std::size_t>(i)];
            const std::string source=t.orbital.kind==NboOrbitalKind::Canonical
                ?canonical_mo_source_label(canonical,t.orbital.index)
                :std::string(nbo_orbital_kind_name(t.orbital.kind))+" "+std::to_string(t.orbital.index+1);
            ImGui::Text("%s / %s; c=%+.8g",source.c_str(),nbo_spin_name(t.orbital.spin),t.coefficient);
        }
        ImGui::TreePop();
    }
    const bool spin_open=s && s->spatial_spin && ImGui::TreeNode("α / β");
    if(s && s->spatial_spin)validation::item("details.inspection.spin-tree");
    if(spin_open){
        for(const auto& channel:s->spatial_spin->channels){
            ImGui::TextUnformatted(nbo_spin_name(channel.spin));
            if(channel.energy_hartree)ImGui::Text("%s: %s",tr(Text::Energy,language),
                format_energy(*channel.energy_hartree,state.energy_unit,6).c_str());
            if(channel.occupation)ImGui::Text("%s: %.8g",tr(Text::Occupation,language),*channel.occupation);
        }
        ImGui::TreePop();
    }
}
} // namespace cov::ui
