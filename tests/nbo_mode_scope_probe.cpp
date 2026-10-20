#include "cov/nbo_ui.hpp"
#include "cov/orbital_ui.hpp"
#include "cov/wavefunction_io.hpp"
#include <imgui.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

namespace {
void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
}

// Replay the actual parent/sidebar draw path, including its central cache.
// No second implementation of the display filter is used by this observer.
int main(int argc,char** argv) {try {
    require(argc==4,"Usage: cov_nbo_mode_scope_probe canonical.fchk analysis_directory output.json");
    const auto wf=cov::parse_wavefunction(argv[1]);
    const auto identity=cov::nbo_canonical_fingerprint(wf);
    const auto found=cov::discover_nbo_inputs({argv[2]});
    require(found.candidates.size()==1,"Expected one NBO analysis");
    const auto integration=cov::read_nbo_integration(wf,found.candidates.front());
    const auto route=cov::route_chemistry(wf,&integration);
    cov::ui::NboUIState nbo;nbo.integration=&integration;nbo.routed=&route;
    cov::ui::OrbitalUIState state;state.nbo_ui=&nbo;
    cov::ui::OrbitalUIActions actions;
    require(cov::ui::prepare_nbo_aomo_state(nbo.aomo,integration,wf),"AOMO unavailable");
    IMGUI_CHECKVERSION();ImGui::CreateContext();
    auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={1600,1000};io.DeltaTime=1.f/60;
    require(cov::ui::configure_fonts(),"Font setup failed");
    std::size_t selected=0;
    for(std::size_t i=0;i<wf.orbitals.size();++i)if(wf.orbitals[i].occupation>0)selected=i;
    const auto draw=[&]() {
        ImGui::NewFrame();ImGui::SetNextWindowPos({0,0});ImGui::SetNextWindowSize({1550,960});
        ImGui::Begin("scope replay",nullptr,ImGuiWindowFlags_NoSavedSettings);
        cov::ui::draw_energy_diagram(wf,selected,state,cov::ui::Language::ChineseSimplified,1,actions);
        ImGui::End();ImGui::Render();
        require(bool(nbo.aomo.drawn_snapshot),"No final display snapshot");
        require(bool(actions.drawn_diagram),"No central snapshot");
    };
    std::ofstream out(argv[3]);out<<"{\"frames\":[";bool comma=false;
    using Preset=cov::ui::NboAomoPreset;
    std::map<Preset,std::set<std::string>> identities;
    for(const auto preset:{Preset::Teaching,Preset::Research,Preset::Full,Preset::Research,Preset::Teaching}) {
        cov::ui::apply_nbo_aomo_preset(nbo.aomo,preset);draw();draw();
        const auto frozen_view=nbo.aomo.drawn_snapshot;
        const auto& view=*frozen_view;
        std::set<std::string> ids;std::set<std::size_t> central;
        std::size_t side=0,visible_edges=0;
        for(const auto& node:view.nodes) {
            require(ids.insert(node.id).second,"Duplicate displayed object identity");
            if(node.canonical_index)central.insert(*node.canonical_index);else if(!node.group_header)++side;
        }
        for(const auto& edge:view.edges)if(edge.visible)++visible_edges;
        if(identities.contains(preset))require(identities[preset]==ids,"Preset depends on preceding mode");
        else identities[preset]=ids;
        if(preset==Preset::Full)require(central.size()==wf.orbitals.size(),"Full mode omits canonical members");
        if(comma)out<<',';comma=true;
        out<<"{\"preset\":"<<int(preset)<<",\"central_count\":"<<central.size()<<",\"side_count\":"<<side
           <<",\"visible_edges\":"<<visible_edges<<",\"central_members\":[";
        bool member_comma=false;for(auto index:central){if(member_comma)out<<',';member_comma=true;out<<index;}
        out<<"],\"side_groups\":[";bool first=true;
        for(const auto& node:view.nodes)if(!node.canonical_index&&!node.group_header) {
            if(!first)out<<',';first=false;
            out<<"{\"id\":"<<std::quoted(node.id)<<",\"label\":"<<std::quoted(node.label)
               <<",\"group\":"<<std::quoted(node.display_group_id)<<",\"energy\":";
            if(node.energy_hartree)out<<*node.energy_hartree;else out<<"null";
            out<<'}';
        }
        out<<"],\"default_dash_segments\":"<<nbo.aomo.last_drawn_dash_segments
           <<",\"unclipped_dash_segments\":"<<nbo.aomo.last_unclipped_dash_segments
           <<",\"connection_draw_ms\":"<<nbo.aomo.last_connection_draw_ms;
        if(preset==Preset::Full) {
            nbo.aomo.all_connections=true;++nbo.aomo.revision;draw();draw();
            const auto& expanded=*nbo.aomo.drawn_snapshot;
            require(expanded.central_mo_indices==view.central_mo_indices,"All connections changed source membership");
            std::size_t all=0;for(const auto& edge:expanded.edges) {
                if(edge.visible)++all;
                if(std::abs(edge.coefficient)>expanded.numerical_zero_bound)
                    require(edge.visible,"All-connections control omitted an actual nonzero link");
            }
            out<<",\"all_connections_count\":"<<all<<",\"all_connections_dash_segments\":"<<nbo.aomo.last_drawn_dash_segments
               <<",\"all_connections_unclipped_segments\":"<<nbo.aomo.last_unclipped_dash_segments
               <<",\"all_connections_draw_ms\":"<<nbo.aomo.last_connection_draw_ms;
            nbo.aomo.all_connections=false;++nbo.aomo.revision;draw();draw();
        }
        out<<'}';
        const auto before=central;selected=selected?selected-1:1;draw();
        std::set<std::size_t> after;
        for(const auto& node:nbo.aomo.drawn_snapshot->nodes)if(node.canonical_index)after.insert(*node.canonical_index);
        require(before==after,"Inspecting another MO changed central filter membership");
    }
    // Rendering settings never change scientific identities or numerical
    // values. Exercise Gaussian AO and back through the same parent draw path.
    cov::ui::apply_nbo_aomo_preset(nbo.aomo,Preset::Teaching);draw();draw();
    const auto baseline=*nbo.aomo.drawn_snapshot;
    nbo.aomo.show_atom_numbers=!nbo.aomo.show_atom_numbers;
    nbo.aomo.show_fragment_numbers=!nbo.aomo.show_fragment_numbers;
    nbo.aomo.number_ignore_h=!nbo.aomo.number_ignore_h;++nbo.aomo.revision;draw();
    require(baseline.central_mo_indices==nbo.aomo.drawn_snapshot->central_mo_indices,"Number toggle changed source membership");
    for(const auto& before:baseline.nodes) {
        const auto found=std::find_if(nbo.aomo.drawn_snapshot->nodes.begin(),nbo.aomo.drawn_snapshot->nodes.end(),
            [&](const auto& after){return after.id==before.id;});
        require(found!=nbo.aomo.drawn_snapshot->nodes.end(),"Number toggle changed object identity");
        require(found->energy_hartree==before.energy_hartree && found->occupation==before.occupation,"Number toggle changed scientific values");
    }
    nbo.aomo.basis_kind=cov::NboOrbitalKind::GaussianAO;++nbo.aomo.revision;draw();draw();
    require(nbo.aomo.drawn_snapshot->central_mo_indices==baseline.central_mo_indices,"Gaussian AO changed central filter");
    require(std::any_of(nbo.aomo.drawn_snapshot->nodes.begin(),nbo.aomo.drawn_snapshot->nodes.end(),
        [](const auto& node){return node.lane!=cov::ui::NboAomoLane::Centre;}),"Gaussian AO sides absent");
    nbo.aomo.basis_kind=cov::NboOrbitalKind::NAO;++nbo.aomo.revision;draw();draw();
    require(nbo.aomo.drawn_snapshot->central_mo_indices==baseline.central_mo_indices,"NAO restoration changed filter");

    // Renderer-only fixture: a verified weak near-energy relation. Scientific
    // acceptance of weak classification is tested independently by pi probes.
    auto fixture_data=actions.drawn_diagram->data;
    auto degenerate=std::find_if(fixture_data.levels.begin(),fixture_data.levels.end(),[&](const auto& level){
        return level.member_indices.size()>=2 && level.member_spin_counterparts.empty() &&
            std::abs(wf.orbitals[level.member_indices[0]].energy_hartree-wf.orbitals[level.member_indices[1]].energy_hartree)<1e-5;});
    bool weak_fixture=false;
    if(degenerate!=fixture_data.levels.end()) {
        const auto a=degenerate->member_indices[0],b=degenerate->member_indices[1];
        fixture_data.pi_interactions.clear();
        for(auto& level:fixture_data.levels){level.homo=level.lumo=false;level.sigma_fraction=0;}
        cov::PiInteractionDescriptor relation;relation.lower_orbitals={a};relation.upper_orbitals={b};
        auto evidence=std::make_shared<cov::PiPartnerAssessment>();evidence->accepted=true;evidence->weak=true;
        evidence->channel.display_calibration.validated=true;evidence->channel.display_calibration.energy_budget_ev=.025;
        relation.orbital_evidence=evidence;fixture_data.pi_interactions.push_back(relation);
        const cov::MODiagramViewSnapshot fixture{fixture_data,actions.drawn_diagram->options};
        const auto draw_fixture=[&](){ImGui::NewFrame();ImGui::Begin("weak container fixture");
            cov::ui::draw_nbo_aomo_diagram(nbo.aomo,integration,wf,fixture,cov::ui::Language::ChineseSimplified,1);
            ImGui::End();ImGui::Render();};
        draw_fixture();
        const auto& nodes=nbo.aomo.drawn_snapshot->nodes;
        const auto folded=std::find_if(nodes.begin(),nodes.end(),[](const auto& node){return node.weak_display_container;});
        require(folded!=nodes.end(),"Weak near-energy display container absent");
        require(folded->member_canonical_indices==std::vector<std::size_t>{a,b},"Weak container lost source members");
        require(folded->occupation==wf.orbitals[a].occupation+wf.orbitals[b].occupation,"Weak container lost occupation");
        require(!folded->occupation_on_bar&&folded->occupation_label.empty(),"Folded electrons still drawn");
        const auto id=folded->id;nbo.aomo.expanded_weak_groups.insert(id);++nbo.aomo.revision;draw_fixture();
        for(auto member:{a,b})require(std::any_of(nbo.aomo.drawn_snapshot->nodes.begin(),nbo.aomo.drawn_snapshot->nodes.end(),
            [&](const auto& node){return node.canonical_index==member;}),"Expanding weak container did not restore actual source");
        weak_fixture=true;
    }
    require(cov::nbo_canonical_fingerprint(wf)==identity,"Source wavefunction changed");
    out<<"],\"source_preserved\":true,\"mode_order_invariant\":true,\"numbers_and_gaussian_ao_checked\":true,\"weak_renderer_fixture\":"<<(weak_fixture?"true":"false")<<"}";
    ImGui::DestroyContext();std::cout<<"mode scope replay passed\n";return 0;
} catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
