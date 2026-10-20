#include "cov/chemistry_route.hpp"
#include "cov/orbital_inspection_ui.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/wavefunction_io.hpp"
#include <imgui.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Compact observation of the production draw path for external case replays.
// This is not a scientific oracle: raw-file checks and ordinary-window checks
// remain separate from these exported layout measurements.
int main(int argc,char** argv) {try {
    if(argc<4)throw std::runtime_error("Usage: cov_nbo_layout_probe canonical.fchk analysis_directory output.json [export_base]");
    // Routing must depend on the displayed object, not its target or first source term.
    cov::ActiveOrbitalView active;active.kind=cov::ActiveOrbitalKind::Inspection;
    active.label="verified display";active.source_label="producer";
    active.display_name_evidence="independent proof";
    cov::NboOrbitalSelection typed;typed.mode=cov::NboSelectionMode::Combination;
    typed.semantic_kind="salc";typed.terms={{{cov::NboOrbitalKind::NAO,cov::NboSpin::Total,0},1}};
    active.selection=typed;
    if(!cov::ui::uses_inspection_details(active))throw std::runtime_error("SALC routed to canonical details");
    active.selection->mode=cov::NboSelectionMode::WeightedComponent;
    active.selection->target_canonical_index=0;
    active.selection->terms[0].orbital.kind=cov::NboOrbitalKind::Canonical;
    if(!cov::ui::uses_inspection_details(active))throw std::runtime_error("Component inherited target details");
    const auto copied=cov::ui::inspection_copy_metadata(active);
    if(copied!=cov::serialize_active_orbital_view_json(active) || copied.find("independent proof")==std::string::npos)
        throw std::runtime_error("Inspection copy lost display identity or selection");
    active.selection->mode=cov::NboSelectionMode::Orbital;
    if(cov::ui::uses_inspection_details(active))throw std::runtime_error("Canonical source lost canonical details");
    active.selection.reset();active.kind=cov::ActiveOrbitalKind::NboSet;
    if(!cov::ui::uses_inspection_details(active))throw std::runtime_error("NBO set inherited canonical details");
    const auto start=std::chrono::steady_clock::now();
    auto wf=cov::parse_wavefunction(argv[1]);
    const auto original=wf;
    const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1)throw std::runtime_error("Expected one analysis candidate");
    const auto integration=cov::read_nbo_integration(wf,found.candidates.front());
    const auto route=cov::route_chemistry(wf,&integration);
    cov::ui::NboAomoUIState state;
    if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf))
        throw std::runtime_error("Unified graph unavailable: "+state.status);
    cov::MODiagramOptions options;
    options.routed=&route;options.routed_identity=integration.id;
    for(std::size_t i=0;i<wf.orbitals.size();++i)
        if(wf.orbitals[i].occupation>0)options.selected_index=i;
    if(state.names)for(const auto& name:state.names->canonical)
        options.display_names.push_back(name.label);
    std::cerr<<"Building central diagram\n";
    const auto data=cov::build_mo_diagram_data(wf,options);
    const auto graph=cov::make_mo_diagram_view_snapshot(data,options,options.selected_index,"case-layout-probe");
    IMGUI_CHECKVERSION();ImGui::CreateContext();
    auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={1280,1800};io.DeltaTime=1.0f/60;
    std::cerr<<"Loading fonts\n";
    if(!cov::ui::configure_fonts())throw std::runtime_error("Font setup failed");
    const auto draw=[&](cov::ui::Language language) {
        ImGui::NewFrame();ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({1200,1700},ImGuiCond_Always);
        ImGui::Begin("case replay",nullptr,ImGuiWindowFlags_NoSavedSettings);
        const bool ok=cov::ui::draw_nbo_aomo_diagram(state,integration,wf,graph,language,1.0f);
        ImGui::End();ImGui::Render();
        if(!ok||!state.drawn_snapshot)throw std::runtime_error("Draw did not provide a snapshot");
    };
    std::cerr<<"Drawing unified snapshot\n";
    draw(cov::ui::Language::English);
    bool preserved=wf.ao_overlap==original.ao_overlap&&wf.orbitals.size()==original.orbitals.size();
    for(std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& a=wf.orbitals[i];const auto& b=original.orbitals[i];
        preserved=preserved&&a.coefficients==b.coefficients&&a.gaussian_source_coefficients==b.gaussian_source_coefficients
            &&a.energy_hartree==b.energy_hartree&&a.occupation==b.occupation&&a.spin==b.spin
            &&a.source_orbital_index==b.source_orbital_index;
    }
    const auto& view=*state.drawn_snapshot;
    std::filesystem::path output(argv[3]);
    if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output);out<<std::setprecision(17)<<std::boolalpha;
    out<<"{\"canonical_preserved\":"<<preserved<<",\"canonical_count\":"<<wf.orbitals.size()
       <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
       <<",\"width\":"<<view.canvas_width<<",\"height\":"<<view.canvas_height
       <<",\"numeric_span\":"<<view.numeric_span<<",\"hidden_class_count\":"<<view.hidden_class_count
       <<",\"energies\":[";
    bool comma=false;
    for(const auto& e:state.salc_model->energies) {
        if(comma)out<<',';comma=true;
        out<<"{\"spin\":"<<std::quoted(cov::nbo_spin_name(e.spin))<<",\"status\":"<<std::quoted(e.status)
           <<",\"available\":"<<e.available<<",\"canonical_same_operator\":"<<e.canonical_same_operator
           <<",\"printed_operator_verified\":"<<e.printed_operator_verified<<'}';
    }
    out<<"],\"nodes\":[";comma=false;
    for(const auto& n:view.nodes) {
        if(comma)out<<',';comma=true;
        out<<"{\"id\":"<<std::quoted(n.id)<<",\"label\":"<<std::quoted(n.label)
           <<",\"lane\":"<<static_cast<int>(n.lane)<<",\"x\":"<<n.x<<",\"y\":"<<n.y
           <<",\"width\":"<<n.width<<",\"height\":"<<n.height
           <<",\"label_x\":"<<n.label_x<<",\"label_y\":"<<n.label_y
           <<",\"label_width\":"<<n.label_width<<",\"label_height\":"<<n.label_height
           <<",\"quantitative\":"<<n.quantitative_energy<<",\"group_header\":"<<n.group_header
           <<",\"group\":"<<std::quoted(n.display_group_id)<<",\"members\":"<<n.shell_member_count
           <<",\"shell_label\":"<<std::quoted(n.shell_label)<<",\"offset\":"<<n.display_offset_y
           <<",\"energy\":";
        if(n.energy_hartree)out<<*n.energy_hartree;else out<<"null";
        out<<",\"display_energy\":";
        if(n.display_energy_hartree)out<<*n.display_energy_hartree;else out<<"null";
        out<<",\"salc_index\":";
        if(n.salc_index)out<<*n.salc_index;else out<<"null";
        out<<'}';
    }
    out<<"],\"captions\":[";comma=false;
    for(const auto& c:view.captions){if(comma)out<<',';comma=true;out<<std::quoted(c.text);}
    out<<"]}";out.close();if(!out)throw std::runtime_error("Output write failed");
    if(argc>4) {
        const auto simplified=state.salc_model;
        if(simplified->spin_averaged){
            const auto side=std::find_if(simplified->orbitals.begin(),simplified->orbitals.end(),[](const auto& o){return o.spatial_spin.has_value();});
            state.selection=cov::nbo_salc_selection(*simplified,std::distance(simplified->orbitals.begin(),side));
            state.focused_canonical_index=options.selected_index;
            state.preset=cov::ui::NboAomoPreset::Research;
            if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf)||state.salc_model->spin_averaged||state.selection||!state.pending_selection||state.pending_selection->semantic_kind!="canonical")
                throw std::runtime_error("Detailed spin switch left a stale spatial selection");
            state.selection=state.pending_selection;state.pending_selection.reset();draw(cov::ui::Language::English);
            const auto detailed=cov::ui::export_nbo_aomo_bundle(*state.drawn_snapshot,integration,std::string(argv[4])+"-detailed");
            if(!detailed.json||!detailed.csv||!detailed.png||!detailed.svg)throw std::runtime_error("Detailed export failed: "+detailed.error);
            state.preset=cov::ui::NboAomoPreset::Teaching;
            if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf)||state.salc_model!=simplified)
                throw std::runtime_error("Simplified mode failed to restore immutable cached model");
            state.selection=cov::nbo_salc_selection(*simplified,std::distance(simplified->orbitals.begin(),side));
            if(!cov::make_nbo_selection_view(integration,wf,*state.selection).available)
                throw std::runtime_error("Restored spatial selection unavailable");
            state.pending_selection=state.selection;state.selection.reset();state.preset=cov::ui::NboAomoPreset::Research;
            if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf)||!state.pending_selection||state.pending_selection->spatial_spin)
                throw std::runtime_error("Queued spatial selection survived spin-mode change");
            state.pending_selection=cov::nbo_salc_selection(*state.salc_model,0);state.selection=state.pending_selection;
            state.selected_side_node_id="previous-model";
            state.salc_model=std::make_shared<const cov::NboSalcModel>(*state.source_salc_model);
            if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf)||state.selection||!state.selected_side_node_id.empty()||!state.pending_selection||state.pending_selection->semantic_kind!="canonical")
                throw std::runtime_error("Explicit side-model replacement retained old identities");
            state.preset=cov::ui::NboAomoPreset::Teaching;
            if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf))throw std::runtime_error("Mode restoration failed");
            state.pending_selection.reset();state.selection=cov::nbo_salc_selection(*state.salc_model,0);
        }
        const cov::ui::Language languages[]={cov::ui::Language::English,cov::ui::Language::ChineseSimplified,
            cov::ui::Language::Japanese,cov::ui::Language::French};
        for(std::size_t i=0;i<4;++i) {
            draw(languages[i]);
            const auto exported=cov::ui::export_nbo_aomo_bundle(*state.drawn_snapshot,integration,
                std::filesystem::path(std::string(argv[4])+"-"+std::to_string(i)));
            if(!exported.svg||!exported.png||!exported.json||!exported.csv)
                throw std::runtime_error("Bundle export failed: "+exported.error);
        }
        // Rejected/absent attachment must clear cached models and selections.
        auto rejected=integration;for(auto& c:rejected.capabilities)if(c.key=="aomo")c.state=cov::NboCapabilityState::Rejected;
        if(cov::ui::prepare_nbo_aomo_state(state,rejected,wf)||state.selection||state.salc_model||state.drawn_snapshot)
            throw std::runtime_error("Rejected attachment retained stale derived state");
        cov::NboIntegration absent;absent.id="no-nbo";
        if(cov::ui::prepare_nbo_aomo_state(state,absent,wf)||state.pending_selection||state.source_salc_model)
            throw std::runtime_error("No-NBO transition retained source state");
        if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf)||state.salc_model->spin_averaged!=simplified->spin_averaged||state.selection)
            throw std::runtime_error("Reattachment failed to establish fresh mode and identity");
        std::cerr<<"Mode, selection, rejected/no-NBO and reattachment checks passed\n";
    }
    ImGui::DestroyContext();return preserved?0:2;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
