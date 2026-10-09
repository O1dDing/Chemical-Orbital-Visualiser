#include "cov/mo_diagram.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_ui.hpp"

#include <imgui.h>

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace {
void require(bool ok, const char* message) {
    if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}
std::string read(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);
    return {std::istreambuf_iterator<char>(in),{}};
}
std::size_t count(const std::string& value,const std::string& needle) {
    std::size_t result=0,pos=0;
    while ((pos=value.find(needle,pos))!=std::string::npos) { ++result; pos+=needle.size(); }
    return result;
}
}

int main() {
    // A display contract fixture, not a physical electronic-state reference.
    cov::MODiagramData data;
    const double energy[]={-0.20,-0.19999,-0.18,-0.17999,-3.0};
    for (std::size_t i=0;i<5;++i) {
        cov::OrbitalMetadata raw;
        raw.orbital_index=i; raw.raw_mo_number=i+1;
        raw.energy_hartree=energy[i]; raw.occupation=1;
        raw.spin=i>=2?cov::Spin::Beta:cov::Spin::Alpha;
        data.metadata.push_back(raw);
    }
    data.annotations.resize(5);
    data.selection.included_indices={0,1};
    cov::MODiagramLevel level;
    level.metadata=data.metadata[0]; level.metadata.degeneracy_size=2;
    level.member_indices={0,1}; level.member_spin_counterparts={2,3};
    level.member_electrons={{1,1},{1,0}};
    level.layout_energy_hartree=-0.199995;
    level.energy_spread_hartree=0.00001;
    data.levels={level};
    data.energy_transform=cov::build_energy_transform({-0.4,0.1},cov::EnergyAxisMode::Linear);
    cov::MODiagramOptions options;
    options.selected_index=5; options.energy_axis_mode=cov::EnergyAxisMode::Linear;
    options.energy_unit=cov::EnergyUnit::ElectronVolt;
    options.include_hidden_in_metadata=false;
    options.display_names={"1a₁ [alpha]","1b₂ [alpha]","1a₁ [beta]","1b₂ [beta]"};
    options.figure_title="分子轨道能级";
    options.axis_title="能量";
    auto snapshot=cov::make_mo_diagram_view_snapshot(data,options,3,"interactive-canvas");
    const auto members=cov::mo_diagram_member_views(snapshot.data,snapshot.data.levels[0]);
    require(members.size()==2 && !members[0].selected && members[1].selected &&
        members[1].inspected_orbital_index==3 && members[1].orbital_index==1,
        "beta inspection must select exactly its spatial member while retaining both identities");
    require(!snapshot.data.metadata[1].selected && snapshot.data.metadata[3].selected &&
        snapshot.data.metadata[3].energy_hartree==energy[3] &&
        snapshot.options.selected_index==5,"inspection must not alter raw metadata or the row-selection anchor");
    data.metadata[3].energy_hartree=99;
    data.levels.clear(); options.energy_unit=cov::EnergyUnit::Hartree;
    require(snapshot.data.metadata[3].energy_hartree==energy[3] &&
        snapshot.data.levels.size()==1 && snapshot.options.energy_unit==cov::EnergyUnit::ElectronVolt,
        "an export snapshot must own its frame values despite subsequent UI/cache mutations");

    const auto root=std::filesystem::temp_directory_path()/
        ("cov_snapshot_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto result=cov::export_mo_diagram_bundle(snapshot,root/"captured");
    require(result.svg&&result.png&&result.json&&result.csv,"snapshot export failed");
    const auto svg=read(result.svg_path),png=read(result.png_path),json=read(result.json_path),csv=read(result.csv_path);
    require(svg.find("1a₁ [alpha]")!=std::string::npos &&
        svg.find("1b₂ [alpha]")!=std::string::npos && svg.find("分子轨道能级")!=std::string::npos &&
        svg.find("能量 (eV)")!=std::string::npos && svg.find("support (not probability)")==std::string::npos,
        "figure must preserve supplied orbital names, language and axis units without scores");
    // Picture-only export neither creates analysis files nor touches old ones.
    const auto pictures=cov::export_mo_diagram_bundle(snapshot,root/"pictures",cov::DiagramExportContent::Images);
    require(pictures.svg&&pictures.png&&!pictures.json&&!pictures.csv&&
        !std::filesystem::exists(pictures.json_path)&&!std::filesystem::exists(pictures.csv_path),
        "picture export must produce only PNG and SVG");
    require(read(pictures.svg_path)==svg&&read(pictures.png_path)==png,
        "picture-only export must preserve the exact rendered figure");
    std::filesystem::create_directory(pictures.json_path);
    {std::ofstream out(pictures.csv_path);out<<"untouched analysis";}
    const auto images_again=cov::export_mo_diagram_bundle(snapshot,root/"pictures",cov::DiagramExportContent::Images);
    require(images_again.svg&&images_again.png&&std::filesystem::is_directory(pictures.json_path)&&
        read(pictures.csv_path)=="untouched analysis","unrequested analysis paths must be ignored");
    const auto analysis=cov::export_mo_diagram_bundle(snapshot,root/"analysis",cov::DiagramExportContent::AnalysisData);
    require(analysis.json&&analysis.csv&&!analysis.svg&&!analysis.png&&
        !std::filesystem::exists(analysis.svg_path)&&!std::filesystem::exists(analysis.png_path)&&
        read(analysis.json_path)==json&&read(analysis.csv_path)==csv,
        "explicit analysis export must preserve numerical output without rendering images");
    const auto& id=snapshot.data.view->id;
    require(svg.find(id)!=std::string::npos && png.find(id)!=std::string::npos &&
        json.find(id)!=std::string::npos && csv.find(id)!=std::string::npos,
        "all four files must identify the same snapshot");
    require(count(svg,"class=\"mo-member\"")==2 && count(svg,"data-selected=\"true\"")==1 &&
        svg.find("data-inspected-orbital-index=\"3\"")!=std::string::npos,
        "SVG must preserve real member count and beta inspection");
    require(json.find("\"spin_counterpart\":3")!=std::string::npos &&
        json.find("\"index\": 3")!=std::string::npos &&
        json.find("\"index\": 4")==std::string::npos &&
        json.find("\"all_members_energy_min_hartree\":-0.2")!=std::string::npos &&
        json.find("\"all_members_energy_max_hartree\":-0.17999")!=std::string::npos,
        "compact machine output must retain counterpart metadata and the actual energy range");
    auto hidden=cov::make_mo_diagram_view_snapshot(snapshot.data,snapshot.options,4);
    require(!hidden.data.levels[0].metadata.selected && hidden.data.metadata[4].selected,
        "inspecting a hidden MO must not manufacture a visible row");
    require(cov::write_mo_diagram_json(hidden.data,hidden.options,root/"hidden.json"),"hidden export failed");
    require(read(root/"hidden.json").find("\"index\": 4")!=std::string::npos,
        "the inspected hidden MO must retain raw metadata even with compact metadata export");
    auto invalid=cov::make_mo_diagram_view_snapshot(snapshot.data,snapshot.options,999);
    require(!invalid.data.view->inspected_orbital_index,"invalid inspection must not become a false index");
    auto raw=level; raw.member_indices.clear(); raw.member_spin_counterparts.clear();
    raw.metadata.degeneracy_size=3;
    require(cov::mo_diagram_member_views(snapshot.data,raw).size()==1,
        "a degeneracy-size field alone must not invent canonical members");

    IMGUI_CHECKVERSION(); ImGui::CreateContext();
    auto& io=ImGui::GetIO(); io.IniFilename=nullptr;
    io.DisplaySize={1280,1800}; io.DeltaTime=1.0f/60;
    io.Fonts->AddFontDefault(); io.Fonts->Build();
    cov::Wavefunction wf; wf.atoms.resize(1); wf.atoms[0].atomic_number=1;
    for (int i=0;i<8;++i) {
        cov::MolecularOrbital mo;
        mo.energy_hartree=-0.5+0.11*(i%4); mo.occupation=i%4<2?1.0f:0.0f;
        mo.spin=i>=4?cov::Spin::Beta:cov::Spin::Alpha;
        wf.orbitals.push_back(mo);
    }
    cov::ui::OrbitalUIState state;
    cov::ui::OrbitalUIActions actions;
    ImVec2 linear_target;
    auto frame=[&](std::size_t selected) {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({1200,1700},ImGuiCond_Always);
        ImGui::Begin("snapshot controls",nullptr,ImGuiWindowFlags_NoSavedSettings);
        const auto start=ImGui::GetCursorScreenPos();
        // The scope combo occupies one actual ImGui control row before the
        // energy radios. Derive its vertical advance from the current style.
        linear_target={start.x+ImGui::CalcTextSize(cov::ui::tr(cov::ui::Text::EnergyScale,
            cov::ui::Language::English)).x+ImGui::GetStyle().ItemSpacing.x+8,
            start.y+ImGui::GetFrameHeightWithSpacing()+ImGui::GetFrameHeight()*0.5f};
        cov::ui::draw_energy_diagram(wf,selected,state,cov::ui::Language::English,1.0f,actions);
        ImGui::End(); ImGui::Render();
        require(actions.drawn_diagram!=nullptr,"live UI must provide its drawn snapshot");
    };
    frame(0); const auto first=actions.drawn_diagram;
    frame(0); require(first==actions.drawn_diagram,"unchanged view should reuse its immutable snapshot");
    frame(7); const auto beta=actions.drawn_diagram;
    require(beta!=first && beta->data.view->inspected_orbital_index==7 &&
        first->data.view->inspected_orbital_index==0,"frame selection must not mutate earlier snapshots");
    // Exercise the actual ImGui radio button. The same click frame must have
    // consistent widget mode, transform and export options.
    io.AddMousePosEvent(linear_target.x,linear_target.y); frame(7);
    io.AddMouseButtonEvent(0,true); frame(7);
    io.AddMouseButtonEvent(0,false); frame(7);
    require(state.energy_axis_mode==cov::EnergyAxisMode::Linear &&
        actions.drawn_diagram->options.energy_axis_mode==state.energy_axis_mode &&
        actions.drawn_diagram->data.energy_transform.mode==state.energy_axis_mode,
        "actual axis click must freeze the new transform in the same drawn frame");
    state.energy_unit=cov::EnergyUnit::ElectronVolt; frame(7);
    require(actions.drawn_diagram->options.energy_unit==cov::EnergyUnit::ElectronVolt,
        "energy-unit change must produce a matching export snapshot");

    // Open-shell alpha/beta rows may have the same producer AO/NAO number.
    // The live graph and its exports must show distinct, legible spin labels.
    wf.atoms[0].symbol="Cr";wf.atoms[0].atomic_number=24;
    cov::NboIntegration integration;
    integration.id="open-shell-heavy-atom-display-fixture";
    integration.capabilities.push_back({"aomo",cov::NboCapabilityState::Available,
        "Display contract fixture"});
    for(const auto kind:{cov::NboOrbitalKind::NAO,cov::NboOrbitalKind::GaussianAO})
        for(const auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}) {
            cov::NboOrbitalDescriptor orbital;
            orbital.ref={kind,spin,0};
            orbital.id=std::string(cov::nbo_orbital_kind_name(kind))+":"+
                cov::nbo_spin_name(spin)+":0";
            // The visible element must come from wf, not this opaque producer
            // label (nor an atomic-number placeholder such as Z241).
            orbital.label="3d Val";
            orbital.atoms={0};
            integration.orbitals.push_back(orbital);
            cov::NboMoLink link;
            link.orbital=orbital.ref;
            link.canonical_index=spin==cov::NboSpin::Alpha?0:2;
            link.coefficient=0.5;
            integration.links.push_back(link);
        }
    for(const auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}) {
        cov::NboNao nao;
        nao.id=1;nao.atom=1;nao.symbol="Cr";nao.angular="dxy";nao.spin=spin;
        integration.dataset.naos.push_back(nao);
    }
    cov::ui::NboAomoUIState aomo_state;
    auto draw_aomo=[&](cov::NboOrbitalKind kind,const cov::MODiagramViewSnapshot& graph) {
        aomo_state.basis_kind=kind;
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({1200,1700},ImGuiCond_Always);
        ImGui::Begin("open-shell labels",nullptr,ImGuiWindowFlags_NoSavedSettings);
        const bool drawn=cov::ui::draw_nbo_aomo_diagram(aomo_state,integration,wf,graph,
            cov::ui::Language::English,1.0f);
        ImGui::End();ImGui::Render();
        require(drawn && aomo_state.drawn_snapshot,
            "display fixture must draw the whole AO/NAO graph");
        return aomo_state.drawn_snapshot;
    };
    auto aomo_frame=[&](cov::NboOrbitalKind kind) {
        const auto drawn=draw_aomo(kind,snapshot);
        const auto& view=*drawn;
        const std::string stem=kind==cov::NboOrbitalKind::NAO?"Cr1 dxy / NAO 1":"Cr1 / AO 1";
        // Expanded orbitals carry their atom labels themselves. An atom
        // disclosure header belongs only to the folded state of the new graph.
        for(const auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}) {
            const auto label=stem+" ["+cov::nbo_spin_name(spin)+"]";
            bool found=false;
            for(const auto& node:view.nodes)
                if(node.orbital && node.orbital->kind==kind && node.orbital->spin==spin &&
                   node.label==label)found=true;
            require(found,"live open-shell AO/NAO node must identify its spin");
        }
        const auto base=root/(kind==cov::NboOrbitalKind::NAO?"nao-spin":"ao-spin");
        const auto exported=cov::ui::export_nbo_aomo_bundle(view,integration,base);
        require(exported.json && exported.svg && exported.png,
            "open-shell AO/NAO graph export failed");
        const auto json=read(exported.json_path),svg=read(exported.svg_path);
        require(json.find(stem+" [alpha]")!=std::string::npos &&
            json.find(stem+" [beta]")!=std::string::npos &&
            svg.find(stem+" [alpha]")!=std::string::npos &&
            svg.find(stem+" [beta]")!=std::string::npos &&
            json.find("Z241")==std::string::npos && svg.find("Z241")==std::string::npos,
            "JSON and SVG must preserve the heavy-atom symbol and both visible spin labels");
    };
    aomo_frame(cov::NboOrbitalKind::NAO);
    aomo_frame(cov::NboOrbitalKind::GaussianAO);
    // Name evidence is shared once, not copied for every coefficient edge.
    // This fixture includes a hidden canonical target and both spin channels.
    {
        auto compact_view=*draw_aomo(cov::NboOrbitalKind::NAO,snapshot);
        auto names=std::make_shared<cov::ui::NboAomoNames>();
        names->canonical.resize(3);
        for(auto& name:names->canonical){
            name.irrep="A";name.label="source name";
            cov::ui::NboIrrepComponent component;component.irrep="A";
            component.source_coefficients.assign(4096,0.12345678901234567);
            name.components.push_back(std::move(component));
        }
        compact_view.names=names;
        auto dense=integration;
        for(std::size_t i=0;i<50;++i){auto link=integration.links[i%integration.links.size()];
            link.coefficient=-0.375;link.source.path="producer,quoted\"source";
            link.source.line_begin=17;link.source.block="retained block";
            dense.links.push_back(link);}
        const auto exported=cov::ui::export_nbo_aomo_bundle(compact_view,dense,root/std::filesystem::path(u8"共用证据"));
        require(exported.csv&&exported.json,"shared name export failed");
        const auto pictures=cov::ui::export_nbo_aomo_bundle(compact_view,dense,root/"dense-pictures",cov::DiagramExportContent::Images);
        require(pictures.svg&&pictures.png&&!pictures.csv&&!pictures.json&&
            !std::filesystem::exists(pictures.csv_path)&&!std::filesystem::exists(pictures.json_path)&&
            read(pictures.svg_path)==read(exported.svg_path)&&read(pictures.png_path)==read(exported.png_path),
            "dense NBO picture export must preserve figures and skip all coefficient payloads");
        std::filesystem::create_directory(pictures.json_path);
        {std::ofstream out(pictures.csv_path);out<<"untouched analysis";}
        const auto repeated=cov::ui::export_nbo_aomo_bundle(compact_view,dense,root/"dense-pictures",cov::DiagramExportContent::Images);
        require(repeated.svg&&repeated.png&&std::filesystem::is_directory(pictures.json_path)&&
            read(pictures.csv_path)=="untouched analysis","NBO picture export must not rewrite older analysis");
        const auto data_only=cov::ui::export_nbo_aomo_bundle(compact_view,dense,root/"dense-data",cov::DiagramExportContent::AnalysisData);
        require(data_only.csv&&data_only.json&&!data_only.png&&!data_only.svg&&
            !std::filesystem::exists(data_only.svg_path)&&!std::filesystem::exists(data_only.png_path)&&
            read(data_only.json_path).find("\"object_table_file\":\"dense-data.aomo.csv\"")!=std::string::npos &&
            read(data_only.json_path).find("\"object_table_scope\":\"current_display_objects\"")!=std::string::npos &&
            read(data_only.json_path).find("\"full_numerical_scope\"")==std::string::npos,
            "NBO analysis must be independently exportable with complete source data");
        const auto csv=read(exported.csv_path),json=read(exported.json_path);
        require(csv.size()<150000&&csv.find("source_coefficients")==std::string::npos,
            "CSV must not repeat projected coefficient arrays on each link");
        require(json.find("cov_aomo_unified_view_v5")!=std::string::npos &&
            json.find("\"name_evidence_ref\"")!=std::string::npos,
            "current names must remain once in the shared scoped JSON table");
        require(count(csv,"\n")==compact_view.nodes.size()+1 &&
            csv.find("producer,quoted")==std::string::npos,
            "current CSV must have exactly one row per displayed object, not all raw source links");
        require(csv.find("name_evidence_ref_json,object_evidence_ref_json")!=std::string::npos &&
            csv.find("共用证据.aomo.json")!=std::string::npos &&
            json.find("current display; keys are original global indices")!=std::string::npos,
            "current objects must retain source-index references and the UTF-8 companion filename");

    }
    aomo_state.collapsed_atoms.insert(0);
    const auto folded_atom=draw_aomo(cov::NboOrbitalKind::NAO,snapshot);
    bool found_header=false;
    for(const auto& node:folded_atom->nodes)
        if(node.group_header && node.id=="atom:0" && node.label=="+ Cr1 (2)")found_header=true;
    require(found_header,"folded heavy-atom header must retain the canonical element symbol and real count");
    const auto folded_export=cov::ui::export_nbo_aomo_bundle(*folded_atom,integration,root/"folded-heavy-atom");
    require(folded_export.json && folded_export.svg &&
        read(folded_export.json_path).find("+ Cr1 (2)")!=std::string::npos &&
        read(folded_export.svg_path).find("+ Cr1 (2)")!=std::string::npos,
        "folded header must remain visible in JSON and SVG");
    aomo_state.collapsed_atoms.clear();

    // Exercise the public live-UI path with actual MODiagram snapshots whose
    // unmatched spin slots use the producer's sentinel. Valid counterparts
    // must survive, a one-member row must not turn into a fictitious set, and
    // the exported attention list must equal the real visible canonical set.
    const auto no_counterpart=std::numeric_limits<std::size_t>::max();
    auto sentinel_data=snapshot.data;
    wf.orbitals[3].occupation=1.0;
    // This fixture declares source-energy mode. Its selected energy metadata
    // must describe the wavefunction actually drawn, rather than the unrelated
    // early snapshot energies. Spin metadata deliberately remains stale to
    // exercise the separate immutable canonical-spin identity contract.
    for(std::size_t i=0;i<sentinel_data.metadata.size();++i) {
        sentinel_data.metadata[i].energy_hartree=wf.orbitals[i].energy_hartree;
        sentinel_data.metadata[i].occupation=wf.orbitals[i].occupation;
    }
    auto single=level;single.metadata=sentinel_data.metadata[0];
    single.member_indices={0};single.member_spin_counterparts={no_counterpart};
    single.layout_energy_hartree=wf.orbitals[0].energy_hartree;
    auto pair=level;pair.metadata=sentinel_data.metadata[1];
    pair.member_indices={1,2};pair.member_spin_counterparts={no_counterpart,no_counterpart};
    pair.layout_energy_hartree=(wf.orbitals[1].energy_hartree+wf.orbitals[2].energy_hartree)*0.5;
    auto matched=level;matched.metadata=sentinel_data.metadata[3];
    matched.member_indices={3};matched.member_spin_counterparts={4};
    matched.layout_energy_hartree=wf.orbitals[3].energy_hartree;
    auto invalid_row=level;invalid_row.member_indices={wf.orbitals.size()};
    invalid_row.member_spin_counterparts={no_counterpart};
    sentinel_data.levels={single,pair,matched,invalid_row};
    sentinel_data.selection.included_indices={0,1,2,3,4};
    const auto sentinel_snapshot=cov::make_mo_diagram_view_snapshot(
        sentinel_data,snapshot.options,3,"unmatched-spin-display");
    const auto sentinel_view=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    const std::vector<std::size_t> expected_attention{0,1,2,3,4};
    require(sentinel_view->central_mo_indices==expected_attention,
        "unmatched or invalid indices must not enter attention, while every real counterpart remains");
    for(const auto& tick:sentinel_view->energy_ticks) {
        const auto label=cov::format_energy(tick.energy_hartree,sentinel_view->energy_unit,3);
        const float tick_right=8.0f+ImGui::GetFont()->CalcTextSizeA(
            sentinel_view->label_font_size,10000.0f,0.0f,label.c_str()).x;
        require(sentinel_view->lane_x[0]>=tick_right+16.0f,
            "native energy tick text must leave a measured gutter before orbital bars and masks");
    }
    std::set<std::size_t> shown_canonical;
    std::size_t group_count=0;
    for(const auto& node:sentinel_view->nodes) {
        if(node.canonical_index)shown_canonical.insert(*node.canonical_index);
        if(node.lane!=cov::ui::NboAomoLane::Centre || !node.group_header)continue;
        ++group_count;
    }
    require(group_count==0 && shown_canonical==std::set<std::size_t>({0,1,2,3,4}),
        "all five real canonical members must remain visible without Set buttons");
    const auto central_geometry=[&](std::size_t index)->const cov::ui::NboAomoNode& {
        for(const auto& node:sentinel_view->nodes)
            if(node.canonical_index==index)return node;
        std::cerr<<"missing central member "<<index<<'\n';std::exit(1);
    };
    const auto& centre_pair_a=central_geometry(1);
    const auto& centre_pair_b=central_geometry(2);
    const float centre_mid=sentinel_view->canvas_width*0.5f;
    require(std::abs((std::min(centre_pair_a.x,centre_pair_b.x)+
            std::max(centre_pair_a.x+centre_pair_a.width,
                     centre_pair_b.x+centre_pair_b.width))*0.5f-centre_mid)<0.08f,
        "a degenerate canonical row bar envelope must centre on MO lane");
    auto isolated_data=sentinel_snapshot.data;
    isolated_data.levels={single};isolated_data.selection.included_indices={0};
    const auto isolated_graph=cov::make_mo_diagram_view_snapshot(
        isolated_data,sentinel_snapshot.options,0,"isolated-central-fixture");
    const auto isolated_view=draw_aomo(cov::NboOrbitalKind::NAO,isolated_graph);
    const auto isolated=std::find_if(isolated_view->nodes.begin(),isolated_view->nodes.end(),
        [](const auto& node){return node.canonical_index==0;});
    require(isolated!=isolated_view->nodes.end() &&
        std::abs(isolated->x+isolated->width*0.5f-
                 isolated_view->canvas_width*0.5f)<0.08f,
        "an isolated canonical bar must centre on the full canvas");
    // The real wavefunction, not the older display metadata or NBO links,
    // owns canonical spin: index 4 is beta, while indices 0-3 are alpha.
    // Indices 0 and 4 also share an energy and must remain legible side by side.
    for(const auto& node:sentinel_view->nodes) {
        if(!node.canonical_index)continue;
        const auto index=*node.canonical_index;
        require(wf.basis_count==0 && wf.orbitals[index].coefficients.empty() &&
            wf.orbitals[index].symmetry.empty() &&
            wf.orbitals[index].source_orbital_index==std::numeric_limits<std::size_t>::max(),
            "the unknown-symmetry fixture must retain absent physical and source-index evidence");
        require(node.individual_label==std::string("MO [list] ")+std::to_string(index+1)+
            (index==4?" [beta]":" [alpha]") &&
            node.id=="canonical_mo:"+std::to_string(index) &&
            node.detail.find("original MO "+std::to_string(index+1))!=std::string::npos &&
            !node.symmetry_name_verified && node.symmetry_irrep.empty() && node.symmetry_ordinal==0,
            "source-number fallback must preserve unknown physical symmetry and each actual spin");
        const float right=std::max(node.label_x+node.label_width,
            node.occupation_x+node.occupation_width);
        for(const auto& other:sentinel_view->nodes) {
            if(!other.canonical_index || *other.canonical_index<=index)continue;
            const bool vertical_overlap=node.label_y<other.label_y+other.label_height &&
                other.label_y<node.label_y+node.label_height;
            const float other_right=std::max(other.label_x+other.label_width,
                other.occupation_x+other.occupation_width);
            require(!vertical_overlap || right<=other.label_x || other_right<=node.label_x,
                "expanded alpha/beta canonical labels and occupations must not overlap");
        }
    }
    const auto central_node=[&](std::size_t index)->const cov::ui::NboAomoNode& {
        for(const auto& node:sentinel_view->nodes)if(node.canonical_index==index)return node;
        std::abort();
    };
    const auto& mean_a=central_node(1);const auto& mean_b=central_node(2);
    require(mean_a.energy_hartree==wf.orbitals[1].energy_hartree &&
        mean_b.energy_hartree==wf.orbitals[2].energy_hartree &&
        std::abs(*mean_a.display_energy_hartree-
            (wf.orbitals[1].energy_hartree+wf.orbitals[2].energy_hartree)*0.5)<1e-12 &&
        mean_a.display_energy_hartree==mean_b.display_energy_hartree && mean_a.y==mean_b.y,
        "same-spin degenerate members must share a display mean without modifying either raw eigenvalue");
    require(central_node(3).spatial_pair_id==central_node(4).spatial_pair_id &&
        !central_node(3).spatial_pair_id.empty() &&
        central_node(3).spatial_pair_label.find("↑↓")!=std::string::npos &&
        central_node(3).x==central_node(4).x,
        "matched alpha/beta canonical members must form one readable spatial pair with separate bars");
    require(central_node(3).display_energy_hartree==central_node(3).energy_hartree &&
        central_node(4).display_energy_hartree==central_node(4).energy_hartree &&
        central_node(3).display_group_id!=central_node(4).display_group_id,
        "a matched opposite-spin spatial counterpart is not an energy-degenerate partner");
    // A second display-contract fixture explicitly selects common-operator
    // expectations. These synthetic values are not an electronic-state or
    // operator-qualification reference; those have separate scientific tests.
    // Here their difference from the source energies proves that AOMO uses
    // the selected snapshot values throughout layout, semantics and exports.
    auto common_data=sentinel_snapshot.data;
    common_data.using_ro_common_energy=true;
    common_data.ro_common_energy.available=true;
    common_data.ro_common_energy.status="display_contract_fixture";
    common_data.ro_common_energy.restricted_open_shell.verified=true;
    common_data.ro_common_energy.restricted_open_shell.method="ROHF";
    const double selected_energies[]={.4,.15,.25,.5,.7};
    for(std::size_t i=0;i<5;++i)common_data.metadata[i].energy_hartree=selected_energies[i];
    const auto common_graph=cov::make_mo_diagram_view_snapshot(
        common_data,sentinel_snapshot.options,1,"common-operator-display-fixture");
    const auto common_view=draw_aomo(cov::NboOrbitalKind::NAO,common_graph);
    const auto common_node=[&](std::size_t index)->const cov::ui::NboAomoNode& {
        for(const auto& node:common_view->nodes)if(node.canonical_index==index)return node;
        std::abort();
    };
    const auto& common_a=common_node(1);const auto& common_b=common_node(2);
    const double selected_mean=(selected_energies[1]+selected_energies[2])*0.5;
    require(common_view->using_ro_common_energy&&
        common_view->display_energy_definition==common_data.ro_common_energy.operator_semantics&&
        common_a.energy_hartree==wf.orbitals[1].energy_hartree&&
        common_b.energy_hartree==wf.orbitals[2].energy_hartree&&
        common_a.occupation==wf.orbitals[1].occupation&&
        std::abs(*common_a.display_energy_hartree-selected_mean)<1e-12&&
        common_a.display_energy_hartree==common_b.display_energy_hartree&&common_a.y==common_b.y&&
        common_a.energy_semantics=="source RO effective energy"&&
        common_a.display_energy_semantics.find("spin-average operator expectations")!=std::string::npos&&
        common_a.display_energy_semantics.find("not canonical eigenvalues")!=std::string::npos,
        "common display mode must use selected expectations while retaining actual source values and occupations");
    const double selected_coordinate=cov::energy_display_coordinate(selected_mean,common_view->energy_transform);
    const double selected_y=common_view->numeric_top+(1-(selected_coordinate-common_view->axis_coordinate_min)/
        (common_view->axis_coordinate_max-common_view->axis_coordinate_min))*common_view->numeric_span;
    require(std::abs(common_a.y-selected_y)<1e-4&&common_view->energy_transform.knots.front().energy_hartree>=.15,
        "common-mode node positions and axis knots must use the selected expectation definition");
    const auto common_export=cov::ui::export_nbo_aomo_bundle(*common_view,integration,
        root/"common-energy-contract",cov::DiagramExportContent::AnalysisData);
    const auto common_json=read(common_export.json_path);
    require(common_export.json&&common_export.csv&&
        common_json.find("\"using_ro_common_energy\":true")!=std::string::npos&&
        common_json.find("\"display_energy_definition\":\""+common_view->display_energy_definition+"\"")!=std::string::npos&&
        common_json.find("source RO effective energy")!=std::string::npos&&
        common_json.find("spin-average operator expectations")!=std::string::npos,
        "common-energy exports must retain source and selected energy definitions");
    const auto restored_view=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    const auto restored_a=std::find_if(restored_view->nodes.begin(),restored_view->nodes.end(),
        [](const auto& node){return node.canonical_index==1;});
    require(!restored_view->using_ro_common_energy&&restored_view->id!=common_view->id&&
        restored_a!=restored_view->nodes.end()&&restored_a->display_energy_hartree==mean_a.display_energy_hartree&&
        common_a.display_energy_hartree==selected_mean&&mean_a.energy_hartree==wf.orbitals[1].energy_hartree,
        "switching energy definitions must rebuild the view while preserving both frozen snapshots and source energies");
    for(const auto& node:sentinel_view->nodes)if(node.canonical_index && node.occupation && *node.occupation>0) {
        const float line=node.y+node.height*0.5f;
        require(node.occupation_on_bar && node.occupation_y<line &&
            node.occupation_y+node.occupation_height>line && node.occupation_x>=node.x &&
            node.occupation_x+node.occupation_width<=node.x+node.width,
            "occupied canonical electron arrows must cross their own orbital bar");
    }
    const auto sentinel_export=cov::ui::export_nbo_aomo_bundle(*sentinel_view,integration,root/"unmatched-spin");
    require(sentinel_export.json && sentinel_export.svg && sentinel_export.png,
        "unmatched-spin display export failed");
    const auto sentinel_json=read(sentinel_export.json_path);
    const auto sentinel_svg=read(sentinel_export.svg_path);
    require(sentinel_json.find(std::to_string(no_counterpart))==std::string::npos &&
        sentinel_json.find("\"central_mo_indices\":[0,1,2,3,4]")!=std::string::npos &&
        sentinel_svg.find("degenerate:")==std::string::npos &&
        sentinel_svg.find("id=\"canonical_mo:4\"")!=std::string::npos,
        "exports must exclude sentinel/false-group identities without losing the valid spin counterpart");
    require(sentinel_json.find("MO [list] 1 [alpha]")!=std::string::npos &&
        sentinel_json.find("MO [list] 5 [beta]")!=std::string::npos &&
        sentinel_json.find("\"symmetry_name_verified\":false")!=std::string::npos &&
        sentinel_json.find("\"symmetry_irrep\":\"\"")!=std::string::npos &&
        sentinel_json.find("\"symmetry_ordinal\":0")!=std::string::npos &&
        sentinel_svg.find("data-individual-label=\"MO [list] 1 [alpha]\"")!=std::string::npos &&
        sentinel_svg.find("data-individual-label=\"MO [list] 5 [beta]\"")!=std::string::npos &&
        sentinel_svg.find("id=\"canonical_mo:4\"")!=std::string::npos,
        "JSON and SVG must preserve unknown symmetry, spin and exact canonical identity");
    require(sentinel_json.find("spatial_pair_id")!=std::string::npos &&
        sentinel_svg.find("MO [list] 4 [alpha] / MO [list] 5 [beta]")!=std::string::npos,
        "exports must preserve the shared pair cue and both individual member identities");
    for(auto& orbital:wf.orbitals)orbital.spin=cov::Spin::Alpha;
    aomo_state.names.reset(); // The fixture mutates an otherwise immutable attachment.
    ++aomo_state.revision;
    const auto closed_shell=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    for(const auto& node:closed_shell->nodes)
        if(node.canonical_index)require(node.label=="MO [list] "+std::to_string(*node.canonical_index+1) &&
            node.id=="canonical_mo:"+std::to_string(*node.canonical_index) &&
            !node.symmetry_name_verified && node.symmetry_irrep.empty() && node.symmetry_ordinal==0,
            "a canonical set without beta orbitals must retain unknown symmetry and true member identity");
    // A fixed local-shell fixture with a deliberately Rydberg-dominated virtual
    // MO exercises the filter boundary, compact stack, and raw/display split.
    auto shell_model=std::make_shared<cov::NboSalcModel>();
    shell_model->available=true;shell_model->dataset_id=integration.id;
    shell_model->fragments.push_back({"atom:0","Cr1",{0},0});
    for(std::size_t i=0;i<5;++i) {
        cov::NboSalcOrbital orbital;
        orbital.id="fixture-shell-"+std::to_string(i);orbital.fragment_id="atom:0";
        orbital.type=i==4?"Ryd(4f)":"Val(3p)";
        orbital.angular=i==0?"px":i==1?"py":"pz";
        orbital.spin=i==3?cov::NboSpin::Beta:cov::NboSpin::Alpha;
        orbital.label="raw member "+std::to_string(i);
        orbital.subspace_id=i==4?"extra":i==3?"beta-p":"alpha-p";
        orbital.atoms={0};orbital.energy_hartree=-0.3-0.1*static_cast<double>(i);
        orbital.occupation=i==0?0.03:i==1?1.2:1.0;
        orbital.terms.push_back({{cov::NboOrbitalKind::NAO,orbital.spin,i},1.0});
        shell_model->orbitals.push_back(orbital);
        shell_model->links.push_back({i,2,i==4?0.99:0.02,i==4?0.9801:0.0004});
    }
    aomo_state.salc_model=shell_model;aomo_state.show_rydberg=false;
    aomo_state.focused_canonical_index=2;
    aomo_state.last_inspected=sentinel_snapshot.data.view->inspected_orbital_index;
    ++aomo_state.revision;
    const auto shell_view=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    require(shell_view->focused_canonical_index==2 && !shell_view->auto_rydberg_expanded &&
        !shell_view->show_rydberg,"virtual inspection must preserve the user's Rydberg filter");
    std::vector<const cov::ui::NboAomoNode*> alpha_p;
    for(const auto& node:shell_view->nodes)if(node.salc_index) {
        require(*node.salc_index!=4,"Rydberg-dominated virtual MO must not reveal a hidden higher shell");
        if(*node.salc_index<3)alpha_p.push_back(&node);
        if(*node.salc_index<3)require(node.individual_label.find("3p")!=std::string::npos&&node.individual_label.find("2p")==std::string::npos,
            "NAO component hover/export must retain its producer principal shell");
        if(*node.salc_index<2)require(node.occupation_label.empty() && !node.occupation_on_bar,
            "fractional side occupations must not be rounded into electron arrows");
        if(*node.salc_index==3)require(node.shell_member_count==1 &&
            node.display_energy_hartree==node.energy_hartree,
            "opposite-spin atomic shells must not be averaged together");
    }
    require(alpha_p.size()==3,"all three real p members must remain independently represented");
    std::size_t shared_labels=0;double offset_sum=0;
    for(const auto* node:alpha_p) {
        shared_labels+=!node->label.empty();offset_sum+=node->display_offset_y;
        require(node->x==alpha_p.front()->x && node->shell_member_count==3 &&
            std::abs(*node->display_energy_hartree+0.4)<1e-12 &&
            node->energy_hartree==shell_model->orbitals[*node->salc_index].energy_hartree &&
            std::abs((node->y-node->display_offset_y)-
                (alpha_p.front()->y-alpha_p.front()->display_offset_y))<1e-4,
            "p members must share x and mean-energy anchor while retaining their individual raw energy");
    }
    require(shared_labels==1 && alpha_p.front()->label=="Cr1 3p [alpha]" &&
        std::abs(offset_sum)<1e-6 && alpha_p[0]->y<alpha_p[1]->y && alpha_p[1]->y<alpha_p[2]->y,
        "atomic shell must show one shared shell label and a symmetric ordered compact stack");
    const auto shell_export=cov::ui::export_nbo_aomo_bundle(*shell_view,integration,root/"shell-display");
    require(shell_export.json && shell_export.csv && shell_export.svg && shell_export.png &&
        read(shell_export.json_path).find("\"display_energy_hartree\"")!=std::string::npos &&
        read(shell_export.csv_path).find("source_display_energy_hartree")!=std::string::npos &&
        read(shell_export.svg_path).find("data-display-group=")!=std::string::npos &&
        read(shell_export.svg_path).find("class=\"electron-arrows\"")!=std::string::npos &&
        read(shell_export.svg_path).find("stroke-dasharray=\"5 5\"")!=std::string::npos,
        "all export formats must carry the mean/stack snapshot and on-bar arrow rendering");
    // Crowded but distinct side expectations must reserve vertical typography
    // on the shared reversible axis. Identical energies keep independent
    // centred columns; neither case changes producer energies or coefficients.
    auto right_atom=wf.atoms.front();right_atom.symbol="F";right_atom.atomic_number=9;
    wf.atoms.push_back(right_atom);
    auto crowded=std::make_shared<cov::NboSalcModel>(*shell_model);
    crowded->fragments.push_back({"right-fixture","right",{1},1});
    const auto append_side=[&](const char* fragment,const char* type,
                               const char* angular,double energy) {
        auto orbital=crowded->orbitals.front();
        orbital.id="crowded-"+std::to_string(crowded->orbitals.size());
        orbital.fragment_id=fragment;orbital.subspace_id=orbital.id;
        if(orbital.fragment_id=="right-fixture")orbital.atoms={1};
        orbital.type=type;orbital.angular=angular;orbital.energy_hartree=energy;
        orbital.spin=cov::NboSpin::Alpha;orbital.occupation=1.0;
        crowded->orbitals.push_back(orbital);
        return crowded->orbitals.size()-1;
    };
    const auto left_near=append_side("atom:0","Val(4s)","4s",-0.3995);
    const auto right_a=append_side("right-fixture","Val(2p)","px",-0.4000);
    const auto right_b=append_side("right-fixture","Val(2p)","py",-0.3998);
    const auto right_c=append_side("right-fixture","Val(2s)","2s",-0.3992);
    const auto right_same=append_side("right-fixture","Val(4s)","4s",-0.3992);
    const auto right_next=append_side("right-fixture","Val(3s)","3s",-0.399199);
    auto adaptive_options=sentinel_snapshot.options;
    adaptive_options.energy_axis_mode=cov::EnergyAxisMode::NonlinearFocus;
    const auto adaptive_graph=cov::make_mo_diagram_view_snapshot(
        sentinel_snapshot.data,adaptive_options,3,"adaptive-side-fixture");
    aomo_state.salc_model=crowded;aomo_state.names.reset();++aomo_state.revision;
    const auto crowded_view=draw_aomo(cov::NboOrbitalKind::NAO,adaptive_graph);
    const auto side_node=[&](std::size_t index)->const cov::ui::NboAomoNode& {
        for(const auto& node:crowded_view->nodes)
            if(node.salc_index==index)return node;
        std::cerr<<"missing crowded side node "<<index<<'\n';std::exit(1);
    };
    const auto& p_a=side_node(right_a);
    const auto& p_b=side_node(right_b);
    const auto& s_c=side_node(right_c);
    const auto& s_same=side_node(right_same);
    const auto& s_next=side_node(right_next);
    const auto& left=side_node(left_near);
    require(p_a.x==p_b.x && p_a.shell_member_count==2 &&
        std::abs((p_a.y-p_a.display_offset_y)-
                 (p_b.y-p_b.display_offset_y))<1e-4,
        "verified same-shell side members retain one aligned symmetric stack");
    require(s_c.x!=s_same.x && s_c.display_energy_hartree==s_same.display_energy_hartree,
        "coincident independent side groups retain separate centred fallback columns");
    const float p_top=std::min(p_a.y,p_b.y);
    const float p_bottom=std::max(p_a.y,p_b.y)+22.0f;
    const float p_c_gap=std::max(p_top-(s_c.y+22.0f),s_c.y-p_bottom);
    if(p_c_gap<11.9f)
        std::cerr<<"p/c physical rows: p_top="<<p_top<<" p_bottom="<<p_bottom
            <<" c_top="<<s_c.y<<" c_bottom="<<s_c.y+22.0f
            <<" gap="<<p_c_gap<<" next_top="<<s_next.y<<'\n';
    require(p_c_gap>=11.9f,
        "crowded distinct side bar envelopes must retain readable separation");
    const float bar=crowded_view->orbital_bar_width;
    const float right_row_centre=(std::min(s_c.x,s_same.x)+
        std::max(s_c.x,s_same.x)+bar)*0.5f;
    const float physical_row_gap=std::abs(s_next.y-s_c.y)-22.0f;
    require(physical_row_gap>=11.9f &&
        std::abs((s_next.x+bar*0.5f)-right_row_centre)<0.08f,
        "distinct side bar rows retain readable group spacing and a centred anchor");
    const float mo_mid=crowded_view->lane_x[1]+crowded_view->lane_width[1]*0.5f;
    require(std::abs(side_node(0).x-left.x)<0.08f &&
        std::abs(right_row_centre-(p_a.x+bar*0.5f))<0.08f &&
        std::abs(((left.x+bar*0.5f)+(p_a.x+bar*0.5f))*0.5f-mo_mid)<0.08f &&
        std::abs(crowded_view->canvas_width*0.5f-mo_mid)<0.08f &&
        crowded_view->lane_width[0]==crowded_view->lane_width[2],
        "physical side bar rows and reserved envelopes must mirror about the MO centre");
    require(left.energy_hartree==crowded->orbitals[left_near].energy_hartree &&
        s_c.energy_hartree==crowded->orbitals[right_c].energy_hartree &&
        p_a.energy_hartree==crowded->orbitals[right_a].energy_hartree,
        "adaptive layout must not change any raw side expectation");
    const auto& knots=crowded_view->energy_transform.knots;
    for(std::size_t i=1;i<knots.size();++i)
        require(knots[i].energy_hartree>knots[i-1].energy_hartree &&
            knots[i].coordinate>knots[i-1].coordinate,
            "adaptive shared transform must be strictly monotone");
    for(const auto& node:crowded_view->nodes)if(node.quantitative_energy && node.energy_hartree) {
        const auto coordinate=cov::energy_display_coordinate(*node.energy_hartree,
            crowded_view->energy_transform);
        require(std::abs(cov::energy_from_display_coordinate(coordinate,
            crowded_view->energy_transform)-*node.energy_hartree)<1e-9,
            "adaptive shared transform must invert every raw node energy");
    }
    require(crowded_view->energy_tick_semantics=="adaptive-nonlinear-amber" &&
        crowded_view->energy_tick_screen_rgb!=
            std::array<std::uint8_t,3>{139,157,178},
        "nonlinear numeric ticks must declare a distinct amber palette");
    const auto crowded_export=cov::ui::export_nbo_aomo_bundle(*crowded_view,integration,
        root/"adaptive-crowded");
    require(crowded_export.json && crowded_export.svg && crowded_export.png &&
        read(crowded_export.json_path).find("adaptive-nonlinear-amber")!=std::string::npos &&
        read(crowded_export.svg_path).find("#e1ad59")!=std::string::npos,
        "frozen nonlinear tick color and semantics must reach exported artifacts");
    // A dense side must grow the scrollable nonlinear axis, even when that
    // exceeds the old global height cap. Distinct energies below 1e-7 Ha are
    // still distinct rows; grouping is governed by verified orbital identity.
    auto dense=std::make_shared<cov::NboSalcModel>();
    dense->available=true;dense->dataset_id=integration.id;
    dense->fragments.push_back({"dense-left","left",{0},0});
    dense->fragments.push_back({"dense-right","right",{1},1});
    for(std::size_t i=0;i<80;++i)for(unsigned side=0;side<2;++side) {
        auto orbital=crowded->orbitals.front();
        orbital.id="dense-"+std::to_string(side)+":"+std::to_string(i);
        orbital.fragment_id=side?"dense-right":"dense-left";
        orbital.subspace_id=orbital.id;orbital.atoms={side};
        orbital.type="Val("+std::to_string(i+1)+"s)";
        orbital.spin=cov::NboSpin::Alpha;orbital.energy_hartree=-0.4+i*2e-8;
        dense->orbitals.push_back(orbital);
    }
    aomo_state.salc_model=dense;++aomo_state.revision;
    const auto dense_view=draw_aomo(cov::NboOrbitalKind::NAO,adaptive_graph);
    std::array<std::vector<const cov::ui::NboAomoNode*>,2> dense_sides;
    for(const auto& node:dense_view->nodes)if(node.salc_index)
        dense_sides[node.lane==cov::ui::NboAomoLane::Right?1:0].push_back(&node);
    for(auto& side:dense_sides) {
        require(side.size()==80,"dense layout dropped a real side orbital");
        std::sort(side.begin(),side.end(),[](auto a,auto b){return a->y<b->y;});
        for(std::size_t i=0;i<side.size();++i) {
            require(std::abs(side[i]->x-side.front()->x)<0.08f,
                "distinct dense side rows must stay in one column");
            require(side[i]->energy_hartree==dense->orbitals[*side[i]->salc_index].energy_hartree,
                "dense layout altered a source expectation");
            if(i)require(side[i]->y-side[i-1]->y>=33.9f,
                "dense side group spacing collapsed to satisfy a height cap");
        }
    }
    require(std::abs((dense_sides[0][0]->x+dense_sides[1][0]->x+
        dense_view->orbital_bar_width)*0.5f-dense_view->canvas_width*0.5f)<0.08f,
        "single side columns must mirror around the central MO column");
    // Verified irrep partners do not lose their shared stack when the global
    // occurrence number is unknown. Distinct occurrences/scopes stay separate.
    auto partner_model=std::make_shared<cov::NboSalcModel>();
    partner_model->available=true;partner_model->dataset_id=integration.id;
    partner_model->fragments={{"partners-left","left",{0,1},0},
        {"partners-right","right",{0,1},1}};
    auto partner_names=std::make_shared<cov::ui::NboAomoNames>();
    for(std::size_t block=0;block<3;++block) {
        const std::string space="verified-copy-"+std::to_string(block);
        {
            cov::NboSalcSubspace subspace;subspace.id=space;
            subspace.symmetry_verified=true;subspace.spin=cov::NboSpin::Alpha;
            subspace.dimension=subspace.irrep_dimension=3;subspace.multiplicity=1;
            subspace.energy_degeneracy_verified=true;
            partner_model->subspaces.push_back(subspace);
        }
        for(std::size_t member=0;member<3;++member) {
            auto orbital=crowded->orbitals.front();
            orbital.id="partner-"+std::to_string(block)+":"+std::to_string(member);
            orbital.fragment_id=block<2?"partners-left":"partners-right";
            orbital.subspace_id=space;orbital.atoms={0,1};orbital.type="Val(2p)";
            orbital.energy_hartree=-0.42+0.04*block+1e-6*member;
            orbital.spin=cov::NboSpin::Alpha;orbital.symmetry_adapted=true;
            orbital.partner_index=member;orbital.partner_dimension=3;
            partner_model->orbitals.push_back(orbital);
            cov::ui::NboAomoName name;name.verified=true;name.irrep="T1u";
            name.label="t₁u [SALC "+std::to_string(partner_names->salc.size()+1)+"] [alpha]";name.ordinal=0;
            name.partner_block_id="block-"+std::to_string(block);
            name.partner_block_size=3;partner_names->salc.push_back(name);
        }
    }
    aomo_state.salc_model=partner_model;aomo_state.names=partner_names;
    aomo_state.names_model=partner_model.get();++aomo_state.revision;
    const auto partner_view=draw_aomo(cov::NboOrbitalKind::NAO,adaptive_graph);
    std::set<std::string> partner_ids,display_groups;
    std::array<std::vector<const cov::ui::NboAomoNode*>,3> partners;
    for(const auto& n:partner_view->nodes)if(n.salc_index) {
        partners[*n.salc_index/3].push_back(&n);partner_ids.insert(n.id);
        display_groups.insert(n.display_group_id);
    }
    require(partner_ids.size()==9&&display_groups.size()==3,
        "same-label partner blocks must preserve all independent identities and scopes");
    for(std::size_t block=0;block<3;++block) {
        require(partners[block].size()==3,"a verified triplet lost a real member");
        std::size_t labels=0;
        for(const auto* n:partners[block]) {
            labels+=!n->label.empty();
            const auto display_ordinal=block<2?block+1:1;
            require(n->shell_member_count==3&&n->symmetry_ordinal==display_ordinal&&
                n->shell_label==std::to_string(display_ordinal)+"t₁u [alpha]"&&n->x==partners[block][0]->x&&
                std::abs(*n->display_energy_hartree-(-0.42+0.04*block+1e-6))<1e-12&&
                n->energy_hartree==partner_model->orbitals[*n->salc_index].energy_hartree&&
                partner_names->salc[*n->salc_index].ordinal==0&&
                partner_view->names->salc[*n->salc_index].complete_set_ordinal==0,
                "visible occurrence counters must preserve proved SALC stacks, unresolved complete-set ordinals and raw values");
            if(!n->label.empty())require(block<2?n->label_x+n->label_width<n->x:
                n->label_x>n->x+n->width,"SALC shared labels must face outward on each side");
        }
        require(labels==1&&partners[block][0]->y<partners[block][1]->y&&
            partners[block][1]->y<partners[block][2]->y,
            "verified partners need one label and separate ordered interaction targets");
    }
    partner_model->orbitals[2].type="Ryd(4f)";++aomo_state.revision;
    const auto partial_partners=draw_aomo(cov::NboOrbitalKind::NAO,adaptive_graph);
    for(const auto& n:partial_partners->nodes)if(n.salc_index&&*n.salc_index<2)
        require(n.shell_member_count==1&&n.display_group_id.rfind("verified-salc:",0)!=0,
            "an incomplete filtered occurrence cannot claim a complete SALC partner stack");
    auto unknown_names=std::make_shared<cov::ui::NboAomoNames>(*partner_names);
    unknown_names->salc[0]=cov::ui::NboAomoName{};
    aomo_state.names=unknown_names;++aomo_state.revision;
    const auto unknown_view=draw_aomo(cov::NboOrbitalKind::NAO,adaptive_graph);
    for(const auto& n:unknown_view->nodes)if(n.salc_index==0)
        require(!n.symmetry_name_verified&&!n.label.empty()&&n.individual_label.find("SALC 1")!=std::string::npos,
            "an unresolved SALC name must retain its real source identity instead of a blank line");
    aomo_state.salc_model=crowded;aomo_state.names.reset();++aomo_state.revision;
    adaptive_options.energy_axis_mode=cov::EnergyAxisMode::Linear;
    const auto linear_graph=cov::make_mo_diagram_view_snapshot(
        sentinel_snapshot.data,adaptive_options,3,"linear-side-fixture");
    ++aomo_state.revision;
    const auto linear_side=draw_aomo(cov::NboOrbitalKind::NAO,linear_graph);
    require(linear_side->energy_tick_semantics=="linear-neutral" &&
        linear_side->energy_tick_screen_rgb==
            std::array<std::uint8_t,3>{139,157,178} &&
        linear_side->nodes.size()==crowded_view->nodes.size(),
        "explicit linear mode must retain neutral ticks and the same real nodes");
    for(std::size_t i=0;i<linear_side->nodes.size();++i)
        require(linear_side->nodes[i].id==crowded_view->nodes[i].id &&
            linear_side->nodes[i].energy_hartree==crowded_view->nodes[i].energy_hartree,
            "axis mode must not alter node identity or raw energy");
    wf.atoms.pop_back();
    aomo_state.salc_model=shell_model;aomo_state.names.reset();++aomo_state.revision;
    aomo_state.show_rydberg=true;++aomo_state.revision;
    const auto explicit_extra=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    bool found_extra=false;
    for(const auto& node:explicit_extra->nodes)if(node.salc_index==4)found_extra=true;
    require(found_extra && explicit_extra->central_mo_indices==shell_view->central_mo_indices,
        "an explicit extra-shell toggle may reveal Rydberg without changing canonical attention");
    // A producer valence frame alone does not reveal arbitrary Ryd shells.
    // Only a shell actually needed by retained MO composition is shown;
    // the producer's Ryd class stays unchanged.
    auto metal_model=std::make_shared<cov::NboSalcModel>(*shell_model);
    for(const auto& [id,type]:std::vector<std::pair<std::size_t,std::string>>{
        {5,"Ryd(4p)"},{6,"Ryd(5p)"}}) {
        auto orbital=metal_model->orbitals.front();
        orbital.id="frame-test-"+std::to_string(id);
        orbital.type=type;orbital.subspace_id=orbital.id;
        orbital.terms={{{cov::NboOrbitalKind::NAO,cov::NboSpin::Alpha,id},1.0}};
        metal_model->orbitals.push_back(orbital);
    }
    for(const auto& type:{"Val(3d)","Val(4s)"}) {
        cov::NboNao nao;nao.id=integration.dataset.naos.size()+1;
        nao.atom=1;nao.spin=cov::NboSpin::Alpha;nao.type=type;
        integration.dataset.naos.push_back(nao);
    }
    aomo_state.salc_model=metal_model;aomo_state.names.reset();
    aomo_state.show_rydberg=false;++aomo_state.revision;
    const auto unsupported_frame=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    for(const auto& node:unsupported_frame->nodes)
        require(node.salc_index!=5 && node.salc_index!=6,
            "a valence frame without current MO composition must not reveal Ryd p shells");
    cov::NboMoDecomposition p_need;
    p_need.canonical_index=2;p_need.available=true;p_need.spin=cov::NboSpin::Alpha;
    cov::NboNaoContribution p_row;
    p_row.atom=1;p_row.nao_id=5;p_row.type="Ryd(4p)";p_row.weight=0.5;
    p_need.rows.push_back(p_row);
    integration.dataset.mo_decompositions.push_back(p_need);
    ++aomo_state.revision;
    const auto frame_view=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    bool shown_4p=false,shown_5p=false;
    for(const auto& node:frame_view->nodes)if(node.salc_index) {
        shown_4p|=*node.salc_index==5;
        shown_5p|=*node.salc_index==6;
    }
    require(shown_4p && !shown_5p && !frame_view->show_rydberg &&
        metal_model->orbitals[5].type=="Ryd(4p)",
        "default metal frame must keep evidenced Ryd 4p without opening other Ryd shells");
    wf.atoms[0].atomic_number=1;
    aomo_state.hide_h_orbitals=true;++aomo_state.revision;
    const auto hidden_h=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    require(hidden_h->hidden_h_count>0 &&
        std::none_of(hidden_h->nodes.begin(),hidden_h->nodes.end(),
            [](const auto& node){return node.salc_index.has_value();}) &&
        hidden_h->central_mo_indices==frame_view->central_mo_indices,
        "H display toggle must hide pure-H side orbitals without changing central identities");
    aomo_state.hide_h_orbitals=false;++aomo_state.revision;
    const auto restored_h=draw_aomo(cov::NboOrbitalKind::NAO,sentinel_snapshot);
    require(std::any_of(restored_h->nodes.begin(),restored_h->nodes.end(),
        [](const auto& node){return node.salc_index.has_value();}),
        "clearing the H toggle must restore the same underlying side orbitals");
    ImGui::DestroyContext();
    // Only this uniquely created test directory is removed.
    std::filesystem::remove_all(root);
    std::cout << "mo_diagram_snapshot_smoke ok\n";
}
