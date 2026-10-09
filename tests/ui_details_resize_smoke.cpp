#include "cov/orbital_ui.hpp"
#include "cov/orbital_ui_text.hpp"
#include "cov/ui_window_layers.hpp"

#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
void require(bool value, const char* message) {
    if (!value) { std::cerr << message << '\n'; std::exit(1); }
}
ImGuiWindow* details() {
    auto* window = ImGui::FindWindowByName("###cov.orbital.details");
    require(window && window->Active && !window->Hidden, "actual details window missing");
    return window;
}
void inspect_bounds(const char* scenario) {
    const auto* window = details();
    const auto* viewport = ImGui::GetMainViewport();
    const auto pos = window->Pos, size = window->Size;
    const auto start = viewport->WorkPos, extent = viewport->WorkSize;
    std::cout << scenario << ": position=" << pos.x << ',' << pos.y
              << " size=" << size.x << ',' << size.y
              << " viewport=" << extent.x << ',' << extent.y << '\n';
    require(pos.x >= start.x && pos.y >= start.y &&
            pos.x + size.x <= start.x + extent.x + 1 &&
            pos.y + size.y <= start.y + extent.y + 1,
            "details title/body/close area escapes actual viewport");
    require(window->TitleBarRect().GetWidth() > 80,
            "details title and close affordance must remain usable");
}
}

int main() {
    // A display fixture. It makes no claim about a physical electronic state.
    IMGUI_CHECKVERSION(); ImGui::CreateContext();
    auto& io = ImGui::GetIO(); io.IniFilename = nullptr;
    io.DisplaySize = {3413,1392}; io.DeltaTime = 1.0f/60;
    io.Fonts->AddFontDefault(); io.Fonts->Build();
    cov::Wavefunction wf; wf.atoms.resize(1); wf.atoms[0].atomic_number = 1;
    for (int i = 0; i < 8; ++i) {
        cov::MolecularOrbital mo;
        mo.energy_hartree = -0.5 + 0.11*(i%4);
        mo.occupation = i%4 < 2 ? 1.0f : 0.0f;
        mo.spin = i >= 4 ? cov::Spin::Beta : cov::Spin::Alpha;
        wf.orbitals.push_back(mo);
    }
    cov::ui::OrbitalUIState state;
    state.hide_ligand_centred_intermediates = false;
    state.filter.mode = cov::OrbitalFilterMode::All;
    state.show_diagram_details = true;
    cov::ui::OrbitalUIActions actions;
    auto language = cov::ui::Language::English;
    auto frame = [&]() {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({0,0},ImGuiCond_Always);
        ImGui::SetNextWindowSize({500,std::max(200.0f,io.DisplaySize.y)},ImGuiCond_Always);
        ImGui::Begin("details test controls",nullptr,ImGuiWindowFlags_NoSavedSettings |
                     cov::ui::background_panel_flags);
        cov::ui::draw_energy_diagram(wf,7,state,language,1,actions);
        ImGui::End(); ImGui::Render();
        require(actions.drawn_diagram &&
                actions.drawn_diagram->data.view->inspected_orbital_index == 7,
                "window interaction changed the inspected source MO");
    };
    frame(); frame(); frame(); inspect_bounds("maximized-open");
    const auto snapshot = actions.drawn_diagram;
    auto drag = [&](ImVec2 from, ImVec2 to) {
        io.AddMousePosEvent(from.x,from.y); frame();
        io.AddMouseButtonEvent(0,true); frame();
        io.AddMousePosEvent(to.x,to.y); frame();
        io.AddMouseButtonEvent(0,false); frame(); frame();
    };
    auto* window = details();
    const auto before_drag = window->Pos;
    const ImVec2 title(before_drag.x+100,before_drag.y+10);
    drag(title,{title.x+400,title.y+50});
    require(details()->Pos.x > before_drag.x+300, "normal title dragging stopped working");
    const auto dragged = details()->Pos; frame(); frame();
    require(details()->Pos.x == dragged.x && details()->Pos.y == dragged.y,
            "valid user position must not be recentered each frame");
    const auto before_size = details()->Size;
    const ImVec2 grip(dragged.x+before_size.x-3,dragged.y+before_size.y-3);
    drag(grip,{grip.x-100,grip.y-100});
    require(details()->Size.x < before_size.x-50 && details()->Size.y < before_size.y-50,
            "normal user resizing stopped working");
    inspect_bounds("moved-and-resized");

    // Reproduce the persistent loss of the floating details behind the left
    // panel, including opening details again after that panel receives focus.
    ImGui::SetWindowPos("###cov.orbital.details",{100,100}); frame(); frame();
    auto assert_above_panel = [&]() {
        const auto* panel = ImGui::FindWindowByName("details test controls");
        int panel_order=-1, details_order=-1;
        for(int i=0;i<GImGui->Windows.Size;++i) {
            if(GImGui->Windows[i]==panel) panel_order=i;
            if(GImGui->Windows[i]==details()) details_order=i;
        }
        require(details_order>panel_order,"base panel permanently covered floating details");
    };
    io.AddMousePosEvent(25,40); frame();
    io.AddMouseButtonEvent(0,true); frame();
    io.AddMouseButtonEvent(0,false); frame(); frame(); assert_above_panel();
    state.focus_diagram_details=true; frame(); frame(); assert_above_panel();
    require(!state.focus_diagram_details,"details must not steal focus every frame");

    for (const ImVec2 size : {ImVec2(1003,658),ImVec2(960,1600),ImVec2(640,360),
                              ImVec2(320,240),ImVec2(2100,1250),ImVec2(3413,1392)}) {
        io.DisplaySize = size; frame(); inspect_bounds("open-across-resize");
        frame(); inspect_bounds("settled-after-resize");
        require(actions.drawn_diagram == snapshot, "resize changed scientific view snapshot");
    }
    auto* original_window = details();
    for (const auto next : {cov::ui::Language::ChineseSimplified,
                            cov::ui::Language::Japanese,cov::ui::Language::French}) {
        language = next; frame(); frame(); inspect_bounds("language-change");
        require(details() == original_window, "localisation created a different details window");
    }
    ImGui::SetWindowCollapsed("###cov.orbital.details",true); frame();
    io.DisplaySize = {640,360}; frame(); frame(); inspect_bounds("collapsed-resize");
    ImGui::SetWindowCollapsed("###cov.orbital.details",false); frame(); frame();
    inspect_bounds("expanded-after-resize");

    // Click the actual first content button using the rendered window's cursor origin.
    const auto start = details()->DC.CursorStartPos;
    io.AddMousePosEvent(start.x+20,start.y+ImGui::GetFrameHeight()/2); frame();
    io.AddMouseButtonEvent(0,true); frame();
    io.AddMouseButtonEvent(0,false); frame();
    require(!state.show_diagram_details, "visible close button did not close the details");
    state.show_diagram_details = true; frame(); frame(); inspect_bounds("reopened-small");

    // Old saved expansion state must not resurrect the removed long chemistry
    // panel. Adding pair evidence must not add sidebar text or alter its height.
    state.show_diagram_details = false;
    io.DisplaySize = {2100,1250};
    auto& chemistry = wf.orbitals[7].chemistry;
    chemistry.available = true;
    cov::OrbitalPairInteraction pair;
    pair.atom_a_label = "Cl2"; pair.atom_b_label = "Cl5";
    pair.channel = {.sigma=.64,.pi=.36,.delta=0,.phi=0,.undetermined=0,
                    .dominant=cov::OrbitalAngularFamily::Sigma,
                    .status=cov::ChemistryStatus::Percentages};
    pair.bonding.dominant = cov::OrbitalBondingRole::Nonbonding;
    pair.bonding.status = cov::ChemistryStatus::Determined;
    pair.occupied_overlap_contribution = .0003; pair.total_mayer_index = .1247;
    auto card_height = [&](float width, int pairs) {
        chemistry.interactions.assign(pairs,pair);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({0,0},ImGuiCond_Always);
        ImGui::SetNextWindowSize({width,1200},ImGuiCond_Always);
        ImGui::Begin("chemistry wrap test",nullptr,ImGuiWindowFlags_NoSavedSettings);
        ImGui::PushTextWrapPos(0);
        // The ordinary sidebar now starts collapsed; this regression measures
        // wrapping after the user opens the same real chemistry section.
        const std::string chemistry_header=std::string(cov::ui::orbital_tr(
            cov::ui::OrbitalText::SelectedMOChemistry,language))+"##diagram.chemistry";
        ImGui::GetStateStorage()->SetInt(ImGui::GetID(chemistry_header.c_str()),1);
        cov::ui::draw_energy_diagram(wf,7,state,language,1,actions);
        const float height = ImGui::GetCursorPosY();
        ImGui::PopTextWrapPos(); ImGui::End(); ImGui::Render();
        return height;
    };
    for (const auto next : {cov::ui::Language::English,cov::ui::Language::ChineseSimplified,
                            cov::ui::Language::Japanese,cov::ui::Language::French}) {
        language = next;
        for (const float width : {500.0f,350.0f}) {
            card_height(width,1); card_height(width,1);
            const float one = card_height(width,1);
            const float two = card_height(width,2);
            const float extra_lines = (two-one)/ImGui::GetTextLineHeightWithSpacing();
            std::cout << "pair wrapping: width=" << width << " language=" << int(language)
                      << " extra lines=" << extra_lines << '\n';
            require(std::abs(extra_lines)<0.01f,
                    "removed long chemistry panel is still being rendered");
        }
    }
    ImGui::DestroyContext();
    std::cout << "ui_details_resize_smoke ok\n";
}
