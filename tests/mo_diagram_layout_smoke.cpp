/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#include "cov/mo_diagram_layout.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/orbital_ui.hpp"
#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {
void require(bool value,const char* message) {
    if (!value) { std::cerr<<message<<'\n'; std::exit(1); }
}
void inspect_packing(const std::vector<cov::DiagramRowFootprint>& rows,double width) {
    const auto packed=cov::layout_diagram_lanes(rows,width,8,2);
    for (std::size_t a=0;a<rows.size();++a) {
        require(packed.centre_x[a]-rows[a].width/2>=-1e-9 &&
            packed.centre_x[a]+rows[a].width/2<=packed.width+1e-9,"row outside packed canvas");
        for (std::size_t b=a+1;b<rows.size();++b) {
            const bool overlap=rows[a].y-rows[a].above<=rows[b].y+rows[b].below+2 &&
                rows[b].y-rows[b].above<=rows[a].y+rows[a].below+2;
            if (overlap) require(std::abs(packed.centre_x[a]-packed.centre_x[b])>=
                (rows[a].width+rows[b].width)/2+8-1e-9,"distinct row footprints overlap");
        }
    }
}
}

int main() {
    // Display geometry controls, not synthetic evidence of a physical state.
    std::vector<cov::DiagramRowFootprint> rows{{0,14,14,24},{0.02,14,14,117},
        {0.03,14,14,55},{100,14,35,200},{135,14,14,40}};
    inspect_packing(rows,1); inspect_packing(rows,1000);
    const auto base=cov::layout_diagram_lanes(rows,300,8,2);
    auto shifted=rows;
    for (auto& row:shifted) row.y+=1234;
    const auto translated=cov::layout_diagram_lanes(shifted,300,8,2);
    require(base.centre_x==translated.centre_x,"translation changed horizontal packing");
    for (auto& row:shifted) { row.y*=2;row.above*=2;row.below*=2;row.width*=2; }
    const auto scaled=cov::layout_diagram_lanes(shifted,600,16,4);
    for (std::size_t i=0;i<rows.size();++i) require(std::abs(scaled.centre_x[i]-2*base.centre_x[i])<1e-9,
        "DPI scaling changed lane assignment");
    std::vector<cov::DiagramRowFootprint> chain;
    for (int i=0;i<100;++i) chain.push_back({20.0*i,11,11,40});
    inspect_packing(chain,1);
    require(cov::layout_diagram_lanes(chain,1,8,2).width==88,"ended lanes were not reused");

    IMGUI_CHECKVERSION();ImGui::CreateContext();
    auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={1000,1000};io.DeltaTime=1.0f/60;
    io.Fonts->AddFontDefault();io.Fonts->Build();
    cov::Wavefunction wf;wf.atoms.resize(1);wf.atoms[0].atomic_number=1;
    for (const auto energy:{-0.30,-0.2999862,0.20}) {
        cov::MolecularOrbital mo;mo.energy_hartree=energy;mo.occupation=energy<0?2.0f:0.0f;
        wf.orbitals.push_back(mo);
    }
    cov::ui::OrbitalUIState state;state.energy_axis_mode=cov::EnergyAxisMode::Linear;
    state.hide_ligand_centred_intermediates=false;state.filter.mode=cov::OrbitalFilterMode::All;
    state.diagram_neighbourhood=32;
    cov::ui::OrbitalUIActions actions;
    auto frame=[&]() {
        actions={};ImGui::NewFrame();
        ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({370,950},ImGuiCond_Always);
        ImGui::Begin("physical energy lanes",nullptr,ImGuiWindowFlags_NoSavedSettings);
        cov::ui::draw_energy_diagram(wf,0,state,cov::ui::Language::English,1,actions);
        ImGui::End();ImGui::Render();
    };
    auto child=[]() -> ImGuiWindow* {
        for (auto* window:ImGui::GetCurrentContext()->Windows) {
            if (std::string(window->Name).find("##energy_diagram_scroll")!=std::string::npos) return window;
        }
        return nullptr;
    };
    frame();frame();
    require(actions.drawn_diagram && actions.drawn_diagram->data.levels.size()==3,
        "close but distinct energies must not be merged into a fake degeneracy");
    auto* window=child();require(window!=nullptr,"scrollable energy canvas missing");
    // Locate the actual rendered line vertices, independently of the lane
    // planner, then click those coordinates through normal ImGui mouse input.
    auto stroke_centre=[&](ImU32 colour) {
        ImVec2 lo(1e9f,1e9f),hi(-1e9f,-1e9f);
        for (const auto& vertex:window->DrawList->VtxBuffer) if (vertex.col==colour) {
            lo.x=std::min(lo.x,vertex.pos.x);lo.y=std::min(lo.y,vertex.pos.y);
            hi.x=std::max(hi.x,vertex.pos.x);hi.y=std::max(hi.y,vertex.pos.y);
        }
        require(lo.x<=hi.x,"actual level stroke vertices missing");
        return ImVec2((lo.x+hi.x)/2,(lo.y+hi.y)/2);
    };
    const auto first=stroke_centre(IM_COL32(196,210,227,255));
    const auto second=stroke_centre(IM_COL32(79,210,157,255));
    require(std::abs(first.y-second.y)<0.1 && std::abs(first.x-second.x)>40,
        "linear energy must be preserved while nearby members remain individually accessible");
    for (const auto index:{0,1}) {
        const auto p=index==0?first:second;
        io.AddMousePosEvent(p.x,p.y);frame();
        io.AddMouseButtonEvent(0,true);frame();
        require(actions.select_orbital && *actions.select_orbital==static_cast<std::size_t>(index),
            "a nearby row cannot be individually selected at its actual drawn position");
        io.AddMouseButtonEvent(0,false);frame();
    }
    wf.orbitals.clear();
    for (int i=0;i<21;++i) {
        cov::MolecularOrbital mo;mo.energy_hartree=i<20?-0.3+i*0.00002:0.2;mo.occupation=i<20?2.0f:0.0f;
        wf.orbitals.push_back(mo);
    }
    state.diagram_cache={};frame();frame();
    window=child();require(window->ScrollMax.x>0,"dense diagram must retain full-width scrollable lanes");
    const auto bar=ImGui::GetWindowScrollbarRect(window,ImGuiAxis_X).GetCenter();
    io.AddMousePosEvent(bar.x,bar.y);frame();
    io.AddMouseWheelEvent(-3,0);frame();frame();
    require(window->Scroll.x>0,"real horizontal-wheel input cannot reach overflowing lanes");
    ImGui::DestroyContext();
    std::cout<<"mo_diagram_layout_smoke ok\n";
}
