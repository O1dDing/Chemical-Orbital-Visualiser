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

#pragma once
#include "cov/scientific_text.hpp"
#include <imgui.h>
#include <limits>

namespace cov::ui {
inline ScientificTextLayout scientific_text_layout(const std::string& text,float size,
        bool point_group=false,bool atom_indices=false) {
    return layout_scientific_text(text,size,[&](const std::string& run,float scale) {
        if(!ImGui::GetCurrentContext() || !ImGui::GetFont()) {
            std::size_t glyphs=0;for(unsigned char c:run)if((c&0xc0)!=0x80)++glyphs;
            return static_cast<float>(glyphs)*size*.6f*scale;
        }
        return ImGui::GetFont()->CalcTextSizeA(size*scale,
            std::numeric_limits<float>::max(),0,run.c_str()).x;
    },point_group,atom_indices);
}
inline void draw_scientific_text(ImDrawList* draw,float size,ImVec2 position,
                                ImU32 colour,const std::string& text,
                                bool point_group=false,bool atom_indices=false) {
    const auto layout=scientific_text_layout(text,size,point_group,atom_indices);
    for(const auto& run:layout.runs)
        draw->AddText(ImGui::GetFont(),size*run.scale,
            ImVec2(position.x+run.x,position.y+run.y),colour,run.text.c_str());
}
inline void scientific_text(const std::string& text,
                            ImU32 colour=ImGui::GetColorU32(ImGuiCol_Text),
                            bool wrap=false,bool point_group=false,bool atom_indices=false) {
    const auto origin=ImGui::GetCursorScreenPos();
    const float size=ImGui::GetFontSize(),line=ImGui::GetTextLineHeight();
    const float limit=wrap?std::max(1.0f,ImGui::GetContentRegionAvail().x):
        std::numeric_limits<float>::max();
    float x=0,y=0,width=0;
    // Wrap whole words; never separate an orbital's scripts from its base.
    for(std::size_t at=0;at<text.size();) {
        const auto end=text.find_first_of(" \n",at);
        const auto next=end==std::string::npos?text.size():end;
        const auto word=text.substr(at,next-at);
        const auto layout=scientific_text_layout(word,size,point_group,atom_indices);
        if(x>0 && x+layout.width>limit){x=0;y+=line;}
        draw_scientific_text(ImGui::GetWindowDrawList(),size,
            ImVec2(origin.x+x,origin.y+y),colour,word,point_group,atom_indices);
        x+=layout.width;width=std::max(width,x);
        if(next==text.size())break;
        if(text[next]=='\n'){x=0;y+=line;}
        else x+=ImGui::CalcTextSize(" ").x;
        at=next+1;
    }
    ImGui::Dummy(ImVec2(width,y+line*1.08f));
}
inline bool scientific_selectable(const std::string& label,bool selected=false,
                                  ImGuiSelectableFlags flags=0,bool atom_indices=false) {
    const auto origin=ImGui::GetCursorScreenPos();
    const auto colour=ImGui::GetColorU32(ImGuiCol_Text);
    // Retain ImGui's original identity (including ##/### suffixes), hit testing
    // and keyboard navigation. Only replace its flat visible text.
    ImGui::PushStyleColor(ImGuiCol_Text,IM_COL32(0,0,0,0));
    const bool clicked=ImGui::Selectable(label.c_str(),selected,flags,
        ImVec2(0,ImGui::GetTextLineHeight()*1.08f));
    ImGui::PopStyleColor();
    draw_scientific_text(ImGui::GetWindowDrawList(),ImGui::GetFontSize(),
        origin,colour,label.substr(0,label.find("##")),false,atom_indices);
    return clicked;
}
inline bool scientific_combo(const char* id,const std::string& preview) {
    auto* draw=ImGui::GetWindowDrawList();
    const auto origin=ImGui::GetCursorScreenPos();
    const auto padding=ImGui::GetStyle().FramePadding;
    const float size=ImGui::GetFontSize(),width=ImGui::CalcItemWidth();
    const float height=ImGui::GetFrameHeight();
    const auto colour=ImGui::GetColorU32(ImGuiCol_Text);
    const bool open=ImGui::BeginCombo(id,"");
    draw->PushClipRect(origin,ImVec2(origin.x+width-height,origin.y+height),true);
    draw_scientific_text(draw,size,ImVec2(origin.x+padding.x,origin.y+padding.y),colour,preview);
    draw->PopClipRect();
    return open;
}
} // namespace cov::ui
