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

#include "cov/validation_navigation.hpp"
#include <imgui_internal.h>
#include <algorithm>
#include <array>
#include <cmath>

namespace cov::validation {
namespace {
ImRect desktop() { return {ImVec2(1,1),ImVec2(ImGui::GetIO().DisplaySize.x-1,ImGui::GetIO().DisplaySize.y-1)}; }
ImRect visible_target(const NavigationTarget& target) {
    ImRect result(target.lo,target.hi);
    result.ClipWith(ImRect(target.clip_lo,target.clip_hi));
    result.ClipWith(desktop());
    return result;
}
bool usable(const ImRect& rect) { return rect.GetWidth()>1 && rect.GetHeight()>1; }
ImRect visible_window(const ImGuiWindow* window) {
    auto result=window->InnerRect;
    result.ClipWith(window->ClipRect);result.ClipWith(desktop());
    return result;
}
ImVec2 wheel_position(ImGuiWindow* window,ImGuiAxis axis) {
    // A scrollbar belongs to its window and cannot be intercepted by a
    // nested child. When the bar is clipped, wheel over visible content.
    ImRect area;
    if (axis==ImGuiAxis_X?window->ScrollbarX:window->ScrollbarY) {
        area=ImGui::GetWindowScrollbarRect(window,axis);
        area.ClipWith(window->OuterRectClipped);area.ClipWith(desktop());
        for (auto* parent=window->ParentWindow;parent;parent=parent->ParentWindow)
            area.ClipWith(parent->ClipRect);
    }
    if (!usable(area)) area=visible_window(window);
    if (!usable(area)) return ImVec2(-1,-1);
    auto point=area.GetCenter();
    const auto& context=*ImGui::GetCurrentContext();
    if (context.WheelingWindow && context.WheelingWindow!=window) {
        const auto ref=context.WheelingWindowRefMousePos;
        const float limit=ImGui::GetIO().MouseDragThreshold+2;
        if (std::hypot(point.x-ref.x,point.y-ref.y)<=limit) {
            const std::array<ImVec2,4> corners{{
                {area.Min.x+1,area.Min.y+1},{area.Max.x-1,area.Min.y+1},
                {area.Min.x+1,area.Max.y-1},{area.Max.x-1,area.Max.y-1}}};
            for (const auto candidate:corners)
                if (std::hypot(candidate.x-ref.x,candidate.y-ref.y)>
                    std::hypot(point.x-ref.x,point.y-ref.y)) point=candidate;
        }
    }
    // AddMousePosEvent quantizes to integer pixels. Compare and emit that
    // same position, or a half-pixel centre would request Move forever.
    return {std::max(std::ceil(area.Min.x),std::floor(point.x)),
            std::max(std::ceil(area.Min.y),std::floor(point.y))};
}
}
bool navigation_target_visible(const NavigationTarget& target) { return usable(visible_target(target)); }
NavigationStep plan_navigation(const NavigationTarget& target, bool reveal_entire_item) {
    const auto visible=visible_target(target);
    auto clip=ImRect(target.clip_lo,target.clip_hi);clip.ClipWith(desktop());
    bool ready=usable(visible);
    if (reveal_entire_item) {
        for (const auto axis:{ImGuiAxis_Y,ImGuiAxis_X}) {
            const float span=target.hi[axis]-target.lo[axis];
            if (span<=clip.Max[axis]-clip.Min[axis] &&
                (target.lo[axis]<clip.Min[axis] || target.hi[axis]>clip.Max[axis])) ready=false;
        }
    }
    if (ready) return {NavigationKind::Ready,visible.GetCenter()};
    const ImVec2 point((target.lo.x+target.hi.x)*0.5f,(target.lo.y+target.hi.y)*0.5f);
    // Scroll only a viewport that actually excludes this target on this axis.
    // A horizontal miss must not unconditionally realign vertical ancestors.
    for (auto* window=target.window;window;window=window->ParentWindow) {
        if (window->Flags&(ImGuiWindowFlags_NoScrollWithMouse|ImGuiWindowFlags_NoMouseInputs)) continue;
        for (const auto axis:{ImGuiAxis_Y,ImGuiAxis_X}) {
            if (window->ScrollMax[axis]<=0) continue;
            const float lo=window->InnerRect.Min[axis]+4,hi=window->InnerRect.Max[axis]-4;
            if (hi<=lo) continue;
            float delta=point[axis]-std::clamp(point[axis],lo,hi);
            if (reveal_entire_item && target.hi[axis]-target.lo[axis]<=hi-lo) {
                if (target.lo[axis]<lo) delta=target.lo[axis]-lo;
                else if (target.hi[axis]>hi) delta=target.hi[axis]-hi;
                else delta=0;
            }
            if (std::abs(delta)<0.5f) continue;
            const float desired=std::clamp(window->Scroll[axis]+delta,0.0f,window->ScrollMax[axis]);
            if (std::abs(desired-window->Scroll[axis])<0.5f) continue;
            const auto mouse=wheel_position(window,axis);
            if (mouse.x<1 || mouse.y<1) continue;
            // Match the installed ImGui wheel pixel scale; only the emitted
            // wheel event changes its scroll state. No SetScrollX/Y call here.
            const float span=window->InnerRect.Max[axis]-window->InnerRect.Min[axis];
            const float pixels=std::floor(std::min((axis==ImGuiAxis_X?2.0f:5.0f)*window->CalcFontSize(),span*0.67f));
            if (pixels<1) continue;
            // ImGui truncates scroll coordinates to whole pixels. A half-pixel
            // correction otherwise repeats forever without moving, preventing
            // navigation from reaching the next scrollable ancestor.
            const float wheel_delta=std::copysign(std::ceil(std::abs(delta)),delta);
            ImVec2 wheel{};wheel[axis]=std::clamp(-wheel_delta/pixels,-3.0f,3.0f);
            const auto actual_mouse=ImGui::GetIO().MousePos;
            const bool needs_move=std::abs(actual_mouse.x-mouse.x)>0.01f || std::abs(actual_mouse.y-mouse.y)>0.01f;
            return {needs_move?NavigationKind::Move:NavigationKind::Wheel,mouse,wheel,window};
        }
    }
    return {};
}
}
