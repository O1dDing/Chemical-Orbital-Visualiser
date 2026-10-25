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
#include <imgui.h>

struct ImGuiWindow;
namespace cov::validation {
struct NavigationTarget {
    ImVec2 lo, hi;
    ImGuiWindow* window = nullptr;
    ImVec2 clip_lo, clip_hi;
};
enum class NavigationKind { Ready, Move, Wheel, Unreachable };
struct NavigationStep {
    NavigationKind kind = NavigationKind::Unreachable;
    ImVec2 mouse{}, wheel{};
    ImGuiWindow* scrolling_window = nullptr;
};
// Observe current ImGui geometry; propose only ordinary pointer/wheel input.
// A pointer move occupies its own frame, so ImGui's input trickling cannot
// leave our wheel queued for the native driver's next ClearEventsQueue().
// Hover evidence may request the entire item on axes where it fits the clip.
// Ordinary clicks retain their partial-visibility behavior.
[[nodiscard]] NavigationStep plan_navigation(const NavigationTarget& target, bool reveal_entire_item = false);
[[nodiscard]] bool navigation_target_visible(const NavigationTarget& target);
}
