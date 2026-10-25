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

#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <utility>

namespace cov {

MODiagramViewSnapshot make_mo_diagram_view_snapshot(
    const MODiagramData& source, const MODiagramOptions& options,
    std::optional<std::size_t> inspected, std::string origin) {
    MODiagramData data=source;
    if (inspected && std::none_of(data.metadata.begin(),data.metadata.end(),
        [&](const auto& item) { return item.orbital_index==*inspected; })) {
        inspected.reset();
    }
    static std::atomic<std::uint64_t> serial{0};
    const auto tick=std::chrono::steady_clock::now().time_since_epoch().count();
    data.view=MODiagramViewContext{
        "view-"+std::to_string(tick)+"-"+std::to_string(++serial),
        std::move(origin),options.selected_index,inspected};
    for (auto& item:data.metadata) {
        item.selected=inspected && item.orbital_index==*inspected;
    }
    for (auto& level:data.levels) {
        level.metadata.selected=inspected &&
            mo_diagram_level_covers_orbital(level,*inspected);
    }
    return {std::move(data),options};
}

std::vector<MODiagramMemberView> mo_diagram_member_views(
    const MODiagramData& data, const MODiagramLevel& level) {
    const std::size_t count=std::max<std::size_t>(1,level.member_indices.size());
    std::vector<MODiagramMemberView> members;
    members.reserve(count);
    for (std::size_t i=0;i<count;++i) {
        MODiagramMemberView member;
        member.orbital_index=level.member_indices.empty()
            ?level.metadata.orbital_index:level.member_indices[i];
        member.inspected_orbital_index=member.orbital_index;
        member.electrons=i<level.member_electrons.size()
            ?level.member_electrons[i]:level.electrons;
        if (i<level.member_spin_counterparts.size()) {
            const auto candidate=level.member_spin_counterparts[i];
            if (std::any_of(data.metadata.begin(),data.metadata.end(),
                [&](const auto& item) { return item.orbital_index==candidate; })) {
                member.spin_counterpart=candidate;
            }
        }
        std::optional<std::size_t> selected;
        if (data.view) selected=data.view->inspected_orbital_index;
        else {
            for (const auto& item:data.metadata) {
                if (item.selected) { selected=item.orbital_index; break; }
            }
            if (!selected && level.metadata.selected) {
                selected=level.metadata.orbital_index;
            }
        }
        if (selected && (*selected==member.orbital_index ||
                         member.spin_counterpart==selected)) {
            member.selected=true;
            member.inspected_orbital_index=*selected;
        }
        members.push_back(member);
    }
    return members;
}

std::optional<std::size_t> mo_diagram_row_for_orbital(
    const MODiagramData& data, const std::size_t index) noexcept {
    for (std::size_t row=0;row<data.levels.size();++row) {
        if (mo_diagram_level_covers_orbital(data.levels[row],index)) return row;
    }
    return std::nullopt;
}

} // namespace cov
