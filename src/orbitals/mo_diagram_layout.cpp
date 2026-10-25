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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cov {

DiagramLaneLayout layout_diagram_lanes(
    const std::span<const DiagramRowFootprint> rows,
    const double minimum_width,
    const double horizontal_gap,
    const double vertical_gap) {
    const auto valid=[](double x) { return std::isfinite(x) && x>=0.0; };
    if (!valid(minimum_width) || !valid(horizontal_gap) || !valid(vertical_gap)) {
        throw std::invalid_argument("Invalid diagram lane dimensions");
    }
    for (const auto& row:rows) {
        if (!std::isfinite(row.y) || !valid(row.above) ||
            !valid(row.below) || !valid(row.width)) {
            throw std::invalid_argument("Invalid diagram row footprint");
        }
    }
    DiagramLaneLayout result;
    result.width=minimum_width;
    result.centre_x.resize(rows.size());
    std::vector<std::size_t> order(rows.size());
    std::iota(order.begin(),order.end(),0);
    const auto top=[&](std::size_t i) { return rows[i].y-rows[i].above; };
    const auto bottom=[&](std::size_t i) { return rows[i].y+rows[i].below; };
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b) { return top(a)<top(b); });
    struct Cluster { std::vector<std::size_t> ids; double width=0.0; };
    std::vector<Cluster> clusters;
    for (std::size_t begin=0;begin<order.size();) {
        std::size_t end=begin+1;
        double cluster_bottom=bottom(order[begin]);
        while (end<order.size() && top(order[end])<=cluster_bottom+vertical_gap) {
            cluster_bottom=std::max(cluster_bottom,bottom(order[end]));
            ++end;
        }
        struct Lane { double bottom=0.0; double width=0.0; };
        std::vector<Lane> lanes;
        std::vector<std::size_t> membership;
        Cluster cluster;
        for (auto position=begin;position<end;++position) {
            const auto id=order[position];
            auto lane=std::find_if(lanes.begin(),lanes.end(),[&](const auto& candidate) {
                return candidate.bottom+vertical_gap<top(id);
            });
            const auto index=static_cast<std::size_t>(lane-lanes.begin());
            if (lane==lanes.end()) lanes.push_back({bottom(id),rows[id].width});
            else { lane->bottom=bottom(id); lane->width=std::max(lane->width,rows[id].width); }
            membership.push_back(index);
            cluster.ids.push_back(id);
        }
        std::vector<double> centres;
        for (const auto& lane:lanes) {
            if (!centres.empty()) cluster.width+=horizontal_gap;
            centres.push_back(cluster.width+0.5*lane.width);
            cluster.width+=lane.width;
        }
        for (std::size_t i=0;i<cluster.ids.size();++i) {
            result.centre_x[cluster.ids[i]]=centres[membership[i]];
        }
        result.width=std::max(result.width,cluster.width);
        clusters.push_back(std::move(cluster));
        begin=end;
    }
    for (const auto& cluster:clusters) {
        for (const auto id:cluster.ids) result.centre_x[id]+=0.5*(result.width-cluster.width);
    }
    return result;
}

} // namespace cov
