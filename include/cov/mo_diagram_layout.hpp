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

#include <span>
#include <vector>

namespace cov {

// A row's actual stroke, arrow, label or hit-target footprint in display units.
// Layout never changes y or interprets pixel proximity as physical degeneracy.
struct DiagramRowFootprint {
    double y = 0.0;
    double above = 0.0;
    double below = 0.0;
    double width = 0.0;
};

struct DiagramLaneLayout {
    std::vector<double> centre_x;
    double width = 0.0;
};

// Overlapping vertical footprints receive disjoint horizontal lanes. Separate
// vertical clusters are centred independently; lanes can be reused after a
// preceding footprint ends. The canvas grows when its minimum width cannot
// contain every lane, rather than shrinking strokes or hiding members.
[[nodiscard]] DiagramLaneLayout layout_diagram_lanes(
    std::span<const DiagramRowFootprint> rows,
    double minimum_width,
    double horizontal_gap,
    double vertical_gap = 0.0);

} // namespace cov
