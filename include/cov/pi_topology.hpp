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

#include "cov/model.hpp"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace cov {

// Geometry only proposes compatible p directions. Electronic assignment is
// deliberately left to S-metric active-subspace analysis: a geometrical
// candidate must never be presented as proof of a pi bond by itself.
struct PiTopologyOptions {
    double linear_rank_tolerance = 0.08;
    double planar_rank_tolerance = 0.12;
    // Slight pyramidalisation specifically contributed by attached hydrogens
    // must reduce confidence, not sever an otherwise electronically complete
    // amide/resonance p network. Heavy-atom tetrahedral/spiro centres retain
    // the stricter ordinary threshold.
    double hydrogen_pyramidal_rank_tolerance = 0.16;
    double common_direction_tolerance = 0.12;
    double adjacent_direction_cosine = 0.90;
    double cyclic_global_coherence = 0.86;
};

struct OrientedPiNetwork {
    std::vector<std::uint32_t> atoms;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    std::array<double, 3> representative_direction{0.0, 0.0, 1.0};
    double orientation_coherence = 0.0;
    double confidence = 0.0;
    bool cyclic = false;
    bool globally_oriented = false;
};

// bonded_pairs is the strong molecular graph, not a renderer distance graph.
// Hydrogens contribute to the local hybridisation frame but are never p
// centres. Separate orthogonal channels sharing an sp centre remain separate.
[[nodiscard]] std::vector<OrientedPiNetwork> infer_oriented_pi_networks(
    const Wavefunction& wavefunction,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& bonded_pairs,
    const PiTopologyOptions& options = {});

} // namespace cov
