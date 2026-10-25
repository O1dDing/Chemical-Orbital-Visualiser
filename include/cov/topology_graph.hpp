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
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cov {
// A graph witness, not a bond detector or a complete chemical classification.
// Each closed simple cycle contains hub and they share no other vertex.
struct CycleUnionWitness {
    std::uint32_t hub=0;
    std::array<std::vector<std::uint32_t>,2> cycles;
    bool additional_connection_without_hub=false;
};
std::optional<CycleUnionWitness> find_cycle_union_at_hub(
    std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    std::uint32_t hub);
bool verify_cycle_union(std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    const CycleUnionWitness& witness);

// A stronger, explicitly scoped association with two electronic channels.
// Each ring has at least two atoms exclusive to its channel, has no atom of
// the other channel except hub, and may include intervening non-channel atoms.
// Failure to associate is not proof that the molecule has no ring union.
struct ChannelCycleSearchResult {
    std::optional<CycleUnionWitness> witness;
    bool unscoped_cycle_union_available=false;
    std::array<std::size_t,2> exclusive_channel_atoms{};
    std::uint64_t qualifying_first_cycles_examined=0;
};
ChannelCycleSearchResult find_channel_cycle_union_at_hub(
    std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    std::uint32_t hub,
    const std::array<std::vector<std::uint32_t>,2>& channels);
}
