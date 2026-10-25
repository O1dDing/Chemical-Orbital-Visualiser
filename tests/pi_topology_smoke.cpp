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

#include "cov/pi_topology.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

namespace {

cov::Atom atom(const std::uint32_t z, const double x, const double y,
               const double z_coordinate) {
    cov::Atom result;
    result.atomic_number = z;
    result.x = x;
    result.y = y;
    result.z = z_coordinate;
    return result;
}

bool same_atoms(const cov::OrientedPiNetwork& network,
                const std::set<std::uint32_t>& expected) {
    return std::set<std::uint32_t>(network.atoms.begin(), network.atoms.end()) ==
           expected;
}

int acetylene_two_channels() {
    cov::Wavefunction wf;
    wf.atoms = {
        atom(6, -1.2, 0.0, 0.0), atom(6, 1.2, 0.0, 0.0),
        atom(1, -3.2, 0.0, 0.0), atom(1, 3.2, 0.0, 0.0)};
    const std::vector<std::pair<std::uint32_t, std::uint32_t>> bonds = {
        {0, 1}, {0, 2}, {1, 3}};
    const auto networks = cov::infer_oriented_pi_networks(wf, bonds);
    const auto count = std::count_if(networks.begin(), networks.end(),
        [](const auto& network) { return same_atoms(network, {0, 1}); });
    if (count != 2) {
        std::cerr << "acetylene must have two independent perpendicular channels\n";
        return 1;
    }
    const auto first = std::find_if(networks.begin(), networks.end(),
        [](const auto& network) { return same_atoms(network, {0, 1}); });
    const auto second = std::find_if(std::next(first), networks.end(),
        [](const auto& network) { return same_atoms(network, {0, 1}); });
    double axis_dot = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        axis_dot += first->representative_direction[i] *
                    second->representative_direction[i];
    }
    if (std::abs(axis_dot) > 1.0e-6) {
        std::cerr << "acetylene channels are not perpendicular\n";
        return 2;
    }
    return 0;
}

int allene_channels_do_not_merge() {
    cov::Wavefunction wf;
    // C=C=C lies on x. The left CH2 plane is xy (p along z); the right CH2
    // plane is xz (p along y). The central sp carbon admits both directions.
    wf.atoms = {
        atom(6, -2.4, 0.0, 0.0), atom(6, 0.0, 0.0, 0.0),
        atom(6, 2.4, 0.0, 0.0),
        atom(1, -3.1, 1.4, 0.0), atom(1, -3.1, -1.4, 0.0),
        atom(1, 3.1, 0.0, 1.4), atom(1, 3.1, 0.0, -1.4)};
    const std::vector<std::pair<std::uint32_t, std::uint32_t>> bonds = {
        {0, 1}, {1, 2}, {0, 3}, {0, 4}, {2, 5}, {2, 6}};
    const auto networks = cov::infer_oriented_pi_networks(wf, bonds);
    bool left = false;
    bool right = false;
    bool incorrectly_merged = false;
    for (const auto& network : networks) {
        left = left || same_atoms(network, {0, 1});
        right = right || same_atoms(network, {1, 2});
        incorrectly_merged = incorrectly_merged || same_atoms(network, {0, 1, 2});
    }
    if (!left || !right || incorrectly_merged) {
        std::cerr << "allene orthogonal pi channels were merged or lost\n";
        return 3;
    }
    return 0;
}

int planar_ring_and_tub_negative() {
    cov::Wavefunction planar;
    cov::Wavefunction tub;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> bonds;
    constexpr double pi = 3.14159265358979323846;
    for (std::uint32_t i = 0; i < 8; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 8.0;
        planar.atoms.push_back(atom(6, 2.6 * std::cos(angle),
                                    2.6 * std::sin(angle), 0.0));
        tub.atoms.push_back(atom(6, 2.6 * std::cos(angle),
                                 2.6 * std::sin(angle),
                                 (i % 2 == 0) ? 0.95 : -0.95));
    }
    for (std::uint32_t i = 0; i < 8; ++i) {
        const auto next = (i + 1) % 8;
        bonds.emplace_back(i, next);
        const double angle = 2.0 * pi * static_cast<double>(i) / 8.0;
        planar.atoms.push_back(atom(1, 4.6 * std::cos(angle),
                                    4.6 * std::sin(angle), 0.0));
        tub.atoms.push_back(atom(1, 4.6 * std::cos(angle),
                                 4.6 * std::sin(angle),
                                 (i % 2 == 0) ? 1.25 : -1.25));
        bonds.emplace_back(i, 8 + i);
    }
    const auto planar_networks = cov::infer_oriented_pi_networks(planar, bonds);
    const auto tub_networks = cov::infer_oriented_pi_networks(tub, bonds);
    const bool planar_cycle = std::any_of(
        planar_networks.begin(), planar_networks.end(),
        [](const auto& network) {
            return network.cyclic && network.atoms.size() == 8;
        });
    const bool tub_cycle = std::any_of(
        tub_networks.begin(), tub_networks.end(),
        [](const auto& network) {
            return network.cyclic && network.atoms.size() == 8;
        });
    if (!planar_cycle || tub_cycle) {
        std::cerr << "planar/tub cyclic pi topology negative control failed\n";
        return 4;
    }
    return 0;
}

} // namespace

int main() {
    if (const auto code = acetylene_two_channels(); code != 0) return code;
    if (const auto code = allene_channels_do_not_merge(); code != 0) return code;
    if (const auto code = planar_ring_and_tub_negative(); code != 0) return code;
    std::cout << "oriented pi topology smoke test passed\n";
    return 0;
}
