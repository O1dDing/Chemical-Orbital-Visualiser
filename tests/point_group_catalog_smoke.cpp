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

#include "cov/point_group_catalog.hpp"

#include <array>
#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {

bool has_irrep(const cov::PointGroupDefinition& group,
               std::string_view label,
               std::size_t dimension) {
    for (const auto& irrep : group.irreps) {
        if (irrep.label == label && irrep.dimension == dimension) return true;
    }
    return false;
}

bool validate_group(std::string_view name) {
    const auto* group = cov::find_point_group(name);
    if (!group) {
        std::cerr << "point-group catalog is missing " << name << '\n';
        return false;
    }
    const auto spd = cov::decompose_metal_spd(name);
    if (!spd) {
        std::cerr << "metal s/p/d decomposition is missing for " << name << '\n';
        return false;
    }
    if (cov::irrep_total_dimension(spd->s) != 1u ||
        cov::irrep_total_dimension(spd->p) != 3u ||
        cov::irrep_total_dimension(spd->d) != 5u ||
        spd->total_dimension() != 9u) {
        std::cerr << name << ": s/p/d dimensions are not 1/3/5\n";
        return false;
    }
    for (const auto* decomposition : {&spd->s, &spd->p, &spd->d}) {
        for (const auto& copy : *decomposition) {
            if (!has_irrep(*group, copy.label, copy.dimension)) {
                std::cerr << name << ": decomposition references unknown irrep "
                          << copy.label << '\n';
                return false;
            }
        }
    }
    return true;
}

bool expect_copy(std::string_view group,
                 cov::MetalAOShell shell,
                 std::string_view label,
                 std::size_t expected_multiplicity,
                 std::size_t expected_last_index) {
    const auto decomposition = cov::decompose_metal_ao_shell(group, shell);
    if (!decomposition ||
        cov::irrep_multiplicity(*decomposition, label) != expected_multiplicity) {
        std::cerr << group << ": multiplicity for " << label << " was lost\n";
        return false;
    }
    std::size_t last_index = 0;
    for (const auto& copy : *decomposition) {
        if (copy.label == label) last_index = copy.copy_index;
    }
    if (last_index != expected_last_index) {
        std::cerr << group << ": copy ordinals for " << label
                  << " are not preserved\n";
        return false;
    }
    return true;
}

template <std::size_t N>
bool expect_labels(std::string_view group,
                   cov::MetalAOShell shell,
                   const std::array<std::string_view, N>& expected) {
    const auto decomposition = cov::decompose_metal_ao_shell(group, shell);
    if (!decomposition || decomposition->size() != expected.size()) {
        std::cerr << group << ": wrong number of irrep copies\n";
        return false;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if ((*decomposition)[i].label != expected[i] ||
            (*decomposition)[i].basis_functions.empty()) {
            std::cerr << group << ": wrong irrep-copy sequence at " << i << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main() {
    constexpr std::array required_groups{
        "Dinfh", "C2v", "D3h", "C3v", "Td", "D4h",
        "C4v", "Oh", "D5h", "D5d", "D4d", "D2d"};
    constexpr std::array subgroup_extensions{
        "C1", "Cs", "Ci", "C2", "C2h", "D2", "D2h", "D3d"};

    for (const auto group : required_groups) {
        if (!validate_group(group)) return EXIT_FAILURE;
    }
    for (const auto group : subgroup_extensions) {
        if (!validate_group(group)) return EXIT_FAILURE;
    }

    const auto* oh = cov::find_point_group("O_h");
    const auto* linear = cov::find_point_group("D-infinity-h");
    const auto* gaussian_linear = cov::find_point_group("D*H");
    if (cov::point_group_catalog().size() != 20u ||
        !oh || oh->canonical_name != "Oh" || !oh->centrosymmetric ||
        !linear || linear->canonical_name != "Dinfh" || !linear->linear ||
        gaussian_linear != linear || linear->order != 0u ||
        cov::find_point_group("not-a-group")) {
        std::cerr << "point-group alias/canonical lookup regression\n";
        return EXIT_FAILURE;
    }

    if (!expect_copy("C2v", cov::MetalAOShell::D, "A1", 2u, 2u) ||
        !expect_copy("C3v", cov::MetalAOShell::D, "E", 2u, 2u) ||
        !expect_copy("D3d", cov::MetalAOShell::D, "Eg", 2u, 2u) ||
        !expect_copy("C1", cov::MetalAOShell::P, "A", 3u, 3u) ||
        !expect_copy("Ci", cov::MetalAOShell::D, "Ag", 5u, 5u)) {
        return EXIT_FAILURE;
    }

    if (!expect_labels("Dinfh", cov::MetalAOShell::P,
                       std::array<std::string_view, 2>{"Sigma_u+", "Pi_u"}) ||
        !expect_labels("C2v", cov::MetalAOShell::D,
                       std::array<std::string_view, 5>{"A1", "A1", "A2", "B1", "B2"}) ||
        !expect_labels("D3h", cov::MetalAOShell::D,
                       std::array<std::string_view, 3>{"A1'", "E'", "E''"}) ||
        !expect_labels("C3v", cov::MetalAOShell::D,
                       std::array<std::string_view, 3>{"A1", "E", "E"}) ||
        !expect_labels("Td", cov::MetalAOShell::D,
                       std::array<std::string_view, 2>{"E", "T2"}) ||
        !expect_labels("D4h", cov::MetalAOShell::D,
                       std::array<std::string_view, 4>{"A1g", "B1g", "B2g", "Eg"}) ||
        !expect_labels("C4v", cov::MetalAOShell::D,
                       std::array<std::string_view, 4>{"A1", "B1", "B2", "E"}) ||
        !expect_labels("Oh", cov::MetalAOShell::D,
                       std::array<std::string_view, 2>{"Eg", "T2g"}) ||
        !expect_labels("D5h", cov::MetalAOShell::D,
                       std::array<std::string_view, 3>{"A1'", "E2'", "E1''"}) ||
        !expect_labels("D5d", cov::MetalAOShell::D,
                       std::array<std::string_view, 3>{"A1g", "E2g", "E1g"}) ||
        !expect_labels("D4d", cov::MetalAOShell::D,
                       std::array<std::string_view, 3>{"A1", "E2", "E3"}) ||
        !expect_labels("D2d", cov::MetalAOShell::D,
                       std::array<std::string_view, 4>{"A1", "B1", "B2", "E"})) {
        return EXIT_FAILURE;
    }

    const auto d4d = cov::decompose_metal_spd("D4d");
    const auto d5h = cov::decompose_metal_spd("D5h");
    const auto dinfh = cov::decompose_metal_spd("Dinfh");
    if (!d4d || d4d->p[0].label != "B2" || d4d->p[1].label != "E1" ||
        d4d->d[0].label != "A1" || d4d->d[1].label != "E2" ||
        d4d->d[2].label != "E3" ||
        !d5h || d5h->p[0].label != "A2''" || d5h->d[1].label != "E2'" ||
        !dinfh || dinfh->p[0].label != "Sigma_u+" ||
        dinfh->p[1].label != "Pi_u" || dinfh->d[2].label != "Delta_g") {
        std::cerr << "high/linear point-group s/p/d decomposition regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "point-group catalog smoke test passed\n";
    return EXIT_SUCCESS;
}
