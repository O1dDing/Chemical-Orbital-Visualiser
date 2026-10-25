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

#include "cov/coordination_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string_view>
#include <vector>

namespace {

using cov::CoordinationDirection;

constexpr std::array<double, 9> kTestRotation{
    0.36, -0.48, 0.80,
    0.80, 0.60, 0.00,
    -0.48, 0.64, 0.60,
};

CoordinationDirection rotate(const std::array<double, 9>& matrix,
                             const CoordinationDirection& value) {
    return {
        matrix[0] * value.x + matrix[1] * value.y + matrix[2] * value.z,
        matrix[3] * value.x + matrix[4] * value.y + matrix[5] * value.z,
        matrix[6] * value.x + matrix[7] * value.y + matrix[8] * value.z,
    };
}

CoordinationDirection scale(const CoordinationDirection& value,
                            const double factor) {
    return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(const CoordinationDirection& left,
           const CoordinationDirection& right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

double norm(const CoordinationDirection& value) {
    return std::sqrt(dot(value, value));
}

CoordinationDirection unit(const CoordinationDirection& value) {
    return scale(value, 1.0 / norm(value));
}

bool validate_catalogue() {
    const std::array<std::pair<std::string_view, std::string_view>, 26> expected{{
        {"L-2", "Dinfh"}, {"A-2", "C2v"},
        {"TP-3", "D3h"}, {"TPY-3", "C3v"}, {"TS-3", "C2v"},
        {"T-4", "Td"}, {"SP-4", "D4h"}, {"SS-4", "C2v"},
        {"vTBPY-4", "C3v"},
        {"TBPY-5", "D3h"}, {"SPY-5", "C4v"},
        {"OC-6", "Oh"}, {"TPR-6", "D3h"},
        {"PBPY-7", "D5h"}, {"COC-7", "C3v"}, {"CTPR-7", "C2v"},
        {"SAPR-8", "D4d"}, {"TDD-8", "D2d"}, {"BTPR-8", "C2v"},
        {"CSAPR-9", "C4v"}, {"TCTPR-9", "D3h"},
        {"PPR-10", "D5h"}, {"PAPR-10", "D5d"},
        {"BCSAPR-10", "D4d"}, {"SPC-10", "C2v"}, {"TD-10", "C2v"},
    }};

    const auto& catalog = cov::coordination_geometry_catalog();
    if (catalog.size() != expected.size()) {
        std::cerr << "catalogue size is " << catalog.size() << ", expected 26\n";
        return false;
    }
    std::set<cov::GeometryId> ids;
    std::set<std::string_view> machine_ids;
    for (std::size_t i = 0; i < catalog.size(); ++i) {
        const auto& item = catalog[i];
        if (item.machine_id != expected[i].first ||
            item.point_group != expected[i].second ||
            item.coordination_number < 2 || item.coordination_number > 10 ||
            item.id == cov::GeometryId::Unknown || item.name.empty() ||
            item.reference_note.empty() || !ids.insert(item.id).second ||
            !machine_ids.insert(item.machine_id).second ||
            cov::coordination_geometry_descriptor(item.id) != &item ||
            cov::coordination_geometry_descriptor(item.machine_id) != &item) {
            std::cerr << "invalid catalogue entry at index " << i << '\n';
            return false;
        }
        if (!item.matchable()) {
            std::cerr << "catalogue reference unavailable for "
                      << item.machine_id << '\n';
            return false;
        }
    }
    if (cov::coordination_geometry_descriptor("not-a-shape") != nullptr ||
        cov::coordination_geometry_descriptor(cov::GeometryId::Unknown) != nullptr) {
        std::cerr << "catalogue pending/lookup contract failed\n";
        return false;
    }
    for (std::size_t cn = 2; cn <= 10; ++cn) {
        const auto entries = cov::coordination_geometries_for_cn(cn);
        if (entries.empty() || std::any_of(entries.begin(), entries.end(),
            [&](const auto* item) { return item->coordination_number != cn; })) {
            std::cerr << "CN=" << cn << " catalogue partition failed\n";
            return false;
        }
    }
    return true;
}

bool validate_all_matchable_references() {
    for (const auto& geometry : cov::coordination_geometry_catalog()) {
        if (!geometry.matchable()) continue;
        std::vector<CoordinationDirection> input;
        input.reserve(geometry.reference_directions.size());
        // Reverse ordering, rotate rigidly, and vary radii.  Angular identity,
        // assignment semantics, rigid-frame direction and radial_cv are tested
        // in one pass.
        for (std::size_t i = 0; i < geometry.reference_directions.size(); ++i) {
            const std::size_t reference_index =
                geometry.reference_directions.size() - 1 - i;
            const double radius = 2.0 + 0.07 * static_cast<double>(i);
            input.push_back(scale(
                rotate(kTestRotation, geometry.reference_directions[reference_index]),
                radius));
        }

        const auto match = cov::match_coordination_geometry(input);
        if (!match.best || match.best->id != geometry.id || !match.accepted ||
            match.ambiguous || match.best->angular_rms > 2.0e-6 ||
            match.best->radial_cv <= 0.0 ||
            match.best->assignment.size() != input.size()) {
            std::cerr << "ideal transformed reference did not recover "
                      << geometry.machine_id;
            if (match.best) {
                const auto* got = cov::coordination_geometry_descriptor(match.best->id);
                std::cerr << " (got " << (got ? got->machine_id : "unknown")
                          << ", rms=" << match.best->angular_rms << ')';
            }
            std::cerr << '\n';
            return false;
        }

        std::set<std::size_t> assigned;
        for (std::size_t i = 0; i < input.size(); ++i) {
            const std::size_t reference_index = match.best->assignment[i];
            if (reference_index >= input.size() ||
                !assigned.insert(reference_index).second) {
                std::cerr << "assignment is not input->reference permutation for "
                          << geometry.machine_id << '\n';
                return false;
            }
            const auto reconstructed = rotate(
                match.best->rotation_reference_to_input,
                geometry.reference_directions[reference_index]);
            if (dot(unit(input[i]), unit(reconstructed)) < 1.0 - 1.0e-9) {
                std::cerr << "reference->input rotation failed for "
                          << geometry.machine_id << '\n';
                return false;
            }
        }

        auto perturbed = input;
        for (std::size_t i = 0; i < perturbed.size(); ++i) {
            const double sign = (i % 2 == 0) ? 1.0 : -1.0;
            perturbed[i].x += 0.012 * sign;
            perturbed[i].y -= 0.007 * sign;
            perturbed[i].z += 0.005 * static_cast<double>(i % 3) - 0.005;
        }
        const auto noisy_match = cov::match_coordination_geometry(perturbed);
        if (!noisy_match.best || noisy_match.best->id != geometry.id ||
            !noisy_match.accepted || noisy_match.ambiguous) {
            std::cerr << "small angular perturbation changed "
                      << geometry.machine_id << " classification\n";
            return false;
        }
    }
    return true;
}

cov::Wavefunction shell_wavefunction(
    const std::vector<CoordinationDirection>& directions,
    const double mayer,
    const double distance = 3.0) {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(directions.size() + 1);
    wavefunction.atoms[0].symbol = "Ce";
    wavefunction.atoms[0].atomic_number = 58;
    for (std::size_t i = 0; i < directions.size(); ++i) {
        const auto position = scale(directions[i], distance);
        auto& atom = wavefunction.atoms[i + 1];
        atom.symbol = "F";
        atom.atomic_number = 9;
        atom.x = position.x;
        atom.y = position.y;
        atom.z = position.z;
        wavefunction.bond_orders.push_back({
            0u, static_cast<std::uint32_t>(i + 1), mayer,
            cov::DataProvenance::Derived,
        });
    }
    return wavefunction;
}

bool validate_shell_extraction() {
    const auto* pbpy = cov::coordination_geometry_descriptor("PBPY-7");
    const auto* bcsapr = cov::coordination_geometry_descriptor("BCSAPR-10");
    if (!pbpy || !bcsapr) return false;

    auto seven = shell_wavefunction(pbpy->reference_directions, 0.20);
    const std::size_t first_remote = seven.atoms.size();
    seven.atoms.resize(first_remote + pbpy->reference_directions.size());
    for (std::size_t i = 0; i < pbpy->reference_directions.size(); ++i) {
        const auto position = scale(pbpy->reference_directions[i], 6.0);
        auto& atom = seven.atoms[first_remote + i];
        atom = {"O", 8, position.x, position.y, position.z, 8.0};
        seven.bond_orders.push_back({
            0u, static_cast<std::uint32_t>(first_remote + i), 0.35,
            cov::DataProvenance::Derived,
        });
    }
    const auto seven_shell = cov::extract_coordination_shell(seven, 0);
    const auto seven_match = cov::analyse_coordination_shell(seven_shell);
    if (seven_shell.contacts.size() != 7 ||
        seven_shell.radial_candidate_count != 7 ||
        !seven_match.best || seven_match.best->id != pbpy->id ||
        std::any_of(seven_shell.contacts.begin(), seven_shell.contacts.end(),
                    [&](const auto& item) { return item.atom_index >= first_remote; })) {
        std::cerr << "CN=7 radial shadowing/no-six-truncation failed\n";
        return false;
    }

    auto ten = shell_wavefunction(bcsapr->reference_directions, 0.025);
    const auto ten_shell = cov::extract_coordination_shell(ten, 0);
    const auto ten_match = cov::analyse_coordination_shell(ten_shell);
    if (ten_shell.contacts.size() != 10 || !ten_shell.used_low_mayer_retry ||
        !ten_match.best || ten_match.best->id != bcsapr->id) {
        std::cerr << "CN=10 low-Mayer retry/no-six-truncation failed\n";
        return false;
    }

    std::vector<CoordinationDirection> eleven_directions;
    constexpr double pi = 3.14159265358979323846;
    for (std::size_t i = 0; i < 11; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 11.0;
        eleven_directions.push_back({std::cos(angle), std::sin(angle), 0.15});
    }
    auto eleven = shell_wavefunction(eleven_directions, 0.20);
    const auto eleven_shell = cov::extract_coordination_shell(eleven, 0);
    if (eleven_shell.contacts.size() != 11 ||
        cov::analyse_coordination_shell(eleven_shell).best.has_value()) {
        std::cerr << "shell extractor silently cropped an >10 diagnostic shell\n";
        return false;
    }

    auto negative = shell_wavefunction({}, 0.0);
    negative.atoms.push_back({"F", 9, 2.0, 0.0, 0.0, 9.0});
    negative.atoms.push_back({"F", 9, 0.0, 0.0, 8.0, 9.0});
    negative.bond_orders.push_back({
        0u, 1u, -0.40, cov::DataProvenance::Derived,
    });
    negative.bond_orders.push_back({
        0u, 2u, -0.80, cov::DataProvenance::Derived,
    });
    if (!cov::extract_coordination_shell(negative, 0).contacts.empty()) {
        std::cerr << "negative Mayer evidence became a coordination neighbour\n";
        return false;
    }
    return true;
}

bool validate_diagnostics() {
    const auto* bcsapr = cov::coordination_geometry_descriptor("BCSAPR-10");
    const auto* planar = cov::coordination_geometry_descriptor("TP-3");
    if (!bcsapr || !planar) return false;
    const auto ten = cov::match_coordination_geometry(bcsapr->reference_directions);
    if (!ten.unavailable_candidates.empty() || !ten.best ||
        !ten.runner_up.has_value()) {
        std::cerr << "CN=10 best/runner-up catalogue diagnostics failed\n";
        return false;
    }

    cov::GeometryMatchOptions forced_ambiguity;
    forced_ambiguity.ambiguity_margin = 10.0;
    const auto ambiguous = cov::match_coordination_geometry(
        planar->reference_directions, forced_ambiguity);
    if (!ambiguous.best || !ambiguous.runner_up || !ambiguous.ambiguous) {
        std::cerr << "best/runner-up/ambiguous contract failed\n";
        return false;
    }

    std::vector<CoordinationDirection> invalid{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    if (cov::match_coordination_geometry(invalid).best.has_value()) {
        std::cerr << "zero vector should be rejected\n";
        return false;
    }
    return true;
}

} // namespace

int main() {
    if (!validate_catalogue() ||
        !validate_all_matchable_references() ||
        !validate_shell_extraction() ||
        !validate_diagnostics()) {
        return 1;
    }
    std::cout << "coordination geometry smoke tests passed\n";
    return 0;
}
