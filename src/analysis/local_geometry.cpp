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

#include "cov/local_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cov {
namespace {

double geometry_confidence(const GeometryMatch& match) {
    if (!match.accepted || match.ambiguous || !match.best) return 0.0;
    constexpr double accepted_rms = 0.40;
    const double fit = std::clamp(
        1.0 - match.best->angular_rms / accepted_rms, 0.0, 1.0);
    double separation = 1.0;
    if (match.runner_up) {
        separation = std::clamp(
            (match.runner_up->angular_rms - match.best->angular_rms) / 0.12,
            0.0, 1.0);
    }
    const double radial_penalty = 0.15 * std::clamp(
        match.best->radial_cv / 0.80, 0.0, 1.0);
    return std::clamp(
        fit * (0.75 + 0.25 * separation) * (1.0 - radial_penalty),
        0.0, 1.0);
}

LocalGeometryEnvironment make_environment(
    const CoordinationShell& shell,
    const GeometryMatch& match) {
    LocalGeometryEnvironment result;
    result.centre_atom = shell.centre_atom;
    result.ambiguous = match.ambiguous;
    if (!match.accepted || !match.best) return result;
    result.geometry_id = match.best->id;
    result.angular_rms = match.best->angular_rms;
    result.shape_measure = match.best->shape_measure;
    result.radial_cv = match.best->radial_cv;
    result.confidence = geometry_confidence(match);
    result.neighbour_atoms.reserve(shell.contacts.size());
    for (const auto& contact : shell.contacts) {
        result.neighbour_atoms.push_back(contact.atom_index);
    }
    std::sort(result.neighbour_atoms.begin(), result.neighbour_atoms.end());
    return result;
}

} // namespace

std::vector<LocalGeometryEnvironment>
analyse_local_molecular_geometries(const Wavefunction& wavefunction) {
    std::vector<LocalGeometryEnvironment> result;
    result.reserve(wavefunction.atoms.size());
    for (std::size_t centre = 0; centre < wavefunction.atoms.size(); ++centre) {
        const CoordinationShell shell = extract_coordination_shell(
            wavefunction, centre);
        if (shell.contacts.size() < 2u || shell.contacts.size() > 10u) continue;
        const GeometryMatch match = analyse_coordination_shell(shell);
        auto environment = make_environment(shell, match);
        if (environment.available()) result.push_back(std::move(environment));
    }
    return result;
}

std::optional<LocalGeometryEnvironment>
principal_local_molecular_geometry(const Wavefunction& wavefunction) {
    auto environments = analyse_local_molecular_geometries(wavefunction);
    if (environments.empty()) return std::nullopt;
    std::sort(environments.begin(), environments.end(),
              [](const auto& left, const auto& right) {
                  if (left.coordination_number() != right.coordination_number()) {
                      return left.coordination_number() > right.coordination_number();
                  }
                  if (std::abs(left.confidence - right.confidence) > 1.0e-12) {
                      return left.confidence > right.confidence;
                  }
                  return left.angular_rms < right.angular_rms;
              });
    if (environments.size() > 1u) {
        const auto& first = environments[0];
        const auto& second = environments[1];
        const bool same_rank =
            first.coordination_number() == second.coordination_number() &&
            std::abs(first.confidence - second.confidence) < 0.035 &&
            std::abs(first.angular_rms - second.angular_rms) < 0.035;
        if (same_rank) return std::nullopt;
    }
    return environments.front();
}

} // namespace cov
