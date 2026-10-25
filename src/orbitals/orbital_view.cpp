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

#include "cov/orbital_view.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>

namespace cov {
namespace {

constexpr double kHartreeToEV = 27.211386245988;
constexpr double kHartreeToJPerMol = 2625499.6394798255;
constexpr double kJoulePerCalorie = 4.184;

bool occupied(const MolecularOrbital& orbital, const double threshold) noexcept {
    return static_cast<double>(orbital.occupation) > threshold;
}

bool has_beta_orbitals(const std::vector<MolecularOrbital>& orbitals) noexcept {
    return std::any_of(orbitals.begin(), orbitals.end(), [](const MolecularOrbital& mo) {
        return mo.spin == Spin::Beta;
    });
}

std::optional<std::size_t> highest_occupied_for_spin(
    const std::vector<MolecularOrbital>& orbitals,
    const Spin spin,
    const double threshold) {
    std::optional<std::size_t> result;
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (orbitals[i].spin == spin && occupied(orbitals[i], threshold) &&
            (!result || orbitals[i].energy_hartree>
                            orbitals[*result].energy_hartree)) {
            result = i;
        }
    }
    return result;
}

std::optional<std::size_t> lowest_virtual_for_spin(
    const std::vector<MolecularOrbital>& orbitals,
    const Spin spin,
    const double threshold) {
    std::optional<std::size_t> result;
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (orbitals[i].spin == spin && !occupied(orbitals[i], threshold) &&
            (!result || orbitals[i].energy_hartree<
                            orbitals[*result].energy_hartree)) {
            result=i;
        }
    }
    return result;
}

bool is_meaningful_symmetry(const std::string& value) {
    if (value.empty()) return false;
    if (value == "?" || value == "-" || value == "none" || value == "None") return false;
    return true;
}

bool compatible_symmetry(const MolecularOrbital& a,
                         const MolecularOrbital& b,
                         const DegeneracySettings& settings) {
    if (!settings.require_compatible_symmetry) return true;
    if (!is_meaningful_symmetry(a.symmetry) || !is_meaningful_symmetry(b.symmetry)) return true;
    return a.symmetry == b.symmetry;
}

} // namespace

double convert_hartree(const double energy_hartree, const EnergyUnit unit) noexcept {
    switch (unit) {
        case EnergyUnit::ElectronVolt:
            return energy_hartree * kHartreeToEV;
        case EnergyUnit::JoulePerMol:
            return energy_hartree * kHartreeToJPerMol;
        case EnergyUnit::KilojoulePerMol:
            return energy_hartree * kHartreeToJPerMol / 1000.0;
        case EnergyUnit::CaloriePerMol:
            return energy_hartree * kHartreeToJPerMol / kJoulePerCalorie;
        case EnergyUnit::KilocaloriePerMol:
            return energy_hartree * kHartreeToJPerMol / (1000.0 * kJoulePerCalorie);
        default:
            return energy_hartree;
    }
}

const char* energy_unit_symbol(const EnergyUnit unit) noexcept {
    switch (unit) {
        case EnergyUnit::ElectronVolt: return "eV";
        case EnergyUnit::JoulePerMol: return "J/mol";
        case EnergyUnit::KilojoulePerMol: return "kJ/mol";
        case EnergyUnit::CaloriePerMol: return "cal/mol";
        case EnergyUnit::KilocaloriePerMol: return "kcal/mol";
        default: return "Ha";
    }
}

std::string format_energy(const double energy_hartree,
                          const EnergyUnit unit,
                          const int precision) {
    std::ostringstream out;
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setprecision(std::max(0, precision))
        << convert_hartree(energy_hartree, unit)
        << ' ' << energy_unit_symbol(unit);
    return out.str();
}

FrontierOrbitals find_frontier_orbitals(
    const std::vector<MolecularOrbital>& orbitals,
    const double occupation_threshold) {
    FrontierOrbitals result;
    if (orbitals.empty()) return result;

    result.separate_spin_sets = has_beta_orbitals(orbitals);

    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (occupied(orbitals[i], occupation_threshold)) {
            if (!result.homo || orbitals[i].energy_hartree>
                                orbitals[*result.homo].energy_hartree) {
                result.homo=i;
            }
        } else if (!result.lumo || orbitals[i].energy_hartree<
                                  orbitals[*result.lumo].energy_hartree) {
            result.lumo=i;
        }
    }

    result.alpha_homo = highest_occupied_for_spin(orbitals, Spin::Alpha, occupation_threshold);
    result.alpha_lumo = lowest_virtual_for_spin(orbitals, Spin::Alpha, occupation_threshold);
    if (result.separate_spin_sets) {
        result.beta_homo = highest_occupied_for_spin(orbitals, Spin::Beta, occupation_threshold);
        result.beta_lumo = lowest_virtual_for_spin(orbitals, Spin::Beta, occupation_threshold);
    }
    return result;
}

OrbitalRegion classify_orbital_region(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept {
    if (!occupied(orbital, settings.occupation_threshold)) {
        return OrbitalRegion::Virtual;
    }
    return orbital.energy_hartree < settings.core_energy_cutoff_hartree
               ? OrbitalRegion::Core
               : OrbitalRegion::Valence;
}

bool orbital_visible(const std::vector<MolecularOrbital>& orbitals,
                     const std::size_t index,
                     const FrontierOrbitals& frontier,
                     const OrbitalFilterSettings& settings) noexcept {
    if (index >= orbitals.size()) return false;
    const auto& mo = orbitals[index];
    const OrbitalRegion region = classify_orbital_region(mo, settings);

    switch (settings.mode) {
        case OrbitalFilterMode::All:
            return true;
        case OrbitalFilterMode::Occupied:
            return region != OrbitalRegion::Virtual;
        case OrbitalFilterMode::Virtual:
            return region == OrbitalRegion::Virtual;
        case OrbitalFilterMode::Core:
            return region == OrbitalRegion::Core;
        case OrbitalFilterMode::Valence:
            return region == OrbitalRegion::Valence;
        case OrbitalFilterMode::AutoReasonable:
        default:
            if (region != OrbitalRegion::Virtual) return true;
            if (!frontier.lumo.has_value() || *frontier.lumo >= orbitals.size()) return true;
            return mo.energy_hartree <=
                   orbitals[*frontier.lumo].energy_hartree + settings.virtual_window_hartree;
    }
}

std::vector<std::size_t> visible_orbital_indices(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) {
    std::vector<std::size_t> indices;
    indices.reserve(orbitals.size());
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (orbital_visible(orbitals, i, frontier, settings)) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::string degeneracy_suffix(std::size_t member_index) {
    std::string suffix;
    do {
        const std::size_t digit = member_index % 26;
        suffix.push_back(static_cast<char>('a' + digit));
        member_index = member_index / 26;
        if (member_index == 0) break;
        --member_index;
    } while (true);
    std::reverse(suffix.begin(), suffix.end());
    return suffix;
}

std::vector<OrbitalLabel> build_orbital_labels(
    const std::vector<MolecularOrbital>& orbitals,
    const DegeneracySettings& settings) {
    std::vector<OrbitalLabel> labels(orbitals.size());
    for (std::size_t i = 0; i < labels.size(); ++i) {
        labels[i].orbital_index = i;
        labels[i].raw_mo_number = i + 1;
        labels[i].group_base_number = i + 1;
        labels[i].display_label = std::to_string(i + 1);
    }

    if (orbitals.empty()) return labels;
    const double tolerance = std::max(0.0, settings.tolerance_hartree);

    std::size_t begin = 0;
    while (begin < orbitals.size()) {
        std::size_t end = begin + 1;
        while (end < orbitals.size()) {
            if (settings.maximum_group_size > 0u &&
                end - begin >= settings.maximum_group_size) {
                break;
            }
            const auto& first = orbitals[begin];
            const auto& current = orbitals[end];
            if (settings.require_same_spin && current.spin != first.spin) break;
            if (!compatible_symmetry(first, current, settings)) break;
            if (std::abs(current.energy_hartree - first.energy_hartree) > tolerance) break;
            ++end;
        }

        const std::size_t count = end - begin;
        if (count > 1) {
            const std::size_t base = begin + 1;
            for (std::size_t i = begin; i < end; ++i) {
                const std::size_t member = i - begin;
                labels[i].group_base_number = base;
                labels[i].group_member_index = member;
                labels[i].group_size = count;
                labels[i].display_label =
                    std::to_string(base) + "-" + degeneracy_suffix(member);
            }
        }
        begin = end;
    }

    return labels;
}

SymmetrySummary analyse_symmetry_labels(
    const std::vector<MolecularOrbital>& orbitals) {
    SymmetrySummary result;
    if (orbitals.empty()) return result;

    std::set<std::string> distinct;
    std::size_t labelled = 0;
    for (const auto& mo : orbitals) {
        if (is_meaningful_symmetry(mo.symmetry)) {
            ++labelled;
            distinct.insert(mo.symmetry);
        }
    }

    result.labelled_fraction = static_cast<double>(labelled) /
                               static_cast<double>(orbitals.size());
    result.has_meaningful_labels = labelled > 0;
    result.distinct_labels.assign(distinct.begin(), distinct.end());
    return result;
}

DiagramPlan choose_diagram_plan(const Wavefunction& wavefunction,
                                const std::size_t complex_atom_threshold,
                                const double symmetry_coverage_threshold) {
    DiagramPlan plan;
    plan.complex_system = wavefunction.atoms.size() > complex_atom_threshold ||
                          wavefunction.orbitals.size() > 80;

    const SymmetrySummary symmetry = analyse_symmetry_labels(wavefunction.orbitals);
    if (!plan.complex_system) {
        plan.classification = DiagramClassification::Simple;
        plan.machine_reason = "simple-system: direct orbital ordering";
        return plan;
    }

    if (symmetry.labelled_fraction >= symmetry_coverage_threshold) {
        plan.classification = DiagramClassification::SymmetryGrouped;
        plan.machine_reason = "complex-system: producer symmetry labels available";
        return plan;
    }

    plan.classification = DiagramClassification::SalcUnavailable;
    plan.strict_salc_available = false;
    plan.machine_reason =
        "complex-system: strict SALC requires point-group/basis transformation data not present";
    return plan;
}

ElectronGlyphs electron_glyphs_for_orbital(
    const MolecularOrbital& orbital,
    const bool separate_spin_sets) noexcept {
    ElectronGlyphs result;
    const double occupation = std::max(0.0, static_cast<double>(orbital.occupation));

    if (!separate_spin_sets) {
        if (occupation > 0.25) result.alpha = 1;
        if (occupation > 1.25) result.beta = 1;
        return result;
    }

    if (occupation > 0.25) {
        if (orbital.spin == Spin::Beta) result.beta = 1;
        else result.alpha = 1;
    }
    return result;
}

std::vector<OrbitalMetadata> build_orbital_metadata(
    const Wavefunction& wavefunction,
    const std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter) {
    const FrontierOrbitals frontier =
        find_frontier_orbitals(wavefunction.orbitals, filter.occupation_threshold);
    const auto labels = build_orbital_labels(wavefunction.orbitals, degeneracy);

    std::vector<OrbitalMetadata> result;
    result.reserve(wavefunction.orbitals.size());
    for (std::size_t i = 0; i < wavefunction.orbitals.size(); ++i) {
        const auto& mo = wavefunction.orbitals[i];
        OrbitalMetadata item;
        item.orbital_index = i;
        item.raw_mo_number = i + 1;
        item.display_label = labels[i].display_label;
        item.energy_hartree = mo.energy_hartree;
        item.occupation = mo.occupation;
        item.spin = mo.spin;
        item.symmetry = mo.symmetry;
        item.degeneracy_size = labels[i].group_size;
        item.region = classify_orbital_region(mo, filter);
        item.visible = orbital_visible(wavefunction.orbitals, i, frontier, filter);
        item.selected = i == selected_index;
        result.push_back(std::move(item));
    }
    return result;
}

} // namespace cov
