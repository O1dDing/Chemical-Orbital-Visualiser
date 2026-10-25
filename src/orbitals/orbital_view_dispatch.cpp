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
#include "cov/ligand_field.hpp"
#include "cov/point_group_catalog.hpp"

#include <algorithm>
#include <cctype>
#include <string>

namespace cov {

OrbitalRegion classify_orbital_region_legacy(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept;
bool orbital_visible_legacy(
    const std::vector<MolecularOrbital>& orbitals,
    std::size_t index,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) noexcept;
std::vector<std::size_t> visible_orbital_indices_legacy(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings);
std::vector<OrbitalMetadata> build_orbital_metadata_legacy(
    const Wavefunction& wavefunction,
    std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter);

namespace {

bool chemistry_available(const std::vector<MolecularOrbital>& orbitals) {
    return std::any_of(orbitals.begin(),orbitals.end(),
        [](const MolecularOrbital& orbital) {
            return orbital.chemistry.available;
        });
}

bool occupied(const MolecularOrbital& orbital,const double threshold) {
    return static_cast<double>(orbital.occupation)>threshold;
}

} // namespace

DegeneracySettings point_group_limited_degeneracy(
    const Wavefunction& wavefunction,
    DegeneracySettings settings) {
    const auto* definition=find_point_group(wavefunction.point_group_detected);
    std::size_t group_maximum=0u;
    if (definition==nullptr || definition->irreps.empty()) {
        // Gaussian writes the heteronuclear linear group as C*V.  It is not a
        // finite catalogue entry, but its Sigma/Pi/Delta irreps have the same
        // one- or two-dimensional ceiling as Dinfh.  Applying that ceiling is
        // enough to keep two nearby Pi copies from becoming a fake quartet.
        std::string compact;
        compact.reserve(wavefunction.point_group_detected.size());
        for (const unsigned char character:wavefunction.point_group_detected) {
            if (!std::isspace(character) && character!='_' && character!='-') {
                compact.push_back(static_cast<char>(std::toupper(character)));
            }
        }
        if (compact=="C*V" || compact=="CINFV") {
            group_maximum=2u;
        }
    } else {
        for (const auto& irrep:definition->irreps) {
            group_maximum=std::max<std::size_t>(
                group_maximum,irrep.dimension);
        }
    }
    if (group_maximum==0u) return settings;
    settings.maximum_group_size=settings.maximum_group_size==0u
        ?group_maximum
        :std::min(settings.maximum_group_size,group_maximum);
    return settings;
}

bool confidently_deep_core_orbital(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept {
    if (!orbital.chemistry.available ||
        !occupied(orbital,settings.occupation_threshold) ||
        orbital.energy_hartree>=settings.core_energy_cutoff_hartree) {
        return false;
    }

    // A large, dominant projection onto the explicit deep-core reference is
    // required.  Semicore, unresolved/diffuse occupied, and SOMO levels remain
    // visible even when they were not selected by the minimal valence rank.
    const auto& chemistry=orbital.chemistry;
    const double competing=chemistry.semicore_weight+
                           chemistry.valence_weight+
                           chemistry.unresolved_weight;
    return chemistry.deep_core_weight>=0.70 &&
           chemistry.deep_core_weight>=competing+0.10;
}

OrbitalRegion classify_orbital_region(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept {
    if (!orbital.chemistry.available) {
        return classify_orbital_region_legacy(orbital,settings);
    }
    if (confidently_deep_core_orbital(orbital,settings)) {
        return OrbitalRegion::Core;
    }
    if (orbital.chemistry.valence_manifold) return OrbitalRegion::Valence;
    if (!occupied(orbital,settings.occupation_threshold)) {
        return OrbitalRegion::Virtual;
    }
    return OrbitalRegion::Valence;
}

bool orbital_visible(
    const std::vector<MolecularOrbital>& orbitals,
    const std::size_t index,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) noexcept {
    if (index>=orbitals.size()) return false;
    if (!chemistry_available(orbitals)) {
        return orbital_visible_legacy(orbitals,index,frontier,settings);
    }

    const auto& orbital=orbitals[index];
    const auto region=classify_orbital_region(orbital,settings);
    switch (settings.mode) {
        case OrbitalFilterMode::All:
            return true;
        case OrbitalFilterMode::Occupied:
            return occupied(orbital,settings.occupation_threshold);
        case OrbitalFilterMode::Virtual:
            return !occupied(orbital,settings.occupation_threshold);
        case OrbitalFilterMode::Core:
            return region==OrbitalRegion::Core;
        case OrbitalFilterMode::Valence:
        case OrbitalFilterMode::AutoReasonable:
        default:
            if (region==OrbitalRegion::Valence) return true;
            return region==OrbitalRegion::Virtual &&
                   orbital.chemistry.valence_manifold;
    }
}

std::vector<std::size_t> visible_orbital_indices(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) {
    std::vector<std::size_t> result;
    result.reserve(orbitals.size());
    for (std::size_t i=0;i<orbitals.size();++i) {
        if (orbital_visible(orbitals,i,frontier,settings)) {
            result.push_back(i);
        }
    }
    return result;
}

std::vector<OrbitalMetadata> build_orbital_metadata(
    const Wavefunction& wavefunction,
    const std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter) {
    const auto effective_degeneracy=point_group_limited_degeneracy(
        wavefunction,degeneracy);
    auto result=build_orbital_metadata_legacy(
        wavefunction,selected_index,effective_degeneracy,filter);
    const auto frontier=find_frontier_orbitals(
        wavefunction.orbitals,filter.occupation_threshold);
    for (std::size_t i=0;i<result.size();++i) {
        result[i].symmetry=wavefunction.orbitals[i].symmetry;
        result[i].symmetry_view=molecular_orbital_symmetry(wavefunction,i);
        result[i].molecular_member_symmetries={result[i].symmetry_view};
        result[i].region=classify_orbital_region(
            wavefunction.orbitals[i],filter);
        result[i].visible=orbital_visible(
            wavefunction.orbitals,i,frontier,filter);
    }


    // A positive-ion calculation can leave the chemically important vacancy
    // just outside a minimal-rank valence projection.  Preserve the complete
    // lowest virtual degeneracy group for each available spin channel.  This
    // is driven only by producer charge, occupations, energies and degeneracy;
    // no molecule identity is involved.
    if (wavefunction.charge>0 &&
        wavefunction.charge_provenance!=DataProvenance::Unavailable &&
        (filter.mode==OrbitalFilterMode::AutoReasonable ||
         filter.mode==OrbitalFilterMode::Valence)) {
        const auto labels=build_orbital_labels(
            wavefunction.orbitals,effective_degeneracy);
        const auto retain_group=[&](const std::optional<std::size_t>& seed) {
            if (!seed || *seed>=labels.size()) return;
            const std::size_t begin=*seed-labels[*seed].group_member_index;
            const std::size_t end=std::min(
                result.size(),begin+labels[*seed].group_size);
            for (std::size_t i=begin;i<end;++i) result[i].visible=true;
        };
        retain_group(frontier.alpha_lumo);
        if (frontier.separate_spin_sets) retain_group(frontier.beta_lumo);
    }
    apply_local_ligand_field_symmetry(wavefunction,result);
    return result;
}

} // namespace cov
