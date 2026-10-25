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
#include <iostream>
#include <string>
#include <vector>

namespace {

bool near(const double a, const double b, const double tol) {
    return std::abs(a - b) <= tol;
}

cov::MolecularOrbital mo(const double e,
                         const float occ,
                         const cov::Spin spin = cov::Spin::Alpha,
                         std::string sym = {}) {
    cov::MolecularOrbital value;
    value.energy_hartree = e;
    value.occupation = occ;
    value.spin = spin;
    value.symmetry = std::move(sym);
    return value;
}

} // namespace

int main() {
    if (!near(cov::convert_hartree(1.0, cov::EnergyUnit::ElectronVolt),
              27.211386245988, 1.0e-12)) {
        std::cerr << "Hartree -> eV conversion failed\n";
        return 1;
    }
    if (!near(cov::convert_hartree(1.0, cov::EnergyUnit::KilojoulePerMol),
              2625.4996394798255, 1.0e-9)) {
        std::cerr << "Hartree -> kJ/mol conversion failed\n";
        return 2;
    }
    if (!near(cov::convert_hartree(1.0, cov::EnergyUnit::KilocaloriePerMol),
              627.5094740630558, 1.0e-9)) {
        std::cerr << "Hartree -> kcal/mol conversion failed\n";
        return 3;
    }

    std::vector<cov::MolecularOrbital> orbitals = {
        mo(-2.2, 2.0f, cov::Spin::Alpha, "a1'"),
        mo(-0.236093, 2.0f, cov::Spin::Alpha, "a2''"),
        mo(-0.06804765, 2.0f, cov::Spin::Alpha, "e1''"),
        mo(-0.06804758, 2.0f, cov::Spin::Alpha, "e1''"),
        mo(0.0200, 0.0f, cov::Spin::Alpha, "e2''"),
        mo(0.020000005, 0.0f, cov::Spin::Alpha, "e2''"),
        mo(2.0, 0.0f, cov::Spin::Alpha, "a1'"),
    };

    cov::DegeneracySettings degeneracy;
    degeneracy.tolerance_hartree = 1.0e-5;
    const auto labels = cov::build_orbital_labels(orbitals, degeneracy);
    if (labels.size() != orbitals.size()) return 4;
    if (labels[2].display_label != "3-a" || labels[3].display_label != "3-b") {
        std::cerr << "occupied degeneracy labels failed\n";
        return 5;
    }
    if (labels[4].display_label != "5-a" || labels[5].display_label != "5-b") {
        std::cerr << "virtual degeneracy labels failed\n";
        return 6;
    }
    if (labels[6].display_label != "7") {
        std::cerr << "single-level label failed\n";
        return 7;
    }

    const std::vector<cov::MolecularOrbital> accidental = {
        mo(-0.1, 2.0f, cov::Spin::Alpha, "a1"),
        mo(-0.100000001, 2.0f, cov::Spin::Alpha, "b2"),
    };
    const auto accidental_labels = cov::build_orbital_labels(accidental, degeneracy);
    if (accidental_labels[0].group_size != 1 || accidental_labels[1].group_size != 1) {
        std::cerr << "different producer symmetry labels were incorrectly grouped\n";
        return 8;
    }
    degeneracy.require_compatible_symmetry = false;
    const auto energy_only_labels = cov::build_orbital_labels(accidental, degeneracy);
    if (energy_only_labels[0].display_label != "1-a" ||
        energy_only_labels[1].display_label != "1-b") {
        std::cerr << "explicit energy-only degeneracy mode failed\n";
        return 9;
    }
    degeneracy.require_compatible_symmetry = true;

    const std::vector<cov::MolecularOrbital> adjacent_doublets = {
        mo(-1.0000000, 2.0f, cov::Spin::Alpha, ""),
        mo(-0.9999999, 2.0f, cov::Spin::Alpha, ""),
        mo(-0.9999992, 2.0f, cov::Spin::Alpha, ""),
        mo(-0.9999991, 2.0f, cov::Spin::Alpha, ""),
    };
    degeneracy.maximum_group_size = 2u;
    const auto capped_labels = cov::build_orbital_labels(
        adjacent_doublets, degeneracy);
    if (capped_labels[0].group_size != 2u ||
        capped_labels[1].group_base_number != 1u ||
        capped_labels[2].group_size != 2u ||
        capped_labels[2].group_base_number != 3u) {
        std::cerr << "point-group degeneracy ceiling failed\n";
        return 18;
    }
    degeneracy.maximum_group_size = 0u;

    cov::Wavefunction linear_wavefunction;
    linear_wavefunction.point_group_detected="D*H";
    if (cov::point_group_limited_degeneracy(
            linear_wavefunction,degeneracy).maximum_group_size!=2u) {
        std::cerr << "Gaussian D*H degeneracy ceiling failed\n";
        return 19;
    }
    linear_wavefunction.point_group_detected="C*V";
    if (cov::point_group_limited_degeneracy(
            linear_wavefunction,degeneracy).maximum_group_size!=2u) {
        std::cerr << "Gaussian C*V degeneracy ceiling failed\n";
        return 20;
    }
    cov::Wavefunction octahedral_wavefunction;
    octahedral_wavefunction.point_group_detected="Oh";
    degeneracy.maximum_group_size=6u;
    if (cov::point_group_limited_degeneracy(
            octahedral_wavefunction,degeneracy).maximum_group_size!=3u) {
        std::cerr<<"point-group ceiling did not constrain a looser user cap\n";
        return 21;
    }
    degeneracy.maximum_group_size=2u;
    if (cov::point_group_limited_degeneracy(
            octahedral_wavefunction,degeneracy).maximum_group_size!=2u) {
        std::cerr<<"point-group ceiling did not preserve a stricter user cap\n";
        return 22;
    }
    degeneracy.maximum_group_size=0u;

    const auto frontier = cov::find_frontier_orbitals(orbitals);
    if (!frontier.homo || *frontier.homo != 3 || !frontier.lumo || *frontier.lumo != 4) {
        std::cerr << "frontier detection failed\n";
        return 10;
    }

    cov::OrbitalFilterSettings filter;
    filter.mode = cov::OrbitalFilterMode::AutoReasonable;
    filter.virtual_window_hartree = 1.5;
    const auto visible = cov::visible_orbital_indices(orbitals, frontier, filter);
    if (visible.size() != 6 || visible.back() != 5) {
        std::cerr << "high-energy virtual filter failed\n";
        return 11;
    }

    if (cov::classify_orbital_region(orbitals[0], filter) != cov::OrbitalRegion::Core ||
        cov::classify_orbital_region(orbitals[2], filter) != cov::OrbitalRegion::Valence ||
        cov::classify_orbital_region(orbitals[4], filter) != cov::OrbitalRegion::Virtual) {
        std::cerr << "region classification failed\n";
        return 12;
    }

    cov::Wavefunction wf;
    wf.atoms.resize(20);
    wf.orbitals = orbitals;
    const auto plan = cov::choose_diagram_plan(wf);
    if (plan.classification != cov::DiagramClassification::SymmetryGrouped ||
        plan.strict_salc_available) {
        std::cerr << "complex symmetry grouping plan failed\n";
        return 13;
    }

    wf.orbitals[0].symmetry.clear();
    wf.orbitals[1].symmetry.clear();
    wf.orbitals[2].symmetry.clear();
    wf.orbitals[3].symmetry.clear();
    wf.orbitals[4].symmetry.clear();
    const auto fallback = cov::choose_diagram_plan(wf);
    if (fallback.classification != cov::DiagramClassification::SalcUnavailable ||
        fallback.strict_salc_available) {
        std::cerr << "SALC honesty fallback failed\n";
        return 14;
    }

    const auto restricted_pair = cov::electron_glyphs_for_orbital(mo(-0.5, 2.0f), false);
    if (restricted_pair.alpha != 1 || restricted_pair.beta != 1) {
        std::cerr << "restricted electron glyph failed\n";
        return 15;
    }

    const auto alpha = cov::electron_glyphs_for_orbital(mo(-0.5, 1.0f, cov::Spin::Alpha), true);
    const auto beta = cov::electron_glyphs_for_orbital(mo(-0.5, 1.0f, cov::Spin::Beta), true);
    if (alpha.alpha != 1 || alpha.beta != 0 || beta.alpha != 0 || beta.beta != 1) {
        std::cerr << "spin-resolved electron glyph failed\n";
        return 16;
    }

    // Gaussian UHF data is stored as a complete alpha block followed by a
    // complete beta block.  Global and spin-resolved frontier orbitals must
    // be selected by energy, not by whichever occupied/virtual item happens
    // to occur last/first in that storage order.
    const std::vector<cov::MolecularOrbital> unrestricted = {
        mo(-0.70, 1.0f, cov::Spin::Alpha),
        mo(-0.20, 1.0f, cov::Spin::Alpha),
        mo( 0.30, 0.0f, cov::Spin::Alpha),
        mo(-0.80, 1.0f, cov::Spin::Beta),
        mo(-0.35, 1.0f, cov::Spin::Beta),
        mo( 0.10, 0.0f, cov::Spin::Beta),
    };
    const auto unrestricted_frontier=cov::find_frontier_orbitals(unrestricted);
    if (!unrestricted_frontier.homo || *unrestricted_frontier.homo!=1u ||
        !unrestricted_frontier.lumo || *unrestricted_frontier.lumo!=5u ||
        !unrestricted_frontier.alpha_homo ||
        *unrestricted_frontier.alpha_homo!=1u ||
        !unrestricted_frontier.alpha_lumo ||
        *unrestricted_frontier.alpha_lumo!=2u ||
        !unrestricted_frontier.beta_homo ||
        *unrestricted_frontier.beta_homo!=4u ||
        !unrestricted_frontier.beta_lumo ||
        *unrestricted_frontier.beta_lumo!=5u) {
        std::cerr << "UHF frontier detection depended on block order\n";
        return 17;
    }
    auto reordered=unrestricted;
    std::rotate(reordered.begin(),reordered.begin()+3,reordered.end());
    const auto reordered_frontier=cov::find_frontier_orbitals(reordered);
    if (!reordered_frontier.homo || !reordered_frontier.lumo ||
        !near(reordered[*reordered_frontier.homo].energy_hartree,-0.20,1.0e-12) ||
        !near(reordered[*reordered_frontier.lumo].energy_hartree,0.10,1.0e-12)) {
        std::cerr << "global UHF frontier changed after block reordering\n";
        return 18;
    }

    std::cout << "orbital_view_smoke ok\n";
    return 0;
}
