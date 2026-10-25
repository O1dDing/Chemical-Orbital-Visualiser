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

#include "cov/bond_analysis.hpp"
#include "cov/density.hpp"
#include "cov/overlap.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

bool close(double a, double b, double eps = 1.0e-6) {
    return std::abs(a - b) <= eps;
}

cov::Wavefunction make_h2() {
    cov::Wavefunction wf;
    wf.atoms = {{"H", 1, 0.0, 0.0, -0.7}, {"H", 1, 0.0, 0.0, 0.7}};
    wf.basis_count = 2;
    wf.shells = {
        {0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0},
    };

    constexpr double s = 0.2;
    const double cb = 1.0 / std::sqrt(2.0 * (1.0 + s));
    const double ca = 1.0 / std::sqrt(2.0 * (1.0 - s));

    cov::MolecularOrbital bonding;
    bonding.energy_hartree = -0.5;
    bonding.occupation = 2.0f;
    bonding.spin = cov::Spin::Alpha;
    bonding.coefficients = {static_cast<float>(cb), static_cast<float>(cb)};

    cov::MolecularOrbital antibonding;
    antibonding.energy_hartree = 0.2;
    antibonding.occupation = 0.0f;
    antibonding.spin = cov::Spin::Alpha;
    antibonding.coefficients = {static_cast<float>(ca), static_cast<float>(-ca)};

    wf.orbitals = {bonding, antibonding};
    wf.orbital_occupation_model = cov::OrbitalOccupationModel::SharedIntegerDeterminant;
    wf.orbital_occupation_model_provenance = cov::DataProvenance::Derived;
    for (auto& mo : wf.orbitals) {
        mo.occupation_provenance = cov::DataProvenance::Derived;
        mo.spin_provenance = cov::DataProvenance::Derived;
    }
    wf.alpha_electrons = 1;
    wf.beta_electrons = 1;
    return wf;
}

cov::Wavefunction make_three_center() {
    cov::Wavefunction wf;
    wf.atoms = {
        {"H", 1, -1.0, 0.0, 0.0},
        {"H", 1, 0.0, 0.0, 0.0},
        {"H", 1, 1.0, 0.0, 0.0},
    };
    wf.basis_count = 3;
    wf.shells = {
        {0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0},
        {2, 0, 0, 2, 0, 0},
    };

    const double r3 = std::sqrt(3.0);
    const double r2 = std::sqrt(2.0);
    const double r6 = std::sqrt(6.0);

    cov::MolecularOrbital a;
    a.occupation = 2.0f;
    a.coefficients = {
        static_cast<float>(1.0 / r3),
        static_cast<float>(1.0 / r3),
        static_cast<float>(1.0 / r3),
    };

    cov::MolecularOrbital b;
    b.occupation = 0.0f;
    b.coefficients = {
        static_cast<float>(1.0 / r2),
        static_cast<float>(-1.0 / r2),
        0.0f,
    };

    cov::MolecularOrbital c;
    c.occupation = 0.0f;
    c.coefficients = {
        static_cast<float>(1.0 / r6),
        static_cast<float>(1.0 / r6),
        static_cast<float>(-2.0 / r6),
    };

    wf.orbitals = {a, b, c};
    wf.orbital_occupation_model = cov::OrbitalOccupationModel::SharedIntegerDeterminant;
    wf.orbital_occupation_model_provenance = cov::DataProvenance::Derived;
    for (auto& mo : wf.orbitals) {
        mo.occupation_provenance = cov::DataProvenance::Derived;
        mo.spin_provenance = cov::DataProvenance::Derived;
    }
    wf.alpha_electrons = 1;
    wf.beta_electrons = 1;
    return wf;
}

} // namespace

int main() {
    {
        auto wf = make_h2();
        const auto overlap = cov::derive_ao_overlap_from_mos(wf);
        if (!overlap.available() || overlap.matrix.size() != 4u ||
            !close(overlap.matrix[0], 1.0) ||
            !close(overlap.matrix[1], 0.2) ||
            !close(overlap.matrix[2], 0.2) ||
            !close(overlap.matrix[3], 1.0) ||
            overlap.max_orthonormality_error > 1.0e-6) {
            std::cerr << "AO overlap recovery regression\n";
            return EXIT_FAILURE;
        }

        wf.ao_overlap = overlap.matrix;
        wf.ao_overlap_provenance = cov::DataProvenance::Derived;
        cov::establish_density(wf);
        cov::derive_bond_and_multicentre_analysis(wf);

        if (wf.bond_orders.size() != 1u ||
            !close(wf.bond_orders[0].mayer_order, 1.0, 2.0e-5) ||
            wf.bond_order_provenance != cov::DataProvenance::Derived) {
            std::cerr << "Mayer H2 bond-order regression\n";
            return EXIT_FAILURE;
        }
    }

    {
        auto wf = make_three_center();
        const auto overlap = cov::derive_ao_overlap_from_mos(wf);
        if (!overlap.available()) {
            std::cerr << "three-centre overlap recovery regression\n";
            return EXIT_FAILURE;
        }
        wf.ao_overlap = overlap.matrix;
        wf.ao_overlap_provenance = cov::DataProvenance::Derived;
        cov::establish_density(wf);
        cov::derive_bond_and_multicentre_analysis(wf);

        if (wf.multicentre_candidates.size() != 1u ||
            wf.multicentre_candidates[0].orbital_index != 0u ||
            wf.multicentre_candidates[0].atoms.size() != 3u) {
            std::cerr << "multicentre participation regression\n";
            return EXIT_FAILURE;
        }
        for (double weight : wf.multicentre_candidates[0].participation) {
            if (!close(weight, 1.0 / 3.0, 2.0e-5)) {
                std::cerr << "multicentre participation weight regression\n";
                return EXIT_FAILURE;
            }
        }
    }

    std::cout << "wavefunction analysis smoke test passed\n";
    return EXIT_SUCCESS;
}
