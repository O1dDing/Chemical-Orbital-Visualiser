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
#include "cov/mo_diagram.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

cov::Atom atom(const char* symbol, int z, double x, double y, double zc) {
    cov::Atom result;
    result.symbol = symbol;
    result.atomic_number = z;
    result.nuclear_charge = static_cast<double>(z);
    result.x = x;
    result.y = y;
    result.z = zc;
    return result;
}

void connect(cov::Wavefunction& wf, std::size_t a, std::size_t b,
             double order = 0.35) {
    cov::BondOrderRecord record;
    record.atom_a = static_cast<std::uint32_t>(a);
    record.atom_b = static_cast<std::uint32_t>(b);
    record.mayer_order = order;
    record.provenance = cov::DataProvenance::Derived;
    wf.bond_orders.push_back(record);
    wf.bond_order_provenance = cov::DataProvenance::Derived;
}

cov::Wavefunction octahedral_sbf6() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Sb", 51, 0.0, 0.0, 0.0));
    constexpr double r = 3.6;
    const std::array<std::array<double, 3>, 6> directions{{
        {{r, 0.0, 0.0}}, {{-r, 0.0, 0.0}},
        {{0.0, r, 0.0}}, {{0.0, -r, 0.0}},
        {{0.0, 0.0, r}}, {{0.0, 0.0, -r}},
    }};
    for (const auto& p : directions) {
        wf.atoms.push_back(atom("F", 9, p[0], p[1], p[2]));
        connect(wf, 0u, wf.atoms.size() - 1u);
    }
    return wf;
}

cov::Wavefunction linear_hf2() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("H", 1, 0.0, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, -2.2, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, 2.2, 0.0, 0.0));
    connect(wf, 0u, 1u, 0.45);
    connect(wf, 0u, 2u, 0.45);
    return wf;
}

} // namespace

int main() {
    auto sb_wavefunction=octahedral_sbf6();
    sb_wavefunction.source=cov::WavefunctionSource::Fchk;
    sb_wavefunction.charge=-1;
    sb_wavefunction.multiplicity=1u;
    sb_wavefunction.alpha_electrons=35u;
    sb_wavefunction.beta_electrons=35u;
    sb_wavefunction.charge_provenance=cov::DataProvenance::Producer;
    sb_wavefunction.multiplicity_provenance=cov::DataProvenance::Producer;
    sb_wavefunction.electron_counts_provenance=cov::DataProvenance::Producer;
    sb_wavefunction.scf_convergence=cov::ScfConvergenceStatus::Converged;
    sb_wavefunction.scf_convergence_provenance=cov::DataProvenance::Producer;
    sb_wavefunction.stability=cov::WavefunctionStabilityStatus::Stable;
    sb_wavefunction.stability_provenance=cov::DataProvenance::Producer;
    sb_wavefunction.stability_detail="integration fixture stable";
    sb_wavefunction.source_title="SbF6 integration fixture";
    sb_wavefunction.source_route="#p integration";
    sb_wavefunction.enrichment_source="fixture.log";
    const auto sb = cov::principal_local_molecular_geometry(sb_wavefunction);
    if (!sb || sb->centre_atom != 0u ||
        sb->geometry_id != cov::GeometryId::Octahedral6 ||
        sb->coordination_number() != 6u) {
        std::cerr << "main-group octahedral geometry regression\n";
        return 1;
    }
    const auto diagram=cov::build_mo_diagram_data(
        sb_wavefunction,cov::MODiagramOptions{});
    if (diagram.local_geometries.empty() ||
        diagram.local_geometries.front().centre_atom!=0u ||
        diagram.local_geometries.front().geometry_id!="OC-6" ||
        diagram.local_geometries.front().point_group!="Oh" ||
        diagram.electronic_state.source!=cov::WavefunctionSource::Fchk ||
        diagram.electronic_state.charge!=-1 ||
        diagram.electronic_state.multiplicity!=1u ||
        diagram.electronic_state.alpha_electrons!=35u ||
        diagram.electronic_state.beta_electrons!=35u ||
        diagram.electronic_state.scf_convergence!=
            cov::ScfConvergenceStatus::Converged ||
        diagram.electronic_state.stability!=
            cov::WavefunctionStabilityStatus::Stable ||
        diagram.electronic_state.stability_detail!="integration fixture stable" ||
        diagram.electronic_state.source_title!="SbF6 integration fixture" ||
        diagram.electronic_state.source_route!="#p integration" ||
        diagram.electronic_state.enrichment_source!="fixture.log") {
        std::cerr << "local geometry or electronic state was not wired into diagram metadata\n";
        return 1;
    }

    const auto hf2 = cov::principal_local_molecular_geometry(linear_hf2());
    if (!hf2 || hf2->centre_atom != 0u ||
        hf2->geometry_id != cov::GeometryId::Linear2 ||
        hf2->coordination_number() != 2u) {
        std::cerr << "linear shared-proton geometry regression\n";
        return 1;
    }

    // Geometry is structural: neither result carries or invents a ligand-field
    // label, and no element/name special case is involved.
    std::cout << "local molecular geometry smoke test passed\n";
    return 0;
}
