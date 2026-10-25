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

#include "cov/density.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/wavefunction_io.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

// Instrument the actual production parser/postprocessor and reconstruction.
// There is no alternate density algorithm or expected result in this program.
namespace {
template<class Range> void numbers(const Range& values) {
    std::cout << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) std::cout << ',';
        first = false;
        cov::numerical_json::number(std::cout, value);
    }
    std::cout << ']';
}

template<class Calculate> void reconstructed(Calculate calculate) {
    try {
        const auto values = calculate();
        std::cout << "{\"returned\":true,\"packed\":";
        numbers(values);
        std::cout << '}';
    } catch (const std::exception& error) {
        std::cout << "{\"returned\":false,\"error\":";
        cov::numerical_json::string(std::cout, error.what());
        std::cout << '}';
    }
}
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    try {
        cov::WavefunctionParseOptions options;
        options.max_atoms = 1000;
        options.keep_density = true;
        options.reconstruct_density_if_missing = true;
        const auto wf = cov::parse_wavefunction(argv[1], options);
        std::cout << std::setprecision(17)
                  << "{\"parse_status\":\"returned\",\"nbasis\":" << wf.basis_count
                  << ",\"source\":" << static_cast<int>(wf.source)
                  << ",\"alpha_electrons\":" << wf.alpha_electrons
                  << ",\"beta_electrons\":" << wf.beta_electrons
                  << ",\"electron_counts_provenance\":" << static_cast<int>(wf.electron_counts_provenance)
                  << ",\"total_density_provenance\":" << static_cast<int>(wf.total_density_provenance)
                  << ",\"spin_density_provenance\":" << static_cast<int>(wf.spin_density_provenance)
                  << ",\"total_density_packed\":";
        numbers(wf.total_density_packed);
        std::cout << ",\"spin_density_packed\":";
        numbers(wf.spin_density_packed);
        const auto metadata=[&](const char* field,const cov::DensityMatrixDiagnostics& info) {
            std::cout << ",\"" << field << "_density_status\":";
            cov::numerical_json::string(std::cout,cov::numerical_status_name(info.status));
            std::cout << ",\"" << field << "_density_reason\":";
            cov::numerical_json::string(std::cout,info.detail);
            std::cout << ",\"" << field << "_density_model\":";
            cov::numerical_json::string(std::cout,cov::orbital_occupation_model_name(info.occupation_model));
        };
        metadata("total",wf.total_density_diagnostics);
        metadata("spin",wf.spin_density_diagnostics);
        std::cout << ",\"density_evidence\":";
        cov::write_density_evidence_json(std::cout,wf);
        std::cout << ",\"overlap\":";
        numbers(wf.ao_overlap);
        std::cout << ",\"numerical_diagnostics\":";
        cov::write_numerical_diagnostics_json(std::cout, wf);
        std::cout << ",\"bond_order_provenance\":" << static_cast<int>(wf.bond_order_provenance)
                  << ",\"bond_orders\":[";
        bool first = true;
        for (const auto& bond : wf.bond_orders) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"atoms\":[" << bond.atom_a << ',' << bond.atom_b << "],\"mayer\":";
            cov::numerical_json::number(std::cout, bond.mayer_order);
            std::cout << ",\"provenance\":" << static_cast<int>(bond.provenance) << '}';
        }
        std::cout << "],\"orbitals\":[";
        first = true;
        for (const auto& mo : wf.orbitals) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"occupation\":";
            cov::numerical_json::number(std::cout, mo.occupation);
            std::cout << ",\"occupation_provenance\":" << static_cast<int>(mo.occupation_provenance)
                      << ",\"spin\":" << static_cast<int>(mo.spin)
                      << ",\"spin_provenance\":" << static_cast<int>(mo.spin_provenance)
                      << ",\"source_orbital_index\":";
            if(mo.source_orbital_index==std::numeric_limits<std::size_t>::max())std::cout<<"null";
            else std::cout<<mo.source_orbital_index;
            std::cout << ",\"spin_source_text\":";
            cov::numerical_json::string(std::cout,mo.spin_source_text);
            std::cout << ",\"coefficients\":";
            numbers(mo.coefficients);
            std::cout << '}';
        }
        std::cout << "],\"production_reconstructed_total\":";
        reconstructed([&] { return cov::reconstruct_total_density_packed(wf); });
        std::cout << ",\"production_reconstructed_spin\":";
        reconstructed([&] { return cov::reconstruct_spin_density_packed(wf); });
        std::cout << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cout << "{\"parse_status\":\"rejected\",\"error\":";
        cov::numerical_json::string(std::cout, error.what());
        std::cout << "}\n";
        return 2;
    }
}
