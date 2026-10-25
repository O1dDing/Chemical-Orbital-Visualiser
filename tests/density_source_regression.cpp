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
#include "cov/numerical_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using cov::DataProvenance;
using cov::NumericalStatus;
using cov::OrbitalOccupationModel;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void matrix(const std::vector<double>& actual, const std::vector<double>& expected) {
    require(actual.size()==expected.size(), "Density dimensions differ");
    for (std::size_t i=0;i<actual.size();++i)
        require(std::isfinite(actual[i]) && std::abs(actual[i]-expected[i])<2e-14,
                "Density matrix differs from independent two-AO algebra");
}
cov::MolecularOrbital orbital(double occupation, cov::Spin spin,
                             std::vector<double> coefficients, std::size_t source_index) {
    cov::MolecularOrbital mo;
    mo.occupation=occupation;mo.spin=spin;mo.coefficients=std::move(coefficients);
    mo.occupation_provenance=mo.spin_provenance=DataProvenance::Producer;
    mo.source_orbital_index=source_index;
    return mo;
}
cov::Wavefunction canonical() {
    cov::Wavefunction wf;
    wf.source=cov::WavefunctionSource::Fchk;wf.basis_count=2;
    wf.orbital_occupation_model=OrbitalOccupationModel::CanonicalShared;
    wf.orbital_occupation_model_provenance=DataProvenance::Derived;
    wf.alpha_electrons=2;wf.beta_electrons=1;
    wf.electron_counts_provenance=DataProvenance::Producer;
    wf.ao_overlap={1,0,0,1};
    wf.orbitals={orbital(2,cov::Spin::Alpha,{1,0},0),orbital(1,cov::Spin::Alpha,{0,1},1)};
    return wf;
}
}

int main() {
    try {
        // A reordered canonical list retains the original occupied directions.
        auto wf=canonical();
        cov::establish_density(wf);
        matrix(wf.total_density_packed,{2,0,1});matrix(wf.spin_density_packed,{0,0,1});
        require(wf.total_density_diagnostics.metric_trace==3,"Total metric trace changed");
        std::reverse(wf.orbitals.begin(),wf.orbitals.end());
        for(auto& mo:wf.orbitals)for(auto& c:mo.coefficients)c=-c;
        cov::establish_density(wf);
        matrix(wf.total_density_packed,{2,0,1});matrix(wf.spin_density_packed,{0,0,1});

        // Truncating an occupied direction cannot fabricate a full matrix.
        wf.orbitals.erase(wf.orbitals.begin());
        cov::establish_density(wf);
        require(wf.total_density_packed.empty() && wf.spin_density_packed.empty(),
                "Truncated occupied space looks complete");
        require(wf.total_density_diagnostics.status==NumericalStatus::MissingInput,
                "Missing occupied space lacks a status");

        // Complete stored matrices are independent of that missing MO direction.
        wf.total_density_packed={2,0,1};wf.spin_density_packed={0,0,1};
        wf.total_density_provenance=wf.spin_density_provenance=DataProvenance::Producer;
        cov::establish_density(wf);
        matrix(wf.total_density_packed,{2,0,1});matrix(wf.spin_density_packed,{0,0,1});
        require(wf.mo_density_reconstruction.total.packed.empty() &&
                wf.mo_density_reconstruction.spin.packed.empty(),"Raw reconstruction conceals missing MOs");
        require(wf.total_density_provenance==DataProvenance::Producer,"Stored source overwritten");

        // A rejected producer field remains rejected on repeated postprocessing.
        wf=canonical();wf.total_density_packed={2,std::numeric_limits<double>::quiet_NaN(),1};
        wf.total_density_provenance=DataProvenance::Producer;
        cov::establish_density(wf);cov::establish_density(wf);
        require(wf.total_density_packed.empty() && wf.total_density_provenance==DataProvenance::Unavailable,
                "Bad producer matrix was hidden by a derived fallback");
        require(wf.total_density_diagnostics.status==NumericalStatus::InvalidInput,
                "Bad producer matrix lost its reason");
        require(wf.mo_density_reconstruction.total.available(),"Independent reconstruction was lost");

        // Equal spin counts may occupy different directions and have nonzero Q.
        wf={};wf.source=cov::WavefunctionSource::Molden;wf.basis_count=2;
        wf.ao_overlap={1,0,0,1};
        wf.orbitals={orbital(1,cov::Spin::Alpha,{1,0},0),orbital(1,cov::Spin::Beta,{0,1},0)};
        cov::derive_molden_electron_counts(wf);cov::establish_density(wf);
        require(wf.alpha_electrons==1 && wf.beta_electrons==1,"Explicit spin counts wrong");
        matrix(wf.total_density_packed,{1,0,1});matrix(wf.spin_density_packed,{1,0,-1});

        // Fractional alpha/beta occupations remain meaningful without integer counters.
        wf.orbitals[0].occupation=.25;wf.orbitals[1].occupation=.75;
        cov::derive_molden_electron_counts(wf);cov::establish_density(wf);
        require(wf.electron_counts_provenance==DataProvenance::Unavailable,"Fractional counters rounded");
        matrix(wf.total_density_packed,{.25,0,.75});matrix(wf.spin_density_packed,{.25,0,-.75});

        // Missing spin does not erase known total density or invent a zero Q.
        wf={};wf.source=cov::WavefunctionSource::Molden;wf.basis_count=2;
        wf.orbitals={orbital(1,cov::Spin::Alpha,{1,0},0)};
        wf.ao_overlap={1,0,0,1};
        wf.atoms={{"H",1,0,0,0},{"H",1,0,0,2}};
        wf.shells={{0,0,0,0,0,0},{1,0,0,1,0,0}};
        wf.orbitals[0].spin_provenance=DataProvenance::Unavailable;
        cov::derive_molden_electron_counts(wf);cov::establish_density(wf);
        matrix(wf.total_density_packed,{1,0,0});
        require(wf.spin_density_packed.empty() && wf.spin_density_diagnostics.status==NumericalStatus::MissingInput,
                "Unknown spin became a usable density");
        require(wf.orbital_occupation_model==OrbitalOccupationModel::Unspecified,"Incomplete metadata named a determinant");
        cov::derive_bond_and_multicentre_analysis(wf);
        require(wf.bond_orders.empty() && wf.bond_order_provenance==DataProvenance::Unavailable,
                "Full Mayer analysis omitted unknown spin density");

        // Missing Occup retains the actual input absence instead of default zero.
        wf.orbitals[0].occupation_provenance=DataProvenance::Unavailable;
        cov::establish_density(wf);
        require(wf.total_density_packed.empty() && wf.total_density_diagnostics.status==NumericalStatus::MissingInput,
                "Missing Occup was treated as zero");
        std::ostringstream out;cov::write_density_evidence_json(out,wf);
        require(out.str().find("missing-input")!=std::string::npos &&
                out.str().find("\"metric_trace\":null")!=std::string::npos,
                "Density evidence lost explicit missing values");
        std::cout << "Density source, identity, partial information, and spin controls passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';return 1;
    }
}
