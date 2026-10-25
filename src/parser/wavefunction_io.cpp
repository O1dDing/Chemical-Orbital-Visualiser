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

#include "cov/wavefunction_io.hpp"

#include "cov/bond_analysis.hpp"
#include "cov/density.hpp"
#include "cov/fchk_overlap.hpp"
#include "cov/fchk_nuclear_charge.hpp"
#include "cov/fchk_parser.hpp"
#include "cov/formchk.hpp"
#include "cov/gaussian_log.hpp"
#include "cov/molden_parser.hpp"
#include "cov/orbital_symmetry.hpp"
#include "cov/orbital_chemistry.hpp"
#include "cov/overlap.hpp"
#include "cov/symmetry.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

namespace cov {
namespace {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool looks_like_molden(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string line;
    for (int i = 0; i < 8 && std::getline(in, line); ++i) {
        const auto folded = lower(line);
        if (folded.find("[molden format]") != std::string::npos) return true;
    }
    return false;
}

void postprocess_wavefunction(Wavefunction& wf,
                              const WavefunctionParseOptions& options) {
    for (auto& mo : wf.orbitals) {
        // Format parsers, not this shared post-processing stage, know whether
        // an optional per-MO field was actually present.  In particular, a
        // missing Molden Occup= line must remain unavailable rather than being
        // promoted from its numeric zero default to producer data.
        if (wf.source == WavefunctionSource::Fchk &&
                   mo.occupation_provenance == DataProvenance::Unavailable) {
            mo.occupation_provenance = DataProvenance::Derived;
        }
    }

    // Producer-reported Gaussian point-group information wins when a sibling
    // log/out enrichment supplied it. FCHK-only and Molden-only inputs receive a
    // geometry-derived point group from COV's own operation/permutation engine.
    derive_point_group_from_geometry(wf);

    if (options.keep_density &&
        wf.mo_density_reconstruction.total.diagnostics.status == NumericalStatus::NotComputed) {
        establish_density(wf, options.reconstruct_density_if_missing);
    }

    // Independent AO integrals support rectangular/truncated MO blocks and
    // expose inconsistent producer matrices rather than fitting S to C.
    establish_ao_metric(wf);
    update_density_metric_diagnostics(wf);
    const bool usable_orbitals=orbital_metric_usable(wf);

    // FCHK normally lacks per-MO irreps. Once a validated geometry operation set
    // and AO overlap are available, derive irreps from the transformed AO/MO
    // representation. Producer labels remain immutable and always take priority.
    if (usable_orbitals) {
        (void)derive_orbital_symmetry(wf);
    }

    if (usable_orbitals && wf.total_density_diagnostics.status == NumericalStatus::Available &&
        wf.spin_density_diagnostics.status == NumericalStatus::Available) {
        derive_bond_and_multicentre_analysis(wf);
    }

    // Chemistry is derived last so it can consume the validated AO metric,
    // density-level Mayer indices, multicentre assignments and MO symmetry.
    if (usable_orbitals) {
        derive_orbital_chemistry(wf);
    }
}

} // namespace

Wavefunction parse_wavefunction(const std::filesystem::path& path,
                                const WavefunctionParseOptions& options) {
    const std::string filename = lower(path.filename().string());
    const std::string extension = lower(path.extension().string());

    Wavefunction wf;
    if (extension == ".fchk" || extension == ".fch" || extension == ".chk") {
        FchkParseOptions fchk;
        fchk.max_atoms = options.max_atoms;
        fchk.require_orbitals = options.require_orbitals;
        fchk.keep_density = options.keep_density;
        fchk.reconstruct_density_if_missing = options.reconstruct_density_if_missing;
        if (extension == ".chk") {
            wf = parse_gaussian_chk_via_formchk(path, fchk);
        } else {
            wf = parse_fchk(path, fchk);
            (void)enrich_fchk_nuclear_charges_from_file(wf, path);
            (void)enrich_fchk_overlap_from_file(wf, path);
        }
    } else if (extension == ".molden" ||
               filename.ends_with(".molden.input") || filename.ends_with(".molden.inp") ||
               looks_like_molden(path)) {
        MoldenParseOptions molden;
        molden.max_atoms = options.max_atoms;
        molden.require_orbitals = options.require_orbitals;
        wf = parse_molden_format(path, molden);
        wf.source = WavefunctionSource::Molden;
    } else {
        throw std::runtime_error(
            "Unsupported wavefunction input. Use Gaussian .fchk/.fch/.chk or Molden .molden files.");
    }

    std::filesystem::path enrichment = options.gaussian_log_path;
    if (enrichment.empty() && options.auto_enrich_gaussian_log) {
        enrichment = find_sibling_gaussian_log(path);
    }
    if (!enrichment.empty()) {
        (void)enrich_from_gaussian_log(wf, enrichment);
    }

    postprocess_wavefunction(wf, options);
    return wf;
}

} // namespace cov
