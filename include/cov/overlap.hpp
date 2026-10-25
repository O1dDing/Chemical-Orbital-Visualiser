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

#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <vector>

namespace cov {

struct OverlapDerivationResult {
    std::vector<double> matrix; // row-major basis_count x basis_count
    double max_orthonormality_error = std::numeric_limits<double>::quiet_NaN();
    double pivot_ratio = std::numeric_limits<double>::quiet_NaN();

    [[nodiscard]] bool available() const noexcept { return !matrix.empty(); }
};

struct BasisOverlapResult {
    std::vector<double> matrix;
    NumericalStatus status = NumericalStatus::NotComputed;
    std::string detail;
    [[nodiscard]] bool available() const noexcept {
        return status==NumericalStatus::Available && !matrix.empty();
    }
};

// Uses contracted, individually normalized primitive functions in COV's actual
// Cartesian/real-spherical convention. No MO coefficients or fit enter S.
[[nodiscard]] BasisOverlapResult derive_ao_overlap_from_basis(const Wavefunction& wavefunction);

// Full spectrum and all entries of each possibly rectangular C^T S C block.
// Rank is reported, never truncated; absent MO blocks retain MissingInput.
[[nodiscard]] AoMetricDiagnostics inspect_ao_metric(
    const Wavefunction& wavefunction,const std::vector<double>& matrix);

// Establish the production AO metric and retain independent diagnostics of any
// producer matrix. Invalid coefficients are reported, never renormalized.
void establish_ao_metric(Wavefunction& wavefunction);

[[nodiscard]] bool orbital_metric_usable(const Wavefunction& wavefunction) noexcept;

// Recover the AO overlap matrix from a complete orthonormal MO coefficient
// matrix C using C^T S C = I, hence S = C^{-T} C^{-1}. This is format-agnostic
// as a diagnostic comparison only, never as the production AO integral source.
// A complete, nonsingular spin block with exactly basis_count
// orbitals is required; otherwise an unavailable result is returned.
[[nodiscard]] OverlapDerivationResult derive_ao_overlap_from_mos(
    const Wavefunction& wavefunction,
    double relative_pivot_tolerance = 1.0e-11,
    // 1024 keeps the dense recovery workspace comfortably below 128 MiB
    // while covering large, chemically relevant multi-ring/ECP validation
    // systems. Callers analysing still larger matrices may opt in explicitly.
    std::size_t maximum_basis = 1024);

} // namespace cov
