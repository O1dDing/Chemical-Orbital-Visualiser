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

#include <vector>

namespace cov {

// Common production result for both parsers and subsequent analyses. It
// distinguishes a known zero matrix from missing, invalid or failed data.
DensityReconstruction reconstruct_density(const Wavefunction& wavefunction);

// Preserve valid producer matrices, retain their failures, and optionally use
// the shared MO reconstruction when a producer matrix was absent.
void establish_density(Wavefunction& wavefunction, bool reconstruct_missing = true);

// Update trace(P S) and trace(Q S) after the independent AO metric is available.
// Matrix availability itself does not require an overlap matrix.
void update_density_metric_diagnostics(Wavefunction& wavefunction);

// Molden has no mandatory global counters. Derive integer counters only if its
// explicit spin data or declared shared integer model establishes them; real
// fractional occupation sums remain in the density diagnostics.
void derive_molden_electron_counts(Wavefunction& wavefunction);

// Lower-triangular packed AO density, index(i,j) = i*(i+1)/2 + j for i >= j.
// This is the ordinary spin-summed AO density P = C n C^T reconstructed from
// the orbitals already stored in Wavefunction. It intentionally does not
// require an overlap matrix; S enters population/bond-index analyses later.
std::vector<double> reconstruct_total_density_packed(const Wavefunction& wavefunction);

// Reconstruct P_alpha - P_beta when the orbital model carries enough spin/
// occupation information. Restricted open-shell files are handled from the
// alpha/beta electron counts and canonical orbital ordering; unrestricted files
// use their explicit alpha/beta orbital blocks.
std::vector<double> reconstruct_spin_density_packed(const Wavefunction& wavefunction);

} // namespace cov
