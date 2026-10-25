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
#include <filesystem>
#include <string>

namespace cov {

struct GaussianLogEnrichmentResult {
    bool opened = false;
    std::size_t symmetry_labels_applied = 0;
    bool point_group_applied = false;
    bool spin_squared_applied = false;
    bool scf_diagnostic_applied = false;
    bool stability_diagnostic_applied = false;
    std::string warning;

    [[nodiscard]] bool applied() const noexcept {
        return symmetry_labels_applied != 0 || point_group_applied ||
               spin_squared_applied || scf_diagnostic_applied ||
               stability_diagnostic_applied;
    }
};

// Apply producer-reported metadata from a Gaussian .log/.out sidecar. At this
// stage the enrichment covers final orbital symmetry/point-group metadata and
// explicit producer diagnostics (SCF convergence, stability, and S**2 before
// and after annihilation). It never substitutes the log for the wavefunction
// coefficients carried by FCHK/Molden and never infers a diagnostic merely
// from route keywords.
[[nodiscard]] GaussianLogEnrichmentResult enrich_from_gaussian_log(
    Wavefunction& wavefunction,
    const std::filesystem::path& path);

// Look for a producer output with the same stem. Preference is an adjacent
// .log/.out, followed by an equivalent relative path below a parallel .log
// tree when the wavefunction resides below a .chk directory. Returns an empty
// path when no companion output exists.
[[nodiscard]] std::filesystem::path find_sibling_gaussian_log(
    const std::filesystem::path& wavefunction_path);

} // namespace cov
