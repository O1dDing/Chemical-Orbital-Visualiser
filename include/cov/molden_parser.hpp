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

namespace cov {

struct MoldenParseOptions {
    std::size_t max_atoms = 100;
    bool require_orbitals = true;
};

// The historical implementation is compiled with a source-local macro that
// renames parse_molden -> parse_molden_format. Avoid declaring both names in
// that translation unit, because the macro would otherwise collapse them into
// duplicate declarations with duplicate default arguments.
#ifdef parse_molden
Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options = {});
#else
// Strict Molden-format parser. Most callers should use parse_wavefunction();
// this symbol remains public for format-specific regression tests.
Wavefunction parse_molden_format(const std::filesystem::path& path,
                                 const MoldenParseOptions& options = {});

// Compatibility entry point retained for existing code. It now dispatches
// through the unified input layer, so legacy viewer code can open FCHK without
// losing the existing Molden API surface.
Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options = {});
#endif

} // namespace cov
