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

#include <filesystem>
#include <string>
#include "cov/ui.hpp"

namespace cov {

struct FileDialogResult {
    // False when no native picker is available; the caller may offer its own path input.
    bool supported = true;
    bool cancelled = false;
    std::filesystem::path path;
    std::string error;

    [[nodiscard]] bool selected() const noexcept {
        return supported && !cancelled && error.empty() && !path.empty();
    }
};

[[nodiscard]] FileDialogResult open_wavefunction_file_dialog(
    ui::Language language = ui::Language::English, bool nbo_input = false);

// Compatibility alias for the PR #2 viewer call-site. FCHK is now included in
// the dialog and the function may be removed once all callers use the generic name.
[[nodiscard]] inline FileDialogResult open_molden_file_dialog() {
    return open_wavefunction_file_dialog();
}

} // namespace cov
