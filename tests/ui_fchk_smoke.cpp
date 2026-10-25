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

#include "cov/ui.hpp"

#include <cstring>
#include <cstdlib>
#include <iostream>

int main() {
    using cov::ui::Language;
    using cov::ui::Text;

    const char* en_file=cov::ui::tr(Text::MoldenPath,Language::English);
    const char* en_idle=cov::ui::tr(Text::IdleHint,Language::English);
    const char* en_mo=cov::ui::tr(Text::MoldenMO,Language::English);
    const char* zh_file=cov::ui::tr(Text::MoldenPath,Language::ChineseSimplified);
    const char* en_tracking=cov::ui::tr(Text::FrameTracking,Language::English);
    const char* zh_fallback=cov::ui::tr(
        Text::ConservativeFallback,Language::ChineseSimplified);

    if (!en_file || std::strcmp(en_file,"Calculation file or folder")!=0 ||
        !en_idle || std::strstr(en_idle,"FCHK")==nullptr || std::strstr(en_idle,"NBO")==nullptr ||
        !en_mo || std::strcmp(en_mo,"Source MO")!=0 ||
        !zh_file || std::strcmp(zh_file,"计算文件或文件夹")!=0 ||
        !en_tracking || std::strcmp(en_tracking,"Frame continuity") != 0 ||
        !zh_fallback || std::strcmp(zh_fallback,"已使用简化匹配") != 0) {
        std::cerr << "Calculation-input UI wording regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "Calculation-input UI smoke test passed\n";
    return EXIT_SUCCESS;
}
