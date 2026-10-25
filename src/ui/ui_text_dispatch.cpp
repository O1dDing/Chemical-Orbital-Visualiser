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
#include "cov/orbital_ui_text.hpp"

namespace cov::ui {

const char* tr_legacy(Text key, Language language) noexcept;

const char* tr(Text key, Language language) noexcept {
    switch (key) {
        case Text::MoldenPath:
            return orbital_tr(OrbitalText::CalculationInput,language);
        case Text::IdleHint:
            return orbital_tr(OrbitalText::SupportedInputs,language);
        case Text::MoldenMO:
            return orbital_tr(OrbitalText::SourceMONumber,language);
        case Text::RawMO:
            switch (language) {
                case Language::ChineseSimplified: return "源 MO";
                case Language::Japanese: return "入力 MO";
                case Language::French: return "MO source";
                default: return "Source MO";
            }
        case Text::ClassificationSource:
            return orbital_tr(OrbitalText::LabelSource,language);
        default:
            return tr_legacy(key,language);
    }
}

} // namespace cov::ui
