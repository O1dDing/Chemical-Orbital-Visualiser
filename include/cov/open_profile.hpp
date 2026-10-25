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
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace cov {
// Optional wall-clock checkpoints. They do not alter parsing or UI state.
class OpenProfile {
    bool enabled_=std::getenv("COV_PROFILE_OPEN")!=nullptr;
    std::chrono::steady_clock::time_point previous_=std::chrono::steady_clock::now();
public:
    void stage(const char* label) {
        if(!enabled_)return;
        const auto now=std::chrono::steady_clock::now();
        std::fprintf(stderr,"COV open detail: %s seconds=%.6f\n",label,
            std::chrono::duration<double>(now-previous_).count());
        std::fflush(stderr);previous_=now;
    }
};
}
