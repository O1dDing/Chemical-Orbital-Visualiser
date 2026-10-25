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

#include "cov/compute_abi.h"

#include <algorithm>
#include <cstring>

namespace {
int fake_context = 1;

void message(char* error, size_t capacity, const char* text) noexcept {
    if (!error || !capacity) return;
    const size_t count = std::min(capacity - 1, std::strlen(text));
    std::memcpy(error, text, count);
    error[count] = 0;
}
void* create(int device_index, char* error, size_t capacity) noexcept {
    if (device_index != -1 && device_index != 0) {
        message(error, capacity, "fake HIP device index unavailable");
        return nullptr;
    }
    message(error, capacity, "");
    return &fake_context;
}
void destroy(void*) noexcept {}
int set_terms(void*, const CovGaussianTerm*, uint32_t, char* error, size_t capacity) noexcept {
    message(error, capacity, "");
    return 1;
}
int evaluate(void*, const CovGridRequest*, float*, char* error, size_t capacity) noexcept {
    message(error, capacity, "fake HIP evaluate failure");
    return 0;
}
const char* device_name(void*) noexcept { return "fake HIP GPU"; }
const CovComputeApi api = {COV_COMPUTE_ABI_VERSION, "hip", create, destroy,
                           set_terms, evaluate, device_name};
}

extern "C" COV_COMPUTE_EXPORT const CovComputeApi* cov_get_compute_api(uint32_t version) {
    return version == COV_COMPUTE_ABI_VERSION ? &api : nullptr;
}
