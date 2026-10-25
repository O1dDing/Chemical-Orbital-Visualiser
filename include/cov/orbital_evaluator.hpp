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
#include <memory>
#include <string>

namespace cov {

struct ComputeOptions {
    std::string backend = "auto";
    int device_index = -1;
};

// Only the UI thread touches OpenGL. Portable compute runs in a cancellable
// worker and publishes a complete grid before the texture is replaced.
class OrbitalEvaluator {
public:
    OrbitalEvaluator(const Wavefunction&, ComputeOptions options = {});
    ~OrbitalEvaluator();
    OrbitalEvaluator(const OrbitalEvaluator&) = delete;
    OrbitalEvaluator& operator=(const OrbitalEvaluator&) = delete;
    void attach_gl_texture(unsigned int texture);
    void detach_gl_texture();
    void evaluate(std::size_t mo, const GridBox&, int nx, int ny, int nz);
    void begin_evaluate(std::size_t mo, const GridBox&, int nx, int ny, int nz);
    bool poll();
    void cancel();
    bool busy() const noexcept;
    bool ready() const noexcept;
    const char* device_name() const noexcept;
    double last_kernel_ms() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cov
