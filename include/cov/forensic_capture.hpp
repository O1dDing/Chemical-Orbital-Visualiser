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
#include <cstddef>
#include <string>

struct ImDrawData;
struct ImFontAtlas;

namespace cov::validation {

// Limits bound evidence size/work. A truncated capture explicitly reports it.
struct RenderedFrameCaptureLimits {
    bool include_glyph_records = false; // Runs/counts remain enabled.
    std::size_t max_glyphs = 12000;
    std::size_t max_runs = 4000;
    std::size_t max_windows = 128;
    // Dense full-basis diagrams submit many non-text line triangles before
    // the floating details window. Keep the scan bounded without starving
    // that final window; glyph/text/JSON budgets remain independent.
    std::size_t max_index_elements = 100000000;
    std::size_t max_text_bytes = 131072;
    std::size_t max_json_bytes = 8388608;
};

// Call immediately after ImGui::Render(), before the next NewFrame() or atlas
// change, with that frame's draw data and its font atlas. No UI hooks required.
// Decodes axis-aligned indexed ImGui font quads, including raw AddText calls.
// Reports geometry/UV evidence, not pixel OCR or proof of raster legibility.
// Missing input characters cannot be recovered: the actual fallback glyph is
// reported. Spaces are geometry-inferred; line/run boundaries are heuristic.
// Never persists texture pixels, complete meshes, or arbitrary GPU buffers.
[[nodiscard]] std::string capture_rendered_frame_json(
    const ImDrawData* draw_data, ImFontAtlas* atlas);
[[nodiscard]] std::string capture_rendered_frame_json(
    const ImDrawData* draw_data, ImFontAtlas* atlas,
    const RenderedFrameCaptureLimits& limits);

} // namespace cov::validation
