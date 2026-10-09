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
