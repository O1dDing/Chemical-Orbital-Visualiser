#pragma once
#include <algorithm>
#include <cmath>

namespace cov {
struct ViewerRect {
    float x = 0, y = 0, width = 0, height = 0;
    [[nodiscard]] bool contains(float px, float py) const noexcept {
        return px >= x && py >= y && px < x + width && py < y + height;
    }
};
struct FramebufferViewport {
    int x = 0, y = 0, width = 0, height = 0; // OpenGL bottom-left origin, physical pixels.
};
struct ViewerLayout {
    ViewerRect controls;
    ViewerRect scene; // Window coordinates, top-left origin.
    FramebufferViewport framebuffer;
    float window_width = 0, window_height = 0;
    int framebuffer_width = 0, framebuffer_height = 0;
};

inline ViewerLayout viewer_layout(float width, float height, int fb_width, int fb_height,
                                  float ui_scale, bool wide_diagram=false) {
    ViewerLayout layout;
    layout.window_width = width; layout.window_height = height;
    layout.framebuffer_width = fb_width; layout.framebuffer_height = fb_height;
    if (!(width > 0 && height > 0 && fb_width > 0 && fb_height > 0) ||
        !std::isfinite(width) || !std::isfinite(height) || !std::isfinite(ui_scale)) return layout;
    const float scale = std::max(0.5f, ui_scale);
    // Match ImGui's integral window position/size, and keep both regions inside
    // the client area. The scrollable panel never consumes the scene region.
    const float margin = std::floor(std::min(14.0f * scale, std::min(width, height) / 8.0f));
    const float usable_width = std::max(0.0f, width - 3.0f * margin);
    const float preferred = wide_diagram ? std::min(1100.0f*scale,width*0.60f)
        : std::min(540.0f * scale, std::max(370.0f, width * 0.46f));
    const float reserved_scene = std::min(320.0f * scale, usable_width * 0.5f);
    const float panel_width = std::floor(std::min(preferred, usable_width - reserved_scene));
    const float panel_height = std::floor(std::max(0.0f, height - 2.0f * margin));
    layout.controls = {margin, margin, panel_width, panel_height};
    const float left = 2.0f * margin + panel_width;
    const float right = width - margin;
    const float top = margin;
    const float bottom = height - margin;
    const double sx = static_cast<double>(fb_width) / width;
    const double sy = static_cast<double>(fb_height) / height;
    const int x0 = std::clamp(static_cast<int>(std::ceil(left * sx)), 0, fb_width);
    const int x1 = std::clamp(static_cast<int>(std::floor(right * sx)), x0, fb_width);
    const int y0 = std::clamp(static_cast<int>(std::ceil(top * sy)), 0, fb_height);
    const int y1 = std::clamp(static_cast<int>(std::floor(bottom * sy)), y0, fb_height);
    layout.framebuffer = {x0, fb_height - y1, x1 - x0, y1 - y0};
    // Hit testing uses the actual inward-rounded pixel viewport, transformed
    // back to logical coordinates; it does not assume DPI is exactly one.
    layout.scene = {static_cast<float>(x0 / sx), static_cast<float>(y0 / sy),
                    static_cast<float>((x1 - x0) / sx), static_cast<float>((y1 - y0) / sy)};
    return layout;
}

struct SceneProjection { float aspect = 1.0f, tan_half_vertical_fov = 0.0f; };
inline SceneProjection scene_projection(int width, int height, float fov_degrees) {
    const float aspect = height > 0 ? static_cast<float>(width) / height : 1.0f;
    constexpr float pi = 3.14159265358979323846f;
    // The stated field of view belongs to the shorter viewport dimension.
    // Narrow windows therefore retain the same horizontal field of view.
    return {aspect, std::tan(fov_degrees * pi / 360.0f) / std::min(1.0f, std::max(aspect, 1.0e-6f))};
}
} // namespace cov
