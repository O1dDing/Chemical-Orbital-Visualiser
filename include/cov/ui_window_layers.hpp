#pragma once

#include <imgui.h>

namespace cov::ui {
// Fixed application panels may receive input without covering floating detail
// windows. Popup/menu ordering is still owned by ImGui.
inline constexpr ImGuiWindowFlags background_panel_flags =
    ImGuiWindowFlags_NoBringToFrontOnFocus;
}
