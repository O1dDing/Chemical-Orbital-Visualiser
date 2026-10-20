#include "cov/orbital_ui.hpp"

#include "cov/mo_diagram.hpp"
#include "cov/nbo_aomo_text.hpp"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

namespace cov::ui {
namespace {

const char* filter_name(const OrbitalFilterMode mode, const Language language) {
    switch (mode) {
        case OrbitalFilterMode::All: return tr(Text::FilterAll, language);
        case OrbitalFilterMode::Occupied: return tr(Text::FilterOccupied, language);
        case OrbitalFilterMode::Virtual: return tr(Text::FilterVirtual, language);
        case OrbitalFilterMode::Core: return tr(Text::FilterCore, language);
        case OrbitalFilterMode::Valence: return tr(Text::FilterValence, language);
        default: return tr(Text::FilterAuto, language);
    }
}

const char* region_name(const OrbitalRegion region, const Language language) {
    switch (region) {
        case OrbitalRegion::Core: return tr(Text::Core, language);
        case OrbitalRegion::Valence: return tr(Text::Valence, language);
        default: return tr(Text::Virtual, language);
    }
}

const char* spin_name(const Spin spin, const Language language) {
    return tr(spin == Spin::Beta ? Text::Beta : Text::Alpha, language);
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool search_matches(const OrbitalMetadata& item, const char* query) {
    if (!query || !*query) return true;
    const std::string needle = lower_ascii(query);
    std::ostringstream text;
    text << item.raw_mo_number << ' ' << item.display_label << ' '
         << item.symmetry << ' ' << spin_name(item.spin, Language::English) << ' '
         << region_name(item.region, Language::English);
    return lower_ascii(text.str()).find(needle) != std::string::npos;
}

std::optional<std::size_t> previous_occupied(const std::vector<MolecularOrbital>& orbitals,
                                             std::size_t from,
                                             const Spin spin,
                                             const double threshold) {
    if (orbitals.empty()) return std::nullopt;
    if (from > orbitals.size()) from = orbitals.size();
    while (from > 0) {
        --from;
        if (orbitals[from].spin == spin && orbitals[from].occupation > threshold) return from;
    }
    return std::nullopt;
}

std::optional<std::size_t> next_virtual(const std::vector<MolecularOrbital>& orbitals,
                                        std::size_t from,
                                        const Spin spin,
                                        const double threshold) {
    for (std::size_t i = from + 1; i < orbitals.size(); ++i) {
        if (orbitals[i].spin == spin && orbitals[i].occupation <= threshold) return i;
    }
    return std::nullopt;
}

std::optional<std::size_t> frontier_for_selected_spin(const FrontierOrbitals& frontier,
                                                      const Spin spin,
                                                      const bool lumo) {
    if (!frontier.separate_spin_sets) return lumo ? frontier.lumo : frontier.homo;
    if (spin == Spin::Beta) return lumo ? frontier.beta_lumo : frontier.beta_homo;
    return lumo ? frontier.alpha_lumo : frontier.alpha_homo;
}

void quick_nav_button(const char* label,
                      const std::optional<std::size_t> target,
                      OrbitalUIActions& actions,
                      const float width) {
    if (!target) ImGui::BeginDisabled();
    if (ImGui::Button(label, ImVec2(width, 0.0f)) && target) {
        actions.select_orbital = *target;
    }
    if (!target) ImGui::EndDisabled();
}

bool begin_filter_combo(OrbitalUIState& state, const Language language) {
    bool changed = false;
    if (ImGui::BeginCombo("##orbital_filter", filter_name(state.filter.mode, language))) {
        constexpr OrbitalFilterMode modes[] = {
            OrbitalFilterMode::AutoReasonable,
            OrbitalFilterMode::All,
            OrbitalFilterMode::Occupied,
            OrbitalFilterMode::Virtual,
            OrbitalFilterMode::Core,
            OrbitalFilterMode::Valence,
        };
        for (const OrbitalFilterMode mode : modes) {
            const bool selected = state.filter.mode == mode;
            if (ImGui::Selectable(filter_name(mode, language), selected)) {
                state.filter.mode = mode;
                changed = true;
            }
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    return changed;
}

void energy_unit_combo(OrbitalUIState& state) {
    constexpr EnergyUnit units[] = {
        EnergyUnit::Hartree,
        EnergyUnit::ElectronVolt,
        EnergyUnit::JoulePerMol,
        EnergyUnit::KilojoulePerMol,
        EnergyUnit::CaloriePerMol,
        EnergyUnit::KilocaloriePerMol,
    };
    if (ImGui::BeginCombo("##energy_unit", energy_unit_symbol(state.energy_unit))) {
        for (const EnergyUnit unit : units) {
            const bool selected = state.energy_unit == unit;
            if (ImGui::Selectable(energy_unit_symbol(unit), selected)) state.energy_unit = unit;
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

std::string selected_metadata_text(const OrbitalMetadata& item,
                                   const EnergyUnit unit) {
    std::ostringstream out;
    out << "label=" << item.display_label
        << "; raw_mo=" << item.raw_mo_number
        << "; internal_index=" << item.orbital_index
        << "; energy_hartree=" << item.energy_hartree
        << "; energy_display=" << convert_hartree(item.energy_hartree, unit)
        << ' ' << energy_unit_symbol(unit)
        << "; occupation=" << item.occupation
        << "; spin=" << (item.spin == Spin::Beta ? "Beta" : "Alpha")
        << "; symmetry=" << item.symmetry
        << "; degeneracy_size=" << item.degeneracy_size;
    return out.str();
}

struct DiagramPoint {
    const MODiagramLevel* level = nullptr;
    ImVec2 left{};
    ImVec2 right{};
    float y = 0.0f;
};

float map_energy_y(const double energy,
                   const double e_min,
                   const double e_max,
                   const double centre,
                   const double threshold,
                   const EnergyAxisMode mode,
                   const float top,
                   const float bottom) {
    if (mode == EnergyAxisMode::Linear) {
        const double t = (energy - e_min) / std::max(1.0e-12, e_max - e_min);
        return bottom - static_cast<float>(t) * (bottom - top);
    }
    const double warped_min = std::asinh((e_min - centre) / threshold);
    const double warped_max = std::asinh((e_max - centre) / threshold);
    const double warped = std::asinh((energy - centre) / threshold);
    const double t = (warped - warped_min) /
                     std::max(1.0e-12, warped_max - warped_min);
    return bottom - static_cast<float>(t) * (bottom - top);
}

double energy_at_warp_fraction(const double fraction,
                               const double e_min,
                               const double e_max,
                               const double centre,
                               const double threshold,
                               const EnergyAxisMode mode) {
    if (mode == EnergyAxisMode::Linear) return e_min + fraction * (e_max - e_min);
    const double warped_min = std::asinh((e_min - centre) / threshold);
    const double warped_max = std::asinh((e_max - centre) / threshold);
    return centre + threshold *
                        std::sinh(warped_min + fraction * (warped_max - warped_min));
}

int suffix_offset(const std::string& label, const std::size_t group_size) {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    if (suffix < 'a' || suffix > 'z') return 0;
    return static_cast<int>(std::lround(
        (static_cast<int>(suffix - 'a') -
         (static_cast<int>(group_size) - 1) * 0.5) * 90.0));
}

void draw_arrow(ImDrawList* draw,
                const ImVec2 start,
                const bool up,
                const ImU32 colour) {
    const float dy = up ? -13.0f : 13.0f;
    const ImVec2 tip(start.x, start.y + dy);
    draw->AddLine(start, tip, colour, 1.6f);
    const float head = up ? 4.0f : -4.0f;
    draw->AddTriangleFilled(tip,
                            ImVec2(tip.x - 3.5f, tip.y + head),
                            ImVec2(tip.x + 3.5f, tip.y + head),
                            colour);
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          const std::size_t selected_index,
                          OrbitalUIState& state,
                          const Language language,
                          const float ui_scale,
                          OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }

    state.degeneracy.tolerance_hartree =
        std::clamp(state.degeneracy.tolerance_hartree, 1.0e-9, 1.0e-2);
    state.filter.virtual_window_hartree =
        std::clamp(state.filter.virtual_window_hartree, 0.01, 20.0);

    const FrontierOrbitals frontier = find_frontier_orbitals(
        wavefunction.orbitals, state.filter.occupation_threshold);
    const auto metadata = build_orbital_metadata(
        wavefunction, selected_index, state.degeneracy, state.filter);

    const Spin selected_spin = selected_index < wavefunction.orbitals.size()
                                   ? wavefunction.orbitals[selected_index].spin
                                   : Spin::Alpha;
    const auto homo = frontier_for_selected_spin(frontier, selected_spin, false);
    const auto lumo = frontier_for_selected_spin(frontier, selected_spin, true);
    const auto homo_minus_1 = homo
        ? previous_occupied(wavefunction.orbitals, *homo, selected_spin,
                            state.filter.occupation_threshold)
        : std::nullopt;
    const auto lumo_plus_1 = lumo
        ? next_virtual(wavefunction.orbitals, *lumo, selected_spin,
                       state.filter.occupation_threshold)
        : std::nullopt;

    const float gap = ImGui::GetStyle().ItemSpacing.x;
    const float button_width = std::max(52.0f * ui_scale,
        (ImGui::GetContentRegionAvail().x - 3.0f * gap) * 0.25f);
    quick_nav_button(tr(Text::HOMOMinus1, language), homo_minus_1, actions, button_width);
    ImGui::SameLine();
    quick_nav_button(tr(Text::HOMO, language), homo, actions, button_width);
    ImGui::SameLine();
    quick_nav_button(tr(Text::LUMO, language), lumo, actions, button_width);
    ImGui::SameLine();
    quick_nav_button(tr(Text::LUMOPlus1, language), lumo_plus_1, actions, button_width);

    ImGui::TextDisabled("%s", tr(Text::Search, language));
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputTextWithHint("##orbital_search", tr(Text::Search, language),
                             state.search.data(), state.search.size());

    if (ImGui::BeginTable("##orbital_filter_controls", 2,
                          ImGuiTableFlags_SizingStretchSame |
                          ImGuiTableFlags_NoSavedSettings)) {
        ImGui::TableNextColumn();
        ImGui::TextDisabled("%s", tr(Text::Filter, language));
        ImGui::SetNextItemWidth(-1.0f);
        begin_filter_combo(state, language);
        ImGui::TableNextColumn();
        ImGui::TextDisabled("%s", tr(Text::EnergyUnit, language));
        ImGui::SetNextItemWidth(-1.0f);
        energy_unit_combo(state);
        ImGui::EndTable();
    }

    if (state.filter.mode == OrbitalFilterMode::AutoReasonable) {
        float window = static_cast<float>(state.filter.virtual_window_hartree);
        ImGui::TextDisabled("%s (Ha)", tr(Text::HighVirtualWindow, language));
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::SliderFloat("##virtual_window", &window, 0.05f, 5.0f, "%.2f")) {
            state.filter.virtual_window_hartree = window;
        }
    }

    double tolerance = state.degeneracy.tolerance_hartree;
    ImGui::TextDisabled("%s (Ha)", tr(Text::DegeneracyTolerance, language));
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputDouble("##degeneracy_tolerance", &tolerance,
                           1.0e-6, 1.0e-5, "%.2e")) {
        state.degeneracy.tolerance_hartree = std::clamp(tolerance, 1.0e-9, 1.0e-2);
    }
    ImGui::Checkbox(tr(Text::GroupedLabels, language), &state.grouped_labels);

    std::vector<std::size_t> candidates;
    candidates.reserve(metadata.size());
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        if (metadata[i].visible && search_matches(metadata[i], state.search.data())) {
            candidates.push_back(i);
        }
    }
    ImGui::SameLine();
    ImGui::TextDisabled("%s: %zu / %zu", tr(Text::VisibleOrbitals, language),
                        candidates.size(), metadata.size());

    const float table_height = 270.0f * ui_scale;
    if (ImGui::BeginTable("##orbital_table", 4,
                          ImGuiTableFlags_RowBg |
                          ImGuiTableFlags_BordersInnerH |
                          ImGuiTableFlags_ScrollY |
                          ImGuiTableFlags_SizingStretchProp |
                          ImGuiTableFlags_NoSavedSettings,
                          ImVec2(0.0f, table_height))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_WidthFixed, 76.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Energy, language), ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn(tr(Text::Occupation, language), ImGuiTableColumnFlags_WidthFixed, 64.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Symmetry, language), ImGuiTableColumnFlags_WidthFixed, 72.0f * ui_scale);
        ImGui::TableHeadersRow();

        ImGuiListClipper clipper;
        clipper.Begin(static_cast<int>(candidates.size()));
        while (clipper.Step()) {
            for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                const std::size_t index = candidates[static_cast<std::size_t>(row)];
                const auto& item = metadata[index];
                const std::string label = state.grouped_labels
                                              ? item.display_label
                                              : std::to_string(item.raw_mo_number);
                ImGui::PushID(static_cast<int>(index));
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                const bool selected = index == selected_index;
                if (ImGui::Selectable(label.c_str(), selected,
                                      ImGuiSelectableFlags_SpanAllColumns |
                                      ImGuiSelectableFlags_AllowOverlap)) {
                    actions.select_orbital = index;
                }
                const bool hovered = ImGui::IsItemHovered();

                ImGui::TableNextColumn();
                const std::string energy = format_energy(item.energy_hartree,
                                                         state.energy_unit,
                                                         state.energy_unit == EnergyUnit::JoulePerMol ? 1 : 5);
                ImGui::TextUnformatted(energy.c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.2f", item.occupation);
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(item.symmetry.empty() ? "—" : item.symmetry.c_str());

                if (hovered) {
                    ImGui::BeginTooltip();
                    ImGui::Text("%s  %s", tr(Text::MoldenMO, language), label.c_str());
                    ImGui::Text("%s: %zu", tr(Text::RawMO, language), item.raw_mo_number);
                    ImGui::Text("%s: %zu", tr(Text::InternalIndex, language), item.orbital_index);
                    ImGui::Text("%s: %s", tr(Text::Energy, language), energy.c_str());
                    ImGui::Text("Ha: %.10f", item.energy_hartree);
                    ImGui::Text("eV: %.8f", convert_hartree(item.energy_hartree, EnergyUnit::ElectronVolt));
                    ImGui::Text("%s: %.3f", tr(Text::Occupation, language), item.occupation);
                    ImGui::Text("%s: %s", tr(Text::Spin, language), spin_name(item.spin, language));
                    ImGui::Text("%s: %s", tr(Text::Symmetry, language),
                                item.symmetry.empty() ? tr(Text::NoneValue, language) : item.symmetry.c_str());
                    ImGui::Text("%s: %zu", tr(Text::DegenerateSet, language), item.degeneracy_size);
                    ImGui::Text("%s: %s", tr(Text::Region, language), region_name(item.region, language));
                    ImGui::EndTooltip();
                }
                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    if (selected_index < metadata.size()) {
        const auto& selected = metadata[selected_index];
        if (ImGui::Button(tr(Text::CopyMetadata, language))) {
            const std::string text = selected_metadata_text(selected, state.energy_unit);
            ImGui::SetClipboardText(text.c_str());
        }
    }
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         const std::size_t selected_index,
                         OrbitalUIState& state,
                         const Language language,
                         const float ui_scale,
                         OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }

    MODiagramOptions options;
    options.energy_unit = state.energy_unit;
    options.energy_axis_mode = state.energy_axis_mode;
    options.degeneracy = state.degeneracy;
    options.filter = state.filter;
    options.selected_index = selected_index;
    options.neighbourhood = static_cast<std::size_t>(std::max(2, state.diagram_neighbourhood));
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);

    const char* plan_label = tr(Text::SimpleDiagram, language);
    ImGui::TextDisabled("%s", plan_label);
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s", data.plan.machine_reason.c_str());
    }

    ImGui::TextDisabled("%s", tr(Text::EnergyScale, language));
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::LinearEnergyScale, language),
                           state.energy_axis_mode == EnergyAxisMode::Linear)) {
        state.energy_axis_mode = EnergyAxisMode::Linear;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::NonlinearFocus, language),
                           state.energy_axis_mode == EnergyAxisMode::NonlinearFocus)) {
        state.energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    }

    ImGui::SetNextItemWidth(155.0f * ui_scale);
    ImGui::SliderInt("##diagram_neighbourhood", &state.diagram_neighbourhood, 3, 32,
                     "%d", ImGuiSliderFlags_AlwaysClamp);
    ImGui::SameLine();
    ImGui::TextDisabled("%s", tr(Text::AroundSelected, language));

    const float height = 300.0f * ui_scale;
    ImGui::InvisibleButton("##energy_diagram_canvas", ImVec2(-1.0f, height),
                           ImGuiButtonFlags_MouseButtonLeft);
    const ImVec2 p0 = ImGui::GetItemRectMin();
    const ImVec2 p1 = ImGui::GetItemRectMax();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(p0, p1, IM_COL32(12, 18, 27, 235), 7.0f * ui_scale);
    draw->AddRect(p0, p1, IM_COL32(43, 58, 77, 220), 7.0f * ui_scale);

    if (!data.levels.empty()) {
        double e_min = data.levels.front().metadata.energy_hartree;
        double e_max = e_min;
        for (const auto& level : data.levels) {
            e_min = std::min(e_min, level.metadata.energy_hartree);
            e_max = std::max(e_max, level.metadata.energy_hartree);
        }
        if (std::abs(e_max - e_min) < 1.0e-10) {
            e_min -= 0.5;
            e_max += 0.5;
        } else {
            const double pad = (e_max - e_min) * 0.08;
            e_min -= pad;
            e_max += pad;
        }
        double energy_centre = (e_min + e_max) * 0.5;
        for (const auto& level : data.levels) {
            if (level.metadata.selected) {
                energy_centre = level.metadata.energy_hartree;
                break;
            }
        }
        const double energy_threshold = std::max(1.0e-5, (e_max - e_min) * 0.04);

        const float top = p0.y + 30.0f * ui_scale;
        const float bottom = p1.y - 20.0f * ui_scale;
        const float left = p0.x + 48.0f * ui_scale;
        const float right = p1.x - 42.0f * ui_scale;
        const float lane_span = right - left;
        const float axis_x = left - 22.0f * ui_scale;
        draw->AddLine(ImVec2(axis_x, bottom), ImVec2(axis_x, top),
                      IM_COL32(129, 148, 171, 255), 1.4f * ui_scale);
        draw->AddTriangleFilled(ImVec2(axis_x, top - 5.0f * ui_scale),
                                ImVec2(axis_x - 4.0f * ui_scale, top + 3.0f * ui_scale),
                                ImVec2(axis_x + 4.0f * ui_scale, top + 3.0f * ui_scale),
                                IM_COL32(129, 148, 171, 255));
        draw->AddText(ImVec2(p0.x + 6.0f * ui_scale, p0.y + 6.0f * ui_scale),
                      IM_COL32(129, 148, 171, 255), tr(Text::Energy, language));
        draw->AddText(ImVec2(p0.x + 62.0f * ui_scale, p0.y + 6.0f * ui_scale),
                      IM_COL32(100, 118, 140, 255),
                      state.energy_axis_mode == EnergyAxisMode::Linear
                          ? tr(Text::LinearEnergyScale, language)
                          : tr(Text::NonlinearEnergyScale, language));
        for (const double fraction : {0.0, 0.5, 1.0}) {
            const double tick_energy = energy_at_warp_fraction(
                fraction, e_min, e_max, energy_centre, energy_threshold,
                state.energy_axis_mode);
            const float tick_y = map_energy_y(tick_energy, e_min, e_max,
                                              energy_centre, energy_threshold,
                                              state.energy_axis_mode,
                                              top, bottom);
            draw->AddLine(ImVec2(axis_x - 4.0f * ui_scale, tick_y),
                          ImVec2(axis_x + 4.0f * ui_scale, tick_y),
                          IM_COL32(100, 118, 140, 255), 1.0f);
            const std::string tick = format_energy(tick_energy, state.energy_unit, 3);
            draw->AddText(ImVec2(axis_x + 6.0f * ui_scale,
                                 tick_y - ImGui::GetTextLineHeight() * 0.5f),
                          IM_COL32(100, 118, 140, 255), tick.c_str());
        }

        std::vector<DiagramPoint> points;
        points.reserve(data.levels.size());
        for (const auto& level : data.levels) {
            float cx = left + 0.5f * lane_span;
            cx += static_cast<float>(suffix_offset(level.metadata.display_label,
                                                   level.metadata.degeneracy_size)) * ui_scale;
            const float y = map_energy_y(level.metadata.energy_hartree,
                                         e_min, e_max, energy_centre,
                                         energy_threshold, state.energy_axis_mode,
                                         top, bottom);
            const float half = std::min(32.0f * ui_scale, lane_span * 0.30f);

            ImU32 colour = IM_COL32(196, 210, 227, 255);
            if (level.metadata.region == OrbitalRegion::Virtual) colour = IM_COL32(126, 143, 165, 255);
            if (level.homo) colour = IM_COL32(79, 210, 157, 255);
            if (level.lumo) colour = IM_COL32(228, 176, 82, 255);
            if (level.metadata.selected) colour = IM_COL32(92, 151, 255, 255);
            draw->AddLine(ImVec2(cx - half, y), ImVec2(cx + half, y), colour,
                          level.metadata.selected ? 3.2f : 1.8f);

            if (level.electrons.alpha > 0) {
                draw_arrow(draw, ImVec2(cx - 7.0f * ui_scale, y - 2.0f), true,
                           IM_COL32(234, 242, 252, 255));
            }
            if (level.electrons.beta > 0) {
                draw_arrow(draw, ImVec2(cx + 7.0f * ui_scale, y - 2.0f), false,
                           IM_COL32(234, 242, 252, 255));
            }

            if (level.metadata.selected || data.levels.size() <= 20) {
                const std::string label = state.grouped_labels
                                              ? level.metadata.display_label
                                              : std::to_string(level.metadata.raw_mo_number);
                draw->AddText(ImVec2(cx + half + 5.0f * ui_scale,
                                     y - ImGui::GetTextLineHeight() * 0.5f),
                              colour, label.c_str());
            }
            points.push_back({&level, ImVec2(cx - half, y), ImVec2(cx + half, y), y});
        }

        if (ImGui::IsItemHovered()) {
            const ImVec2 mouse = ImGui::GetIO().MousePos;
            const DiagramPoint* nearest = nullptr;
            float best = 1.0e9f;
            for (const auto& point : points) {
                const float xmin = point.left.x - 10.0f * ui_scale;
                const float xmax = point.right.x + 70.0f * ui_scale;
                if (mouse.x < xmin || mouse.x > xmax) continue;
                const float distance = std::abs(mouse.y - point.y);
                if (distance < best && distance <= 7.0f * ui_scale) {
                    best = distance;
                    nearest = &point;
                }
            }
            if (nearest && nearest->level) {
                const auto& level = *nearest->level;
                ImGui::BeginTooltip();
                const std::string label = state.grouped_labels
                                              ? level.metadata.display_label
                                              : std::to_string(level.metadata.raw_mo_number);
                ImGui::Text("MO %s", label.c_str());
                ImGui::Text("%s: %zu", tr(Text::RawMO, language), level.metadata.raw_mo_number);
                const std::string energy = format_energy(level.metadata.energy_hartree,
                                                         state.energy_unit,
                                                         state.energy_unit == EnergyUnit::JoulePerMol ? 1 : 6);
                ImGui::Text("%s: %s", tr(Text::Energy, language), energy.c_str());
                ImGui::Text("Ha: %.10f", level.metadata.energy_hartree);
                ImGui::Text("eV: %.8f", convert_hartree(level.metadata.energy_hartree,
                                                        EnergyUnit::ElectronVolt));
                ImGui::Text("%s: %.3f", tr(Text::Occupation, language), level.metadata.occupation);
                ImGui::Text("%s: %s", tr(Text::Spin, language),
                            spin_name(level.metadata.spin, language));
                ImGui::Text("%s: %s", tr(Text::Symmetry, language),
                            level.metadata.symmetry_view.label.empty() ? tr(Text::NoneValue, language)
                                                            : level.metadata.symmetry_view.label.c_str());
                ImGui::Text("%s: %zu", tr(Text::DegenerateSet, language),
                            level.metadata.degeneracy_size);
                ImGui::Text("%s: %s", tr(Text::OrbitalFamily, language),
                            level.annotation.family == "unavailable"
                                ? tr(Text::NoneValue, language)
                                : level.annotation.family.c_str());
                ImGui::Text("%s: %s", tr(Text::BondingClassLabel, language),
                            level.annotation.bonding_class == BondingClass::Unclassified
                                ? tr(Text::NoneValue, language)
                                : bonding_class_name(level.annotation.bonding_class));
                ImGui::EndTooltip();
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    actions.select_orbital = level.metadata.orbital_index;
                }
            }
        }
    }

    if (ImGui::Button(tr(Text::ExportBundle, language), ImVec2(-1.0f, 0.0f))) {
        actions.export_diagram = true;
    }
    if (ImGui::CollapsingHeader(aomo_text(language,"Analysis data (advanced)"))) {
        if (ImGui::Button(aomo_text(language,"Export analysis data"))) actions.export_analysis=true;
    }
}

} // namespace cov::ui
