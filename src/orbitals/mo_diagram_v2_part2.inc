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

            return std::pair{centres, electrons};
        }
    }
    return std::nullopt;
}

[[maybe_unused]] std::string annotation_summary(const OrbitalAnnotation& annotation) {
    std::string out;
    if (annotation.family != "unavailable") out += annotation.family;
    if (annotation.bonding_class != BondingClass::Unclassified) {
        if (!out.empty()) out += " / ";
        out += bonding_class_name(annotation.bonding_class);
    }
    if (annotation.multicentre.available) {
        if (!out.empty()) out += " / ";
        out += annotation.multicentre.label;
    }
    if (annotation.delocalised_pi.available) {
        if (!out.empty()) out += " / ";
        out += annotation.delocalised_pi.label;
    }
    return out;
}

int group_member_index(const std::string& label) noexcept {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    if (suffix < 'a' || suffix > 'z') return 0;
    return static_cast<int>(suffix - 'a');
}

std::size_t group_base_number(const OrbitalMetadata& metadata) noexcept {
    const int member = group_member_index(metadata.display_label);
    if (metadata.degeneracy_size <= 1 || member < 0) return metadata.raw_mo_number;
    const auto m = static_cast<std::size_t>(member);
    return metadata.raw_mo_number >= m ? metadata.raw_mo_number - m
                                       : metadata.raw_mo_number;
}

double group_layout_energy(const std::vector<OrbitalMetadata>& metadata,
                           const std::size_t index) {
    if (index >= metadata.size()) return 0.0;
    const auto& item = metadata[index];
    if (item.degeneracy_size <= 1) return item.energy_hartree;
    const std::size_t base = group_base_number(item);
    if (base == 0) return item.energy_hartree;
    const std::size_t begin = base - 1;
    const std::size_t end = std::min(metadata.size(), begin + item.degeneracy_size);
    if (begin >= end) return item.energy_hartree;
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i) sum += metadata[i].energy_hartree;
    return sum / static_cast<double>(end - begin);
}

struct PositionedLevel {
    const MODiagramLevel* level = nullptr;
    double x = 0.0;
    double y = 0.0;
    double label_y = 0.0;
};

double coordinate_min(const EnergyTransform& transform) noexcept {
    if (transform.knots.empty()) return 0.0;
    return transform.knots.front().coordinate;
}

double coordinate_max(const EnergyTransform& transform) noexcept {
    if (transform.knots.empty()) return 1.0;
    return transform.knots.back().coordinate;
}

double map_energy_y(const double layout_energy_hartree,
                    const EnergyTransform& transform,
                    const int height,
                    const int top = 110,
                    const int bottom = 70) {
    const double c = energy_display_coordinate(layout_energy_hartree, transform);
    const double c0 = coordinate_min(transform);
    const double c1 = coordinate_max(transform);
    const double t = (c - c0) / std::max(kCoordinateEpsilon, c1 - c0);
    return static_cast<double>(height - bottom) -
           t * static_cast<double>(height - top - bottom);
}

double energy_at_axis_fraction(const double fraction,
                               const EnergyTransform& transform) noexcept {
    const double c0 = coordinate_min(transform);
    const double c1 = coordinate_max(transform);
    return energy_from_display_coordinate(c0 + fraction * (c1 - c0), transform);
}

std::string orbital_type_plain(const OrbitalAnnotation& annotation);

std::vector<PositionedLevel> position_levels(const MODiagramData& data,
                                             int& width,
                                             const int height) {
    std::vector<PositionedLevel> positioned;
    positioned.reserve(data.levels.size());
    std::vector<DiagramRowFootprint> footprints;
    for (const auto& level:data.levels) {
        const auto members=mo_diagram_member_views(data,level);
        const double strokes=members.empty()?0.0:24.0+31.0*static_cast<double>(members.size()-1);
        // Reserve label space as well as strokes. PNG uses a six-pixel bitmap
        // advance; SVG has the same labels at ten pixels, with a conservative
        // eight-pixel advance. Full glyph validation remains a separate check.
        double label_width=8.0*static_cast<double>(orbital_type_plain(level.annotation).size());
        if (level.annotation.delocalised_pi.available) {
            const auto& pi=level.annotation.delocalised_pi;
            const std::string text="PI "+std::to_string(pi.participating_atoms)+" "+
                std::to_string(static_cast<int>(std::lround(pi.participating_electrons)))+
                " "+compact_pi_topology_suffix(pi);
            label_width=std::max(label_width,8.0*static_cast<double>(text.size()));
        }
        footprints.push_back({map_energy_y(level.layout_energy_hartree,data.energy_transform,height),
            18.0,level.annotation.delocalised_pi.available?35.0:23.0,
            std::max(strokes,label_width)+8.0});
    }
    // Equal or merely nearby *display* footprints need separate lanes, while
    // their exact physical y coordinates and canonical members stay intact.
    constexpr double margin=180.0;
    const auto lanes=layout_diagram_lanes(footprints,std::max(1.0,width-2.0*margin),18.0,2.0);
    width=static_cast<int>(std::ceil(lanes.width+2.0*margin));
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        const double y=footprints[i].y;
        positioned.push_back({&level,margin+lanes.centre_x[i],y,y-7.0});
    }

    // Only text labels move for collision avoidance. The actual level line y
    // remains exactly the energy-axis coordinate.
    std::vector<std::size_t> ids(positioned.size());
    std::iota(ids.begin(), ids.end(), 0);
    std::sort(ids.begin(), ids.end(), [&](const std::size_t a, const std::size_t b) {
        if (positioned[a].label_y != positioned[b].label_y) {
            return positioned[a].label_y < positioned[b].label_y;
        }
        return positioned[a].x < positioned[b].x;
    });
    const double spacing = data.levels.size() > 28 ? 23.0 : 29.0;
    const double top_label = 108.0;
    const double bottom_label = static_cast<double>(height - 48);
    double previous = top_label - spacing;
    for (const std::size_t id : ids) {
        positioned[id].label_y = std::max(positioned[id].label_y, previous + spacing);
        previous = positioned[id].label_y;
    }
    const double overflow = previous - bottom_label;
    if (overflow > 0.0) {
        for (const std::size_t id : ids) positioned[id].label_y -= overflow;
    }
    double next = bottom_label + spacing;
    for (auto it = ids.rbegin(); it != ids.rend(); ++it) {
        positioned[*it].label_y = std::min(positioned[*it].label_y, next - spacing);
        next = positioned[*it].label_y;
    }
    return positioned;
}

struct Canvas {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgba;

    Canvas(const int w, const int h)
        : width(std::max(1, w)),
          height(std::max(1, h)),
          rgba(static_cast<std::size_t>(std::max(1, w)) *
               static_cast<std::size_t>(std::max(1, h)) * 4u, 255u) {
        clear(250, 251, 253, 255);
    }

    void clear(const std::uint8_t r, const std::uint8_t g,
               const std::uint8_t b, const std::uint8_t a) {
        for (std::size_t i = 0; i < rgba.size(); i += 4) {
            rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = a;
        }
    }

    void pixel(const int x, const int y,
               const std::uint8_t r, const std::uint8_t g,
