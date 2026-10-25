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

#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace cov {

namespace {

constexpr double kCoordinateEpsilon = 1.0e-12;

bool occupied(const MolecularOrbital& orbital, const double threshold) noexcept {
    return static_cast<double>(orbital.occupation) > threshold;
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) {
        if (c >= 'A' && c <= 'Z') return static_cast<char>(c - 'A' + 'a');
        return static_cast<char>(c);
    });
    return value;
}

bool contains_any(const std::string& text,
                  const std::initializer_list<std::string_view> needles) {
    for (const auto needle : needles) {
        if (text.find(needle) != std::string::npos) return true;
    }
    return false;
}

const char* spin_name(const Spin spin) noexcept {
    return spin == Spin::Beta ? "Beta" : "Alpha";
}

const char* region_name(const OrbitalRegion region) noexcept {
    switch (region) {
        case OrbitalRegion::Core: return "core";
        case OrbitalRegion::Valence: return "valence";
        default: return "virtual";
    }
}

std::string json_escape(std::string_view value) {
    std::ostringstream out;
    for (const unsigned char c : value) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(c) << std::dec << std::setfill(' ');
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\n\r") == std::string::npos) return value;
    std::string escaped = "\"";
    for (const char c : value) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
    }
    escaped += '"';
    return escaped;
}

std::string xml_escape(std::string_view value) {
    std::string out;
    out.reserve(value.size());
    for (const char c : value) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += c; break;
        }
    }
    return out;
}

std::size_t effective_max_levels(const MODiagramOptions& options) noexcept {
    if (options.max_levels > 0) return std::max<std::size_t>(4, options.max_levels);
    return std::clamp<std::size_t>(options.neighbourhood * 2u + 1u, 10u, 48u);
}

std::vector<std::size_t> expand_degenerate_sets(
    const std::vector<std::size_t>& initial,
    const std::vector<OrbitalLabel>& labels) {
    std::set<std::size_t> expanded(initial.begin(), initial.end());
    for (const std::size_t index : initial) {
        if (index >= labels.size() || labels[index].group_size <= 1) continue;
        const std::size_t begin = labels[index].group_base_number - 1;
        const std::size_t end = std::min(labels.size(), begin + labels[index].group_size);
        for (std::size_t i = begin; i < end; ++i) expanded.insert(i);
    }
    return {expanded.begin(), expanded.end()};
}

double median_positive_gap(const std::vector<double>& values) {
    if (values.size() < 2) return 1.0e-3;
    std::vector<double> gaps;
    gaps.reserve(values.size() - 1);
    for (std::size_t i = 1; i < values.size(); ++i) {
        const double gap = values[i] - values[i - 1];
        if (gap > kCoordinateEpsilon) gaps.push_back(gap);
    }
    if (gaps.empty()) return 1.0e-3;
    const std::size_t middle = gaps.size() / 2;
    std::nth_element(gaps.begin(),
                     gaps.begin() + static_cast<std::ptrdiff_t>(middle),
                     gaps.end());
    double result = gaps[middle];
    if (gaps.size() % 2 == 0) {
        const auto max_lower = *std::max_element(
            gaps.begin(), gaps.begin() + static_cast<std::ptrdiff_t>(middle));
        result = 0.5 * (result + max_lower);
    }
    return std::max(result, 1.0e-8);
}

double interpolate(const double x,
                   const double x0,
                   const double x1,
                   const double y0,
                   const double y1) noexcept {
    const double dx = x1 - x0;
    if (std::abs(dx) <= kCoordinateEpsilon) return y0;
    const double t = (x - x0) / dx;
    return y0 + t * (y1 - y0);
}

std::optional<std::pair<std::size_t, std::size_t>> parse_ncme_token(
    const std::string& text) {
    // Accept explicit producer text such as "3c2e", "4c-2e", "3c_4e".
    // Do not infer this descriptor from energy, occupation or a molecule name.
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') continue;
        std::size_t j = i;
        std::size_t centres = 0;
        while (j < text.size() && text[j] >= '0' && text[j] <= '9') {
            centres = centres * 10u + static_cast<std::size_t>(text[j] - '0');
            ++j;
        }
        if (j >= text.size() || text[j] != 'c') continue;
        ++j;
        while (j < text.size() && (text[j] == '-' || text[j] == '_' || text[j] == ' ')) ++j;
        if (j >= text.size() || text[j] < '0' || text[j] > '9') continue;
        std::size_t electrons = 0;
        while (j < text.size() && text[j] >= '0' && text[j] <= '9') {
            electrons = electrons * 10u + static_cast<std::size_t>(text[j] - '0');
            ++j;
        }
        if (j < text.size() && text[j] == 'e' && centres >= 2 && electrons >= 1) {
