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

#include "cov/point_group_catalog.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <string>
#include <utility>

namespace cov {
namespace {

using PGF = PointGroupFamily;
using Shell = MetalAOShell;

template <typename T, std::size_t N>
constexpr std::span<const T> as_span(const std::array<T, N>& values) noexcept {
    return {values.data(), values.size()};
}

constexpr std::array<std::string_view, 1> c1_aliases{"1"};
constexpr std::array<std::string_view, 1> cs_aliases{"C_s"};
constexpr std::array<std::string_view, 1> ci_aliases{"C_i"};
constexpr std::array<std::string_view, 1> c2_aliases{"C_2"};
constexpr std::array<std::string_view, 1> c2h_aliases{"C_2h"};
constexpr std::array<std::string_view, 1> c2v_aliases{"C_2v"};
constexpr std::array<std::string_view, 1> c3v_aliases{"C_3v"};
constexpr std::array<std::string_view, 1> c4v_aliases{"C_4v"};
constexpr std::array<std::string_view, 1> d2_aliases{"D_2"};
constexpr std::array<std::string_view, 1> d2h_aliases{"D_2h"};
constexpr std::array<std::string_view, 1> d2d_aliases{"D_2d"};
constexpr std::array<std::string_view, 1> d3d_aliases{"D_3d"};
constexpr std::array<std::string_view, 1> d3h_aliases{"D_3h"};
constexpr std::array<std::string_view, 1> d4h_aliases{"D_4h"};
constexpr std::array<std::string_view, 1> d4d_aliases{"D_4d"};
constexpr std::array<std::string_view, 1> d5h_aliases{"D_5h"};
constexpr std::array<std::string_view, 1> d5d_aliases{"D_5d"};
constexpr std::array<std::string_view, 2> td_aliases{"T_d", "tetrahedral"};
constexpr std::array<std::string_view, 2> oh_aliases{"O_h", "octahedral"};
constexpr std::array<std::string_view, 6> dinfh_aliases{
    "D_inf_h", "D-infinity-h", "D_infinity_h", "D∞h",
    "linear-centrosymmetric", "D*H"};

constexpr std::array<IrrepDefinition, 1> c1_irreps{{{"A", 1}}};
constexpr std::array<IrrepDefinition, 2> cs_irreps{{{"A'", 1}, {"A''", 1}}};
constexpr std::array<IrrepDefinition, 2> ci_irreps{{{"Ag", 1}, {"Au", 1}}};
constexpr std::array<IrrepDefinition, 2> c2_irreps{{{"A", 1}, {"B", 1}}};
constexpr std::array<IrrepDefinition, 4> c2h_irreps{{
    {"Ag", 1}, {"Bg", 1}, {"Au", 1}, {"Bu", 1}}};
constexpr std::array<IrrepDefinition, 4> c2v_irreps{{
    {"A1", 1}, {"A2", 1}, {"B1", 1}, {"B2", 1}}};
constexpr std::array<IrrepDefinition, 3> c3v_irreps{{
    {"A1", 1}, {"A2", 1}, {"E", 2}}};
constexpr std::array<IrrepDefinition, 5> c4v_irreps{{
    {"A1", 1}, {"A2", 1}, {"B1", 1}, {"B2", 1}, {"E", 2}}};
constexpr std::array<IrrepDefinition, 4> d2_irreps{{
    {"A", 1}, {"B1", 1}, {"B2", 1}, {"B3", 1}}};
constexpr std::array<IrrepDefinition, 8> d2h_irreps{{
    {"Ag", 1}, {"B1g", 1}, {"B2g", 1}, {"B3g", 1},
    {"Au", 1}, {"B1u", 1}, {"B2u", 1}, {"B3u", 1}}};
constexpr std::array<IrrepDefinition, 5> d2d_irreps{{
    {"A1", 1}, {"A2", 1}, {"B1", 1}, {"B2", 1}, {"E", 2}}};
constexpr std::array<IrrepDefinition, 6> d3d_irreps{{
    {"A1g", 1}, {"A2g", 1}, {"Eg", 2},
    {"A1u", 1}, {"A2u", 1}, {"Eu", 2}}};
constexpr std::array<IrrepDefinition, 6> d3h_irreps{{
    {"A1'", 1}, {"A2'", 1}, {"E'", 2},
    {"A1''", 1}, {"A2''", 1}, {"E''", 2}}};
constexpr std::array<IrrepDefinition, 10> d4h_irreps{{
    {"A1g", 1}, {"A2g", 1}, {"B1g", 1}, {"B2g", 1}, {"Eg", 2},
    {"A1u", 1}, {"A2u", 1}, {"B1u", 1}, {"B2u", 1}, {"Eu", 2}}};
constexpr std::array<IrrepDefinition, 7> d4d_irreps{{
    {"A1", 1}, {"A2", 1}, {"B1", 1}, {"B2", 1},
    {"E1", 2}, {"E2", 2}, {"E3", 2}}};
constexpr std::array<IrrepDefinition, 8> d5h_irreps{{
    {"A1'", 1}, {"A2'", 1}, {"E1'", 2}, {"E2'", 2},
    {"A1''", 1}, {"A2''", 1}, {"E1''", 2}, {"E2''", 2}}};
constexpr std::array<IrrepDefinition, 8> d5d_irreps{{
    {"A1g", 1}, {"A2g", 1}, {"E1g", 2}, {"E2g", 2},
    {"A1u", 1}, {"A2u", 1}, {"E1u", 2}, {"E2u", 2}}};
constexpr std::array<IrrepDefinition, 5> td_irreps{{
    {"A1", 1}, {"A2", 1}, {"E", 2}, {"T1", 3}, {"T2", 3}}};
constexpr std::array<IrrepDefinition, 10> oh_irreps{{
    {"A1g", 1}, {"A2g", 1}, {"Eg", 2}, {"T1g", 3}, {"T2g", 3},
    {"A1u", 1}, {"A2u", 1}, {"Eu", 2}, {"T1u", 3}, {"T2u", 3}}};
constexpr std::array<IrrepDefinition, 8> dinfh_irreps{{
    {"Sigma_g+", 1}, {"Sigma_g-", 1}, {"Sigma_u+", 1}, {"Sigma_u-", 1},
    {"Pi_g", 2}, {"Pi_u", 2}, {"Delta_g", 2}, {"Delta_u", 2}}};

const std::array<PointGroupDefinition, 20> catalog{{
    {"C1", PGF::Trivial, 1, false, false, as_span(c1_aliases), as_span(c1_irreps)},
    {"Cs", PGF::Reflection, 2, false, false, as_span(cs_aliases), as_span(cs_irreps)},
    {"Ci", PGF::Reflection, 2, false, true, as_span(ci_aliases), as_span(ci_irreps)},
    {"C2", PGF::Cyclic, 2, false, false, as_span(c2_aliases), as_span(c2_irreps)},
    {"C2h", PGF::Cyclic, 4, false, true, as_span(c2h_aliases), as_span(c2h_irreps)},
    {"C2v", PGF::Cyclic, 4, false, false, as_span(c2v_aliases), as_span(c2v_irreps)},
    {"C3v", PGF::Cyclic, 6, false, false, as_span(c3v_aliases), as_span(c3v_irreps)},
    {"C4v", PGF::Cyclic, 8, false, false, as_span(c4v_aliases), as_span(c4v_irreps)},
    {"D2", PGF::Dihedral, 4, false, false, as_span(d2_aliases), as_span(d2_irreps)},
    {"D2h", PGF::Dihedral, 8, false, true, as_span(d2h_aliases), as_span(d2h_irreps)},
    {"D2d", PGF::Dihedral, 8, false, false, as_span(d2d_aliases), as_span(d2d_irreps)},
    {"D3d", PGF::Dihedral, 12, false, true, as_span(d3d_aliases), as_span(d3d_irreps)},
    {"D3h", PGF::Dihedral, 12, false, false, as_span(d3h_aliases), as_span(d3h_irreps)},
    {"D4h", PGF::Dihedral, 16, false, true, as_span(d4h_aliases), as_span(d4h_irreps)},
    {"D4d", PGF::Dihedral, 16, false, false, as_span(d4d_aliases), as_span(d4d_irreps)},
    {"D5h", PGF::Dihedral, 20, false, false, as_span(d5h_aliases), as_span(d5h_irreps)},
    {"D5d", PGF::Dihedral, 20, false, true, as_span(d5d_aliases), as_span(d5d_irreps)},
    {"Td", PGF::Polyhedral, 24, false, false, as_span(td_aliases), as_span(td_irreps)},
    {"Oh", PGF::Polyhedral, 48, false, true, as_span(oh_aliases), as_span(oh_irreps)},
    {"Dinfh", PGF::Linear, 0, true, true, as_span(dinfh_aliases), as_span(dinfh_irreps)},
}};

struct CopySpec {
    std::string_view label;
    std::uint8_t dimension;
    std::string_view basis;
};

struct SPDRecord {
    std::string_view group;
    std::span<const CopySpec> s;
    std::span<const CopySpec> p;
    std::span<const CopySpec> d;
};

#define COV_COPY(LABEL, DIMENSION, BASIS) CopySpec{LABEL, DIMENSION, BASIS}

constexpr std::array<CopySpec, 1> s_a{{COV_COPY("A", 1, "s")}};
constexpr std::array<CopySpec, 1> s_ap{{COV_COPY("A'", 1, "s")}};
constexpr std::array<CopySpec, 1> s_ag{{COV_COPY("Ag", 1, "s")}};
constexpr std::array<CopySpec, 1> s_a1{{COV_COPY("A1", 1, "s")}};
constexpr std::array<CopySpec, 1> s_a1p{{COV_COPY("A1'", 1, "s")}};
constexpr std::array<CopySpec, 1> s_a1g{{COV_COPY("A1g", 1, "s")}};
constexpr std::array<CopySpec, 1> s_sigma_g{{COV_COPY("Sigma_g+", 1, "s")}};

constexpr std::array<CopySpec, 3> c1_p{{
    COV_COPY("A", 1, "p_z"), COV_COPY("A", 1, "p_x"), COV_COPY("A", 1, "p_y")}};
constexpr std::array<CopySpec, 5> c1_d{{
    COV_COPY("A", 1, "d_z2"), COV_COPY("A", 1, "d_x2-y2"),
    COV_COPY("A", 1, "d_xy"), COV_COPY("A", 1, "d_xz"),
    COV_COPY("A", 1, "d_yz")}};

// Cs convention: the mirror plane is xy.
constexpr std::array<CopySpec, 3> cs_p{{
    COV_COPY("A'", 1, "p_x"), COV_COPY("A'", 1, "p_y"),
    COV_COPY("A''", 1, "p_z")}};
constexpr std::array<CopySpec, 5> cs_d{{
    COV_COPY("A'", 1, "d_z2"), COV_COPY("A'", 1, "d_x2-y2"),
    COV_COPY("A'", 1, "d_xy"), COV_COPY("A''", 1, "d_xz"),
    COV_COPY("A''", 1, "d_yz")}};

constexpr std::array<CopySpec, 3> ci_p{{
    COV_COPY("Au", 1, "p_z"), COV_COPY("Au", 1, "p_x"),
    COV_COPY("Au", 1, "p_y")}};
constexpr std::array<CopySpec, 5> ci_d{{
    COV_COPY("Ag", 1, "d_z2"), COV_COPY("Ag", 1, "d_x2-y2"),
    COV_COPY("Ag", 1, "d_xy"), COV_COPY("Ag", 1, "d_xz"),
    COV_COPY("Ag", 1, "d_yz")}};

// Cyclic/dihedral conventions use z as the principal axis.
constexpr std::array<CopySpec, 3> c2_p{{
    COV_COPY("A", 1, "p_z"), COV_COPY("B", 1, "p_x"),
    COV_COPY("B", 1, "p_y")}};
constexpr std::array<CopySpec, 5> c2_d{{
    COV_COPY("A", 1, "d_z2"), COV_COPY("A", 1, "d_x2-y2"),
    COV_COPY("A", 1, "d_xy"), COV_COPY("B", 1, "d_xz"),
    COV_COPY("B", 1, "d_yz")}};

constexpr std::array<CopySpec, 3> c2h_p{{
    COV_COPY("Au", 1, "p_z"), COV_COPY("Bu", 1, "p_x"),
    COV_COPY("Bu", 1, "p_y")}};
constexpr std::array<CopySpec, 5> c2h_d{{
    COV_COPY("Ag", 1, "d_z2"), COV_COPY("Ag", 1, "d_x2-y2"),
    COV_COPY("Ag", 1, "d_xy"), COV_COPY("Bg", 1, "d_xz"),
    COV_COPY("Bg", 1, "d_yz")}};

constexpr std::array<CopySpec, 3> c2v_p{{
    COV_COPY("A1", 1, "p_z"), COV_COPY("B1", 1, "p_x"),
    COV_COPY("B2", 1, "p_y")}};
constexpr std::array<CopySpec, 5> c2v_d{{
    COV_COPY("A1", 1, "d_z2"), COV_COPY("A1", 1, "d_x2-y2"),
    COV_COPY("A2", 1, "d_xy"), COV_COPY("B1", 1, "d_xz"),
    COV_COPY("B2", 1, "d_yz")}};

constexpr std::array<CopySpec, 3> d2_p{{
    COV_COPY("B1", 1, "p_z"), COV_COPY("B2", 1, "p_y"),
    COV_COPY("B3", 1, "p_x")}};
constexpr std::array<CopySpec, 5> d2_d{{
    COV_COPY("A", 1, "d_z2"), COV_COPY("A", 1, "d_x2-y2"),
    COV_COPY("B1", 1, "d_xy"), COV_COPY("B2", 1, "d_xz"),
    COV_COPY("B3", 1, "d_yz")}};

constexpr std::array<CopySpec, 3> d2h_p{{
    COV_COPY("B1u", 1, "p_z"), COV_COPY("B2u", 1, "p_y"),
    COV_COPY("B3u", 1, "p_x")}};
constexpr std::array<CopySpec, 5> d2h_d{{
    COV_COPY("Ag", 1, "d_z2"), COV_COPY("Ag", 1, "d_x2-y2"),
    COV_COPY("B1g", 1, "d_xy"), COV_COPY("B2g", 1, "d_xz"),
    COV_COPY("B3g", 1, "d_yz")}};

constexpr std::array<CopySpec, 2> c3v_p{{
    COV_COPY("A1", 1, "p_z"), COV_COPY("E", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> c3v_d{{
    COV_COPY("A1", 1, "d_z2"), COV_COPY("E", 2, "d_xz,d_yz"),
    COV_COPY("E", 2, "d_x2-y2,d_xy")}};

constexpr std::array<CopySpec, 2> c4v_p{{
    COV_COPY("A1", 1, "p_z"), COV_COPY("E", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 4> c4v_d{{
    COV_COPY("A1", 1, "d_z2"), COV_COPY("B1", 1, "d_x2-y2"),
    COV_COPY("B2", 1, "d_xy"), COV_COPY("E", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d2d_p{{
    COV_COPY("B2", 1, "p_z"), COV_COPY("E", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 4> d2d_d{{
    COV_COPY("A1", 1, "d_z2"), COV_COPY("B1", 1, "d_x2-y2"),
    COV_COPY("B2", 1, "d_xy"), COV_COPY("E", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d3d_p{{
    COV_COPY("A2u", 1, "p_z"), COV_COPY("Eu", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> d3d_d{{
    COV_COPY("A1g", 1, "d_z2"), COV_COPY("Eg", 2, "d_xz,d_yz"),
    COV_COPY("Eg", 2, "d_x2-y2,d_xy")}};

constexpr std::array<CopySpec, 2> d3h_p{{
    COV_COPY("A2''", 1, "p_z"), COV_COPY("E'", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> d3h_d{{
    COV_COPY("A1'", 1, "d_z2"), COV_COPY("E'", 2, "d_x2-y2,d_xy"),
    COV_COPY("E''", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d4h_p{{
    COV_COPY("A2u", 1, "p_z"), COV_COPY("Eu", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 4> d4h_d{{
    COV_COPY("A1g", 1, "d_z2"), COV_COPY("B1g", 1, "d_x2-y2"),
    COV_COPY("B2g", 1, "d_xy"), COV_COPY("Eg", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d4d_p{{
    COV_COPY("B2", 1, "p_z"), COV_COPY("E1", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> d4d_d{{
    COV_COPY("A1", 1, "d_z2"), COV_COPY("E2", 2, "d_x2-y2,d_xy"),
    COV_COPY("E3", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d5h_p{{
    COV_COPY("A2''", 1, "p_z"), COV_COPY("E1'", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> d5h_d{{
    COV_COPY("A1'", 1, "d_z2"), COV_COPY("E2'", 2, "d_x2-y2,d_xy"),
    COV_COPY("E1''", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> d5d_p{{
    COV_COPY("A2u", 1, "p_z"), COV_COPY("E1u", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> d5d_d{{
    COV_COPY("A1g", 1, "d_z2"), COV_COPY("E2g", 2, "d_x2-y2,d_xy"),
    COV_COPY("E1g", 2, "d_xz,d_yz")}};

constexpr std::array<CopySpec, 1> td_p{{COV_COPY("T2", 3, "p_x,p_y,p_z")}};
constexpr std::array<CopySpec, 2> td_d{{
    COV_COPY("E", 2, "d_z2,d_x2-y2"),
    COV_COPY("T2", 3, "d_xy,d_xz,d_yz")}};

constexpr std::array<CopySpec, 1> oh_p{{COV_COPY("T1u", 3, "p_x,p_y,p_z")}};
constexpr std::array<CopySpec, 2> oh_d{{
    COV_COPY("Eg", 2, "d_z2,d_x2-y2"),
    COV_COPY("T2g", 3, "d_xy,d_xz,d_yz")}};

constexpr std::array<CopySpec, 2> dinfh_p{{
    COV_COPY("Sigma_u+", 1, "p_z"), COV_COPY("Pi_u", 2, "p_x,p_y")}};
constexpr std::array<CopySpec, 3> dinfh_d{{
    COV_COPY("Sigma_g+", 1, "d_z2"), COV_COPY("Pi_g", 2, "d_xz,d_yz"),
    COV_COPY("Delta_g", 2, "d_x2-y2,d_xy")}};

const std::array<SPDRecord, 20> spd_records{{
    {"C1", as_span(s_a), as_span(c1_p), as_span(c1_d)},
    {"Cs", as_span(s_ap), as_span(cs_p), as_span(cs_d)},
    {"Ci", as_span(s_ag), as_span(ci_p), as_span(ci_d)},
    {"C2", as_span(s_a), as_span(c2_p), as_span(c2_d)},
    {"C2h", as_span(s_ag), as_span(c2h_p), as_span(c2h_d)},
    {"C2v", as_span(s_a1), as_span(c2v_p), as_span(c2v_d)},
    {"C3v", as_span(s_a1), as_span(c3v_p), as_span(c3v_d)},
    {"C4v", as_span(s_a1), as_span(c4v_p), as_span(c4v_d)},
    {"D2", as_span(s_a), as_span(d2_p), as_span(d2_d)},
    {"D2h", as_span(s_ag), as_span(d2h_p), as_span(d2h_d)},
    {"D2d", as_span(s_a1), as_span(d2d_p), as_span(d2d_d)},
    {"D3d", as_span(s_a1g), as_span(d3d_p), as_span(d3d_d)},
    {"D3h", as_span(s_a1p), as_span(d3h_p), as_span(d3h_d)},
    {"D4h", as_span(s_a1g), as_span(d4h_p), as_span(d4h_d)},
    {"D4d", as_span(s_a1), as_span(d4d_p), as_span(d4d_d)},
    {"D5h", as_span(s_a1p), as_span(d5h_p), as_span(d5h_d)},
    {"D5d", as_span(s_a1g), as_span(d5d_p), as_span(d5d_d)},
    {"Td", as_span(s_a1), as_span(td_p), as_span(td_d)},
    {"Oh", as_span(s_a1g), as_span(oh_p), as_span(oh_d)},
    {"Dinfh", as_span(s_sigma_g), as_span(dinfh_p), as_span(dinfh_d)},
}};

#undef COV_COPY

std::string normalise_group_name(std::string_view name) {
    std::string result;
    result.reserve(name.size());
    for (const unsigned char ch : name) {
        if (std::isspace(ch) || ch == '_' || ch == '{' || ch == '}' ||
            ch == '^' || ch == '-') {
            continue;
        }
        result.push_back(static_cast<char>(std::tolower(ch)));
    }
    return result;
}

const SPDRecord* find_spd_record(std::string_view canonical_name) noexcept {
    const auto found = std::find_if(
        spd_records.begin(), spd_records.end(),
        [canonical_name](const SPDRecord& record) {
            return record.group == canonical_name;
        });
    return found == spd_records.end() ? nullptr : &*found;
}

IrrepDecomposition materialise(std::span<const CopySpec> specs, Shell shell) {
    IrrepDecomposition result;
    result.reserve(specs.size());
    for (const auto& spec : specs) {
        std::uint8_t copy_index = 1;
        for (const auto& previous : result) {
            if (previous.label == spec.label) ++copy_index;
        }
        result.push_back({spec.label, spec.dimension, shell, copy_index, spec.basis});
    }
    return result;
}

} // namespace

std::size_t MetalSPDDecomposition::total_dimension() const noexcept {
    return irrep_total_dimension(s) + irrep_total_dimension(p) +
           irrep_total_dimension(d);
}

std::span<const PointGroupDefinition> point_group_catalog() noexcept {
    return catalog;
}

const PointGroupDefinition* find_point_group(std::string_view name) noexcept {
    const std::string wanted = normalise_group_name(name);
    for (const auto& group : catalog) {
        if (normalise_group_name(group.canonical_name) == wanted) return &group;
        for (const auto alias : group.aliases) {
            if (normalise_group_name(alias) == wanted) return &group;
        }
    }
    return nullptr;
}

std::optional<IrrepDecomposition> decompose_metal_ao_shell(
    std::string_view point_group,
    MetalAOShell shell) {
    const auto* group = find_point_group(point_group);
    if (!group) return std::nullopt;
    const auto* record = find_spd_record(group->canonical_name);
    if (!record) return std::nullopt;
    switch (shell) {
        case Shell::S: return materialise(record->s, shell);
        case Shell::P: return materialise(record->p, shell);
        case Shell::D: return materialise(record->d, shell);
    }
    return std::nullopt;
}

std::optional<MetalSPDDecomposition> decompose_metal_spd(
    std::string_view point_group) {
    const auto s = decompose_metal_ao_shell(point_group, Shell::S);
    const auto p = decompose_metal_ao_shell(point_group, Shell::P);
    const auto d = decompose_metal_ao_shell(point_group, Shell::D);
    if (!s || !p || !d) return std::nullopt;
    return MetalSPDDecomposition{*s, *p, *d};
}

std::size_t irrep_multiplicity(
    std::span<const IrrepCopy> decomposition,
    std::string_view label) noexcept {
    return static_cast<std::size_t>(std::count_if(
        decomposition.begin(), decomposition.end(),
        [label](const IrrepCopy& copy) { return copy.label == label; }));
}

std::size_t irrep_total_dimension(
    std::span<const IrrepCopy> decomposition) noexcept {
    std::size_t result = 0;
    for (const auto& copy : decomposition) result += copy.dimension;
    return result;
}

} // namespace cov
