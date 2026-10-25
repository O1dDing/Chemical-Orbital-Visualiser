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

#include "cov/version.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace cov::export_attribution {

namespace detail {
inline std::string xml_escape(const std::string_view text) {
    std::string result;
    for (const char ch : text) {
        switch (ch) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '\"': result += "&quot;"; break;
        case '\'': result += "&apos;"; break;
        default: result += ch;
        }
    }
    return result;
}

inline void put_u32(std::vector<std::uint8_t>& out, const std::uint32_t value) {
    for (int shift = 24; shift >= 0; shift -= 8)
        out.push_back(static_cast<std::uint8_t>(value >> shift));
}

// PNG iTXt stores uncompressed UTF-8 text. The fixed ASCII keywords are valid
// PNG keywords; the author fields explicitly describe the generating software.
inline void append_itxt(std::vector<std::uint8_t>& out,
                        const std::string_view keyword, const std::string_view value) {
    const auto payload_size = keyword.size() + 5u + value.size();
    put_u32(out, static_cast<std::uint32_t>(payload_size));
    const auto begin = out.size();
    for (const char ch : std::string_view("iTXt")) out.push_back(static_cast<std::uint8_t>(ch));
    for (const char ch : keyword) out.push_back(static_cast<std::uint8_t>(ch));
    // Keyword terminator, compression flag/method, empty language tag and
    // empty translated keyword, followed by the UTF-8 text.
    out.insert(out.end(), 5u, 0u);
    for (const char ch : value) out.push_back(static_cast<std::uint8_t>(ch));
    std::uint32_t crc = 0xffffffffu;
    for (auto i = begin; i < out.size(); ++i) {
        crc ^= out[i];
        for (int bit = 0; bit < 8; ++bit)
            crc = (crc >> 1u) ^ ((crc & 1u) ? 0xedb88320u : 0u);
    }
    put_u32(out, ~crc);
}
} // namespace detail

inline std::string software_name() {
    return std::string(cov::display_name) + " " + std::string(cov::version);
}

// Call after IDAT and before IEND. This only appends ancillary chunks and does
// not touch the pixel stream, its compression, or any existing PNG chunks.
inline void append_png_text_chunks(std::vector<std::uint8_t>& png) {
    detail::append_itxt(png, "Software", software_name());
    detail::append_itxt(png, "Software Version", cov::version);
    detail::append_itxt(png, "Software Author", cov::software_author);
    detail::append_itxt(png, "Software URL", cov::project_url);
}

// XMP identifies the generating tool; DOAP describes that software project and
// its developer/release. No artwork creator or artwork copyright is asserted.
inline std::string svg_metadata() {
    const auto name = detail::xml_escape(cov::display_name);
    const auto version = detail::xml_escape(cov::version);
    const auto author = detail::xml_escape(cov::software_author);
    const auto project = detail::xml_escape(cov::project_url);
    return "<metadata><rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" "
        "xmlns:xmp=\"http://ns.adobe.com/xap/1.0/\" "
        "xmlns:doap=\"http://usefulinc.com/ns/doap#\" "
        "xmlns:foaf=\"http://xmlns.com/foaf/0.1/\">"
        "<rdf:Description rdf:about=\"\"><xmp:CreatorTool>" +
        detail::xml_escape(software_name()) + "</xmp:CreatorTool></rdf:Description>"
        "<doap:Project rdf:about=\"" + project + "\"><doap:name>" + name +
        "</doap:name><doap:homepage rdf:resource=\"" + project + "\"/>"
        "<doap:developer><foaf:Person><foaf:name>" + author +
        "</foaf:name></foaf:Person></doap:developer>"
        "<doap:release><doap:Version><doap:revision>" + version +
        "</doap:revision></doap:Version></doap:release></doap:Project>"
        "</rdf:RDF></metadata>";
}

} // namespace cov::export_attribution
