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

               const std::uint8_t b, const std::uint8_t a = 255) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        const std::size_t i =
            (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
             static_cast<std::size_t>(x)) * 4u;
        rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = a;
    }

    void line(int x0, int y0, const int x1, const int y1,
              const std::uint8_t r, const std::uint8_t g,
              const std::uint8_t b, const int thickness = 1) {
        const int dx = std::abs(x1 - x0);
        const int sx = x0 < x1 ? 1 : -1;
        const int dy = -std::abs(y1 - y0);
        const int sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        for (;;) {
            const int radius = std::max(0, thickness / 2);
            for (int oy = -radius; oy <= radius; ++oy) {
                for (int ox = -radius; ox <= radius; ++ox) {
                    pixel(x0 + ox, y0 + oy, r, g, b);
                }
            }
            if (x0 == x1 && y0 == y1) break;
            const int e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }
};

std::array<std::uint8_t, 7> glyph_rows(char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    switch (c) {
        case '0': return {14,17,19,21,25,17,14};
        case '1': return {4,12,4,4,4,4,14};
        case '2': return {14,17,1,2,4,8,31};
        case '3': return {30,1,1,14,1,1,30};
        case '4': return {2,6,10,18,31,2,2};
        case '5': return {31,16,16,30,1,1,30};
        case '6': return {14,16,16,30,17,17,14};
        case '7': return {31,1,2,4,8,8,8};
        case '8': return {14,17,17,14,17,17,14};
        case '9': return {14,17,17,15,1,1,14};
        case 'A': return {14,17,17,31,17,17,17};
        case 'B': return {30,17,17,30,17,17,30};
        case 'C': return {14,17,16,16,16,17,14};
        case 'D': return {30,17,17,17,17,17,30};
        case 'E': return {31,16,16,30,16,16,31};
        case 'F': return {31,16,16,30,16,16,16};
        case 'G': return {14,17,16,23,17,17,15};
        case 'H': return {17,17,17,31,17,17,17};
        case 'I': return {14,4,4,4,4,4,14};
        case 'J': return {7,2,2,2,2,18,12};
        case 'K': return {17,18,20,24,20,18,17};
        case 'L': return {16,16,16,16,16,16,31};
        case 'M': return {17,27,21,21,17,17,17};
        case 'N': return {17,25,21,19,17,17,17};
        case 'O': return {14,17,17,17,17,17,14};
        case 'P': return {30,17,17,30,16,16,16};
        case 'Q': return {14,17,17,17,21,18,13};
        case 'R': return {30,17,17,30,20,18,17};
        case 'S': return {15,16,16,14,1,1,30};
        case 'T': return {31,4,4,4,4,4,4};
        case 'U': return {17,17,17,17,17,17,14};
        case 'V': return {17,17,17,17,17,10,4};
        case 'W': return {17,17,17,21,21,21,10};
        case 'X': return {17,17,10,4,10,17,17};
        case 'Y': return {17,17,10,4,4,4,4};
        case 'Z': return {31,1,2,4,8,16,31};
        case '-': return {0,0,0,31,0,0,0};
        case '.': return {0,0,0,0,0,12,12};
        case '/': return {1,2,2,4,8,8,16};
        case '+': return {0,4,4,31,4,4,0};
        case ':': return {0,12,12,0,12,12,0};
        case '?': return {14,17,1,2,4,0,4};
        case '_': return {0,0,0,0,0,0,31};
        case '*': return {0,21,14,31,14,21,0};
        case '=': return {0,31,0,31,0,0,0};
        default: return {0,0,0,0,0,0,0};
    }
}

void draw_text(Canvas& canvas, int x, const int y, const std::string& text,
               const std::uint8_t r = 32, const std::uint8_t g = 41,
               const std::uint8_t b = 55, const int scale = 2) {
    for (const char c : text) {
        if (c == ' ') { x += 4 * scale; continue; }
        const auto rows = glyph_rows(c);
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if ((rows[static_cast<std::size_t>(row)] & (1u << (4 - col))) == 0) continue;
                for (int yy = 0; yy < scale; ++yy) {
                    for (int xx = 0; xx < scale; ++xx) {
                        canvas.pixel(x + col * scale + xx,
                                     y + row * scale + yy, r, g, b);
                    }
                }
            }
        }
        x += 6 * scale;
    }
}

void draw_electron(Canvas& canvas, const int x, const int level_y,
                   const bool up, const std::uint8_t r = 28,
                   const std::uint8_t g = 38, const std::uint8_t b = 52) {
    // Textbook convention: both spin arrows occupy the space above the energy
    // level. The beta arrow points down *towards* the level instead of below it,
    // leaving the orbital-type and symmetry labels underneath completely clear.
    constexpr int kArrowLength = 11;
    const int shaft_top = level_y - kArrowLength;
    const int tip = up ? shaft_top : level_y;
    const int tail = up ? level_y : shaft_top;
    canvas.line(x, tail, x, tip, r, g, b, 2);
    const int dy = up ? 4 : -4;
    canvas.line(x, tip, x - 4, tip + dy, r, g, b, 2);
    canvas.line(x, tip, x + 4, tip + dy, r, g, b, 2);
}

std::uint32_t crc32(const std::uint8_t* data, const std::size_t size) {
    std::uint32_t crc = 0xffffffffu;
    for (std::size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) {
            const std::uint32_t mask = 0u - (crc & 1u);
            crc = (crc >> 1u) ^ (0xedb88320u & mask);
        }
    }
    return ~crc;
}

std::uint32_t adler32(const std::uint8_t* data, const std::size_t size) {
    constexpr std::uint32_t mod = 65521u;
    std::uint32_t a = 1u, b = 0u;
    for (std::size_t i = 0; i < size; ++i) {
        a = (a + data[i]) % mod;
        b = (b + a) % mod;
    }
    return (b << 16u) | a;
}

void push_u32_be(std::vector<std::uint8_t>& out, const std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>((value >> 24u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 16u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 8u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>(value & 0xffu));
}

void append_png_chunk(std::vector<std::uint8_t>& png,
                      const char type[4],
                      const std::vector<std::uint8_t>& payload) {
    push_u32_be(png, static_cast<std::uint32_t>(payload.size()));
    const std::size_t crc_begin = png.size();
    png.insert(png.end(), type, type + 4);
    png.insert(png.end(), payload.begin(), payload.end());
    push_u32_be(png, crc32(png.data() + crc_begin, 4 + payload.size()));
}

bool write_png_rgba(const Canvas& canvas,
                    const std::filesystem::path& path,
                    std::string* error,
                    const std::string& view_id = {}) {
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(canvas.height) *
                (1u + static_cast<std::size_t>(canvas.width) * 4u));
    for (int y = 0; y < canvas.height; ++y) {
        raw.push_back(0);
        const std::size_t offset = static_cast<std::size_t>(y) *
                                   static_cast<std::size_t>(canvas.width) * 4u;
        raw.insert(raw.end(),
                   canvas.rgba.begin() + static_cast<std::ptrdiff_t>(offset),
                   canvas.rgba.begin() + static_cast<std::ptrdiff_t>(
                       offset + static_cast<std::size_t>(canvas.width) * 4u));
    }

    std::vector<std::uint8_t> zlib{0x78, 0x01};
    std::size_t pos = 0;
    while (pos < raw.size()) {
        const std::size_t block = std::min<std::size_t>(65535u, raw.size() - pos);
        const bool final = pos + block == raw.size();
        zlib.push_back(final ? 0x01 : 0x00);
        const std::uint16_t len = static_cast<std::uint16_t>(block);
        const std::uint16_t nlen = static_cast<std::uint16_t>(~len);
        zlib.push_back(static_cast<std::uint8_t>(len & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>((len >> 8u) & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>(nlen & 0xffu));
