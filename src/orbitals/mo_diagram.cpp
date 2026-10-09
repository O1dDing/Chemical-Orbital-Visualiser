#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string_view>
#include <vector>

namespace cov {

const char* energy_axis_mode_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "linear" : "nonlinear-focus";
}

const char* energy_transform_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "identity" : "asinh-focus-v1";
}

double energy_display_coordinate(const double energy_hartree,
                                 const EnergyTransform& transform) noexcept {
    if (transform.mode == EnergyAxisMode::Linear) return energy_hartree;
    const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
    return std::asinh((energy_hartree - transform.focus_hartree) / scale);
}

double energy_from_display_coordinate(const double coordinate,
                                      const EnergyTransform& transform) noexcept {
    if (transform.mode == EnergyAxisMode::Linear) return coordinate;
    const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
    return transform.focus_hartree + scale * std::sinh(coordinate);
}
namespace {

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

const char* family_symbol(const std::string& family) noexcept {
    if (family == "sigma") return "σ";
    if (family == "pi") return "π";
    if (family == "delta") return "δ";
    if (family == "phi") return "φ";
    return "";
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

struct PositionedLevel {
    const MODiagramLevel* level = nullptr;
    double x = 0.0;
    double y = 0.0;
    double label_y = 0.0;
};

struct EnergyRange {
    double min = -1.0;
    double max = 1.0;
};

struct EnergyScale {
    EnergyRange range;
    EnergyTransform transform;
    double warped_min = 0.0;
    double warped_max = 1.0;
};

double warp_energy(const double energy, const EnergyScale& scale) noexcept {
    return energy_display_coordinate(energy, scale.transform);
}

EnergyScale energy_scale(const MODiagramData& data,
                         const EnergyRange& range,
                         const MODiagramOptions& options) {
    EnergyScale scale;
    scale.range = range;
    scale.transform.mode = options.energy_axis_mode;
    scale.transform.focus_hartree = (range.min + range.max) * 0.5;
    for (const auto& level : data.levels) {
        if (level.metadata.selected) {
            scale.transform.focus_hartree = level.metadata.energy_hartree;
            break;
        }
    }
    scale.transform.scale_hartree = std::max(1.0e-5, (range.max - range.min) * 0.04);
    scale.warped_min = warp_energy(range.min, scale);
    scale.warped_max = warp_energy(range.max, scale);
    return scale;
}

double energy_at_fraction(const double fraction, const EnergyScale& scale) noexcept {
    const double warped = scale.warped_min + fraction * (scale.warped_max - scale.warped_min);
    return energy_from_display_coordinate(warped, scale.transform);
}

EnergyRange energy_range(const MODiagramData& data) {
    EnergyRange range;
    if (data.levels.empty()) return range;
    range.min = data.levels.front().metadata.energy_hartree;
    range.max = range.min;
    for (const auto& level : data.levels) {
        range.min = std::min(range.min, level.metadata.energy_hartree);
        range.max = std::max(range.max, level.metadata.energy_hartree);
    }
    if (std::abs(range.max - range.min) < 1.0e-10) {
        range.min -= 0.5;
        range.max += 0.5;
    } else {
        const double pad = std::max(0.02, (range.max - range.min) * 0.09);
        range.min -= pad;
        range.max += pad;
    }
    return range;
}

double map_energy_y(const double energy,
                    const EnergyScale& scale,
                    const int height,
                    const int top = 105,
                    const int bottom = 70) {
    const double t = (warp_energy(energy, scale) - scale.warped_min) /
                     std::max(1.0e-12, scale.warped_max - scale.warped_min);
    return static_cast<double>(height - bottom) -
           t * static_cast<double>(height - top - bottom);
}

int group_member_index(const std::string& label) noexcept {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    if (suffix < 'a' || suffix > 'z') return 0;
    return static_cast<int>(suffix - 'a');
}

std::vector<PositionedLevel> position_levels(const MODiagramData& data,
                                             const MODiagramOptions& options,
                                             const int width,
                                             const int height) {
    std::vector<PositionedLevel> positioned;
    positioned.reserve(data.levels.size());
    const EnergyRange range = energy_range(data);
    const EnergyScale scale = energy_scale(data, range, options);
    const double centre = static_cast<double>(width) * 0.50;

    for (const auto& level : data.levels) {
        const std::size_t n = std::max<std::size_t>(1, level.metadata.degeneracy_size);
        const int member = group_member_index(level.metadata.display_label);
        const double spacing = std::clamp(static_cast<double>(width) * 0.105, 92.0, 138.0);
        const double offset = n <= 1
                                  ? 0.0
                                  : (static_cast<double>(member) -
                                     (static_cast<double>(n) - 1.0) * 0.5) * spacing;
        const double y = map_energy_y(level.metadata.energy_hartree, scale, height);
        positioned.push_back({&level, centre + offset, y, y - 7.0});
    }

    // Keep physical level y quantitative. Labels share a dedicated outer
    // column and are spaced globally, so labels can never collide with a
    // neighbouring degeneracy lane or masquerade as shifted energy levels.
    std::vector<std::size_t> ids(positioned.size());
    std::iota(ids.begin(), ids.end(), 0);
    std::sort(ids.begin(), ids.end(), [&](const std::size_t a, const std::size_t b) {
        if (positioned[a].label_y != positioned[b].label_y) {
            return positioned[a].label_y < positioned[b].label_y;
        }
        return positioned[a].x < positioned[b].x;
    });
    constexpr double kLabelSpacing = 32.0;
    constexpr double kTopLabel = 104.0;
    const double bottom_label = static_cast<double>(height - 48);
    double previous = kTopLabel - kLabelSpacing;
    for (const std::size_t id : ids) {
        positioned[id].label_y = std::max(positioned[id].label_y, previous + kLabelSpacing);
        previous = positioned[id].label_y;
    }
    const double overflow = previous - bottom_label;
    if (overflow > 0.0) {
        for (const std::size_t id : ids) positioned[id].label_y -= overflow;
    }
    double next = bottom_label + kLabelSpacing;
    for (auto it = ids.rbegin(); it != ids.rend(); ++it) {
        positioned[*it].label_y = std::min(positioned[*it].label_y, next - kLabelSpacing);
        next = positioned[*it].label_y;
    }
    return positioned;
}

struct Canvas {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgba;

    Canvas(const int w, const int h)
        : width(std::max(1, w)), height(std::max(1, h)),
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
               const std::uint8_t b, const std::uint8_t a = 255) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        const std::size_t i = (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
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
                for (int ox = -radius; ox <= radius; ++ox) pixel(x0 + ox, y0 + oy, r, g, b);
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

void draw_electron(Canvas& canvas, const int x, const int y,
                   const bool up, const std::uint8_t r = 28,
                   const std::uint8_t g = 38, const std::uint8_t b = 52) {
    const int tip = up ? y - 17 : y + 17;
    canvas.line(x, y, x, tip, r, g, b, 2);
    const int dy = up ? 5 : -5;
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
                    std::string* error) {
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(canvas.height) *
                (1u + static_cast<std::size_t>(canvas.width) * 4u));
    for (int y = 0; y < canvas.height; ++y) {
        raw.push_back(0);
        const std::size_t offset = static_cast<std::size_t>(y) *
                                   static_cast<std::size_t>(canvas.width) * 4u;
        raw.insert(raw.end(), canvas.rgba.begin() + static_cast<std::ptrdiff_t>(offset),
                   canvas.rgba.begin() + static_cast<std::ptrdiff_t>(offset +
                       static_cast<std::size_t>(canvas.width) * 4u));
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
        zlib.push_back(static_cast<std::uint8_t>((nlen >> 8u) & 0xffu));
        zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(pos),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos + block));
        pos += block;
    }
    push_u32_be(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> png = {137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.width));
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.height));
    ihdr.push_back(8); ihdr.push_back(6); ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    append_png_chunk(png, "IHDR", ihdr);
    append_png_chunk(png, "IDAT", zlib);
    append_png_chunk(png, "IEND", {});

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        if (error) *error = "Unable to open PNG output path";
        return false;
    }
    file.write(reinterpret_cast<const char*>(png.data()), static_cast<std::streamsize>(png.size()));
    if (!file) {
        if (error) *error = "Failed while writing PNG";
        return false;
    }
    return true;
}

std::string annotation_summary(const OrbitalAnnotation& annotation) {
    std::string out;
    if (annotation.family != "unavailable") out += annotation.family;
    if (annotation.bonding_class != BondingClass::Unclassified) {
        if (!out.empty()) out += " / ";
        out += bonding_class_name(annotation.bonding_class);
    }
    return out;
}

} // namespace

const char* annotation_source_name(const AnnotationSource source) noexcept {
    switch (source) {
        case AnnotationSource::Direct: return "direct";
        case AnnotationSource::ParsedLabel: return "parsed-label";
        case AnnotationSource::Heuristic: return "heuristic";
        default: return "unavailable";
    }
}

const char* bonding_class_name(const BondingClass value) noexcept {
    switch (value) {
        case BondingClass::Bonding: return "bonding";
        case BondingClass::Nonbonding: return "nonbonding";
        case BondingClass::Antibonding: return "antibonding";
        case BondingClass::Mixed: return "mixed";
        default: return "unclassified";
    }
}

OrbitalAnnotation annotate_orbital(const MolecularOrbital& orbital) {
    OrbitalAnnotation result;
    const std::string label = lower_ascii(orbital.symmetry);

    if (contains_any(label, {"sigma", "σ"})) {
        result.family = "sigma";
        result.family_source = AnnotationSource::ParsedLabel;
        result.family_confidence = 1.0;
    } else if (contains_any(label, {"delta", "δ"})) {
        result.family = "delta";
        result.family_source = AnnotationSource::ParsedLabel;
        result.family_confidence = 1.0;
    } else if (contains_any(label, {"phi", "φ"})) {
        result.family = "phi";
        result.family_source = AnnotationSource::ParsedLabel;
        result.family_confidence = 1.0;
    } else if (contains_any(label, {"pi", "π"})) {
        result.family = "pi";
        result.family_source = AnnotationSource::ParsedLabel;
        result.family_confidence = 1.0;
    }

    if (contains_any(label, {"antibond", "anti-bond", "sigma*", "pi*", "delta*"})) {
        result.bonding_class = BondingClass::Antibonding;
        result.bonding_source = AnnotationSource::ParsedLabel;
        result.bonding_confidence = 0.95;
    } else if (contains_any(label, {"nonbond", "non-bond", " n.b.", " nb"})) {
        result.bonding_class = BondingClass::Nonbonding;
        result.bonding_source = AnnotationSource::ParsedLabel;
        result.bonding_confidence = 0.95;
    } else if (contains_any(label, {"bonding", "bonding-mo"})) {
        result.bonding_class = BondingClass::Bonding;
        result.bonding_source = AnnotationSource::ParsedLabel;
        result.bonding_confidence = 0.95;
    }

    // Deliberately do not infer bonding/nonbonding/antibonding from energy or
    // occupancy alone: standard Molden MO data are insufficient for that.
    result.heuristic = result.family_source == AnnotationSource::Heuristic ||
                       result.bonding_source == AnnotationSource::Heuristic;
    return result;
}

DiagramSelectionPlan build_valence_selection_plan(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata) {
    DiagramSelectionPlan plan;
    if (wavefunction.orbitals.empty() || metadata.empty()) {
        plan.summary = "valence-auto: no orbitals";
        return plan;
    }

    const FrontierOrbitals frontier = find_frontier_orbitals(
        wavefunction.orbitals, options.filter.occupation_threshold);
    const std::size_t max_levels = effective_max_levels(options);
    const std::size_t max_virtual = std::min(options.max_virtual_levels,
                                             std::max<std::size_t>(2, max_levels / 2));

    std::vector<std::size_t> occupied_valence;
    std::vector<std::size_t> frontier_virtual;
    for (std::size_t i = 0; i < wavefunction.orbitals.size(); ++i) {
        const auto& mo = wavefunction.orbitals[i];
        if (occupied(mo, options.filter.occupation_threshold)) {
            if (metadata[i].region == OrbitalRegion::Valence) occupied_valence.push_back(i);
            continue;
        }
        if (metadata[i].region != OrbitalRegion::Virtual) continue;
        if (frontier.lumo && *frontier.lumo < wavefunction.orbitals.size()) {
            const double cutoff = wavefunction.orbitals[*frontier.lumo].energy_hartree +
                                  options.filter.virtual_window_hartree;
            if (mo.energy_hartree <= cutoff) frontier_virtual.push_back(i);
        } else {
            frontier_virtual.push_back(i);
        }
    }

    const std::size_t virtual_keep = std::min(max_virtual, frontier_virtual.size());
    const std::size_t occupied_budget = max_levels > virtual_keep ? max_levels - virtual_keep : 1;
    const std::size_t occ_begin = occupied_valence.size() > occupied_budget
                                      ? occupied_valence.size() - occupied_budget
                                      : 0;

    std::vector<std::size_t> chosen;
    chosen.insert(chosen.end(), occupied_valence.begin() + static_cast<std::ptrdiff_t>(occ_begin),
                  occupied_valence.end());
    chosen.insert(chosen.end(), frontier_virtual.begin(),
                  frontier_virtual.begin() + static_cast<std::ptrdiff_t>(virtual_keep));

    if (options.selected_index < wavefunction.orbitals.size() &&
        metadata[options.selected_index].region != OrbitalRegion::Core) {
        chosen.push_back(options.selected_index);
    }

    std::sort(chosen.begin(), chosen.end());
    chosen.erase(std::unique(chosen.begin(), chosen.end()), chosen.end());
    const auto labels = build_orbital_labels(wavefunction.orbitals, options.degeneracy);
    chosen = expand_degenerate_sets(chosen, labels);

    plan.included_indices = std::move(chosen);
    plan.hidden_count = metadata.size() > plan.included_indices.size()
                            ? metadata.size() - plan.included_indices.size()
                            : 0;
    plan.valence_occupied_count = 0;
    plan.frontier_virtual_count = 0;
    for (const std::size_t index : plan.included_indices) {
        if (index >= metadata.size()) continue;
        if (metadata[index].region == OrbitalRegion::Valence &&
            occupied(wavefunction.orbitals[index], options.filter.occupation_threshold)) {
            ++plan.valence_occupied_count;
        } else if (metadata[index].region == OrbitalRegion::Virtual) {
            ++plan.frontier_virtual_count;
        }
    }

    std::ostringstream summary;
    summary << "valence-auto: " << plan.included_indices.size() << '/' << metadata.size()
            << " levels; occupied-valence=" << plan.valence_occupied_count
            << "; frontier-virtual=" << plan.frontier_virtual_count
            << "; core/high-virtual hidden non-destructively";
    plan.summary = summary.str();
    return plan;
}

MODiagramData build_mo_diagram_data(const Wavefunction& wavefunction,
                                    const MODiagramOptions& options) {
    MODiagramData data;
    data.mode = MODiagramMode::ValenceCentral;
    data.plan.classification = DiagramClassification::Simple;
    data.plan.complex_system = wavefunction.atoms.size() > 12 || wavefunction.orbitals.size() > 80;
    data.plan.strict_salc_available = false;
    data.plan.machine_reason =
        "valence-central: direct Molden MO data; strict SALC reconstruction is not claimed";
    data.frontier = find_frontier_orbitals(wavefunction.orbitals,
                                           options.filter.occupation_threshold);
    data.metadata = build_orbital_metadata(wavefunction,
                                           options.selected_index,
                                           options.degeneracy,
                                           options.filter);
    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& mo : wavefunction.orbitals) data.annotations.push_back(annotate_orbital(mo));

    data.selection = build_valence_selection_plan(wavefunction, options, data.metadata);
    data.levels.reserve(data.selection.included_indices.size());
    for (const std::size_t index : data.selection.included_indices) {
        if (index >= data.metadata.size()) continue;
        MODiagramLevel level;
        level.metadata = data.metadata[index];
        level.annotation = data.annotations[index];
        level.electrons = electron_glyphs_for_orbital(
            wavefunction.orbitals[index], data.frontier.separate_spin_sets);
        level.homo = data.frontier.homo && *data.frontier.homo == index;
        level.lumo = data.frontier.lumo && *data.frontier.lumo == index;
        data.levels.push_back(std::move(level));
    }
    return data;
}

bool write_mo_diagram_svg(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open SVG output path";
        return false;
    }

    const int width = std::max(760, options.width);
    const int height = std::max(620, options.height);
    const EnergyRange range = energy_range(data);
    const EnergyScale scale = energy_scale(data, range, options);
    const auto positioned = position_levels(data, options, width, height);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << ' ' << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#fbfcfe\"/>\n";
    out << "<text x=\"36\" y=\"40\" fill=\"#182231\" font-family=\"sans-serif\" font-size=\"24\" font-weight=\"600\">Valence MO diagram</text>\n";
    out << "<text x=\"36\" y=\"64\" fill=\"#68778b\" font-family=\"sans-serif\" font-size=\"12\">"
        << xml_escape(data.selection.summary) << "</text>\n";
    out << "<text x=\"36\" y=\"82\" fill=\"#7f8c9c\" font-family=\"sans-serif\" font-size=\"11\">"
        << "Molden-derived MO energies/occupations/symmetry; orbital family/class is shown only when explicitly supported by labels." << "</text>\n";

    const int axis_x = 70;
    out << "<line x1=\"" << axis_x << "\" y1=\"105\" x2=\"" << axis_x
        << "\" y2=\"" << height - 62 << "\" stroke=\"#2d3745\" stroke-width=\"2\"/>\n";
    out << "<path d=\"M " << axis_x << " 96 l -5 10 h 10 z\" fill=\"#2d3745\"/>\n";
    out << "<text x=\"24\" y=\"120\" fill=\"#4d5a6c\" font-family=\"sans-serif\" font-size=\"12\" transform=\"rotate(-90 24 120)\">Energy · "
        << (options.energy_axis_mode == EnergyAxisMode::Linear ? "linear" : "nonlinear focus (asinh)")
        << " ("
        << xml_escape(energy_unit_symbol(options.energy_unit)) << ")</text>\n";

    for (const double fraction : {0.0, 0.5, 1.0}) {
        const double e = energy_at_fraction(fraction, scale);
        const double y = map_energy_y(e, scale, height);
        out << "<line x1=\"64\" y1=\"" << y << "\" x2=\"76\" y2=\"" << y
            << "\" stroke=\"#7c8795\" stroke-width=\"1\"/>\n";
        out << "<text x=\"80\" y=\"" << y + 4 << "\" fill=\"#7c8795\" font-family=\"sans-serif\" font-size=\"10\">"
            << xml_escape(format_energy(e, options.energy_unit,
                                        options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 3))
            << "</text>\n";
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        const char* colour = level.metadata.selected ? "#286de0" :
                             level.homo ? "#138a62" :
                             level.lumo ? "#b56a00" :
                             level.metadata.region == OrbitalRegion::Virtual ? "#778496" : "#29384c";
        const double x0 = item.x - 45.0;
        const double x1 = item.x + 45.0;
        out << "<line x1=\"" << x0 << "\" y1=\"" << item.y << "\" x2=\"" << x1
            << "\" y2=\"" << item.y << "\" stroke=\"" << colour << "\" stroke-width=\""
            << (level.metadata.selected ? 4 : 2.2) << "\"/>\n";

        if (level.electrons.alpha > 0) {
            out << "<text x=\"" << item.x - 14 << "\" y=\"" << item.y - 5
                << "\" fill=\"#151d28\" font-family=\"sans-serif\" font-size=\"19\">↑</text>\n";
        }
        if (level.electrons.beta > 0) {
            out << "<text x=\"" << item.x + 3 << "\" y=\"" << item.y - 5
                << "\" fill=\"#151d28\" font-family=\"sans-serif\" font-size=\"19\">↓</text>\n";
        }

        const double label_x = static_cast<double>(width) * 0.72;
        if (std::abs(item.label_y - item.y) > 5.0) {
            out << "<line x1=\"" << x1 + 2 << "\" y1=\"" << item.y << "\" x2=\""
                << label_x - 6 << "\" y2=\"" << item.label_y << "\" stroke=\"#b6c0cc\" stroke-width=\"1\"/>\n";
        }
        out << "<text x=\"" << label_x << "\" y=\"" << item.label_y
            << "\" fill=\"" << colour << "\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"600\">"
            << xml_escape(level.metadata.display_label) << "</text>\n";

        std::string secondary;
        if (!level.metadata.symmetry_view.label.empty()) secondary += level.metadata.symmetry_view.label;
        const char* symbol = family_symbol(level.annotation.family);
        if (*symbol) {
            if (!secondary.empty()) secondary += " · ";
            secondary += symbol;
        }
        if (level.annotation.bonding_class != BondingClass::Unclassified) {
            if (!secondary.empty()) secondary += " · ";
            secondary += bonding_class_name(level.annotation.bonding_class);
        }
        if (!secondary.empty()) {
            out << "<text x=\"" << label_x << "\" y=\"" << item.label_y + 14
                << "\" fill=\"#68778b\" font-family=\"sans-serif\" font-size=\"10\">"
                << xml_escape(secondary) << "</text>\n";
        }
        out << "<text x=\"" << label_x << "\" y=\"" << item.label_y + 27
            << "\" fill=\"#8793a3\" font-family=\"sans-serif\" font-size=\"9\">"
            << xml_escape(format_energy(level.metadata.energy_hartree, options.energy_unit,
                                        options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 4))
            << " · raw MO " << level.metadata.raw_mo_number << "</text>\n";
    }

    out << "</svg>\n";
    if (!out) {
        if (error) *error = "Failed while writing SVG";
        return false;
    }
    return true;
}

bool write_mo_diagram_png(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    const int width = std::max(760, options.width);
    const int height = std::max(620, options.height);
    Canvas canvas(width, height);
    const EnergyRange range = energy_range(data);
    const EnergyScale scale = energy_scale(data, range, options);
    const auto positioned = position_levels(data, options, width, height);

    draw_text(canvas, 28, 22, "VALENCE MO DIAGRAM", 24, 34, 49, 3);
    draw_text(canvas, 28, 50,
              std::to_string(data.levels.size()) + " OF " + std::to_string(data.metadata.size()) +
                  " ORBITALS SHOWN",
              101, 116, 136, 2);

    const int axis_x = 70;
    canvas.line(axis_x, 96, axis_x, height - 60, 45, 55, 70, 2);
    canvas.line(axis_x, 96, axis_x - 5, 106, 45, 55, 70, 2);
    canvas.line(axis_x, 96, axis_x + 5, 106, 45, 55, 70, 2);
    draw_text(canvas, 18, 96,
              options.energy_axis_mode == EnergyAxisMode::Linear ? "E LINEAR" : "E FOCUS ASINH",
              45, 55, 70, 2);

    for (const double fraction : {0.0, 0.5, 1.0}) {
        const double e = energy_at_fraction(fraction, scale);
        const int y = static_cast<int>(std::lround(map_energy_y(e, scale, height)));
        canvas.line(axis_x - 5, y, axis_x + 5, y, 124, 135, 149, 1);
        draw_text(canvas, 80, y - 5, format_energy(e, options.energy_unit, 3),
                  101, 116, 136, 1);
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        std::array<std::uint8_t, 3> colour{41, 56, 76};
        if (level.metadata.region == OrbitalRegion::Virtual) colour = {119, 132, 150};
        if (level.homo) colour = {19, 138, 98};
        if (level.lumo) colour = {181, 106, 0};
        if (level.metadata.selected) colour = {40, 109, 224};
        const int x = static_cast<int>(std::lround(item.x));
        const int y = static_cast<int>(std::lround(item.y));
        canvas.line(x - 45, y, x + 45, y, colour[0], colour[1], colour[2],
                    level.metadata.selected ? 4 : 2);
        if (level.electrons.alpha > 0) draw_electron(canvas, x - 9, y - 2, true);
        if (level.electrons.beta > 0) draw_electron(canvas, x + 9, y - 2, false);

        const int label_y = static_cast<int>(std::lround(item.label_y));
        const int label_x = static_cast<int>(std::lround(static_cast<double>(width) * 0.72));
        if (std::abs(label_y - y) > 5 || label_x > x + 58) {
            canvas.line(x + 47, y, label_x - 8, label_y, 184, 194, 206, 1);
        }
        draw_text(canvas, label_x, label_y - 7, level.metadata.display_label,
                  colour[0], colour[1], colour[2], 2);
        const std::string annotation = annotation_summary(level.annotation);
        if (!annotation.empty()) {
            draw_text(canvas, x + 58, label_y + 9, annotation, 101, 116, 136, 1);
        }
    }
    return write_png_rgba(canvas, path, error);
}

bool write_mo_diagram_json(const MODiagramData& data,
                           const MODiagramOptions& options,
                           const std::filesystem::path& path,
                           std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open JSON output path";
        return false;
    }
    const std::set<std::size_t> included(data.selection.included_indices.begin(),
                                         data.selection.included_indices.end());
    const EnergyScale axis = energy_scale(data, energy_range(data), options);
    out << "{\n";
    out << "  \"schema\": \"cov.mo-diagram.v2\",\n";
    out << "  \"mode\": \"valence-central\",\n";
    out << "  \"energy_axis_mode\": \"" << energy_axis_mode_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_name\": \"" << energy_transform_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_parameters\": {\"focus_hartree\": " << std::setprecision(15)
        << axis.transform.focus_hartree << ", \"scale_hartree\": "
        << axis.transform.scale_hartree << "},\n";
    out << "  \"selection_summary\": \"" << json_escape(data.selection.summary) << "\",\n";
    out << "  \"shown_orbitals\": " << data.selection.included_indices.size() << ",\n";
    out << "  \"hidden_orbitals\": " << data.selection.hidden_count << ",\n";
    out << "  \"energy_unit\": \"" << energy_unit_symbol(options.energy_unit) << "\",\n";
    out << "  \"degeneracy_tolerance_hartree\": " << std::setprecision(12)
        << options.degeneracy.tolerance_hartree << ",\n";
    out << "  \"strict_salc_claimed\": false,\n";
    out << "  \"orbitals\": [\n";
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        const auto& item = data.metadata[i];
        const auto& annotation = data.annotations[i];
        out << "    {\"index\": " << item.orbital_index
            << ", \"raw_mo\": " << item.raw_mo_number
            << ", \"label\": \"" << json_escape(item.display_label) << "\""
            << ", \"source_energy_hartree\": " << std::setprecision(15) << item.energy_hartree
            << ", \"energy_hartree\": " << item.energy_hartree
            << ", \"energy_display\": " << convert_hartree(item.energy_hartree, options.energy_unit)
            << ", \"display_coordinate\": "
            << energy_display_coordinate(item.energy_hartree, axis.transform)
            << ", \"occupation\": " << item.occupation
            << ", \"spin\": \"" << spin_name(item.spin) << "\""
            << ", \"symmetry\": \"" << json_escape(item.symmetry) << "\""
            << ", \"degeneracy_size\": " << item.degeneracy_size
            << ", \"region\": \"" << region_name(item.region) << "\""
            << ", \"included_in_diagram\": " << (included.count(i) ? "true" : "false")
            << ", \"visible\": " << (item.visible ? "true" : "false")
            << ", \"selected\": " << (item.selected ? "true" : "false")
            << ", \"orbital_family\": \"" << json_escape(annotation.family) << "\""
            << ", \"family_source\": \"" << annotation_source_name(annotation.family_source) << "\""
            << ", \"family_confidence\": " << annotation.family_confidence
            << ", \"bonding_class\": \"" << bonding_class_name(annotation.bonding_class) << "\""
            << ", \"bonding_class_source\": \"" << annotation_source_name(annotation.bonding_source) << "\""
            << ", \"bonding_confidence\": " << annotation.bonding_confidence
            << ", \"heuristic\": " << (annotation.heuristic ? "true" : "false") << '}';
        if (i + 1 != data.metadata.size()) out << ',';
        out << '\n';
    }
    out << "  ]\n}\n";
    if (!out) {
        if (error) *error = "Failed while writing JSON";
        return false;
    }
    return true;
}

bool write_mo_diagram_csv(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open CSV output path";
        return false;
    }
    const std::set<std::size_t> included(data.selection.included_indices.begin(),
                                         data.selection.included_indices.end());
    const EnergyScale axis = energy_scale(data, energy_range(data), options);
    out << "internal_index,raw_mo,display_label,source_energy_hartree,energy_display,energy_unit,"
           "energy_axis_mode,transform_name,focus_hartree,scale_hartree,display_coordinate,"
           "occupation,spin,symmetry,degeneracy_size,region,included_in_diagram,visible,selected,"
           "orbital_family,family_source,family_confidence,bonding_class,bonding_class_source,"
           "bonding_confidence,heuristic\n";
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        const auto& item = data.metadata[i];
        const auto& annotation = data.annotations[i];
        out << item.orbital_index << ',' << item.raw_mo_number << ','
            << csv_escape(item.display_label) << ',' << std::setprecision(15)
            << item.energy_hartree << ','
            << convert_hartree(item.energy_hartree, options.energy_unit) << ','
            << energy_unit_symbol(options.energy_unit) << ','
            << energy_axis_mode_name(options.energy_axis_mode) << ','
            << energy_transform_name(options.energy_axis_mode) << ','
            << axis.transform.focus_hartree << ',' << axis.transform.scale_hartree << ','
            << energy_display_coordinate(item.energy_hartree, axis.transform) << ','
            << item.occupation << ','
            << spin_name(item.spin) << ',' << csv_escape(item.symmetry) << ','
            << item.degeneracy_size << ',' << region_name(item.region) << ','
            << (included.count(i) ? "true" : "false") << ','
            << (item.visible ? "true" : "false") << ','
            << (item.selected ? "true" : "false") << ','
            << csv_escape(annotation.family) << ','
            << annotation_source_name(annotation.family_source) << ','
            << annotation.family_confidence << ','
            << bonding_class_name(annotation.bonding_class) << ','
            << annotation_source_name(annotation.bonding_source) << ','
            << annotation.bonding_confidence << ','
            << (annotation.heuristic ? "true" : "false") << '\n';
    }
    if (!out) {
        if (error) *error = "Failed while writing CSV";
        return false;
    }
    return true;
}

MODiagramExportResult export_mo_diagram_bundle(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::filesystem::path& base_path) {
    MODiagramExportResult result;
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);
    std::filesystem::path base = base_path;
    if (base.has_extension()) base.replace_extension();
    std::string error;

    result.svg = write_mo_diagram_svg(data, options, base.string() + ".mo.svg", &error);
    if (!result.svg && result.error.empty()) result.error = error;
    error.clear();
    result.png = write_mo_diagram_png(data, options, base.string() + ".mo.png", &error);
    if (!result.png && result.error.empty()) result.error = error;
    error.clear();
    result.json = write_mo_diagram_json(data, options, base.string() + ".mo.json", &error);
    if (!result.json && result.error.empty()) result.error = error;
    error.clear();
    result.csv = write_mo_diagram_csv(data, options, base.string() + ".mo.csv", &error);
    if (!result.csv && result.error.empty()) result.error = error;
    return result;
}

} // namespace cov
