        zlib.push_back(static_cast<std::uint8_t>((nlen >> 8u) & 0xffu));
        zlib.insert(zlib.end(),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos + block));
        pos += block;
    }
    push_u32_be(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> png = {137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.width));
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.height));
    ihdr.push_back(8); ihdr.push_back(6); ihdr.push_back(0);
    ihdr.push_back(0); ihdr.push_back(0);
    append_png_chunk(png, "IHDR", ihdr);
    if (!view_id.empty()) {
        const std::string key="cov.view.snapshot";
        std::vector<std::uint8_t> payload(key.begin(),key.end());
        payload.push_back(0);
        payload.insert(payload.end(),view_id.begin(),view_id.end());
        append_png_chunk(png,"tEXt",payload);
    }
    append_png_chunk(png, "IDAT", zlib);
    append_png_chunk(png, "IEND", {});

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        if (error) *error = "Unable to open PNG output path";
        return false;
    }
    file.write(reinterpret_cast<const char*>(png.data()),
               static_cast<std::streamsize>(png.size()));
    if (!file) {
        if (error) *error = "Failed while writing PNG";
        return false;
    }
    return true;
}

std::string trim_symmetry(std::string value) {
    while (!value.empty() && (value.front() == ' ' || value.front() == '(')) value.erase(value.begin());
    while (!value.empty() && (value.back() == ' ' || value.back() == ')')) value.pop_back();
    return value;
}

std::string subscript_digits(std::string_view value) {
    std::string out;
    for (const char c : value) {
        switch (c) {
            case '0': out += "₀"; break;
            case '1': out += "₁"; break;
            case '2': out += "₂"; break;
            case '3': out += "₃"; break;
            case '4': out += "₄"; break;
            case '5': out += "₅"; break;
            case '6': out += "₆"; break;
            case '7': out += "₇"; break;
            case '8': out += "₈"; break;
            case '9': out += "₉"; break;
            default: out.push_back(c); break;
        }
    }
    return out;
}

} // namespace

const char* energy_axis_mode_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "linear" : "nonlinear-focus";
}

const char* energy_transform_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "identity" : "adaptive-log-gap-v3";
}

EnergyTransform build_energy_transform(const std::vector<double>& energies_hartree,
                                       const EnergyAxisMode mode,
                                       const double minimum_gap_weight) {
    EnergyTransform transform;
    transform.mode = mode;
    transform.minimum_gap_weight = std::max(0.0, minimum_gap_weight);

    std::vector<double> energies;
    energies.reserve(energies_hartree.size());
    for (const double energy : energies_hartree) {
        if (std::isfinite(energy)) energies.push_back(energy);
    }
    if (energies.empty()) {
        transform.knots = {{0.0, 0.0}, {1.0, 1.0}};
        transform.focus_hartree = 0.5;
        transform.scale_hartree = 1.0;
        return transform;
    }

    std::sort(energies.begin(), energies.end());
    energies.erase(std::unique(energies.begin(), energies.end(),
                               [](const double a, const double b) {
                                   return std::abs(a - b) <= kCoordinateEpsilon;
                               }),
                   energies.end());

    transform.focus_hartree = energies[energies.size() / 2];
    transform.scale_hartree = median_positive_gap(energies);

    if (energies.size() == 1) {
        transform.knots.push_back({energies.front(), 0.0});
        return transform;
    }

    if (mode == EnergyAxisMode::Linear) {
        transform.knots.reserve(energies.size());
        for (const double energy : energies) transform.knots.push_back({energy, energy});
        return transform;
    }

    const std::size_t gap_count = energies.size() - 1;
    std::vector<double> raw_weights(gap_count, 0.0);
    double raw_sum = 0.0;
    for (std::size_t i = 0; i < gap_count; ++i) {
        const double gap = energies[i + 1] - energies[i];
        const double raw = std::log1p(gap / std::max(1.0e-12, transform.scale_hartree));
        raw_weights[i] = std::max(kCoordinateEpsilon, raw);
        raw_sum += raw_weights[i];
    }

    // Give every distinct adjacent energy gap an explicit minimum share of
    // the axis, then distribute the remainder according to a logarithmic gap
    // weight. This is strictly monotone, reversible, and specifically avoids
    // collapsing crowded but non-degenerate valence levels. The 0.90 cap
    // leaves at least ten percent of the axis for genuine gap-size variation.
    const double maximum_floor = 0.90 / static_cast<double>(gap_count);
    const double floor_fraction = std::clamp(transform.minimum_gap_weight,
                                             0.0,
                                             maximum_floor);
    const double remaining = std::max(0.0,
        1.0 - floor_fraction * static_cast<double>(gap_count));

    transform.knots.reserve(energies.size());
    double coordinate = 0.0;
    transform.knots.push_back({energies.front(), coordinate});
    for (std::size_t i = 0; i < gap_count; ++i) {
        const double proportional = raw_sum > 0.0
            ? remaining * raw_weights[i] / raw_sum
            : remaining / static_cast<double>(gap_count);
        coordinate += floor_fraction + proportional;
        transform.knots.push_back({energies[i + 1], coordinate});
    }
    // Avoid tiny floating drift so downstream pixel metrics can treat the
    // adaptive coordinate span as exactly [0, 1].
    transform.knots.back().coordinate = 1.0;
    return transform;
}

double energy_display_coordinate(const double energy_hartree,
                                 const EnergyTransform& transform) noexcept {
    if (!std::isfinite(energy_hartree)) return 0.0;
    if (transform.knots.empty()) {
        if (transform.mode == EnergyAxisMode::Linear) return energy_hartree;
        const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
        return std::asinh((energy_hartree - transform.focus_hartree) / scale);
    }
    if (transform.knots.size() == 1) return transform.knots.front().coordinate;
    if (energy_hartree <= transform.knots.front().energy_hartree) {
        const auto& a = transform.knots[0];
        const auto& b = transform.knots[1];
        return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                           a.coordinate, b.coordinate);
    }
    if (energy_hartree >= transform.knots.back().energy_hartree) {
        const auto& a = transform.knots[transform.knots.size() - 2];
        const auto& b = transform.knots.back();
        return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                           a.coordinate, b.coordinate);
    }
    const auto it = std::upper_bound(
        transform.knots.begin(), transform.knots.end(), energy_hartree,
        [](const double value, const EnergyAxisKnot& knot) {
            return value < knot.energy_hartree;
        });
    const std::size_t hi = static_cast<std::size_t>(std::distance(transform.knots.begin(), it));
    const std::size_t lo = hi - 1;
    const auto& a = transform.knots[lo];
    const auto& b = transform.knots[hi];
    return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                       a.coordinate, b.coordinate);
}

double energy_from_display_coordinate(const double coordinate,
                                      const EnergyTransform& transform) noexcept {
    if (!std::isfinite(coordinate)) return 0.0;
    if (transform.knots.empty()) {
        if (transform.mode == EnergyAxisMode::Linear) return coordinate;
        const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
        return transform.focus_hartree + scale * std::sinh(coordinate);
    }
    if (transform.knots.size() == 1) return transform.knots.front().energy_hartree;
    if (coordinate <= transform.knots.front().coordinate) {
        const auto& a = transform.knots[0];
        const auto& b = transform.knots[1];
        return interpolate(coordinate, a.coordinate, b.coordinate,
                           a.energy_hartree, b.energy_hartree);
    }
    if (coordinate >= transform.knots.back().coordinate) {
        const auto& a = transform.knots[transform.knots.size() - 2];
        const auto& b = transform.knots.back();
        return interpolate(coordinate, a.coordinate, b.coordinate,
                           a.energy_hartree, b.energy_hartree);
    }
    const auto it = std::upper_bound(
        transform.knots.begin(), transform.knots.end(), coordinate,
        [](const double value, const EnergyAxisKnot& knot) {
            return value < knot.coordinate;
        });
    const std::size_t hi = static_cast<std::size_t>(std::distance(transform.knots.begin(), it));
    const std::size_t lo = hi - 1;
    const auto& a = transform.knots[lo];
    const auto& b = transform.knots[hi];
    return interpolate(coordinate, a.coordinate, b.coordinate,
                       a.energy_hartree, b.energy_hartree);
}

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
        default: return "N/A";
    }
}

SymmetryNotation parse_symmetry_notation(std::string_view raw_view) {
    SymmetryNotation result;
    result.raw = std::string(raw_view);
    std::string value = trim_symmetry(result.raw);
    if (value.empty() || value == "?" || value == "none" || value == "None") return result;

    const std::string lower = lower_ascii(value);
    // Gaussian abbreviates the angular-momentum symbol itself: SG is Sigma,
    // while SGG/SGU carry an additional inversion suffix. Preserve the raw
    // producer text; only its presentation is normalised. Greek irreps are
    // uppercase, distinct from the sigma/pi bonding-channel names.
    struct LinearAlias { std::string_view prefix; const char* glyph; };
    constexpr LinearAlias linear_aliases[] = {
        {"sigma","Σ"},{"sg","Σ"},{"Σ","Σ"},{"σ","Σ"},
        {"delta","Δ"},{"dlta","Δ"},{"dlt","Δ"},{"Δ","Δ"},{"δ","Δ"},
        {"phi","Φ"},{"Φ","Φ"},{"φ","Φ"},
        {"gamma","Γ"},{"Γ","Γ"},{"γ","Γ"},
        {"pi","Π"},{"Π","Π"},{"π","Π"},
    };
    for (const auto& alias : linear_aliases) {
        if (!lower.starts_with(alias.prefix)) continue;
        std::string suffix=lower.substr(alias.prefix.size());
        for (std::size_t at=suffix.find("−");at!=std::string::npos;at=suffix.find("−")) {
            suffix.replace(at,std::string_view("−").size(),"-");
        }
        bool valid=true;
        std::string parity,reflection;
        for (const unsigned char c:suffix) {
            if (std::isspace(c) || c=='_' || c=='.' || c=='^' || c=='{' || c=='}') continue;
            if ((c=='g' || c=='u') && parity.empty()) parity.push_back(static_cast<char>(c));
            else if ((c=='+' || c=='-') && reflection.empty()) reflection=c=='+'?"+":"−";
            else { valid=false;break; }
        }
        if (valid) {
            result.base=alias.glyph;
            result.subscript=parity;
            result.superscript=reflection;
        } else {
            // Unknown longer text remains literal; a recognised prefix is not
            // enough evidence to invent a valid irrep or an inversion suffix.
            result.base=value;
        }
        return result;
    }
    std::size_t pos = 0;
    {
        while (pos < value.size()) {
            const char c = value[pos];
            if ((c >= '0' && c <= '9') || c == '_' || c == '\'' || c == '"' ||
                c == '*' || c == 'g' || c == 'G' || c == 'u' || c == 'U') break;
            result.base.push_back(c);
            ++pos;
        }
    }
    if (result.base.empty() && pos < value.size()) {
        result.base.push_back(value[pos]);
        ++pos;
    }

    bool explicit_subscript = false;
    for (; pos < value.size(); ++pos) {
        const char c = value[pos];
        if (c == '_') { explicit_subscript = true; continue; }
        if (c >= '0' && c <= '9') { result.subscript.push_back(c); continue; }
        if (c == 'g' || c == 'G' || c == 'u' || c == 'U') {
            result.subscript.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
            continue;
        }
        if (c == '\'') { result.superscript += "′"; continue; }
        if (c == '"') { result.superscript += "″"; continue; }
        if (c == '*') { result.superscript += "*"; continue; }
        if (explicit_subscript) result.subscript.push_back(c);
        else result.base.push_back(c);
    }
    return result;
}

std::string point_group_display(std::string_view raw) {
    if(raw=="Dinfh")return "D∞h";
    if(raw=="Cinfv")return "C∞v";
    return std::string(raw);
}

std::string format_symmetry_unicode(std::string_view raw) {
    const SymmetryNotation notation = parse_symmetry_notation(raw);
    if (notation.base.empty()) return notation.raw.empty() ? std::string{"N/A"} : notation.raw;
    // Unicode has robust numeric subscripts but no universal subscript glyph
    // for g. GUI rich rendering uses parse_symmetry_notation() to lower the
    // entire subscript; this string is only the plain-text fallback.
    std::string digits;
    std::string letters;
    for (const char c : notation.subscript) {
        if (c >= '0' && c <= '9') digits.push_back(c);
        else letters.push_back(c);
    }
    return notation.base + subscript_digits(digits) + letters + notation.superscript;
}
