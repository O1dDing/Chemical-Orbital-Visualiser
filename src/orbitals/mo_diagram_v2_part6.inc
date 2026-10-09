    data.plan.classification = DiagramClassification::Simple;
    data.plan.complex_system =
        wavefunction.atoms.size() > 12 || wavefunction.orbitals.size() > 80;
    data.plan.strict_salc_available = false;
    data.plan.machine_reason =
        "valence-central: direct Molden MO data; strict SALC reconstruction is not claimed";
    data.frontier =
        find_frontier_orbitals(wavefunction.orbitals,
                               options.filter.occupation_threshold);
    data.metadata =
        build_orbital_metadata(wavefunction,
                               options.selected_index,
                               options.degeneracy,
                               options.filter);

    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& mo : wavefunction.orbitals) {
        data.annotations.push_back(annotate_orbital(mo));
    }

    data.selection =
        build_valence_selection_plan(wavefunction, options, data.metadata);
    data.levels.reserve(data.selection.included_indices.size());
    std::vector<double> axis_energies;
    axis_energies.reserve(data.selection.included_indices.size());

    for (const std::size_t index : data.selection.included_indices) {
        if (index >= data.metadata.size()) continue;
        MODiagramLevel level;
        level.metadata = data.metadata[index];
        level.annotation = data.annotations[index];
        level.electrons =
            electron_glyphs_for_orbital(
                wavefunction.orbitals[index],
                data.frontier.separate_spin_sets);
        level.member_indices={index};
        level.member_electrons={level.electrons};
        level.total_occupation=wavefunction.orbitals[index].occupation;
        level.chemistry=wavefunction.orbitals[index].chemistry;
        level.homo = data.frontier.homo && *data.frontier.homo == index;
        level.lumo = data.frontier.lumo && *data.frontier.lumo == index;
        level.layout_energy_hartree =
            group_layout_energy(data.metadata, index);
        axis_energies.push_back(level.layout_energy_hartree);
        data.levels.push_back(std::move(level));
    }

    data.energy_transform =
        build_energy_transform(axis_energies,
                               options.energy_axis_mode,
                               options.nonlinear_minimum_gap_weight);
    return data;
}

namespace {

constexpr const char* kSvgFontFamily =
    "Segoe UI, Microsoft YaHei UI, Yu Gothic UI, Noto Sans, sans-serif";

std::string family_symbol_svg(const std::string& family) {
    if (family == "sigma") return "σ";
    if (family == "pi") return "π";
    if (family == "delta") return "δ";
    if (family == "phi") return "φ";
    return "N/A";
}

std::string family_plain(const std::string& family) {
    if (family == "sigma") return "SIGMA";
    if (family == "pi") return "PI";
    if (family == "delta") return "DELTA";
    if (family == "phi") return "PHI";
    return "N/A";
}

std::string bonding_plain(const BondingClass value) {
    switch (value) {
        case BondingClass::Bonding: return "BONDING";
        case BondingClass::Nonbonding: return "NONBONDING";
        case BondingClass::Antibonding: return "ANTIBONDING";
        case BondingClass::Mixed: return "MIXED";
        default: return "N/A";
    }
}

std::string orbital_type_plain(const OrbitalAnnotation& annotation) {
    const std::string family = family_plain(annotation.family);
    const std::string bonding = bonding_plain(annotation.bonding_class);
    std::string out = family == "N/A" && bonding == "N/A"
                          ? "N/A"
                          : family + " | " + bonding;
    if (annotation.multicentre.available) out += " | " + annotation.multicentre.label;
    if (annotation.delocalised_pi.available) out += " | LARGE-PI";
    return out;
}

std::string orbital_type_svg(const OrbitalAnnotation& annotation) {
    const std::string family = family_symbol_svg(annotation.family);
    const std::string bonding = annotation.bonding_class == BondingClass::Unclassified
                                    ? "N/A"
                                    : bonding_class_name(annotation.bonding_class);
    std::string out = family == "N/A" && bonding == "N/A"
                          ? "N/A"
                          : family + " · " + bonding;
    if (annotation.multicentre.available) {
        out += " · " + annotation.multicentre.label;
    }
    return out;
}

void write_svg_electron(std::ostream& out,
                        const double x,
                        const double level_anchor_y,
                        const bool up) {
    // Textbook convention: both spin arrows occupy the space above the energy
    // line. Alpha points upward from the near-line anchor. Beta starts above
    // the line and points downward *towards* the same anchor, so neither arrow
    // can collide with the orbital-type/symmetry text below the level.
    constexpr double length = 11.0;
    constexpr double head_dx = 3.2;
    constexpr double head_dy = 3.5;
    const double tail = up ? level_anchor_y : level_anchor_y - length;
    const double tip = up ? level_anchor_y - length : level_anchor_y;
    const double wing_y = up ? tip + head_dy : tip - head_dy;
    out << "<line x1=\"" << x << "\" y1=\"" << tail
        << "\" x2=\"" << x << "\" y2=\"" << tip
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
    out << "<line x1=\"" << x << "\" y1=\"" << tip
        << "\" x2=\"" << x - head_dx << "\" y2=\"" << wing_y
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
    out << "<line x1=\"" << x << "\" y1=\"" << tip
        << "\" x2=\"" << x + head_dx << "\" y2=\"" << wing_y
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
}

const char* svg_pi_topology_id(const DelocalisedPiTopology topology) noexcept {
    return pi_topology_name(topology);
}

void write_svg_delocalised_pi(std::ostream& out,
                              const double x,
                              const double y,
                              const DelocalisedPiDescriptor& pi,
                              const char* colour) {
    if (!pi.available) return;
    const double base_x = x - 4.0;
    out << "<g data-pi-topology-available=\""
        << (pi.topology_available ? "true" : "false")
        << "\" data-pi-topology=\"";
    if (pi.topology_available) out << svg_pi_topology_id(pi.topology);
    out << "\" data-pi-orientation-channels=\"";
    if (pi.topology_available) out << pi.orientation_channels;
    out << "\" data-pi-cyclic=\"";
    if (pi.topology_available) out << (pi.cyclic_topology ? "true" : "false");
    out << "\">\n<title>delocalised pi; ";
    if (pi.topology_available) {
        out << svg_pi_topology_id(pi.topology) << "; "
            << pi.orientation_channels << " orientation channel(s); "
            << (pi.cyclic_topology ? "cyclic pi channel" : "no cyclic pi channel claim");
    } else {
        out << "topology unavailable";
    }
    out << "</title>\n";
    std::ostringstream graph_evidence;
    write_pi_topology_graph_json(graph_evidence,pi.topology_graph);
    out << "<metadata class=\"pi-topology-graph\">"
        << xml_escape(graph_evidence.str()) << "</metadata>\n";
    out << "<text x=\"" << base_x << "\" y=\"" << y
        << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"10\">Π</text>\n";
    if (pi.participating_atoms > 0) {
        out << "<text x=\"" << base_x + 7.0 << "\" y=\"" << y - 5.0
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"7\">"
            << pi.participating_atoms << "</text>\n";
    }
    if (pi.participating_electrons > 0.0) {
        out << "<text x=\"" << base_x + 7.0 << "\" y=\"" << y + 4.0
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"7\">"
            << static_cast<int>(std::lround(pi.participating_electrons)) << "</text>\n";
    }
    const std::string topology_suffix=compact_pi_topology_suffix(pi);
    if (!topology_suffix.empty()) {
        out << "<text x=\"" << base_x + 16.0 << "\" y=\"" << y
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily
            << "\" font-size=\"7\">" << topology_suffix << "</text>\n";
    }
    out << "</g>\n";
}

struct ExportPiColour {
    const char* hex = "#e25ddc";
    std::uint8_t r = 226;
    std::uint8_t g = 93;
    std::uint8_t b = 220;
};

ExportPiColour export_pi_colour(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return {"#f2a340", 242, 163, 64};
        case PiInteractionKind::Acceptor: return {"#46ceda", 70, 206, 218};
        case PiInteractionKind::WeakNearNonbonding:
            return {"#ae9abe", 174, 154, 190};
        default: return {};
    }
}

const char* export_pi_kind_id(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return "donor";
        case PiInteractionKind::Acceptor: return "acceptor";
        case PiInteractionKind::WeakNearNonbonding: return "weak-near-nonbonding";
        default: return "coupled";
    }
}

struct ExportPiBracket {
    double top = 0.0;
    double bottom = 0.0;
};

std::optional<ExportPiBracket> export_pi_bracket(
    const PiInteractionDescriptor& interaction,
    const std::vector<PositionedLevel>& positioned) {
    if (interaction.retained_level >= positioned.size()) return std::nullopt;

    ExportPiBracket bracket;
    if ((interaction.gap_kind==OrbitalEnergyGapKind::PiPartner &&
         interaction.kind == PiInteractionKind::WeakNearNonbonding) ||
        !interaction.lower_visible || !interaction.upper_visible ||
        interaction.lower_level >= positioned.size() ||
        interaction.upper_level >= positioned.size()) {
        const double y = positioned[interaction.retained_level].y;
        bracket.top = y - 7.0;
        bracket.bottom = y + 7.0;
        return bracket;
    }

    bracket.top = std::min(positioned[interaction.lower_level].y,
                           positioned[interaction.upper_level].y) - 5.0;
    bracket.bottom = std::max(positioned[interaction.lower_level].y,
                              positioned[interaction.upper_level].y) + 5.0;
    if (bracket.bottom - bracket.top < 14.0) {
        const double middle = 0.5 * (bracket.top + bracket.bottom);
        bracket.top = middle - 7.0;
        bracket.bottom = middle + 7.0;
    }
    return bracket;
}

} // namespace

namespace {
std::string export_level_name(const MODiagramData& data, const MODiagramLevel& level,
                              const MODiagramOptions& options) {
    std::string result;
    for(const auto& member : mo_diagram_member_views(data,level)) {
        if(!result.empty()) result += ", ";
        const auto i=member.orbital_index;
        result += i<options.display_names.size() ? options.display_names[i] :
            "MO [list] "+std::to_string(i+1);
    }
    return result;
}
std::string export_axis_title(const MODiagramOptions& options) {
    return (options.axis_title.empty() ?
        std::string(options.energy_axis_mode==EnergyAxisMode::Linear ? "Energy" : "Nonlinear energy axis") :
        options.axis_title) + " (" + energy_unit_symbol(options.energy_unit) + ")";
}
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

    int width = std::max(760, options.width);
    const int height = std::max(720, options.height);
    const auto positioned = position_levels(data, width, height);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
         << "\" height=\"" << height << "\" viewBox=\"0 0 "
         << width << ' ' << height << "\" data-diagram-row-count=\""
         << data.levels.size() << "\" data-pi-interaction-count=\""
         << data.pi_interactions.size() << "\" data-crystal-field-gap-count=\""
         << data.crystal_field_gaps.size() << "\" data-diagram-mode=\""
         << mo_diagram_mode_name(data.mode) << "\" data-view-snapshot=\""
         << (data.view?xml_escape(data.view->id):"") << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#fbfcfe\"/>\n";
    out << "<text x=\"36\" y=\"40\" fill=\"#182231\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"24\" font-weight=\"600\">"
        << xml_escape(options.figure_title.empty()?mo_diagram_mode_title(data.mode):options.figure_title) << "</text>\n";

    const int axis_x = 70;
    out << "<line x1=\"" << axis_x << "\" y1=\"110\" x2=\"" << axis_x
        << "\" y2=\"" << height - 62 << "\" stroke=\"#2d3745\" stroke-width=\"2\"/>\n";
    out << "<path d=\"M " << axis_x << " 98 l -5 10 h 10 z\" fill=\"#2d3745\"/>\n";
    out << "<text x=\"24\" y=\"340\" fill=\"#4d5a6c\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"12\" transform=\"rotate(-90 24 340)\">"
        << xml_escape(export_axis_title(options)) << "</text>\n";

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_axis_fraction(fraction, data.energy_transform);
        const double y = map_energy_y(energy, data.energy_transform, height);
        out << "<line x1=\"64\" y1=\"" << y << "\" x2=\"76\" y2=\"" << y
            << "\" stroke=\"#7c8795\" stroke-width=\"1\"/>\n";
        out << "<text x=\"80\" y=\"" << y + 4
            << "\" fill=\"#7c8795\" font-family=\"" << kSvgFontFamily << "\" font-size=\"10\">"
            << xml_escape(format_energy(energy, options.energy_unit,
                   options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 4))
            << "</text>\n";
    }

    for (std::size_t row_index = 0; row_index < positioned.size(); ++row_index) {
        const auto& item = positioned[row_index];
        const auto& level = *item.level;
        const char* colour = level.homo ? "#138a62" :
                             level.lumo ? "#b56a00" :
                             level.metadata.region == OrbitalRegion::Virtual ? "#778496" : "#29384c";
        const auto members=mo_diagram_member_views(data,level);
        out << "<g class=\"mo-level\" data-row-index=\"" << row_index
            << "\" data-symmetry=\"" << xml_escape(level.metadata.symmetry_view.label)
            << "\" data-symmetry-origin=\"" << orbital_symmetry_origin_name(level.metadata.symmetry_view.origin)
            << "\" data-symmetry-explanation=\"" << xml_escape(orbital_symmetry_json(level.metadata.symmetry_view))
            << "\" data-member-count=\"" << members.size() << "\">\n";
        out << "<title>" << xml_escape(orbital_symmetry_compact_text(level.metadata.symmetry_view)) << "</title>\n";
        for (std::size_t member_index=0;member_index<members.size();++member_index) {
            const auto& member=members[member_index];
            const double x=item.x+(static_cast<double>(member_index)-
                0.5*static_cast<double>(members.size()-1))*31.0;
            out << "<line class=\"mo-member\" data-orbital-index=\"" << member.orbital_index
                << "\" data-current-name=\"" << xml_escape(member.orbital_index<options.display_names.size()?options.display_names[member.orbital_index]:"")
                << "\" data-current-irrep=\"" << xml_escape(member.orbital_index<options.display_irreps.size()?options.display_irreps[member.orbital_index]:"?")
                << "\" data-current-point-group=\"" << xml_escape(member.orbital_index<options.display_point_groups.size()?point_group_display(options.display_point_groups[member.orbital_index]):"")
                << "\" data-inspected-orbital-index=\"" << member.inspected_orbital_index
                << "\" data-selected=\"" << (member.selected?"true":"false")
                << "\" x1=\"" << x-12.0 << "\" y1=\"" << item.y
                << "\" x2=\"" << x+12.0 << "\" y2=\"" << item.y
                << "\" stroke=\"" << (member.selected?"#286de0":colour)
                << "\" stroke-width=\"" << (member.selected?4.0:2.2) << "\"/>\n";
            if (member.electrons.alpha>0) write_svg_electron(out,x-5.0,item.y-2.0,true);
            if (member.electrons.beta>0) write_svg_electron(out,x+5.0,item.y-2.0,false);
        }
        out << "</g>\n";

        // Names come from the same canonical naming table as the live viewer.
        out << "<text x=\"" << item.x << "\" y=\"" << item.y + 16
            << "\" text-anchor=\"middle\" fill=\"#435166\" font-family=\"" << kSvgFontFamily << "\" "
               "font-size=\"10\">" << xml_escape(export_level_name(data,level,options)) << "</text>\n";
        if (level.annotation.delocalised_pi.available) {
            write_svg_delocalised_pi(out, item.x, item.y + 29,
                                     level.annotation.delocalised_pi, "#5b6575");
        }
    }

    out << "<g id=\"orbital-energy-gaps\">\n";
    std::size_t rendered_interactions = 0;
    for (const bool cf:{false,true}) {
    const auto& gaps=cf?data.crystal_field_gaps:data.pi_interactions;
    out << "<g id=\"" << (cf?"crystal-field-gaps":"pi-interactions") << "\">\n";
    for (std::size_t interaction_index = 0; interaction_index < gaps.size(); ++interaction_index) {
        const auto& interaction = gaps[interaction_index];
        const char* symbol=cf?"ΔCF":"Δπ";
        const auto bracket = export_pi_bracket(interaction, positioned);
        if (!bracket) continue;
        const ExportPiColour colour = export_pi_colour(interaction.kind);
        const double x = static_cast<double>(width) - 26.0 -
                         static_cast<double>(rendered_interactions) * 48.0;
        constexpr double cap = 7.0;
        const double middle = 0.5 * (bracket->top + bracket->bottom);
        out << "<g class=\"" << (cf?"crystal-field-gap":"pi-interaction") << "\" data-interaction-index=\""
            << interaction_index << "\" data-gap-kind=\"" << orbital_energy_gap_kind_name(interaction.gap_kind)
            << "\" data-kind=\"" << (cf?"weak-crystal-field":export_pi_kind_id(interaction.kind))
            << "\" data-evidence=\"" << xml_escape(orbital_energy_gap_json(interaction,options.energy_unit)) << "\">\n";
        out << "<title>" << xml_escape(cf?"Local crystal-field separation":pi_interaction_kind_name(interaction.kind));
        if (!interaction.symmetry.empty()) {
            out << " · " << xml_escape(interaction.symmetry);
        }
        out << " · " << symbol << "=" << std::setprecision(10)
            << interaction.splitting_hartree << " Ha</title>\n";
        out << "<line x1=\"" << x << "\" y1=\"" << bracket->top
            << "\" x2=\"" << x << "\" y2=\"" << bracket->bottom
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<line x1=\"" << x - cap << "\" y1=\"" << bracket->top
            << "\" x2=\"" << x << "\" y2=\"" << bracket->top
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<line x1=\"" << x - cap << "\" y1=\"" << bracket->bottom
            << "\" x2=\"" << x << "\" y2=\"" << bracket->bottom
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<text x=\"" << x - cap - 5.0 << "\" y=\"" << middle + 4.0
            << "\" text-anchor=\"end\" fill=\"" << colour.hex
            << "\" font-family=\"" << kSvgFontFamily
            << "\" font-size=\"13\">" << symbol << "</text>\n";
        out << "</g>\n";
        ++rendered_interactions;
    }
    out << "</g>\n";
    }
    out << "</g>\n";

    out << "</svg>\n";
    if (!out) {
        if (error) *error = "Failed while writing SVG";
        return false;
    }
    return true;
}
