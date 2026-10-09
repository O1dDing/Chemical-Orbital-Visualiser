namespace {

int bitmap_text_width(const std::string& text, const int scale) {
    int width = 0;
    for (const char c : text) width += (c == ' ' ? 4 : 6) * scale;
    return width;
}

void draw_delta_pi_png(Canvas& canvas,
                       const int right_x,
                       const int centre_y,
                       const ExportPiColour colour, const bool crystal_field=false) {
    // Raster exports use vector strokes for the same compact Δπ symbol as the
    // interactive canvas.  This avoids depending on a platform bitmap font for
    // Greek glyphs and keeps the symbol legible at all supported export sizes.
    const int delta_x = right_x - 25;
    canvas.line(delta_x, centre_y + 5, delta_x + 5, centre_y - 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(delta_x + 5, centre_y - 5, delta_x + 10, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(delta_x, centre_y + 5, delta_x + 10, centre_y + 5,
                colour.r, colour.g, colour.b, 2);

    if (crystal_field) {
        draw_text(canvas,delta_x+14,centre_y-4,"CF",colour.r,colour.g,colour.b,1);
        return;
    }
    const int pi_x = delta_x + 14;
    canvas.line(pi_x, centre_y - 4, pi_x + 10, centre_y - 4,
                colour.r, colour.g, colour.b, 2);
    canvas.line(pi_x + 2, centre_y - 4, pi_x + 2, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(pi_x + 8, centre_y - 4, pi_x + 8, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
}

void draw_pi_interactions_png(Canvas& canvas,
                              const MODiagramData& data,
                              const std::vector<PositionedLevel>& positioned) {
    std::size_t rendered = 0;
    for (const auto* gap : orbital_energy_gaps(data)) {
        const auto& interaction=*gap;
        const auto bracket = export_pi_bracket(interaction, positioned);
        if (!bracket) continue;
        const ExportPiColour colour = export_pi_colour(interaction.kind);
        const int x = canvas.width - 26 - static_cast<int>(rendered) * 48;
        const int top = static_cast<int>(std::lround(bracket->top));
        const int bottom = static_cast<int>(std::lround(bracket->bottom));
        constexpr int cap = 7;
        canvas.line(x, top, x, bottom, colour.r, colour.g, colour.b, 2);
        canvas.line(x - cap, top, x, top, colour.r, colour.g, colour.b, 2);
        canvas.line(x - cap, bottom, x, bottom, colour.r, colour.g, colour.b, 2);
        draw_delta_pi_png(canvas, x - cap - 5, (top + bottom) / 2, colour,
                          interaction.gap_kind==OrbitalEnergyGapKind::CrystalField);
        ++rendered;
    }
}

} // namespace

bool write_mo_diagram_png(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    int width = std::max(760, options.width);
    const int height = std::max(720, options.height);
    const auto positioned = position_levels(data, width, height);
    Canvas canvas(width, height);

    const auto text=[&](int x,int y,const std::string& value,
                       std::uint8_t r,std::uint8_t g,std::uint8_t b,int scale) {
        if(options.raster_text) options.raster_text(canvas.rgba,width,height,x,y,value,r,g,b,scale);
        else draw_text(canvas,x,y,value,r,g,b,scale);
    };
    const auto text_width=[&](const std::string& value,int scale) {
        return options.raster_text_width ? options.raster_text_width(value,scale) : bitmap_text_width(value,scale);
    };
    text(36,28,options.figure_title.empty()?mo_diagram_mode_title(data.mode):options.figure_title,24,34,49,2);
    text(36,60,export_axis_title(options),92,106,126,1);

    const int axis_x = 70;
    canvas.line(axis_x, 110, axis_x, height - 62, 45, 55, 69, 2);
    canvas.line(axis_x, 110, axis_x - 5, 120, 45, 55, 69, 2);
    canvas.line(axis_x, 110, axis_x + 5, 120, 45, 55, 69, 2);

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_axis_fraction(fraction, data.energy_transform);
        const int y = static_cast<int>(std::lround(
            map_energy_y(energy, data.energy_transform, height)));
        canvas.line(axis_x - 5, y, axis_x + 5, y, 116, 130, 148, 1);
        const std::string tick = format_energy(
            energy, options.energy_unit,
            options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 3);
        text(axis_x + 10, y - 6, tick, 105, 119, 138, 1);
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        std::uint8_t r = 41, g = 56, b = 76;
        if (level.metadata.region == OrbitalRegion::Virtual) { r = 119; g = 132; b = 150; }
        if (level.homo) { r = 19; g = 138; b = 98; }
        if (level.lumo) { r = 181; g = 106; b = 0; }
        const int y = static_cast<int>(std::lround(item.y));
        const auto members=mo_diagram_member_views(data,level);
        for (std::size_t member_index=0;member_index<members.size();++member_index) {
            const auto& member=members[member_index];
            const int x=static_cast<int>(std::lround(item.x+
                (static_cast<double>(member_index)-
                 0.5*static_cast<double>(members.size()-1))*31.0));
            canvas.line(x-12,y,x+12,y,member.selected?40:r,
                member.selected?109:g,member.selected?224:b,member.selected?4:2);
            if (member.electrons.alpha>0) draw_electron(canvas,x-5,y-2,true);
            if (member.electrons.beta>0) draw_electron(canvas,x+5,y-2,false);
        }

        const std::string name=export_level_name(data,level,options);
        text(static_cast<int>(item.x)-text_width(name,1)/2,y+9,name,67,81,102,1);
        if (level.annotation.delocalised_pi.available) {
            std::ostringstream pi;
            pi << "PI";
            if (level.annotation.delocalised_pi.participating_atoms > 0) {
                pi << "^" << level.annotation.delocalised_pi.participating_atoms;
            }
            if (level.annotation.delocalised_pi.participating_electrons > 0.0) {
                pi << "_" << static_cast<int>(std::lround(
                    level.annotation.delocalised_pi.participating_electrons));
            }
            const std::string topology_suffix=compact_pi_topology_suffix(
                level.annotation.delocalised_pi);
            if (!topology_suffix.empty()) {
                pi << " " << topology_suffix;
            }
            const std::string value = pi.str();
            text(static_cast<int>(item.x)-text_width(value,1)/2,
                 y+25,value,91,101,117,1);
        }
    }

    draw_pi_interactions_png(canvas, data, positioned);

    return write_png_rgba(canvas,path,error,data.view?data.view->id:std::string{});
}
