namespace {

const char* machine_bonding_class(const BondingClass value) noexcept {
    switch (value) {
        case BondingClass::Bonding: return "bonding";
        case BondingClass::Nonbonding: return "nonbonding";
        case BondingClass::Antibonding: return "antibonding";
        case BondingClass::Mixed: return "mixed";
        default: return "unclassified";
    }
}

const char* machine_pi_topology(
    const DelocalisedPiTopology value) noexcept {
    return pi_topology_name(value);
}

const char* machine_scf_status(const ScfConvergenceStatus value) noexcept {
    switch (value) {
        case ScfConvergenceStatus::Converged: return "converged";
        case ScfConvergenceStatus::Failed: return "failed";
        default: return "unavailable";
    }
}

const char* machine_stability_status(
    const WavefunctionStabilityStatus value) noexcept {
    switch (value) {
        case WavefunctionStabilityStatus::Stable: return "stable";
        case WavefunctionStabilityStatus::Unstable: return "unstable";
        default: return "unavailable";
    }
}

void write_json_index_array(std::ostream& out,
                            const std::vector<std::size_t>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i];
    }
    out << ']';
}

void write_json_pi_channels(
    std::ostream& out,
    const std::vector<PiOrientationChannel>& channels) {
    out << '[';
    for (std::size_t index=0;index<channels.size();++index) {
        if (index) out << ',';
        const auto& channel=channels[index];
        out << "{\"atoms\":[";
        for (std::size_t atom=0;atom<channel.atoms.size();++atom) {
            if (atom) out << ',';
            out << channel.atoms[atom];
        }
        out << "],\"direction\":[" << std::setprecision(17)
            << channel.direction[0] << ',' << channel.direction[1] << ','
            << channel.direction[2] << "],\"coherence\":"
            << channel.coherence << ",\"cyclic\":"
            << (channel.cyclic ? "true" : "false") << '}';
    }
    out << ']';
}

std::string csv_index_list(const std::vector<std::size_t>& values) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ';';
        out << values[i];
    }
    return out.str();
}

std::string csv_pi_interaction_summary(const MODiagramData& data) {
    std::ostringstream out;
    out << std::setprecision(17);
    for (std::size_t i = 0; i < data.pi_interactions.size(); ++i) {
        if (i) out << " || ";
        const auto& interaction = data.pi_interactions[i];
        out << "kind=" << export_pi_kind_id(interaction.kind)
            << ";symmetry=" << interaction.symmetry
            << ";lower_level=" << interaction.lower_level
            << ";upper_level=" << interaction.upper_level
            << ";retained_level=" << interaction.retained_level
            << ";lower_orbitals=" << csv_index_list(interaction.lower_orbitals)
            << ";upper_orbitals=" << csv_index_list(interaction.upper_orbitals)
            << ";splitting_hartree=" << interaction.splitting_hartree
            << ";confidence=" << interaction.confidence
            << ";lower_visible=" << (interaction.lower_visible ? 1 : 0)
            << ";upper_visible=" << (interaction.upper_visible ? 1 : 0);
    }
    return out.str();
}

std::string csv_pi_channel_summary(const DelocalisedPiDescriptor& descriptor) {
    if (!descriptor.topology_available) return {};
    std::ostringstream out;
    out << std::setprecision(17);
    for (std::size_t index=0;index<descriptor.orientation_channel_details.size();++index) {
        if (index) out << " || ";
        const auto& channel=descriptor.orientation_channel_details[index];
        out << "atoms=";
        for (std::size_t atom=0;atom<channel.atoms.size();++atom) {
            if (atom) out << ';';
            out << channel.atoms[atom];
        }
        out << ";direction=" << channel.direction[0] << ';'
            << channel.direction[1] << ';' << channel.direction[2]
            << ";coherence=" << channel.coherence
            << ";cyclic=" << (channel.cyclic ? 1 : 0);
    }
    return out.str();
}

std::string csv_local_geometry_summary(const MODiagramData& data) {
    std::ostringstream out;
    out << std::setprecision(17);
    for (std::size_t index=0;index<data.local_geometries.size();++index) {
        if (index) out << " || ";
        const auto& geometry=data.local_geometries[index];
        out << "centre=" << geometry.centre_atom
            << ";geometry_id=" << geometry.geometry_id
            << ";name=" << geometry.geometry_name
            << ";point_group=" << geometry.point_group
            << ";neighbours=" << csv_index_list(geometry.neighbour_atoms)
            << ";confidence=" << geometry.confidence
            << ";angular_rms=" << geometry.angular_rms
            << ";shape_measure=" << geometry.shape_measure
            << ";radial_cv=" << geometry.radial_cv;
    }
    return out.str();
}

std::string csv_real_list(const std::vector<double>& values) {
    std::ostringstream out;
    out << std::setprecision(17);
    for (std::size_t index=0;index<values.size();++index) {
        if (index) out << ';';
        out << values[index];
    }
    return out.str();
}

} // namespace

bool write_mo_diagram_json(const MODiagramData& data,
                           const MODiagramOptions& options,
                           const std::filesystem::path& path,
                           std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open JSON output path";
        return false;
    }

    out << "{\n";
    out << "  \"schema\": \"cov.mo-diagram.v2\",\n";
    out << "  \"view_snapshot\": ";
    if (!data.view) out << "null";
    else {
        out << "{\"id\":\"" << json_escape(data.view->id)
            << "\",\"origin\":\"" << json_escape(data.view->origin)
            << "\",\"selection_anchor\":" << data.view->selection_anchor
            << ",\"inspected_orbital_index\":";
        if (data.view->inspected_orbital_index) out << *data.view->inspected_orbital_index;
        else out << "null";
        out << '}';
    }
    out << ",\n";
    out << "  \"mode\": \"" << mo_diagram_mode_name(data.mode) << "\",\n";
    out << "  \"strict_salc_claimed\": false,\n";
    out << "  \"energy_unit\": \"" << json_escape(energy_unit_symbol(options.energy_unit)) << "\",\n";
    out << "  \"energy_axis_mode\": \"" << energy_axis_mode_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_name\": \"" << energy_transform_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_parameters\": {\n";
    out << "    \"minimum_gap_fraction\": " << std::setprecision(12)
        << data.energy_transform.minimum_gap_weight << ",\n";
    out << "    \"scale_hartree\": " << data.energy_transform.scale_hartree << ",\n";
    out << "    \"focus_hartree\": " << data.energy_transform.focus_hartree << ",\n";
    out << "    \"knots\": [";
    for (std::size_t i = 0; i < data.energy_transform.knots.size(); ++i) {
        if (i) out << ',';
        const auto& knot = data.energy_transform.knots[i];
        out << "{\"energy_hartree\":" << knot.energy_hartree
            << ",\"coordinate\":" << knot.coordinate << '}';
    }
    out << "]\n  },\n";
    out << "  \"selection_summary\": \"" << json_escape(data.selection.summary) << "\",\n";
    out << "  \"shown_orbitals\": " << data.selection.included_indices.size() << ",\n";
    out << "  \"diagram_row_count\": " << data.levels.size() << ",\n";
    out << "  \"hidden_orbitals\": " << data.selection.hidden_count << ",\n";
    out << "  \"valence_occupied_count\": " << data.selection.valence_occupied_count << ",\n";
    out << "  \"frontier_virtual_count\": " << data.selection.frontier_virtual_count << ",\n";
    const auto& state=data.electronic_state;
    out << "  \"electronic_state\": {\n";
    out << "    \"wavefunction_source\": \""
        << wavefunction_source_name(state.source) << "\",\n";
    out << "    \"charge\": ";
    if (state.charge_provenance==DataProvenance::Unavailable) out << "null";
    else out << state.charge;
    out << ",\n    \"charge_provenance\": \""
        << data_provenance_name(state.charge_provenance) << "\",\n";
    out << "    \"multiplicity\": ";
    if (state.multiplicity_provenance==DataProvenance::Unavailable) out << "null";
    else out << state.multiplicity;
    out << ",\n    \"multiplicity_provenance\": \""
        << data_provenance_name(state.multiplicity_provenance) << "\",\n";
    out << "    \"alpha_electrons\": ";
    if (state.electron_counts_provenance==DataProvenance::Unavailable) out << "null";
    else out << state.alpha_electrons;
    out << ",\n    \"beta_electrons\": ";
    if (state.electron_counts_provenance==DataProvenance::Unavailable) out << "null";
    else out << state.beta_electrons;
    out << ",\n    \"electron_counts_provenance\": \""
        << data_provenance_name(state.electron_counts_provenance) << "\",\n";
    out << "    \"scf_convergence\": \""
        << machine_scf_status(state.scf_convergence) << "\",\n";
    out << "    \"scf_convergence_provenance\": \""
        << data_provenance_name(state.scf_convergence_provenance) << "\",\n";
    out << "    \"wavefunction_stability\": \""
        << machine_stability_status(state.stability) << "\",\n";
    out << "    \"stability_provenance\": \""
        << data_provenance_name(state.stability_provenance) << "\",\n";
    out << "    \"stability_detail\": \""
        << json_escape(state.stability_detail) << "\",\n";
    out << "    \"spin_squared_before_annihilation\": ";
    if (state.spin_squared_provenance==DataProvenance::Unavailable) out << "null";
    else out << std::setprecision(17) << state.spin_squared_before;
    out << ",\n    \"spin_squared_after_annihilation\": ";
    if (state.spin_squared_provenance==DataProvenance::Unavailable) out << "null";
    else out << std::setprecision(17) << state.spin_squared_after;
    out << ",\n    \"spin_squared_provenance\": \""
        << data_provenance_name(state.spin_squared_provenance) << "\",\n";
    out << "    \"atomic_partial_charge_scheme\": \""
        << json_escape(state.atomic_partial_charge_scheme) << "\",\n";
    out << "    \"atomic_partial_charge_provenance\": \""
        << data_provenance_name(state.atomic_partial_charge_provenance) << "\",\n";
    out << "    \"atomic_partial_charges\": ";
    if (state.atomic_partial_charge_provenance==DataProvenance::Unavailable) {
        out << "null";
    } else {
        out << '[';
        for (std::size_t index=0;index<state.atomic_partial_charges.size();++index) {
            if (index) out << ',';
            out << std::setprecision(17) << state.atomic_partial_charges[index];
        }
        out << ']';
    }
    out << ",\n    \"point_group_detected\": \""
        << json_escape(state.point_group_detected) << "\",\n";
    out << "    \"point_group_used\": \""
        << json_escape(state.point_group_used) << "\",\n";
    out << "    \"point_group_provenance\": \""
        << data_provenance_name(state.point_group_provenance) << "\",\n";
    out << "    \"point_group_detected_provenance\": \""
        << data_provenance_name(state.point_group_detected_provenance) << "\",\n";
    out << "    \"point_group_used_provenance\": \""
        << data_provenance_name(state.point_group_used_provenance) << "\",\n";
    out << "    \"source_title\": \"" << json_escape(state.source_title)
        << "\",\n    \"source_route\": \"" << json_escape(state.source_route)
        << "\",\n    \"enrichment_source\": \""
        << json_escape(state.enrichment_source) << "\"\n";
    out << "  },\n";
    out << "  \"local_coordination\": {\n";
    out << "    \"available\": "
        << (!data.ligand_field_geometry_id.empty() &&
            data.ligand_field_coordination_number > 0 ? "true" : "false") << ",\n";
    out << "    \"point_group\": \""
        << json_escape(data.ligand_field_point_group) << "\",\n";
    out << "    \"geometry_id\": \""
        << json_escape(data.ligand_field_geometry_id) << "\",\n";
    out << "    \"geometry_name\": \""
        << json_escape(data.ligand_field_geometry_name) << "\",\n";
    out << "    \"coordination_number\": "
        << data.ligand_field_coordination_number << ",\n";
    out << "    \"metal_atom_index\": " << data.ligand_field_metal_atom << ",\n";
    out << "    \"ligand_atom_indices\": ";
    write_json_index_array(out, data.ligand_field_ligand_atoms);
    out << ",\n";
    out << "    \"confidence\": " << std::setprecision(17)
        << data.ligand_field_confidence << ",\n";
    out << "    \"angular_rms\": " << data.ligand_field_angular_rms << ",\n";
    out << "    \"shape_measure\": " << data.ligand_field_shape_measure << ",\n";
    out << "    \"radial_cv\": " << data.ligand_field_radial_cv << "\n";
    out << "  },\n";
    out << "  \"local_molecular_geometries\": [\n";
    for (std::size_t index=0;index<data.local_geometries.size();++index) {
        const auto& geometry=data.local_geometries[index];
        out << "    {\"centre_atom_index\":" << geometry.centre_atom
            << ",\"geometry_id\":\"" << json_escape(geometry.geometry_id)
            << "\",\"geometry_name\":\"" << json_escape(geometry.geometry_name)
            << "\",\"point_group\":\"" << json_escape(geometry.point_group)
            << "\",\"neighbour_atom_indices\":";
        write_json_index_array(out,geometry.neighbour_atoms);
        out << ",\"confidence\":" << std::setprecision(17)
            << geometry.confidence << ",\"angular_rms\":"
            << geometry.angular_rms << ",\"shape_measure\":"
            << geometry.shape_measure << ",\"radial_cv\":"
            << geometry.radial_cv << '}'
            << (index+1u<data.local_geometries.size()?",":"") << "\n";
    }
    out << "  ],\n";
    out << "  \"pi_partner_candidates\": " << pi_partner_candidates_json(data.pi_partner_candidates) << ",\n";
    out << "  \"sigma_framework\": " << mo_sigma_framework_json(data.sigma_framework) << ",\n";
    out << "  \"using_ro_common_energy\": " << (data.using_ro_common_energy?"true":"false") << ",\n";
    out << "  \"ro_common_energy\": " << serialize_nbo_ro_common_energy_json(data.ro_common_energy) << ",\n";
    out << "  \"pi_field_response\": " << pi_field_response_analysis_json(data.pi_field_response) << ",\n";
    out << "  \"composition_scope\": {\"applicable\":" << (data.composition_scope.applicable?"true":"false")
        << ",\"centres\":";
    write_json_index_array(out,data.composition_scope.centre_atoms);
    out << ",\"ligands\":";write_json_index_array(out,data.composition_scope.ligand_atoms);
    out << ",\"other_atoms\":";write_json_index_array(out,data.composition_scope.other_atoms);
    out << ",\"reason\":\"" << json_escape(data.composition_scope.detail) << "\"},\n";
    out << "  \"pi_mode_networks\": [";
    for(std::size_t i=0;i<data.pi_mode_networks.size();++i) {
        if(i)out << ',';
        out << pi_mode_network_assessment_json(data.pi_mode_networks[i]);
    }
    out << "],\n";
    out << "  \"pi_interactions\": " << orbital_energy_gap_array_json(data.pi_interactions,options.energy_unit) << ",\n";
    out << "  \"crystal_field_gaps\": " << orbital_energy_gap_array_json(data.crystal_field_gaps,options.energy_unit) << ",\n";
    out << "  \"delocalised_pi_families\": [\n";
    std::set<std::string> exported_pi_families;
    bool first_pi_family=true;
    for (const auto& annotation:data.annotations) {
        const auto& pi=annotation.delocalised_pi;
        if (!pi.available || pi.family_id.empty() ||
            !exported_pi_families.insert(pi.family_id).second) {
            continue;
        }
        if (!first_pi_family) out << ",\n";
        first_pi_family=false;
        out << "    {\"family_id\":\"" << json_escape(pi.family_id)
            << "\",\"atom_indices\":";
        write_json_index_array(out,pi.atom_indices);
        out << ",\"orbital_indices\":";
        write_json_index_array(out,pi.orbital_indices);
        out << ",\"electron_count\":" << std::setprecision(17)
            << pi.participating_electrons
            << ",\"topology_available\":"
            << (pi.topology_available ? "true" : "false")
            << ",\"topology\":";
        if (pi.topology_available) {
            out << '"' << machine_pi_topology(pi.topology) << '"';
        } else {
            out << "null";
        }
        out << ",\"orientation_channels\":";
        if (pi.topology_available) out << pi.orientation_channels;
        else out << "null";
        out << ",\"cyclic_topology\":";
        if (pi.topology_available) out << (pi.cyclic_topology ? "true" : "false");
        else out << "null";
        out << ",\"channel_details\":";
        if (pi.topology_available) {
            write_json_pi_channels(out,pi.orientation_channel_details);
        } else {
            out << "null";
        }
        out << ",\"topology_graph\":";
        write_pi_topology_graph_json(out,pi.topology_graph);
        out << '}';
    }
    out << "\n  ],\n";
    out << "  \"diagram_rows\": [\n";
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        out << "    {\"row_index\":" << i
            << ",\"layout_energy_hartree\":" << std::setprecision(17)
            << level.layout_energy_hartree
            << ",\"symmetry\":\"" << json_escape(level.metadata.symmetry_view.label)
            << "\",\"member_orbitals\":";
        if (level.member_indices.empty()) {
            write_json_index_array(out, std::vector<std::size_t>{level.metadata.orbital_index});
        } else {
            write_json_index_array(out, level.member_indices);
        }
        const auto members=mo_diagram_member_views(data,level);
        out << ",\"symmetry_explanation\":" << orbital_symmetry_json(level.metadata.symmetry_view);
        out << ",\"composition\":" << mo_group_composition_json(level.composition)
            << ",\"display_decision\":" << mo_group_display_decision_json(level.display_decision);
        out << ",\"representative_metadata_index\":" << level.metadata.orbital_index
            << ",\"layout_energy_meaning\":\"row-layout-value-not-an-individual-MO-energy\""
            << ",\"energy_spread_hartree\":" << level.energy_spread_hartree
            << ",\"total_occupation\":" << level.total_occupation
            << ",\"members\":[";
        std::optional<double> minimum_energy,maximum_energy;
        for (std::size_t member_index=0;member_index<members.size();++member_index) {
            if (member_index) out << ',';
            const auto& member=members[member_index];
            out << "{\"orbital_index\":" << member.orbital_index
                << ",\"spin_counterpart\":";
            if (member.spin_counterpart) out << *member.spin_counterpart;
            else out << "null";
            out << ",\"inspected_orbital_index\":" << member.inspected_orbital_index
                << ",\"selected\":" << (member.selected?"true":"false")
                << ",\"alpha_arrows\":" << member.electrons.alpha
                << ",\"beta_arrows\":" << member.electrons.beta << '}';
            for (const auto& raw:data.metadata) {
                if (raw.orbital_index!=member.orbital_index &&
                    member.spin_counterpart!=raw.orbital_index) continue;
                if (!std::isfinite(raw.energy_hartree)) continue;
                minimum_energy=minimum_energy?std::min(*minimum_energy,raw.energy_hartree):raw.energy_hartree;
                maximum_energy=maximum_energy?std::max(*maximum_energy,raw.energy_hartree):raw.energy_hartree;
            }
        }
        out << "],\"all_members_energy_min_hartree\":";
        if (minimum_energy) out << *minimum_energy;
        else out << "null";
        out << ",\"all_members_energy_max_hartree\":";
        if (maximum_energy) out << *maximum_energy;
        else out << "null";
        out << '}' << (i + 1 < data.levels.size() ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"orbitals\": [\n";

    bool first = true;
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        const auto row=mo_diagram_row_for_orbital(data,i);
        const bool inspected=data.view && data.view->inspected_orbital_index==i;
        if (!options.include_hidden_in_metadata && !row && !inspected) continue;
        if (!first) out << ",\n";
        first = false;
        const auto& item = data.metadata[i];
        const auto& annotation = i < data.annotations.size() ? data.annotations[i] : OrbitalAnnotation{};
        const double coordinate = energy_display_coordinate(
            row?data.levels[*row].layout_energy_hartree:data.metadata[i].energy_hartree,
            data.energy_transform);
        out << "    {\n";
        out << "      \"index\": " << item.orbital_index << ",\n";
        out << "      \"raw_mo\": " << item.raw_mo_number << ",\n";
        out << "      \"label\": \"" << json_escape(i<options.display_names.size()?options.display_names[i]:item.display_label) << "\",\n";
        out << "      \"current_display_irrep\": \"" << json_escape(i<options.display_irreps.size()?options.display_irreps[i]:"?") << "\",\n";
        out << "      \"current_display_point_group\": \"" << json_escape(i<options.display_point_groups.size()?point_group_display(options.display_point_groups[i]):"") << "\",\n";
        out << "      \"energy_hartree\": " << std::setprecision(17) << item.energy_hartree << ",\n";
        out << "      \"source_energy_hartree\": " << item.energy_hartree << ",\n";
        out << "      \"energy_display\": " << convert_hartree(item.energy_hartree, options.energy_unit) << ",\n";
        out << "      \"display_coordinate\": " << coordinate << ",\n";
        out << "      \"diagram_row_index\": ";
        if (row) out << *row;
        else out << "null";
        out << ",\n      \"coordinate_scope\": \"" << (row?"row-layout":"unplotted-energy-map") << "\",\n";
        out << "      \"occupation\": " << item.occupation << ",\n";
        out << "      \"spin\": \"" << spin_name(item.spin) << "\",\n";
        out << "      \"symmetry\": \"" << json_escape(item.symmetry_view.label) << "\",\n";
        out << "      \"symmetry_display\": \"" << json_escape(orbital_symmetry_compact_text(item.symmetry_view)) << "\",\n";
        out << "      \"source_symmetry\": \"" << json_escape(item.symmetry) << "\",\n";
        out << "      \"symmetry_explanation\": " << orbital_symmetry_json(item.symmetry_view) << ",\n";
        out << "      \"source_symmetry_explanation\": ";
        out << (item.molecular_member_symmetries.empty() ? "null" :
            orbital_symmetry_json(item.molecular_member_symmetries.front())) << ",\n";
        out << "      \"degeneracy_size\": " << item.degeneracy_size << ",\n";
        out << "      \"region\": \"" << region_name(item.region) << "\",\n";
        out << "      \"included_in_diagram\": " << (row ? "true" : "false") << ",\n";
        out << "      \"visible\": " << (item.visible ? "true" : "false") << ",\n";
        out << "      \"selected\": " << (item.selected ? "true" : "false") << ",\n";
        out << "      \"orbital_family\": \"" << json_escape(annotation.family) << "\",\n";
        out << "      \"family_source\": \"" << annotation_source_name(annotation.family_source) << "\",\n";
        out << "      \"family_confidence\": " << annotation.family_confidence << ",\n";
        out << "      \"bonding_class\": \"" << machine_bonding_class(annotation.bonding_class) << "\",\n";
        out << "      \"bonding_class_source\": \"" << annotation_source_name(annotation.bonding_source) << "\",\n";
        out << "      \"bonding_confidence\": " << annotation.bonding_confidence << ",\n";
        out << "      \"multicentre_available\": " << (annotation.multicentre.available ? "true" : "false") << ",\n";
        out << "      \"multicentre_centres\": " << annotation.multicentre.centres << ",\n";
        out << "      \"multicentre_electrons\": " << annotation.multicentre.electrons << ",\n";
        out << "      \"multicentre_label\": \"" << json_escape(annotation.multicentre.label) << "\",\n";
        out << "      \"multicentre_channel_count\": "
            << annotation.multicentre.channel_count << ",\n";
        out << "      \"multicentre_source_subspace_id\": \""
            << json_escape(annotation.multicentre.source_subspace_id) << "\",\n";
        out << "      \"multicentre_source_subspace_electron_count\": "
            << annotation.multicentre.source_subspace_electron_count << ",\n";
        out << "      \"multicentre_source\": \"" << annotation_source_name(annotation.multicentre.source) << "\",\n";
        out << "      \"multicentre_confidence\": " << annotation.multicentre.confidence << ",\n";
        out << "      \"multicentre_heuristic\": " << (annotation.multicentre.heuristic ? "true" : "false") << ",\n";
        out << "      \"delocalised_pi_available\": " << (annotation.delocalised_pi.available ? "true" : "false") << ",\n";
        out << "      \"delocalised_pi_participating_atoms\": " << annotation.delocalised_pi.participating_atoms << ",\n";
        out << "      \"delocalised_pi_participating_electrons\": " << annotation.delocalised_pi.participating_electrons << ",\n";
        out << "      \"delocalised_pi_topology_available\": "
            << (annotation.delocalised_pi.topology_available ? "true" : "false") << ",\n";
        out << "      \"delocalised_pi_topology\": ";
        if (annotation.delocalised_pi.topology_available) {
            out << '"' << machine_pi_topology(
                annotation.delocalised_pi.topology) << '"';
        } else {
            out << "null";
        }
        out << ",\n";
        out << "      \"delocalised_pi_orientation_channels\": ";
        if (annotation.delocalised_pi.topology_available) {
            out << annotation.delocalised_pi.orientation_channels;
        } else {
            out << "null";
        }
        out << ",\n      \"delocalised_pi_cyclic_topology\": ";
        if (annotation.delocalised_pi.topology_available) {
            out << (annotation.delocalised_pi.cyclic_topology ? "true" : "false");
        } else {
            out << "null";
        }
        out << ",\n      \"delocalised_pi_orientation_channel_details\": ";
        if (annotation.delocalised_pi.topology_available) {
            write_json_pi_channels(
                 out,annotation.delocalised_pi.orientation_channel_details);
        } else {
            out << "null";
        }
        out << ",\n";
        out << "      \"delocalised_pi_topology_graph\": ";
        write_pi_topology_graph_json(out,annotation.delocalised_pi.topology_graph);
        out << ",\n";
        out << "      \"delocalised_pi_label\": \"" << json_escape(annotation.delocalised_pi.label) << "\",\n";
        out << "      \"delocalised_pi_source\": \"" << annotation_source_name(annotation.delocalised_pi.source) << "\",\n";
        out << "      \"delocalised_pi_confidence\": " << annotation.delocalised_pi.confidence << ",\n";
        out << "      \"delocalised_pi_heuristic\": " << (annotation.delocalised_pi.heuristic ? "true" : "false") << ",\n";
        out << "      \"heuristic\": " << (annotation.heuristic ? "true" : "false") << "\n";
        out << "    }";
    }
    out << "\n  ]\n}\n";
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
    out << "diagram_mode,index,raw_mo,label,energy_hartree,energy_display,energy_unit,display_coordinate,"
           "energy_axis_mode,transform_name,symmetry,symmetry_display,source_symmetry,symmetry_explanation_json,occupation,spin,"
           "degeneracy_size,region,included_in_diagram,visible,selected,orbital_family,"
           "family_source,family_confidence,bonding_class,bonding_source,bonding_confidence,"
           "multicentre_available,multicentre_centres,multicentre_electrons,multicentre_label,"
           "multicentre_source,multicentre_confidence,multicentre_heuristic,"
           "delocalised_pi_available,delocalised_pi_participating_atoms,"
           "delocalised_pi_participating_electrons,delocalised_pi_label,delocalised_pi_source,"
           "delocalised_pi_confidence,delocalised_pi_heuristic,heuristic,"
           "local_point_group,local_geometry_id,local_geometry_name,local_coordination_number,"
           "local_metal_atom_index,local_ligand_atom_indices,local_confidence,local_angular_rms,"
           "local_shape_measure,local_radial_cv,diagram_row_count,pi_interaction_count,"
           "pi_interactions,delocalised_pi_topology_available,delocalised_pi_topology,"
           "delocalised_pi_orientation_channels,delocalised_pi_cyclic_topology,"
           "delocalised_pi_orientation_channel_details,local_geometry_count,"
           "local_geometries,wavefunction_source,charge,charge_provenance,"
           "multiplicity,multiplicity_provenance,alpha_electrons,beta_electrons,"
           "electron_counts_provenance,scf_convergence,scf_convergence_provenance,"
           "wavefunction_stability,stability_provenance,spin_squared_before,"
           "spin_squared_after,spin_squared_provenance,atomic_partial_charge_scheme,"
           "atomic_partial_charge_provenance,atomic_partial_charges,point_group_detected,"
           "point_group_used,point_group_provenance,point_group_detected_provenance,point_group_used_provenance,multicentre_channel_count,"
           "multicentre_source_subspace_id,multicentre_source_subspace_electron_count,"
           "stability_detail,source_title,source_route,enrichment_source,"
           "view_snapshot_id,view_origin,view_selection_anchor,view_inspected_orbital_index,"
           "diagram_row_index,coordinate_scope,row_member_orbitals,row_spin_counterparts,delocalised_pi_topology_graph,crystal_field_gap_count,crystal_field_gaps_json,pi_partner_evidence_json,current_display_irrep,current_display_point_group,pi_partner_candidates_json\n";
    const std::string ligand_atoms = csv_index_list(data.ligand_field_ligand_atoms);
    const std::string pi_interactions = csv_pi_interaction_summary(data);
    const std::string local_geometries = csv_local_geometry_summary(data);
    const auto& state=data.electronic_state;
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        const auto row=mo_diagram_row_for_orbital(data,i);
        const bool inspected=data.view && data.view->inspected_orbital_index==i;
        if (!options.include_hidden_in_metadata && !row && !inspected) continue;
        const auto& item = data.metadata[i];
        const auto& annotation = i < data.annotations.size() ? data.annotations[i] : OrbitalAnnotation{};
        const double coordinate = energy_display_coordinate(
            row?data.levels[*row].layout_energy_hartree:data.metadata[i].energy_hartree,
            data.energy_transform);
        out << mo_diagram_mode_name(data.mode) << ','
            << item.orbital_index << ',' << item.raw_mo_number << ','
            << csv_escape(i<options.display_names.size()?options.display_names[i]:item.display_label) << ',' << std::setprecision(17) << item.energy_hartree << ','
            << convert_hartree(item.energy_hartree, options.energy_unit) << ','
            << csv_escape(energy_unit_symbol(options.energy_unit)) << ',' << coordinate << ','
            << energy_axis_mode_name(options.energy_axis_mode) << ','
            << energy_transform_name(options.energy_axis_mode) << ','
            << csv_escape(item.symmetry_view.label) << ',' << csv_escape(orbital_symmetry_compact_text(item.symmetry_view)) << ','
            << csv_escape(item.symmetry) << ',' << csv_escape(orbital_symmetry_json(item.symmetry_view)) << ','
            << item.occupation << ',' << spin_name(item.spin) << ',' << item.degeneracy_size << ','
            << region_name(item.region) << ',' << (row ? 1 : 0) << ','
            << (item.visible ? 1 : 0) << ',' << (item.selected ? 1 : 0) << ','
            << csv_escape(annotation.family) << ',' << annotation_source_name(annotation.family_source) << ','
            << annotation.family_confidence << ',' << machine_bonding_class(annotation.bonding_class) << ','
            << annotation_source_name(annotation.bonding_source) << ',' << annotation.bonding_confidence << ','
            << (annotation.multicentre.available ? 1 : 0) << ',' << annotation.multicentre.centres << ','
            << annotation.multicentre.electrons << ',' << csv_escape(annotation.multicentre.label) << ','
            << annotation_source_name(annotation.multicentre.source) << ',' << annotation.multicentre.confidence << ','
            << (annotation.multicentre.heuristic ? 1 : 0) << ','
            << (annotation.delocalised_pi.available ? 1 : 0) << ','
            << annotation.delocalised_pi.participating_atoms << ','
            << annotation.delocalised_pi.participating_electrons << ','
            << csv_escape(annotation.delocalised_pi.label) << ','
            << annotation_source_name(annotation.delocalised_pi.source) << ','
            << annotation.delocalised_pi.confidence << ','
            << (annotation.delocalised_pi.heuristic ? 1 : 0) << ','
            << (annotation.heuristic ? 1 : 0) << ','
            << csv_escape(data.ligand_field_point_group) << ','
            << csv_escape(data.ligand_field_geometry_id) << ','
            << csv_escape(data.ligand_field_geometry_name) << ','
            << data.ligand_field_coordination_number << ','
            << data.ligand_field_metal_atom << ','
            << csv_escape(ligand_atoms) << ','
            << data.ligand_field_confidence << ','
            << data.ligand_field_angular_rms << ','
            << data.ligand_field_shape_measure << ','
            << data.ligand_field_radial_cv << ','
            << data.levels.size() << ','
            << data.pi_interactions.size() << ','
            << csv_escape(pi_interactions) << ','
            << (annotation.delocalised_pi.topology_available ? 1 : 0) << ',';
        if (annotation.delocalised_pi.topology_available) {
            out << machine_pi_topology(annotation.delocalised_pi.topology);
        }
        out << ',';
        if (annotation.delocalised_pi.topology_available) {
            out << annotation.delocalised_pi.orientation_channels << ','
                << (annotation.delocalised_pi.cyclic_topology ? 1 : 0);
        } else {
            out << ',';
        }
        out << ',' << csv_escape(csv_pi_channel_summary(annotation.delocalised_pi))
            << ',' << data.local_geometries.size() << ','
            << csv_escape(local_geometries) << ','
            << wavefunction_source_name(state.source) << ',';
        if (state.charge_provenance!=DataProvenance::Unavailable) out << state.charge;
        out << ',' << data_provenance_name(state.charge_provenance) << ',';
        if (state.multiplicity_provenance!=DataProvenance::Unavailable) {
            out << state.multiplicity;
        }
        out << ',' << data_provenance_name(state.multiplicity_provenance) << ',';
        if (state.electron_counts_provenance!=DataProvenance::Unavailable) {
            out << state.alpha_electrons;
        }
        out << ',';
        if (state.electron_counts_provenance!=DataProvenance::Unavailable) {
            out << state.beta_electrons;
        }
        out << ',' << data_provenance_name(state.electron_counts_provenance) << ','
            << machine_scf_status(state.scf_convergence) << ','
            << data_provenance_name(state.scf_convergence_provenance) << ','
            << machine_stability_status(state.stability) << ','
            << data_provenance_name(state.stability_provenance) << ',';
        if (state.spin_squared_provenance!=DataProvenance::Unavailable) {
            out << state.spin_squared_before;
        }
        out << ',';
        if (state.spin_squared_provenance!=DataProvenance::Unavailable) {
            out << state.spin_squared_after;
        }
        out << ',' << data_provenance_name(state.spin_squared_provenance) << ','
            << csv_escape(state.atomic_partial_charge_scheme) << ','
            << data_provenance_name(state.atomic_partial_charge_provenance) << ',';
        if (state.atomic_partial_charge_provenance!=DataProvenance::Unavailable) {
            out << csv_escape(csv_real_list(state.atomic_partial_charges));
        }
        out << ',' << csv_escape(state.point_group_detected) << ','
            << csv_escape(state.point_group_used) << ','
            << data_provenance_name(state.point_group_provenance) << ','
            << data_provenance_name(state.point_group_detected_provenance) << ','
            << data_provenance_name(state.point_group_used_provenance) << ','
            << annotation.multicentre.channel_count << ','
            << csv_escape(annotation.multicentre.source_subspace_id) << ','
            << annotation.multicentre.source_subspace_electron_count << ','
            << csv_escape(state.stability_detail) << ','
            << csv_escape(state.source_title) << ','
            << csv_escape(state.source_route) << ','
            << csv_escape(state.enrichment_source) << ',';
        if (data.view) {
            out << csv_escape(data.view->id) << ',' << csv_escape(data.view->origin)
                << ',' << data.view->selection_anchor << ',';
            if (data.view->inspected_orbital_index) out << *data.view->inspected_orbital_index;
        } else out << ",,,";
        out << ',';
        if (row) out << *row;
        out << ',' << (row?"row-layout":"unplotted-energy-map") << ',';
        if (row) {
            std::vector<std::size_t> primary,counterparts;
            for (const auto& member:mo_diagram_member_views(data,data.levels[*row])) {
                primary.push_back(member.orbital_index);
                if (member.spin_counterpart) counterparts.push_back(*member.spin_counterpart);
            }
            out << csv_escape(csv_index_list(primary)) << ',' << csv_escape(csv_index_list(counterparts));
        } else out << ',';
        std::ostringstream topology_evidence;
        write_pi_topology_graph_json(topology_evidence,annotation.delocalised_pi.topology_graph);
        out << ',' << csv_escape(topology_evidence.str());
        out << ',' << data.crystal_field_gaps.size()
            << ',' << csv_escape(orbital_energy_gap_array_json(data.crystal_field_gaps,options.energy_unit))
            << ',' << csv_escape(orbital_energy_gap_array_json(data.pi_interactions,options.energy_unit))
            << ',' << csv_escape(i<options.display_irreps.size()?options.display_irreps[i]:"?")
            << ',' << csv_escape(i<options.display_point_groups.size()?point_group_display(options.display_point_groups[i]):"")
            << ',' << csv_escape(pi_partner_candidates_json(data.pi_partner_candidates));
        out << '\n';
    }
    if (!out) {
        if (error) *error = "Failed while writing CSV";
        return false;
    }
    return true;
}

MODiagramExportResult export_mo_diagram_bundle(const MODiagramViewSnapshot& snapshot,
                                               const std::filesystem::path& base_path,
                                               const DiagramExportContent content) {
    static std::mutex bundle_mutex;
    static std::atomic<std::uint64_t> generation{0u};
    const std::lock_guard lock(bundle_mutex);

    MODiagramExportResult result;
    const MODiagramData& data=snapshot.data;
    const MODiagramOptions& options=snapshot.options;
    std::filesystem::path normalised_base=base_path;
    if (normalised_base.has_extension()) normalised_base.replace_extension();
    result.svg_path=normalised_base.string()+".mo.svg";
    result.png_path=normalised_base.string()+".mo.png";
    result.json_path=normalised_base.string()+".mo.json";
    result.csv_path=normalised_base.string()+".mo.csv";
    const std::array<std::filesystem::path,4> targets{
        result.svg_path,result.png_path,result.json_path,result.csv_path};

    const bool images=content!=DiagramExportContent::AnalysisData;
    const bool analysis=content!=DiagramExportContent::Images;
    const std::array<bool,4> enabled{images,images,analysis,analysis};
    std::error_code ec;
    for (std::size_t i=0;i<targets.size();++i) {
        if(!enabled[i]) continue;
        const auto& target=targets[i];
        if (std::filesystem::exists(target,ec) &&
            !std::filesystem::is_regular_file(target,ec)) {
            result.error="Export target exists and is not a regular file: "+
                target.string();
            return result;
        }
        ec.clear();
    }

    const auto time_tag=static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    const std::string suffix=".covtmp."+std::to_string(time_tag)+"."+
        std::to_string(generation.fetch_add(1u,std::memory_order_relaxed));
    std::array<std::filesystem::path,4> temporary{};
    std::array<std::filesystem::path,4> backup{};
    for (std::size_t i=0;i<targets.size();++i) {
        temporary[i]=targets[i].string()+suffix;
        backup[i]=targets[i].string()+suffix+".bak";
    }
    const auto cleanup_temporaries=[&] {
        for (const auto& path:temporary) {
            ec.clear();
            std::filesystem::remove(path,ec);
        }
    };
    std::string error;
    if ((images && (!write_mo_diagram_svg(data,options,temporary[0],&error) ||
                    !write_mo_diagram_png(data,options,temporary[1],&error))) ||
        (analysis && (!write_mo_diagram_json(data,options,temporary[2],&error) ||
                      !write_mo_diagram_csv(data,options,temporary[3],&error)))) {
        cleanup_temporaries();
        result.error=error;
        return result;
    }

    std::array<bool,4> backed_up{};
    std::array<bool,4> installed{};
    for (std::size_t i=0;i<targets.size();++i) {
        if(!enabled[i]) continue;
        ec.clear();
        if (std::filesystem::exists(targets[i],ec)) {
            std::filesystem::rename(targets[i],backup[i],ec);
            if (ec) {
                result.error="Unable to stage existing export: "+ec.message();
                break;
            }
            backed_up[i]=true;
        }
        ec.clear();
        std::filesystem::rename(temporary[i],targets[i],ec);
        if (ec) {
            result.error="Unable to install export bundle: "+ec.message();
            break;
        }
        installed[i]=true;
    }
    if (!result.error.empty()) {
        for (std::size_t reverse=targets.size();reverse>0u;--reverse) {
            const std::size_t i=reverse-1u;
            if (installed[i]) {
                ec.clear();
                std::filesystem::remove(targets[i],ec);
            }
            if (backed_up[i]) {
                ec.clear();
                std::filesystem::rename(backup[i],targets[i],ec);
            }
        }
        cleanup_temporaries();
        return result;
    }
    for (std::size_t i=0;i<targets.size();++i) {
        if (!backed_up[i]) continue;
        ec.clear();
        std::filesystem::remove(backup[i],ec);
    }
    result.svg=result.png=images;
    result.json=result.csv=analysis;
    return result;
}

MODiagramExportResult export_mo_diagram_bundle(const Wavefunction& wavefunction,
                                               const MODiagramOptions& options,
                                               const std::filesystem::path& base_path) {
    const auto data=build_mo_diagram_data(wavefunction,options);
    return export_mo_diagram_bundle(make_mo_diagram_view_snapshot(
        data,options,options.selected_index,"computed-export"),base_path);
}

} // namespace cov
