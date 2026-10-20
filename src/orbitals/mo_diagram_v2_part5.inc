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

    if (const auto ncme = parse_ncme_token(label)) {
        result.multicentre.available = true;
        result.multicentre.centres = ncme->first;
        result.multicentre.electrons = static_cast<double>(ncme->second);
        result.multicentre.source = AnnotationSource::ParsedLabel;
        result.multicentre.confidence = 1.0;
        result.multicentre.heuristic = false;
        result.multicentre.label =
            std::to_string(ncme->first) + "c" + std::to_string(ncme->second) + "e";
    }

    if (result.family == "pi" &&
        contains_any(label, {"delocal", "deloc", "large-pi", "large_pi"})) {
        result.delocalised_pi.available = true;
        result.delocalised_pi.source = AnnotationSource::ParsedLabel;
        result.delocalised_pi.confidence = 1.0;
        result.delocalised_pi.heuristic = false;
        if (result.multicentre.available) {
            result.delocalised_pi.participating_atoms = result.multicentre.centres;
            result.delocalised_pi.participating_electrons = result.multicentre.electrons;
        }
        std::ostringstream descriptor;
        descriptor << "delocalised-pi";
        if (result.delocalised_pi.participating_atoms > 0) {
            descriptor << '-' << result.delocalised_pi.participating_atoms << 'c';
        }
        if (result.delocalised_pi.participating_electrons > 0.0) {
            descriptor << '-' << result.delocalised_pi.participating_electrons << 'e';
        }
        result.delocalised_pi.label = descriptor.str();
    }

    // No Cp-, benzene, occupancy or energy special cases. If producer text
    // does not support a classification, the UI deliberately shows N/A and
    // machine metadata records an unavailable source.
    result.heuristic =
        result.family_source == AnnotationSource::Heuristic ||
        result.bonding_source == AnnotationSource::Heuristic ||
        result.multicentre.heuristic ||
        result.delocalised_pi.heuristic;
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

    const FrontierOrbitals frontier =
        find_frontier_orbitals(wavefunction.orbitals,
                               options.filter.occupation_threshold);
    const std::size_t max_levels = effective_max_levels(options);
    const std::size_t max_virtual = std::min(
        options.max_virtual_levels,
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
    chosen.insert(chosen.end(),
                  occupied_valence.begin() + static_cast<std::ptrdiff_t>(occ_begin),
                  occupied_valence.end());
    chosen.insert(chosen.end(), frontier_virtual.begin(),
                  frontier_virtual.begin() + static_cast<std::ptrdiff_t>(virtual_keep));

    if (options.selected_index < wavefunction.orbitals.size() &&
        metadata[options.selected_index].region != OrbitalRegion::Core) {
        chosen.push_back(options.selected_index);
    }

    std::sort(chosen.begin(), chosen.end());
    if(options.aomo_scope>=2) {
        chosen.clear();
        for(std::size_t i=0;i<metadata.size();++i)
            if(options.aomo_scope==3 || options.show_core_background ||
               metadata[i].region!=OrbitalRegion::Core)chosen.push_back(i);
    }
    chosen.erase(std::unique(chosen.begin(), chosen.end()), chosen.end());
    const auto labels = build_orbital_labels(
        wavefunction.orbitals,
        point_group_limited_degeneracy(wavefunction, options.degeneracy));
    chosen = expand_degenerate_sets(chosen, labels);

    plan.included_indices = std::move(chosen);
    plan.hidden_count = metadata.size() > plan.included_indices.size()
                            ? metadata.size() - plan.included_indices.size()
                            : 0;
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
