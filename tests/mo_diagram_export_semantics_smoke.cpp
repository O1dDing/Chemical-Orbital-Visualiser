#include "cov/mo_diagram.hpp"

#include <array>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace {

std::string read_all(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

std::size_t count_token(const std::string& text, const std::string& token) {
    std::size_t count = 0;
    std::size_t offset = 0;
    while ((offset = text.find(token, offset)) != std::string::npos) {
        ++count;
        offset += token.size();
    }
    return count;
}

bool number_after(const std::string& text,const std::string& marker,double expected) {
    const auto pos=text.find(marker);
    if (pos==std::string::npos) return false;
    try { return std::abs(std::stod(text.substr(pos+marker.size()))-expected)<1e-12; }
    catch (...) { return false; }
}

bool numeric_array_after(const std::string& text,const std::string& marker,
                         const std::vector<double>& expected) {
    const auto pos=text.find(marker);
    if (pos==std::string::npos) return false;
    auto begin=text.find('[',pos+marker.size());
    if (begin==std::string::npos) return false;
    ++begin;
    for (std::size_t i=0;i<expected.size();++i) {
        try {
            std::size_t used=0;
            if (std::abs(std::stod(text.substr(begin),&used)-expected[i])>=1e-12) return false;
            begin+=used;
            const auto delimiter=text.find_first_not_of(" \t\r\n",begin);
            if (delimiter==std::string::npos || text[delimiter]!=(i+1==expected.size()?']':',')) return false;
            begin=delimiter+1;
        } catch (...) { return false; }
    }
    return true;
}

cov::OrbitalMetadata metadata(const std::size_t index,
                              const double energy,
                              const char* symmetry) {
    cov::OrbitalMetadata result;
    result.orbital_index = index;
    result.raw_mo_number = index + 1;
    result.display_label = std::to_string(index + 1);
    result.energy_hartree = energy;
    result.occupation = index == 0 ? 2.0f : 0.0f;
    result.symmetry = symmetry;
    // Direct export fixture: model a supplied label with unresolved axes.
    // The exporter consumes the scoped value, while preserving the raw label.
    result.symmetry_view.label = symmetry;
    result.symmetry_view.origin = cov::OrbitalSymmetryOrigin::Producer;
    result.symmetry_view.orbital_indices = {index};
    result.symmetry_view.source_path = "synthetic-export-fixture";
    result.region = index == 0 ? cov::OrbitalRegion::Valence
                               : cov::OrbitalRegion::Virtual;
    return result;
}

} // namespace

int main() {
    cov::MODiagramData data;
    data.selection.included_indices = {0, 1};
    data.selection.summary = "compact validation rows";
    data.metadata = {metadata(0, -0.20, "T2g"),
                     metadata(1, 0.10, "T2g*")};
    data.annotations.resize(2);
    data.annotations[0].family = "pi";
    data.annotations[1].family = "pi";
    auto& derived_pi=data.annotations[0].delocalised_pi;
    derived_pi.available=true;
    derived_pi.participating_atoms=6u;
    derived_pi.participating_electrons=6.0;
    derived_pi.atom_indices={0,1,2,3,4,5};
    derived_pi.orbital_indices={0};
    derived_pi.family_id="pi-family-derived";
    derived_pi.topology_available=true;
    derived_pi.topology=cov::DelocalisedPiTopology::HapticMetal;
    derived_pi.orientation_channels=2u;
    derived_pi.cyclic_topology=true;
    derived_pi.orientation_channel_details={
        {{0,1,2,3,4,5},{0.0,0.0,1.0},0.98,true},
        {{0,1},{1.0,0.0,0.0},0.94,false},
    };
    derived_pi.label="PI6_6";
    derived_pi.source=cov::AnnotationSource::Derived;
    if (cov::compact_pi_topology_suffix(derived_pi)!="2ch H") {
        std::cerr<<"shared compact pi topology label mismatch\n";
        return 11;
    }
    auto& parsed_pi=data.annotations[1].delocalised_pi;
    parsed_pi.available=true;
    parsed_pi.family_id="pi-family-parsed-only";
    parsed_pi.label="delocalised-pi";
    parsed_pi.source=cov::AnnotationSource::ParsedLabel;
    auto& shared_multicentre=data.annotations[0].multicentre;
    shared_multicentre.available=true;
    shared_multicentre.label="2x3c2e (shared 4c/4e source)";
    shared_multicentre.centres=4u;
    shared_multicentre.electrons=4.0;
    shared_multicentre.atom_indices={0,1,2,3};
    shared_multicentre.channel_count=2u;
    shared_multicentre.source_subspace_id="mc-source-bridge-pair";
    shared_multicentre.source_subspace_electron_count=4.0;

    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        cov::MODiagramLevel level;
        level.metadata = data.metadata[i];
        level.annotation = data.annotations[i];
        level.layout_energy_hartree = data.metadata[i].energy_hartree;
        level.member_indices = {i};
        level.electrons = i == 0 ? cov::ElectronGlyphs{1, 1}
                                 : cov::ElectronGlyphs{};
        data.levels.push_back(level);
    }
    cov::MODiagramLevel collapsed_udft_level;
    collapsed_udft_level.member_indices={2u,3u};
    collapsed_udft_level.member_spin_counterparts={6u,7u};
    if (!cov::mo_diagram_level_covers_orbital(collapsed_udft_level,2u) ||
        !cov::mo_diagram_level_covers_orbital(collapsed_udft_level,7u) ||
        cov::mo_diagram_level_covers_orbital(collapsed_udft_level,8u)) {
        std::cerr<<"collapsed UDFT pi-family row coverage mismatch\n";
        return 12;
    }
    data.energy_transform = cov::build_energy_transform(
        {-0.20, 0.10}, cov::EnergyAxisMode::NonlinearFocus, 0.055);

    data.ligand_field_point_group = "Oh";
    data.ligand_field_geometry_id = "OC-6";
    data.ligand_field_geometry_name = "Octahedron";
    data.ligand_field_coordination_number = 6;
    data.ligand_field_metal_atom = 4;
    data.ligand_field_ligand_atoms = {0, 1, 2, 3, 5, 6};
    data.ligand_field_confidence = 0.93;
    data.ligand_field_angular_rms = 0.0125;
    data.ligand_field_shape_measure = 0.018;
    data.ligand_field_radial_cv = 0.021;
    cov::LocalGeometryDiagramDescriptor local_geometry;
    local_geometry.centre_atom=4u;
    local_geometry.geometry_id="OC-6";
    local_geometry.geometry_name="Octahedron";
    local_geometry.point_group="Oh";
    local_geometry.neighbour_atoms={0,1,2,3,5,6};
    local_geometry.confidence=0.93;
    local_geometry.angular_rms=0.0125;
    local_geometry.shape_measure=0.018;
    local_geometry.radial_cv=0.021;
    data.local_geometries.push_back(local_geometry);
    data.electronic_state.source=cov::WavefunctionSource::Fchk;
    data.electronic_state.charge=-1;
    data.electronic_state.multiplicity=2u;
    data.electronic_state.alpha_electrons=5u;
    data.electronic_state.beta_electrons=4u;
    data.electronic_state.charge_provenance=cov::DataProvenance::Producer;
    data.electronic_state.multiplicity_provenance=cov::DataProvenance::Producer;
    data.electronic_state.electron_counts_provenance=cov::DataProvenance::Producer;
    data.electronic_state.scf_convergence=cov::ScfConvergenceStatus::Converged;
    data.electronic_state.scf_convergence_provenance=cov::DataProvenance::Producer;
    data.electronic_state.stability=cov::WavefunctionStabilityStatus::Stable;
    data.electronic_state.stability_provenance=cov::DataProvenance::Producer;
    data.electronic_state.stability_detail="stable internal test";
    data.electronic_state.spin_squared_before=0.76;
    data.electronic_state.spin_squared_after=0.75;
    data.electronic_state.spin_squared_provenance=cov::DataProvenance::Producer;
    data.electronic_state.atomic_partial_charge_scheme="Mulliken";
    data.electronic_state.atomic_partial_charges={-0.2,0.2};
    data.electronic_state.atomic_partial_charge_provenance=
        cov::DataProvenance::Producer;
    data.electronic_state.point_group_detected="Oh";
    data.electronic_state.point_group_used="Oh";
    data.electronic_state.point_group_provenance=cov::DataProvenance::Producer;
    data.electronic_state.source_title="export semantics fixture";
    data.electronic_state.source_route="#p test route";
    data.electronic_state.enrichment_source="fixture.log";

    cov::PiInteractionDescriptor pi;
    pi.lower_level = 0;
    pi.upper_level = 1;
    pi.retained_level = 0;
    pi.lower_orbitals = {0, 2, 4};
    pi.upper_orbitals = {1, 3, 5};
    pi.symmetry = "T2g";
    pi.kind = cov::PiInteractionKind::Acceptor;
    pi.splitting_hartree = 0.30;
    pi.lower_energy_hartree = -0.20;
    pi.upper_energy_hartree = 0.10;
    pi.confidence = 0.87;
    data.pi_interactions.push_back(pi);

    cov::MODiagramOptions options;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 760;
    options.height = 720;
    options.include_hidden_in_metadata = false;

    // CPU and CUDA build trees may execute this contract concurrently.  Keep
    // their transient artifacts disjoint so one process cannot replace an SVG
    // while the other is reading it.
    const auto build_id=std::hash<std::string>{}(
        std::filesystem::absolute(std::filesystem::current_path()).string());
    const auto base = std::filesystem::temp_directory_path() /
        ("cov_mo_diagram_export_semantics_smoke_"+
         std::to_string(build_id));
    const auto svg_path = std::filesystem::path(base.string() + ".svg");
    const auto png_path = std::filesystem::path(base.string() + ".png");
    const auto json_path = std::filesystem::path(base.string() + ".json");
    const auto csv_path = std::filesystem::path(base.string() + ".csv");
    std::string error;
    if (!cov::write_mo_diagram_svg(data, options, svg_path, &error) ||
        !cov::write_mo_diagram_png(data, options, png_path, &error) ||
        !cov::write_mo_diagram_json(data, options, json_path, &error) ||
        !cov::write_mo_diagram_csv(data, options, csv_path, &error)) {
        std::cerr << "export failed: " << error << '\n';
        return 1;
    }

    const std::string svg = read_all(svg_path);
    if (svg.find("data-diagram-row-count=\"2\"") == std::string::npos ||
        svg.find("data-diagram-mode=\"valence-central\"")==std::string::npos ||
        svg.find("Valence MO diagram")==std::string::npos ||
        svg.find("data-symmetry=\"T2g\"")==std::string::npos ||
        svg.find(">T</text>")!=std::string::npos ||
        count_token(svg, "class=\"mo-level\"") != data.levels.size() ||
        count_token(svg, "class=\"pi-interaction\"") != 1 ||
        svg.find("data-kind=\"acceptor\"") == std::string::npos ||
        svg.find("data-pi-topology-available=\"true\"") == std::string::npos ||
        svg.find("data-pi-topology=\"haptic-metal\"") == std::string::npos ||
        svg.find("data-pi-orientation-channels=\"2\"") == std::string::npos ||
        svg.find("data-pi-cyclic=\"true\"") == std::string::npos ||
        svg.find("data-pi-topology-available=\"false\"") == std::string::npos ||
        svg.find("splitting; approximately") != std::string::npos) {
        std::cerr << "SVG compact-row or pi-marker semantics mismatch\n";
        return 2;
    }

    // Screen/SVG semantics must preserve the actual chemistry topology, not
    // collapse every family to the coarser cyclic/non-cyclic flag.
    const std::array<std::pair<cov::DelocalisedPiTopology,const char*>,4>
        topology_contracts{{
            {cov::DelocalisedPiTopology::HapticMetal,"haptic-metal"},
            {cov::DelocalisedPiTopology::Spiro,"spiro"},
            {cov::DelocalisedPiTopology::SymmetryDirectSum,
             "symmetry-direct-sum"},
            {cov::DelocalisedPiTopology::MultiChannel,"multi-channel"},
        }};
    for (const auto& [topology,machine]:topology_contracts) {
        data.levels[0].annotation.delocalised_pi.topology=topology;
        if (!cov::write_mo_diagram_svg(data,options,svg_path,&error)) {
            std::cerr<<"topology-specific SVG export failed: "<<error<<'\n';
            return 9;
        }
        const auto topology_svg=read_all(svg_path);
        if (topology_svg.find(
                std::string("data-pi-topology=\"")+machine+"\"")==
                std::string::npos ||
            topology_svg.find(std::string("<title>delocalised pi; ")+machine+
                              ";")==std::string::npos) {
            std::cerr<<"SVG lost concrete pi topology "<<machine<<'\n';
            return 10;
        }
    }
    data.levels[0].annotation.delocalised_pi.topology=
        cov::DelocalisedPiTopology::HapticMetal;

    const std::string png = read_all(png_path);
    const std::string acceptor_rgb{
        static_cast<char>(70), static_cast<char>(206),
        static_cast<char>(218), static_cast<char>(255)};
    if (png.find(acceptor_rgb) == std::string::npos) {
        std::cerr << "PNG acceptor side marker colour missing\n";
        return 3;
    }

    const std::string json = read_all(json_path);
    if (json.find("\"mode\": \"valence-central\"")==std::string::npos ||
        json.find("\"local_coordination\": {") == std::string::npos ||
        json.find("\"geometry_id\": \"OC-6\"") == std::string::npos ||
        json.find("\"geometry_name\": \"Octahedron\"") == std::string::npos ||
        json.find("\"coordination_number\": 6") == std::string::npos ||
        json.find("\"angular_rms\": 0.0125") == std::string::npos ||
        !number_after(json,"\"shape_measure\":",0.018) ||
        json.find("\"radial_cv\": 0.021") == std::string::npos ||
        json.find("\"local_molecular_geometries\": [") == std::string::npos ||
        json.find("\"centre_atom_index\":4") == std::string::npos ||
        json.find("\"neighbour_atom_indices\":[0,1,2,3,5,6]") ==
            std::string::npos ||
        json.find("\"electronic_state\": {") == std::string::npos ||
        json.find("\"charge\": -1") == std::string::npos ||
        json.find("\"multiplicity\": 2") == std::string::npos ||
        json.find("\"alpha_electrons\": 5") == std::string::npos ||
        json.find("\"beta_electrons\": 4") == std::string::npos ||
        json.find("\"scf_convergence\": \"converged\"") == std::string::npos ||
        json.find("\"wavefunction_stability\": \"stable\"") == std::string::npos ||
        !numeric_array_after(json,"\"atomic_partial_charges\":",{-0.2,0.2}) ||
        json.find("\"pi_interactions\": [") == std::string::npos ||
        json.find("\"kind\": \"acceptor\"") == std::string::npos ||
        !number_after(json,"\"splitting_hartree\":",0.3) ||
        json.find("\"lower_orbitals\": [0,2,4]") == std::string::npos ||
        json.find("\"family_id\":\"pi-family-derived\"") == std::string::npos ||
        json.find("\"topology_available\":true") == std::string::npos ||
        json.find("\"topology\":\"haptic-metal\"") == std::string::npos ||
        json.find("\"orientation_channels\":2") == std::string::npos ||
        json.find("\"direction\":[0,0,1]") == std::string::npos ||
        json.find("\"delocalised_pi_topology_available\": false") == std::string::npos ||
        json.find("\"delocalised_pi_topology\": null") == std::string::npos ||
        json.find("\"delocalised_pi_orientation_channels\": null") == std::string::npos ||
        json.find("\"delocalised_pi_cyclic_topology\": null") == std::string::npos ||
        json.find("\"multicentre_channel_count\": 2") == std::string::npos ||
        json.find("\"multicentre_source_subspace_id\": \"mc-source-bridge-pair\"") ==
            std::string::npos ||
        json.find("\"multicentre_source_subspace_electron_count\": 4") ==
            std::string::npos ||
        json.find("\"diagram_row_count\": 2") == std::string::npos) {
        std::cerr << "JSON local coordination or pi metadata missing\n";
        return 4;
    }

    const std::string csv = read_all(csv_path);
    if (csv.find("diagram_mode,index,raw_mo")==std::string::npos ||
        csv.find("valence-central,0,1,")==std::string::npos ||
        csv.find("local_geometry_id,local_geometry_name,local_coordination_number") ==
            std::string::npos ||
        csv.find("local_angular_rms,local_shape_measure,local_radial_cv") ==
            std::string::npos ||
        csv.find("diagram_row_count,pi_interaction_count,pi_interactions") ==
            std::string::npos ||
        csv.find("delocalised_pi_participating_electrons,delocalised_pi_label") ==
            std::string::npos ||
        csv.find("pi_interactions,delocalised_pi_topology_available,"
                 "delocalised_pi_topology,") ==
            std::string::npos ||
        csv.find("delocalised_pi_orientation_channel_details") ==
            std::string::npos ||
        csv.find("local_geometry_count,local_geometries") == std::string::npos ||
        csv.find("centre=4;geometry_id=OC-6") == std::string::npos ||
        csv.find("wavefunction_source,charge,charge_provenance") ==
            std::string::npos ||
        csv.find("stability_detail,source_title,source_route,enrichment_source") ==
            std::string::npos ||
        csv.find(",FCHK,-1,producer,2,producer,5,4,producer,converged,producer,stable,") ==
            std::string::npos ||
        csv.find(",Oh,OC-6,Octahedron,6,") == std::string::npos ||
        csv.find("kind=acceptor;symmetry=T2g") == std::string::npos ||
        !number_after(csv,"splitting_hartree=",0.3)) {
        std::cerr << "CSV local coordination or pi metadata missing\n";
        return 5;
    }
    if (csv.find("multicentre_channel_count,multicentre_source_subspace_id,") ==
            std::string::npos ||
        csv.find(",2,mc-source-bridge-pair,4,stable internal test,"
                 "export semantics fixture,#p test route,fixture.log,") ==
            std::string::npos) {
        std::cerr << "CSV shared multicentre source metadata missing\n";
        return 6;
    }

    struct ModeContract {
        cov::MODiagramMode mode;
        const char* machine;
        const char* title;
    };
    const std::array<ModeContract,2> compact_modes{{
        {cov::MODiagramMode::DelocalisedPiFamilyOnly,
         "delocalised-pi-family-only","Delocalised pi MO diagram"},
        {cov::MODiagramMode::MulticentreActiveSpaceOnly,
         "multicentre-active-space-only","Multicentre active-space MO diagram"},
    }};
    for (const auto& contract:compact_modes) {
        data.mode=contract.mode;
        if (!cov::write_mo_diagram_svg(data,options,svg_path,&error) ||
            !cov::write_mo_diagram_json(data,options,json_path,&error) ||
            !cov::write_mo_diagram_csv(data,options,csv_path,&error)) {
            std::cerr<<"mode-specific export failed: "<<error<<'\n';
            return 7;
        }
        const auto mode_svg=read_all(svg_path);
        const auto mode_json=read_all(json_path);
        const auto mode_csv=read_all(csv_path);
        if (mode_svg.find(std::string("data-diagram-mode=\"")+contract.machine+
                          "\"")==std::string::npos ||
            mode_svg.find(contract.title)==std::string::npos ||
            mode_json.find(std::string("\"mode\": \"")+contract.machine+
                           "\"")==std::string::npos ||
            mode_csv.find(std::string("diagram_mode,index,raw_mo"))==
                std::string::npos ||
            mode_csv.find(std::string(contract.machine)+",0,1,")==
                std::string::npos) {
            std::cerr<<"mode-specific export labels diverged\n";
            return 8;
        }
    }

    // Canonical display identity must be independent of the row's local
    // label and the original producer spelling; no relabelling by dimension.
    options.display_names={"1πg","2πg"};
    options.display_irreps={"πg","πg"};
    options.display_point_groups={"Dinfh","Dinfh"};
    const auto raw_before=data.metadata[0].symmetry;
    const auto scoped_before=data.metadata[0].symmetry_view.label;
    if (!cov::write_mo_diagram_json(data,options,json_path,&error) ||
        !cov::write_mo_diagram_csv(data,options,csv_path,&error) ||
        !cov::write_mo_diagram_svg(data,options,svg_path,&error)) return 12;
    const auto named_json=read_all(json_path);
    const auto named_csv=read_all(csv_path);
    const auto named_svg=read_all(svg_path);
    if(named_json.find("\"label\": \"1πg\"")==std::string::npos ||
       named_json.find("\"current_display_irrep\": \"πg\"")==std::string::npos ||
       named_json.find("\"current_display_point_group\": \"D∞h\"")==std::string::npos ||
       named_json.find("\"source_symmetry\": \"T2g\"")==std::string::npos ||
       named_csv.find("current_display_irrep,current_display_point_group")==std::string::npos ||
       named_csv.find("πg,D∞h")==std::string::npos ||
       named_svg.find("1πg")==std::string::npos ||
       named_svg.find("data-symmetry=\"T2g\"")==std::string::npos ||
       data.metadata[0].symmetry!=raw_before || data.metadata[0].symmetry_view.label!=scoped_before ||
       cov::point_group_display("Cinfv")!="C∞v" || cov::point_group_display("D5h")!="D5h") {
        std::cerr<<"current canonical display and raw/local evidence were conflated\n";
        return 13;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(csv_path, ec);
    std::cout << "mo_diagram_export_semantics_smoke ok\n";
    return 0;
}
