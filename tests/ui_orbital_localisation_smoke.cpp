#include "cov/orbital_ui_text.hpp"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>

namespace {

using cov::AnnotationSource;
using cov::BondingClass;
using cov::DataProvenance;
using cov::MODiagramData;
using cov::MODiagramLevel;
using cov::MODiagramMode;
using cov::OrbitalBondingRole;
using cov::PiInteractionDescriptor;
using cov::PiInteractionKind;
using cov::WavefunctionSource;
using cov::ui::Language;
using cov::ui::OrbitalText;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "ui_orbital_localisation_smoke: " << message << '\n';
    std::exit(1);
}

void require(const bool condition, const std::string& message) {
    if (!condition) fail(message);
}

void require_contains(const std::string_view value,
                      const std::string_view needle,
                      const std::string& context) {
    require(value.find(needle) != std::string_view::npos,
            context + " missing '" + std::string(needle) + "'");
}

} // namespace

int main() {
    constexpr std::array languages{
        Language::English,
        Language::ChineseSimplified,
        Language::Japanese,
        Language::French,
    };

    for (const auto language : languages) {
        for (std::size_t key = 0;
             key < static_cast<std::size_t>(OrbitalText::Count); ++key) {
            const char* text = cov::ui::orbital_tr(
                static_cast<OrbitalText>(key), language);
            require(text != nullptr && *text != '\0',
                    "empty orbital localisation entry");
        }
        require(std::string_view(cov::ui::orbital_ui_glyph_seed(language)).size() > 100u,
                "glyph seed is unexpectedly short");
        for (const auto key : {OrbitalText::OrbitalDetails,
                               OrbitalText::CloseOrbitalDetails,
                               OrbitalText::OrbitalDetailsHint,
                               OrbitalText::OrbitalDetailsOutsideDiagram,
                               OrbitalText::OrbitalDetailsDataScope,
                               OrbitalText::LevelGroupContainsMOs}) {
            const std::string_view text = cov::ui::orbital_tr(key, language);
            require_contains(cov::ui::orbital_ui_glyph_seed(language), text,
                             "orbital details glyph seed");
            if (language != Language::English) {
                require(text != cov::ui::orbital_tr(key, Language::English),
                        "orbital details still uses the English fallback");
            }
        }
        require_contains(cov::ui::orbital_tr(OrbitalText::LevelGroupContainsMOs, language),
                         "%zu", "level group member count format");
        const std::string_view representative=cov::ui::orbital_tr(OrbitalText::GroupRepresentativeData,language);
        require_contains(representative,"%s","level group representative name format");
        require(representative.find("%zu")==std::string_view::npos,
                "representative name must not use an integer format");
        std::array<char,512> representative_text{};
        const auto length=std::snprintf(representative_text.data(),representative_text.size(),
            representative.data(),"2a1 [beta]");
        require(length>0 && static_cast<std::size_t>(length)<representative_text.size(),
                "representative name formatting failed");
        require_contains(representative_text.data(),"2a1 [beta]","representative name and spin");
    }

    // A group can have a measured zero overlap from opposite signed member
    // contributions. Missing member evidence must not be averaged as a zero,
    // and populated numeric fields cannot establish applicability by themselves.
    cov::Wavefunction group_source;
    group_source.atoms.resize(2);
    group_source.atoms[0].atomic_number=24;
    group_source.atoms[1].atomic_number=6;
    group_source.orbitals.resize(2);
    for (std::size_t i=0;i<2;++i) {
        auto& chemistry=group_source.orbitals[i].chemistry;
        chemistry.available=true;
        chemistry.ao_contributions.push_back({0,3,2,"3d",0.2});
        chemistry.ao_contributions.push_back({1,2,1,"2p",0.8});
        cov::OrbitalPairInteraction pair;
        pair.atom_a=0;pair.atom_b=1;
        pair.overlap_character=i==0?0.02:-0.02;
        pair.channel.status=cov::ChemistryStatus::Determined;
        pair.channel.pi=1.0;
        pair.channel.undetermined=0.0;
        chemistry.interactions.push_back(pair);
    }
    cov::MODiagramData group_data;
    group_data.ligand_field_point_group="declared local frame";
    group_data.ligand_field_metal_atom=0;
    group_data.ligand_field_ligand_atoms={1};
    cov::MODiagramLevel group_level;
    group_level.member_indices={0,1};
    group_level.metal_d_weight=0.2;
    group_level.ligand_p_weight=0.8;
    group_level.pi_fraction=1.0;
    auto availability=cov::metal_ligand_detail_availability(group_source,group_data,group_level);
    require(availability.scope==cov::ChemistryStatus::Determined && availability.populations &&
            availability.overlap && availability.channels && group_level.metal_ligand_overlap==0.0,
            "measured group zero was confused with missing evidence");
    auto incomplete_source=group_source;
    incomplete_source.orbitals[1].chemistry.interactions.clear();
    availability=cov::metal_ligand_detail_availability(incomplete_source,group_data,group_level);
    require(availability.populations && !availability.overlap && !availability.channels,
            "a missing member interaction was presented as measured zero");
    incomplete_source.orbitals[1].chemistry.available=false;
    availability=cov::metal_ligand_detail_availability(incomplete_source,group_data,group_level);
    require(!availability.populations && !availability.overlap && !availability.channels,
            "a partially unavailable group was presented as complete");
    auto invalid_level=group_level;
    invalid_level.metal_ligand_overlap=std::numeric_limits<double>::quiet_NaN();
    require(!cov::metal_ligand_detail_availability(group_source,group_data,invalid_level).overlap,
            "a non-finite overlap was presented as available");
    auto missing_scope=group_data;
    missing_scope.ligand_field_ligand_atoms.clear();
    require(cov::metal_ligand_detail_availability(group_source,missing_scope,group_level).scope==
                cov::ChemistryStatus::Unavailable,"missing local context was declared not applicable");
    auto organic_source=group_source;
    organic_source.atoms[0].atomic_number=6;
    availability=cov::metal_ligand_detail_availability(organic_source,group_data,group_level);
    require(availability.scope==cov::ChemistryStatus::NotApplicable && !availability.populations &&
            !availability.overlap && !availability.channels,
            "organic group displayed metal-ligand defaults");
    auto atomic_source=group_source;
    atomic_source.atoms.resize(1);
    require(cov::metal_ligand_detail_availability(atomic_source,group_data,group_level).scope==
                cov::ChemistryStatus::NotApplicable,"isolated atom acquired a metal-ligand scope");

    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::CoordinationGeometry,
                Language::ChineseSimplified)) == "配位几何",
            "Simplified Chinese coordination label mismatch");
    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::PiWeakNearNonbonding,
                Language::Japanese)).find("弱場") != std::string_view::npos,
            "Japanese weak-field label is not localised");
    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::GaussianEnrichmentAttached,
                Language::French)).find("Enrichissement") != std::string_view::npos,
            "French enrichment label is not localised");

    require(std::string_view(cov::ui::localised_wavefunction_source(
                WavefunctionSource::Fchk, Language::ChineseSimplified)) == "FCHK",
            "FCHK scientific identifier changed");
    require(std::string_view(cov::ui::localised_data_provenance(
                DataProvenance::Derived, Language::ChineseSimplified)) == "推导数据",
            "provenance value is not localised");
    require(std::string_view(cov::ui::localised_annotation_source(
                AnnotationSource::ParsedLabel, Language::Japanese)) == "ラベル解析",
            "annotation source is not localised");
    require(std::string_view(cov::ui::localised_bonding_class(
                BondingClass::Antibonding, Language::French)) == "antiliante",
            "bonding class is not localised");
    require(std::string_view(cov::ui::localised_orbital_bonding_role(
                OrbitalBondingRole::Nonbonding,
                Language::ChineseSimplified)) == "非键",
            "orbital bonding role is not localised");
    require(std::string_view(cov::ui::localised_pi_interaction_kind(
                PiInteractionKind::Acceptor,
                Language::Japanese)).find("逆供与") != std::string_view::npos,
            "pi interaction is not localised");

    require(cov::ui::localised_geometry_name(
                "OC-6", "Octahedral", Language::ChineseSimplified) == "八面体形",
            "OC-6 geometry is not localised");
    require(cov::ui::localised_geometry_name(
                "SPC-10", "Sphenocorona", Language::French) == "sphénocouronne",
            "SPC-10 geometry is not localised");
    require(cov::ui::localised_geometry_name(
                "future-id", "Future geometry", Language::Japanese) == "Future geometry",
            "unknown geometry fallback changed");

    const std::string method =
        "COV FCHK S-metric minimal atomic-reference projection";
    require_contains(cov::ui::localised_chemistry_method(
                         method, Language::ChineseSimplified),
                     "S 度量", "Chinese chemistry method");
    require_contains(cov::ui::localised_chemistry_note(
                         "No stable atom-pair interaction frame; chemistry remains UND",
                         Language::French),
                     "indéterminée", "French chemistry note");
    require(cov::ui::localised_chemistry_note(
                "future backend note", Language::ChineseSimplified) ==
                "future backend note",
            "unknown backend note fallback changed");

    MODiagramData data;
    data.mode = MODiagramMode::ValenceCentral;
    data.metadata.resize(12u);
    data.levels.resize(5u);
    data.levels[2].raw_data_fallback = true;
    data.selection.valence_occupied_count = 3u;
    data.selection.frontier_virtual_count = 2u;
    data.selection.hidden_count = 7u;
    data.selection.protected_overflow_count = 1u;
    data.ligand_field_point_group = "Oh";
    data.ligand_field_geometry_id = "OC-6";
    data.ligand_field_coordination_number = 6u;
    data.pi_interactions.resize(2u);
    data.spin_counterpart_pair_count = 4u;
    data.spin_counterpart_unmatched_visible = 1u;

    const std::string zh_summary =
        cov::ui::localised_diagram_selection_summary(
            data, Language::ChineseSimplified);
    const std::string ja_summary =
        cov::ui::localised_diagram_selection_summary(data, Language::Japanese);
    const std::string fr_summary =
        cov::ui::localised_diagram_selection_summary(data, Language::French);
    require_contains(zh_summary, "MO 图摘要", "Chinese diagram summary");
    require_contains(zh_summary, "Oh", "Chinese diagram summary");
    require_contains(ja_summary, "表示準位", "Japanese diagram summary");
    require_contains(fr_summary, "niveaux visibles", "French diagram summary");

    require_contains(cov::ui::orbital_ui_glyph_seed(
                         Language::ChineseSimplified),
                     "斯芬诺冠形", "Chinese glyph seed");
    require_contains(cov::ui::orbital_ui_glyph_seed(Language::Japanese),
                     "逆供与", "Japanese glyph seed");
    require_contains(cov::ui::orbital_ui_glyph_seed(Language::French),
                     "tétradécaédrique", "French glyph seed");

    std::cout << "ui_orbital_localisation_smoke: four-language orbital UI catalogue verified\n";
    return 0;
}
