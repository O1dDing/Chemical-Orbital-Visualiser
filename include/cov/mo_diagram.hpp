#pragma once
#include "cov/orbital_group_bonding.hpp"

#include "cov/model.hpp"
#include "cov/orbital_view.hpp"
#include "cov/pi_pair_evidence.hpp"
#include "cov/mo_sigma_framework.hpp"
#include "cov/nbo_spin_average.hpp"
#include "cov/pi_field_response.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace cov {
struct RoutedAnalysis;

enum class MODiagramMode {
    ValenceCentral = 0,
    DelocalisedPiFamilyOnly,
    MulticentreActiveSpaceOnly,
};

// Stable machine and human labels shared by the live view and every export
// format.  The built diagram is authoritative because a requested compact
// mode may legitimately fall back to the valence view when no active space is
// supported by the wavefunction.
[[nodiscard]] const char* mo_diagram_mode_name(MODiagramMode mode) noexcept;
[[nodiscard]] const char* mo_diagram_mode_title(MODiagramMode mode) noexcept;

// Compact main-group conjugated systems are clearer as their complete
// delocalised-pi active space.  Transition-metal ligand-field cases retain the
// valence-central diagram even if a ligand happens to contain a pi family.
[[nodiscard]] MODiagramMode preferred_compact_mo_diagram_mode(
    const Wavefunction& wavefunction,
    bool compact) noexcept;

enum class EnergyAxisMode {
    Linear = 0,
    NonlinearFocus,
};

struct EnergyAxisKnot {
    double energy_hartree = 0.0;
    double coordinate = 0.0;
};

struct EnergyTransform {
    EnergyAxisMode mode = EnergyAxisMode::NonlinearFocus;
    double focus_hartree = 0.0;
    double scale_hartree = 1.0e-4;
    double minimum_gap_weight = 0.070;
    std::vector<EnergyAxisKnot> knots;
};

[[nodiscard]] const char* energy_axis_mode_name(EnergyAxisMode mode) noexcept;
[[nodiscard]] const char* energy_transform_name(EnergyAxisMode mode) noexcept;
[[nodiscard]] EnergyTransform build_energy_transform(
    const std::vector<double>& energies_hartree,
    EnergyAxisMode mode,
    double minimum_gap_weight = 0.070);
[[nodiscard]] double energy_display_coordinate(
    double energy_hartree,
    const EnergyTransform& transform) noexcept;
[[nodiscard]] double energy_from_display_coordinate(
    double coordinate,
    const EnergyTransform& transform) noexcept;

enum class AnnotationSource {
    Direct,
    ParsedLabel,
    Derived,
    Heuristic,
    Unavailable,
};

enum class BondingClass {
    Unclassified,
    Bonding,
    Nonbonding,
    Antibonding,
    Mixed,
};

struct MulticentreDescriptor {
    bool available = false;
    std::size_t centres = 0;
    double electrons = 0.0;
    std::vector<std::size_t> atom_indices;
    std::string label;
    std::size_t channel_count = 0;
    std::string source_subspace_id;
    double source_subspace_electron_count = 0.0;
    AnnotationSource source = AnnotationSource::Unavailable;
    double confidence = 0.0;
    bool heuristic = false;
};

struct DelocalisedPiDescriptor {
    bool available = false;
    std::size_t participating_atoms = 0;
    double participating_electrons = 0.0;
    std::vector<std::size_t> atom_indices;
    std::vector<std::size_t> orbital_indices;
    std::string family_id;
    // `available` says that a delocalised-pi family is known.  Topology is a
    // separate claim: parsed producer labels can establish the family without
    // establishing channel directions or cyclicity.
    bool topology_available = false;
    DelocalisedPiTopology topology = DelocalisedPiTopology::Unknown;
    std::size_t orientation_channels = 0;
    bool cyclic_topology = false;
    std::vector<PiOrientationChannel> orientation_channel_details;
    PiTopologyGraphEvidence topology_graph;
    std::string label;
    AnnotationSource source = AnnotationSource::Unavailable;
    double confidence = 0.0;
    bool heuristic = false;
};

// The same compact topology suffix is used by the live diagram and bitmap/
// vector exports. An empty suffix means that the family is known but a
// concrete topology is not.
[[nodiscard]] inline const char* compact_pi_topology_code(
    const DelocalisedPiTopology topology) noexcept {
    switch (topology) {
        case DelocalisedPiTopology::Path: return "P";
        case DelocalisedPiTopology::Cycle: return "C";
        case DelocalisedPiTopology::BranchedResonance: return "B";
        case DelocalisedPiTopology::Spiro: return "S";
        case DelocalisedPiTopology::HapticMetal: return "H";
        case DelocalisedPiTopology::SymmetryDirectSum: return "D";
        case DelocalisedPiTopology::MultiChannel: return "M";
        default: return "?";
    }
}
[[nodiscard]] inline std::string compact_pi_topology_suffix(
    const DelocalisedPiDescriptor& descriptor) {
    if (!descriptor.topology_available) return {};
    return std::to_string(descriptor.orientation_channels)+"ch "+
           compact_pi_topology_code(descriptor.topology);
}

struct OrbitalAnnotation {
    // Machine-friendly canonical family names: sigma / pi / delta / phi / unavailable.
    std::string family = "unavailable";
    BondingClass bonding_class = BondingClass::Unclassified;
    AnnotationSource family_source = AnnotationSource::Unavailable;
    AnnotationSource bonding_source = AnnotationSource::Unavailable;
    double family_confidence = 0.0;
    double bonding_confidence = 0.0;
    MulticentreDescriptor multicentre;
    DelocalisedPiDescriptor delocalised_pi;
    bool heuristic = false;
};

struct SymmetryNotation {
    std::string base;
    std::string subscript;
    std::string superscript;
    std::string raw;
};

[[nodiscard]] SymmetryNotation parse_symmetry_notation(std::string_view raw);
[[nodiscard]] std::string format_symmetry_unicode(std::string_view raw);
// Display only: preserve machine group keys and never convert group irreps.
[[nodiscard]] std::string point_group_display(std::string_view raw);

struct DiagramSelectionPlan {
    std::vector<std::size_t> included_indices;
    std::size_t hidden_count = 0;
    std::size_t valence_occupied_count = 0;
    std::size_t frontier_virtual_count = 0;
    std::size_t final_group_count=0, final_member_count=0;
    bool counts_are_final=false;
    // Complete protected manifolds may legitimately exceed the visual row
    // target. This records only that unavoidable excess; it is never a
    // licence to keep arbitrary unprotected rows.
    std::size_t protected_overflow_count = 0;
    std::string summary;
};

struct MODiagramOptions {
    // Presentation supplied by the viewer; no dependency on its font library.
    std::vector<std::string> display_names;
    // Full canonical MO labels supplied by the same verified naming table.
    // Separate from row/local symmetry_view and immutable producer labels.
    std::vector<std::string> display_irreps;
    std::vector<std::string> display_point_groups;
    std::string figure_title;
    std::string axis_title;
    std::function<void(std::vector<std::uint8_t>&,int,int,int,int,const std::string&,
                       std::uint8_t,std::uint8_t,std::uint8_t,int)> raster_text;
    std::function<int(const std::string&,int)> raster_text_width;
    const RoutedAnalysis* routed = nullptr; // immutable, same canonical fingerprint
    const NboIntegration* nbo_source=nullptr; // immutable same-source sigma projector
    const NboSalcModel* source_salc_model=nullptr; // verified unaveraged source
    const NboRoCommonEnergyModel* ro_common_energy=nullptr; // immutable prepared cache
    const PiFieldResponseAnalysis* pi_field_response=nullptr; // same wavefunction/attachment
    bool use_ro_common_energy=true;
    // Internal view override. Source wavefunction and its fingerprint never change.
    std::vector<double> canonical_display_energies;
    std::string routed_identity; // cache generation; changes on every reattachment
    MODiagramMode mode = MODiagramMode::ValenceCentral;
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    EnergyAxisMode energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    std::size_t selected_index = 0;

    std::size_t neighbourhood = 12;
    std::size_t max_levels = 0;
    std::size_t max_virtual_levels = 10;
    bool hide_ligand_centred_intermediates = false;
    // 0 preserves the non-NBO path; 1/2/3 are brief/research/full AOMO scopes.
    // These describe display membership, never a different electronic state.
    unsigned aomo_scope=0;
    bool show_core_background=false, show_fragment_background=false;
    std::vector<std::size_t> display_centre_atoms;

    double nonlinear_minimum_gap_weight = 0.070;

    // A small same-symmetry pi splitting is shown as an approximately
    // nonbonding weak-coupling level instead of forcing a donor/acceptor
    // assignment.  The retained member is the one with the larger metal
    // valence contribution, so the reduced ligand-field diagram keeps the
    // chemically relevant d-level rather than an arbitrary canonical MO.
    double weak_pi_split_hartree = 0.020;
    double weak_metal_ligand_overlap = 0.025;

    int width = 1200;
    int height = 900;
    bool include_hidden_in_metadata = true;
    double weak_crystal_field_split_hartree = 0.020;
    double weak_crystal_field_overlap = 0.025;
};

enum class PiInteractionKind {
    Donor,
    Acceptor,
    Coupled,
    WeakNearNonbonding,
};

enum class OrbitalEnergyGapKind { PiPartner, CrystalField };

struct OrbitalEnergyGapDescriptor {
    std::size_t lower_level = 0;
    std::size_t upper_level = 0;
    std::vector<std::size_t> lower_orbitals;
    std::vector<std::size_t> upper_orbitals;
    std::string symmetry;
    PiInteractionKind kind = PiInteractionKind::Coupled;
    double splitting_hartree = 0.0;
    double confidence = 0.0;
    // Reading priority only: minimum absolute cross-Fock trace of the two
    // complete uniform-sign endpoint groups, in hartree; not a bond energy.
    double display_strength_hartree = std::numeric_limits<double>::quiet_NaN();
    bool lower_visible = true;
    bool upper_visible = true;
    std::size_t retained_level = 0;
    std::shared_ptr<const PiPartnerAssessment> orbital_evidence;
    // Equivalent physical projector scopes share one displayed relation;
    // retain each original channel reference for reproducible diagnostics.
    std::vector<std::string> equivalent_channel_ids;
    OrbitalEnergyGapKind gap_kind = OrbitalEnergyGapKind::PiPartner;
    OrbitalSymmetryExplanation lower_symmetry_scope;
    OrbitalSymmetryExplanation upper_symmetry_scope;
    double lower_energy_hartree = std::numeric_limits<double>::quiet_NaN();
    double upper_energy_hartree = std::numeric_limits<double>::quiet_NaN();
    double lower_energy_spread_hartree = std::numeric_limits<double>::quiet_NaN();
    double upper_energy_spread_hartree = std::numeric_limits<double>::quiet_NaN();
    double weak_split_threshold_hartree = std::numeric_limits<double>::quiet_NaN();
    double weak_overlap_threshold = std::numeric_limits<double>::quiet_NaN();
    std::shared_ptr<const WeakCrystalFieldAssessment> crystal_field_evidence;
};

// Source compatibility for callers that construct actual pi partner records.
using PiInteractionDescriptor = OrbitalEnergyGapDescriptor;
[[nodiscard]] const char* orbital_energy_gap_kind_name(OrbitalEnergyGapKind kind) noexcept;
[[nodiscard]] std::string orbital_energy_gap_json(
    const OrbitalEnergyGapDescriptor& gap, EnergyUnit unit = EnergyUnit::Hartree);
[[nodiscard]] std::string orbital_energy_gap_array_json(
    const std::vector<OrbitalEnergyGapDescriptor>& gaps, EnergyUnit unit = EnergyUnit::Hartree);

[[nodiscard]] const char* pi_interaction_kind_name(
    PiInteractionKind kind) noexcept;

// Exclusive buckets in the original complete NAO norm, averaged over the
// actual canonical members. Missing data is never interpreted as a zero.
struct MOGroupCompositionLedger {
    bool available=false, complete=false;
    std::string status="unavailable", detail, source="unavailable";
    std::size_t member_count=0;
    double weight_sum=0, normalization_error=0;
    double centre_current_s=0, centre_current_p=0, centre_current_d=0,
           centre_current_f=0, centre_other=0;
    double ligand_valence=0, ligand_other=0, other_atoms=0, core=0, unresolved=0;
    // Subtotals of ligand_valence, not additional exclusive buckets.
    double ligand_valence_s=0, ligand_valence_p=0;
};
struct MOCurrentRadialShell {
    std::size_t atom=0; // canonical zero-based identity
    int n=0, l=0;
    std::string evidence;
};
// Chemical scope is independent of whether a numerical NAO partition exists.
// A user-selected nonmetal centre never becomes a metal, and disconnected
// counterions never become ligands merely by being outside the centre bucket.
struct MOCompositionScope {
    bool applicable=false;
    std::vector<std::size_t> centre_atoms,ligand_atoms,other_atoms;
    std::string detail;
};
[[nodiscard]] MOCompositionScope mo_composition_scope(
    const Wavefunction&,const RoutedAnalysis*,const std::vector<std::size_t>& centres);
struct MOGroupDisplayDecision {
    bool included=false, energy_window=false, major_relation=false, frontier=false;
    double coverage=0; // complete-norm current-centre coverage, not conditional
    double sigma_coverage=0; // overlapping verified donor projector, full norm
    std::vector<std::string> reason_codes;
};
struct MODiagramGroupAudit {
    std::vector<std::size_t> member_indices, member_spin_counterparts;
    double energy_hartree=0, total_occupation=0;
    MOGroupCompositionLedger composition;
    MOGroupDisplayDecision display_decision;
};
// Current shells are frozen from source shell identities before MO selection.
[[nodiscard]] std::vector<MOCurrentRadialShell> mo_current_radial_shells(
    const Wavefunction&, const RoutedAnalysis&, const std::vector<std::size_t>& centres);
[[nodiscard]] MOGroupCompositionLedger mo_group_composition_ledger(
    const Wavefunction&, const RoutedAnalysis*, const std::vector<std::size_t>& members,
    const std::vector<std::size_t>& centres, const std::vector<MOCurrentRadialShell>& shells);

struct MODiagramLevel {
    OrbitalMetadata metadata;
    OrbitalAnnotation annotation;
    OrbitalChemistry chemistry;
    ElectronGlyphs electrons;
    bool homo = false;
    bool lumo = false;
    double layout_energy_hartree = 0.0;

    // One row represents a complete degenerate subspace.  Group-level
    // quantities are averaged per member and therefore remain invariant to a
    // rotation of canonical MOs inside an exactly degenerate subspace.
    std::vector<std::size_t> member_indices;
    std::vector<ElectronGlyphs> member_electrons;
    // For a UDFT spatial-row view, each majority-spin member can retain the
    // matched minority-spin canonical MO.  This keeps selection highlighting
    // and exact member metadata available without drawing a duplicate row.
    std::vector<std::size_t> member_spin_counterparts;
    double energy_spread_hartree = 0.0;
    double total_occupation = 0.0;
    double metal_s_weight = 0.0;
    double metal_p_weight = 0.0;
    double metal_d_weight = 0.0;
    // p population on the atoms directly coordinated to the selected metal.
    // ligand_p_weight may additionally include the rest of those ligand
    // fragments (for example O in CO or N in CN).
    double direct_ligand_p_weight = 0.0;
    double ligand_p_weight = 0.0;
    double sigma_fraction = 0.0;
    double pi_fraction = 0.0;
    double metal_ligand_overlap = 0.0;
    bool raw_data_fallback = false;
    bool approximate_nonbonding = false;
    std::vector<OrbitalGroupBondingResult> bonding_scopes;
    MOGroupCompositionLedger composition;
    MOGroupDisplayDecision display_decision;
};

// One visible row may represent an exactly-degenerate canonical-MO set and,
// after UDFT spatial collapse, its matched opposite-spin counterparts.
[[nodiscard]] inline bool mo_diagram_level_covers_orbital(
    const MODiagramLevel& level,
    const std::size_t orbital_index) noexcept {
    if (std::find(level.member_indices.begin(),level.member_indices.end(),
                  orbital_index)!=level.member_indices.end()) {
        return true;
    }
    if (std::find(level.member_spin_counterparts.begin(),
                  level.member_spin_counterparts.end(),orbital_index)!=
            level.member_spin_counterparts.end()) {
        return true;
    }
    return level.member_indices.empty() &&
           level.metadata.orbital_index==orbital_index;
}

struct LocalGeometryDiagramDescriptor {
    std::size_t centre_atom = 0;
    std::string geometry_id;
    std::string geometry_name;
    std::string point_group;
    std::vector<std::size_t> neighbour_atoms;
    double confidence = 0.0;
    double angular_rms = 0.0;
    double shape_measure = 0.0;
    double radial_cv = 0.0;
};

struct ElectronicStateDiagramMetadata {
    WavefunctionSource source = WavefunctionSource::Unknown;
    std::int32_t charge = 0;
    std::uint32_t multiplicity = 0;
    std::uint32_t alpha_electrons = 0;
    std::uint32_t beta_electrons = 0;
    DataProvenance charge_provenance = DataProvenance::Unavailable;
    DataProvenance multiplicity_provenance = DataProvenance::Unavailable;
    DataProvenance electron_counts_provenance = DataProvenance::Unavailable;
    ScfConvergenceStatus scf_convergence = ScfConvergenceStatus::Unavailable;
    DataProvenance scf_convergence_provenance = DataProvenance::Unavailable;
    WavefunctionStabilityStatus stability =
        WavefunctionStabilityStatus::Unavailable;
    DataProvenance stability_provenance = DataProvenance::Unavailable;
    std::string stability_detail;
    double spin_squared_before = 0.0;
    double spin_squared_after = 0.0;
    DataProvenance spin_squared_provenance = DataProvenance::Unavailable;
    std::vector<double> atomic_partial_charges;
    std::string atomic_partial_charge_scheme;
    DataProvenance atomic_partial_charge_provenance =
        DataProvenance::Unavailable;
    std::string point_group_detected;
    std::string point_group_used;
    DataProvenance point_group_provenance = DataProvenance::Unavailable;
    DataProvenance point_group_detected_provenance = DataProvenance::Unavailable;
    DataProvenance point_group_used_provenance = DataProvenance::Unavailable;
    std::string source_title;
    std::string source_route;
    std::string enrichment_source;
};

struct MODiagramViewContext {
    // Correlates one immutable view with its exports; artifact hashes identify
    // the bytes separately. The selection anchor controls row construction,
    // while inspection can select a hidden MO or an opposite-spin counterpart.
    std::string id;
    std::string origin;
    std::size_t selection_anchor = 0;
    std::optional<std::size_t> inspected_orbital_index;
};

struct MODiagramData {
    std::optional<MODiagramViewContext> view;
    DiagramPlan plan;
    FrontierOrbitals frontier;
    DiagramSelectionPlan selection;
    std::vector<MODiagramLevel> levels;
    std::vector<OrbitalMetadata> metadata;
    std::vector<OrbitalAnnotation> annotations;
    std::vector<PiInteractionDescriptor> pi_interactions;
    // Includes rejected counterparts with their original member identities.
    std::vector<PiPartnerAssessment> pi_partner_candidates;
    // One scientific object per actual mode; canonical groups are nodes, not
    // an arbitrary Cartesian product of purported two-level counterparts.
    std::vector<PiModeNetworkAssessment> pi_mode_networks;
    MOSigmaFramework sigma_framework;
    NboRoCommonEnergyModel ro_common_energy;
    bool using_ro_common_energy=false;
    PiFieldResponseAnalysis pi_field_response;
    MOCompositionScope composition_scope;
    std::vector<MOCurrentRadialShell> current_radial_shells;
    // Complete source groups, including folded and opposite-spin groups.
    std::vector<MODiagramGroupAudit> group_audit;
    ElectronicStateDiagramMetadata electronic_state;
    // Every unambiguous Mayer-supported CN2--CN10 centre, including main-group
    // and non-metal centres.  This is structural metadata and never creates a
    // transition-metal ligand-field diagram by itself.
    std::vector<LocalGeometryDiagramDescriptor> local_geometries;
    MODiagramMode mode = MODiagramMode::ValenceCentral;
    EnergyTransform energy_transform;

    // Local first-shell ligand-field information is kept separate from the
    // full molecular point group.  For an unrestricted calculation the
    // diagram can use majority-spin spatial representatives while retaining
    // the paired alpha/beta occupation and member metadata.
    std::string ligand_field_point_group;
    std::string ligand_field_geometry_id;
    std::string ligand_field_geometry_name;
    std::size_t ligand_field_coordination_number = 0;
    std::size_t ligand_field_metal_atom = 0;
    std::vector<std::size_t> ligand_field_ligand_atoms;
    double ligand_field_confidence = 0.0;
    double ligand_field_angular_rms = 0.0;
    double ligand_field_shape_measure = 0.0;
    double ligand_field_radial_cv = 0.0;
    bool spin_counterparts_collapsed = false;
    bool spin_counterparts_partial = false;
    std::size_t spin_counterpart_pair_count = 0;
    std::size_t spin_counterpart_unmatched_visible = 0;
    // Different-irrep local d-level gaps never contribute to pi partner counts.
    std::vector<OrbitalEnergyGapDescriptor> crystal_field_gaps;
};

[[nodiscard]] std::vector<const OrbitalEnergyGapDescriptor*> orbital_energy_gaps(
    const MODiagramData& data);
[[nodiscard]] std::string pi_partner_candidates_json(
    const std::vector<PiPartnerAssessment>& candidates);

struct MODiagramViewSnapshot {
    const MODiagramData data;
    const MODiagramOptions options;
};

struct MODiagramMemberView {
    std::size_t orbital_index = 0;
    std::optional<std::size_t> spin_counterpart;
    ElectronGlyphs electrons;
    bool selected = false;
    std::size_t inspected_orbital_index = 0;
};

[[nodiscard]] MODiagramViewSnapshot make_mo_diagram_view_snapshot(
    const MODiagramData& data,
    const MODiagramOptions& options,
    std::optional<std::size_t> inspected_orbital_index,
    std::string origin = "explicit-data");
// A mark always has a real canonical MO identity. A degeneracy count alone
// cannot manufacture adjacent MO identities or extra electron arrows.
[[nodiscard]] std::vector<MODiagramMemberView> mo_diagram_member_views(
    const MODiagramData& data, const MODiagramLevel& level);
[[nodiscard]] std::optional<std::size_t> mo_diagram_row_for_orbital(
    const MODiagramData& data, std::size_t orbital_index) noexcept;

[[nodiscard]] const char* annotation_source_name(AnnotationSource source) noexcept;
[[nodiscard]] const char* bonding_class_name(BondingClass value) noexcept;
[[nodiscard]] OrbitalAnnotation annotate_orbital(const MolecularOrbital& orbital);
[[nodiscard]] DiagramSelectionPlan build_valence_selection_plan(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata);
[[nodiscard]] MODiagramData build_mo_diagram_data(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options);

// Availability of the local metal/first-shell fields shown for a row's
// representative members. A default numeric zero is not availability evidence.
struct MetalLigandDetailAvailability {
    ChemistryStatus scope = ChemistryStatus::Unavailable;
    bool composition = false;
    bool populations = false;
    bool overlap = false;
    bool channels = false;
};
[[nodiscard]] MetalLigandDetailAvailability metal_ligand_detail_availability(
    const Wavefunction& wavefunction, const MODiagramData& data,
    const MODiagramLevel& level);

enum class DiagramExportContent { Images, AnalysisData, All };

struct MODiagramExportResult {
    bool svg = false;
    bool png = false;
    bool json = false;
    bool csv = false;
    std::filesystem::path svg_path;
    std::filesystem::path png_path;
    std::filesystem::path json_path;
    std::filesystem::path csv_path;
    std::string error;
};

[[nodiscard]] MODiagramExportResult export_mo_diagram_bundle(
    const MODiagramViewSnapshot& snapshot,
    const std::filesystem::path& base_path,
    DiagramExportContent content = DiagramExportContent::All);
// Explicit computed-export API for noninteractive callers. The live UI must
// pass its already drawn snapshot through the overload above.
[[nodiscard]] MODiagramExportResult export_mo_diagram_bundle(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::filesystem::path& base_path);
[[nodiscard]] bool write_mo_diagram_svg(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_png(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_json(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_csv(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);

} // namespace cov
