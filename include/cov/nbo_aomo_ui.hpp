#pragma once

#include "cov/nbo_integration.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/ui.hpp"
#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cov::ui {

std::string nbo_aomo_glyph_seed(Language language);

enum class NboAomoPreset { Teaching, Research, Full };
enum class NboAomoLane { Left, Centre, Right };

struct NboAomoFragmentGroup {
    std::size_t id=0;
    std::set<std::size_t> atoms; // zero-based canonical atom identities
    bool suggested=false;
};

struct NboAomoNode {
    std::string id,label,detail,energy_semantics;
    std::string symmetry_irrep,name_detail;
    std::size_t symmetry_ordinal=0;
    std::size_t symmetry_multiplicity=1;
    bool symmetry_name_verified=false;
    BondingClass bonding_class=BondingClass::Unclassified;
    std::string bonding_scope_status;
    std::optional<NboOrbitalRef> orbital;
    std::optional<std::size_t> canonical_index;
    std::optional<std::size_t> fragment_group_id;
    std::optional<std::size_t> salc_index;
    std::optional<NboSpatialSpinInfo> spatial_spin;
    std::string subspace_id;
    std::vector<std::size_t> atoms;
    std::vector<std::size_t> member_canonical_indices;
    std::vector<double> member_energies_hartree,member_display_energies_hartree,member_occupations;
    bool weak_display_container=false;
    std::optional<double> energy_hartree,occupation;
    // Source energies above remain unchanged. Display energies use the selected
    // operator definition, then group means and symmetric shell-stack offsets.
    std::optional<double> display_energy_hartree;
    std::string display_energy_semantics="individual",display_group_id;
    std::string individual_label,shell_label;
    std::string spatial_pair_id,spatial_pair_label;
    std::size_t shell_member_index=0,shell_member_count=1;
    float display_offset_y=0;
    std::optional<double> metric_norm2;
    NboAomoLane lane=NboAomoLane::Left;
    // Quantitative side levels use a verified same-operator expectation value.
    // A missing value belongs in the labelled non-quantitative side band.
    bool quantitative_energy=false;
    float x=0,y=0,width=160,height=24,label_y=0;
    // All renderers consume these frozen text rectangles (top-left + size).
    // They are separate from the true energy line and its click target.
    float label_x=0,label_width=0,label_height=14;
    float occupation_x=0,occupation_y=0,occupation_width=0,occupation_height=14;
    std::string occupation_label;
    bool occupation_on_bar=false;
    bool group_header=false,available=true,composition_available=true;
};
struct NboAomoEdge {
    std::string id;
    std::size_t source_node=0,target_node=0;
    double coefficient=0;
    std::optional<double> weight;
    std::optional<double> projection_strength_nonadditive;
    NboSource source;
    std::optional<std::size_t> salc_link_index;
    bool visible=false; // all verified links remain in JSON/CSV and numeric detail
};
struct NboAomoEnergyTick { double energy_hartree=0; float y=0; };
struct NboAomoCaption {
    std::string role,text;
    float x=0,y=0,width=0,height=0;
};
struct NboAomoViewSnapshot {
    NboRoCommonEnergyModel ro_common_energy;
    bool using_ro_common_energy=false;
    std::string display_energy_definition="source_canonical_energy";
    PiFieldResponseAnalysis pi_field_response;
    MOSigmaFramework sigma_framework;
    std::vector<MOCurrentRadialShell> current_radial_shells;
    std::vector<MODiagramGroupAudit> group_audit;
    DiagramSelectionPlan final_selection;
    Language language=Language::English;
    std::string id,integration_id,mo_snapshot_id,capability_status,capability_detail;
    std::string mo_energy_axis_mode,mo_energy_axis_detail;
    std::string display_energy_unit;
    EnergyUnit energy_unit=EnergyUnit::Hartree;
    NboOrbitalKind basis_kind=NboOrbitalKind::NAO;
    std::size_t focused_canonical_index=0;
    std::optional<double> focused_projection_weight,focused_projection_residual_norm;
    std::vector<std::size_t> central_mo_indices;
    std::vector<NboAomoNode> nodes;
    std::vector<NboAomoEdge> edges;
    std::vector<NboAomoEnergyTick> energy_ticks;
    // A single frozen palette drives screen, SVG, PNG and JSON axis labels.
    std::array<std::uint8_t,3> energy_tick_screen_rgb{139,157,178};
    std::array<std::uint8_t,3> energy_tick_export_rgb{139,157,178};
    std::string energy_tick_semantics="linear-neutral";
    std::vector<NboAomoCaption> captions;
    std::vector<PiPartnerAssessment> pi_partner_candidates;
    std::vector<PiModeNetworkAssessment> pi_mode_networks;
    std::vector<PiInteractionDescriptor> pi_interactions;
    std::vector<OrbitalGroupBondingResult> bonding_groups;
    std::optional<NboOrbitalSelection> selection;
    std::optional<ActiveOrbitalView> active_view;
    std::vector<NboAomoFragmentGroup> fragment_groups;
    std::vector<std::string> sum_component_ids;
    std::size_t hidden_basis_count=0,hidden_mo_count=0;
    std::size_t hidden_class_count=0;
    std::size_t hidden_h_count=0;
    NboAomoPreset preset=NboAomoPreset::Teaching;
    bool overview=false;
    bool all_connections=false;
    bool illustrative_side_layout=false;
    bool paper_export=false;
    bool auto_rydberg_expanded=false;
    std::optional<double> focused_display_weight,focused_hidden_weight;
    std::optional<double> focused_hidden_core_weight,focused_hidden_rydberg_weight,
        focused_hidden_valence_weight,focused_hidden_group_weight,focused_hidden_h_weight;
    double numerical_zero_bound=0;
    std::string numerical_zero_reason;
    std::string selected_side_node_id;
    std::shared_ptr<const NboSalcModel> salc_model;
    std::shared_ptr<const NboSalcModel> source_salc_model;
    std::shared_ptr<const NboAomoNames> names;
    std::string name_ordinal_scope;
    std::size_t hidden_numeric_zero_count=0,hidden_readability_count=0;
    std::size_t hidden_group_count=0;
    float canvas_width=0,canvas_height=0;
    float qualitative_band_y=0;
    EnergyTransform energy_transform;
    double axis_coordinate_min=0,axis_coordinate_max=1;
    float numeric_top=0,numeric_span=0;
    float label_font_size=14;
    float orbital_bar_width=62,orbital_bar_stroke_width=1.5f;
    float footer_y=0;
    std::array<float,3> lane_x{},lane_width{};
    bool show_core=false,show_rydberg=false,hide_h_orbitals=false;
    bool show_fragment_background=false;
    bool show_atom_numbers=true,show_fragment_numbers=true,number_ignore_h=false;
    float zoom=1,pan_x=0,pan_y=0;
};

struct NboAomoUIState {
    Language language=Language::English;
    bool initial_fit_pending=true;
    float initial_fit_viewport_width=0;
    int initial_fit_stable_frames=0;
    NboOrbitalKind basis_kind=NboOrbitalKind::NAO;
    float zoom=1,pan_x=0,pan_y=0;
    std::set<std::size_t> collapsed_atoms; // zero-based canonical
    std::set<std::size_t> expanded_atoms;
    bool show_core=false,show_rydberg=false,hide_h_orbitals=false;
    bool show_fragment_background=false;
    bool show_atom_numbers=true,show_fragment_numbers=true,number_ignore_h=false;
    NboAomoPreset preset=NboAomoPreset::Teaching;
    bool overview=true;
    bool all_connections=false;
    std::size_t last_drawn_dash_segments=0,last_unclipped_dash_segments=0;
    double last_connection_draw_ms=0;
    bool illustrative_side_layout=false;
    bool paper_export=false;
    std::set<std::string> expanded_subspaces,collapsed_subspaces;
    std::set<std::string> expanded_weak_groups;
    std::set<std::string> expanded_fragments;
    std::set<std::size_t> expanded_user_fragments;
    std::string selected_side_node_id;
    std::set<std::size_t> draft_fragment_atoms;
    std::optional<std::size_t> editing_fragment_id;
    std::set<std::string> sum_component_ids;
    std::optional<std::size_t> sum_canonical_index;
    NboOrbitalKind sum_basis_kind=NboOrbitalKind::NAO;
    std::vector<std::string> ambiguous_edge_ids;
    std::vector<std::string> ambiguous_node_ids;
    std::vector<NboAomoFragmentGroup> fragment_groups;
    bool suggested_fragments_initialized=false;
    std::size_t next_fragment_id=1;
    std::optional<std::size_t> focused_canonical_index;
    std::optional<std::size_t> last_inspected;
    std::optional<NboOrbitalSelection> pending_selection;
    std::optional<NboOrbitalSelection> selection; // last applied 3D identity; root updates
    std::optional<ActiveOrbitalView> active_view; // same presentation identity as details/export
    bool export_requested=false;
    DiagramExportContent export_content=DiagramExportContent::Images;
    bool show_full_numeric=false;
    std::string source_id,status,export_status;
    std::array<char,2048> export_path{};
    std::shared_ptr<const NboAomoViewSnapshot> drawn_snapshot;
    std::shared_ptr<const NboSalcModel> salc_model;
    std::shared_ptr<const NboSalcModel> source_salc_model,spin_averaged_salc_model;
    std::shared_ptr<const NboRoCommonEnergyModel> common_energy_model;
    std::shared_ptr<const NboAomoNames> source_names,spin_averaged_names;
    std::uint64_t revision=0;
    std::shared_ptr<const NboAomoNames> names;
    const NboSalcModel* names_model=nullptr;
    std::shared_ptr<const NboAomoNames> filtered_names,filtered_names_source;
    std::vector<std::size_t> filtered_canonical_indices,filtered_salc_indices;
    std::vector<double> filtered_canonical_display_energies;
    std::string filtered_name_scope;
};

// Returns true only when the whole AO/NAO--MO view is scientifically available.
// Otherwise the caller should draw its existing Gaussian-only MO diagram.
void apply_nbo_aomo_preset(NboAomoUIState& state,NboAomoPreset preset);
bool prepare_nbo_aomo_state(NboAomoUIState& state,const NboIntegration& integration,
                           const Wavefunction& canonical);
bool draw_nbo_aomo_diagram(NboAomoUIState& state,const NboIntegration& integration,
                           const Wavefunction& canonical,const MODiagramViewSnapshot& diagram,
                           Language language,float ui_scale);

struct NboAomoExportResult {
    bool json=false,csv=false,svg=false,png=false;
    std::filesystem::path json_path,csv_path,svg_path,png_path;
    std::string error;
};
NboAomoExportResult export_nbo_aomo_bundle(const NboAomoViewSnapshot& snapshot,
    const NboIntegration& integration,const std::filesystem::path& base,
    DiagramExportContent content=DiagramExportContent::All);

// Picks an occupied, non-core NBO associated with the current central MO when
// matrix evidence exists. No selection is made when that association is absent.
std::optional<std::size_t> suggest_initial_nbo_index(const NboIntegration& integration,
    const MODiagramViewSnapshot* diagram);

} // namespace cov::ui
