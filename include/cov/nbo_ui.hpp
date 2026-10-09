#pragma once

#include "cov/nbo.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/ui.hpp"
#include <array>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <set>
#include <vector>

namespace cov::ui {

struct NboFocusAtomGroup {
    std::size_t id = 0;
    std::set<std::size_t> atoms; // explicit producer one-based atom IDs
};

// Inspection controls only. These never change the central canonical MO diagram.
struct NboFocusUIState {
    std::optional<std::size_t> canonical_index;
    std::set<std::size_t> visible_atoms; // NBO producer one-based atom IDs
    std::set<std::size_t> ligand_atoms;  // unsaved group editor only
    std::vector<NboFocusAtomGroup> ligand_groups;
    std::size_t next_group_id = 1;
    std::set<std::string> hidden_shells;
    std::set<std::string> explicitly_included_shells;
    bool group_by_atom = false;
    bool group_ligands_by_l = false;
    bool show_core = false;
    bool show_rydberg = false;
    bool show_full_details = false;
    std::string diagram_id;
    std::optional<std::size_t> last_inspected;
    std::optional<std::size_t> pending_canonical_selection;
};

struct NboUIState {
    std::array<char, 2048> path{};
    std::array<char, 2048> archive47{};
    std::array<char, 2048> aonbo{};
    std::array<char, 2048> nbomo{};
    std::array<char, 2048> naomo{};
    std::array<char, 2048> aonao{};
    std::array<char, 2048> naonbo{};
    std::array<char, 2048> export_path{};
    std::optional<NboDataset> dataset;
    const NboIntegration* integration = nullptr; // main-owned, clear on reload
    const RoutedAnalysis* routed = nullptr; // main-owned, rebuilt on reload/reattach
    std::optional<NboInputDiscovery> input_discovery;
    std::optional<std::size_t> pending_candidate;
    std::string input_status;
    bool show_advanced_inputs = false;
    NboAomoUIState aomo;
    std::set<std::size_t> selected_atoms; // zero-based canonical atoms
    std::optional<std::size_t> selected_structure; // integration.structure index
    std::optional<NboOrbitalRef> inspected_nlmo;
    std::optional<NboOrbitalRef> inspected_nho;
    std::optional<NboOrbitalRef> nho_sum_owner;
    std::set<std::size_t> nho_sum_nao_indices; // selected zero-based NAO columns
    std::string inspected_dataset_id;
    int atom_colour_mode = 0; // 0 element, 1 NPA charge, 2 spin density
    bool show_bond_indices = false;
    bool show_e2 = false;
    bool show_lewis_skeleton = false;
    std::string error;
    std::string export_status;
    std::size_t selected_orbital = std::numeric_limits<std::size_t>::max();
    int analysis_segment = -1;
    double contribution_threshold = 0.01;
    NboFocusUIState focus;
};

struct NboUIActions {
    bool choose_input = false;
    bool attach = false;
    bool canonical_set = false;
    bool nbo_set = false;
    bool export_bundle = false;
    std::optional<std::size_t> selected_orbital;
};

NboUIActions draw_nbo_panel(NboUIState& state, Language language,
                            bool canonical_loaded, bool nbo_renderable,
                            ActiveOrbitalKind active_kind, float scale,
                            const Wavefunction* canonical, std::size_t canonical_index,
                            const MODiagramViewSnapshot* diagram);
std::string nbo_glyph_seed(Language language);
std::string nbo_structure_display_label(const NboStructureEvidence& record,
                                        const Wavefunction* canonical,
                                        Language language);

void draw_nbo_focus_view(NboFocusUIState& focus, const NboDataset& dataset,
                         const Wavefunction& canonical,
                         const MODiagramViewSnapshot& diagram,
                         Language language, float scale);

// Exports one immutable dataset snapshot and a matched SVG/PNG report view.
// The report uses occupation only; diagonal Fock values retain their own label.
void export_nbo_bundle(const NboDataset& dataset, std::size_t selected_orbital,
                       const std::filesystem::path& base,
                       const Wavefunction* canonical, std::size_t canonical_index,
                       double contribution_threshold,
                       const ActiveOrbitalView& active_view,
                       const NboIntegration* integration,
                       const MODiagramViewSnapshot* diagram,
                       const NboFocusUIState* focus,
                       Language language = Language::English);

} // namespace cov::ui
