#pragma once

#include "cov/nbo_salc.hpp"
#include <string>
#include <vector>
#include <memory>
#include <limits>

namespace cov::ui {
struct NboIrrepContent {
    std::string irrep;
    std::size_t dimension=0,multiplicity=0;
};
struct NboIrrepComponent {
    std::string irrep;
    std::size_t dimension=0;
    // Squared S norm divided by the unchanged source column's squared S norm.
    // This is symmetry composition, not an occupation or a probability.
    double weight=0;
    // Unnormalised projected component in component_source_members order. Summing
    // the reconstructed components recovers the unchanged source column.
    std::vector<double> source_coefficients;
};
struct NboAomoName {
    std::string label, irrep;
    // Group actually used for the verified character comparison; never a
    // dimension/energy guess or a group inferred from the label spelling.
    std::string point_group;
    // Zero means that a complete, unambiguous ordinal could not be established.
    std::size_t ordinal=0;
    std::size_t view_row_ordinal=0, visible_partner_count=0;
    bool approximate_dominant_label=false;
    std::string dominant_irrep;
    double dominant_weight=0;
    // Irrep evidence only: true with ordinal==0 means the symmetry is known
    // while occurrence order / repeated-copy membership remains unresolved.
    bool verified=false;
    // Number of irreducible copies in the measured containing span.
    // This does not assign an individual copy or occurrence ordinal.
    std::size_t representation_multiplicity=1;
    std::string detail;
    // A verified single irrep occurrence, independent of whether its position
    // among other occurrences can be numbered. Empty means unproved partners.
    std::string partner_block_id;
    std::size_t partner_block_size=0;
    // Failed source-copy and ordering checks are distinct from missing input.
    // These diagnostics never promote a projection to a verified partner block.
    std::string partner_status, partner_failure_quantity;
    std::optional<double> partner_failure_value, partner_failure_limit;
    // Scientific containing span is independent of the current display filter.
    std::vector<std::size_t> containing_members;
    std::vector<NboIrrepContent> containing_irreps;
    std::string status,ordinal_scope;
    std::size_t complete_set_ordinal=0;
    std::string ordinal_status,complete_set_ordinal_status;
    std::vector<std::size_t> ordinal_blocking_members;
    // Conservative possible complete-set ranks, not a chosen physical order.
    // Absent if the relevant irrep counts or energy bounds are unavailable.
    std::optional<std::size_t> complete_set_ordinal_lower,complete_set_ordinal_upper;
    std::optional<double> projection_residual;
    bool decomposition_verified=false;
    std::string decomposition_status;
    // S-norm amplitude and source-normalised scalar errors, respectively.
    double decomposition_reconstruction_residual=std::numeric_limits<double>::quiet_NaN();
    double decomposition_orthogonality_error=std::numeric_limits<double>::quiet_NaN();
    double decomposition_weight_sum_error=std::numeric_limits<double>::quiet_NaN();
    std::vector<NboIrrepComponent> components;
    // Separate from the original containing span: a non-closed local SALC can
    // have a global composition reconstructed in the full canonical basis.
    std::string component_source_kind;
    std::vector<std::size_t> component_source_members;
};
struct NboAomoNames {
    std::vector<NboAomoName> canonical, salc;
};
std::string serialize_orbital_name_json(const NboAomoName&);
std::string serialize_orbital_names_json(const NboAomoNames&);
// Presentation-only evidence cache. Call once per immutable attachment, never
// per frame. No coefficients, source labels, orbital ordering or fields change.
NboAomoNames build_nbo_aomo_names(const Wavefunction&,const NboIntegration&,
                                 const NboSalcModel*);
// Canonical naming is available without an NBO attachment. The cache is for
// immutable loaded wavefunctions; it never changes producer data or coefficients.
std::shared_ptr<const NboAomoNames> canonical_mo_names(const Wavefunction&);
void invalidate_canonical_mo_names_cache();
// Display ordinals count every visible certified occurrence, including partial
// views; the energy ordering key always uses its full certified partner block.
// An explicit display-energy list uses immutable canonical indices and must
// cover the entire block. Missing/nonfinite entries never fall back to a
// different energy definition. Omitting the list retains source-energy order.
// Source IDs, complete-set ordinals, evidence and unresolved copies are retained.
NboAomoNames nbo_aomo_names_for_view(const Wavefunction&,const NboAomoNames&,
    const std::vector<std::size_t>& canonical_indices,
    const std::vector<std::size_t>& salc_indices,const NboSalcModel*,
    const std::string& scope,const std::vector<double>* display_energies=nullptr);
// Source identity is spin-block based; unavailable source indices are explicitly
// identified as list positions. Internal indices remain the selection addresses.
std::string canonical_mo_source_label(const Wavefunction&,std::size_t);
std::string canonical_mo_display_label(const Wavefunction&,std::size_t,
                                     const NboAomoName* = nullptr);
std::string orbital_irrep_display_label(const NboAomoName&);
std::string canonical_mo_current_irrep(const Wavefunction&,std::size_t,
                                     const NboAomoName* = nullptr);
} // namespace cov::ui
