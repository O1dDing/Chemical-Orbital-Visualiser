#pragma once
#include "cov/nbo_integration.hpp"
#include "cov/pi_pair_evidence.hpp"
#include <utility>

namespace cov {
// All indices are zero-based except literal producer NAO IDs. The local
// projectors, not a chosen SVD vector gauge, define each coupled channel.
struct NboPiCanonicalGroup {
    std::vector<std::size_t> members;
    NboSpin spin=NboSpin::Total;
    double centre_weight=0, ligand_weight=0;
    std::optional<double> occupation_per_mo;
    double cross_fock_min_hartree=0, cross_fock_max_hartree=0,
           cross_fock_mean_hartree=0;
    std::string character="unresolved";
    // Tr(C_group^T P_ordered-role C_group)/rank, using unique complete NBO
    // vectors (not printed row counts). The two directions remain separate.
    double ligand_to_centre_donor_weight=0, ligand_to_centre_acceptor_weight=0;
    double centre_to_ligand_donor_weight=0, centre_to_ligand_acceptor_weight=0;
};
struct NboPiAngularEvidence {
    std::size_t atom=0;
    std::string method="AO S-metric local angular projector in bond-axis frame";
    std::size_t rank=2;
    double p_metric_min_eigenvalue=0, max_sigma_leakage=0,
           centre_projection=0, partition_residual=0;
    std::array<double,3> pi_generalized_eigenvalues{};
};
struct NboPiDirectionProjection {
    std::string e2_id;
    std::size_t donor_nbo_id=0, acceptor_nbo_id=0;
    double donor_ligand_weight=0, donor_centre_weight=0,
           acceptor_ligand_weight=0, acceptor_centre_weight=0;
    NboSource source;
    std::string method, status, reason;
    double donor_occupation=0, acceptor_occupation=0;
    double energy_gap_hartree=0, fock_hartree=0, e2_hartree=0;
    double actual_occupation_second_order_hartree=0, perturbation_occupation=0;
    std::string e2_reference="unavailable";
    bool recovered_from_unprinted_fock=false;
    std::string direction="unresolved";
};
struct NboPiModeGroup {
    std::size_t group_index=0;
    double centre_weight=0, ligand_weight=0;
    double cross_fock_min_hartree=0, cross_fock_max_hartree=0,
           cross_fock_mean_hartree=0;
    // Mode-side projections of the whole canonical group. Their gauge is
    // shared within this mode; only contracted norms/traces are meaningful.
    NboMatrix centre_coordinates, ligand_coordinates;
};
struct NboPiModeEdgeGroup {
    std::size_t group_index=0;
    // Original orthonormal NBO-role weights, used without renormalisation.
    double donor_weight=0, acceptor_weight=0;
    // Weights after projecting that very same source edge into this mode.
    double donor_mode_weight=0, acceptor_mode_weight=0;
};
struct NboPiModeEdge {
    std::string id, direction;
    std::size_t donor_nbo_id=0, acceptor_nbo_id=0;
    double projected_fock_hartree=0;
    std::vector<NboPiModeEdgeGroup> groups;
};
struct NboPiModeShell {
    std::string shell, angular;
    double fraction=0;
};
struct NboPiCouplingMode {
    std::string id, centre_space_kind, ligand_space_kind;
    std::size_t rank=0;
    double singular_min_hartree=0, singular_max_hartree=0;
    std::vector<double> singular_values_hartree;
    double numerical_coverage_bound=1;
    NboMatrix centre_projector_basis, ligand_projector_basis;
    std::vector<NboPiModeShell> centre_shells;
    std::vector<NboPiModeGroup> groups;
    std::vector<NboPiModeEdge> ordered_edges;
};
struct NboPiCoupling {
    std::string id, channel="pi", direction="unresolved", direction_evidence;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> centre_atoms, ligand_atoms;
    std::vector<std::size_t> centre_nao_ids, ligand_nao_ids;
    std::vector<double> singular_values_hartree;
    std::size_t coupled_rank=0;
    double centre_onsite_hartree=0, ligand_onsite_hartree=0;
    std::optional<double> centre_occupation, ligand_occupation;
    std::array<double,2> centre_onsite_range_hartree{}, ligand_onsite_range_hartree{};
    std::optional<std::array<double,2>> centre_occupation_range, ligand_occupation_range;
    std::string occupation_status="insufficient_evidence", occupation_reason;
    double operator_max_error_hartree=0, nao_orthogonality_error=0;
    double angular_leakage=0, minimum_centre_projection=0;
    std::optional<double> direction_donor_weight, direction_acceptor_weight;
    std::vector<NboPiAngularEvidence> angular_projector_evidence;
    std::vector<NboPiDirectionProjection> direction_projection_evidence;
    std::vector<NboPiCanonicalGroup> groups;
    // Empty for the broad atom-pi projector. Nonempty families retain full
    // ligand support and independently assigned internal pi/antipi character.
    std::string ligand_family;
    std::string centre_family;
    std::string ligand_space_kind="transverse-p-basis";
    std::vector<std::size_t> ligand_family_nbo_ids;
    bool localized_family_verified=false;
    std::string operator_kind="canonical-same-operator";
    double canonical_operator_residual_hartree=0;
    double operator_validation_tolerance_hartree=2e-5;
    bool canonical_members_are_verified_shared_spatial=false;
    bool direction_verified=false;
    double direction_mapping_error_bound=1;
    std::string direction_reference="ordered-localized-NBO/same-spin-Fock";
    // Full orthonormal projector columns in the verified NAO reference.
    // Scientific scope comparisons use these, never names or atom lists.
    NboMatrix centre_projector_basis, ligand_projector_basis;
    PiFrozenOperatorAssessment frozen_operator;
    std::vector<PiFrozenSpectralGroup> frozen_groups;
    NboSource source;
    std::vector<NboPiCouplingMode> modes;
};
// Evaluate common modes and linked source edges, never the independent union
// of all donors and all acceptors. Source members are immutable complete groups.
PiModePairAssessment assess_pi_mode_pair(const NboPiCoupling&,
    const std::vector<std::size_t>& lower_members,
    const std::vector<std::size_t>& upper_members,
    const PiEndpointSymmetryEvidence& lower_symmetry={},
    const PiEndpointSymmetryEvidence& upper_symmetry={});
// One mode is one network of complete source groups, never an N-squared list
// of pair energy gaps. A two-endpoint descriptor is a separately verified
// dominant reduction; the full network remains available in the source mode.
PiModeNetworkAssessment assess_pi_mode_network(const NboPiCoupling&,
    const NboPiCouplingMode&);
// All matrices use the same orthonormal NAO reference. Optional source-column
// mapping connects immutable global canonical indices to this spin's columns.
// Recomputes modes without changing the parent channel, source or canonical MOs.
void populate_pi_coupling_modes(NboPiCoupling&,const NboMatrix& fock,
    const NboMatrix& canonical_columns,const NboMatrix* localized_columns=nullptr,
    const std::vector<std::size_t>& canonical_source_columns={});
struct NboPiCouplingAnalysis {
    std::string status="insufficient_evidence", reason;
    std::vector<NboPiCoupling> couplings;
    std::vector<NboSource> evidence;
};
// Same-operator, rotation-covariant local p-pi projections. No canonical or
// producer orbital is modified. A failed operator gate yields no pseudozero.
NboPiCouplingAnalysis analyse_nbo_pi_couplings(const Wavefunction& canonical,
    const NboIntegration& integration,
    const std::vector<std::pair<std::size_t,std::size_t>>& strong_connectivity);
// Frozen Fock removal in one orthonormal reference. canonical_columns must
// diagonalize fock; groups contain complete source column indices. This is
// a display sensitivity diagnostic, not a self-consistent deletion energy.
PiFrozenOperatorAssessment assess_pi_frozen_operator(const NboMatrix& fock,
    const NboMatrix& centre_basis,const NboMatrix& ligand_basis,
    const NboMatrix& canonical_columns,
    const std::vector<std::vector<std::size_t>>& groups,
    const std::vector<double>& occupations,double operator_error_hartree,
    std::vector<PiFrozenSpectralGroup>* per_group=nullptr);
// source_members must be fixed by the display domain before evaluating this
// diagnostic. Every touched degenerate spectral cluster must be complete.
PiFrozenOperatorAssessment assess_pi_frozen_display_scope(const NboPiCoupling&,
    const std::vector<std::size_t>& source_members);
// A versioned display-error calibration for structurally defined channel
// families. Eligibility never comes from a molecule/ligand name or MO number.
PiDisplayCalibration pi_channel_display_calibration(const NboPiCoupling&);
} // namespace cov
