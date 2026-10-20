#pragma once
#include <limits>
#include <string>
#include <vector>

namespace cov {
enum class LigandPiPrior { Unresolved, SigmaOnly, Donor, Acceptor, Ambiguous };
enum class PiPairDirection { Unresolved, Donor, Acceptor, Coupled, WeakNearNonbonding };
struct PiFrozenOperatorAssessment {
    bool available=false, tracking_verified=false, occupation_boundary_preserved=false;
    std::string reason;
    double max_energy_shift_hartree=std::numeric_limits<double>::quiet_NaN();
    double full_spectrum_max_energy_shift_hartree=std::numeric_limits<double>::quiet_NaN();
    double max_group_width_change_hartree=std::numeric_limits<double>::quiet_NaN();
    double frontier_gap_change_hartree=std::numeric_limits<double>::quiet_NaN();
    double max_subspace_sin2=std::numeric_limits<double>::quiet_NaN();
    double minimum_external_gap_hartree=std::numeric_limits<double>::quiet_NaN();
    double removal_norm_hartree=std::numeric_limits<double>::quiet_NaN();
    double numerical_error_bound_hartree=std::numeric_limits<double>::quiet_NaN();
};
struct PiFrozenSpectralGroup {
    std::vector<std::size_t> members;
    PiFrozenOperatorAssessment assessment;
};
struct PiDisplayCalibration {
    // Set only for a scope/method family that passed the recorded held-out scan.
    bool validated=false;
    std::string version, family;
    double energy_budget_ev=0.025, subspace_sin2_budget=0.02;
};
[[nodiscard]] bool pi_display_negligible(const PiFrozenOperatorAssessment&,
    const PiDisplayCalibration&);
struct PiEndpointNboWeights {
    double ligand_to_centre_donor=0, ligand_to_centre_acceptor=0;
    double centre_to_ligand_donor=0, centre_to_ligand_acceptor=0;
    double occupation=std::numeric_limits<double>::quiet_NaN();
};
struct PiEndpointDirectionAssessment {
    std::string direction="unresolved", reason;
    bool verified=false, ligand_to_centre_supported=false, centre_to_ligand_supported=false;
    double ligand_to_centre_coverage=0, centre_to_ligand_coverage=0;
};
// Only set verified for a common, certified operation/atom/spin domain.
// Equal point-group names alone are not a common representation domain.
struct PiEndpointSymmetryEvidence {
    bool verified=false;
    std::string scope_id, irrep;
};
struct PiModePairAssessment {
    bool verified=false, direction_verified=false;
    bool ordinary_display_eligible=false, multi_group_relation=false;
    bool two_endpoint_relation=false;
    std::string direction="unresolved", channel_family, reason;
    std::vector<std::string> shared_mode_ids, matched_edge_ids;
    double lower_coverage=0, upper_coverage=0;
    double lower_role_coverage=0, upper_role_coverage=0;
    double lower_cross_fock_mean_hartree=0, upper_cross_fock_mean_hartree=0;
    double shared_fragment_contraction_norm_hartree=0;
    double numerical_coverage_bound=1;
    // A majority of each complete endpoint's original norm must belong to
    // the declared common mode space for a primary explanation. This is a
    // semantic majority rule, not a chemical-energy significance threshold.
    double primary_coverage_floor=0.5;
};
struct PiModeNetworkNode {
    std::size_t group_index=0;
    std::vector<std::size_t> members;
    double centre_weight=0, ligand_weight=0;
    double donor_role_coverage=0, acceptor_role_coverage=0;
    std::string character;
    bool primary=false;
};
struct PiModeNetworkAssessment {
    std::string channel_id, mode_id, spin, channel_family, ligand_space_kind;
    std::string direction="unresolved", reason;
    bool verified=false, direction_verified=false, ordinary_display_eligible=false;
    bool two_endpoint_relation=false;
    double numerical_coverage_bound=1;
    double primary_coverage_floor=0.5;
    std::vector<PiModeNetworkNode> nodes;
    std::vector<std::string> matched_edge_ids;
};
[[nodiscard]] std::string pi_mode_network_assessment_json(const PiModeNetworkAssessment&);
// Mapping support is tested against its propagated numerical error, not a
// fitted composition percentage. The caller supplies complete canonical groups.
[[nodiscard]] PiEndpointDirectionAssessment assess_pi_endpoint_direction(
    const PiEndpointNboWeights& lower,const PiEndpointNboWeights& upper,
    const std::string& scope_direction,bool scope_direction_verified,
    bool occupations_verified,double full_occupation,double mapping_error_bound);
struct PiPartnerComponents {
    double energy_hartree=0;
    double metal_d=0;
    double ligand_p=0;
    double pi_fraction=0;
    double metal_ligand_overlap=0;
};
// A counterpart is a relation within one independently verified local Fock
// channel. Composition contrast alone cannot establish this relation.
struct PiPartnerChannelEvidence {
    std::string channel_id, canonical_fingerprint, spin;
    std::vector<std::size_t> lower_members, upper_members;
    bool same_operator_verified=false, occupations_verified=false;
    bool complete_membership_verified=false;
    std::string operator_kind;
    std::string direction="unresolved";
    std::string lower_character, upper_character;
    double lower_cross_fock_max_hartree=0, upper_cross_fock_min_hartree=0;
    double operator_error_hartree=0;
    double operator_tolerance_hartree=0;
    bool direction_verified=false;
    std::string direction_reference;
    PiFrozenOperatorAssessment frozen_operator;
    PiDisplayCalibration display_calibration;
    PiModePairAssessment mode_assessment;
};
struct PiPartnerAssessment {
    bool input_valid=false;
    bool accepted=false;
    PiPairDirection direction=PiPairDirection::Unresolved;
    LigandPiPrior prior=LigandPiPrior::Unresolved;
    // Verified channels rank by minimum endpoint cross-Fock magnitude (hartree).
    // Verified ranking is endpoint cross-Fock support; support is not a probability.
    double ranking_score=std::numeric_limits<double>::quiet_NaN();
    double support_score=std::numeric_limits<double>::quiet_NaN();
    double splitting_hartree=std::numeric_limits<double>::quiet_NaN();
    double composition_contrast=std::numeric_limits<double>::quiet_NaN();
    bool complementary_composition=false;
    bool opposite_overlap_signs=false;
    bool weak=false;
    std::string prior_relation="undetermined";
    std::string detail;
    PiPartnerComponents lower;
    PiPartnerComponents upper;
    PiPartnerChannelEvidence channel;
};
struct WeakCrystalFieldAssessment {
    bool input_valid=false;
    bool accepted=false;
    double splitting_hartree=std::numeric_limits<double>::quiet_NaN();
    double support_score=std::numeric_limits<double>::quiet_NaN();
    double split_threshold_hartree=std::numeric_limits<double>::quiet_NaN();
    double overlap_threshold=std::numeric_limits<double>::quiet_NaN();
    PiPartnerComponents first;
    PiPartnerComponents second;
    std::string detail;
};
// Numerical weak-field screen only. Callers must independently establish the
// compatible local axes, spin, d-shell irreps and all orbital identities.
[[nodiscard]] WeakCrystalFieldAssessment assess_weak_crystal_field(
    const PiPartnerComponents& first,const PiPartnerComponents& second,
    double split_threshold_hartree=0.020,double overlap_threshold=0.025);
[[nodiscard]] std::string weak_crystal_field_assessment_json(const WeakCrystalFieldAssessment&);
[[nodiscard]] PiPartnerAssessment assess_pi_partner(
    const PiPartnerComponents& lower,const PiPartnerComponents& upper,LigandPiPrior,
    double weak_split_hartree=0.020,double weak_overlap=0.025);
[[nodiscard]] PiPartnerAssessment assess_pi_channel_partner(
    const PiPartnerComponents& lower,const PiPartnerComponents& upper,
    LigandPiPrior,const PiPartnerChannelEvidence&);
[[nodiscard]] const char* ligand_pi_prior_name(LigandPiPrior) noexcept;
[[nodiscard]] const char* pi_pair_direction_name(PiPairDirection) noexcept;
[[nodiscard]] std::string pi_partner_assessment_json(const PiPartnerAssessment&);
}
