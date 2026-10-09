#include "cov/pi_pair_evidence.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace cov {
bool pi_display_negligible(const PiFrozenOperatorAssessment& a,const PiDisplayCalibration& c) {
    const auto nonnegative=[](double x){return std::isfinite(x)&&x>=0;};
    if(!c.validated||c.version.empty()||c.family.empty()||!a.available||!a.tracking_verified||
       !a.occupation_boundary_preserved||!nonnegative(c.energy_budget_ev)||c.energy_budget_ev==0||
       !nonnegative(c.subspace_sin2_budget)||c.subspace_sin2_budget>=0.5||
       !nonnegative(a.numerical_error_bound_hartree)||!nonnegative(a.max_energy_shift_hartree)||
       !nonnegative(a.max_group_width_change_hartree)||!nonnegative(a.frontier_gap_change_hartree)||
       !nonnegative(a.max_subspace_sin2))return false;
    if(std::isnan(a.minimum_external_gap_hartree)||a.minimum_external_gap_hartree<=0)return false;
    // Bound uncertainty in the tracked spectral projector as well as energy.
    const double angle_error=2*a.numerical_error_bound_hartree/a.minimum_external_gap_hartree;
    const double sin2_upper=std::pow(std::min(1.0,std::sqrt(a.max_subspace_sin2)+angle_error),2);
    const double budget=c.energy_budget_ev/27.211386245988;
    return a.max_energy_shift_hartree+a.numerical_error_bound_hartree<=budget&&
        a.max_group_width_change_hartree+2*a.numerical_error_bound_hartree<=budget&&
        a.frontier_gap_change_hartree+2*a.numerical_error_bound_hartree<=budget&&
        sin2_upper<=c.subspace_sin2_budget;
}

PiEndpointDirectionAssessment assess_pi_endpoint_direction(
    const PiEndpointNboWeights& a,const PiEndpointNboWeights& b,
    const std::string& scope,bool scope_verified,bool occupations_verified,
    double full_occupation,double error) {
    PiEndpointDirectionAssessment out;
    const auto valid=[](double x){return std::isfinite(x)&&x>=0&&x<=1+1e-5;};
    if(!std::isfinite(error)||error<0||error>=1||!std::isfinite(full_occupation)||full_occupation<=0||
       !valid(a.ligand_to_centre_donor)||!valid(a.ligand_to_centre_acceptor)||
       !valid(a.centre_to_ligand_donor)||!valid(a.centre_to_ligand_acceptor)||
       !valid(b.ligand_to_centre_donor)||!valid(b.ligand_to_centre_acceptor)||
       !valid(b.centre_to_ligand_donor)||!valid(b.centre_to_ligand_acceptor)) {
        out.reason="ordered-NBO-group-mapping-unavailable";return out;
    }
    if(occupations_verified&&std::isfinite(a.occupation)&&std::isfinite(b.occupation)&&
       std::abs(a.occupation-full_occupation)<=1e-6&&std::abs(b.occupation-full_occupation)<=1e-6){
        out.direction="occupied_space_mixing";out.verified=true;
        out.reason="both-complete-canonical-endpoints-occupied-no-electron-transfer-arrow";return out;
    }
    if(!scope_verified){out.reason="ordered-NBO-scope-direction-unavailable";return out;}
    out.ligand_to_centre_coverage=std::max(std::min(a.ligand_to_centre_donor,b.ligand_to_centre_acceptor),
        std::min(b.ligand_to_centre_donor,a.ligand_to_centre_acceptor));
    out.centre_to_ligand_coverage=std::max(std::min(a.centre_to_ligand_donor,b.centre_to_ligand_acceptor),
        std::min(b.centre_to_ligand_donor,a.centre_to_ligand_acceptor));
    const double numerical_gate=std::max(1e-8,error);
    out.ligand_to_centre_supported=(scope=="ligand_to_centre"||scope=="bidirectional")&&out.ligand_to_centre_coverage>numerical_gate;
    out.centre_to_ligand_supported=(scope=="centre_to_ligand"||scope=="bidirectional")&&out.centre_to_ligand_coverage>numerical_gate;
    out.verified=out.ligand_to_centre_supported||out.centre_to_ligand_supported;
    if(out.ligand_to_centre_supported&&out.centre_to_ligand_supported)out.direction="bidirectional";
    else if(out.ligand_to_centre_supported)out.direction="ligand_to_centre";
    else if(out.centre_to_ligand_supported)out.direction="centre_to_ligand";
    out.reason=out.verified?"ordered-donor-and-acceptor-spans-map-to-distinct-complete-endpoints-above-numerical-error":
        "this-endpoint-pair-does-not-capture-both-ordered-NBO-roles-above-mapping-error";
    return out;
}
const char* ligand_pi_prior_name(LigandPiPrior p) noexcept {
    switch(p) {
        case LigandPiPrior::SigmaOnly:return "sigma-only";
        case LigandPiPrior::Donor:return "donor";
        case LigandPiPrior::Acceptor:return "acceptor";
        case LigandPiPrior::Ambiguous:return "ambiguous";
        default:return "unresolved";
    }
}
const char* pi_pair_direction_name(PiPairDirection p) noexcept {
    switch(p) {
        case PiPairDirection::Donor:return "donor";
        case PiPairDirection::Acceptor:return "acceptor";
        case PiPairDirection::Coupled:return "coupled";
        case PiPairDirection::WeakNearNonbonding:return "weak-near-nonbonding";
        default:return "unresolved";
    }
}
PiPartnerAssessment assess_pi_partner(const PiPartnerComponents& a,const PiPartnerComponents& b,
    LigandPiPrior prior,double weak_split,double weak_overlap) {
    PiPartnerAssessment result;result.prior=prior;result.lower=a;result.upper=b;
    const auto valid=[](const auto& c) {
        return std::isfinite(c.energy_hartree) && std::isfinite(c.metal_d) && c.metal_d>=0 &&
               std::isfinite(c.ligand_p) && c.ligand_p>=0 && std::isfinite(c.pi_fraction) &&
               c.pi_fraction>=0 && std::isfinite(c.metal_ligand_overlap);
    };
    result.input_valid=valid(a)&&valid(b)&&std::isfinite(weak_split)&&weak_split>=0&&
        std::isfinite(weak_overlap)&&weak_overlap>=0&&b.energy_hartree>=a.energy_hartree;
    if(!result.input_valid){result.detail="invalid-input";return result;}
    result.splitting_hartree=b.energy_hartree-a.energy_hartree;
    const double low=a.ligand_p-a.metal_d,high=b.ligand_p-b.metal_d;
    const double donor=low-high;
    result.composition_contrast=donor;
    result.complementary_composition=low*high<0&&std::abs(low)>=.05&&std::abs(high)>=.05;
    result.opposite_overlap_signs=a.metal_ligand_overlap*b.metal_ligand_overlap<0&&
        std::abs(a.metal_ligand_overlap)>=weak_overlap&&std::abs(b.metal_ligand_overlap)>=weak_overlap;
    // Composition, gap and overlap remain observations. Without the actual
    // same-source operator and a calibrated frozen-removal diagnostic they
    // cannot establish either a directional channel or a negligible one.
    result.weak=false;
    result.detail="composition-only-candidate-needs-verified-local-channel-and-calibrated-sensitivity";
    return result;
}

PiPartnerAssessment assess_pi_channel_partner(const PiPartnerComponents& a,
    const PiPartnerComponents& b,LigandPiPrior prior,const PiPartnerChannelEvidence& channel) {
    // Reuse finite-input validation and descriptive quantities, not its
    // composition-only acceptance decision or heuristic direction.
    auto result=assess_pi_partner(a,b,prior);
    result.channel=channel;result.weak=false;result.accepted=false;result.direction=PiPairDirection::Unresolved;
    result.prior_relation="undetermined";
    if(!result.input_valid)return result;
    if(channel.channel_id.empty()||channel.canonical_fingerprint.empty()||
       channel.spin.empty()||channel.lower_members.empty()||channel.upper_members.empty()||
       !channel.same_operator_verified||!channel.complete_membership_verified||channel.operator_kind.empty()||
       !std::isfinite(channel.operator_error_hartree)||
       channel.operator_error_hartree<0||!std::isfinite(channel.operator_tolerance_hartree)||
       channel.operator_tolerance_hartree<=0||channel.operator_error_hartree>channel.operator_tolerance_hartree||
       !std::isfinite(channel.lower_cross_fock_max_hartree)||
       !std::isfinite(channel.upper_cross_fock_min_hartree)) {
        result.detail="verified-channel-identity-or-operator-unavailable";return result;
    }
    for(const auto* members:{&channel.lower_members,&channel.upper_members}) {
        auto copy=*members;std::sort(copy.begin(),copy.end());
        if(std::adjacent_find(copy.begin(),copy.end())!=copy.end()) {
            result.detail="duplicate-counterpart-member-identity";return result;
        }
    }
    for(const auto i:channel.lower_members)
        if(std::find(channel.upper_members.begin(),channel.upper_members.end(),i)!=channel.upper_members.end()) {
        result.detail="overlapping-counterpart-members";return result;
    }
    if(!channel.mode_assessment.verified||channel.mode_assessment.shared_mode_ids.empty()) {
        result.detail=channel.mode_assessment.reason.empty()?"verified-common-fragment-mode-unavailable":channel.mode_assessment.reason;return result;
    }
    if(!std::isfinite(channel.mode_assessment.lower_cross_fock_mean_hartree)||
       !std::isfinite(channel.mode_assessment.upper_cross_fock_mean_hartree)||
       !(channel.mode_assessment.lower_cross_fock_mean_hartree<0)||
       !(channel.mode_assessment.upper_cross_fock_mean_hartree>0)) {
        result.detail="common-mode-cross-fock-traces-unresolved";return result;
    }
    result.accepted=true;result.weak=pi_display_negligible(channel.frozen_operator,channel.display_calibration);
    result.direction=PiPairDirection::Coupled;
    result.channel.direction=channel.mode_assessment.direction;
    result.channel.direction_verified=channel.mode_assessment.direction_verified;
    if(channel.mode_assessment.direction_verified&&channel.mode_assessment.direction=="centre_to_ligand")result.direction=PiPairDirection::Acceptor;
    if(channel.mode_assessment.direction_verified&&channel.mode_assessment.direction=="ligand_to_centre")result.direction=PiPairDirection::Donor;
    result.ranking_score=std::min(-channel.mode_assessment.lower_cross_fock_mean_hartree,
                                channel.mode_assessment.upper_cross_fock_mean_hartree);
    result.support_score=std::numeric_limits<double>::quiet_NaN();
    if(result.direction==PiPairDirection::Acceptor||result.direction==PiPairDirection::Donor) {
        if(prior==LigandPiPrior::Acceptor||prior==LigandPiPrior::Donor)
            result.prior_relation=(prior==LigandPiPrior::Acceptor)==(result.direction==PiPairDirection::Acceptor)
                ?"consistent":"contradicted";
    }
    result.detail=channel.mode_assessment.multi_group_relation?
        "shared-mode-multigroup-membership-not-unique-two-level-pair":"shared-verified-fragment-mode-canonical-energy-separation";
    return result;
}
namespace {
void number(std::ostream& out,double v){if(std::isfinite(v))out<<v;else out<<"null";}
void json_string(std::ostream& out,const std::string& value) {
    static constexpr char hex[]="0123456789abcdef";
    out << '"';
    for(unsigned char c:value) {
        if(c=='"'||c=='\\')out << '\\' << static_cast<char>(c);
        else if(c<0x20)out << "\\u00" << hex[c>>4] << hex[c&15];
        else out << static_cast<char>(c);
    }
    out << '"';
}
void components(std::ostream& out,const PiPartnerComponents& v){
    out<<"{\"energy_hartree\":";number(out,v.energy_hartree);
    out<<",\"metal_d\":";number(out,v.metal_d);out<<",\"ligand_p\":";number(out,v.ligand_p);
    out<<",\"pi_fraction\":";number(out,v.pi_fraction);out<<",\"metal_ligand_overlap\":";number(out,v.metal_ligand_overlap);out<<'}';
}
}
std::string pi_partner_assessment_json(const PiPartnerAssessment& v){
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"{\"input_valid\":"<<(v.input_valid?"true":"false")<<",\"accepted\":"<<(v.accepted?"true":"false")
       <<",\"direction\":\""<<pi_pair_direction_name(v.direction)<<"\",\"catalogue_prior\":\""<<ligand_pi_prior_name(v.prior)
       <<"\",\"prior_relation\":";json_string(out,v.prior_relation);out<<",\"detail\":";json_string(out,v.detail);
    out<<",\"ranking_score\":";number(out,v.ranking_score);out<<",\"support_score\":";number(out,v.support_score);
    const bool local_channel=!v.channel.channel_id.empty();
    out<<",\"ranking_measure\":\""<<(local_channel?"minimum-endpoint-common-mode-cross-Fock-trace-per-member":"weak-composition-heuristic")
       <<"\",\"ranking_unit\":\""<<(local_channel?"hartree":"dimensionless")
       <<"\",\"score_meaning\":\""<<(local_channel?"cross-Fock-strength-not-probability":"heuristic-support-not-probability")
       <<"\",\"splitting_hartree\":";number(out,v.splitting_hartree);
    out<<",\"composition_contrast\":";number(out,v.composition_contrast);
    out<<",\"complementary_composition\":"<<(v.complementary_composition?"true":"false")
       <<",\"opposite_overlap_signs\":"<<(v.opposite_overlap_signs?"true":"false")<<",\"weak\":"<<(v.weak?"true":"false")<<",\"lower\":";
    components(out,v.lower);out<<",\"upper\":";components(out,v.upper);
    out<<",\"channel_id\":";json_string(out,v.channel.channel_id);
    out<<",\"canonical_fingerprint\":";json_string(out,v.channel.canonical_fingerprint);
    out<<",\"channel_spin\":";json_string(out,v.channel.spin);
    out<<",\"operator_kind\":";json_string(out,v.channel.operator_kind);
    out<<",\"channel_direction\":";json_string(out,v.channel.direction);
    out<<",\"direction_verified\":"<<(v.channel.direction_verified?"true":"false");
    out<<",\"direction_reference\":";json_string(out,v.channel.direction_reference);
    const auto& mode=v.channel.mode_assessment;
    out<<",\"mode_assessment\":{\"verified\":"<<(mode.verified?"true":"false")
       <<",\"direction_verified\":"<<(mode.direction_verified?"true":"false")
       <<",\"ordinary_display_eligible\":"<<(mode.ordinary_display_eligible?"true":"false")
       <<",\"multi_group_relation\":"<<(mode.multi_group_relation?"true":"false");
    out<<",\"two_endpoint_relation\":"<<(mode.two_endpoint_relation?"true":"false");
    out<<",\"direction\":";json_string(out,mode.direction);
    out<<",\"channel_family\":";json_string(out,mode.channel_family);
    out<<",\"reason\":";json_string(out,mode.reason);
    out<<",\"lower_coverage\":";number(out,mode.lower_coverage);
    out<<",\"upper_coverage\":";number(out,mode.upper_coverage);
    out<<",\"lower_role_coverage\":";number(out,mode.lower_role_coverage);
    out<<",\"upper_role_coverage\":";number(out,mode.upper_role_coverage);
    out<<",\"lower_cross_fock_mean_hartree\":";number(out,mode.lower_cross_fock_mean_hartree);
    out<<",\"upper_cross_fock_mean_hartree\":";number(out,mode.upper_cross_fock_mean_hartree);
    out<<",\"shared_fragment_contraction_norm_hartree\":";number(out,mode.shared_fragment_contraction_norm_hartree);
    out<<",\"numerical_coverage_bound\":";number(out,mode.numerical_coverage_bound);
    out<<",\"primary_coverage_floor\":";number(out,mode.primary_coverage_floor);
    out<<",\"shared_mode_ids\":[";for(std::size_t i=0;i<mode.shared_mode_ids.size();++i){if(i)out<<',';json_string(out,mode.shared_mode_ids[i]);}
    out<<"],\"matched_edge_ids\":[";for(std::size_t i=0;i<mode.matched_edge_ids.size();++i){if(i)out<<',';json_string(out,mode.matched_edge_ids[i]);}out<<"]}";
    out<<",\"display_calibration\":{\"validated\":"<<(v.channel.display_calibration.validated?"true":"false");
    out<<",\"version\":";json_string(out,v.channel.display_calibration.version);
    out<<",\"family\":";json_string(out,v.channel.display_calibration.family);
    out<<",\"energy_budget_ev\":";number(out,v.channel.display_calibration.energy_budget_ev);
    out<<",\"subspace_sin2_budget\":";number(out,v.channel.display_calibration.subspace_sin2_budget);out<<'}';
    const auto& frozen=v.channel.frozen_operator;
    out<<",\"frozen_operator\":{\"available\":"<<(frozen.available?"true":"false")
       <<",\"tracking_verified\":"<<(frozen.tracking_verified?"true":"false")
       <<",\"occupation_boundary_preserved\":"<<(frozen.occupation_boundary_preserved?"true":"false");
    out<<",\"reason\":";json_string(out,frozen.reason);
    out<<",\"max_energy_shift_hartree\":";number(out,frozen.max_energy_shift_hartree);
    out<<",\"full_spectrum_max_energy_shift_hartree\":";number(out,frozen.full_spectrum_max_energy_shift_hartree);
    out<<",\"max_group_width_change_hartree\":";number(out,frozen.max_group_width_change_hartree);
    out<<",\"frontier_gap_change_hartree\":";number(out,frozen.frontier_gap_change_hartree);
    out<<",\"max_subspace_sin2\":";number(out,frozen.max_subspace_sin2);
    out<<",\"minimum_external_gap_hartree\":";number(out,frozen.minimum_external_gap_hartree);
    out<<",\"removal_norm_hartree\":";number(out,frozen.removal_norm_hartree);
    out<<",\"numerical_error_bound_hartree\":";number(out,frozen.numerical_error_bound_hartree);out<<'}';

    out<<",\"lower_character\":";json_string(out,v.channel.lower_character);
    out<<",\"upper_character\":";json_string(out,v.channel.upper_character);
    out<<",\"lower_cross_fock_max_hartree\":";number(out,v.channel.lower_cross_fock_max_hartree);
    out<<",\"upper_cross_fock_min_hartree\":";number(out,v.channel.upper_cross_fock_min_hartree);
    out<<",\"operator_error_hartree\":";number(out,v.channel.operator_error_hartree);
    out<<",\"operator_tolerance_hartree\":";number(out,v.channel.operator_tolerance_hartree);
    out<<",\"lower_members\":[";for(std::size_t i=0;i<v.channel.lower_members.size();++i){if(i)out<<',';out<<v.channel.lower_members[i];}
    out<<"],\"upper_members\":[";for(std::size_t i=0;i<v.channel.upper_members.size();++i){if(i)out<<',';out<<v.channel.upper_members[i];}
    out<<"],\"same_operator_verified\":"<<(v.channel.same_operator_verified?"true":"false")
       <<",\"occupations_verified\":"<<(v.channel.occupations_verified?"true":"false")
       <<",\"complete_membership_verified\":"<<(v.channel.complete_membership_verified?"true":"false")
       <<",\"energy_semantics\":\"canonical-group-energy-separation\"}";return out.str();
}
WeakCrystalFieldAssessment assess_weak_crystal_field(const PiPartnerComponents& a,
    const PiPartnerComponents& b,double split_threshold,double overlap_threshold) {
    WeakCrystalFieldAssessment result;result.first=a;result.second=b;
    result.split_threshold_hartree=split_threshold;result.overlap_threshold=overlap_threshold;
    const auto valid=[](const auto& c) {
        return std::isfinite(c.energy_hartree)&&std::isfinite(c.metal_d)&&c.metal_d>=0&&
            std::isfinite(c.metal_ligand_overlap);
    };
    result.input_valid=valid(a)&&valid(b)&&std::isfinite(split_threshold)&&split_threshold>=0&&
        std::isfinite(overlap_threshold)&&overlap_threshold>=0;
    if(!result.input_valid){result.detail="invalid-input";return result;}
    result.splitting_hartree=std::abs(b.energy_hartree-a.energy_hartree);
    if(a.metal_d<.60||b.metal_d<.60||result.splitting_hartree>split_threshold||
        std::abs(a.metal_ligand_overlap)>overlap_threshold||
        std::abs(b.metal_ligand_overlap)>overlap_threshold) {
        result.detail="weak-local-d-evidence-insufficient";return result;
    }
    result.support_score=std::clamp(split_threshold>0?
        1-result.splitting_hartree/split_threshold:1,0.0,1.0);
    result.accepted=result.support_score>=.20;
    result.detail=result.accepted?"weak-local-d-separation":"weak-gap-support-insufficient";
    return result;
}
std::string weak_crystal_field_assessment_json(const WeakCrystalFieldAssessment& v) {
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"{\"input_valid\":"<<(v.input_valid?"true":"false")<<",\"accepted\":"<<(v.accepted?"true":"false");
    out<<",\"splitting_hartree\":";number(out,v.splitting_hartree);
    out<<",\"support_score\":";number(out,v.support_score);
    out<<",\"score_meaning\":\"heuristic-support-not-probability\",\"split_threshold_hartree\":";
    number(out,v.split_threshold_hartree);out<<",\"overlap_threshold\":";number(out,v.overlap_threshold);
    out<<",\"detail\":";json_string(out,v.detail);out<<",\"first\":";components(out,v.first);
    out<<",\"second\":";components(out,v.second);out<<'}';return out.str();
}
std::string pi_mode_network_assessment_json(const PiModeNetworkAssessment& v) {
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"{\"channel_id\":";json_string(out,v.channel_id);
    out<<",\"mode_id\":";json_string(out,v.mode_id);
    out<<",\"spin\":";json_string(out,v.spin);
    out<<",\"channel_family\":";json_string(out,v.channel_family);
    out<<",\"ligand_space_kind\":";json_string(out,v.ligand_space_kind);
    out<<",\"direction\":";json_string(out,v.direction);
    out<<",\"verified\":"<<(v.verified?"true":"false")
       <<",\"direction_verified\":"<<(v.direction_verified?"true":"false")
       <<",\"ordinary_display_eligible\":"<<(v.ordinary_display_eligible?"true":"false")
       <<",\"two_endpoint_relation\":"<<(v.two_endpoint_relation?"true":"false");
    out<<",\"reason\":";json_string(out,v.reason);
    out<<",\"numerical_coverage_bound\":";number(out,v.numerical_coverage_bound);
    out<<",\"primary_coverage_floor\":";number(out,v.primary_coverage_floor);
    out<<",\"weight_definition\":\"trace-of-complete-source-group-projector-divided-by-original-group-rank; no-renormalization\"";
    out<<",\"nodes\":[";for(std::size_t i=0;i<v.nodes.size();++i){const auto& node=v.nodes[i];if(i)out<<',';
        out<<"{\"group_index\":"<<node.group_index<<",\"members\":[";
        for(std::size_t j=0;j<node.members.size();++j){if(j)out<<',';out<<node.members[j];}
        out<<"],\"centre_weight\":";number(out,node.centre_weight);
        out<<",\"ligand_weight\":";number(out,node.ligand_weight);
        out<<",\"donor_role_coverage\":";number(out,node.donor_role_coverage);
        out<<",\"acceptor_role_coverage\":";number(out,node.acceptor_role_coverage);
        out<<",\"character\":";json_string(out,node.character);
        out<<",\"primary\":"<<(node.primary?"true":"false")<<'}';
    }
    out<<"],\"matched_edge_ids\":[";for(std::size_t i=0;i<v.matched_edge_ids.size();++i){if(i)out<<',';json_string(out,v.matched_edge_ids[i]);}
    out<<"]}";return out.str();
}
}
