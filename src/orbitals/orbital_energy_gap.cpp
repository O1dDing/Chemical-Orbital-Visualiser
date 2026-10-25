/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#include "cov/mo_diagram.hpp"
#include "cov/orbital_symmetry_scope.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>

namespace cov {
std::string pi_partner_candidates_json(const std::vector<PiPartnerAssessment>& candidates) {
    std::ostringstream out;out << '[';
    for(std::size_t i=0;i<candidates.size();++i) {
        if(i)out << ',';
        out << pi_partner_assessment_json(candidates[i]);
    }
    out << ']';return out.str();
}
namespace {
void number(std::ostream& out, double value) {
    if (std::isfinite(value)) out << value;
    else out << "null";
}
void quoted(std::ostream& out, const std::string& text) {
    static constexpr char hex[]="0123456789abcdef";
    out << '"';
    for (unsigned char c:text) {
        if (c=='"' || c=='\\') out << '\\' << static_cast<char>(c);
        else if (c<0x20) out << "\\u00" << hex[c>>4] << hex[c&15];
        else out << static_cast<char>(c);
    }
    out << '"';
}
void indices(std::ostream& out, const std::vector<std::size_t>& values) {
    out << '[';
    for (std::size_t i=0;i<values.size();++i) {if(i) out << ',';out << values[i];}
    out << ']';
}
const char* pi_kind_id(PiInteractionKind kind) {
    switch(kind) {
        case PiInteractionKind::Donor:return "donor";
        case PiInteractionKind::Acceptor:return "acceptor";
        case PiInteractionKind::WeakNearNonbonding:return "weak-near-nonbonding";
        default:return "coupled";
    }
}
}
const char* orbital_energy_gap_kind_name(OrbitalEnergyGapKind kind) noexcept {
    return kind==OrbitalEnergyGapKind::CrystalField?"crystal-field":"pi-partner";
}
std::vector<const OrbitalEnergyGapDescriptor*> orbital_energy_gaps(const MODiagramData& data) {
    std::vector<const OrbitalEnergyGapDescriptor*> result;
    result.reserve(data.pi_interactions.size()+data.crystal_field_gaps.size());
    for(const auto& gap:data.pi_interactions)result.push_back(&gap);
    for(const auto& gap:data.crystal_field_gaps)result.push_back(&gap);
    return result;
}
std::string orbital_energy_gap_json(const OrbitalEnergyGapDescriptor& gap, EnergyUnit unit) {
    const bool cf=gap.gap_kind==OrbitalEnergyGapKind::CrystalField;
    std::ostringstream out;out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "{\"gap_kind\":\"" << orbital_energy_gap_kind_name(gap.gap_kind)
        << "\",\"kind\": \"" << (cf?"weak-crystal-field":pi_kind_id(gap.kind))
        << "\",\"description\":";
    quoted(out,cf?"Local d-manifold crystal-field separation":pi_interaction_kind_name(gap.kind));
    out << ",\"interpretation_scope\":\"local-metal-ligand-subspace\",\"symmetry\":";
    quoted(out,gap.symmetry);
    out << ",\"lower_level\":" << gap.lower_level << ",\"upper_level\":" << gap.upper_level
        << ",\"retained_level\":" << gap.retained_level << ",\"lower_orbitals\": ";
    indices(out,gap.lower_orbitals);out << ",\"upper_orbitals\":";indices(out,gap.upper_orbitals);
    out << ",\"lower_visible\":" << (gap.lower_visible?"true":"false")
        << ",\"upper_visible\":" << (gap.upper_visible?"true":"false");
    out << ",\"splitting_hartree\":";number(out,gap.splitting_hartree);
    out << ",\"splitting_display\":";number(out,convert_hartree(gap.splitting_hartree,unit));
    out << ",\"display_unit\":";quoted(out,energy_unit_symbol(unit));
    out << ",\"confidence\":";number(out,gap.confidence);
    const bool verified_channel=gap.orbital_evidence&&gap.orbital_evidence->channel.same_operator_verified;
    out << ",\"score_meaning\":";quoted(out,verified_channel?
        "minimum-endpoint-cross-Fock-magnitude-not-probability":"heuristic-support-not-probability");
    out << ",\"score_unit\":";quoted(out,verified_channel?"hartree":"dimensionless");
    out << ",\"ranking_score\":";
    if(verified_channel)number(out,gap.orbital_evidence->ranking_score);else out<<"null";
    out << ",\"display_strength_hartree\":";number(out,gap.display_strength_hartree);
    out << ",\"display_strength_definition\":";
    quoted(out,"minimum-absolute-cross-Fock-trace-of-complete-endpoint-groups; display-priority-only");
    out << ",\"equivalent_channel_ids\":[";
    for(std::size_t i=0;i<gap.equivalent_channel_ids.size();++i){if(i)out<<',';quoted(out,gap.equivalent_channel_ids[i]);}
    out << ']';
    out << ",\"lower_energy_hartree\":";number(out,gap.lower_energy_hartree);
    out << ",\"upper_energy_hartree\":";number(out,gap.upper_energy_hartree);
    out << ",\"lower_energy_spread_hartree\":";number(out,gap.lower_energy_spread_hartree);
    out << ",\"upper_energy_spread_hartree\":";number(out,gap.upper_energy_spread_hartree);
    out << ",\"weak_split_threshold_hartree\":";number(out,gap.weak_split_threshold_hartree);
    out << ",\"weak_overlap_threshold\":";number(out,gap.weak_overlap_threshold);
    out << ",\"lower_symmetry_scope\":" << orbital_symmetry_json(gap.lower_symmetry_scope)
        << ",\"upper_symmetry_scope\":" << orbital_symmetry_json(gap.upper_symmetry_scope);
    out << ",\"orbital_evidence\":";
    if(gap.orbital_evidence && !cf)out << pi_partner_assessment_json(*gap.orbital_evidence);
    else out << "null";
    out << ",\"crystal_field_evidence\":";
    if(gap.crystal_field_evidence && cf)out << weak_crystal_field_assessment_json(*gap.crystal_field_evidence);
    else out << "null";
    if(gap.crystal_field_evidence && cf)
        out << ",\"crystal_field_component_order\":\"first=lower;second=upper\"";
    out << '}';return out.str();
}
std::string orbital_energy_gap_array_json(const std::vector<OrbitalEnergyGapDescriptor>& gaps, EnergyUnit unit) {
    std::ostringstream out;out << '[';
    for(std::size_t i=0;i<gaps.size();++i){if(i)out << ',';out << orbital_energy_gap_json(gaps[i],unit);}
    out << ']';return out.str();
}
}
