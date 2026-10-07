#pragma once

#include "cov/mo_diagram.hpp"
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <locale>
#include <sstream>
#include <string>
#include <utility>

namespace cov {
namespace mo_display_json_detail {
inline std::string quote(const std::string& value) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << '"';
    for (const unsigned char c:value) {
        switch(c) {
        case '"':out << "\\\"";break;
        case '\\':out << "\\\\";break;
        case '\n':out << "\\n";break;
        case '\r':out << "\\r";break;
        case '\t':out << "\\t";break;
        default:
            if(c<32)out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
            else out << c;
        }
    }
    out << '"';
    return out.str();
}
inline void number(std::ostream& out,double value) {
    if(std::isfinite(value))out << std::setprecision(17) << value;
    else out << "null";
}
}

// The same complete-group ledger is used by the view, copy/export and native
// display capture. Subitems never enter the mutually exclusive total twice.
inline std::string mo_group_composition_json(const MOGroupCompositionLedger& c) {
    using namespace mo_display_json_detail;
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "{\"available\":" << (c.available?"true":"false")
        << ",\"complete\":" << (c.complete?"true":"false")
        << ",\"status\":" << quote(c.status) << ",\"detail\":" << quote(c.detail)
        << ",\"source\":" << quote(c.source) << ",\"member_count\":" << c.member_count;
    for(const auto& [key,value]:std::initializer_list<std::pair<const char*,double>>{
        {"weight_sum",c.weight_sum},{"normalization_error",c.normalization_error},
        {"centre_current_s",c.centre_current_s},{"centre_current_p",c.centre_current_p},
        {"centre_current_d",c.centre_current_d},{"centre_current_f",c.centre_current_f},
        {"centre_other",c.centre_other},{"ligand_valence",c.ligand_valence},
        {"ligand_other",c.ligand_other},{"other_atoms",c.other_atoms},{"core",c.core},{"unresolved",c.unresolved},
        {"ligand_valence_s",c.ligand_valence_s},{"ligand_valence_p",c.ligand_valence_p}}) {
        out << ',' << quote(key) << ':'; number(out,value);
    }
    out << '}';return out.str();
}
inline std::string mo_group_display_decision_json(const MOGroupDisplayDecision& d) {
    using namespace mo_display_json_detail;
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "{\"included\":" << (d.included?"true":"false")
        << ",\"energy_window\":" << (d.energy_window?"true":"false")
        << ",\"major_relation\":" << (d.major_relation?"true":"false")
        << ",\"frontier\":" << (d.frontier?"true":"false")
        << ",\"coverage\":";number(out,d.coverage);
    out << ",\"sigma_coverage\":";number(out,d.sigma_coverage);
    out << ",\"reason_codes\":[";
    for(std::size_t i=0;i<d.reason_codes.size();++i){if(i)out<<',';out<<quote(d.reason_codes[i]);}
    out << "]}";return out.str();
}
inline std::string mo_sigma_framework_json(const MOSigmaFramework& s) {
    using namespace mo_display_json_detail;
    std::ostringstream out;out.imbue(std::locale::classic());
    out << "{\"available\":" << (s.available?"true":"false")
        << ",\"status\":" << quote(s.status) << ",\"detail\":" << quote(s.detail)
        << ",\"source\":" << quote(s.source) << ",\"rank\":" << s.rank
        << ",\"alpha_source_rank\":" << s.alpha_source_rank
        << ",\"beta_source_rank\":" << s.beta_source_rank
        << ",\"total_source_rank\":" << s.total_source_rank
        << ",\"shared_spatial_average\":" << (s.shared_spatial_average?"true":"false")
        << ",\"display_member_budget\":" << s.display_member_budget
        << ",\"retained_members\":" << s.retained_members;
    for(const auto& [key,value]:std::initializer_list<std::pair<const char*,double>>{
        {"source_orthogonality_error",s.source_orthogonality_error},
        {"source_canonical_closure_error",s.source_canonical_closure_error},{"occupied_trace",s.occupied_trace},
        {"retained_occupied_trace",s.retained_occupied_trace},{"trace_target",s.trace_target}}) {
        out << ',' << quote(key) << ':';number(out,value);
    }
    out << ",\"sources\":[";
    for(std::size_t i=0;i<s.sources.size();++i) {
        if(i)out << ',';const auto& q=s.sources[i];
        out << "{\"source_index\":" << q.orbital.index << ",\"spin\":" << quote(nbo_spin_name(q.orbital.spin))
            << ",\"centre\":" << q.centre << ",\"ligand\":" << q.ligand << ",\"kind\":" << quote(q.kind);
        out << ",\"occupation\":";number(out,q.occupation);
        out << ",\"ligand_valence_sp_weight\":";number(out,q.ligand_valence_sp_weight);
        out << ",\"axis_sigma_fraction\":";number(out,q.axis_sigma_fraction);
        out << '}';
    }
    out << "]}";return out.str();
}
} // namespace cov
