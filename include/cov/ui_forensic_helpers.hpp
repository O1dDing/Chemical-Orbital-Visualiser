#pragma once

#include "cov/validation.hpp"
#include "cov/chemistry_route.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <locale>
#include <optional>
#include <sstream>

namespace cov::ui::forensic {
inline std::string number(double value) {
    if(!std::isfinite(value))return "null";
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<std::setprecision(17)<<value;return out.str();
}
inline std::string number(const std::optional<double>& value) {
    return value?number(*value):"null";
}
inline std::string index(const std::optional<std::size_t>& value) {
    return value?std::to_string(*value):"null";
}
inline const char* boolean(bool value) {return value?"true":"false";}
inline std::string rect(ImVec2 lo,ImVec2 hi) {
    return "["+number(lo.x)+","+number(lo.y)+","+number(hi.x)+","+number(hi.y)+"]";
}
template<class Range> inline std::string indices(const Range& values) {
    std::string out="[";bool first=true;
    for(const auto value:values){if(!first)out+=',';first=false;out+=std::to_string(value);}
    return out+"]";
}
inline std::string ref(const NboOrbitalRef& orbital) {
    return "{\"kind\":"+validation::quote(nbo_orbital_kind_name(orbital.kind))+
        ",\"spin\":"+validation::quote(nbo_spin_name(orbital.spin))+
        ",\"index\":"+std::to_string(orbital.index)+"}";
}
// Selection identity is small and bounded: never serialize coefficient vectors,
// spin-density/Fock matrices, full subspaces, or display-name provenance here.
inline std::string selection(const std::optional<NboOrbitalSelection>& value) {
    if(!value)return "null";
    const auto& s=*value;
    std::string out="{\"source_id\":"+validation::quote(s.source_id)+
        ",\"dataset_id\":"+validation::quote(s.dataset_id)+
        ",\"group_id\":"+validation::quote(s.group_id)+
        ",\"semantic_kind\":"+validation::quote(s.semantic_kind)+
        ",\"mode\":"+std::to_string(static_cast<int>(s.mode))+
        ",\"target_canonical_index\":"+index(s.target_canonical_index)+
        ",\"term_count\":"+std::to_string(s.terms.size())+",\"terms\":[";
    const auto count=std::min<std::size_t>(128,s.terms.size());
    for(std::size_t i=0;i<count;++i){if(i)out+=',';
        out+="{\"orbital\":"+ref(s.terms[i].orbital)+
            ",\"coefficient\":"+number(s.terms[i].coefficient)+"}";}
    return out+"],\"terms_truncated\":"+boolean(count<s.terms.size())+"}";
}
inline std::string active(const std::optional<ActiveOrbitalView>& value) {
    if(!value)return "null";
    const auto& v=*value;
    return "{\"kind\":"+validation::quote(active_orbital_kind_name(v.kind))+
        ",\"source_id\":"+validation::quote(v.source_id)+
        ",\"group_id\":"+validation::quote(v.group_id)+
        ",\"semantic_kind\":"+validation::quote(v.semantic_kind)+
        ",\"spin\":"+validation::quote(nbo_spin_name(v.spin))+
        ",\"source_spin\":"+validation::quote(nbo_spin_name(v.source_spin))+
        ",\"canonical_index\":"+index(v.canonical_index)+
        ",\"rendered_index\":"+index(v.rendered_index)+
        ",\"selection\":"+selection(v.selection)+"}";
}
} // namespace cov::ui::forensic
