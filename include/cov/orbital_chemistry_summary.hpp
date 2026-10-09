#pragma once

#include "cov/model.hpp"

#include <cmath>
#include <sstream>
#include <string>

namespace cov {

// Shared, evidence-preserving text for the compact diagram hover and the
// detailed orbital panel. Fractions describe this canonical orbital's
// analysis; they are not probabilities or whole-bond populations.
inline std::string orbital_channel_fraction_summary(const OrbitalChannelDistribution& value) {
    if(value.status!=ChemistryStatus::Percentages)return {};
    std::ostringstream out;
    out<<" [σ "<<std::lround(100.0*value.sigma)
       <<"% · π "<<std::lround(100.0*value.pi)
       <<"% · δ "<<std::lround(100.0*value.delta)
       <<"% · φ "<<std::lround(100.0*value.phi)<<'%';
    if(value.undetermined>0.005)
        out<<" · UND "<<std::lround(100.0*value.undetermined)<<'%';
    return out.str()+']';
}

inline std::string orbital_bonding_fraction_summary(const OrbitalBondingDistribution& value,
    const char* bonding,const char* antibonding,const char* nonbonding) {
    if(value.status!=ChemistryStatus::Percentages)return {};
    std::ostringstream out;
    out<<" ["<<bonding<<' '<<std::lround(100.0*value.bonding)
       <<"% · "<<antibonding<<' '<<std::lround(100.0*value.antibonding)
       <<"% · "<<nonbonding<<' '<<std::lround(100.0*value.nonbonding)<<'%';
    if(value.undetermined>0.005)
        out<<" · UND "<<std::lround(100.0*value.undetermined)<<'%';
    return out.str()+']';
}

} // namespace cov
