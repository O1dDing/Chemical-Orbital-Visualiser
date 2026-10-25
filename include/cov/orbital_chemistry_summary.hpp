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
