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

#include "cov/pi_pair_evidence.hpp"
#include <iostream>
#include <limits>
int main(){
    using namespace cov;
    for(int scenario=0;scenario<7;++scenario)for(int prior=0;prior<5;++prior){
        PiPartnerComponents a{-.5,.1,.8,.9,.1},b{0,.8,.1,.9,-.1};
        if(scenario==1){a.metal_d=.8;a.ligand_p=.1;b.metal_d=.1;b.ligand_p=.8;}
        if(scenario==2){b.energy_hartree=-.495;a.metal_ligand_overlap=.001;b.metal_ligand_overlap=.002;}
        if(scenario==3)a.pi_fraction=.59;
        if(scenario==4){b.metal_d=.1;b.ligand_p=.8;}
        if(scenario==5)b.energy_hartree=2;
        if(scenario==6)a.pi_fraction=std::numeric_limits<double>::quiet_NaN();
        std::cout<<"{\"scenario\":"<<scenario<<",\"prior\":"<<prior<<",\"assessment\":"
            <<pi_partner_assessment_json(assess_pi_partner(a,b,static_cast<LigandPiPrior>(prior)))<<"}\n";
    }
}
