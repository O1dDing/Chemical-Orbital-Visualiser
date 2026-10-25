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

#include "cov/local_angular_generators.hpp"
#include <iomanip>
#include <iostream>
int main() {
    std::cout<<std::setprecision(17);
    int id=0,L=0,pure=0,l=0;
    while(std::cin>>id>>L>>pure>>l) {
        std::array<double,9> rotation{};
        for(auto& v:rotation)if(!(std::cin>>v))return 2;
        const auto x=cov::local_angular_generators(L,pure!=0,l,rotation);
        std::cout<<"{\"id\":"<<id<<",\"available\":"<<(x.available?"true":"false")
            <<",\"basis_count\":"<<x.basis_count<<",\"columns\":"<<x.columns
            <<",\"rotation_error\":"<<x.rotation_orthogonality_error
            <<",\"polynomial_error\":"<<x.polynomial_reconstruction_error<<",\"coefficients\":[";
        bool first=true;for(double v:x.coefficients){if(!first)std::cout<<',';first=false;std::cout<<v;}
        std::cout<<"]}\n";
    }
    return std::cin.eof()?0:2;
}
