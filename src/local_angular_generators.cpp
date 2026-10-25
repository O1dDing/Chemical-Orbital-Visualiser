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
#include "cov/ao_angular_basis.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace cov {
namespace {
using Powers = std::array<int,3>;
using Polynomial = std::map<Powers,double>;
using ao_angular::cartesian_components;
using ao_angular::odd_factorial;
Polynomial solid_harmonic(int l,int component) {
    Polynomial result;
    for (const auto& term:ao_angular::solid_harmonic(l,component))
        result.emplace(term.powers,term.coefficient);
    return result;
}
Polynomial multiply(const Polynomial& a,const Polynomial& b) {
    Polynomial out;
    for(const auto& [p,x]:a)for(const auto& [q,y]:b)
        out[{p[0]+q[0],p[1]+q[1],p[2]+q[2]}]+=x*y;
    return out;
}
Polynomial rotated(const Polynomial& p,const std::array<double,9>& r) {
    std::array<Polynomial,3> coordinates;
    for(int local=0;local<3;++local)for(int input=0;input<3;++input) {
        Powers power{};power[input]=1;
        coordinates[local][power]=r[3*input+local];
    }
    Polynomial out;
    for(const auto& [powers,coefficient]:p) {
        Polynomial term{{Powers{0,0,0},coefficient}};
        for(int axis=0;axis<3;++axis)for(int n=0;n<powers[axis];++n)
            term=multiply(term,coordinates[axis]);
        for(const auto& [power,value]:term)out[power]+=value;
    }
    return out;
}
double sphere_inner_product(const Polynomial& a,const Polynomial& b) {
    double sum=0;
    for(const auto& [p,x]:a)for(const auto& [q,y]:b) {
        const Powers powers{p[0]+q[0],p[1]+q[1],p[2]+q[2]};
        if(powers[0]%2 || powers[1]%2 || powers[2]%2)continue;
        const double moment=2*std::tgamma((powers[0]+1)*.5)*std::tgamma((powers[1]+1)*.5)*
            std::tgamma((powers[2]+1)*.5)/std::tgamma((powers[0]+powers[1]+powers[2]+3)*.5);
        sum+=x*y*moment;
    }
    return sum;
}
double coefficient_scale(const Powers& p) {
    return std::sqrt(odd_factorial(2*p[0]-1)*odd_factorial(2*p[1]-1)*odd_factorial(2*p[2]-1));
}
}
LocalAngularGenerators local_angular_generators(int degree,bool pure,int angular,
    const std::array<double,9>& rotation) {
    LocalAngularGenerators out;
    out.shell_degree=degree;out.angular_degree=angular;out.pure_shell=pure;
    if(degree<0 || degree>4 || angular<0 || angular>degree || (degree-angular)%2 || (pure&&angular!=degree)) {
        out.detail="Requested angular subspace is absent or outside the supported shell definition";
        return out;
    }
    for(double v:rotation)if(!std::isfinite(v)) {out.detail="Nonfinite local frame";return out;}
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        double dot=0;for(int k=0;k<3;++k)dot+=rotation[3*k+i]*rotation[3*k+j];
        out.rotation_orthogonality_error=std::max(out.rotation_orthogonality_error,std::abs(dot-(i==j?1.0:0.0)));
    }
    if(out.rotation_orthogonality_error>1e-10) {out.detail="Local frame is not orthogonal";return out;}
    const auto powers=cartesian_components(degree);
    out.basis_count=pure?static_cast<std::size_t>(2*degree+1):powers.size();
    out.columns=static_cast<std::size_t>(2*angular+1);
    out.coefficients.assign(out.basis_count*out.columns,0);
    std::vector<Polynomial> basis;
    for(std::size_t i=0;i<out.basis_count;++i) {
        if(pure)basis.push_back(solid_harmonic(degree,static_cast<int>(i)));
        else basis.push_back(Polynomial{{powers[i],1/coefficient_scale(powers[i])}});
    }
    const Polynomial radius_squared{{Powers{2,0,0},1},{Powers{0,2,0},1},{Powers{0,0,2},1}};
    for(std::size_t column=0;column<out.columns;++column) {
        auto target=rotated(solid_harmonic(angular,static_cast<int>(column)),rotation);
        for(int k=0;k<(degree-angular)/2;++k)target=multiply(target,radius_squared);
        Polynomial reconstruction;
        for(std::size_t row=0;row<out.basis_count;++row) {
            double value=0;
            if(pure)value=sphere_inner_product(basis[row],target);
            else {
                const auto found=target.find(powers[row]);
                if(found!=target.end())value=found->second*coefficient_scale(powers[row]);
            }
            out.coefficients[row*out.columns+column]=value;
            for(const auto& [p,c]:basis[row])reconstruction[p]+=value*c;
        }
        for(const auto& [p,c]:target)reconstruction[p]-=c;
        for(const auto& [p,c]:reconstruction)out.polynomial_reconstruction_error=std::max(out.polynomial_reconstruction_error,std::abs(c));
    }
    out.available=std::isfinite(out.polynomial_reconstruction_error) && out.polynomial_reconstruction_error<=1e-10;
    out.detail=out.available?"Local angular generators in COV AO order and normalization; radial envelope remains that of each source shell"
                            :"Angular polynomial reconstruction is unresolved";
    return out;
}
}
