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

#include "cov/local_orbital_symmetry.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <span>
#include <string_view>
#include <vector>

namespace {

using Mat3 = std::array<double, 9>;

void require(const bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

constexpr double kInvSqrt2 = 0.70710678118654752440;
constexpr double kInvSqrt6 = 0.40824829046386301637;
constexpr std::array<Mat3, 5> d_basis{{
    Mat3{-kInvSqrt6, 0.0, 0.0, 0.0, -kInvSqrt6, 0.0, 0.0, 0.0, 2.0*kInvSqrt6},
    Mat3{0.0, 0.0, kInvSqrt2, 0.0, 0.0, 0.0, kInvSqrt2, 0.0, 0.0},
    Mat3{0.0, 0.0, 0.0, 0.0, 0.0, kInvSqrt2, 0.0, kInvSqrt2, 0.0},
    Mat3{kInvSqrt2, 0.0, 0.0, 0.0, -kInvSqrt2, 0.0, 0.0, 0.0, 0.0},
    Mat3{0.0, kInvSqrt2, 0.0, kInvSqrt2, 0.0, 0.0, 0.0, 0.0, 0.0},
}};

Mat3 multiply(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (std::size_t r=0;r<3;++r) for (std::size_t c=0;c<3;++c)
        for (std::size_t k=0;k<3;++k) out[3*r+c]+=a[3*r+k]*b[3*k+c];
    return out;
}

Mat3 transpose(const Mat3& a) {
    return {a[0],a[3],a[6],a[1],a[4],a[7],a[2],a[5],a[8]};
}

double dot(const Mat3& a,const Mat3& b) {
    double out=0.0; for (std::size_t i=0;i<9;++i) out+=a[i]*b[i]; return out;
}

Mat3 rotation_xyz(const double x,const double y,const double z) {
    const double cx=std::cos(x), sx=std::sin(x);
    const double cy=std::cos(y), sy=std::sin(y);
    const double cz=std::cos(z), sz=std::sin(z);
    const Mat3 rx{1,0,0,0,cx,-sx,0,sx,cx};
    const Mat3 ry{cy,0,sy,0,1,0,-sy,0,cy};
    const Mat3 rz{cz,-sz,0,sz,cz,0,0,0,1};
    return multiply(rz,multiply(ry,rx));
}

cov::Wavefunction base_wavefunction(const bool pure_d=true) {
    cov::Wavefunction wf;
    wf.atoms.push_back({"M",24,0,0,0});
    wf.shells.push_back({0,0,0,0,0,0});
    wf.shells.push_back({0,0,0,1,1,1});
    wf.shells.push_back({0,0,0,4,2,static_cast<std::uint8_t>(pure_d)});
    wf.basis_count=pure_d?9u:10u;
    // Controlled angular metric: pure functions are orthonormal. Cartesian d
    // diagonal functions have overlap 1/3; their scalar radial direction is
    // independent of the separate s reference in this mathematical fixture.
    wf.ao_overlap.assign(wf.basis_count*wf.basis_count,0);
    for(std::size_t i=0;i<wf.basis_count;++i)wf.ao_overlap[i*wf.basis_count+i]=1;
    if(!pure_d)for(std::size_t i=4;i<7;++i)for(std::size_t j=4;j<7;++j)
        if(i!=j)wf.ao_overlap[i*wf.basis_count+j]=1.0/3.0;
    return wf;
}

std::vector<double> rotated_p_coefficients(const Mat3& rotation,
                                          const std::array<double,3>& local) {
    std::array<double,3> input{};
    for (std::size_t r=0;r<3;++r) for (std::size_t c=0;c<3;++c)
        input[r]+=rotation[3*r+c]*local[c];
    std::vector<double> coefficients(9,0.0f);
    coefficients[1]=static_cast<float>(input[2]);
    coefficients[2]=static_cast<float>(-input[0]);
    coefficients[3]=static_cast<float>(-input[1]);
    return coefficients;
}

std::vector<double> rotated_pure_d_coefficients(const Mat3& rotation,
                                                const std::size_t component) {
    const Mat3 input=multiply(multiply(rotation,d_basis[component]),transpose(rotation));
    std::vector<double> coefficients(9,0.0f);
    constexpr std::array<double,5> signs{1.0,-1.0,-1.0,1.0,1.0};
    for (std::size_t i=0;i<5;++i)
        coefficients[4+i]=static_cast<float>(signs[i]*dot(input,d_basis[i]));
    return coefficients;
}

std::vector<double> rotated_cartesian_d_coefficients(const Mat3& rotation,
                                                     const std::size_t component) {
    const Mat3 input=multiply(multiply(rotation,d_basis[component]),transpose(rotation));
    constexpr double sqrt3=1.7320508075688772935;
    std::vector<double> coefficients(10,0.0f);
    coefficients[4]=static_cast<float>(sqrt3*input[0]);
    coefficients[5]=static_cast<float>(sqrt3*input[4]);
    coefficients[6]=static_cast<float>(sqrt3*input[8]);
    coefficients[7]=static_cast<float>(2.0*input[1]);
    coefficients[8]=static_cast<float>(2.0*input[2]);
    coefficients[9]=static_cast<float>(2.0*input[5]);
    return coefficients;
}

std::size_t add_orbital(cov::Wavefunction& wf,std::vector<double> coefficients) {
    cov::MolecularOrbital mo;
    mo.coefficients=std::move(coefficients);
    wf.orbitals.push_back(std::move(mo));
    return wf.orbitals.size()-1;
}

void expect(const cov::Wavefunction& wf,
            std::initializer_list<std::size_t> indices,
            std::string_view group,const Mat3& rotation,
            std::string_view label,cov::MetalAOShell shell,
            const std::uint8_t copy=1) {
    const std::vector<std::size_t> owned(indices);
    const auto assignment=cov::classify_local_metal_irrep(
        wf,owned,0,group,rotation);
    require(assignment.has_value(),"expected local irrep assignment");
    require(assignment->label==label,"unexpected local irrep label");
    require(assignment->shell==shell,"unexpected dominant AO shell");
    require(assignment->copy_index==copy,"irrep copy index was merged or changed");
    require(assignment->confidence>0.99999,"unexpectedly weak assignment confidence");
    require(assignment->source==cov::LocalIrrepSource::MetricAngularProjection && assignment->projection,
            "metric result must retain its source and scoped evidence");
}

} // namespace

int main() {
    const Mat3 r1=rotation_xyz(0.37,-0.51,0.22);
    const Mat3 r2=rotation_xyz(-0.43,0.29,0.61);

    // D4h: a rotated pure dz2 remains the first d-shell A1g copy.
    {
        auto wf=base_wavefunction();
        add_orbital(wf,rotated_pure_d_coefficients(r1,0));
        expect(wf,{0},"D4h",r1,"A1g",cov::MetalAOShell::D);
    }

    // C2v preserves the two independent A1 copies instead of string-merging.
    {
        auto wf=base_wavefunction();
        add_orbital(wf,rotated_pure_d_coefficients(r2,3));
        expect(wf,{0},"C2v",r2,"A1",cov::MetalAOShell::D,2);
    }

    // A complete rotated D3h E' subspace is aggregated over two orbitals.
    {
        auto wf=base_wavefunction();
        const auto a=add_orbital(wf,rotated_pure_d_coefficients(r1,3));
        const auto b=add_orbital(wf,rotated_pure_d_coefficients(r1,4));
        expect(wf,{a,b},"D3h",r1,"E'",cov::MetalAOShell::D);
    }

    // D4d E3 checks a different degenerate d subspace and orientation.
    {
        auto wf=base_wavefunction();
        add_orbital(wf,rotated_pure_d_coefficients(r2,1));
        expect(wf,{0},"D4d",r2,"E3",cov::MetalAOShell::D);
    }

    // A traceless Cartesian d example remains d; scalar directions are separately retained.
    {
        auto wf=base_wavefunction(false);
        add_orbital(wf,rotated_cartesian_d_coefficients(r1,4));
        expect(wf,{0},"D4h",r1,"B2g",cov::MetalAOShell::D);
    }

    // Pure p uses Gaussian's pz,-px,-py ordering before local rotation.
    {
        auto wf=base_wavefunction();
        add_orbital(wf,rotated_p_coefficients(r2,{0.0,0.0,1.0}));
        expect(wf,{0},"D3h",r2,"A2''",cov::MetalAOShell::P);
    }

    // Equivalent C2v A1 copies may mix without destroying the A1 irrep.
    {
        auto wf=base_wavefunction();
        auto coefficients=rotated_pure_d_coefficients(r1,0);
        const auto second=rotated_pure_d_coefficients(r1,3);
        for (std::size_t i=4;i<9;++i) coefficients[i]+=second[i];
        add_orbital(wf,std::move(coefficients));
        const std::array<std::size_t,1> index{0};
        const auto assignment=cov::classify_local_metal_irrep(
            wf,index,0,"C2v",r1);
        require(assignment && assignment->label=="A1" &&
                    assignment->copy_index==0u,
                "mixed equivalent copies must retain their irrep label");
    }

    // If the dominant p evidence is internally ambiguous, a clean d signal
    // may be used at 2% or above, but not below the transition-metal d floor.
    {
        auto wf=base_wavefunction();
        auto coefficients=rotated_p_coefficients(
            r1,{kInvSqrt2,kInvSqrt2,0.0});
        const auto d=rotated_pure_d_coefficients(r1,0);
        for (std::size_t i=4;i<9;++i)
            coefficients[i]+=static_cast<float>(std::sqrt(0.03))*d[i];
        add_orbital(wf,std::move(coefficients));
        expect(wf,{0},"C2v",r1,"A1",cov::MetalAOShell::D,1);
    }
    {
        auto wf=base_wavefunction();
        auto coefficients=rotated_p_coefficients(
            r1,{kInvSqrt2,kInvSqrt2,0.0});
        const auto d=rotated_pure_d_coefficients(r1,0);
        for (std::size_t i=4;i<9;++i)
            coefficients[i]+=0.1f*d[i];
        add_orbital(wf,std::move(coefficients));
        const std::array<std::size_t,1> index{0};
        require(!cov::classify_local_metal_irrep(wf,index,0,"C2v",r1),
                "sub-threshold fallback d evidence must be ignored");
    }

    // Invalid input is rejected rather than producing a misleading label.
    {
        auto wf=base_wavefunction();
        add_orbital(wf,rotated_pure_d_coefficients(r1,0));
        const std::array<std::size_t,1> index{0};
        Mat3 invalid=r1;
        invalid[0]*=2.0;
        require(!cov::classify_local_metal_irrep(wf,index,0,"D4h",invalid),
                "non-orthogonal frame must be rejected");
        require(!cov::classify_local_metal_irrep(wf,index,1,"D4h",r1),
                "invalid metal atom must be rejected");
        require(!cov::classify_local_metal_irrep(wf,index,0,"missing",r1),
                "unknown point group must be rejected");
    }

    return 0;
}
