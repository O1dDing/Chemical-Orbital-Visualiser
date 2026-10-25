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

// Numerical (not screenshot/UI) comparison of the production CUDA evaluator
// against an independently generated, single-orbital Gaussian cube.
#include "cov/cuda_orbital.hpp"
#include "cov/gl_api.hpp"
#include "cov/wavefunction_io.hpp"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr GLenum texture3d = 0x806F, red = 0x1903, r32f = 0x822E;
struct Cube {
    std::array<int,3> n{};
    std::array<double,3> origin{}, step{};
    std::vector<cov::Atom> atoms;
    std::vector<double> values;
    int orbital = -1;
};
Cube read_cube(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read reference cube");
    std::string line;
    std::getline(input,line); std::getline(input,line);
    Cube cube;
    int natoms;
    input >> natoms >> cube.origin[0] >> cube.origin[1] >> cube.origin[2];
    std::getline(input,line);
    for(int axis=0;axis<3;++axis) {
        double vector[3];
        input >> cube.n[axis] >> vector[0] >> vector[1] >> vector[2];
        std::getline(input,line);
        if(cube.n[axis]<2 || cube.n[axis]>512)
            throw std::runtime_error("Require positive Bohr cube axes, 2..512 points");
        for(int j=0;j<3;++j)
            if(j!=axis && std::abs(vector[j])>1e-12)
                throw std::runtime_error("Reference grid must be axis aligned");
        cube.step[axis]=vector[axis];
        if(cube.step[axis]<=0) throw std::runtime_error("Nonpositive grid step");
    }
    for(int i=0;i<std::abs(natoms);++i) {
        cov::Atom atom;
        input >> atom.atomic_number >> atom.nuclear_charge >> atom.x >> atom.y >> atom.z;
        cube.atoms.push_back(atom);
    }
    if(natoms<0) {
        int orbitals;
        input >> orbitals >> cube.orbital;
        if(orbitals!=1) throw std::runtime_error("Require one orbital per cube");
    } else throw std::runtime_error("Require orbital cube, not density cube");
    const std::size_t count=static_cast<std::size_t>(cube.n[0])*cube.n[1]*cube.n[2];
    cube.values.resize(count);
    for(double& value:cube.values)
        if(!(input>>value) || !std::isfinite(value))
            throw std::runtime_error("Invalid/truncated cube values");
    double extra;
    if(input>>extra) throw std::runtime_error("Extra cube values");
    return cube;
}
void identity(const Cube& cube,const cov::Wavefunction& wf) {
    if(cube.atoms.size()!=wf.atoms.size()) throw std::runtime_error("Atom count mismatch");
    for(std::size_t i=0;i<wf.atoms.size();++i) {
        const auto& a=cube.atoms[i]; const auto& b=wf.atoms[i];
        if(a.atomic_number!=b.atomic_number || std::abs(a.x-b.x)>6e-6 ||
           std::abs(a.y-b.y)>6e-6 || std::abs(a.z-b.z)>6e-6)
            throw std::runtime_error("Cube/FCHK atom identity or coordinate mismatch");
    }
}
}
int main(int argc,char** argv) {
    GLFWwindow* window=nullptr;
    GLuint texture=0;
    try {
        if(argc<4 || argc%2!=0)
            throw std::runtime_error("Usage: cov_cuda_cube_reference file.fch internal_index reference.cube [index cube ...]");
        cov::WavefunctionParseOptions options;
        options.max_atoms=1000;
        options.auto_enrich_gaussian_log=false;
        auto wf=cov::parse_wavefunction(argv[1],options);
        // Opt-in counterfactual to diagnose a phase convention mismatch.
        // The default remains the unmodified production path. Diagnostic
        // results must never be reported as a passing production build.
        const bool phase_diagnostic=std::getenv("COV_REFERENCE_DIAGNOSTIC_ODD_M")!=nullptr;
        if(phase_diagnostic) {
            for(const auto& shell:wf.shells) {
                if(!shell.pure || shell.angular_momentum<2) continue;
                for(unsigned component=1;component<cov::shell_basis_count(shell);++component) {
                    const unsigned m=(component+1)/2;
                    if(m%2==1)
                        for(auto& orbital:wf.orbitals)
                            orbital.coefficients[shell.basis_offset+component]*=-1.0f;
                }
            }
        }
        if(!glfwInit()) throw std::runtime_error("GLFW init failed");
        glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
        window=glfwCreateWindow(16,16,"COV numeric reference context",nullptr,nullptr);
        if(!window) throw std::runtime_error("Hidden GL context failed");
        glfwMakeContextCurrent(window);
        if(!cov::gl::load()) throw std::runtime_error("GL entry points unavailable");
        glGenTextures(1,&texture);
        bool all_pass=true;
        {
            cov::CudaOrbitalEvaluator evaluator(wf);
            for(int arg=2;arg<argc;arg+=2) {
                const std::size_t mo=std::stoull(argv[arg]);
                if(mo>=wf.orbitals.size()) throw std::runtime_error("MO index out of range");
                const Cube cube=read_cube(argv[arg+1]);
                identity(cube,wf);
                cov::GridBox box;
                box.min_x=static_cast<float>(cube.origin[0]);
                box.min_y=static_cast<float>(cube.origin[1]);
                box.min_z=static_cast<float>(cube.origin[2]);
                box.max_x=static_cast<float>(cube.origin[0]+(cube.n[0]-1)*cube.step[0]);
                box.max_y=static_cast<float>(cube.origin[1]+(cube.n[1]-1)*cube.step[1]);
                box.max_z=static_cast<float>(cube.origin[2]+(cube.n[2]-1)*cube.step[2]);
                glBindTexture(texture3d,texture);
                cov::gl::TexImage3D(texture3d,0,r32f,cube.n[0],cube.n[1],cube.n[2],0,red,GL_FLOAT,nullptr);
                if(glGetError()!=GL_NO_ERROR) throw std::runtime_error("GL texture allocation failed");
                evaluator.attach_gl_texture(texture);
                evaluator.evaluate(mo,box,cube.n[0],cube.n[1],cube.n[2]);
                std::vector<float> actual(cube.values.size());
                glBindTexture(texture3d,texture);
                glGetTexImage(texture3d,0,red,GL_FLOAT,actual.data());
                if(glGetError()!=GL_NO_ERROR) throw std::runtime_error("GL texture readback failed");
                evaluator.detach_gl_texture();
                long double dot=0, ref2=0, gpu2=0;
                double max_ref=0;
                for(int x=0;x<cube.n[0];++x) for(int y=0;y<cube.n[1];++y) for(int z=0;z<cube.n[2];++z) {
                    const auto ci=(static_cast<std::size_t>(x)*cube.n[1]+y)*cube.n[2]+z;
                    const auto gi=(static_cast<std::size_t>(z)*cube.n[1]+y)*cube.n[0]+x;
                    const double a=actual[gi],r=cube.values[ci];
                    if(!std::isfinite(a)) throw std::runtime_error("Nonfinite CUDA value");
                    dot+=a*r; ref2+=r*r; gpu2+=a*a; max_ref=std::max(max_ref,std::abs(r));
                }
                if(ref2<1e-20 || gpu2<1e-20) throw std::runtime_error("Grid has insufficient orbital signal");
                const double phase=dot<0?-1:1;
                long double error2=0; double max_error=0;
                for(int x=0;x<cube.n[0];++x) for(int y=0;y<cube.n[1];++y) for(int z=0;z<cube.n[2];++z) {
                    const auto ci=(static_cast<std::size_t>(x)*cube.n[1]+y)*cube.n[2]+z;
                    const auto gi=(static_cast<std::size_t>(z)*cube.n[1]+y)*cube.n[0]+x;
                    const double error=phase*actual[gi]-cube.values[ci];
                    error2+=error*error; max_error=std::max(max_error,std::abs(error));
                }
                const double nrms=std::sqrt(static_cast<double>(error2/ref2));
                const double cosine=std::abs(static_cast<double>(dot/std::sqrt(ref2*gpu2)));
                const bool pass=nrms<=1e-4 && cosine>=1-1e-7 && max_error/max_ref<=1e-3;
                all_pass=all_pass&&pass;
                std::cout<<std::setprecision(15)<<"{\"internal_index\":"<<mo
                    <<",\"phase_diagnostic\":"<<(phase_diagnostic?"true":"false")
                    <<",\"reference_orbital\":"<<cube.orbital<<",\"points\":"<<actual.size()
                    <<",\"nrms\":"<<nrms<<",\"cosine\":"<<cosine<<",\"phase\":"<<phase
                    <<",\"max_abs_error\":"<<max_error<<",\"max_reference\":"<<max_ref
                    <<",\"kernel_ms\":"<<evaluator.last_kernel_ms()
                    <<",\"pass\":"<<(pass?"true":"false")<<"}\n";
            }
        }
        glDeleteTextures(1,&texture); glfwDestroyWindow(window); glfwTerminate();
        return all_pass?0:2;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        if(window) glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }
}
