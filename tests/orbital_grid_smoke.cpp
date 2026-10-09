#include "cov/orbital_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <stop_token>
#include <vector>

namespace {

void require(bool condition, const char* name) {
    if (!condition) throw std::runtime_error(name);
}

void near(double actual, double expected, const char* name) {
    if (std::abs(actual-expected) > 2.0e-5*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(name);
}

template <class Exception, class F> void rejects(F&& function, const char* name) {
    bool rejected = false;
    try { function(); } catch (const Exception&) { rejected = true; }
    require(rejected,name);
}

constexpr double pi = 3.141592653589793238462643383279502884;

double factorial_odd(int n) {
    double result = 1.0;
    for (; n > 1; n -= 2) result *= n;
    return result;
}

double pure_harmonic(int l, int c, double x, double y, double z) {
    const double x2=x*x, y2=y*y, z2=z*z, r2=x2+y2+z2;
    if (l == 0) return 0.28209479177387814;
    if (l == 1) {
        const double xyz[] = {z,-x,-y};
        return 0.4886025119029199*xyz[c];
    }
    if (l == 2) {
        const double values[] = {
            0.31539156525252005*(3*z2-r2),
            -1.0925484305920792*x*z,
            -1.0925484305920792*y*z,
            0.5462742152960396*(x2-y2),
            1.0925484305920792*x*y
        };
        return values[c];
    }
    if (l == 3) {
        const double values[] = {
            0.3731763325901154*z*(5*z2-3*r2),
            -0.4570457994644658*x*(5*z2-r2),
            -0.4570457994644658*y*(5*z2-r2),
            1.445305721320277*z*(x2-y2),
            2.890611442640554*x*y*z,
            -0.5900435899266435*x*(x2-3*y2),
            -0.5900435899266435*y*(3*x2-y2)
        };
        return values[c];
    }
    const double values[] = {
        0.10578554691520431*(35*z2*z2-30*z2*r2+3*r2*r2),
        -0.6690465435572892*x*z*(7*z2-3*r2),
        -0.6690465435572892*y*z*(7*z2-3*r2),
        0.47308734787878004*(x2-y2)*(7*z2-r2),
        0.9461746957575601*x*y*(7*z2-r2),
        -1.7701307697799304*x*z*(x2-3*y2),
        -1.7701307697799304*y*z*(3*x2-y2),
        0.6258357354491761*(x2*x2-6*x2*y2+y2*y2),
        2.5033429417967046*x*y*(x2-y2)
    };
    return values[c];
}

const std::array<std::vector<std::array<int,3>>,5> cartesian_components{{
    {{{0,0,0}}},
    {{{1,0,0},{0,1,0},{0,0,1}}},
    {{{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}}},
    {{{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
      {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}}},
    {{{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
      {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
      {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}}}
}};

cov::Wavefunction one_shell(int l, bool pure, int component) {
    cov::Wavefunction wf;
    cov::Atom atom;
    atom.x = 0.2; atom.y = -0.3; atom.z = 0.1;
    wf.atoms.push_back(atom);
    wf.primitives.push_back({1.25,-0.75});
    cov::Shell shell;
    shell.primitive_count = 1;
    shell.angular_momentum = static_cast<std::uint8_t>(l);
    shell.pure = pure ? 1 : 0;
    wf.shells.push_back(shell);
    wf.basis_count = cov::shell_basis_count(shell);
    cov::MolecularOrbital mo;
    mo.coefficients.assign(wf.basis_count,0.0);
    mo.coefficients[component] = 0.6;
    wf.orbitals.push_back(mo);
    return wf;
}

double expected(int l, bool pure, int component, double x, double y, double z) {
    const double alpha=1.25;
    const double r2=x*x+y*y+z*z;
    const double scale=-0.75*0.6*std::exp(-alpha*r2);
    if (pure) {
        const double norm=std::sqrt(2*std::pow(2*alpha,l+1.5)/std::tgamma(l+1.5));
        return scale*norm*pure_harmonic(l,component,x,y,z);
    }
    const auto powers=cartesian_components[l].at(component);
    const double norm=std::pow(2*alpha/pi,0.75)*std::pow(4*alpha,0.5*l)/
        std::sqrt(factorial_odd(2*powers[0]-1)*
                  factorial_odd(2*powers[1]-1)*
                  factorial_odd(2*powers[2]-1));
    return scale*norm*std::pow(x,powers[0])*
           std::pow(y,powers[1])*std::pow(z,powers[2]);
}

void numeric_components() {
    // Distinct nonzero offsets expose phase, exponent order and mixed powers.
    const double dx=0.43, dy=-0.37, dz=0.59;
    const cov::GridBox box{static_cast<float>(0.2+dx),static_cast<float>(-0.3+dy),
                           static_cast<float>(0.1+dz),3.0f,4.0f,5.0f};
    for (int l=0; l<=4; ++l) for (bool pure : {false,true}) {
        const int count=pure ? 2*l+1 : (l+1)*(l+2)/2;
        if (!pure) require(cartesian_components[l].size()==static_cast<std::size_t>(count),
                           "analytic Cartesian reference component count");
        for (int c=0; c<count; ++c) {
            const auto wf=one_shell(l,pure,c);
            const auto before=wf;
            const auto terms=cov::pack_orbital_terms(wf,0);
            require(!terms.empty(),"packed component has terms");
            require(wf.atoms[0].x == before.atoms[0].x &&
                    wf.primitives[0].coefficient == before.primitives[0].coefficient &&
                    wf.orbitals[0].coefficients == before.orbitals[0].coefficients,
                    "source scientific data stays unchanged");
            float actual=0.0f;
            cov::evaluate_cpu_grid(terms,cov::make_grid_request(box,1,1,1,0,1),
                                   std::span<float>(&actual,1));
            const double reference=expected(l,pure,c,dx,dy,dz);
            if (std::abs(static_cast<double>(actual)-reference) >
                2.0e-5*std::max(1.0,std::abs(reference))) {
                std::ostringstream message;
                message << "s through g component l=" << l << " pure=" << pure
                        << " c=" << c << " actual=" << actual
                        << " expected=" << reference;
                if (!pure) {
                    const auto p=cartesian_components[l].at(c);
                    message << " powers=" << p[0] << ',' << p[1] << ',' << p[2]
                            << " count=" << cartesian_components[l].size();
                }
                throw std::runtime_error(message.str());
            }
            if (pure && l > 0)
                for (const auto& term : terms)
                    require(term.ax+term.ay+term.az == static_cast<std::uint32_t>(l),
                            "pure harmonic degree");
        }
    }
}

void grid_chunks_and_bounds() {
    const auto wf=one_shell(2,true,3);
    const auto terms=cov::pack_orbital_terms(wf,0);
    const cov::GridBox box{-0.2f,0.1f,-0.6f,0.9f,0.1f,0.7f};
    constexpr int nx=4,ny=1,nz=3;
    require(cov::checked_grid_size(nx,ny,nz)==12,"noncube grid size");
    const auto full=cov::make_grid_request(box,nx,ny,nz,0,12);
    std::vector<float> all(12),chunks(12);
    cov::evaluate_cpu_grid(terms,full,all);
    for (std::size_t first=0; first<12; first+=5) {
        const auto count=std::min<std::size_t>(5,12-first);
        const auto part=cov::make_grid_request(box,nx,ny,nz,first,count);
        cov::evaluate_cpu_grid(terms,part,std::span<float>(chunks).subspan(first,count));
    }
    require(all==chunks,"split grid matches full grid");
    for (int z=0; z<nz; ++z) for (int x=0; x<nx; ++x) {
        const double gx=box.min_x+static_cast<float>(x)/3*(box.max_x-box.min_x)-0.2;
        const double gy=box.min_y-(-0.3);
        const double gz=box.min_z+static_cast<float>(z)/2*(box.max_z-box.min_z)-0.1;
        near(all[x+nx*z],expected(2,true,3,gx,gy,gz),"x fastest and dimension one");
    }
    rejects<std::invalid_argument>([]{ cov::checked_grid_size(0,2,2); },"zero dimension");
    rejects<std::length_error>([]{ cov::checked_grid_size(513,513,513); },"bounded volume");
    rejects<std::out_of_range>([&]{ cov::make_grid_request(box,nx,ny,nz,11,2); },"chunk range");
    auto invalid_box=box;
    invalid_box.max_x=std::numeric_limits<float>::infinity();
    rejects<std::invalid_argument>([&]{ cov::make_grid_request(invalid_box,1,1,1,0,1); },"nonfinite box");
    invalid_box=box;
    invalid_box.max_x=invalid_box.min_x-1.0f;
    rejects<std::invalid_argument>([&]{ cov::make_grid_request(invalid_box,1,1,1,0,1); },"reversed box");
    rejects<std::invalid_argument>([&]{ cov::evaluate_cpu_grid(terms,full,std::span<float>(all).first(11)); },"short output");
    cov::stop_source cancelled;
    cancelled.request_stop();
    rejects<std::runtime_error>([&]{ cov::evaluate_cpu_grid(terms,full,all,cancelled.get_token()); },"cancellation");
}

void malformed_wavefunction() {
    const auto valid=one_shell(1,false,0);
    auto check=[&](std::function<void(cov::Wavefunction&)> mutate,const char* name) {
        auto wf=valid;
        mutate(wf);
        rejects<std::exception>([&]{ cov::pack_orbital_terms(wf,0); },name);
    };
    rejects<std::out_of_range>([&]{ cov::pack_orbital_terms(valid,1); },"MO index");
    check([](auto& wf){ wf.basis_count++; },"MO dimension");
    check([](auto& wf){ wf.shells[0].angular_momentum=5; },"shell angular degree");
    check([](auto& wf){ wf.shells[0].pure=2; },"shell representation");
    check([](auto& wf){ wf.shells[0].atom_index=1; },"atom index");
    check([](auto& wf){ wf.shells[0].primitive_offset=1; },"primitive index");
    check([](auto& wf){ wf.shells[0].basis_offset=1; },"basis index");
    check([](auto& wf){ wf.primitives[0].exponent=0; },"zero exponent");
    check([](auto& wf){ wf.primitives[0].exponent=-1; },"negative exponent");
    check([](auto& wf){ wf.primitives[0].exponent=std::numeric_limits<double>::infinity(); },"infinite exponent");
    check([](auto& wf){ wf.primitives[0].coefficient=std::numeric_limits<double>::quiet_NaN(); },"NaN primitive");
    check([](auto& wf){ wf.orbitals[0].coefficients[0]=std::numeric_limits<double>::quiet_NaN(); },"NaN MO");
    check([](auto& wf){ wf.atoms[0].x=std::numeric_limits<double>::max(); },"coordinate float overflow");
    auto terms=cov::pack_orbital_terms(valid,0);
    terms[0].exponent=0;
    float output=0;
    rejects<std::invalid_argument>([&]{ cov::evaluate_cpu_grid(terms,cov::make_grid_request({},1,1,1,0,1),std::span<float>(&output,1)); },"invalid packed term");
}

} // namespace

int main() {
    try {
        numeric_components();
        grid_chunks_and_bounds();
        malformed_wavefunction();
        std::cout << "orbital_grid_smoke: normalized s-p-d-f-g grids verified\n";
    } catch (const std::exception& error) {
        std::cerr << "orbital_grid_smoke: " << error.what() << '\n';
        return 1;
    }
}
