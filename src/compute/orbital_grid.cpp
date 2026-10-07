#include "cov/orbital_grid.hpp"

#include "cov/analysis_threads.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace cov {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::size_t maximum_grid_points = 512u * 512u * 512u;

float packed(double value, const char* name) {
    const float result = static_cast<float>(value);
    if (!std::isfinite(value) || !std::isfinite(result))
        throw std::invalid_argument(std::string("Invalid rendering value: ") + name);
    return result;
}

// A small degree-four polynomial is only constructed when an orbital changes.
// The grid loop receives its nonzero monomials in CovGaussianTerm instead.
struct Polynomial {
    std::array<double, 125> c{};
    explicit Polynomial(double constant = 0.0) { c[0] = constant; }
    static Polynomial variable(int axis) {
        Polynomial p;
        p.c[axis == 0 ? 25 : axis == 1 ? 5 : 1] = 1.0;
        return p;
    }
};

Polynomial operator+(Polynomial a, const Polynomial& b) {
    for (std::size_t i = 0; i < a.c.size(); ++i) a.c[i] += b.c[i];
    return a;
}
Polynomial operator-(Polynomial a, const Polynomial& b) {
    for (std::size_t i = 0; i < a.c.size(); ++i) a.c[i] -= b.c[i];
    return a;
}
Polynomial operator*(double scale, Polynomial p) {
    for (double& coefficient : p.c) coefficient *= scale;
    return p;
}
Polynomial operator*(const Polynomial& a, const Polynomial& b) {
    Polynomial result;
    for (int ax = 0; ax <= 4; ++ax) for (int ay = 0; ay <= 4; ++ay)
        for (int az = 0; az <= 4; ++az) {
            const double lhs = a.c[25 * ax + 5 * ay + az];
            if (lhs == 0.0) continue;
            for (int bx = 0; bx + ax <= 4; ++bx)
                for (int by = 0; by + ay <= 4; ++by)
                    for (int bz = 0; bz + az <= 4; ++bz)
                        result.c[25 * (ax + bx) + 5 * (ay + by) + az + bz] +=
                            lhs * b.c[25 * bx + 5 * by + bz];
        }
    return result;
}

Polynomial harmonic(int l, int index) {
    const auto x = Polynomial::variable(0), y = Polynomial::variable(1);
    const auto z = Polynomial::variable(2);
    const auto x2 = x*x, y2 = y*y, z2 = z*z;
    const auto r2 = x2+y2+z2;
    if (l == 0) return Polynomial(0.28209479177387814);
    if (l == 1) {
        switch (index) {
            case 0: return 0.4886025119029199*z;
            case 1: return -0.4886025119029199*x;
            default: return -0.4886025119029199*y;
        }
    }
    if (l == 2) {
        switch (index) {
            case 0: return 0.31539156525252005*(3.0*z2-r2);
            case 1: return -1.0925484305920792*(x*z);
            case 2: return -1.0925484305920792*(y*z);
            case 3: return 0.5462742152960396*(x2-y2);
            default: return 1.0925484305920792*(x*y);
        }
    }
    if (l == 3) {
        switch (index) {
            case 0: return 0.3731763325901154*(z*(5.0*z2-3.0*r2));
            case 1: return -0.4570457994644658*(x*(5.0*z2-r2));
            case 2: return -0.4570457994644658*(y*(5.0*z2-r2));
            case 3: return 1.445305721320277*(z*(x2-y2));
            case 4: return 2.890611442640554*(x*y*z);
            case 5: return -0.5900435899266435*(x*(x2-3.0*y2));
            default: return -0.5900435899266435*(y*(3.0*x2-y2));
        }
    }
    const auto z4 = z2*z2, r4 = r2*r2;
    switch (index) {
        case 0: return 0.10578554691520431*(35.0*z4-30.0*(z2*r2)+3.0*r4);
        case 1: return -0.6690465435572892*(x*z*(7.0*z2-3.0*r2));
        case 2: return -0.6690465435572892*(y*z*(7.0*z2-3.0*r2));
        case 3: return 0.47308734787878004*((x2-y2)*(7.0*z2-r2));
        case 4: return 0.9461746957575601*(x*y*(7.0*z2-r2));
        case 5: return -1.7701307697799304*(x*z*(x2-3.0*y2));
        case 6: return -1.7701307697799304*(y*z*(3.0*x2-y2));
        case 7: return 0.6258357354491761*(x2*x2-6.0*(x2*y2)+y2*y2);
        default: return 2.5033429417967046*(x*y*(x2-y2));
    }
}

constexpr std::array<std::array<std::uint32_t, 3>, 35> cartesian_order{{
    {0,0,0},
    {1,0,0},{0,1,0},{0,0,1},
    {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1},
    {3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
    {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1},
    {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
    {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
    {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}
}};

int odd_double_factorial(int n) {
    int result = 1;
    for (; n > 1; n -= 2) result *= n;
    return result;
}

double cartesian_norm(double alpha, std::uint32_t ax,
                      std::uint32_t ay, std::uint32_t az) {
    const auto l = ax+ay+az;
    const double denominator = std::sqrt(static_cast<double>(
        odd_double_factorial(2*static_cast<int>(ax)-1) *
        odd_double_factorial(2*static_cast<int>(ay)-1) *
        odd_double_factorial(2*static_cast<int>(az)-1)));
    return std::pow(2.0*alpha/pi, 0.75) *
           std::pow(4.0*alpha, 0.5*static_cast<double>(l)) / denominator;
}

double spherical_norm(double alpha, int l) {
    return std::sqrt(2.0*std::pow(2.0*alpha, l+1.5) /
                     std::tgamma(l+1.5));
}

void validate_request(const CovGridRequest& request) {
    if (request.nx > static_cast<std::uint32_t>(std::numeric_limits<int>::max()) ||
        request.ny > static_cast<std::uint32_t>(std::numeric_limits<int>::max()) ||
        request.nz > static_cast<std::uint32_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Grid dimensions are out of range");
    const auto size = checked_grid_size(static_cast<int>(request.nx),
                                        static_cast<int>(request.ny),
                                        static_cast<int>(request.nz));
    const float bounds[] = {request.min_x,request.min_y,request.min_z,
                            request.max_x,request.max_y,request.max_z};
    for (float bound : bounds)
        if (!std::isfinite(bound)) throw std::invalid_argument("Invalid grid bounds");
    if (request.min_x > request.max_x || request.min_y > request.max_y ||
        request.min_z > request.max_z ||
        !std::isfinite(request.max_x-request.min_x) ||
        !std::isfinite(request.max_y-request.min_y) ||
        !std::isfinite(request.max_z-request.min_z))
        throw std::invalid_argument("Invalid grid bounds");
    if (request.first_point > size || request.point_count > size-request.first_point)
        throw std::out_of_range("Grid chunk is out of range");
}

} // namespace

std::vector<CovGaussianTerm> pack_orbital_terms(const Wavefunction& wf,
                                               std::size_t mo_index) {
    if (mo_index >= wf.orbitals.size()) throw std::out_of_range("MO index out of range");
    const auto& mo = wf.orbitals[mo_index];
    if (mo.coefficients.size() != wf.basis_count)
        throw std::invalid_argument("MO coefficient dimension does not match basis");
    for (double coefficient : mo.coefficients)
        if (!std::isfinite(coefficient)) throw std::invalid_argument("Invalid MO coefficient");

    std::vector<CovGaussianTerm> terms;
    for (const Shell& shell : wf.shells) {
        if (shell.angular_momentum > 4 || shell.pure > 1 || shell.primitive_count == 0)
            throw std::invalid_argument("Invalid s through g shell");
        if (shell.atom_index >= wf.atoms.size() ||
            shell.primitive_offset > wf.primitives.size() ||
            shell.primitive_count > wf.primitives.size()-shell.primitive_offset ||
            shell.basis_offset > wf.basis_count ||
            shell_basis_count(shell) > wf.basis_count-shell.basis_offset)
            throw std::out_of_range("Shell index is out of range");
        const Atom& atom = wf.atoms[shell.atom_index];
        const float cx = packed(atom.x,"atom x"), cy = packed(atom.y,"atom y");
        const float cz = packed(atom.z,"atom z");
        const auto l = static_cast<int>(shell.angular_momentum);
        const auto components = shell_basis_count(shell);
        for (std::uint32_t c = 0; c < components; ++c) {
            const double mo_coefficient = mo.coefficients[shell.basis_offset+c];
            const auto powers = cartesian_order[static_cast<std::size_t>(l*(l+1)*(l+2)/6)+c];
            const Polynomial polynomial = shell.pure ? harmonic(l,static_cast<int>(c)) : Polynomial{};
            for (std::uint32_t p = 0; p < shell.primitive_count; ++p) {
                const Primitive& primitive = wf.primitives[shell.primitive_offset+p];
                if (!std::isfinite(primitive.exponent) || primitive.exponent <= 0.0 ||
                    !std::isfinite(primitive.coefficient))
                    throw std::invalid_argument("Invalid Gaussian primitive");
                const float exponent = packed(primitive.exponent,"primitive exponent");
                if (exponent <= 0.0f) throw std::invalid_argument("Primitive exponent underflows float");
                const double scale = mo_coefficient * primitive.coefficient *
                    (shell.pure ? spherical_norm(primitive.exponent,l) :
                                  cartesian_norm(primitive.exponent,powers[0],powers[1],powers[2]));
                if (!std::isfinite(scale)) throw std::invalid_argument("Invalid normalized coefficient");
                auto append = [&](std::uint32_t ax, std::uint32_t ay,
                                  std::uint32_t az, double factor) {
                    if (terms.size() == std::numeric_limits<std::uint32_t>::max())
                        throw std::length_error("Too many Gaussian terms");
                    terms.push_back({cx,cy,cz,exponent,
                                     packed(scale*factor,"normalized coefficient"),ax,ay,az});
                };
                if (!shell.pure) append(powers[0],powers[1],powers[2],1.0);
                else for (std::uint32_t ax = 0; ax <= 4; ++ax)
                    for (std::uint32_t ay = 0; ay <= 4; ++ay)
                        for (std::uint32_t az = 0; az <= 4; ++az) {
                            const double factor = polynomial.c[25*ax+5*ay+az];
                            if (factor != 0.0) append(ax,ay,az,factor);
                        }
            }
        }
    }
    return terms;
}

std::size_t checked_grid_size(int nx, int ny, int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    const auto x = static_cast<std::size_t>(nx), y = static_cast<std::size_t>(ny);
    const auto z = static_cast<std::size_t>(nz);
    if (x > maximum_grid_points/y || x*y > maximum_grid_points/z)
        throw std::length_error("Grid exceeds the supported size");
    return x*y*z;
}

CovGridRequest make_grid_request(const GridBox& box, int nx, int ny, int nz,
                                 std::size_t first, std::size_t count) {
    const auto size = checked_grid_size(nx,ny,nz);
    if (first > size || count > size-first ||
        count > std::numeric_limits<std::uint32_t>::max())
        throw std::out_of_range("Grid chunk is out of range");
    CovGridRequest request{box.min_x,box.min_y,box.min_z,
                           box.max_x,box.max_y,box.max_z,
                           static_cast<std::uint32_t>(nx),
                           static_cast<std::uint32_t>(ny),
                           static_cast<std::uint32_t>(nz),
                           static_cast<std::uint64_t>(first),
                           static_cast<std::uint32_t>(count)};
    validate_request(request);
    return request;
}

void evaluate_cpu_grid(std::span<const CovGaussianTerm> terms,
                       const CovGridRequest& request, std::span<float> output,
                       cov::stop_token stop) {
    validate_request(request);
    if (output.size() < request.point_count)
        throw std::invalid_argument("Grid output is too small");
    for (const auto& term : terms) {
        if (!std::isfinite(term.cx) || !std::isfinite(term.cy) ||
            !std::isfinite(term.cz) || !std::isfinite(term.exponent) ||
            term.exponent <= 0.0f || !std::isfinite(term.coefficient) ||
            term.ax > 4u || term.ay > 4u || term.az > 4u ||
            term.ax+term.ay+term.az > 4u)
            throw std::invalid_argument("Invalid Gaussian term");
    }
    if (stop.stop_requested()) throw std::runtime_error("Grid evaluation cancelled");
    const auto count = static_cast<std::size_t>(request.point_count);
    if (count == 0) return;
    const auto workers = std::min<std::size_t>(analysis_thread_budget(),count);
    const auto block = (count+workers-1)/workers;
    std::atomic<bool> invalid_result{false};
    auto work = [&](std::size_t begin, std::size_t end) {
        for (std::size_t local = begin; local < end; ++local) {
            if ((local & 63u) == 0 && (stop.stop_requested() || invalid_result.load())) return;
            const std::uint64_t point = request.first_point+local;
            const auto ix = static_cast<std::uint32_t>(point%request.nx);
            const auto iy = static_cast<std::uint32_t>((point/request.nx)%request.ny);
            const auto iz = static_cast<std::uint32_t>(point/(static_cast<std::uint64_t>(request.nx)*request.ny));
            const float tx = request.nx > 1 ? static_cast<float>(ix)/static_cast<float>(request.nx-1) : 0.0f;
            const float ty = request.ny > 1 ? static_cast<float>(iy)/static_cast<float>(request.ny-1) : 0.0f;
            const float tz = request.nz > 1 ? static_cast<float>(iz)/static_cast<float>(request.nz-1) : 0.0f;
            const float x = request.min_x+tx*(request.max_x-request.min_x);
            const float y = request.min_y+ty*(request.max_y-request.min_y);
            const float z = request.min_z+tz*(request.max_z-request.min_z);
            double value = 0.0;
            for (const auto& term : terms) {
                const double dx = static_cast<double>(x)-term.cx;
                const double dy = static_cast<double>(y)-term.cy;
                const double dz = static_cast<double>(z)-term.cz;
                const double radial = std::exp(-static_cast<double>(term.exponent)*(dx*dx+dy*dy+dz*dz));
                if (radial == 0.0) continue;
                auto power = [](double base, std::uint32_t n) {
                    double result = 1.0;
                    for (std::uint32_t i=0; i<n; ++i) result *= base;
                    return result;
                };
                value += static_cast<double>(term.coefficient)*radial*
                         power(dx,term.ax)*power(dy,term.ay)*power(dz,term.az);
            }
            const float result = static_cast<float>(value);
            if (!std::isfinite(value) || !std::isfinite(result)) {
                invalid_result.store(true);
                return;
            }
            output[local] = result;
        }
    };
    {
        std::vector<cov::jthread> threads;
        threads.reserve(workers > 0 ? workers-1 : 0);
        for (std::size_t worker = 1; worker < workers; ++worker) {
            const auto begin = worker*block;
            if (begin < count) threads.emplace_back(work,begin,std::min(count,begin+block));
        }
        work(0,std::min(count,block));
    }
    if (stop.stop_requested()) throw std::runtime_error("Grid evaluation cancelled");
    if (invalid_result.load()) throw std::runtime_error("Grid result is not finite float");
}

} // namespace cov
