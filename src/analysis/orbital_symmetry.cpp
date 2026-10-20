#include "cov/orbital_symmetry.hpp"
#include "cov/d2h_orbital_characters.hpp"

#include "cov/symmetry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>;

bool close_value(const double a, const double b, const double rel = 2.0e-5) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

double dot3(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vec3 transpose_mat_vec(const Mat3& m, const Vec3& v) {
    return {
        m[0]*v[0] + m[3]*v[1] + m[6]*v[2],
        m[1]*v[0] + m[4]*v[1] + m[7]*v[2],
        m[2]*v[0] + m[5]*v[1] + m[8]*v[2],
    };
}

Mat3 mat_mul(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            for (int k = 0; k < 3; ++k) out[3*r+c] += a[3*r+k] * b[3*k+c];
        }
    }
    return out;
}

int double_factorial_odd(int n) {
    if (n <= 0) return 1;
    int value = 1;
    for (int k = n; k > 1; k -= 2) value *= k;
    return value;
}

void cartesian_exponents(const int l, const int index, int& ax, int& ay, int& az) {
    ax = ay = az = 0;
    if (l == 0) return;
    if (l == 1) {
        constexpr int table[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    if (l == 2) {
        constexpr int table[6][3] = {
            {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}
        };
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    if (l == 3) {
        constexpr int table[10][3] = {
            {3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
            {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}
        };
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    constexpr int table[15][3] = {
        {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
        {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
        {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}
    };
    ax=table[index][0]; ay=table[index][1]; az=table[index][2];
}

double powi(double x, int n) {
    double out = 1.0;
    for (int i = 0; i < n; ++i) out *= x;
    return out;
}

double cartesian_angular(const int l, const int index, const Vec3& p) {
    int ax=0, ay=0, az=0;
    cartesian_exponents(l, index, ax, ay, az);
    const double denom = std::sqrt(static_cast<double>(
        double_factorial_odd(2*ax-1) *
        double_factorial_odd(2*ay-1) *
        double_factorial_odd(2*az-1)));
    return powi(p[0], ax) * powi(p[1], ay) * powi(p[2], az) / denom;
}

double real_solid_harmonic(const int l, const int index,
                           const double x, const double y, const double z) {
    const double x2=x*x, y2=y*y, z2=z*z, r2=x2+y2+z2;
    if (l == 0) return 0.28209479177387814;
    if (l == 1) {
        switch (index) {
            case 0: return  0.4886025119029199*z;
            case 1: return -0.4886025119029199*x;
            default:return -0.4886025119029199*y;
        }
    }
    if (l == 2) {
        switch (index) {
            case 0: return  0.31539156525252005*(3.0*z2-r2);
            case 1: return -1.0925484305920792*x*z;
            case 2: return -1.0925484305920792*y*z;
            case 3: return  0.5462742152960396*(x2-y2);
            default:return  1.0925484305920792*x*y;
        }
    }
    if (l == 3) {
        switch (index) {
            case 0: return  0.3731763325901154*z*(5.0*z2-3.0*r2);
            case 1: return -0.4570457994644658*x*(5.0*z2-r2);
            case 2: return -0.4570457994644658*y*(5.0*z2-r2);
            case 3: return  1.445305721320277*z*(x2-y2);
            case 4: return  2.890611442640554*x*y*z;
            case 5: return -0.5900435899266435*x*(x2-3.0*y2);
            default:return -0.5900435899266435*y*(3.0*x2-y2);
        }
    }
    const double z4=z2*z2, r4=r2*r2;
    switch (index) {
        case 0: return 0.10578554691520431*(35.0*z4-30.0*z2*r2+3.0*r4);
        case 1: return -0.6690465435572892*x*z*(7.0*z2-3.0*r2);
        case 2: return -0.6690465435572892*y*z*(7.0*z2-3.0*r2);
        case 3: return 0.47308734787878004*(x2-y2)*(7.0*z2-r2);
        case 4: return 0.9461746957575601*x*y*(7.0*z2-r2);
        case 5: return -1.7701307697799304*x*z*(x2-3.0*y2);
        case 6: return -1.7701307697799304*y*z*(3.0*x2-y2);
        case 7: return 0.6258357354491761*(x2*x2-6.0*x2*y2+y2*y2);
        default:return 2.5033429417967046*x*y*(x2-y2);
    }
}

double angular_value(const Shell& shell, const int component, const Vec3& p) {
    if (shell.pure) {
        return real_solid_harmonic(static_cast<int>(shell.angular_momentum), component,
                                   p[0], p[1], p[2]);
    }
    return cartesian_angular(static_cast<int>(shell.angular_momentum), component, p);
}

bool invert_square(std::vector<double> a, std::vector<double>& inverse, const std::size_t n) {
    inverse.assign(n*n, 0.0);
    for (std::size_t i=0;i<n;++i) inverse[i*n+i]=1.0;
    double scale=0.0;
    for (double v:a) scale=std::max(scale,std::abs(v));
    if (!(scale>0.0)) return false;
    for (std::size_t col=0;col<n;++col) {
        std::size_t pivot=col;
        double best=std::abs(a[col*n+col]);
        for (std::size_t row=col+1;row<n;++row) {
            const double v=std::abs(a[row*n+col]);
            if (v>best) { best=v; pivot=row; }
        }
        if (best <= 1.0e-12*scale || !std::isfinite(best)) return false;
        if (pivot!=col) {
            for (std::size_t j=0;j<n;++j) {
                std::swap(a[col*n+j],a[pivot*n+j]);
                std::swap(inverse[col*n+j],inverse[pivot*n+j]);
            }
        }
        const double d=a[col*n+col];
        for (std::size_t j=0;j<n;++j) { a[col*n+j]/=d; inverse[col*n+j]/=d; }
        for (std::size_t row=0;row<n;++row) {
            if (row==col) continue;
            const double f=a[row*n+col];
            if (f==0.0) continue;
            for (std::size_t j=0;j<n;++j) {
                a[row*n+j]-=f*a[col*n+j];
                inverse[row*n+j]-=f*inverse[col*n+j];
            }
        }
    }
    return true;
}

std::vector<Vec3> collocation_points(const std::size_t count) {
    std::vector<Vec3> points;
    points.reserve(count);
    constexpr double golden = 2.39996322972865332;
    for (std::size_t i=0;i<count;++i) {
        const double t=(static_cast<double>(i)+0.5)/static_cast<double>(count);
        const double z=1.0-2.0*t;
        const double rxy=std::sqrt(std::max(0.0,1.0-z*z));
        const double phi=golden*static_cast<double>(i);
        const double radius=0.72+0.07*static_cast<double>(i%5u);
        points.push_back({radius*rxy*std::cos(phi), radius*rxy*std::sin(phi), radius*z});
    }
    return points;
}

std::vector<double> local_angular_transform(const Shell& shell, const Mat3& matrix) {
    const std::size_t n=shell_basis_count(shell);
    if (n==0) return {};
    if (n==1) return {1.0};
    const std::size_t samples=std::max<std::size_t>(64,6*n);
    const auto points=collocation_points(samples);
    std::vector<double> gram(n*n,0.0), crossm(n*n,0.0);
    for (const auto& p:points) {
        const Vec3 q=transpose_mat_vec(matrix,p);
        std::vector<double> b(n),bt(n);
        for (std::size_t i=0;i<n;++i) {
            b[i]=angular_value(shell,static_cast<int>(i),p);
            bt[i]=angular_value(shell,static_cast<int>(i),q);
        }
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t j=0;j<n;++j) {
                gram[i*n+j]+=b[i]*b[j];
                crossm[i*n+j]+=b[i]*bt[j];
            }
        }
    }
    std::vector<double> inv;
    if (!invert_square(gram,inv,n)) return {};
    std::vector<double> transform(n*n,0.0);
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<n;++j) {
            for (std::size_t k=0;k<n;++k) transform[i*n+j]+=inv[i*n+k]*crossm[k*n+j];
        }
    }

    double max_ref=0.0,max_err=0.0;
    for (const auto& p:points) {
        const Vec3 q=transpose_mat_vec(matrix,p);
        for (std::size_t j=0;j<n;++j) {
            double reconstructed=0.0;
            for (std::size_t i=0;i<n;++i) {
                reconstructed+=angular_value(shell,static_cast<int>(i),p)*transform[i*n+j];
            }
            const double exact=angular_value(shell,static_cast<int>(j),q);
            max_ref=std::max(max_ref,std::abs(exact));
            max_err=std::max(max_err,std::abs(reconstructed-exact));
        }
    }
    if (!std::isfinite(max_err) || max_err > 2.0e-8*std::max(1.0,max_ref)) return {};
    return transform;
}

bool same_shell_shape(const Wavefunction& wf, const Shell& a, const Shell& b) {
    if (a.angular_momentum!=b.angular_momentum || a.pure!=b.pure ||
        a.primitive_count!=b.primitive_count) return false;
    for (std::uint32_t p=0;p<a.primitive_count;++p) {
        const auto& pa=wf.primitives[a.primitive_offset+p];
        const auto& pb=wf.primitives[b.primitive_offset+p];
        if (!close_value(pa.exponent,pb.exponent) || !close_value(pa.coefficient,pb.coefficient)) return false;
    }
    return true;
}

struct OperationMap {
    std::vector<std::size_t> target_shell;
    std::vector<std::vector<double>> local;
    bool valid=false;
};

OperationMap prepare_operation(const Wavefunction& wf, const SymmetryOperation& op) {
    OperationMap out;
    if (op.atom_permutation.size()!=wf.atoms.size()) return out;
    out.target_shell.assign(wf.shells.size(),wf.shells.size());
    out.local.resize(wf.shells.size());
    std::vector<bool> used(wf.shells.size(),false);

    for (std::size_t s=0;s<wf.shells.size();++s) {
        const Shell& source=wf.shells[s];
        if (source.atom_index>=op.atom_permutation.size()) return {};
        const std::size_t target_atom=op.atom_permutation[source.atom_index];
        std::size_t target=wf.shells.size();
        for (std::size_t t=0;t<wf.shells.size();++t) {
            if (used[t] || wf.shells[t].atom_index!=target_atom) continue;
            if (same_shell_shape(wf,source,wf.shells[t])) { target=t; break; }
        }
        if (target==wf.shells.size()) return {};
        const auto local=local_angular_transform(source,op.matrix);
        const std::size_t n=shell_basis_count(source);
        if (local.size()!=n*n || shell_basis_count(wf.shells[target])!=n) return {};
        used[target]=true;
        out.target_shell[s]=target;
        out.local[s]=local;
    }
    out.valid=true;
    return out;
}

std::vector<double> apply_operation(const Wavefunction& wf,
                                    const OperationMap& map,
                                    const std::vector<double>& coefficients) {
    if (!map.valid || coefficients.size()!=wf.basis_count) return {};
    std::vector<double> out(wf.basis_count,0.0);
    for (std::size_t s=0;s<wf.shells.size();++s) {
        const Shell& source=wf.shells[s];
        const Shell& target=wf.shells[map.target_shell[s]];
        const std::size_t n=shell_basis_count(source);
        const auto& t=map.local[s];
        for (std::size_t j=0;j<n;++j) {
            for (std::size_t i=0;i<n;++i) {
                out[target.basis_offset+j]+=t[j*n+i]*static_cast<double>(coefficients[source.basis_offset+i]);
            }
        }
    }
    return out;
}

double metric_inner(const Wavefunction& wf,
                    const std::vector<double>& a,
                    const std::vector<double>& b) {
    const std::size_t n=wf.basis_count;
    if (wf.ao_overlap.size()!=n*n || a.size()!=n || b.size()!=n) return 0.0;
    double value=0.0;
    for (std::size_t mu=0;mu<n;++mu) {
        double sb=0.0;
        for (std::size_t nu=0;nu<n;++nu) sb+=wf.ao_overlap[mu*n+nu]*b[nu];
        value+=static_cast<double>(a[mu])*sb;
    }
    return value;
}

double metric_norm(const Wavefunction& wf, const std::vector<double>& a) {
    const std::size_t n=wf.basis_count;
    if (wf.ao_overlap.size()!=n*n || a.size()!=n) return 0.0;
    double value=0.0;
    for (std::size_t mu=0;mu<n;++mu) {
        double sa=0.0;
        for (std::size_t nu=0;nu<n;++nu) sa+=wf.ao_overlap[mu*n+nu]*a[nu];
        value+=a[mu]*sa;
    }
    return value;
}

struct CharacterEval {
    double character=0.0;
    double retention=0.0;
    bool valid=false;
};

CharacterEval evaluate_character(const Wavefunction& wf,
                                 const SymmetryOperation& op,
                                 const std::vector<std::size_t>& orbitals) {
    CharacterEval result;
    const OperationMap map=prepare_operation(wf,op);
    if (!map.valid || orbitals.empty()) return result;
    double min_retention=1.0;
    double chi=0.0;
    for (std::size_t index:orbitals) {
        if (index>=wf.orbitals.size()) return {};
        const auto transformed=apply_operation(wf,map,wf.orbitals[index].coefficients);
        if (transformed.empty()) return {};
        chi+=metric_inner(wf,wf.orbitals[index].coefficients,transformed);
        double projected=0.0;
        for (std::size_t j:orbitals) {
            const double overlap=metric_inner(wf,wf.orbitals[j].coefficients,transformed);
            projected+=overlap*overlap;
        }
        const double n=metric_norm(wf,transformed);
        if (!(n>1.0e-12) || !std::isfinite(n)) return {};
        min_retention=std::min(min_retention,projected/n);
    }
    result.character=chi;
    result.retention=min_retention;
    result.valid=std::isfinite(chi)&&std::isfinite(min_retention);
    return result;
}

const SymmetryOperation* find_operation(const MolecularSymmetry& symmetry,
                                        const SymmetryOperationKind kind,
                                        const int order) {
    for (const auto& op:symmetry.operations) {
        if (op.kind==kind && op.order==order) return &op;
    }
    return nullptr;
}

const SymmetryOperation* principal_rotation(const MolecularSymmetry& symmetry, int& order) {
    order=1;
    const SymmetryOperation* best=nullptr;
    for (const auto& op:symmetry.operations) {
        if (op.kind==SymmetryOperationKind::ProperRotation && op.order>order) {
            order=op.order; best=&op;
        }
    }
    return best;
}

const SymmetryOperation* perpendicular_c2(const Wavefunction& wf,
                                         const MolecularSymmetry& symmetry, const Vec3& axis) {
    const SymmetryOperation* best=nullptr;
    std::size_t most_on_axis=0;
    for (const auto& op:symmetry.operations) {
        if (op.kind==SymmetryOperationKind::ProperRotation && op.order==2 &&
            std::abs(dot3(op.axis_or_normal,axis))<2.0e-5) {
            std::size_t on_axis=0;
            for (const auto& atom:wf.atoms) {
                const Vec3 r{atom.x-symmetry.centre_bohr[0],
                             atom.y-symmetry.centre_bohr[1],
                             atom.z-symmetry.centre_bohr[2]};
                const double along=dot3(r,op.axis_or_normal);
                const double perpendicular2=std::max(0.0,dot3(r,r)-along*along);
                if (perpendicular2<=symmetry.tolerance_bohr*symmetry.tolerance_bohr) ++on_axis;
            }
            if (!best || on_axis>most_on_axis) { best=&op;most_on_axis=on_axis; }
        }
    }
    return best;
}

const SymmetryOperation* horizontal_reflection(const MolecularSymmetry& symmetry, const Vec3& axis) {
    for (const auto& op:symmetry.operations) {
        if (op.kind==SymmetryOperationKind::Reflection &&
            std::abs(dot3(op.axis_or_normal,axis))>1.0-2.0e-5) return &op;
    }
    return nullptr;
}

SymmetryOperation squared_operation(const SymmetryOperation& op) {
    SymmetryOperation out=op;
    out.kind=SymmetryOperationKind::ProperRotation;
    out.order=2;
    out.power=2;
    out.matrix=mat_mul(op.matrix,op.matrix);
    out.atom_permutation.assign(op.atom_permutation.size(),0u);
    for (std::size_t i=0;i<op.atom_permutation.size();++i) {
        const std::size_t mid=op.atom_permutation[i];
        if (mid>=op.atom_permutation.size()) { out.atom_permutation.clear(); break; }
        out.atom_permutation[i]=op.atom_permutation[mid];
    }
    return out;
}

std::optional<std::string> classify_d2h(const Wavefunction& wf,const MolecularSymmetry& symmetry,const std::vector<std::size_t>& group,const OrbitalSymmetryOptions& options,double& retention,DerivedOrbitalSymmetryAssignment& evidence){
    if(symmetry.point_group!="D2h"||group.size()!=1)return std::nullopt;
    const auto ops=orbital_characters::d2h_operations(wf,symmetry.operations);
    const auto axes=orbital_characters::d2h_axes(wf,ops,symmetry.centre_bohr,symmetry.tolerance_bohr);
    std::vector<double> characters;retention=1.;
    for(const auto& op:ops){const auto value=evaluate_character(wf,op,group);if(!value.valid)return std::nullopt;retention=std::min(retention,value.retention);characters.push_back(value.character);}
    if(retention<options.minimum_subspace_retention)return std::nullopt;
    const auto label=orbital_characters::match_d2h(axes,ops,characters,group.size(),options.character_tolerance);
    if(label.empty())return std::nullopt;
    evidence.centre_bohr=symmetry.centre_bohr;evidence.axes_available=axes.resolved;evidence.axis_convention=axes.detail;
    if(axes.resolved){evidence.principal_axis=axes.xyz[2];evidence.secondary_axis=axes.xyz[0];}
    return label;
}

std::optional<std::string> classify_dnh(const Wavefunction& wf,
                                        const MolecularSymmetry& symmetry,
                                        const std::vector<std::size_t>& group,
                                        const OrbitalSymmetryOptions& options,
                                        double& retention,
                                        DerivedOrbitalSymmetryAssignment& evidence) {
    int n=1;
    const auto* cn=principal_rotation(symmetry,n);
    if (!cn || n<3 || symmetry.point_group != "D"+std::to_string(n)+"h") return std::nullopt;
    const auto* c2=perpendicular_c2(wf,symmetry,cn->axis_or_normal);
    const auto* sh=horizontal_reflection(symmetry,cn->axis_or_normal);
    if (!c2 || !sh) return std::nullopt;
    const auto a=evaluate_character(wf,*cn,group);
    const auto b=evaluate_character(wf,*c2,group);
    const auto h=evaluate_character(wf,*sh,group);
    if (!a.valid||!b.valid||!h.valid) return std::nullopt;
    retention=std::min({a.retention,b.retention,h.retention});
    if (retention<options.minimum_subspace_retention) return std::nullopt;
    const int d=static_cast<int>(group.size());
    if (std::abs(std::abs(h.character)-static_cast<double>(d))>options.character_tolerance) return std::nullopt;
    evidence.centre_bohr=symmetry.centre_bohr;
    evidence.principal_axis=cn->axis_or_normal;
    evidence.secondary_axis=c2->axis_or_normal;
    evidence.axes_available=true;
    evidence.axis_convention="C2 prime is the recorded perpendicular axis with maximal on-axis atom support; B1/B2 are relative to it";
    evidence.maximum_character_error=std::abs(std::abs(h.character)-static_cast<double>(d));
    std::string suffix=h.character>=0.0 ? "'" : "''";
    double inversion_character=0.0;
    if ((n%2)==0) {
        const auto* inv=find_operation(symmetry,SymmetryOperationKind::Inversion,2);
        if (!inv) return std::nullopt;
        const auto parity=evaluate_character(wf,*inv,group);
        if (!parity.valid) return std::nullopt;
        retention=std::min(retention,parity.retention);
        const double error=std::abs(std::abs(parity.character)-static_cast<double>(d));
        if (retention<options.minimum_subspace_retention || error>options.character_tolerance) return std::nullopt;
        evidence.maximum_character_error=std::max(evidence.maximum_character_error,error);
        inversion_character=parity.character;
        suffix=parity.character>=0.0 ? "g" : "u";
    }
    if (d==1) {
        if (std::abs(std::abs(a.character)-1.0)>options.character_tolerance ||
            std::abs(std::abs(b.character)-1.0)>options.character_tolerance) return std::nullopt;
        if (a.character<0.0 && (n%2)!=0) return std::nullopt;
        if ((n%2)==0) {
            const bool half_turn_odd=a.character<0.0 && ((n/2)%2)!=0;
            const bool expected_even=(h.character>=0.0)!=half_turn_odd;
            if ((inversion_character>=0.0)!=expected_even) return std::nullopt;
        }
        evidence.maximum_character_error=std::max({evidence.maximum_character_error,
            std::abs(std::abs(a.character)-1.0),std::abs(std::abs(b.character)-1.0)});
        const std::string family=a.character>=0.0?"A":"B";
        return family+(b.character>=0.0?"1":"2")+suffix;
    }
    if (d==2) {
        int best_k=0;
        double best_error=std::numeric_limits<double>::infinity();
        const int kmax=(n-1)/2;
        for (int k=1;k<=kmax;++k) {
            const double expected=2.0*std::cos(2.0*kPi*static_cast<double>(k)/static_cast<double>(n));
            const double error=std::abs(a.character-expected);
            if (error<best_error) { best_error=error; best_k=k; }
        }
        if (best_k==0 || best_error>options.character_tolerance) return std::nullopt;
        if (std::abs(b.character)>options.character_tolerance) return std::nullopt;
        if ((n%2)==0 && (inversion_character>=0.0)!=((h.character>=0.0)!=((best_k%2)!=0))) return std::nullopt;
        evidence.maximum_character_error=std::max({evidence.maximum_character_error,best_error,std::abs(b.character)});
        return "E"+(kmax==1?std::string{}:std::to_string(best_k))+suffix;
    }
    return std::nullopt;
}

struct TableRow {
    const char* label;
    int dimension;
    std::array<double,5> chars;
    int count;
};

std::optional<std::string> nearest_table(const std::vector<double>& observed,
                                         const int dimension,
                                         const std::vector<TableRow>& table,
                                         const double tolerance) {
    const TableRow* best=nullptr;
    double best_error=std::numeric_limits<double>::infinity();
    for (const auto& row:table) {
        if (row.dimension!=dimension || row.count!=static_cast<int>(observed.size())) continue;
        double max_error=0.0;
        for (int i=0;i<row.count;++i) max_error=std::max(max_error,std::abs(observed[static_cast<std::size_t>(i)]-row.chars[static_cast<std::size_t>(i)]));
        if (max_error<best_error) { best_error=max_error; best=&row; }
    }
    if (!best || best_error>tolerance) return std::nullopt;
    return std::string(best->label);
}

std::optional<std::string> classify_td(const Wavefunction& wf,
                                       const MolecularSymmetry& symmetry,
                                       const std::vector<std::size_t>& group,
                                       const OrbitalSymmetryOptions& options,
                                       double& retention) {
    if (symmetry.point_group!="Td") return std::nullopt;
    const auto* c3=find_operation(symmetry,SymmetryOperationKind::ProperRotation,3);
    const auto* c2=find_operation(symmetry,SymmetryOperationKind::ProperRotation,2);
    const auto* s4=find_operation(symmetry,SymmetryOperationKind::ImproperRotation,4);
    const auto* sd=find_operation(symmetry,SymmetryOperationKind::Reflection,2);
    if (!c3||!c2||!s4||!sd) return std::nullopt;
    const auto a=evaluate_character(wf,*c3,group);
    const auto b=evaluate_character(wf,*c2,group);
    const auto c=evaluate_character(wf,*s4,group);
    const auto d=evaluate_character(wf,*sd,group);
    if (!a.valid||!b.valid||!c.valid||!d.valid) return std::nullopt;
    retention=std::min({a.retention,b.retention,c.retention,d.retention});
    if (retention<options.minimum_subspace_retention) return std::nullopt;
    const std::vector<TableRow> table={
        {"A1",1,{ 1, 1, 1, 1,0},4},
        {"A2",1,{ 1, 1,-1,-1,0},4},
        {"E", 2,{-1, 2, 0, 0,0},4},
        {"T1",3,{ 0,-1, 1,-1,0},4},
        {"T2",3,{ 0,-1,-1, 1,0},4},
    };
    return nearest_table({a.character,b.character,c.character,d.character},
                         static_cast<int>(group.size()),table,options.character_tolerance);
}

std::optional<std::string> classify_oh(const Wavefunction& wf,
                                       const MolecularSymmetry& symmetry,
                                       const std::vector<std::size_t>& group,
                                       const OrbitalSymmetryOptions& options,
                                       double& retention) {
    if (symmetry.point_group!="Oh") return std::nullopt;
    const auto* c3=find_operation(symmetry,SymmetryOperationKind::ProperRotation,3);
    const auto* c4=find_operation(symmetry,SymmetryOperationKind::ProperRotation,4);
    const auto* inv=find_operation(symmetry,SymmetryOperationKind::Inversion,2);
    if (!c3||!c4||!inv) return std::nullopt;
    const SymmetryOperation c2axis=squared_operation(*c4);
    const SymmetryOperation* c2edge=nullptr;
    for (const auto& op:symmetry.operations) {
        if (op.kind!=SymmetryOperationKind::ProperRotation || op.order!=2) continue;
        if (std::abs(dot3(op.axis_or_normal,c4->axis_or_normal))<1.0-2.0e-5) { c2edge=&op; break; }
    }
    if (!c2edge) return std::nullopt;
    const auto a=evaluate_character(wf,*c3,group);
    const auto b=evaluate_character(wf,*c2edge,group);
    const auto c=evaluate_character(wf,*c4,group);
    const auto d=evaluate_character(wf,c2axis,group);
    const auto e=evaluate_character(wf,*inv,group);
    if (!a.valid||!b.valid||!c.valid||!d.valid||!e.valid) return std::nullopt;
    retention=std::min({a.retention,b.retention,c.retention,d.retention,e.retention});
    if (retention<options.minimum_subspace_retention) return std::nullopt;
    const int dim=static_cast<int>(group.size());
    if (std::abs(std::abs(e.character)-static_cast<double>(dim))>options.character_tolerance) return std::nullopt;
    const char parity=e.character>=0.0 ? 'g' : 'u';
    const std::vector<TableRow> table={
        {"A1",1,{ 1, 1, 1, 1,0},4},
        {"A2",1,{ 1,-1,-1, 1,0},4},
        {"E", 2,{-1, 0, 0, 2,0},4},
        {"T1",3,{ 0,-1, 1,-1,0},4},
        {"T2",3,{ 0, 1,-1,-1,0},4},
    };
    const auto base=nearest_table({a.character,b.character,c.character,d.character},
                                  dim,table,options.character_tolerance);
    if (!base) return std::nullopt;
    return *base+parity;
}

std::vector<std::vector<std::size_t>> energy_groups(const Wavefunction& wf, const double tolerance) {
    std::vector<std::vector<std::size_t>> groups;
    std::size_t i=0;
    while (i<wf.orbitals.size()) {
        std::vector<std::size_t> group{i};
        const Spin spin=wf.orbitals[i].spin;
        const double energy=wf.orbitals[i].energy_hartree;
        std::size_t j=i+1;
        while (j<wf.orbitals.size() && wf.orbitals[j].spin==spin &&
               std::abs(wf.orbitals[j].energy_hartree-energy)<=tolerance) {
            group.push_back(j++);
        }
        groups.push_back(std::move(group));
        i=j;
    }
    return groups;
}

} // namespace

std::vector<double> apply_orbital_symmetry_operation(const Wavefunction& wf,
    const SymmetryOperation& operation, const std::vector<double>& coefficients,
    const std::size_t columns) {
    if (!columns || coefficients.size()!=std::size_t(wf.basis_count)*columns) return {};
    for (const auto& shell:wf.shells) if (shell.angular_momentum>4) return {};
    const auto map=prepare_operation(wf,operation);
    if (!map.valid) return {};
    std::vector<double> out(coefficients.size(),0.0);
    for (std::size_t s=0;s<wf.shells.size();++s) {
        const auto& source=wf.shells[s];
        const auto& target=wf.shells[map.target_shell[s]];
        const auto n=shell_basis_count(source);
        for (std::size_t j=0;j<n;++j) for (std::size_t i=0;i<n;++i) {
            const double t=map.local[s][j*n+i];
            for (std::size_t k=0;k<columns;++k)
                out[(target.basis_offset+j)*columns+k]+=t*coefficients[(source.basis_offset+i)*columns+k];
        }
    }
    return out;
}

OrbitalSymmetryResult derive_orbital_symmetry(Wavefunction& wavefunction,
                                               const OrbitalSymmetryOptions& options) {
    OrbitalSymmetryResult result;
    if (wavefunction.orbitals.empty() || wavefunction.basis_count==0 ||
        wavefunction.ao_overlap.size()!=static_cast<std::size_t>(wavefunction.basis_count)*wavefunction.basis_count) {
        return result;
    }
    const MolecularSymmetry symmetry=analyse_molecular_symmetry(wavefunction);
    result.point_group=symmetry.point_group;
    if (!symmetry.available()) return result;

    for (const auto& group:energy_groups(wavefunction,options.degeneracy_tolerance_hartree)) {
        bool producer=false;
        bool already_labelled=false;
        for (std::size_t index:group) {
            producer=producer || wavefunction.orbitals[index].symmetry_provenance==DataProvenance::Producer;
            already_labelled=already_labelled || !wavefunction.orbitals[index].symmetry.empty();
        }
        if (producer) continue;
        if (already_labelled) continue;
        ++result.groups_examined;
        double retention=1.0;
        DerivedOrbitalSymmetryAssignment evidence;
        std::optional<std::string> label;
        if (symmetry.point_group=="D2h") {
            label=classify_d2h(wavefunction,symmetry,group,options,retention,evidence);
        } else if (symmetry.point_group=="Td") {
            label=classify_td(wavefunction,symmetry,group,options,retention);
        } else if (symmetry.point_group=="Oh") {
            label=classify_oh(wavefunction,symmetry,group,options,retention);
        } else if (symmetry.point_group.size()>=3 && symmetry.point_group[0]=='D' &&
                   symmetry.point_group.back()=='h') {
            label=classify_dnh(wavefunction,symmetry,group,options,retention,evidence);
        }
        result.worst_subspace_retention=std::min(result.worst_subspace_retention,retention);
        if (!label) continue;
        evidence.point_group=symmetry.point_group;
        evidence.label=*label;
        evidence.orbital_indices=group;
        evidence.subspace_retention=retention;
        evidence.centre_bohr=symmetry.centre_bohr;
        wavefunction.derived_orbital_symmetry_assignments.push_back(std::move(evidence));
        ++result.groups_labelled;
        for (std::size_t index:group) {
            wavefunction.orbitals[index].symmetry=*label;
            wavefunction.orbitals[index].symmetry_provenance=DataProvenance::Derived;
            ++result.orbitals_labelled;
        }
    }
    return result;
}

} // namespace cov
