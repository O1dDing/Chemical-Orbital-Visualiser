#include "cov/orbital_symmetry.hpp"

#include "cov/symmetry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include <utility>

namespace cov {

// The original finite-point-group implementation is compiled under this
// internal symbol.  Keeping it isolated lets the new COV-native linear-group
// implementation mature without taking a dependency on an external symmetry
// package or destabilising the already validated Dnh/Td paths.
OrbitalSymmetryResult derive_orbital_symmetry_legacy(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options);

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>;

double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    };
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a,a));
}

bool normalize(Vec3& a) {
    const double n = norm(a);
    if (!(n > 1.0e-14) || !std::isfinite(n)) return false;
    for (double& v : a) v /= n;
    return true;
}

Mat3 rotation_matrix(Vec3 axis, const double angle) {
    normalize(axis);
    const double x=axis[0], y=axis[1], z=axis[2];
    const double c=std::cos(angle), s=std::sin(angle), t=1.0-c;
    return {
        t*x*x+c,   t*x*y-s*z, t*x*z+s*y,
        t*x*y+s*z, t*y*y+c,   t*y*z-s*x,
        t*x*z-s*y, t*y*z+s*x, t*z*z+c,
    };
}

Mat3 reflection_matrix(Vec3 normal) {
    normalize(normal);
    const double x=normal[0], y=normal[1], z=normal[2];
    return {
        1.0-2.0*x*x, -2.0*x*y,       -2.0*x*z,
        -2.0*y*x,      1.0-2.0*y*y,  -2.0*y*z,
        -2.0*z*x,     -2.0*z*y,        1.0-2.0*z*z,
    };
}

Vec3 transpose_mat_vec(const Mat3& m, const Vec3& v) {
    return {
        m[0]*v[0] + m[3]*v[1] + m[6]*v[2],
        m[1]*v[0] + m[4]*v[1] + m[7]*v[2],
        m[2]*v[0] + m[5]*v[1] + m[8]*v[2],
    };
}

int double_factorial_odd(int n) {
    if (n <= 0) return 1;
    int value=1;
    for (int k=n;k>1;k-=2) value*=k;
    return value;
}

void cartesian_exponents(const int l, const int index, int& ax, int& ay, int& az) {
    ax=ay=az=0;
    if (l==0) return;
    if (l==1) {
        constexpr int table[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    if (l==2) {
        constexpr int table[6][3]={{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    if (l==3) {
        constexpr int table[10][3]={{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
                                     {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}};
        ax=table[index][0]; ay=table[index][1]; az=table[index][2]; return;
    }
    constexpr int table[15][3]={{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
                                 {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
                                 {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}};
    ax=table[index][0]; ay=table[index][1]; az=table[index][2];
}

double powi(double x, int n) {
    double out=1.0;
    for (int i=0;i<n;++i) out*=x;
    return out;
}

double cartesian_angular(const int l, const int index, const Vec3& p) {
    int ax=0,ay=0,az=0;
    cartesian_exponents(l,index,ax,ay,az);
    const double denom=std::sqrt(static_cast<double>(
        double_factorial_odd(2*ax-1)*double_factorial_odd(2*ay-1)*double_factorial_odd(2*az-1)));
    return powi(p[0],ax)*powi(p[1],ay)*powi(p[2],az)/denom;
}

double real_solid_harmonic(const int l, const int index,
                           const double x, const double y, const double z) {
    const double x2=x*x,y2=y*y,z2=z*z,r2=x2+y2+z2;
    if (l==0) return 0.28209479177387814;
    if (l==1) {
        switch(index) {
            case 0:return 0.4886025119029199*z;
            case 1:return -0.4886025119029199*x;
            default:return -0.4886025119029199*y;
        }
    }
    if (l==2) {
        switch(index) {
            case 0:return 0.31539156525252005*(3.0*z2-r2);
            case 1:return -1.0925484305920792*x*z;
            case 2:return -1.0925484305920792*y*z;
            case 3:return 0.5462742152960396*(x2-y2);
            default:return 1.0925484305920792*x*y;
        }
    }
    if (l==3) {
        switch(index) {
            case 0:return 0.3731763325901154*z*(5.0*z2-3.0*r2);
            case 1:return -0.4570457994644658*x*(5.0*z2-r2);
            case 2:return -0.4570457994644658*y*(5.0*z2-r2);
            case 3:return 1.445305721320277*z*(x2-y2);
            case 4:return 2.890611442640554*x*y*z;
            case 5:return -0.5900435899266435*x*(x2-3.0*y2);
            default:return -0.5900435899266435*y*(3.0*x2-y2);
        }
    }
    const double z4=z2*z2,r4=r2*r2;
    switch(index) {
        case 0:return 0.10578554691520431*(35.0*z4-30.0*z2*r2+3.0*r4);
        case 1:return -0.6690465435572892*x*z*(7.0*z2-3.0*r2);
        case 2:return -0.6690465435572892*y*z*(7.0*z2-3.0*r2);
        case 3:return 0.47308734787878004*(x2-y2)*(7.0*z2-r2);
        case 4:return 0.9461746957575601*x*y*(7.0*z2-r2);
        case 5:return -1.7701307697799304*x*z*(x2-3.0*y2);
        case 6:return -1.7701307697799304*y*z*(3.0*x2-y2);
        case 7:return 0.6258357354491761*(x2*x2-6.0*x2*y2+y2*y2);
        default:return 2.5033429417967046*x*y*(x2-y2);
    }
}

double angular_value(const Shell& shell, const int component, const Vec3& p) {
    if (shell.pure) {
        return real_solid_harmonic(static_cast<int>(shell.angular_momentum),component,p[0],p[1],p[2]);
    }
    return cartesian_angular(static_cast<int>(shell.angular_momentum),component,p);
}

bool invert_square(std::vector<double> a, std::vector<double>& inv, const std::size_t n) {
    inv.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) inv[i*n+i]=1.0;
    double scale=0.0;
    for (double v:a) scale=std::max(scale,std::abs(v));
    if (!(scale>0.0)) return false;
    for (std::size_t c=0;c<n;++c) {
        std::size_t p=c;
        double best=std::abs(a[c*n+c]);
        for (std::size_t r=c+1;r<n;++r) {
            const double v=std::abs(a[r*n+c]);
            if (v>best) {best=v;p=r;}
        }
        if (!std::isfinite(best) || best<=1.0e-12*scale) return false;
        if (p!=c) {
            for (std::size_t j=0;j<n;++j) {
                std::swap(a[c*n+j],a[p*n+j]);
                std::swap(inv[c*n+j],inv[p*n+j]);
            }
        }
        const double d=a[c*n+c];
        for (std::size_t j=0;j<n;++j) {a[c*n+j]/=d;inv[c*n+j]/=d;}
        for (std::size_t r=0;r<n;++r) {
            if (r==c) continue;
            const double f=a[r*n+c];
            if (f==0.0) continue;
            for (std::size_t j=0;j<n;++j) {
                a[r*n+j]-=f*a[c*n+j];
                inv[r*n+j]-=f*inv[c*n+j];
            }
        }
    }
    return true;
}

std::vector<Vec3> collocation_points(const std::size_t count) {
    std::vector<Vec3> points;
    points.reserve(count);
    constexpr double golden=2.39996322972865332;
    for (std::size_t i=0;i<count;++i) {
        const double t=(static_cast<double>(i)+0.5)/static_cast<double>(count);
        const double z=1.0-2.0*t;
        const double rxy=std::sqrt(std::max(0.0,1.0-z*z));
        const double phi=golden*static_cast<double>(i);
        const double radius=0.72+0.07*static_cast<double>(i%5u);
        points.push_back({radius*rxy*std::cos(phi),radius*rxy*std::sin(phi),radius*z});
    }
    return points;
}

std::vector<double> local_transform(const Shell& shell, const Mat3& matrix) {
    const std::size_t n=shell_basis_count(shell);
    if (n==0) return {};
    if (n==1) return {1.0};
    const auto points=collocation_points(std::max<std::size_t>(64,6*n));
    std::vector<double> gram(n*n,0.0),crossm(n*n,0.0);
    for (const auto& p:points) {
        const Vec3 q=transpose_mat_vec(matrix,p);
        std::vector<double> b(n),bt(n);
        for (std::size_t i=0;i<n;++i) {
            b[i]=angular_value(shell,static_cast<int>(i),p);
            bt[i]=angular_value(shell,static_cast<int>(i),q);
        }
        for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j) {
            gram[i*n+j]+=b[i]*b[j];
            crossm[i*n+j]+=b[i]*bt[j];
        }
    }
    std::vector<double> inv;
    if (!invert_square(gram,inv,n)) return {};
    std::vector<double> t(n*n,0.0);
    for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j)
        for (std::size_t k=0;k<n;++k) t[i*n+j]+=inv[i*n+k]*crossm[k*n+j];
    return t;
}

bool same_shell_shape(const Wavefunction& wf, const Shell& a, const Shell& b) {
    if (a.angular_momentum!=b.angular_momentum || a.pure!=b.pure ||
        a.primitive_count!=b.primitive_count) return false;
    for (std::uint32_t p=0;p<a.primitive_count;++p) {
        const auto& pa=wf.primitives[a.primitive_offset+p];
        const auto& pb=wf.primitives[b.primitive_offset+p];
        const double escale=std::max({1.0,std::abs(static_cast<double>(pa.exponent)),std::abs(static_cast<double>(pb.exponent))});
        const double cscale=std::max({1.0,std::abs(static_cast<double>(pa.coefficient)),std::abs(static_cast<double>(pb.coefficient))});
        if (std::abs(static_cast<double>(pa.exponent)-pb.exponent)>2.0e-5*escale ||
            std::abs(static_cast<double>(pa.coefficient)-pb.coefficient)>2.0e-5*cscale) return false;
    }
    return true;
}

struct PreparedOperation {
    std::vector<std::size_t> target_shell;
    std::vector<std::vector<double>> local;
    bool valid=false;
};

PreparedOperation prepare(const Wavefunction& wf, const SymmetryOperation& op) {
    PreparedOperation out;
    // This angular implementation has validated Cartesian/pure transforms
    // through g only. Do not index its 15-entry g table for an h/i shell.
    for(const auto& shell:wf.shells)if(shell.angular_momentum>4)return out;
    if (op.atom_permutation.size()!=wf.atoms.size()) return out;
    out.target_shell.assign(wf.shells.size(),wf.shells.size());
    out.local.resize(wf.shells.size());
    std::vector<bool> used(wf.shells.size(),false);
    for (std::size_t s=0;s<wf.shells.size();++s) {
        const auto& source=wf.shells[s];
        if (source.atom_index>=op.atom_permutation.size()) return {};
        const std::size_t target_atom=op.atom_permutation[source.atom_index];
        std::size_t target=wf.shells.size();
        for (std::size_t t=0;t<wf.shells.size();++t) {
            if (used[t] || wf.shells[t].atom_index!=target_atom) continue;
            if (same_shell_shape(wf,source,wf.shells[t])) {target=t;break;}
        }
        if (target==wf.shells.size()) return {};
        const std::size_t count=shell_basis_count(source);
        const auto local=local_transform(source,op.matrix);
        if (local.size()!=count*count || shell_basis_count(wf.shells[target])!=count) return {};
        used[target]=true;
        out.target_shell[s]=target;
        out.local[s]=local;
    }
    out.valid=true;
    return out;
}

std::vector<double> apply(const Wavefunction& wf,
                          const PreparedOperation& map,
                          const std::vector<double>& coefficients) {
    if (!map.valid || coefficients.size()!=wf.basis_count) return {};
    std::vector<double> out(wf.basis_count,0.0);
    for (std::size_t s=0;s<wf.shells.size();++s) {
        const auto& source=wf.shells[s];
        const auto& target=wf.shells[map.target_shell[s]];
        const std::size_t count=shell_basis_count(source);
        const auto& t=map.local[s];
        for (std::size_t j=0;j<count;++j) for (std::size_t i=0;i<count;++i) {
            out[target.basis_offset+j]+=t[j*count+i]*static_cast<double>(coefficients[source.basis_offset+i]);
        }
    }
    return out;
}

double metric_inner(const Wavefunction& wf,
                    const std::vector<double>& a,
                    const std::vector<double>& b) {
    const std::size_t n=wf.basis_count;
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
    double value=0.0;
    for (std::size_t mu=0;mu<n;++mu) {
        double sa=0.0;
        for (std::size_t nu=0;nu<n;++nu) sa+=wf.ao_overlap[mu*n+nu]*a[nu];
        value+=a[mu]*sa;
    }
    return value;
}

struct Character {
    double value=0.0;
    double retention=0.0;
    bool valid=false;
};

Character character(const Wavefunction& wf,
                    const SymmetryOperation& op,
                    const std::vector<std::size_t>& group) {
    Character result;
    const auto map=prepare(wf,op);
    if (!map.valid || group.empty()) return result;
    result.retention=1.0;
    for (std::size_t i:group) {
        if (i>=wf.orbitals.size()) return {};
        const auto transformed=apply(wf,map,wf.orbitals[i].coefficients);
        if (transformed.empty()) return {};
        result.value+=metric_inner(wf,wf.orbitals[i].coefficients,transformed);
        double projected=0.0;
        for (std::size_t j:group) {
            const double x=metric_inner(wf,wf.orbitals[j].coefficients,transformed);
            projected+=x*x;
        }
        const double n=metric_norm(wf,transformed);
        if (!(n>1.0e-12) || !std::isfinite(n)) return {};
        result.retention=std::min(result.retention,projected/n);
    }
    result.valid=std::isfinite(result.value)&&std::isfinite(result.retention);
    return result;
}

std::vector<std::vector<std::size_t>> energy_groups(const Wavefunction& wf, const double tol) {
    std::vector<std::vector<std::size_t>> groups;
    std::size_t i=0;
    while (i<wf.orbitals.size()) {
        std::vector<std::size_t> group{i};
        const auto spin=wf.orbitals[i].spin;
        const double e=wf.orbitals[i].energy_hartree;
        std::size_t j=i+1;
        while (j<wf.orbitals.size() && wf.orbitals[j].spin==spin &&
               std::abs(wf.orbitals[j].energy_hartree-e)<=tol) group.push_back(j++);
        groups.push_back(std::move(group));
        i=j;
    }
    return groups;
}

bool group_unlabelled(const Wavefunction& wf, const std::vector<std::size_t>& group) {
    for (std::size_t i:group) {
        if (wf.orbitals[i].symmetry_provenance==DataProvenance::Producer ||
            !wf.orbitals[i].symmetry.empty()) return false;
    }
    return true;
}

const SymmetryOperation* find_operation(const MolecularSymmetry& sym,
                                        const SymmetryOperationKind kind,
                                        const int order) {
    for (const auto& op:sym.operations) if (op.kind==kind && op.order==order) return &op;
    return nullptr;
}

std::optional<std::string> classify_oh_fallback(const Wavefunction& wf,
                                                const MolecularSymmetry& sym,
                                                const std::vector<std::size_t>& group,
                                                const OrbitalSymmetryOptions& options,
                                                double& retention) {
    if (sym.point_group!="Oh") return std::nullopt;
    const auto* c3=find_operation(sym,SymmetryOperationKind::ProperRotation,3);
    const auto* c4=find_operation(sym,SymmetryOperationKind::ProperRotation,4);
    const auto* inv=find_operation(sym,SymmetryOperationKind::Inversion,2);
    if (!c3||!c4||!inv) return std::nullopt;
    const auto a=character(wf,*c3,group);
    const auto b=character(wf,*c4,group);
    const auto p=character(wf,*inv,group);
    if (!a.valid||!b.valid||!p.valid) return std::nullopt;
    retention=std::min({a.retention,b.retention,p.retention});
    if (retention<options.minimum_subspace_retention) return std::nullopt;
    const int d=static_cast<int>(group.size());
    if (std::abs(std::abs(p.value)-static_cast<double>(d))>options.character_tolerance) return std::nullopt;
    const char parity=p.value>=0.0?'g':'u';
    std::string base;
    if (d==1 && std::abs(a.value-1.0)<=options.character_tolerance) {
        if (std::abs(b.value-1.0)<=options.character_tolerance) base="A1";
        else if (std::abs(b.value+1.0)<=options.character_tolerance) base="A2";
    } else if (d==2 && std::abs(a.value+1.0)<=options.character_tolerance &&
               std::abs(b.value)<=options.character_tolerance) {
        base="E";
    } else if (d==3 && std::abs(a.value)<=options.character_tolerance) {
        if (std::abs(b.value-1.0)<=options.character_tolerance) base="T1";
        else if (std::abs(b.value+1.0)<=options.character_tolerance) base="T2";
    }
    if (base.empty()) return std::nullopt;
    return base+parity;
}

Vec3 linear_axis(const Wavefunction& wf, const MolecularSymmetry& sym) {
    Vec3 axis{0.0,0.0,1.0};
    double best=0.0;
    for (const auto& atom:wf.atoms) {
        Vec3 v{atom.x-sym.centre_bohr[0],atom.y-sym.centre_bohr[1],atom.z-sym.centre_bohr[2]};
        const double n2=dot(v,v);
        if (n2>best) {best=n2;axis=v;}
    }
    normalize(axis);
    return axis;
}

SymmetryOperation continuous_rotation(const Wavefunction& wf, const Vec3& axis, const double angle) {
    SymmetryOperation op;
    op.kind=SymmetryOperationKind::ProperRotation;
    op.order=0;
    op.power=1;
    op.axis_or_normal=axis;
    op.matrix=rotation_matrix(axis,angle);
    op.atom_permutation.resize(wf.atoms.size());
    for (std::size_t i=0;i<wf.atoms.size();++i) op.atom_permutation[i]=i;
    return op;
}

SymmetryOperation vertical_reflection(const Wavefunction& wf, const Vec3& axis) {
    Vec3 seed=std::abs(axis[0])<0.8?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    Vec3 normal=cross(axis,seed);
    if (!normalize(normal)) normal={0.0,1.0,0.0};
    SymmetryOperation op;
    op.kind=SymmetryOperationKind::Reflection;
    op.order=2;
    op.power=1;
    op.axis_or_normal=normal;
    op.matrix=reflection_matrix(normal);
    op.atom_permutation.resize(wf.atoms.size());
    for (std::size_t i=0;i<wf.atoms.size();++i) op.atom_permutation[i]=i;
    return op;
}

const char* linear_symbol(const int m) {
    switch(m) {
        case 0:return "Σ";
        case 1:return "Π";
        case 2:return "Δ";
        case 3:return "Φ";
        case 4:return "Γ";
        default:return nullptr;
    }
}

std::optional<std::string> classify_linear(const Wavefunction& wf,
                                           const MolecularSymmetry& sym,
                                           const std::vector<std::size_t>& group,
                                           const OrbitalSymmetryOptions& options,
                                           double& retention) {
    const bool dinfh=sym.point_group=="Dinfh";
    const bool cinfv=sym.point_group=="Cinfv";
    if (!dinfh&&!cinfv) return std::nullopt;
    const Vec3 axis=linear_axis(wf,sym);
    constexpr double theta=kPi/5.0; // separates |m|=0..4 for the supported s..g basis.
    const auto rot=continuous_rotation(wf,axis,theta);
    const auto r=character(wf,rot,group);
    if (!r.valid) return std::nullopt;
    retention=r.retention;
    const int d=static_cast<int>(group.size());
    int m=-1;
    if (d==1 && std::abs(r.value-1.0)<=options.character_tolerance) {
        m=0;
    } else if (d==2) {
        double best=std::numeric_limits<double>::infinity();
        for (int candidate=1;candidate<=4;++candidate) {
            const double expected=2.0*std::cos(theta*static_cast<double>(candidate));
            const double error=std::abs(r.value-expected);
            if (error<best) {best=error;m=candidate;}
        }
        if (best>options.character_tolerance) return std::nullopt;
    } else {
        return std::nullopt;
    }
    if (retention<options.minimum_subspace_retention) return std::nullopt;

    std::string label=linear_symbol(m);
    if (label.empty()) return std::nullopt;

    if (dinfh) {
        const auto* inv=find_operation(sym,SymmetryOperationKind::Inversion,2);
        if (!inv) return std::nullopt;
        const auto p=character(wf,*inv,group);
        if (!p.valid) return std::nullopt;
        retention=std::min(retention,p.retention);
        if (retention<options.minimum_subspace_retention ||
            std::abs(std::abs(p.value)-static_cast<double>(d))>options.character_tolerance) return std::nullopt;
        label+=(p.value>=0.0?"g":"u");
    }

    if (m==0) {
        const auto sigma=vertical_reflection(wf,axis);
        const auto s=character(wf,sigma,group);
        if (!s.valid) return std::nullopt;
        retention=std::min(retention,s.retention);
        if (retention<options.minimum_subspace_retention ||
            std::abs(std::abs(s.value)-1.0)>options.character_tolerance) return std::nullopt;
        label+=(s.value>=0.0?"+":"-");
    }
    return label;
}

} // namespace

OrbitalSymmetryResult derive_orbital_symmetry(Wavefunction& wavefunction,
                                               const OrbitalSymmetryOptions& options) {
    OrbitalSymmetryResult result=derive_orbital_symmetry_legacy(wavefunction,options);
    if (wavefunction.orbitals.empty() || wavefunction.basis_count==0 ||
        wavefunction.ao_overlap.size()!=static_cast<std::size_t>(wavefunction.basis_count)*wavefunction.basis_count) {
        return result;
    }

    const MolecularSymmetry sym=analyse_molecular_symmetry(wavefunction);
    if (result.point_group.empty()) result.point_group=sym.point_group;
    if (!sym.available()) return result;

    for (const auto& group:energy_groups(wavefunction,options.degeneracy_tolerance_hartree)) {
        if (!group_unlabelled(wavefunction,group)) continue;
        double retention=1.0;
        std::optional<std::string> label;
        if (sym.point_group=="Oh") {
            label=classify_oh_fallback(wavefunction,sym,group,options,retention);
        } else if (sym.point_group=="Dinfh" || sym.point_group=="Cinfv") {
            label=classify_linear(wavefunction,sym,group,options,retention);
        }
        if (!label) continue;
        DerivedOrbitalSymmetryAssignment evidence;
        evidence.point_group=sym.point_group;
        evidence.label=*label;
        evidence.orbital_indices=group;
        evidence.subspace_retention=retention;
        evidence.centre_bohr=sym.centre_bohr;
        wavefunction.derived_orbital_symmetry_assignments.push_back(std::move(evidence));
        ++result.groups_examined;
        ++result.groups_labelled;
        result.worst_subspace_retention=std::min(result.worst_subspace_retention,retention);
        for (std::size_t i:group) {
            wavefunction.orbitals[i].symmetry=*label;
            wavefunction.orbitals[i].symmetry_provenance=DataProvenance::Derived;
            ++result.orbitals_labelled;
        }
    }
    return result;
}

} // namespace cov
