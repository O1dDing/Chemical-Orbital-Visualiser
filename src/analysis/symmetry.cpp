#include "cov/symmetry.hpp"
#include "cov/point_group_irreps.hpp"
#include <Eigen/SVD>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>;

struct RotationAxis {
    Vec3 axis{};
    int order = 1;
    double error = 0.0;
};

struct ImproperAxis {
    Vec3 axis{};
    int order = 1;
    double error = 0.0;
};

struct ReflectionPlane {
    Vec3 normal{};
    double error = 0.0;
};

Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 scale(const Vec3& a, const double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

double norm2(const Vec3& a) {
    return dot(a, a);
}

double norm(const Vec3& a) {
    return std::sqrt(norm2(a));
}

bool normalize(Vec3& a) {
    const double n = norm(a);
    if (!(n > 1.0e-14) || !std::isfinite(n)) return false;
    a = scale(a, 1.0 / n);
    return true;
}

void canonicalize_axis(Vec3& a) {
    if (!normalize(a)) return;
    for (double value : a) {
        if (std::abs(value) <= 1.0e-12) continue;
        if (value < 0.0) a = scale(a, -1.0);
        break;
    }
}

Vec3 mat_vec(const Mat3& m, const Vec3& v) {
    return {
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    };
}

Mat3 mat_mul(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            double value = 0.0;
            for (int k = 0; k < 3; ++k) value += a[3 * r + k] * b[3 * k + c];
            out[3 * r + c] = value;
        }
    }
    return out;
}

Mat3 identity_matrix() {
    return {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
}

Mat3 inversion_matrix() {
    return {
        -1.0, 0.0, 0.0,
         0.0,-1.0, 0.0,
         0.0, 0.0,-1.0,
    };
}

Mat3 rotation_matrix(Vec3 axis, const double angle) {
    normalize(axis);
    const double x = axis[0];
    const double y = axis[1];
    const double z = axis[2];
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1.0 - c;
    return {
        t*x*x + c,     t*x*y - s*z,   t*x*z + s*y,
        t*x*y + s*z,   t*y*y + c,     t*y*z - s*x,
        t*x*z - s*y,   t*y*z + s*x,   t*z*z + c,
    };
}

Mat3 reflection_matrix(Vec3 normal) {
    normalize(normal);
    const double x = normal[0];
    const double y = normal[1];
    const double z = normal[2];
    return {
        1.0 - 2.0*x*x, -2.0*x*y,       -2.0*x*z,
        -2.0*y*x,       1.0 - 2.0*y*y, -2.0*y*z,
        -2.0*z*x,       -2.0*z*y,       1.0 - 2.0*z*z,
    };
}

Mat3 improper_matrix(const Vec3& axis, const int order) {
    const Mat3 r = rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order));
    const Mat3 h = reflection_matrix(axis);
    return mat_mul(h, r);
}

bool same_axis(const Vec3& a, const Vec3& b, const double cosine_tolerance = 1.0e-7) {
    return std::abs(dot(a, b)) >= 1.0 - cosine_tolerance;
}

void add_axis(std::vector<Vec3>& axes,
              Vec3 axis,
              const std::size_t maximum) {
    if (axes.size() >= maximum) return;
    if (!normalize(axis)) return;
    canonicalize_axis(axis);
    for (const auto& existing : axes) {
        if (same_axis(axis, existing, 2.0e-7)) return;
    }
    axes.push_back(axis);
}

Vec3 molecular_centre(const Wavefunction& wf) {
    Vec3 centre{0.0, 0.0, 0.0};
    double weight_sum = 0.0;
    for (const auto& atom : wf.atoms) {
        const double weight = static_cast<double>(std::max(1, atom.atomic_number));
        centre[0] += weight * atom.x;
        centre[1] += weight * atom.y;
        centre[2] += weight * atom.z;
        weight_sum += weight;
    }
    if (weight_sum > 0.0) centre = scale(centre, 1.0 / weight_sum);
    return centre;
}

std::vector<Vec3> centred_positions(const Wavefunction& wf, const Vec3& centre) {
    std::vector<Vec3> positions;
    positions.reserve(wf.atoms.size());
    for (const auto& atom : wf.atoms) {
        positions.push_back({atom.x - centre[0], atom.y - centre[1], atom.z - centre[2]});
    }
    return positions;
}

bool validate_operation(const Wavefunction& wf,
                        const std::vector<Vec3>& positions,
                        const Mat3& matrix,
                        const double tolerance,
                        std::vector<std::size_t>& permutation,
                        double& max_error) {
    const std::size_t n = wf.atoms.size();
    permutation.assign(n, n);
    std::vector<bool> used(n, false);
    max_error = 0.0;

    const double tolerance2 = tolerance * tolerance;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 transformed = mat_vec(matrix, positions[i]);
        std::size_t best = n;
        double best_d2 = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < n; ++j) {
            if (used[j] || wf.atoms[j].atomic_number != wf.atoms[i].atomic_number) continue;
            const Vec3 delta = sub(transformed, positions[j]);
            const double d2 = norm2(delta);
            if (d2 < best_d2) {
                best_d2 = d2;
                best = j;
            }
        }
        if (best == n || best_d2 > tolerance2 || !std::isfinite(best_d2)) return false;
        used[best] = true;
        permutation[i] = best;
        max_error = std::max(max_error, std::sqrt(best_d2));
    }
    return true;
}

void append_operation(MolecularSymmetry& result,
                      const SymmetryOperationKind kind,
                      const int order,
                      const int power,
                      const Vec3& axis,
                      const Mat3& matrix,
                      std::vector<std::size_t> permutation,
                      const double error) {
    SymmetryOperation op;
    op.kind = kind;
    op.order = order;
    op.power = power;
    op.axis_or_normal = axis;
    op.matrix = matrix;
    op.atom_permutation = std::move(permutation);
    op.max_mapping_error_bohr = error;
    result.operations.push_back(std::move(op));
    result.max_mapping_error_bohr = std::max(result.max_mapping_error_bohr, error);
}

bool is_linear_geometry(const std::vector<Vec3>& positions,
                        const double tolerance,
                        Vec3& axis) {
    std::size_t farthest = positions.size();
    double farthest_norm2 = 0.0;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const double r2 = norm2(positions[i]);
        if (r2 > farthest_norm2) {
            farthest_norm2 = r2;
            farthest = i;
        }
    }
    if (farthest == positions.size() || farthest_norm2 <= tolerance * tolerance) return true;
    axis = positions[farthest];
    normalize(axis);
    for (const auto& p : positions) {
        if (norm(cross(p, axis)) > tolerance) return false;
    }
    canonicalize_axis(axis);
    return true;
}

double determinant(const Mat3&);
Vec3 fixed_axis(const Mat3&);

bool refine_pair_operation(const Wavefunction& w,const std::vector<Vec3>& positions,
                           Mat3& matrix,const double tolerance) {
    const double desired_determinant=determinant(matrix)>0?1.0:-1.0;
    // Rounded anchors can amplify their individual coordinate errors. Use the
    // tentative same-element correspondence to fit ALL atoms. The wider radius
    // is only for proposing correspondences; it never accepts an operation.
    for(int iteration=0;iteration<2;++iteration) {
        std::vector<std::size_t> permutation;double error=0;
        if(!validate_operation(w,positions,matrix,4*tolerance,permutation,error))return false;
        Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();
        for(std::size_t i=0;i<positions.size();++i) {
            const auto& source=positions[i];const auto& target=positions[permutation[i]];
            for(int row=0;row<3;++row)for(int col=0;col<3;++col)covariance(row,col)+=target[row]*source[col];
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> decomposition(covariance,Eigen::ComputeFullU|Eigen::ComputeFullV);
        if(decomposition.info()!=Eigen::Success)return false;
        Eigen::Matrix3d sign=Eigen::Matrix3d::Identity();
        sign(2,2)=desired_determinant*(decomposition.matrixU()*decomposition.matrixV().transpose()).determinant();
        const Eigen::Matrix3d fitted=decomposition.matrixU()*sign*decomposition.matrixV().transpose();
        if(!fitted.allFinite())return false;
        for(int row=0;row<3;++row)for(int col=0;col<3;++col)matrix[3*row+col]=fitted(row,col);
        if(validate_operation(w,positions,matrix,tolerance,permutation,error))return true;
    }
    return false;
}

std::vector<Vec3> candidate_axes(const Wavefunction& w,
                                 const std::vector<Vec3>& positions,
                                 const double tolerance,const bool refine) {
    std::vector<Vec3> axes;
    const std::size_t n=positions.size();
    // An orthogonal symmetry is determined by the images of any two independent
    // centred vectors and the sign of their normal. Enumerating those atom-pair
    // images is complete for a nonlinear geometry; a truncated collection of
    // atom/pair/triangle axes is not. Maximal area improves frame conditioning.
    std::size_t first=n,second=n;double area=0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j) {
        const double candidate=norm2(cross(positions[i],positions[j]));
        if(candidate>area) {area=candidate;first=i;second=j;}
    }
    if(first==n)return axes;
    const double r1=norm(positions[first]),r2=norm(positions[second]);
    Vec3 x=scale(positions[first],1/r1);
    Vec3 y=sub(positions[second],scale(x,dot(x,positions[second])));
    if(!normalize(y))return axes;const Vec3 z=cross(x,y);
    const std::array<Vec3,3> source{x,y,z};
    const double source_cosine=dot(positions[first],positions[second])/(r1*r2);
    // This is only a proposal filter. Every resulting operation still has to
    // map all original atoms within the unchanged recorded geometry tolerance.
    const double angle_proposal_tolerance=2*tolerance*(1/r1+1/r2);
    const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max({1.0,r1,r2});
    for(std::size_t i=0;i<n;++i) {
        if(w.atoms[i].atomic_number!=w.atoms[first].atomic_number)continue;
        const double a=norm(positions[i]);if(std::abs(a-r1)>tolerance+roundoff||a<=roundoff)continue;
        for(std::size_t j=0;j<n;++j) {
            if(i==j||w.atoms[j].atomic_number!=w.atoms[second].atomic_number)continue;
            const double b=norm(positions[j]);if(std::abs(b-r2)>tolerance+roundoff||b<=roundoff)continue;
            if(std::abs(dot(positions[i],positions[j])/(a*b)-source_cosine)>angle_proposal_tolerance)continue;
            Vec3 u=scale(positions[i],1/a),v=sub(positions[j],scale(u,dot(u,positions[j])));
            if(!normalize(v))continue;const Vec3 normal=cross(u,v);
            for(int sign:{1,-1}) {
                const std::array<Vec3,3> target{u,v,scale(normal,double(sign))};
                Mat3 matrix{};
                for(int row=0;row<3;++row)for(int col=0;col<3;++col)for(int k=0;k<3;++k)
                    matrix[3*row+col]+=target[k][row]*source[k][col];
                if(refine) {
                    if(!refine_pair_operation(w,positions,matrix,tolerance))continue;
                } else {
                    std::vector<std::size_t> permutation;double error=0;
                    if(!validate_operation(w,positions,matrix,tolerance,permutation,error))continue;
                }
                // A proper operation fixes its rotation axis; an improper
                // operation negates it. Identity/inversion have no unique axis.
                if(determinant(matrix)<0)for(auto& value:matrix)value=-value;
                const auto axis=fixed_axis(matrix);
                add_axis(axes,axis,std::numeric_limits<std::size_t>::max());
            }
        }
    }
    return axes;
}

int count_axes_with_order_at_least(const std::vector<RotationAxis>& axes, const int order) {
    int count = 0;
    // C5 does not contain C4 or C3. Count the actual rotation subgroups.
    for (const auto& axis : axes) if (axis.order % order == 0) ++count;
    return count;
}

int maximum_rotation_order(const std::vector<RotationAxis>& axes) {
    int order = 1;
    for (const auto& axis : axes) order = std::max(order, axis.order);
    return order;
}

int maximum_improper_order(const std::vector<ImproperAxis>& axes) {
    int order = 1;
    for (const auto& axis : axes) order = std::max(order, axis.order);
    return order;
}

const RotationAxis* principal_rotation_axis(const std::vector<RotationAxis>& axes,
                                            const int order,
                                            const std::vector<ImproperAxis>& improper) {
    const RotationAxis* best=nullptr;int best_improper=1;
    for (const auto& axis : axes) if (axis.order == order) {
        int aligned=1;for(const auto& s:improper)if(same_axis(axis.axis,s.axis))aligned=std::max(aligned,s.order);
        if(!best||aligned>best_improper||(aligned==best_improper&&axis.error<best->error)){
            best=&axis;best_improper=aligned;
        }
    }
    return best;
}

bool has_plane_relation(const std::vector<ReflectionPlane>& planes,
                        const Vec3& axis,
                        const bool normal_parallel_axis,
                        const double angular_tolerance) {
    for (const auto& plane : planes) {
        const double alignment = std::abs(dot(plane.normal, axis));
        if (normal_parallel_axis) {
            if (alignment > 1.0 - angular_tolerance*angular_tolerance) return true;
        } else {
            if (alignment < angular_tolerance) return true;
        }
    }
    return false;
}

int perpendicular_c2_count(const std::vector<RotationAxis>& axes,
                           const Vec3& principal,
                           const double angular_tolerance) {
    int count = 0;
    for (const auto& axis : axes) {
        if (axis.order % 2 != 0) continue;
        if (std::abs(dot(axis.axis, principal)) < angular_tolerance) ++count;
    }
    return count;
}

double matrix_error(const Mat3& a,const Mat3& b) {
    double e=0;for(int k=0;k<9;++k)e=std::max(e,std::abs(a[k]-b[k]));return e;
}
Mat3 transpose(const Mat3& a) {return {a[0],a[3],a[6],a[1],a[4],a[7],a[2],a[5],a[8]};}
Mat3 negative(Mat3 a) {for(auto& x:a)x=-x;return a;}
double determinant(const Mat3& a) {
    return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);
}
Vec3 fixed_axis(const Mat3& a) {
    const Vec3 rows[3]{{a[0]-1,a[1],a[2]},{a[3],a[4]-1,a[5]},{a[6],a[7],a[8]-1}};
    Vec3 axis{};double best=0;for(int i=0;i<3;++i)for(int j=i+1;j<3;++j){
        const auto candidate=cross(rows[i],rows[j]);if(norm2(candidate)>best){best=norm2(candidate);axis=candidate;}}
    canonicalize_axis(axis);return axis;
}

// Candidate axes extracted from rounded coordinates are separately accurate,
// but their floating-point generators need not obey exact group relations.
// Fit a common frame, construct exact finite-group matrices in that frame,
// and validate EACH matrix against the original, untouched atom coordinates.
// This changes the symmetry action only within its recorded geometry tolerance.
bool complete_consistent_operations(MolecularSymmetry& result,const Wavefunction& w,
        const std::vector<Vec3>& positions,const std::vector<RotationAxis>& rotations,
        const std::vector<ImproperAxis>& improper,const std::vector<ReflectionPlane>& planes,
        double angular_tolerance) {
    const auto expected=finite_point_group_order(result.point_group);
    if(!expected)return false;
    std::vector<Mat3> generators;
    const auto best_axis=[&](int n)->const RotationAxis*{
        const RotationAxis* best=nullptr;for(const auto& axis:rotations)
            if(axis.order==n&&(!best||axis.error<best->error))best=&axis;return best;};
    const std::string& group=result.point_group;
    if(group=="C1") {}
    else if(group=="Ci")generators.push_back(inversion_matrix());
    else if(group=="Cs"){if(planes.empty())return false;const auto best=std::min_element(planes.begin(),planes.end(),[](const auto& a,const auto& b){return a.error<b.error;});generators.push_back(reflection_matrix(best->normal));}
    else if(group=="T"||group=="Td"||group=="Th"||group=="O"||group=="Oh"){
        const bool octahedral=group[0]=='O';const auto* first=best_axis(octahedral?4:2);if(!first)return false;
        Vec3 x=first->axis,y{};double best=std::numeric_limits<double>::infinity();
        for(const auto& a:rotations)if(a.order==(octahedral?4:2)&&std::abs(dot(x,a.axis))<angular_tolerance&&a.error<best){best=a.error;y=a.axis;}
        y=sub(y,scale(x,dot(x,y)));if(!normalize(y))return false;const auto z=cross(x,y);
        const Mat3 frame{x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2]};
        const Mat3 cycle{0,0,1,1,0,0,0,1,0},halfturn{1,0,0,0,-1,0,0,0,-1};
        generators.push_back(mat_mul(mat_mul(frame,cycle),transpose(frame)));
        generators.push_back(mat_mul(mat_mul(frame,halfturn),transpose(frame)));
        if(octahedral)generators.push_back(rotation_matrix(z,kPi/2));
        if(group=="Td")generators.push_back(mat_mul(mat_mul(frame,Mat3{0,1,0,1,0,0,0,0,1}),transpose(frame)));
        if(group=="Th"||group=="Oh")generators.push_back(inversion_matrix());
    }else if(group=="I"||group=="Ih"){
        const auto* first=best_axis(5);if(!first)return false;const auto z=first->axis;Vec3 x{};double best=std::numeric_limits<double>::infinity();
        const double cosine=1/std::sqrt(5.0);
        for(const auto& a:rotations)if(a.order==5&&std::abs(std::abs(dot(z,a.axis))-cosine)<angular_tolerance&&a.error<best){
            best=a.error;auto adjacent=a.axis;if(dot(adjacent,z)<0)adjacent=scale(adjacent,-1);x=sub(adjacent,scale(z,dot(adjacent,z)));
        }
        if(!normalize(x))return false;
        // Adjacent vertex axes of the icosahedron have exact cos(theta)=1/sqrt5.
        const auto adjacent=add(scale(x,2/std::sqrt(5.0)),scale(z,cosine));
        generators={rotation_matrix(z,2*kPi/5),rotation_matrix(adjacent,2*kPi/5)};
        if(group=="Ih")generators.push_back(inversion_matrix());
    }else{
        std::size_t k=1;int n=0;while(k<group.size()&&group[k]>='0'&&group[k]<='9'){n=10*n+group[k]-'0';++k;}
        if(n<2)return false;const char family=group[0],suffix=k<group.size()?group[k]:0;
        Vec3 z{};
        if(family=='S'){const ImproperAxis* best=nullptr;for(const auto& a:improper)if(a.order==n&&(!best||a.error<best->error))best=&a;if(!best)return false;z=best->axis;}
        else{const auto* principal=principal_rotation_axis(rotations,n,improper);if(!principal)return false;z=principal->axis;}
        if(family=='S')generators.push_back(improper_matrix(z,n));
        else if(family=='D'&&suffix=='d')generators.push_back(improper_matrix(z,2*n));
        else generators.push_back(rotation_matrix(z,2*kPi/n));
        if(family=='D'){
            Vec3 x{};double best=std::numeric_limits<double>::infinity();
            for(const auto& a:rotations)if(a.order%2==0&&std::abs(dot(z,a.axis))<angular_tolerance&&a.error<best){best=a.error;x=a.axis;}
            x=sub(x,scale(z,dot(z,x)));if(!normalize(x))return false;generators.push_back(rotation_matrix(x,kPi));
        }else if(family=='C'&&suffix=='v'){
            Vec3 normal{};double best=std::numeric_limits<double>::infinity();
            for(const auto& p:planes)if(std::abs(dot(z,p.normal))<angular_tolerance&&p.error<best){best=p.error;normal=p.normal;}
            normal=sub(normal,scale(z,dot(z,normal)));if(!normalize(normal))return false;generators.push_back(reflection_matrix(normal));
        }
        if(suffix=='h')generators.push_back(reflection_matrix(z));
    }
    MolecularSymmetry completed=result;completed.operations.clear();completed.max_mapping_error_bohr=0;
    std::vector<Mat3> matrices{identity_matrix()};
    for(std::size_t i=0;i<matrices.size();++i)for(const auto& g:generators){
        const auto p=mat_mul(matrices[i],g);
        if(std::any_of(matrices.begin(),matrices.end(),[&](const auto& a){return matrix_error(a,p)<1e-8;}))continue;
        if(matrices.size()>=expected)return false;matrices.push_back(p);
    }
    if(matrices.size()!=expected)return false;
    for(const auto& matrix:matrices){
        std::vector<std::size_t> permutation;double error=0;
        if(!validate_operation(w,positions,matrix,result.tolerance_bohr,permutation,error))return false;
        SymmetryOperationKind kind;int order=1,power=1;Vec3 axis{};
        const double det=determinant(matrix),trace=matrix[0]+matrix[4]+matrix[8];
        if(matrix_error(matrix,identity_matrix())<1e-8){kind=SymmetryOperationKind::Identity;power=0;axis={0,0,1};}
        else if(matrix_error(matrix,inversion_matrix())<1e-8){kind=SymmetryOperationKind::Inversion;order=2;}
        else if(det<0&&std::abs(trace-1)<1e-8){kind=SymmetryOperationKind::Reflection;order=2;axis=fixed_axis(negative(matrix));}
        else {
            kind=det>0?SymmetryOperationKind::ProperRotation:SymmetryOperationKind::ImproperRotation;
            axis=fixed_axis(det>0?matrix:negative(matrix));
            // order/power describe the rotational component C_n^k in an
            // improper operation C_n^k sigma_h (S3 has matrix period 6 but
            // rotational order 3), matching the existing metadata contract.
            const auto component=det>0?matrix:mat_mul(matrix,reflection_matrix(axis));
            auto p=component;for(order=1;order<=int(expected);++order){if(matrix_error(p,identity_matrix())<1e-8)break;p=mat_mul(p,component);}
            if(order>int(expected))return false;
            const double rotation_trace=component[0]+component[4]+component[8];
            const double angle=std::acos(std::clamp((rotation_trace-1)/2,-1.0,1.0));power=int(std::round(angle*order/(2*kPi)));
            const auto skew=Vec3{component[7]-component[5],component[2]-component[6],component[3]-component[1]};if(dot(skew,axis)<0)power=order-power;
        }
        append_operation(completed,kind,order,power,axis,matrix,std::move(permutation),error);
    }
    // The real chemical table also checks the exact closure and class counts;
    // failure is an unavailable group, never a fabricated downgrade to C1.
    if(!finite_point_group_irreps(completed).valid)return false;
    result=std::move(completed);return true;
}

} // namespace

static MolecularSymmetry analyse_molecular_symmetry_impl(const Wavefunction& wavefunction,
                                             const SymmetryOptions& options,const bool refine) {
    MolecularSymmetry result;
    if (wavefunction.atoms.empty()) return result;

    result.centre_bohr = molecular_centre(wavefunction);
    const auto positions = centred_positions(wavefunction, result.centre_bohr);
    for (const auto& p : positions) {
        result.molecular_radius_bohr = std::max(result.molecular_radius_bohr, norm(p));
    }
    result.tolerance_bohr = std::max(options.absolute_tolerance_bohr,
                                     options.relative_tolerance * std::max(1.0, result.molecular_radius_bohr));

    std::vector<std::size_t> permutation;
    double error = 0.0;
    if (validate_operation(wavefunction, positions, identity_matrix(), result.tolerance_bohr,
                           permutation, error)) {
        append_operation(result, SymmetryOperationKind::Identity, 1, 0,
                         {0.0, 0.0, 1.0}, identity_matrix(), std::move(permutation), error);
    }

    if (wavefunction.atoms.size() == 1u) {
        result.linear = true;
        result.point_group = "Kh";
        return result;
    }

    Vec3 linear_axis{0.0, 0.0, 1.0};
    result.linear = is_linear_geometry(positions, result.tolerance_bohr, linear_axis);

    permutation.clear();
    error = 0.0;
    if (validate_operation(wavefunction, positions, inversion_matrix(), result.tolerance_bohr,
                           permutation, error)) {
        result.has_inversion = true;
        append_operation(result, SymmetryOperationKind::Inversion, 2, 1,
                         {0.0, 0.0, 0.0}, inversion_matrix(), std::move(permutation), error);
    }

    if (result.linear) {
        result.point_group = result.has_inversion ? "Dinfh" : "Cinfv";
        return result;
    }

    const auto axes = candidate_axes(wavefunction,positions,result.tolerance_bohr,refine);
    std::vector<RotationAxis> rotation_axes;
    std::vector<ImproperAxis> improper_axes;
    std::vector<ReflectionPlane> planes;

    for (const Vec3& axis : axes) {
        int best_order = 1;
        double best_error = 0.0;
        std::vector<std::size_t> best_permutation;
        Mat3 best_matrix = identity_matrix();

        for (int order = std::max(2, options.maximum_rotation_order); order >= 2; --order) {
            const Mat3 matrix = rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order));
            permutation.clear();
            error = 0.0;
            if (validate_operation(wavefunction, positions, matrix, result.tolerance_bohr,
                                   permutation, error)) {
                best_order = order;
                best_error = error;
                best_permutation = permutation;
                best_matrix = matrix;
                break;
            }
        }
        if (best_order >= 2) {
            bool duplicate = false;
            for (const auto& existing : rotation_axes) {
                if (same_axis(axis, existing.axis)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                rotation_axes.push_back({axis, best_order, best_error});
                append_operation(result, SymmetryOperationKind::ProperRotation,
                                 best_order, 1, axis, best_matrix,
                                 std::move(best_permutation), best_error);
            }
        }

        int best_improper = 1;
        double best_improper_error = 0.0;
        std::vector<std::size_t> best_improper_permutation;
        Mat3 best_improper_matrix = identity_matrix();
        for (int order = std::max(3, options.maximum_rotation_order); order >= 3; --order) {
            const Mat3 matrix = improper_matrix(axis, order);
            permutation.clear();
            error = 0.0;
            if (validate_operation(wavefunction, positions, matrix, result.tolerance_bohr,
                                   permutation, error)) {
                best_improper = order;
                best_improper_error = error;
                best_improper_permutation = permutation;
                best_improper_matrix = matrix;
                break;
            }
        }
        if (best_improper >= 3) {
            bool duplicate = false;
            for (const auto& existing : improper_axes) {
                if (same_axis(axis, existing.axis)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                improper_axes.push_back({axis, best_improper, best_improper_error});
                append_operation(result, SymmetryOperationKind::ImproperRotation,
                                 best_improper, 1, axis, best_improper_matrix,
                                 std::move(best_improper_permutation), best_improper_error);
            }
        }

        const Mat3 mirror = reflection_matrix(axis);
        permutation.clear();
        error = 0.0;
        if (validate_operation(wavefunction, positions, mirror, result.tolerance_bohr,
                               permutation, error)) {
            bool duplicate = false;
            for (const auto& existing : planes) {
                if (same_axis(axis, existing.normal)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                planes.push_back({axis, error});
                append_operation(result, SymmetryOperationKind::Reflection,
                                 2, 1, axis, mirror, std::move(permutation), error);
            }
        }
    }

    const int max_order = maximum_rotation_order(rotation_axes);
    const int max_improper = maximum_improper_order(improper_axes);
    const int c5_axes = count_axes_with_order_at_least(rotation_axes, 5);
    const int c4_axes = count_axes_with_order_at_least(rotation_axes, 4);
    const int c3_axes = count_axes_with_order_at_least(rotation_axes, 3);
    const int c2_axes = count_axes_with_order_at_least(rotation_axes, 2);
    const double angular_tolerance=std::max(2.0e-5,
        2*result.tolerance_bohr/std::max(1.0,result.molecular_radius_bohr));

    if (c5_axes >= 6 && c3_axes >= 10) {
        result.point_group = result.has_inversion ? "Ih" : "I";
    } else if (c4_axes >= 3) {
        result.point_group = result.has_inversion ? "Oh" : "O";
    } else if (c3_axes >= 4 && c2_axes >= 3) {
        if (result.has_inversion) result.point_group = "Th";
        else result.point_group = planes.size() >= 6u ? "Td" : "T";
    } else if (max_order >= 2) {
        const RotationAxis* principal = principal_rotation_axis(rotation_axes, max_order,improper_axes);
        if (principal) {
            const bool horizontal = has_plane_relation(planes, principal->axis, true,angular_tolerance);
            const bool vertical = has_plane_relation(planes, principal->axis, false,angular_tolerance);
            const int perpendicular_c2 = perpendicular_c2_count(rotation_axes, principal->axis,angular_tolerance);
            const bool dihedral = perpendicular_c2 >= max_order;
            if (dihedral) {
                if (horizontal) result.point_group = "D" + std::to_string(max_order) + "h";
                else if (vertical) result.point_group = "D" + std::to_string(max_order) + "d";
                else result.point_group = "D" + std::to_string(max_order);
            } else if (horizontal) {
                result.point_group = "C" + std::to_string(max_order) + "h";
            } else if (vertical) {
                result.point_group = "C" + std::to_string(max_order) + "v";
            } else if (max_improper > max_order && max_improper > 2) {
                result.point_group = "S" + std::to_string(max_improper);
            } else {
                result.point_group = "C" + std::to_string(max_order);
            }
        }
    }

    if (result.point_group.empty()) {
        if (result.has_inversion) result.point_group = "Ci";
        else if (!planes.empty()) result.point_group = "Cs";
        else result.point_group = "C1";
    }

    if(!complete_consistent_operations(result,wavefunction,positions,rotation_axes,
            improper_axes,planes,angular_tolerance)){
        result.point_group.clear();result.operations.clear();
    }

    return result;
}

MolecularSymmetry analyse_molecular_symmetry(const Wavefunction& wavefunction,
                                             const SymmetryOptions& options) {
    auto direct=analyse_molecular_symmetry_impl(wavefunction,options,false);
    if(wavefunction.atoms.size()<2||direct.linear)return direct;
    auto fitted=analyse_molecular_symmetry_impl(wavefunction,options,true);
    // Least-squares fitting can find valid evidence that noisy anchors missed,
    // but can also move away from an already valid minimax frame. Keep both
    // independently verified full groups; fitting never discards direct proof.
    if(!direct.available())return fitted;
    if(!fitted.available())return direct;
    if(fitted.operations.size()>direct.operations.size()||
       (fitted.operations.size()==direct.operations.size()&&fitted.max_mapping_error_bohr<direct.max_mapping_error_bohr))return fitted;
    return direct;
}

void derive_point_group_from_geometry(Wavefunction& wavefunction,
                                      const SymmetryOptions& options) {
    if (!wavefunction.point_group_detected.empty() &&
        (wavefunction.point_group_detected_provenance == DataProvenance::Producer ||
         wavefunction.point_group_provenance == DataProvenance::Producer)) return;
    const MolecularSymmetry symmetry = analyse_molecular_symmetry(wavefunction, options);
    if (!symmetry.available()) return;
    wavefunction.point_group_detected = symmetry.point_group;
    // Geometry says nothing about the symmetry settings actually used by a
    // producer calculation. Keep that separately supplied field unchanged.
    if (!wavefunction.point_group_used.empty() &&
        wavefunction.point_group_used_provenance==DataProvenance::Unavailable &&
        wavefunction.point_group_provenance==DataProvenance::Producer)
        wavefunction.point_group_used_provenance=DataProvenance::Producer;
    wavefunction.point_group_provenance = DataProvenance::Derived;
    wavefunction.point_group_detected_provenance = DataProvenance::Derived;
}

} // namespace cov
