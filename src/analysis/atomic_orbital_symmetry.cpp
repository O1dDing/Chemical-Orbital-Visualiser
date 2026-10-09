#include "cov/atomic_orbital_symmetry.hpp"
#include "cov/ao_angular_basis.hpp"
#include "cov/local_angular_generators.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <sstream>

namespace cov {
namespace {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using RowMatrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Powers = ao_angular::Powers;
using Polynomial = std::map<Powers,double>;
constexpr std::array<double,9> identity{1,0,0,0,1,0,0,0,1};
constexpr const char* letters[] = {"s","p","d","f","g"};
Polynomial harmonic(int l,int m) {
    Polynomial p;
    for(const auto& term:ao_angular::solid_harmonic(l,m))p[term.powers]+=term.coefficient;
    return p;
}
// Real infinitesimal rotation: x_b d/dx_c - x_c d/dx_b. Differentiation is
// exact on the stored homogeneous polynomials, with no angular sampling.
Polynomial derivative(const Polynomial& p,int axis) {
    const int b=(axis+1)%3,c=(axis+2)%3;
    Polynomial result;
    for(const auto& [powers,value]:p) {
        if(powers[c]) {auto q=powers;--q[c];++q[b];result[q]+=value*powers[c];}
        if(powers[b]) {auto q=powers;--q[b];++q[c];result[q]-=value*powers[b];}
    }
    return result;
}
double sphere_inner_product(const Polynomial& a,const Polynomial& b) {
    double result=0;
    for(const auto& [p,x]:a)for(const auto& [q,y]:b) {
        const Powers r{p[0]+q[0],p[1]+q[1],p[2]+q[2]};
        if(r[0]%2||r[1]%2||r[2]%2)continue;
        result+=x*y*2*std::tgamma((r[0]+1)*.5)*std::tgamma((r[1]+1)*.5)*
            std::tgamma((r[2]+1)*.5)/std::tgamma((r[0]+r[1]+r[2]+3)*.5);
    }
    return result;
}
std::array<Matrix,3> angular_action(int l) {
    const int d=2*l+1;
    std::array<Matrix,3> result;
    std::vector<Polynomial> harmonics;
    for(int m=0;m<d;++m)harmonics.push_back(harmonic(l,m));
    for(int axis=0;axis<3;++axis) {
        result[axis]=Matrix::Zero(d,d);
        for(int col=0;col<d;++col) {
            const auto dp=derivative(harmonics[col],axis);
            for(int row=0;row<d;++row)result[axis](row,col)=sphere_inner_product(harmonics[row],dp);
        }
    }
    return result;
}
Matrix packed_columns(const std::vector<Vector>& columns,Eigen::Index rows) {
    Matrix result(rows,columns.size());
    for(std::size_t j=0;j<columns.size();++j)result.col(j)=columns[j];
    return result;
}
double norm_squared(const Vector& v,const Matrix& s) {return v.dot(s*v);}
struct Reference {
    Matrix orthonormal;
    std::vector<std::size_t> shells;
};
bool orthonormalise(const Matrix& b,const Matrix& s,double rank_tolerance,Matrix& q) {
    if(!b.cols()) {q.resize(s.rows(),0);return true;}
    const Matrix raw=b.transpose()*s*b;
    const Matrix g=(raw+raw.transpose())*.5;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(g);
    if(solver.info()!=Eigen::Success)return false;
    const auto& values=solver.eigenvalues();
    const double limit=rank_tolerance*std::max(1.0,values.cwiseAbs().maxCoeff());
    if(values.minCoeff() < -limit)return false;
    const Eigen::Index rank=(values.array()>limit).count();
    q.resize(s.rows(),rank);
    Eigen::Index col=0;
    for(Eigen::Index i=0;i<values.size();++i)if(values[i]>limit)
        q.col(col++)=(b*solver.eigenvectors().col(i))/std::sqrt(values[i]);
    return q.allFinite();
}
}

AtomicOrbitalSymmetryResult analyse_atomic_orbital_symmetry(
    const Wavefunction& w,const AtomicOrbitalSymmetryOptions& options) {
    AtomicOrbitalSymmetryResult result;
    const auto fail=[&](const char* status,const std::string& detail) {
        result.available=false;result.status=status;result.detail=detail;
        return result;
    };
    if(w.atoms.size()!=1)return fail("not_single_atom","Atomic angular analysis requires exactly one nuclear centre");
    result.applicable=true;result.symmetry_group="SO(3)";
    if(!std::isfinite(options.metric_tolerance)||options.metric_tolerance<=0||
       !std::isfinite(options.squared_residual_tolerance)||options.squared_residual_tolerance<=0||
       !std::isfinite(options.degeneracy_tolerance_hartree)||options.degeneracy_tolerance_hartree<0||
       !std::isfinite(options.relative_rank_tolerance)||options.relative_rank_tolerance<=0)
        return fail("invalid_options","Atomic analysis tolerances are invalid");
    const std::size_t n=w.basis_count;
    if(!n||n>std::numeric_limits<std::size_t>::max()/n||w.ao_overlap.size()!=n*n||w.shells.empty())
        return fail("atomic_evidence_unavailable","Complete AO shell definitions and overlap metric are required");
    Matrix s=Eigen::Map<const RowMatrix>(w.ao_overlap.data(),n,n);
    if(!s.allFinite()||(s-s.transpose()).cwiseAbs().maxCoeff()>options.metric_tolerance)
        return fail("invalid_ao_metric","AO overlap is nonfinite or nonsymmetric");
    s=(s+s.transpose()).eval()*.5;
    Eigen::SelfAdjointEigenSolver<Matrix> metric_solver(s);
    if(metric_solver.info()!=Eigen::Success)return fail("invalid_ao_metric","AO metric spectral analysis failed");
    const double rank_limit=options.relative_rank_tolerance*std::max(1.0,metric_solver.eigenvalues().cwiseAbs().maxCoeff());
    if(metric_solver.eigenvalues().minCoeff() < -rank_limit)
        return fail("invalid_ao_metric","AO overlap has a negative metric direction");
    result.ao_metric_rank=(metric_solver.eigenvalues().array()>rank_limit).count();
    if(!result.ao_metric_rank)return fail("invalid_ao_metric","AO overlap has no resolved positive direction");

    std::array<std::vector<Vector>,5> reference_columns;
    std::array<Reference,5> references;
    std::array<std::array<Matrix,3>,5> actions;
    for(int l=0;l<=4;++l)actions[l]=angular_action(l);
    std::array<Matrix,3> generators{Matrix::Zero(n,n),Matrix::Zero(n,n),Matrix::Zero(n,n)};
    std::vector<bool> covered(n,false);
    for(std::size_t shell_index=0;shell_index<w.shells.size();++shell_index) {
        const auto& shell=w.shells[shell_index];
        const std::size_t offset=shell.basis_offset,count=shell_basis_count(shell);
        if(shell.atom_index!=0||shell.angular_momentum>4||shell.pure>1||offset>n||count>n-offset)
            return fail("invalid_ao_shells","AO shell range, centre, or angular support is invalid");
        for(std::size_t i=offset;i<offset+count;++i) {
            if(covered[i])return fail("invalid_ao_shells","AO shell ranges overlap");
            covered[i]=true;
        }
        Matrix basis(count,count);
        std::array<Matrix,3> shell_actions{Matrix::Zero(count,count),Matrix::Zero(count,count),Matrix::Zero(count,count)};
        std::size_t begin=0;
        for(int l=shell.angular_momentum;l>=0;l-=2) {
            if(shell.pure&&l!=shell.angular_momentum)continue;
            const auto angular=local_angular_generators(shell.angular_momentum,shell.pure!=0,l,identity);
            if(!angular.available||angular.basis_count!=count)
                return fail("angular_basis_unavailable","Analytic angular reference reconstruction failed");
            const std::size_t d=2*l+1;
            basis.block(0,begin,count,d)=Eigen::Map<const RowMatrix>(angular.coefficients.data(),count,d);
            for(int axis=0;axis<3;++axis)shell_actions[axis].block(begin,begin,d,d)=actions[l][axis];
            for(std::size_t m=0;m<d;++m) {
                Vector column=Vector::Zero(n);
                column.segment(offset,count)=basis.col(begin+m);
                reference_columns[l].push_back(std::move(column));
            }
            references[l].shells.push_back(shell_index);
            begin+=d;
        }
        if(begin!=count)return fail("angular_basis_unavailable","Analytic angular references do not span the shell");
        Eigen::FullPivLU<Matrix> lu(basis);
        if(!lu.isInvertible())return fail("angular_basis_unavailable","Analytic shell reference matrix is singular");
        const Matrix inverse=lu.inverse();
        for(int axis=0;axis<3;++axis)generators[axis].block(offset,offset,count,count)=basis*shell_actions[axis]*inverse;
    }
    if(std::find(covered.begin(),covered.end(),false)!=covered.end())
        return fail("invalid_ao_shells","AO definitions omit a basis direction");
    result.generator_metric_error=0;
    const double scale=std::max(1.0,s.cwiseAbs().maxCoeff());
    for(const auto& j:generators) {
        const Matrix error=j.transpose()*s+s*j;
        result.generator_metric_error=std::max(result.generator_metric_error,error.cwiseAbs().maxCoeff()/scale);
    }
    if(result.generator_metric_error>options.metric_tolerance)
        return fail("atomic_metric_not_rotation_invariant","Actual AO metric does not preserve analytic continuous rotations");
    std::size_t total_rank=0;
    for(int l=0;l<=4;++l) {
        if(!orthonormalise(packed_columns(reference_columns[l],n),s,options.relative_rank_tolerance,references[l].orthonormal))
            return fail("angular_reference_unresolved","An angular reference metric is unresolved");
        const std::size_t rank=references[l].orthonormal.cols();
        result.angular_reference_ranks[l]=rank;total_rank+=rank;
        if(rank%std::size_t(2*l+1))return fail("angular_reference_incomplete","Metric-resolved angular rank is not a complete 2l+1 multiple");
        for(int other=0;other<l;++other)if(rank&&references[other].orthonormal.cols()) {
            const Matrix overlap=references[other].orthonormal.transpose()*s*references[l].orthonormal;
            if(overlap.cwiseAbs().maxCoeff()>options.metric_tolerance)
                return fail("angular_references_not_orthogonal","Different angular momenta are not orthogonal in the actual S metric");
        }
    }
    if(total_rank!=result.ao_metric_rank)
        return fail("angular_reference_incomplete","Angular references do not cover the resolved AO metric rank");
    result.available=true;result.status="atomic_angular_evidence_available";
    result.detail="Analytic SO(3) AO action in actual overlap metric; parity=(-1)^l; radial occurrence is not principal n";
    result.orbitals.resize(w.orbitals.size());
    for(std::size_t index=0;index<w.orbitals.size();++index) {
        const auto& mo=w.orbitals[index];auto& assignment=result.orbitals[index];
        assignment.source_orbital_index=index;assignment.partner_status="partners_unresolved";
        if(mo.coefficients.size()!=n||(mo.spin!=Spin::Alpha&&mo.spin!=Spin::Beta)) {
            assignment.status="invalid_source_orbital";continue;
        }
        const Vector c=Eigen::Map<const Vector>(mo.coefficients.data(),n);
        if(!c.allFinite()) {assignment.status="invalid_source_orbital";continue;}
        const double norm=norm_squared(c,s);assignment.source_squared_s_norm=norm;
        if(!std::isfinite(norm)||norm<=rank_limit) {assignment.status="source_metric_norm_unresolved";continue;}
        Vector reconstructed=Vector::Zero(n);double sum=0;int pure=-1;
        for(int l=0;l<=4;++l) {
            const auto& q=references[l].orthonormal;
            if(!q.cols())continue;
            const Vector projected=q*(q.transpose()*(s*c));
            const double weight=std::max(0.0,norm_squared(projected,s))/norm;
            AtomicAngularComponent component;
            component.angular_momentum=l;component.label=letters[l];component.weight=weight;
            component.ao_coefficients.assign(projected.data(),projected.data()+n);
            component.source_shell_indices=references[l].shells;
            assignment.components.push_back(std::move(component));
            reconstructed+=projected;sum+=weight;
            const double loss=std::max(0.0,norm_squared(c-projected,s))/norm;
            if(loss<=options.squared_residual_tolerance) {pure=l;assignment.projection_residual_squared=loss;}
        }
        assignment.reconstruction_residual_squared=std::max(0.0,norm_squared(c-reconstructed,s))/norm;
        assignment.weight_sum_error=std::abs(sum-1);
        assignment.decomposition_verified=assignment.reconstruction_residual_squared<=options.squared_residual_tolerance&&
            assignment.weight_sum_error<=options.metric_tolerance;
        if(!assignment.decomposition_verified) {assignment.status="angular_decomposition_unresolved";continue;}
        if(std::abs(norm-1)>options.metric_tolerance) {assignment.status="source_metric_not_normalized";continue;}
        if(pure<0) {assignment.status="mixed_angular_momentum";continue;}
        assignment.angular_momentum_verified=true;assignment.angular_momentum=pure;
        assignment.label=letters[pure];assignment.inversion_parity=pure%2?-1:1;
        assignment.status="pure_angular_momentum_verified";
        assignment.detail="Unchanged source column is pure "+assignment.label+" in S; principal n is unavailable";
    }

    // Energies only propose candidate partner spans. Continuous generator
    // leakage, source Gram, pure l, and dimension establish the actual evidence.
    for(int spin=0;spin<2;++spin)for(int l=0;l<=4;++l) {
        std::vector<std::size_t> indices;
        for(std::size_t i=0;i<w.orbitals.size();++i)if(result.orbitals[i].angular_momentum_verified&&
            result.orbitals[i].angular_momentum==l&&int(w.orbitals[i].spin)==spin&&std::isfinite(w.orbitals[i].energy_hartree))indices.push_back(i);
        std::stable_sort(indices.begin(),indices.end(),[&](auto a,auto b){return w.orbitals[a].energy_hartree<w.orbitals[b].energy_hartree;});
        const std::size_t d=2*l+1;
        for(std::size_t begin=0;begin<indices.size();) {
            std::size_t end=begin+1;
            if(l)while(end<indices.size()&&w.orbitals[indices[end]].energy_hartree-w.orbitals[indices[begin]].energy_hartree<=options.degeneracy_tolerance_hartree)++end;
            const std::vector<std::size_t> members(indices.begin()+begin,indices.begin()+end);
            const std::size_t k=members.size();Matrix q(n,k);
            for(std::size_t j=0;j<k;++j)q.col(j)=Eigen::Map<const Vector>(w.orbitals[members[j]].coefficients.data(),n);
            const Matrix gram=q.transpose()*s*q;
            const double gram_error=(gram-Matrix::Identity(k,k)).cwiseAbs().maxCoeff();
            double leakage=0;
            for(const auto& generator:generators) {
                const Matrix transformed=generator*q;
                const Matrix residual=transformed-q*(q.transpose()*s*transformed);
                for(std::size_t j=0;j<k;++j)leakage=std::max(leakage,std::max(0.0,norm_squared(residual.col(j),s)));
            }
            const bool closed=k%d==0&&gram_error<=options.metric_tolerance&&leakage<=options.squared_residual_tolerance;
            for(auto index:members) {
                auto& assignment=result.orbitals[index];assignment.containing_members=members;
                assignment.source_gram_error=gram_error;assignment.generator_leakage_squared=leakage;
                assignment.containing_span_closed=closed;
                if(closed) {
                    assignment.representation_multiplicity=k/d;assignment.partner_block_verified=k==d;
                    assignment.partner_status=k==d?"complete_radial_partner_block_verified":"closed_repeated_radial_copies_unresolved";
                    if(k==d)assignment.partner_block_id="atomic:"+std::to_string(spin)+":"+letters[l]+":"+std::to_string(members.front());
                } else assignment.partner_status=gram_error>options.metric_tolerance?"source_metric_not_orthonormal":"radial_partner_span_not_closed";
            }
            begin=end;
        }
        // Number only fully resolved copies in this spin/l set. This conservative
        // gate cannot hide an unnumbered lower occurrence or an equal-energy tie.
        std::vector<std::size_t> representatives;bool ordering_complete=true;
        // A mixed/unresolved source with l content can contain an additional
        // occurrence. Do not number only the conveniently pure subset. Missing
        // columns/energies cannot establish the complete occurrence order.
        for(std::size_t source=0;source<w.orbitals.size();++source) {
            if(int(w.orbitals[source].spin)!=spin)continue;
            const auto& assignment=result.orbitals[source];
            if(assignment.angular_momentum_verified) {
                if(assignment.angular_momentum==l&&!std::isfinite(w.orbitals[source].energy_hartree))ordering_complete=false;
                continue;
            }
            if(assignment.components.empty())ordering_complete=false;
            for(const auto& component:assignment.components)if(component.angular_momentum==l&&
                component.weight>options.squared_residual_tolerance)ordering_complete=false;
        }
        for(auto index:indices) {
            const auto& assignment=result.orbitals[index];
            if(!assignment.partner_block_verified)ordering_complete=false;
            else if(index==assignment.containing_members.front())representatives.push_back(index);
        }
        // Separate energy blocks are not automatically separate radial copies:
        // duplicate or overlapping source columns must not be counted twice.
        if(!indices.empty()) {
            Matrix all(n,indices.size());
            for(std::size_t i=0;i<indices.size();++i)all.col(i)=Eigen::Map<const Vector>(w.orbitals[indices[i]].coefficients.data(),n);
            const Matrix gram=all.transpose()*s*all;
            if((gram-Matrix::Identity(indices.size(),indices.size())).cwiseAbs().maxCoeff()>options.metric_tolerance)
                ordering_complete=false;
        }
        for(std::size_t i=1;i<representatives.size();++i)if(w.orbitals[representatives[i]].energy_hartree-w.orbitals[representatives[i-1]].energy_hartree<=options.degeneracy_tolerance_hartree)ordering_complete=false;
        if(ordering_complete)for(std::size_t i=0;i<representatives.size();++i)
            for(auto index:result.orbitals[representatives[i]].containing_members)result.orbitals[index].radial_copy_ordinal=i+1;
    }
    return result;
}
} // namespace cov
