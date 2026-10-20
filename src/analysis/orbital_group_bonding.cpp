#include "cov/orbital_group_bonding.hpp"
#include "cov/molecule_style.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

namespace cov { namespace {
using Matrix=std::vector<double>;
struct Eigen {Matrix values,vectors;double residual=0;bool converged=false;};
double maximum(const Matrix& a){double m=0;for(double x:a)m=std::max(m,std::abs(x));return m;}
Matrix identity(std::size_t n){Matrix a(n*n);for(std::size_t i=0;i<n;++i)a[i*n+i]=1;return a;}
Matrix multiply(const Matrix& a,const Matrix& b,std::size_t n){
    Matrix c(n*n);for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<n;++k)
        for(std::size_t j=0;j<n;++j)c[i*n+j]+=a[i*n+k]*b[k*n+j];return c;
}
Matrix transpose(const Matrix& a,std::size_t n){Matrix b(n*n);for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)b[i*n+j]=a[j*n+i];return b;}
Matrix congruence(const Matrix& a,const Matrix& x,std::size_t n){return multiply(transpose(x,n),multiply(a,x,n),n);}
Eigen diagonalize(Matrix a,std::size_t n){
    const Matrix original=a;Eigen e;e.vectors=identity(n);
    const double tolerance=1e-14*std::max(1.0,maximum(a));
    for(unsigned sweep=0;sweep<80;++sweep){
        double off=0;
        for(std::size_t p=0;p<n;++p)for(std::size_t q=p+1;q<n;++q){
            const double apq=a[p*n+q];off=std::max(off,std::abs(apq));if(std::abs(apq)<=tolerance)continue;
            const double theta=.5*std::atan2(2*apq,a[q*n+q]-a[p*n+p]);
            const double c=std::cos(theta),s=std::sin(theta),app=a[p*n+p],aqq=a[q*n+q];
            for(std::size_t k=0;k<n;++k)if(k!=p && k!=q){
                const double kp=a[k*n+p],kq=a[k*n+q];
                a[k*n+p]=a[p*n+k]=c*kp-s*kq;a[k*n+q]=a[q*n+k]=s*kp+c*kq;
            }
            a[p*n+p]=c*c*app-2*c*s*apq+s*s*aqq;a[q*n+q]=s*s*app+2*c*s*apq+c*c*aqq;
            a[p*n+q]=a[q*n+p]=0;
            for(std::size_t k=0;k<n;++k){const double kp=e.vectors[k*n+p],kq=e.vectors[k*n+q];e.vectors[k*n+p]=c*kp-s*kq;e.vectors[k*n+q]=s*kp+c*kq;}
        }
        if(off<=tolerance){e.converged=true;break;}
    }
    e.values.resize(n);for(std::size_t i=0;i<n;++i)e.values[i]=a[i*n+i];
    const auto av=multiply(original,e.vectors,n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)e.residual=std::max(e.residual,std::abs(av[i*n+j]-e.vectors[i*n+j]*e.values[j]));
    e.converged=e.converged && e.residual<1e-10*std::max(1.0,maximum(original));return e;
}
Matrix spectral(const Eigen& e,std::size_t n,bool inverse){
    Matrix r(n*n);for(std::size_t k=0;k<n;++k){const double f=inverse?1/std::sqrt(e.values[k]):std::sqrt(e.values[k]);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)r[i*n+j]+=e.vectors[i*n+k]*f*e.vectors[j*n+k];}return r;
}
bool finite(const Matrix& a){return std::all_of(a.begin(),a.end(),[](double x){return std::isfinite(x);});}
FixedBondingScope scope_from_pairs(const Wavefunction& w,const std::vector<std::pair<std::size_t,std::size_t>>& pairs,const std::string& id,const std::string& source){
    FixedBondingScope s;s.id=id;s.source=source;s.atom_count=w.atoms.size();std::set<std::pair<std::size_t,std::size_t>> unique;
    for(auto [a,b]:pairs){if(a>=s.atom_count || b>=s.atom_count || a==b){++s.invalid_edges;continue;}unique.emplace(std::min(a,b),std::max(a,b));}
    s.edges.assign(unique.begin(),unique.end());s.status=s.invalid_edges?"invalid_scope":"available";return s;
}
} // namespace

FixedBondingScope make_fixed_bonding_scope(const Wavefunction& w){
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    for(const auto& b:analyse_bonds(w))pairs.emplace_back(b.atom_a,b.atom_b);
    return scope_from_pairs(w,pairs,"whole_skeleton","wavefunction_density_geometry_connectivity");
}
FixedBondingScope make_fixed_bonding_scope(const Wavefunction& w,const InteractionGraph& graph){
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    for(const auto& e:graph.edges)if(e.strength==InteractionStrength::StrongConnectivity &&
        (e.kind==InteractionKind::CovalentConnectivity || e.kind==InteractionKind::CoordinationContact))pairs.emplace_back(e.atom_a,e.atom_b);
    return scope_from_pairs(w,pairs,"whole_skeleton","validated_interaction_graph");
}
FixedBondingScope make_fixed_bonding_scope(const Wavefunction& w,const std::vector<std::pair<std::size_t,std::size_t>>& pairs,const std::string& id){
    return scope_from_pairs(w,pairs,id,"explicit_fixed_atom_pairs");
}

OrbitalGroupBondingResult analyse_orbital_group_bonding(const Wavefunction& w,const FixedBondingScope& scope,
    const std::vector<std::size_t>& members,const OrbitalGroupBondingOptions& options,const std::vector<double>& supplied_occupation){
    OrbitalGroupBondingResult r;r.scope_id=scope.id;r.scope_source=scope.source;r.edges=scope.edges;r.source_members=members;r.dimension=members.size();
    const std::size_t g=members.size(),n=w.basis_count;
    const auto fail=[&](OrbitalGroupBondingStatus s,const char* reason){r.status=s;r.detail=reason;return r;};
    if(w.atoms.size()<2)return fail(OrbitalGroupBondingStatus::NotApplicable,"single_atom_has_no_interatomic_scope");
    if(scope.status!="available" || scope.atom_count!=w.atoms.size() || scope.invalid_edges)
        return fail(OrbitalGroupBondingStatus::Unavailable,"invalid_or_unavailable_fixed_scope");
    if(scope.edges.empty())return fail(OrbitalGroupBondingStatus::NotApplicable,"fixed_scope_has_no_supported_atom_pairs");
    if(!g || (options.expected_dimension && options.expected_dimension!=g))
        return fail(OrbitalGroupBondingStatus::IncompleteGroup,"complete_group_dimension_not_supplied");
    if(!n || w.ao_overlap.size()!=n*n || !finite(w.ao_overlap))
        return fail(OrbitalGroupBondingStatus::Unavailable,"complete_finite_ao_metric_unavailable");
    if(!std::isfinite(options.absolute_error) || !std::isfinite(options.relative_error) || !std::isfinite(options.metric_rank_relative_tolerance) ||
        !(options.absolute_error>0) || !(options.relative_error>=0) || !(options.metric_rank_relative_tolerance>0))
        return fail(OrbitalGroupBondingStatus::Unavailable,"invalid_numerical_tolerances");
    std::set<std::size_t> unique;
    Matrix coefficients(n*g);Spin spin=Spin::Alpha;double lo=std::numeric_limits<double>::infinity(),hi=-lo;
    for(std::size_t k=0;k<g;++k){
        const auto index=members[k];
        if(index>=w.orbitals.size() || !unique.insert(index).second)
            return fail(OrbitalGroupBondingStatus::IncompleteGroup,"invalid_or_duplicate_source_member");
        const auto& mo=w.orbitals[index];if(!k)spin=mo.spin;
        if(mo.spin!=spin)return fail(OrbitalGroupBondingStatus::IncompleteGroup,"physical_spin_blocks_must_be_analysed_separately");
        if(mo.coefficients.size()!=n || !finite(mo.coefficients))return fail(OrbitalGroupBondingStatus::Unavailable,"complete_source_column_unavailable");
        for(std::size_t mu=0;mu<n;++mu)coefficients[mu*g+k]=mo.coefficients[mu];
        lo=std::min(lo,mo.energy_hartree);hi=std::max(hi,mo.energy_hartree);
    }
    r.source_energy_span_hartree=hi-lo;
    const auto missing=std::numeric_limits<std::size_t>::max();std::vector<std::size_t> atom_of(n,missing);
    std::vector<std::vector<std::size_t>> atom_basis(w.atoms.size());
    for(const auto& shell:w.shells){
        const auto count=shell_basis_count(shell);
        if(shell.atom_index>=w.atoms.size() || shell.basis_offset+count>n)return fail(OrbitalGroupBondingStatus::Unavailable,"invalid_ao_atom_mapping");
        for(std::size_t k=0;k<count;++k){const auto mu=shell.basis_offset+k;if(atom_of[mu]!=missing)return fail(OrbitalGroupBondingStatus::Unavailable,"overlapping_ao_atom_mapping");atom_of[mu]=shell.atom_index;atom_basis[shell.atom_index].push_back(mu);}
    }
    if(std::find(atom_of.begin(),atom_of.end(),missing)!=atom_of.end())return fail(OrbitalGroupBondingStatus::Unavailable,"incomplete_ao_atom_mapping");
    Matrix sc(n*g);
    for(std::size_t mu=0;mu<n;++mu)for(std::size_t nu=0;nu<n;++nu){
        const double s=w.ao_overlap[mu*n+nu];
        r.hermiticity_error=std::max(r.hermiticity_error,std::abs(s-w.ao_overlap[nu*n+mu]));
        for(std::size_t k=0;k<g;++k)sc[mu*g+k]+=s*coefficients[nu*g+k];
    }
    if(r.hermiticity_error>1e-8*std::max(1.0,maximum(w.ao_overlap)))return fail(OrbitalGroupBondingStatus::Unavailable,"ao_metric_not_hermitian");
    r.source_metric.assign(g*g,0);
    for(std::size_t i=0;i<g;++i)for(std::size_t j=0;j<g;++j)for(std::size_t mu=0;mu<n;++mu)r.source_metric[i*g+j]+=coefficients[mu*g+i]*sc[mu*g+j];
    if(!finite(r.source_metric))return fail(OrbitalGroupBondingStatus::Unavailable,"source_metric_not_finite");
    for(std::size_t i=0;i<g;++i)for(std::size_t j=0;j<g;++j)r.metric_error=std::max(r.metric_error,std::abs(r.source_metric[i*g+j]-(i==j?1.:0.)));
    for(std::size_t i=0;i<g;++i)for(std::size_t j=i+1;j<g;++j)r.source_metric[i*g+j]=r.source_metric[j*g+i]=.5*(r.source_metric[i*g+j]+r.source_metric[j*g+i]);
    const auto metric=diagonalize(r.source_metric,g);
    if(!metric.converged)return fail(OrbitalGroupBondingStatus::Unavailable,"metric_eigensolver_not_converged");
    const double rank_floor=options.metric_rank_relative_tolerance*std::max(1.0,maximum(metric.values));
    for(double value:metric.values)if(value>rank_floor)++r.metric_rank;
    if(r.metric_rank!=g)return fail(OrbitalGroupBondingStatus::IncompleteGroup,"source_group_metric_rank_deficient");
    r.normalization=spectral(metric,g,true);const auto metric_root=spectral(metric,g,false);
    Matrix source_operator(g*g);
    for(const auto& [a,b]:scope.edges){
        if(a>=w.atoms.size() || b>=w.atoms.size() || a==b)return fail(OrbitalGroupBondingStatus::Unavailable,"invalid_fixed_atom_pair");
        for(const auto mu:atom_basis[a])for(const auto nu:atom_basis[b]){
            const double s=.5*(w.ao_overlap[mu*n+nu]+w.ao_overlap[nu*n+mu]);
            for(std::size_t i=0;i<g;++i)for(std::size_t j=0;j<g;++j)
                source_operator[i*g+j]+=s*(coefficients[mu*g+i]*coefficients[nu*g+j]+coefficients[nu*g+i]*coefficients[mu*g+j]);
        }
    }
    r.operator_matrix=congruence(source_operator,r.normalization,g);
    if(!finite(r.operator_matrix))return fail(OrbitalGroupBondingStatus::Unavailable,"group_operator_not_finite");
    for(std::size_t i=0;i<g;++i)for(std::size_t j=i+1;j<g;++j){r.hermiticity_error=std::max(r.hermiticity_error,std::abs(r.operator_matrix[i*g+j]-r.operator_matrix[j*g+i]));r.operator_matrix[i*g+j]=r.operator_matrix[j*g+i]=.5*(r.operator_matrix[i*g+j]+r.operator_matrix[j*g+i]);}
    const auto operator_eigen=diagonalize(r.operator_matrix,g);r.eigensolver_residual=operator_eigen.residual;
    if(!operator_eigen.converged)return fail(OrbitalGroupBondingStatus::Unavailable,"operator_eigensolver_not_converged");
    r.eigenvalues=operator_eigen.values;std::sort(r.eigenvalues.begin(),r.eigenvalues.end());
    r.minimum_eigenvalue=r.eigenvalues.front();r.maximum_eigenvalue=r.eigenvalues.back();
    r.trace=std::accumulate(r.eigenvalues.begin(),r.eigenvalues.end(),0.);
    r.error_bound=options.absolute_error+options.relative_error*maximum(r.eigenvalues)+r.eigensolver_residual+r.hermiticity_error;
    Matrix occupation=supplied_occupation;
    r.occupation_available=!occupation.empty();
    if(occupation.empty()){
        occupation.assign(g*g,0);r.occupation_available=true;
        for(std::size_t k=0;k<g;++k){const auto& mo=w.orbitals[members[k]];occupation[k*g+k]=mo.occupation;
            r.occupation_available&=mo.occupation_provenance!=DataProvenance::Unavailable && std::isfinite(mo.occupation);}
    }
    if(r.occupation_available){
        if(occupation.size()!=g*g || !finite(occupation))return fail(OrbitalGroupBondingStatus::Unavailable,"complete_finite_occupation_matrix_unavailable");
        for(std::size_t i=0;i<g;++i)for(std::size_t j=i+1;j<g;++j)
            if(std::abs(occupation[i*g+j]-occupation[j*g+i])>1e-8)return fail(OrbitalGroupBondingStatus::Unavailable,"occupation_matrix_not_hermitian");
        const auto occ_eigen=diagonalize(occupation,g);
        if(!occ_eigen.converged || *std::min_element(occ_eigen.values.begin(),occ_eigen.values.end())<-1e-8)
            return fail(OrbitalGroupBondingStatus::Unavailable,"occupation_matrix_not_positive_semidefinite");
        // C = (C G^-1/2) G^1/2, so the physical occupation matrix in the
        // normalized copy is G^1/2 f G^1/2, not an unchanged diagonal.
        r.occupation_matrix=congruence(occupation,metric_root,g);
        for(std::size_t i=0;i<g;++i){r.electron_count+=r.occupation_matrix[i*g+i];for(std::size_t j=0;j<g;++j)r.electron_weighted_trace+=r.occupation_matrix[i*g+j]*r.operator_matrix[j*g+i];}
    }
    if(r.minimum_eigenvalue>r.error_bound)r.status=OrbitalGroupBondingStatus::Positive;
    else if(r.maximum_eigenvalue<-r.error_bound)r.status=OrbitalGroupBondingStatus::Negative;
    else if(r.minimum_eigenvalue<-r.error_bound && r.maximum_eigenvalue>r.error_bound)r.status=OrbitalGroupBondingStatus::Mixed;
    else r.status=OrbitalGroupBondingStatus::Unresolved;
    r.detail="fixed_scope_hermitian_overlap_spectrum_not_bond_energy_or_all_local_pair_labels";return r;
}
const char* orbital_group_bonding_status_name(OrbitalGroupBondingStatus s) noexcept{
    switch(s){case OrbitalGroupBondingStatus::Positive:return "positive";case OrbitalGroupBondingStatus::Negative:return "negative";case OrbitalGroupBondingStatus::Mixed:return "mixed";case OrbitalGroupBondingStatus::Unresolved:return "unresolved";case OrbitalGroupBondingStatus::NotApplicable:return "not_applicable";case OrbitalGroupBondingStatus::IncompleteGroup:return "incomplete_group";default:return "unavailable";}
}
} // namespace cov
