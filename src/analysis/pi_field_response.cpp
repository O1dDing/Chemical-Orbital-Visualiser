#include "cov/pi_field_response.hpp"
#include "cov/nbo_salc.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace cov { namespace {
using M=Eigen::MatrixXd;using V=Eigen::VectorXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
// Projection conditioning and population resolution have their own units;
// neither can be compared directly with a Fock error in hartree.
constexpr double minimum_squared_d_projection=1e-10; // dimensionless
constexpr double occupation_resolution_electrons=5e-5; // density verification precision
bool valid(const NboMatrix& a){return a.rows&&a.columns&&a.values.size()==a.rows*a.columns&&
 std::all_of(a.values.begin(),a.values.end(),[](double x){return std::isfinite(x);});}
M mat(const NboMatrix& a){return Eigen::Map<const RM>(a.values.data(),a.rows,a.columns);}
NboMatrix matrix(const M& a,const char* kind){NboMatrix b;b.kind=kind;b.rows=a.rows();b.columns=a.cols();RM r=a;b.values.assign(r.data(),r.data()+r.size());return b;}
double mx(const M& a){return a.size()?a.cwiseAbs().maxCoeff():0;}
M sym(const M& a){return (a+a.transpose())*.5;}
double response_norm(const M& a){if(!a.size())return 0;Eigen::SelfAdjointEigenSolver<M> e(sym(a));return e.info()==Eigen::Success?e.eigenvalues().cwiseAbs().maxCoeff():a.norm();}
M basis(const M& a,double tol){if(!a.cols())return M(a.rows(),0);Eigen::JacobiSVD<M> s(a,Eigen::ComputeThinU);std::size_t r=0;
 const double cut=std::max(tol,tol*s.singularValues()[0]);for(double x:s.singularValues())if(x>cut)++r;return s.matrixU().leftCols(r);}
M joined(const std::vector<M>& parts,std::size_t n,double tol){std::size_t k=0;for(const auto& p:parts)k+=p.cols();M raw(n,k);k=0;for(const auto& p:parts){raw.middleCols(k,p.cols())=p;k+=p.cols();}return basis(raw,tol);}
M cross(const M& f,const M& d,const M& l){if(!l.cols())return M::Zero(f.rows(),f.cols());const M v=d.transpose()*f*l;return d*v*l.transpose()+l*v.transpose()*d.transpose();}
std::vector<Eigen::Index> most_overlap(const M& previous,const M& vectors){V w=(previous.transpose()*vectors).colwise().squaredNorm();std::vector<Eigen::Index> ids(vectors.cols());std::iota(ids.begin(),ids.end(),0);
 std::partial_sort(ids.begin(),ids.begin()+previous.cols(),ids.end(),[&](auto a,auto b){return w[a]>w[b];});ids.resize(previous.cols());std::sort(ids.begin(),ids.end());return ids;}
std::vector<Eigen::Index> complete_d_reference(const M& d,const M& vectors){
 // Rank-revealing column pivoting represents every independent d direction.
 // Ranking individual d weights can select two spectral partners of the
 // same direction and omit another even when a complete d reference exists.
 // Left orthogonal changes of the d basis leave this selection invariant.
 Eigen::ColPivHouseholderQR<M> qr(d.transpose()*vectors);std::vector<Eigen::Index> ids(d.cols());
 for(Eigen::Index i=0;i<d.cols();++i)ids[i]=qr.colsPermutation().indices()[i];std::sort(ids.begin(),ids.end());return ids;
}
M columns(const M& c,const std::vector<Eigen::Index>& ids){M x(c.rows(),ids.size());for(std::size_t k=0;k<ids.size();++k)x.col(k)=c.col(ids[k]);return x;}
V values(const V& v,const std::vector<Eigen::Index>& ids){V x(ids.size());for(std::size_t k=0;k<ids.size();++k)x[k]=v[ids[k]];return x;}
M effective(const M& d,const M& c,const V& e,double& minimum,M* embedding=nullptr){Eigen::JacobiSVD<M> sv(d.transpose()*c,Eigen::ComputeFullU|Eigen::ComputeFullV);
 minimum=sv.singularValues().array().square().minCoeff();const M polar=sv.matrixU()*sv.matrixV().transpose();if(embedding)*embedding=c*polar.transpose();return sym(polar*e.asDiagonal()*polar.transpose());}
struct Path {PiFieldTracking record;M start,end,embedding;};
Path track(const M& f0,const M& k,const M& d,const PiFieldResponseOptions& o,double error){Path out;out.record.rank=d.cols();
 Eigen::SelfAdjointEigenSolver<M> initial(f0);if(initial.info()!=Eigen::Success){out.record.reason="reference-eigensolver-failed";return out;}
 auto ids=complete_d_reference(d,initial.eigenvectors());M start=columns(initial.eigenvectors(),ids);V e0=values(initial.eigenvalues(),ids);
 double dp=0;out.start=effective(d,start,e0,dp,&out.embedding);out.record.minimum_d_projection=dp;
 if(dp<=minimum_squared_d_projection){out.record.reason="reference-d-subspace-projection-singular";return out;}
 out.record.reference_energies_hartree.assign(e0.data(),e0.data()+e0.size());
 if(k.norm()<=error){out.end=out.start;out.record.final_energies_hartree=out.record.reference_energies_hartree;out.record.verified=true;out.record.reason="numerically-zero-removed-block";return out;}
 for(std::size_t steps=o.initial_steps;steps<=o.maximum_steps;steps*=2){M previous=start;double overlap=1,min_d=dp;bool ok=true;V en=e0;M hn=f0;
  for(std::size_t step=1;step<=steps;++step){hn=f0+(double(step)/steps)*k;Eigen::SelfAdjointEigenSolver<M> es(hn);
   if(es.info()!=Eigen::Success){out.record.reason="path-eigensolver-failed";return out;}
   ids=most_overlap(previous,es.eigenvectors());M next=columns(es.eigenvectors(),ids);Eigen::JacobiSVD<M> s(previous.transpose()*next);
   const double ov=s.singularValues().array().square().minCoeff();overlap=std::min(overlap,ov);
   if(ov<o.minimum_step_overlap){ok=false;break;}previous=std::move(next);en=values(es.eigenvalues(),ids);
  }
  out.record.steps=steps;out.record.minimum_step_overlap=overlap;
  if(!ok){if(steps>o.maximum_steps/2)break;continue;}
  double end_dp=0;out.end=effective(d,previous,en,end_dp,&out.embedding);min_d=std::min(min_d,end_dp);out.record.minimum_d_projection=min_d;
  if(end_dp<=minimum_squared_d_projection){out.record.reason="final-d-subspace-projection-singular";return out;}
  out.record.endpoint_eigen_residual_hartree=(hn*previous-previous*en.asDiagonal()).norm();
  out.record.final_energies_hartree.assign(en.data(),en.data()+en.size());out.record.verified=true;out.record.reason="complete-d-subspace-adiabatically-tracked";return out;
 }
 out.record.reason="d-subspace-branch-not-resolved-at-maximum-path-resolution";return out;
}
PiFieldRole classify(const M& delta,double gross,double gate){Eigen::SelfAdjointEigenSolver<M> es(sym(delta));if(es.info()!=Eigen::Success)return PiFieldRole::Unavailable;
 const double lo=es.eigenvalues().minCoeff(),hi=es.eigenvalues().maxCoeff();
 if(lo>=-gate&&hi<=gate)return gross>gate?PiFieldRole::Mixed:PiFieldRole::Negligible;
 if(lo< -gate&&hi>gate)return PiFieldRole::Mixed;
 return hi>gate?PiFieldRole::DonorDominant:PiFieldRole::AcceptorDominant;}
bool usual_mechanisms(const M& occupied,const M& vacant,double gate){Eigen::SelfAdjointEigenSolver<M> a(sym(occupied)),b(sym(vacant));return a.info()==Eigen::Success&&b.info()==Eigen::Success&&a.eigenvalues().minCoeff()>=-gate&&b.eigenvalues().maxCoeff()<=gate;}
const NboMatrix* unique(const std::vector<NboMatrix>& all,const char* name,NboSpin spin){const NboMatrix* found=nullptr;for(const auto& x:all)if(x.kind==name&&x.spin==spin){if(found)return nullptr;found=&x;}return found;}
bool square(const NboMatrix* a,std::size_t n){return a&&valid(*a)&&a->rows==n&&a->columns==n;}
void str(std::ostream& o,const std::string& x){o<<'"';for(char c:x){switch(c){case '"':o<<"\\\"";break;case '\\':o<<"\\\\";break;case '\n':o<<"\\n";break;case '\r':o<<"\\r";break;case '\t':o<<"\\t";break;default:if(static_cast<unsigned char>(c)>=32)o<<c;}}o<<'"';}
template<class T>void arr(std::ostream& o,const std::vector<T>& a){o<<'[';for(std::size_t i=0;i<a.size();++i){if(i)o<<',';o<<a[i];}o<<']';}
void jmatrix(std::ostream& o,const NboMatrix& a){o<<"{\"rows\":"<<a.rows<<",\"columns\":"<<a.columns<<",\"values\":";arr(o,a.values);o<<'}';}
}
const char* pi_field_role_name(PiFieldRole r) noexcept {switch(r){case PiFieldRole::Negligible:return "negligible";case PiFieldRole::DonorDominant:return "donor_dominant";case PiFieldRole::AcceptorDominant:return "acceptor_dominant";case PiFieldRole::Mixed:return "mixed";default:return "unavailable";}}

PiFieldResponse assess_pi_field_response(const PiFieldResponseInput& in,const PiFieldResponseOptions& o){PiFieldResponse out;
 out.canonical_fingerprint=in.canonical_fingerprint;out.operator_semantics=in.operator_semantics;out.scope_semantics=in.scope_semantics;out.centre_atoms=in.centre_atoms;out.ligand_atoms=in.ligand_atoms;out.canonical_indices=in.canonical_indices;out.shared_spatial_verified=in.shared_spatial_verified;
 out.response_resolution_hartree=o.response_resolution_hartree;out.density_verified=in.density_verified;
 out.minimum_group_spectral_fraction=o.minimum_group_spectral_fraction;
 const auto fail=[&](const char* reason){out.reason=reason;return out;};
 if(!valid(in.fock)||in.fock.rows!=in.fock.columns||!valid(in.metal_d_basis)||!valid(in.ligand_pi_basis)||
    in.metal_d_basis.rows!=in.fock.rows||in.ligand_pi_basis.rows!=in.fock.rows||!std::isfinite(in.operator_error_hartree)||in.operator_error_hartree<0||
    !std::isfinite(o.numerical_tolerance)||o.numerical_tolerance<=0||!std::isfinite(o.response_resolution_hartree)||o.response_resolution_hartree<0||
    !o.initial_steps||o.maximum_steps<o.initial_steps||o.minimum_step_overlap<=0||o.minimum_step_overlap>=1||
    !std::isfinite(o.minimum_group_spectral_fraction)||o.minimum_group_spectral_fraction<=0||o.minimum_group_spectral_fraction>1)return fail("invalid-complete-shared-field-input");
 const M f=mat(in.fock),d=mat(in.metal_d_basis),l=mat(in.ligand_pi_basis);const auto n=f.rows(),r=d.cols();
 if(r>=n||mx(f-f.transpose())>2e-5||mx(d.transpose()*d-M::Identity(r,r))>5e-5||mx(l.transpose()*l-M::Identity(l.cols(),l.cols()))>5e-5||mx(d.transpose()*l)>5e-5)return fail("shared-orthonormal-partition-or-hermiticity-failed");
 out.full_dimension=n;out.metal_rank=r;out.ligand_rank=l.cols();out.full_environment_retained=true;
 const double error=std::max(o.numerical_tolerance,double(n)*std::max(in.operator_error_hartree,mx(f-f.transpose())));
 out.numerical_error_bound_hartree=error;const M k=cross(f,d,l),f0=sym(f-k);const auto total=track(f0,k,d,o,error);out.tracking=total.record;
 if(!total.record.verified)return fail(total.record.reason.c_str());
 const M delta=sym(total.end-total.start);M occupied=M::Zero(r,r),vacant=M::Zero(r,r);
 M density;
 if(in.density_verified){if(!valid(in.total_density)||in.total_density.rows!=n||in.total_density.columns!=n)return fail("verified-density-missing");density=sym(mat(in.total_density));
  Eigen::SelfAdjointEigenSolver<M> ds(density);if(ds.info()!=Eigen::Success||ds.eigenvalues().minCoeff()< -5e-5||ds.eigenvalues().maxCoeff()>2.00005)return fail("invalid-spin-summed-density-spectrum");
  out.occupied_d_electrons=(d.transpose()*density*d).trace();
 }
 // The occupation partition only resolves competing field responses. It is
 // never itself a transfer direction, nor a reference fragment CT density.
 if(in.density_verified){Eigen::SelfAdjointEigenSolver<M> ls(sym(l.transpose()*density*l));if(ls.info()!=Eigen::Success)return fail("ligand-density-partition-failed");
  std::vector<Eigen::Index> oi,vi;for(Eigen::Index i=0;i<l.cols();++i)(ls.eigenvalues()[i]>=1?oi:vi).push_back(i);
  for(int which=0;which<2;++which){const auto& ids=which?vi:oi;if(ids.empty())continue;const M part=l*columns(ls.eigenvectors(),ids);const M kp=cross(f,d,part);
   const auto path=track(f0,kp,d,o,error);if(!path.record.verified)return fail("competing-field-partition-tracking-failed");
   if(which)vacant=sym(path.end-path.start);else occupied=sym(path.end-path.start);
  }
 }else {occupied=delta;out.scope_semantics+="; competing-occupation-partition-unavailable";}
 out.d_response=matrix(delta,"d-field-response");out.d_occupied_pi_response=matrix(occupied,"higher-occupation-pi-field-response");out.d_vacant_pi_response=matrix(vacant,"lower-occupation-pi-field-response");
 out.trace_shift_hartree=delta.trace();out.mean_shift_hartree=delta.trace()/r;
 out.gross_response_hartree=std::max(response_norm(occupied),response_norm(vacant));out.nonadditivity_norm_hartree=(delta-occupied-vacant).norm();
 out.role=classify(delta,out.gross_response_hartree,std::max(error,o.response_resolution_hartree));
 const bool mechanisms_ok=!in.density_verified||usual_mechanisms(occupied,vacant,std::max(error,o.response_resolution_hartree));
 if(!mechanisms_ok)out.role=PiFieldRole::Mixed;
 out.occupied_metal_backbond_supported=in.density_verified&&in.occupied_metal_backbond_evidence&&out.occupied_d_electrons>occupation_resolution_electrons;
 if(out.occupied_metal_backbond_supported&&valid(in.occupied_backbond_donor_basis)&&in.occupied_backbond_donor_basis.rows==n){const M b=basis(d.transpose()*mat(in.occupied_backbond_donor_basis),1e-6);out.d_backbond_projector=matrix(b*b.transpose(),"occupied-local-metal-backbond-support");}
 if(valid(in.canonical_columns)&&in.canonical_columns.rows==n&&in.canonical_columns.columns==in.canonical_indices.size()){
  out.canonical_d_coordinates=matrix(d.transpose()*mat(in.canonical_columns),"canonical-in-common-d-reference");
  out.canonical_tracked_d_coordinates=matrix(total.embedding.transpose()*mat(in.canonical_columns),"canonical-in-tracked-d-spectral-subspace");
 }
 out.available=true;out.reason=mechanisms_ok?"joined-field-response-in-full-orthonormal-environment":"joined-field-response-with-inverted-or-competing-occupation-family-mechanisms";return out;
}

PiFieldGroupAssessment assess_pi_field_group(const PiFieldResponse& response,const std::vector<std::size_t>& members){PiFieldGroupAssessment out;out.members=members;out.numerical_error_bound_hartree=response.numerical_error_bound_hartree;
 if(!response.available||members.empty()||!valid(response.canonical_d_coordinates)||!valid(response.canonical_tracked_d_coordinates)||!valid(response.d_response)){out.reason="joined-response-or-complete-group-unavailable";return out;}
 std::set<std::size_t> distinct;const M all=mat(response.canonical_d_coordinates),tracked=mat(response.canonical_tracked_d_coordinates);M c(all.rows(),members.size()),tc(all.rows(),members.size());
 for(std::size_t k=0;k<members.size();++k){if(!distinct.insert(members[k]).second){out.reason="duplicate-canonical-member";return out;}const auto it=std::find(response.canonical_indices.begin(),response.canonical_indices.end(),members[k]);
  if(it==response.canonical_indices.end()){out.reason="canonical-member-outside-shared-reference";return out;}c.col(k)=all.col(it-response.canonical_indices.begin());tc.col(k)=tracked.col(it-response.canonical_indices.begin());}
 out.d_fraction=c.squaredNorm()/members.size();const M q=basis(c,1e-6);out.metal_support_rank=q.cols();out.available=true;
 if(!q.cols()){out.reason="no-resolved-metal-d-support";out.role=PiFieldRole::Negligible;return out;}
 const M h=mat(response.d_response),ho=mat(response.d_occupied_pi_response),hv=mat(response.d_vacant_pi_response);
 const double gate=std::max(response.response_resolution_hartree,response.numerical_error_bound_hartree);
 // Responsive d directions include cancelling mechanisms, not just the net
 // matrix. This excludes sigma-only/eg support without deleting cancellation.
 Eigen::SelfAdjointEigenSolver<M> active(sym(h*h+ho*ho+hv*hv));M projector=M::Zero(h.rows(),h.rows());
 for(Eigen::Index i=0;i<h.rows();++i)if(active.eigenvalues()[i]>gate*gate)projector+=active.eigenvectors().col(i)*active.eigenvectors().col(i).transpose();
 out.responsive_d_fraction=(c.transpose()*projector*c).trace()/members.size();
 out.tracked_d_fraction=tc.squaredNorm()/members.size();out.tracked_responsive_fraction=(tc.transpose()*projector*tc).trace()/members.size();
 const M shift=sym(q.transpose()*h*q);Eigen::SelfAdjointEigenSolver<M> es(shift);out.trace_shift_hartree=shift.trace();out.mean_shift_hartree=shift.trace()/q.cols();out.response_min_hartree=es.eigenvalues().minCoeff();out.response_max_hartree=es.eigenvalues().maxCoeff();
 out.gross_response_hartree=std::max(response_norm(q.transpose()*ho*q),response_norm(q.transpose()*hv*q));out.role=classify(shift,out.gross_response_hartree,gate);
 const bool mechanisms_ok=!response.density_verified||usual_mechanisms(q.transpose()*ho*q,q.transpose()*hv*q,gate);if(!mechanisms_ok)out.role=PiFieldRole::Mixed;
 out.applicable=out.responsive_d_fraction>1e-6&&out.tracked_responsive_fraction>=response.minimum_group_spectral_fraction&&out.role!=PiFieldRole::Negligible;
 out.occupied_metal_backbond_supported=out.applicable&&response.occupied_metal_backbond_supported&&valid(response.d_backbond_projector)&&
  (q.transpose()*mat(response.d_backbond_projector)*q).trace()>1e-6;
 out.reason=out.applicable?(mechanisms_ok?"complete-group-corresponds-to-tracked-responsive-d-subspace":"group-has-inverted-or-competing-occupation-family-mechanisms"):
  out.responsive_d_fraction>1e-6?"group-is-not-a-major-component-of-tracked-responsive-d-levels":"group-has-no-resolved-pi-field-response";return out;
}

PiFieldResponseAnalysis analyse_pi_field_responses(const Wavefunction& w,const NboIntegration& data,std::span<const NboPiCoupling* const> channels,const PiFieldResponseOptions& o){PiFieldResponseAnalysis out;
 if(!data.dataset.association.compatible||data.canonical_fingerprint!=nbo_canonical_fingerprint(w)){out.status="rejected";out.reason="canonical-NBO-association-mismatch";return out;}
 const auto n=static_cast<std::size_t>(w.basis_count);if(!n||!data.dataset.archive){out.reason="same-source-archive-unavailable";return out;}
 const auto& arc=data.dataset.archive->matrices;const auto* s=unique(arc,"OVERLAP",NboSpin::Total);if(!square(s,n)){out.reason="AO-metric-unavailable";return out;}
 const auto ro=verify_nbo_restricted_open_shell(w,data);const bool shared_ro=ro.verified;
 const auto* ft=unique(arc,"FOCK",NboSpin::Total);const auto* fa=unique(arc,"FOCK",NboSpin::Alpha);const auto* fb=unique(arc,"FOCK",NboSpin::Beta);
 const bool closed=square(ft,n)&&!square(fa,n)&&!square(fb,n);
 if(!closed&&(!shared_ro||!square(fa,n)||!square(fb,n))){out.status="unsupported";out.reason="requires-verified-shared-RO-or-restricted-closed-shell-reference";return out;}
 const NboSpin primary=closed?NboSpin::Total:NboSpin::Alpha;
 const auto ao_for=[&](NboSpin spin){const auto* a=unique(data.dataset.matrices,"AONAO",spin);if(!a)a=unique(data.dataset.matrices,"AONAO",NboSpin::Total);return a;};
 const auto* a=ao_for(primary);const auto* u=unique(data.dataset.matrices,"NAOMO",primary);if(!square(a,n)||!square(u,n)){out.reason="complete-common-NAO-transform-unavailable";return out;}
 const M A=mat(*a),S=mat(*s),C=mat(*u);if(mx(A.transpose()*S*A-M::Identity(n,n))>5e-5||mx(C.transpose()*C-M::Identity(n,n))>5e-5){out.status="rejected";out.reason="common-reference-orthogonality-failed";return out;}
 // Reuse the independently associated physical-operator gates also used by
 // the SALC model. An available geometric channel, or a default zero error,
 // is not proof that either physical spin Fock is valid.
 std::optional<NboSalcModel> physical=build_nbo_salc_model(w,data);double physical_error=0;
 for(auto spin:closed?std::vector<NboSpin>{NboSpin::Total}:std::vector<NboSpin>{NboSpin::Alpha,NboSpin::Beta}){
  const auto* raw=spin==NboSpin::Total?ft:spin==NboSpin::Alpha?fa:fb;const auto* as=ao_for(spin);
  const auto ev=std::find_if(physical->energies.begin(),physical->energies.end(),[&](const auto& x){return x.spin==spin&&x.available&&(x.canonical_same_operator||x.printed_operator_verified);});
  const auto op=std::find_if(physical->spin_operators.begin(),physical->spin_operators.end(),[&](const auto& x){return x.spin==spin;});
  if(!square(raw,n)||!square(as,n)||ev==physical->energies.end()||op==physical->spin_operators.end()||op->fock.size()!=n*n||op->basis.size()!=n){out.status="rejected";out.reason="independently-verified-physical-spin-Fock-unavailable";return out;}
  const double antisym=mx(mat(*raw)-mat(*raw).transpose());bool same_basis=true;for(std::size_t i=0;i<n;++i)same_basis=same_basis&&op->basis[i].kind==NboOrbitalKind::NAO&&op->basis[i].spin==spin&&op->basis[i].index==i;
  const M verified=Eigen::Map<const RM>(op->fock.data(),n,n);const double mismatch=mx(mat(*as).transpose()*mat(*raw)*mat(*as)-verified);
  if(!same_basis||!std::isfinite(antisym)||!std::isfinite(mismatch)||antisym>2e-5||mismatch>2e-5){out.status="rejected";out.reason="raw-physical-Fock-hermiticity-or-reference-mismatch";return out;}
  physical_error=std::max({physical_error,antisym,mismatch,ev->hermiticity_error});
 }
 const M F=sym(A.transpose()*(closed?mat(*ft):M((mat(*fa)+mat(*fb))*.5))*A);
 // Source occupations must have already passed the archive density check.
 M density=M::Zero(n,n);bool density_ok=true;
 for(auto spin:closed?std::vector<NboSpin>{NboSpin::Total}:std::vector<NboSpin>{NboSpin::Alpha,NboSpin::Beta}){
  const auto val=std::find_if(data.dataset.nao_validation.begin(),data.dataset.nao_validation.end(),[&](const auto& v){return v.spin==spin&&v.available&&v.occupation_density_verified&&v.canonical_occupations.size()==n;});
  const auto* us=unique(data.dataset.matrices,"NAOMO",spin);const auto* as=ao_for(spin);
  if(!square(us,n)||!square(as,n)){density_ok=false;break;}
  const M transform=A.transpose()*S*mat(*as),cs=transform*mat(*us);
  if(val!=data.dataset.nao_validation.end()){const V occ=Eigen::Map<const V>(val->canonical_occupations.data(),n);density+=cs*occ.asDiagonal()*cs.transpose();}
  else if(shared_ro){
   if(!physical)physical=build_nbo_salc_model(w,data);
   const auto op=std::find_if(physical->spin_operators.begin(),physical->spin_operators.end(),[&](const auto& x){return x.spin==spin;});
   if(op==physical->spin_operators.end()||op->density.size()!=n*n||op->basis.size()!=n){density_ok=false;break;}
   bool basis_ok=true;for(std::size_t i=0;i<n;++i)basis_ok=basis_ok&&op->basis[i].kind==NboOrbitalKind::NAO&&op->basis[i].spin==spin&&op->basis[i].index==i;
   const auto electrons=spin==NboSpin::Beta?w.beta_electrons:w.alpha_electrons;V occ=V::Zero(n);if(electrons>n){density_ok=false;break;}occ.head(electrons).setOnes();
   const M ds=Eigen::Map<const RM>(op->density.data(),n,n),usm=mat(*us);
   if(!basis_ok||mx(usm*occ.asDiagonal()*usm.transpose()-ds)>2e-5){density_ok=false;break;}density+=transform*ds*transform.transpose();
  }else {density_ok=false;break;}
 }
 if(!density_ok){out.reason="verified-total-density-unavailable";return out;}
 // All original common-space canonical columns, with their immutable global
 // indices. No source energy ordering or display ordinal is used here.
 std::vector<std::size_t> global,source;for(std::size_t i=0;i<w.orbitals.size();++i)if((w.orbitals[i].spin!=Spin::Beta||shared_ro)&&w.orbitals[i].source_orbital_index<n){global.push_back(i);source.push_back(w.orbitals[i].source_orbital_index);}
 M beta_columns;if(shared_ro){const auto* ub=unique(data.dataset.matrices,"NAOMO",NboSpin::Beta);const auto* ab=ao_for(NboSpin::Beta);if(square(ub,n)&&square(ab,n))beta_columns=A.transpose()*S*mat(*ab)*mat(*ub);}
 M canonical(n,source.size());for(std::size_t k=0;k<source.size();++k){if(w.orbitals[global[k]].spin==Spin::Beta&&beta_columns.size())canonical.col(k)=beta_columns.col(source[k]);else canonical.col(k)=C.col(source[k]);}
 std::map<std::vector<std::size_t>,std::vector<const NboPiCoupling*>> spheres;
 for(const auto* ch:channels)if(ch&&ch->centre_family=="valence-d"&&!ch->centre_atoms.empty()&&valid(ch->centre_projector_basis)&&valid(ch->ligand_projector_basis)&&
    ch->centre_projector_basis.rows==n&&ch->ligand_projector_basis.rows==n&&ch->spin==primary)spheres[ch->centre_atoms].push_back(ch);
 for(const auto& [centres,parts]:spheres){std::vector<M> ligand_parts,metal_parts,backbond_parts;std::set<std::size_t> ligands;double error=physical_error;bool occupied_edge=false;
  for(const auto* ch:parts){const auto* as=ao_for(ch->spin);if(!square(as,n))continue;const M map=A.transpose()*S*mat(*as);
   metal_parts.push_back(map*mat(ch->centre_projector_basis));ligand_parts.push_back(map*mat(ch->ligand_projector_basis));ligands.insert(ch->ligand_atoms.begin(),ch->ligand_atoms.end());error=std::max(error,ch->operator_max_error_hartree);
  }
  // Backbond qualification may use either real occupied spin, but never
  // affects the shared field calculation or its sign.
  for(const auto* ch:channels)if(ch&&ch->centre_atoms==centres&&ch->centre_family=="valence-d")for(const auto& e:ch->direction_projection_evidence)
   if(e.direction=="centre_to_ligand"&&e.donor_occupation>0.5&&e.acceptor_occupation<0.5&&e.donor_centre_weight>0.5&&e.acceptor_ligand_weight>0.5&&std::isfinite(e.fock_hartree)&&std::abs(e.fock_hartree)>1e-7){
    const auto* nb=unique(data.dataset.matrices,"NAONBO",ch->spin);const auto* as=ao_for(ch->spin);if(square(nb,n)&&square(as,n)&&e.donor_nbo_id&&e.donor_nbo_id<=n){backbond_parts.push_back(A.transpose()*S*mat(*as)*mat(*nb).col(e.donor_nbo_id-1));occupied_edge=true;}
   }
  const M d=joined(metal_parts,n,1e-7);M l=joined(ligand_parts,n,1e-7);if(!d.cols()||!l.cols())continue;l=basis(l-d*(d.transpose()*l),1e-7);
  PiFieldResponseInput in;in.fock=matrix(F,"shared-spin-scalar-Fock-NAO");in.metal_d_basis=matrix(d,"complete-valence-d");in.ligand_pi_basis=matrix(l,"union-complete-geometric-pi-projectors");in.total_density=matrix(density,"spin-summed-density");in.canonical_columns=matrix(canonical,"common-spatial-canonical");in.canonical_indices=global;in.centre_atoms=centres;in.ligand_atoms.assign(ligands.begin(),ligands.end());in.canonical_fingerprint=data.canonical_fingerprint;
  in.shared_spatial_verified=shared_ro||closed;in.density_verified=true;in.operator_error_hartree=error;in.occupied_metal_backbond_evidence=occupied_edge;
  if(!backbond_parts.empty())in.occupied_backbond_donor_basis=matrix(joined(backbond_parts,n,1e-7),"occupied-NBO-backbond-donor-span");
  in.operator_semantics=closed?"restricted-closed-shell-physical-Fock":"shared-RO-(Falpha+Fbeta)/2-frozen-one-electron-operator";
  in.scope_semantics="full-NAO-environment; union of complete valence-d geometric transverse-p and internal ligand pi/pi-star projectors; no direction votes";
  auto response=assess_pi_field_response(in,o);response.id="joined-pi-field";for(auto c:centres)response.id+=":"+std::to_string(c);
  // Existing complete source groups supply membership only, never direction.
  std::set<std::vector<std::size_t>> groups;for(const auto* ch:parts)for(const auto& g:ch->groups)groups.insert(g.members);
  for(const auto& group:groups)response.groups.push_back(assess_pi_field_group(response,group));out.responses.push_back(std::move(response));
 }
 if(out.responses.empty()){out.reason="complete-geometric-d-pi-scope-unavailable";return out;}
 const bool any=std::any_of(out.responses.begin(),out.responses.end(),[](const auto& r){return r.available;});out.status=any?"available":"insufficient_evidence";out.reason=any?"joined-field-responses-computed-in-shared-full-reference":"all-shared-field-responses-failed-with-explicit-diagnostics";return out;
}

std::string pi_field_group_assessment_json(const PiFieldGroupAssessment& x){std::ostringstream o;o<<std::setprecision(17)<<"{\"available\":"<<(x.available?"true":"false")<<",\"applicable\":"<<(x.applicable?"true":"false")<<",\"role\":";str(o,pi_field_role_name(x.role));o<<",\"reason\":";str(o,x.reason);o<<",\"method\":";str(o,x.method);o<<",\"members\":";arr(o,x.members);
 o<<",\"tracked_d_fraction\":"<<x.tracked_d_fraction<<",\"tracked_responsive_fraction\":"<<x.tracked_responsive_fraction;
 o<<",\"metal_support_rank\":"<<x.metal_support_rank<<",\"d_fraction\":"<<x.d_fraction<<",\"responsive_d_fraction\":"<<x.responsive_d_fraction<<",\"mean_shift_hartree\":"<<x.mean_shift_hartree<<",\"trace_shift_hartree\":"<<x.trace_shift_hartree<<",\"response_min_hartree\":"<<x.response_min_hartree<<",\"response_max_hartree\":"<<x.response_max_hartree<<",\"gross_response_hartree\":"<<x.gross_response_hartree<<",\"numerical_error_bound_hartree\":"<<x.numerical_error_bound_hartree<<",\"occupied_metal_backbond_supported\":"<<(x.occupied_metal_backbond_supported?"true":"false")<<'}';return o.str();}
std::string pi_field_response_json(const PiFieldResponse& x,bool matrices){std::ostringstream o;o<<std::setprecision(17)<<"{\"id\":";str(o,x.id);o<<",\"available\":"<<(x.available?"true":"false")<<",\"role\":";str(o,pi_field_role_name(x.role));o<<",\"reason\":";str(o,x.reason);o<<",\"canonical_fingerprint\":";str(o,x.canonical_fingerprint);o<<",\"operator_semantics\":";str(o,x.operator_semantics);o<<",\"scope_semantics\":";str(o,x.scope_semantics);o<<",\"centre_atoms\":";arr(o,x.centre_atoms);o<<",\"ligand_atoms\":";arr(o,x.ligand_atoms);
 o<<",\"minimum_group_spectral_fraction\":"<<x.minimum_group_spectral_fraction;
 o<<",\"full_environment_retained\":"<<(x.full_environment_retained?"true":"false")<<",\"shared_spatial_verified\":"<<(x.shared_spatial_verified?"true":"false")<<",\"full_dimension\":"<<x.full_dimension<<",\"metal_rank\":"<<x.metal_rank<<",\"ligand_rank\":"<<x.ligand_rank<<",\"mean_shift_hartree\":"<<x.mean_shift_hartree<<",\"trace_shift_hartree\":"<<x.trace_shift_hartree<<",\"gross_response_hartree\":"<<x.gross_response_hartree<<",\"nonadditivity_norm_hartree\":"<<x.nonadditivity_norm_hartree<<",\"numerical_error_bound_hartree\":"<<x.numerical_error_bound_hartree<<",\"response_resolution_hartree\":"<<x.response_resolution_hartree<<",\"occupied_d_electrons\":"<<x.occupied_d_electrons<<",\"density_verified\":"<<(x.density_verified?"true":"false")<<",\"occupied_metal_backbond_supported\":"<<(x.occupied_metal_backbond_supported?"true":"false");
 o<<",\"tracking\":{\"verified\":"<<(x.tracking.verified?"true":"false")<<",\"reason\":";str(o,x.tracking.reason);o<<",\"rank\":"<<x.tracking.rank<<",\"steps\":"<<x.tracking.steps<<",\"minimum_step_overlap\":"<<x.tracking.minimum_step_overlap<<",\"minimum_d_projection\":"<<x.tracking.minimum_d_projection<<",\"endpoint_eigen_residual_hartree\":"<<x.tracking.endpoint_eigen_residual_hartree<<",\"reference_energies_hartree\":";arr(o,x.tracking.reference_energies_hartree);o<<",\"final_energies_hartree\":";arr(o,x.tracking.final_energies_hartree);o<<'}';
 if(matrices){o<<",\"d_response\":";jmatrix(o,x.d_response);o<<",\"d_occupied_pi_response\":";jmatrix(o,x.d_occupied_pi_response);o<<",\"d_vacant_pi_response\":";jmatrix(o,x.d_vacant_pi_response);o<<",\"canonical_d_coordinates\":";jmatrix(o,x.canonical_d_coordinates);o<<",\"canonical_indices\":";arr(o,x.canonical_indices);}
 if(matrices){o<<",\"canonical_tracked_d_coordinates\":";jmatrix(o,x.canonical_tracked_d_coordinates);o<<",\"d_backbond_projector\":";jmatrix(o,x.d_backbond_projector);}
 o<<",\"groups\":[";for(std::size_t i=0;i<x.groups.size();++i){if(i)o<<',';o<<pi_field_group_assessment_json(x.groups[i]);}o<<"]}";return o.str();}
std::string pi_field_response_analysis_json(const PiFieldResponseAnalysis& x,bool matrices){std::ostringstream o;o<<"{\"status\":";str(o,x.status);o<<",\"reason\":";str(o,x.reason);o<<",\"responses\":[";for(std::size_t i=0;i<x.responses.size();++i){if(i)o<<',';o<<pi_field_response_json(x.responses[i],matrices);}o<<"]}";return o.str();}
}
