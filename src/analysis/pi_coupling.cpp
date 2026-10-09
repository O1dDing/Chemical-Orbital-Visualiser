#include "cov/pi_coupling.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/local_angular_projection.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <span>

namespace cov { namespace {
using M=Eigen::MatrixXd;
using V=Eigen::VectorXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using V3=Eigen::Vector3d;
M mat(const NboMatrix& x){return Eigen::Map<const RM>(x.values.data(),x.rows,x.columns);}
double maxabs(const M& x){return x.size()?x.cwiseAbs().maxCoeff():0;}
const NboMatrix* unique_matrix(const std::vector<NboMatrix>& rows,const char* kind,NboSpin spin){
    const NboMatrix* out=nullptr;
    for(const auto& row:rows)if(row.kind==kind&&row.spin==spin){if(out)return nullptr;out=&row;}
    return out;
}
bool full(const NboMatrix* x,std::size_t n){
    return x&&x->rows==n&&x->columns==n&&x->values.size()==n*n&&
        std::all_of(x->values.begin(),x->values.end(),[](double v){return std::isfinite(v);});
}
bool valence(const NboNao& x,char angular){
    return x.type.find("Val")!=std::string::npos &&
        x.type.find(angular)!=std::string::npos;
}
int axis_of(const std::string& angular){
    if(angular=="px")return 0;if(angular=="py")return 1;
    if(angular=="pz")return 2;return -1;
}
struct Shells {
    std::map<std::size_t,std::array<std::size_t,3>> p;
    std::map<std::size_t,std::vector<std::size_t>> d;
};
Shells shells(const NboIntegration& integration,NboSpin spin){
    Shells out;std::map<std::size_t,std::array<bool,3>> seen;
    for(const auto& row:integration.dataset.naos){
        if(row.spin!=spin||!row.atom||!row.id)continue;
        const auto atom=row.atom-1;
        if(valence(row,'p')){
            const int axis=axis_of(row.angular);
            if(axis>=0){
                if(seen[atom][axis])out.p.erase(atom);
                else {out.p[atom][axis]=row.id;seen[atom][axis]=true;}
            }
        }
        if(valence(row,'d') && row.angular.size() && row.angular[0]=='d')
            out.d[atom].push_back(row.id);
    }
    for(auto it=out.p.begin();it!=out.p.end();){
        const auto& flags=seen[it->first];
        if(!(flags[0]&&flags[1]&&flags[2]))it=out.p.erase(it);
        else ++it;
    }
    for(auto it=out.d.begin();it!=out.d.end();){
        if(it->second.size()!=5)it=out.d.erase(it);else ++it;
    }
    return out;
}
M d_columns(std::size_t n,const std::vector<std::size_t>& ids){
    M q=M::Zero(n,ids.size());
    for(std::size_t i=0;i<ids.size();++i)q(ids[i]-1,i)=1;
    return q;
}
V3 bond_axis(const Wavefunction& w,std::size_t a,std::size_t b){
    V3 v(w.atoms[b].x-w.atoms[a].x,w.atoms[b].y-w.atoms[a].y,
         w.atoms[b].z-w.atoms[a].z);
    return v.norm()>1e-9?v.normalized():V3::Zero();
}
std::array<double,9> axis_frame(const V3& z){
    V3 guide=std::abs(z.z())<0.85?V3::UnitZ():V3::UnitY();
    const V3 x=guide.cross(z).normalized(),y=z.cross(x);
    return {x.x(),y.x(),z.x(),x.y(),y.y(),z.y(),x.z(),y.z(),z.z()};
}
struct PiFrame{M q;double leakage=0,centre=0,minimum_p=0,residual=0;
    std::array<double,3> spectrum{};};
std::optional<PiFrame> derived_pi_frame(const Wavefunction& w,const NboIntegration& data,
                        std::size_t atom,const std::array<std::size_t,3>& ids,
                        NboSpin spin,const V3& axis,std::string& reason){
    const auto fail=[&](const char* why)->std::optional<PiFrame>{reason=why;return std::nullopt;};
    if(axis.norm()<0.5)return fail("Bond axis cannot be defined");
    const auto n=static_cast<std::size_t>(w.basis_count);
    std::array<const NboOrbitalDescriptor*,3> descriptors{};
    for(int component=0;component<3;++component){
        descriptors[component]=nbo_orbital(data,{NboOrbitalKind::NAO,spin,ids[component]-1});
        if(!descriptors[component]&&spin!=NboSpin::Total)
            descriptors[component]=nbo_orbital(data,{NboOrbitalKind::NAO,
                NboSpin::Total,ids[component]-1});
        if(!descriptors[component]||descriptors[component]->coefficients.size()!=n)
            return fail("NAO descriptor or complete coefficient column is missing");
    }
    Wavefunction view=w;view.orbitals.clear();
    const std::array<std::array<double,3>,6> samples{{
        {1,0,0},{0,1,0},{0,0,1},
        {std::numbers::sqrt2/2,std::numbers::sqrt2/2,0},
        {std::numbers::sqrt2/2,0,std::numbers::sqrt2/2},
        {0,std::numbers::sqrt2/2,std::numbers::sqrt2/2}}};
    for(const auto& sample:samples){
        MolecularOrbital mo;mo.spin=spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        mo.occupation=1;mo.coefficients.assign(n,0);
        for(int component=0;component<3;++component)
            for(std::size_t k=0;k<n;++k)
                mo.coefficients[k]+=sample[component]*
                    descriptors[component]->coefficients[k];
        view.orbitals.push_back(std::move(mo));
    }
    try{
        LocalAngularProjectionWorkspace workspace(view,atom,axis_frame(axis));
        std::array<double,6> sigma{},pi{};
        double centre=1,residual=0;
        const M metric=Eigen::Map<const RM>(w.ao_overlap.data(),n,n);
        for(std::size_t sample=0;sample<6;++sample){
            const auto projection=workspace.project(std::span<const std::size_t>(&sample,1));
            if(projection.status!=MetricSubspaceStatus::Available||
               projection.represented_spin_orbital_rank!=1||
               !std::isfinite(projection.angular_partition_residual)||
               std::abs(projection.angular_partition_residual)>2e-5)
                return fail("Local angular projection rank or partition residual failed");
            const V coeff=Eigen::Map<const V>(view.orbitals[sample].coefficients.data(),n);
            const double norm=coeff.dot(metric*coeff);
            if(!std::isfinite(norm)||norm<0.7||std::abs(norm-1)>5e-5)
                return fail("NAO metric normalization failed");
            sigma[sample]=norm*projection.component_projection_traces[1];
            pi[sample]=norm*(projection.component_projection_traces[2]+
                       projection.component_projection_traces[3]);
            centre=std::min(centre,projection.centre_mean_fraction);
            residual=std::max(residual,std::abs(projection.angular_partition_residual));
        }
        M gs=M::Zero(3,3),gp=M::Zero(3,3);
        for(int i=0;i<3;++i){gs(i,i)=sigma[i];gp(i,i)=sigma[i]+pi[i];}
        const std::array<std::pair<int,int>,3> pairs{{{0,1},{0,2},{1,2}}};
        for(int pair=0;pair<3;++pair){
            const auto [i,j]=pairs[pair];
            const double s=sigma[3+pair]-(sigma[i]+sigma[j])*0.5;
            const double p=sigma[3+pair]+pi[3+pair]-(gp(i,i)+gp(j,j))*0.5;
            gs(i,j)=gs(j,i)=s;gp(i,j)=gp(j,i)=p;
        }
        Eigen::SelfAdjointEigenSolver<M> ps(gp);
        if(ps.info()!=Eigen::Success||ps.eigenvalues()[0]<0.7||centre<0.7)
            return fail("Local p metric or centre projection is below the supported threshold");
        const M gpi=gp-gs;
        Eigen::GeneralizedSelfAdjointEigenSolver<M> angular(gpi,gp);
        if(angular.info()!=Eigen::Success||!angular.eigenvalues().allFinite()||
           angular.eigenvalues()[0]< -2e-4||angular.eigenvalues()[2]>1+2e-4||
           std::abs(angular.eigenvalues()[0])>2e-4||
           std::abs(1-angular.eigenvalues()[1])>2e-4||
           std::abs(1-angular.eigenvalues()[2])>2e-4)
            return fail("Generalized pi projector eigenvalues failed the sigma/pi separation gate");
        const M eig=angular.eigenvectors().rightCols(2);
        Eigen::HouseholderQR<M> qr(eig);
        const M local=qr.householderQ()*M::Identity(3,2);
        const M local_gp=local.transpose()*gp*local;
        const M local_gs=local.transpose()*gs*local;
        Eigen::GeneralizedSelfAdjointEigenSolver<M> verify(local_gs,local_gp);
        if(verify.info()!=Eigen::Success||!verify.eigenvalues().allFinite()||
           verify.eigenvalues().maxCoeff()>2e-4)return fail("Residual sigma leakage exceeds tolerance");
        PiFrame result;result.q=M::Zero(n,2);
        for(int row=0;row<3;++row)result.q.row(ids[row]-1)=local.row(row);
        result.leakage=std::max(std::abs(angular.eigenvalues()[0]),
            std::max(std::abs(1-angular.eigenvalues()[1]),
                     std::abs(1-angular.eigenvalues()[2])));
        result.centre=centre;result.minimum_p=ps.eigenvalues()[0];
        result.residual=residual;
        for(int i=0;i<3;++i)result.spectrum[i]=angular.eigenvalues()[i];
        return result;
    }catch(const std::exception& e){reason=std::string("Local angular projector unavailable: ")+e.what();return std::nullopt;}
}
std::string verified_symmetry(const Wavefunction& w,std::size_t index){
    for(const auto& assignment:w.derived_orbital_symmetry_assignments)
        if(std::find(assignment.orbital_indices.begin(),assignment.orbital_indices.end(),index)!=
           assignment.orbital_indices.end() &&
           std::isfinite(assignment.subspace_retention) && assignment.subspace_retention>0.999)
            return assignment.label;
    return w.orbitals[index].symmetry_provenance!=DataProvenance::Unavailable?
        w.orbitals[index].symmetry:"";
}
std::vector<std::vector<std::size_t>> canonical_groups(const Wavefunction& w,NboSpin spin){
    std::vector<std::size_t> indices;
    for(std::size_t i=0;i<w.orbitals.size();++i)
        if((spin==NboSpin::Beta)==(w.orbitals[i].spin==Spin::Beta))indices.push_back(i);
    std::sort(indices.begin(),indices.end(),[&](auto a,auto b){
        if(w.orbitals[a].energy_hartree!=w.orbitals[b].energy_hartree)
            return w.orbitals[a].energy_hartree<w.orbitals[b].energy_hartree;
        return a<b;
    });
    std::vector<std::vector<std::size_t>> groups;
    for(auto index:indices){
        if(groups.empty()||
           std::abs(w.orbitals[index].energy_hartree-
                    w.orbitals[groups.back().front()].energy_hartree)>1e-6 ||
           verified_symmetry(w,index)!=verified_symmetry(w,groups.back().front()))
            groups.push_back({});
        groups.back().push_back(index);
    }
    return groups;
}
struct E2Convention {bool verified=false,nominal=false;};
E2Convention e2_convention(const NboIntegration& data,const M& fnbo,NboSpin spin){
    E2Convention out;std::size_t checked=0;bool nominal=true,actual=true;
    for(const auto& row:data.dataset.e2){
        if(row.spin!=spin||row.value<0.05||row.units!="kcal/mol"||!row.donor||!row.acceptor||
           row.donor>static_cast<std::size_t>(fnbo.rows())||row.acceptor>static_cast<std::size_t>(fnbo.rows()))continue;
        const auto donor=std::find_if(data.dataset.orbitals.begin(),data.dataset.orbitals.end(),
            [&](const auto& x){return x.spin==spin&&x.id==row.donor;});
        if(donor==data.dataset.orbitals.end())continue;
        const double gap=fnbo(row.acceptor-1,row.acceptor-1)-fnbo(row.donor-1,row.donor-1);
        if(gap<=1e-5)continue;
        const double base=std::pow(fnbo(row.donor-1,row.acceptor-1),2)/gap*627.5094740631;
        nominal=nominal&&std::abs((spin==NboSpin::Total?2.0:1.0)*base-row.value)<=0.0051;
        actual=actual&&std::abs(donor->occupation*base-row.value)<=0.0051;++checked;
    }
    out.verified=checked>0&&(nominal||actual);out.nominal=nominal&&!actual;
    return out;
}
struct E2Support{double donor=0,acceptor=0;bool verified=false;std::vector<NboPiDirectionProjection> rows;};
E2Support directed_e2(const NboIntegration& data,const M& transform,const M& fnbo,
                 const M& qc,const M& ql,const std::vector<std::size_t>& centre,
                 const std::vector<std::size_t>& ligand,NboSpin spin,bool ligand_to_centre,
                 double operator_error,const E2Convention& convention){
    E2Support support;std::set<std::size_t> donor_ids,acceptor_ids;
    const auto& from=ligand_to_centre?ligand:centre;
    const auto& to=ligand_to_centre?centre:ligand;
    const double full_occupation=spin==NboSpin::Total?2.0:1.0;
    const auto on_fragment=[&](const NboOrbitalDescriptor* d,const auto& atoms){
        if(!d||d->atoms.empty())return false;
        // A localized acceptor may be a polarized M-L BD* rather than a
        // one-centre LV. Atom labels locate its source; the following full
        // projector norm establishes which side of this pi scope it occupies.
        const auto contains=[](const auto& list,auto atom){return std::find(list.begin(),list.end(),atom)!=list.end();};
        return std::any_of(d->atoms.begin(),d->atoms.end(),[&](auto atom){return contains(atoms,atom);})&&
            std::all_of(d->atoms.begin(),d->atoms.end(),[&](auto atom){return contains(centre,atom)||contains(ligand,atom);});
    };
    // Search every localized occupied/virtual pair, including unprinted E2 rows.
    for(const auto& donor:data.dataset.orbitals){
        if(donor.spin!=spin||!donor.id||donor.id>static_cast<std::size_t>(transform.cols())||
           (donor.kind!="LP"&&donor.kind!="BD"&&donor.kind!="3C")||
           !std::isfinite(donor.occupation)||donor.occupation<0.5*full_occupation)continue;
        const auto* dd=nbo_orbital(data,{NboOrbitalKind::NBO,spin,donor.id-1});
        if(!on_fragment(dd,from))continue;
        const V bd=transform.col(donor.id-1);
        const double dw=(ligand_to_centre?ql.transpose()*bd:qc.transpose()*bd).squaredNorm();
        if(dw<0.5)continue;
        for(const auto& acceptor:data.dataset.orbitals){
            if(acceptor.spin!=spin||!acceptor.id||acceptor.id>static_cast<std::size_t>(transform.cols())||
               (acceptor.kind.find('*')==std::string::npos&&acceptor.kind!="LV"&&acceptor.kind!="RY")||
               !std::isfinite(acceptor.occupation)||acceptor.occupation>0.5*full_occupation||
               donor.occupation-acceptor.occupation<0.25*full_occupation)continue;
            const auto* ad=nbo_orbital(data,{NboOrbitalKind::NBO,spin,acceptor.id-1});
            if(!on_fragment(ad,to))continue;
            const V ba=transform.col(acceptor.id-1);
            const double aw=(ligand_to_centre?qc.transpose()*ba:ql.transpose()*ba).squaredNorm();
            if(aw<0.5)continue;
            NboPiDirectionProjection evidence;
            evidence.direction=ligand_to_centre?"ligand_to_centre":"centre_to_ligand";
            evidence.donor_nbo_id=donor.id;evidence.acceptor_nbo_id=acceptor.id;
            evidence.donor_ligand_weight=(ql.transpose()*bd).squaredNorm();
            evidence.donor_centre_weight=(qc.transpose()*bd).squaredNorm();
            evidence.acceptor_ligand_weight=(ql.transpose()*ba).squaredNorm();
            evidence.acceptor_centre_weight=(qc.transpose()*ba).squaredNorm();
            evidence.donor_occupation=donor.occupation;evidence.acceptor_occupation=acceptor.occupation;
            evidence.energy_gap_hartree=fnbo(acceptor.id-1,acceptor.id-1)-fnbo(donor.id-1,donor.id-1);
            evidence.fock_hartree=fnbo(donor.id-1,acceptor.id-1);
            evidence.method="ordered-NBO-from-verified-same-spin-Fock";evidence.source=donor.source;
            evidence.e2_id="recovered:"+std::string(nbo_spin_name(spin))+":"+std::to_string(donor.id)+":"+std::to_string(acceptor.id);
            evidence.recovered_from_unprinted_fock=true;
            for(std::size_t i=0;i<data.dataset.e2.size();++i){const auto& printed=data.dataset.e2[i];
                if(printed.spin==spin&&printed.donor==donor.id&&printed.acceptor==acceptor.id){
                    evidence.e2_id="e2:"+std::to_string(i);evidence.source=printed.source;
                    evidence.recovered_from_unprinted_fock=false;break;
                }
            }
            const double noise=std::max(1e-7,10*transform.rows()*operator_error);
            if(std::abs(evidence.fock_hartree)<=noise)continue;
            if(evidence.energy_gap_hartree<=noise){
                evidence.status="nonperturbative";evidence.reason="nonpositive-or-unresolved-NBO-denominator";
            }else{
                evidence.actual_occupation_second_order_hartree=donor.occupation*evidence.fock_hartree*evidence.fock_hartree/evidence.energy_gap_hartree;
                evidence.perturbation_occupation=convention.verified&&convention.nominal?full_occupation:donor.occupation;
                evidence.e2_hartree=evidence.perturbation_occupation*evidence.fock_hartree*evidence.fock_hartree/evidence.energy_gap_hartree;
                evidence.e2_reference=convention.verified?(convention.nominal?
                    "producer-validated-nominal-occupation-NBO-E2":"producer-validated-actual-occupation-NBO-E2"):
                    "actual-occupation-second-order-estimate-producer-normalization-unverified";
                evidence.status=std::abs(evidence.fock_hartree)/evidence.energy_gap_hartree>0.25?"nonperturbative":"available";
                evidence.reason=evidence.status=="available"?"occupied-to-acceptor-pi-support":"large-Fock-to-gap-ratio-use-block-strength";
                donor_ids.insert(donor.id-1);acceptor_ids.insert(acceptor.id-1);support.verified=true;
            }
            support.rows.push_back(std::move(evidence));
        }
    }
    for(auto id:donor_ids){const V column=transform.col(id);
        support.donor+=(ligand_to_centre?ql.transpose()*column:qc.transpose()*column).squaredNorm();}
    for(auto id:acceptor_ids){const V column=transform.col(id);
        support.acceptor+=(ligand_to_centre?qc.transpose()*column:ql.transpose()*column).squaredNorm();}
    support.donor/=std::max<Eigen::Index>(1,ligand_to_centre?ql.cols():qc.cols());
    support.acceptor/=std::max<Eigen::Index>(1,ligand_to_centre?qc.cols():ql.cols());
    return support;
}
void append_group(NboPiCoupling& out,const Wavefunction& w,const M& u,const M& f,
                  const M& qc,const M& ql,const M& pc,const M& pl,
                  const V* occupations,const std::vector<std::size_t>& members){
    M c(u.rows(),members.size());double occ=0;
    for(std::size_t j=0;j<members.size();++j){
        const auto& mo=w.orbitals[members[j]];
        if(mo.source_orbital_index>=static_cast<std::size_t>(u.cols()))return;
        c.col(j)=u.col(mo.source_orbital_index);
        if(occupations)occ+=(*occupations)[mo.source_orbital_index];
    }
    NboPiCanonicalGroup group;
    group.members=members;group.spin=out.spin;
    if(occupations)group.occupation_per_mo=occ/members.size();
    group.centre_weight=(qc.transpose()*c).squaredNorm()/members.size();
    group.ligand_weight=(ql.transpose()*c).squaredNorm()/members.size();
    // Keep complete mapped support, including weakly admixed acceptor groups.
    // A fixed one-percent-on-both-fragments screen loses real weak channels.
    if(group.centre_weight+group.ligand_weight<1e-8)return;
    M cross=c.transpose()*(pc*f*pl+pl*f*pc)*c;
    cross=(cross+cross.transpose()).eval()*0.5;
    Eigen::SelfAdjointEigenSolver<M> eig(cross,Eigen::EigenvaluesOnly);
    if(eig.info()!=Eigen::Success)return;
    group.cross_fock_min_hartree=eig.eigenvalues()[0];
    group.cross_fock_max_hartree=eig.eigenvalues()[eig.eigenvalues().size()-1];
    group.cross_fock_mean_hartree=cross.trace()/members.size();
    const double gate=std::max(1e-5,10*u.rows()*out.operator_max_error_hartree);
    if(group.cross_fock_max_hartree< -gate)group.character="bonding_mixing";
    else if(group.cross_fock_min_hartree>gate)group.character="antibonding_mixing";
    out.groups.push_back(std::move(group));
}
} // namespace

void populate_pi_coupling_modes(NboPiCoupling& channel,const NboMatrix& fock,
    const NboMatrix& canonical_columns,const NboMatrix* localized_columns,
    const std::vector<std::size_t>& source_columns) {
    channel.modes.clear();const auto n=fock.rows;
    const auto rectangular=[&](const NboMatrix& a){return a.rows==n&&a.columns>0&&
        a.values.size()==a.rows*a.columns&&std::all_of(a.values.begin(),a.values.end(),[](double x){return std::isfinite(x);});};
    if(!n||fock.columns!=n||!rectangular(fock)||!rectangular(canonical_columns)||
       !rectangular(channel.centre_projector_basis)||!rectangular(channel.ligand_projector_basis)||
       !std::isfinite(channel.operator_max_error_hartree)||channel.operator_max_error_hartree<0)return;
    const M f=mat(fock),cm=mat(canonical_columns),m=mat(channel.centre_projector_basis),l=mat(channel.ligand_projector_basis);
    if(maxabs(f-f.transpose())>std::max(2e-5,2*channel.operator_max_error_hartree)||
       maxabs(cm.transpose()*cm-M::Identity(cm.cols(),cm.cols()))>5e-5||
       maxabs(m.transpose()*m-M::Identity(m.cols(),m.cols()))>5e-5||
       maxabs(l.transpose()*l-M::Identity(l.cols(),l.cols()))>5e-5||maxabs(m.transpose()*l)>5e-5)return;
    std::vector<M> groups;groups.reserve(channel.groups.size());
    for(const auto& g:channel.groups){
        if(g.members.empty())return;M c(n,g.members.size());
        for(std::size_t j=0;j<g.members.size();++j){const auto member=g.members[j];
            if(!source_columns.empty()&&member>=source_columns.size())return;
            const auto column=source_columns.empty()?member:source_columns[member];
            if(column>=canonical_columns.columns)return;c.col(j)=cm.col(column);
        }
        groups.push_back(std::move(c));
    }
    const M block=m.transpose()*f*l;
    Eigen::JacobiSVD<M> svd(block,Eigen::ComputeThinU|Eigen::ComputeThinV);
    if(svd.info()!=Eigen::Success)return;
    // Entrywise operator uncertainty bounds its spectral norm by n*epsilon.
    // SVD residual is included uniformly; no molecule-dependent clustering.
    const double residual=(block-svd.matrixU()*svd.singularValues().asDiagonal()*svd.matrixV().transpose()).norm();
    const double error=std::max(1e-10,double(n)*channel.operator_max_error_hartree+residual);
    const double cluster_tolerance=std::max(1e-6,4*error);
    std::size_t rank=0;for(double s:svd.singularValues())if(s>cluster_tolerance)++rank;
    const bool localized=localized_columns&&rectangular(*localized_columns)&&localized_columns->columns==n&&
        maxabs(mat(*localized_columns).transpose()*mat(*localized_columns)-M::Identity(n,n))<=5e-5;
    const M nb=localized?mat(*localized_columns):M{};
    const auto matrix=[&](const M& v,const char* kind){NboMatrix out;out.kind=kind;out.spin=channel.spin;
        out.rows=v.rows();out.columns=v.cols();const RM flat=v;out.values.assign(flat.data(),flat.data()+flat.size());out.source=channel.source;return out;};
    for(std::size_t begin=0;begin<rank;){
        std::size_t end=begin+1;
        while(end<rank&&svd.singularValues()[end-1]-svd.singularValues()[end]<=cluster_tolerance)++end;
        const auto width=end-begin;const M qm=m*svd.matrixU().middleCols(begin,width),ql=l*svd.matrixV().middleCols(begin,width);
        const V singular=svd.singularValues().segment(begin,width);
        M cross=qm*singular.asDiagonal()*ql.transpose();cross=(cross+cross.transpose()).eval();
        double gap=svd.singularValues()[end-1];
        if(begin)gap=std::min(gap,svd.singularValues()[begin-1]-svd.singularValues()[begin]);
        if(end<static_cast<std::size_t>(svd.singularValues().size()))gap=std::min(gap,svd.singularValues()[end-1]-svd.singularValues()[end]);
        const double angle=std::min(1.0,2*error/std::max(gap,error));
        const double mapping=channel.direction_mapping_error_bound<1?channel.direction_mapping_error_bound:
            4*double(n)*channel.nao_orthogonality_error;
        NboPiCouplingMode mode;mode.id=channel.id+":mode:"+std::to_string(begin)+"-"+std::to_string(end);
        mode.rank=width;mode.singular_max_hartree=singular[0];mode.singular_min_hartree=singular[width-1];
        mode.singular_values_hartree.assign(singular.data(),singular.data()+singular.size());
        mode.numerical_coverage_bound=std::min(1.0,std::max(1e-8,mapping+2*angle+angle*angle));
        mode.centre_space_kind=channel.centre_family=="valence-d"?"valence-d":"auxiliary-pd";
        mode.ligand_space_kind=channel.ligand_space_kind;
        mode.centre_projector_basis=matrix(qm,"pi-mode-centre-basis");mode.ligand_projector_basis=matrix(ql,"pi-mode-ligand-basis");
        for(std::size_t i=0;i<groups.size();++i){const auto& c=groups[i];NboPiModeGroup g;g.group_index=i;
            g.centre_weight=(qm.transpose()*c).squaredNorm()/c.cols();g.ligand_weight=(ql.transpose()*c).squaredNorm()/c.cols();
            g.centre_coordinates=matrix(qm.transpose()*c,"pi-mode-canonical-centre-coordinates");
            g.ligand_coordinates=matrix(ql.transpose()*c,"pi-mode-canonical-ligand-coordinates");
            const M k=c.transpose()*cross*c;Eigen::SelfAdjointEigenSolver<M> eig((k+k.transpose())*.5,Eigen::EigenvaluesOnly);
            if(eig.info()!=Eigen::Success){channel.modes.clear();return;}
            g.cross_fock_min_hartree=eig.eigenvalues()[0];g.cross_fock_max_hartree=eig.eigenvalues()[eig.eigenvalues().size()-1];
            g.cross_fock_mean_hartree=k.trace()/c.cols();mode.groups.push_back(g);
        }
        if(localized)for(const auto& row:channel.direction_projection_evidence){
            if(!row.donor_nbo_id||!row.acceptor_nbo_id||row.donor_nbo_id>n||row.acceptor_nbo_id>n||
               row.energy_gap_hartree<=0||(row.status!="available"&&row.reason!="large-Fock-to-gap-ratio-use-block-strength")||
               (row.direction!="ligand_to_centre"&&row.direction!="centre_to_ligand"))continue;
            const V d=nb.col(row.donor_nbo_id-1),a=nb.col(row.acceptor_nbo_id-1);
            const bool forward=row.direction=="ligand_to_centre";
            const M& donor_basis=forward?ql:qm;const M& acceptor_basis=forward?qm:ql;
            const V dp=donor_basis*(donor_basis.transpose()*d),ap=acceptor_basis*(acceptor_basis.transpose()*a);
            const double projected=forward?(qm.transpose()*a).dot(singular.cwiseProduct(ql.transpose()*d)):
                (qm.transpose()*d).dot(singular.cwiseProduct(ql.transpose()*a));
            if(std::abs(projected)<=std::max(1e-10,10*double(n)*channel.operator_max_error_hartree))continue;
            NboPiModeEdge edge;edge.id=std::string(nbo_spin_name(channel.spin))+":"+std::to_string(row.donor_nbo_id)+":"+std::to_string(row.acceptor_nbo_id);
            edge.direction=row.direction;edge.donor_nbo_id=row.donor_nbo_id;edge.acceptor_nbo_id=row.acceptor_nbo_id;edge.projected_fock_hartree=projected;
            for(std::size_t i=0;i<groups.size();++i){const auto& c=groups[i];NboPiModeEdgeGroup eg;eg.group_index=i;
                eg.donor_weight=(c.transpose()*d).squaredNorm()/c.cols();eg.acceptor_weight=(c.transpose()*a).squaredNorm()/c.cols();
                eg.donor_mode_weight=(c.transpose()*dp).squaredNorm()/c.cols();eg.acceptor_mode_weight=(c.transpose()*ap).squaredNorm()/c.cols();
                if(eg.donor_mode_weight>1e-12||eg.acceptor_mode_weight>1e-12)edge.groups.push_back(eg);
            }
            mode.ordered_edges.push_back(std::move(edge));
        }
        channel.modes.push_back(std::move(mode));begin=end;
    }
}

PiModePairAssessment assess_pi_mode_pair(const NboPiCoupling& channel,
    const std::vector<std::size_t>& lower,const std::vector<std::size_t>& upper,
    const PiEndpointSymmetryEvidence& ls,const PiEndpointSymmetryEvidence& us) {
    PiModePairAssessment out;std::set<std::string> centre_kinds;
    if(ls.verified&&us.verified&&!ls.scope_id.empty()&&ls.scope_id==us.scope_id&&
       !ls.irrep.empty()&&!us.irrep.empty()&&ls.irrep!=us.irrep){out.reason="incompatible-verified-same-domain-irreps";return out;}
    const auto same=[](auto a,auto b){std::sort(a.begin(),a.end());std::sort(b.begin(),b.end());return a==b;};
    const auto find=[&](const auto& members){return std::find_if(channel.groups.begin(),channel.groups.end(),[&](const auto& g){return same(g.members,members);});};
    const auto lo=find(lower),hi=find(upper);
    if(lower.empty()||upper.empty()||lo==channel.groups.end()||hi==channel.groups.end()||lo==hi){out.reason="complete-mode-endpoint-members-unavailable";return out;}
    const auto li=static_cast<std::size_t>(lo-channel.groups.begin()),ui=static_cast<std::size_t>(hi-channel.groups.begin());
    std::array<std::map<std::size_t,double>,2> low_roles,high_roles;
    std::array<bool,2> directions{};std::set<std::string> edges;
    std::vector<double> common_space_coverage(channel.groups.size(),0);
    out.numerical_coverage_bound=0;
    for(const auto& mode:channel.modes){
        const auto mg=[&](std::size_t i){return std::find_if(mode.groups.begin(),mode.groups.end(),[&](const auto& g){return g.group_index==i;});};
        const auto a=mg(li),b=mg(ui);const double gate=mode.numerical_coverage_bound;
        if(a==mode.groups.end()||b==mode.groups.end()||!std::isfinite(gate)||gate<0||gate>=1)continue;
        const auto valid_coordinates=[&](const NboMatrix& c,std::size_t columns){return c.rows==mode.rank&&c.columns==columns&&
            c.values.size()==c.rows*c.columns&&std::all_of(c.values.begin(),c.values.end(),[](double v){return std::isfinite(v);});};
        if(mode.singular_values_hartree.size()!=mode.rank||!mode.rank||
           !valid_coordinates(a->centre_coordinates,lower.size())||!valid_coordinates(a->ligand_coordinates,lower.size())||
           !valid_coordinates(b->centre_coordinates,upper.size())||!valid_coordinates(b->ligand_coordinates,upper.size()))continue;
        const Eigen::Map<const V> singular(mode.singular_values_hartree.data(),mode.rank);
        // Do not sum these operators before taking norms: the full canonical
        // F_ij is diagonal and can cancel even for a genuine shared channel.
        const M ml=mat(a->centre_coordinates).transpose()*singular.asDiagonal()*mat(b->ligand_coordinates);
        const M lm=mat(a->ligand_coordinates).transpose()*singular.asDiagonal()*mat(b->centre_coordinates);
        const double joint=std::hypot(ml.norm(),lm.norm());
        const double joint_error=std::max(1e-10,gate*mode.singular_max_hartree+
            double(channel.centre_projector_basis.rows)*channel.operator_max_error_hartree);
        if(!std::isfinite(joint)||joint<=joint_error)continue;
        const double energy_gate=std::max(1e-8,10*channel.operator_max_error_hartree);
        if(a->centre_weight+a->ligand_weight<=gate||b->centre_weight+b->ligand_weight<=gate||
           a->cross_fock_mean_hartree>=-energy_gate||a->cross_fock_max_hartree>energy_gate||
           b->cross_fock_mean_hartree<=energy_gate||b->cross_fock_min_hartree< -energy_gate)continue;
        out.shared_mode_ids.push_back(mode.id);out.lower_coverage+=a->centre_weight+a->ligand_weight;
        centre_kinds.insert(mode.centre_space_kind);
        out.upper_coverage+=b->centre_weight+b->ligand_weight;out.numerical_coverage_bound=std::max(out.numerical_coverage_bound,gate);
        out.lower_cross_fock_mean_hartree+=a->cross_fock_mean_hartree;
        out.upper_cross_fock_mean_hartree+=b->cross_fock_mean_hartree;
        out.shared_fragment_contraction_norm_hartree=std::hypot(out.shared_fragment_contraction_norm_hartree,joint);
        for(const auto& g:mode.groups)if(g.group_index<common_space_coverage.size())
            common_space_coverage[g.group_index]+=g.centre_weight+g.ligand_weight;
        const auto participating=std::count_if(mode.groups.begin(),mode.groups.end(),[&](const auto& g){return g.centre_weight+g.ligand_weight>gate;});
        out.multi_group_relation=out.multi_group_relation||participating>2;
        for(const auto& edge:mode.ordered_edges){
            const auto eg=[&](std::size_t i){return std::find_if(edge.groups.begin(),edge.groups.end(),[&](const auto& g){return g.group_index==i;});};
            const auto x=eg(li),y=eg(ui);if(x==edge.groups.end()||y==edge.groups.end())continue;
            const double direct=std::min(x->donor_mode_weight,y->acceptor_mode_weight),reverse=std::min(y->donor_mode_weight,x->acceptor_mode_weight);
            if(std::max(direct,reverse)<=gate)continue;
            const int direction=edge.direction=="ligand_to_centre"?0:edge.direction=="centre_to_ligand"?1:-1;
            if(direction<0)continue;directions[direction]=true;edges.insert(edge.id);
            if(direct>=reverse){low_roles[direction][edge.donor_nbo_id]=x->donor_weight;high_roles[direction][edge.acceptor_nbo_id]=y->acceptor_weight;}
            else {low_roles[direction][edge.acceptor_nbo_id]=x->acceptor_weight;high_roles[direction][edge.donor_nbo_id]=y->donor_weight;}
        }
    }
    out.verified=!out.shared_mode_ids.empty();
    if(!out.verified){out.reason="no-common-resolved-fragment-coupling-mode";return out;}
    out.channel_family=centre_kinds==std::set<std::string>{"valence-d"}?"valence-d-pi":
        centre_kinds==std::set<std::string>{"auxiliary-p"}?"auxiliary-p-pi":
        centre_kinds==std::set<std::string>{"higher-radial-d"}?"higher-radial-pi":"mixed-auxiliary-pi";
    std::vector<std::size_t> major_groups;
    for(std::size_t i=0;i<common_space_coverage.size();++i)
        if(common_space_coverage[i]-out.numerical_coverage_bound>=out.primary_coverage_floor)major_groups.push_back(i);
    out.two_endpoint_relation=major_groups.size()==2&&
        std::find(major_groups.begin(),major_groups.end(),li)!=major_groups.end()&&
        std::find(major_groups.begin(),major_groups.end(),ui)!=major_groups.end();
    out.matched_edge_ids.assign(edges.begin(),edges.end());
    const auto sum=[](const auto& roles){double w=0;for(const auto& [id,value]:roles)w+=value;return w;};
    for(int d=0;d<2;++d)if(directions[d]){out.lower_role_coverage=std::max(out.lower_role_coverage,sum(low_roles[d]));out.upper_role_coverage=std::max(out.upper_role_coverage,sum(high_roles[d]));}
    const double full=channel.spin==NboSpin::Total?2:1;
    if(channel.occupation_status=="available"&&lo->occupation_per_mo&&hi->occupation_per_mo&&
       std::abs(*lo->occupation_per_mo-full)<=1e-6&&std::abs(*hi->occupation_per_mo-full)<=1e-6){
        out.direction="occupied_space_mixing";out.direction_verified=true;
    }else if(directions[0]||directions[1]){out.direction_verified=true;out.direction=directions[0]&&directions[1]?"bidirectional":directions[0]?"ligand_to_centre":"centre_to_ligand";}
    const bool family=channel.centre_family=="valence-d"&&
        (channel.ligand_space_kind=="internal-pi-bonding"||channel.ligand_space_kind=="internal-pi-antibonding"||channel.ligand_space_kind=="localized-pi-lone-pair"||channel.ligand_space_kind=="localized-pi-vacancy");
    out.ordinary_display_eligible=out.two_endpoint_relation&&family&&out.direction_verified&&out.direction!="occupied_space_mixing"&&out.direction!="bidirectional"&&
        std::min(out.lower_coverage,out.upper_coverage)-out.numerical_coverage_bound>=out.primary_coverage_floor&&
        std::min(out.lower_role_coverage,out.upper_role_coverage)-out.numerical_coverage_bound>=out.primary_coverage_floor;
    out.reason=out.ordinary_display_eligible?"same-mode-same-ordered-edge-majority-valence-d-pi":
        out.direction_verified?"same-mode-source-direction-secondary-or-outside-valence-d-primary-scope":"same-mode-mixing-without-linked-ordered-direction";
    return out;
}

PiModeNetworkAssessment assess_pi_mode_network(const NboPiCoupling& channel,
    const NboPiCouplingMode& mode) {
    PiModeNetworkAssessment out;out.channel_id=channel.id;out.mode_id=mode.id;
    out.spin=nbo_spin_name(channel.spin);out.ligand_space_kind=mode.ligand_space_kind;
    out.channel_family=mode.centre_space_kind=="valence-d"?"valence-d-pi":
        mode.centre_space_kind=="auxiliary-p"?"auxiliary-p-pi":
        mode.centre_space_kind=="higher-radial-d"?"higher-radial-pi":"mixed-auxiliary-pi";
    out.numerical_coverage_bound=mode.numerical_coverage_bound;
    if(channel.id.empty()||mode.id.empty()||!mode.rank||mode.groups.empty()||
       !std::isfinite(mode.numerical_coverage_bound)||mode.numerical_coverage_bound<0||mode.numerical_coverage_bound>=1){
        out.reason="resolved-fragment-mode-unavailable";return out;
    }
    const double gate=mode.numerical_coverage_bound;
    std::vector<std::size_t> primary;
    std::map<std::size_t,std::size_t> node_index;
    for(const auto& g:mode.groups){
        if(g.group_index>=channel.groups.size()||g.centre_weight+g.ligand_weight<=gate)continue;
        const auto& source=channel.groups[g.group_index];if(source.members.empty())continue;
        PiModeNetworkNode node;node.group_index=g.group_index;node.members=source.members;
        node.centre_weight=g.centre_weight;node.ligand_weight=g.ligand_weight;
        const double energy_gate=std::max(1e-8,10*channel.operator_max_error_hartree);
        node.character=g.cross_fock_mean_hartree< -energy_gate&&g.cross_fock_max_hartree<=energy_gate?"bonding_mixing":
            g.cross_fock_mean_hartree>energy_gate&&g.cross_fock_min_hartree>= -energy_gate?"antibonding_mixing":"mixed_or_unresolved";
        node.primary=node.centre_weight+node.ligand_weight-gate>=out.primary_coverage_floor;
        if(node.primary)primary.push_back(g.group_index);
        node_index[g.group_index]=out.nodes.size();out.nodes.push_back(std::move(node));
    }
    out.verified=!out.nodes.empty();if(!out.verified){out.reason="mode-without-resolved-canonical-support";return out;}
    std::vector<std::map<std::size_t,double>> donor_roles(out.nodes.size()),acceptor_roles(out.nodes.size());
    bool forward=false,reverse=false,major_edge=false;
    std::set<std::string> edge_ids;
    for(const auto& edge:mode.ordered_edges){
        double donor_span=0,acceptor_span=0,donor_mode=0,acceptor_mode=0;
        double donor_full_mode=0,acceptor_full_mode=0;
        for(const auto& g:edge.groups){const auto found=node_index.find(g.group_index);if(found==node_index.end())continue;
            const auto index=found->second;const auto& node=out.nodes[index];
            donor_roles[index][edge.donor_nbo_id]=g.donor_weight;
            acceptor_roles[index][edge.acceptor_nbo_id]=g.acceptor_weight;
            donor_full_mode+=double(node.members.size())*g.donor_mode_weight;
            acceptor_full_mode+=double(node.members.size())*g.acceptor_mode_weight;
            if(node.primary){const double rank=double(node.members.size());
                donor_span+=rank*g.donor_weight;acceptor_span+=rank*g.acceptor_weight;
                donor_mode+=rank*g.donor_mode_weight;acceptor_mode+=rank*g.acceptor_mode_weight;}
        }
        // The same source edge must map both roles into the retained primary
        // network span. No donor/acceptor union creates a fictitious edge.
        if(std::min(donor_full_mode,acceptor_full_mode)<=gate)continue;
        edge_ids.insert(edge.id);forward=forward||edge.direction=="ligand_to_centre";
        reverse=reverse||edge.direction=="centre_to_ligand";
        major_edge=major_edge||(std::min(donor_mode,acceptor_mode)>gate&&std::min(donor_span,acceptor_span)-gate>=out.primary_coverage_floor);
    }
    for(std::size_t i=0;i<out.nodes.size();++i){
        for(const auto& [id,weight]:donor_roles[i])out.nodes[i].donor_role_coverage+=weight;
        for(const auto& [id,weight]:acceptor_roles[i])out.nodes[i].acceptor_role_coverage+=weight;
    }
    out.matched_edge_ids.assign(edge_ids.begin(),edge_ids.end());
    const double full=channel.spin==NboSpin::Total?2:1;
    const bool occupied=primary.size()>1&&channel.occupation_status=="available"&&
        std::all_of(primary.begin(),primary.end(),[&](auto i){return channel.groups[i].occupation_per_mo&&
            std::abs(*channel.groups[i].occupation_per_mo-full)<=1e-6;});
    if(occupied){out.direction="occupied_space_mixing";out.direction_verified=true;}
    else if(forward||reverse){out.direction_verified=true;out.direction=forward&&reverse?"bidirectional":forward?"ligand_to_centre":"centre_to_ligand";}
    if(primary.size()==2){const auto& a=channel.groups[primary[0]].members;const auto& b=channel.groups[primary[1]].members;
        const auto pair=assess_pi_mode_pair(channel,a,b);const auto inverse=assess_pi_mode_pair(channel,b,a);
        out.two_endpoint_relation=(pair.verified&&pair.two_endpoint_relation)||(inverse.verified&&inverse.two_endpoint_relation);
    }
    const bool family=mode.centre_space_kind=="valence-d"&&
        (mode.ligand_space_kind=="internal-pi-bonding"||mode.ligand_space_kind=="internal-pi-antibonding"||
         mode.ligand_space_kind=="localized-pi-lone-pair"||mode.ligand_space_kind=="localized-pi-vacancy");
    out.ordinary_display_eligible=family&&major_edge&&primary.size()>1&&out.direction_verified&&
        out.direction!="occupied_space_mixing"&&out.direction!="bidirectional";
    out.reason=out.two_endpoint_relation?"resolved-mode-network-with-dominant-two-endpoint-reduction":
        out.ordinary_display_eligible?"major-linked-direction-in-multigroup-mode-network":"resolved-mode-network-research-scope";
    return out;
}

PiFrozenOperatorAssessment assess_pi_frozen_operator(const NboMatrix& fock,
    const NboMatrix& centre_basis,const NboMatrix& ligand_basis,
    const NboMatrix& canonical_columns,
    const std::vector<std::vector<std::size_t>>& groups,
    const std::vector<double>& occupations,double operator_error,std::vector<PiFrozenSpectralGroup>* per_group) {
    PiFrozenOperatorAssessment out;
    if(per_group)per_group->clear();
    std::vector<PiFrozenSpectralGroup> saved_groups;
    const auto fail=[&](const char* reason){out.reason=reason;return out;};
    const std::size_t n=fock.rows;
    const auto rectangular=[&](const NboMatrix& a){return a.rows==n&&a.columns>0&&
        a.columns<=n&&a.values.size()==n*a.columns&&
        std::all_of(a.values.begin(),a.values.end(),[](double x){return std::isfinite(x);});};
    if(!n||!full(&fock,n)||!full(&canonical_columns,n)||!rectangular(centre_basis)||
       !rectangular(ligand_basis)||groups.empty()||!std::isfinite(operator_error)||operator_error<0)
        return fail("incomplete-frozen-operator-input");
    const M f=mat(fock),c=mat(canonical_columns),qm=mat(centre_basis),ql=mat(ligand_basis);
    if(maxabs(f-f.transpose())>2e-5||maxabs(c.transpose()*c-M::Identity(n,n))>5e-5||
       maxabs(qm.transpose()*qm-M::Identity(qm.cols(),qm.cols()))>5e-5||
       maxabs(ql.transpose()*ql-M::Identity(ql.cols(),ql.cols()))>5e-5||
       maxabs(qm.transpose()*ql)>1e-8)return fail("nonorthogonal-frozen-operator-scope");
    const M diagonal=c.transpose()*f*c;
    const V energies=diagonal.diagonal();
    if(maxabs(diagonal-M(energies.asDiagonal()))>std::max(2e-5,operator_error))
        return fail("canonical-members-do-not-diagonalize-this-operator");
    const M cross=qm*(qm.transpose()*f*ql)*ql.transpose();
    const M removal=cross+cross.transpose();
    Eigen::SelfAdjointEigenSolver<M> off((f+f.transpose())*0.5-removal);
    if(off.info()!=Eigen::Success||!off.eigenvalues().allFinite())return fail("frozen-operator-eigensolver-failed");
    Eigen::JacobiSVD<M> coupling(qm.transpose()*f*ql);
    out.removal_norm_hartree=coupling.singularValues()[0];
    // The complete residual is available: its Frobenius norm bounds its
    // spectral norm without the very loose n*maximum-entry estimate. The
    // second term bounds the eigenvalue distortion from C^T C differing from I.
    const double residual_bound=(f-c*energies.asDiagonal()*c.transpose()).norm();
    const double metric_bound=energies.cwiseAbs().maxCoeff()*(c.transpose()*c-M::Identity(n,n)).norm();
    out.numerical_error_bound_hartree=2*std::max(operator_error,residual_bound+metric_bound);
    std::vector<std::size_t> order(n),position(n);std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return energies[a]<energies[b];});
    for(std::size_t i=0;i<n;++i)position[order[i]]=i;
    out.max_energy_shift_hartree=0;out.max_group_width_change_hartree=0;
    out.max_subspace_sin2=0;out.minimum_external_gap_hartree=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<n;++i)out.max_energy_shift_hartree=std::max(out.max_energy_shift_hartree,
        std::abs(energies[order[i]]-off.eigenvalues()[i]));
    out.full_spectrum_max_energy_shift_hartree=out.max_energy_shift_hartree;
    bool tracked=true;std::set<std::size_t> covered;
    for(const auto& group:groups) {
        if(group.empty())return fail("empty-spectral-group");
        bool group_tracked=true;
        std::vector<std::size_t> ranks;
        for(auto index:group){if(index>=n||!covered.insert(index).second)return fail("invalid-spectral-group-members");ranks.push_back(position[index]);}
        std::sort(ranks.begin(),ranks.end());
        for(std::size_t k=1;k<ranks.size();++k)if(ranks[k]!=ranks[k-1]+1)group_tracked=false;
        const auto first=ranks.front(),last=ranks.back();
        double gap=std::numeric_limits<double>::infinity();
        if(first)gap=std::min({gap,energies[order[first]]-energies[order[first-1]],off.eigenvalues()[first]-off.eigenvalues()[first-1]});
        if(last+1<n)gap=std::min({gap,energies[order[last+1]]-energies[order[last]],off.eigenvalues()[last+1]-off.eigenvalues()[last]});
        out.minimum_external_gap_hartree=std::min(out.minimum_external_gap_hartree,gap);
        if(gap<=std::max(1e-7,2*out.numerical_error_bound_hartree))group_tracked=false;
        M original(n,group.size()),changed(n,group.size());
        for(std::size_t j=0;j<group.size();++j){original.col(j)=c.col(group[j]);changed.col(j)=off.eigenvectors().col(ranks[j]);}
        Eigen::JacobiSVD<M> overlap(original.transpose()*changed);
        const double sin2=std::clamp(1-std::pow(overlap.singularValues().minCoeff(),2),0.0,1.0);
        out.max_subspace_sin2=std::max(out.max_subspace_sin2,sin2);
        if(sin2>=0.5)group_tracked=false;
        const double width=energies[order[last]]-energies[order[first]];
        out.max_group_width_change_hartree=std::max(out.max_group_width_change_hartree,
            std::abs(off.eigenvalues()[last]-off.eigenvalues()[first]-width));
        PiFrozenSpectralGroup saved;saved.members=group;auto& diagnostic=saved.assessment;
        diagnostic=out;diagnostic.available=true;diagnostic.tracking_verified=group_tracked;
        diagnostic.max_energy_shift_hartree=0;
        for(auto r:ranks)diagnostic.max_energy_shift_hartree=std::max(diagnostic.max_energy_shift_hartree,
            std::abs(energies[order[r]]-off.eigenvalues()[r]));
        diagnostic.max_group_width_change_hartree=std::abs(off.eigenvalues()[last]-off.eigenvalues()[first]-width);
        diagnostic.max_subspace_sin2=sin2;diagnostic.minimum_external_gap_hartree=gap;
        saved_groups.push_back(std::move(saved));tracked=tracked&&group_tracked;
    }
    out.frontier_gap_change_hartree=0;out.occupation_boundary_preserved=occupations.size()==n;
    if(out.occupation_boundary_preserved) {
        for(std::size_t i=0;i<n;++i)if(!std::isfinite(occupations[i])||occupations[i]<0)out.occupation_boundary_preserved=false;
        for(std::size_t i=1;i<n&&out.occupation_boundary_preserved;++i) {
            if(std::abs(occupations[order[i]]-occupations[order[i-1]])<1e-6)continue;
            const double oldgap=energies[order[i]]-energies[order[i-1]],newgap=off.eigenvalues()[i]-off.eigenvalues()[i-1];
            out.frontier_gap_change_hartree=std::max(out.frontier_gap_change_hartree,std::abs(newgap-oldgap));
            if(std::min(oldgap,newgap)<=std::max(1e-7,2*out.numerical_error_bound_hartree))out.occupation_boundary_preserved=false;
            M occupied(n,i);for(std::size_t j=0;j<i;++j)occupied.col(j)=c.col(order[j]);
            Eigen::JacobiSVD<M> overlap(occupied.transpose()*off.eigenvectors().leftCols(i));
            if(overlap.singularValues().minCoeff()*overlap.singularValues().minCoeff()<=0.5)out.occupation_boundary_preserved=false;
        }
    }
    if(per_group){
        for(auto& group:saved_groups){
            group.assessment.occupation_boundary_preserved=out.occupation_boundary_preserved;
            group.assessment.frontier_gap_change_hartree=out.frontier_gap_change_hartree;
            group.assessment.reason=group.assessment.tracking_verified&&out.occupation_boundary_preserved?
                "frozen-complete-source-spectral-cluster":"cluster-tracking-or-global-occupation-boundary-not-certified";
        }
        *per_group=std::move(saved_groups);
    }
    out.available=true;out.tracking_verified=tracked;
    out.reason=tracked&&out.occupation_boundary_preserved?"frozen-same-operator-complete-spectral-groups":
        "spectral-tracking-or-occupation-boundary-not-certified";
    return out;
}

PiFrozenOperatorAssessment assess_pi_frozen_display_scope(const NboPiCoupling& channel,
    const std::vector<std::size_t>& members){
    auto out=channel.frozen_operator;out.tracking_verified=false;
    if(!out.available||members.empty()){out.available=false;out.reason="stable-display-domain-unavailable";return out;}
    const std::set<std::size_t> requested(members.begin(),members.end());
    if(requested.size()!=members.size()){out.available=false;out.reason="duplicate-stable-display-member";return out;}
    std::set<std::size_t> covered;
    out.max_energy_shift_hartree=out.max_group_width_change_hartree=out.max_subspace_sin2=0;
    out.minimum_external_gap_hartree=std::numeric_limits<double>::infinity();out.tracking_verified=true;
    for(const auto& group:channel.frozen_groups){
        const auto count=std::count_if(group.members.begin(),group.members.end(),[&](auto i){return requested.count(i);});
        if(!count)continue;
        if(static_cast<std::size_t>(count)!=group.members.size()){
            out.available=false;out.tracking_verified=false;out.reason="stable-domain-cuts-a-degenerate-spectral-cluster";return out;
        }
        const auto& a=group.assessment;
        out.tracking_verified=out.tracking_verified&&a.tracking_verified;
        out.max_energy_shift_hartree=std::max(out.max_energy_shift_hartree,a.max_energy_shift_hartree);
        out.max_group_width_change_hartree=std::max(out.max_group_width_change_hartree,a.max_group_width_change_hartree);
        out.max_subspace_sin2=std::max(out.max_subspace_sin2,a.max_subspace_sin2);
        out.minimum_external_gap_hartree=std::min(out.minimum_external_gap_hartree,a.minimum_external_gap_hartree);
        covered.insert(group.members.begin(),group.members.end());
    }
    if(covered!=requested){out.available=false;out.tracking_verified=false;out.reason="stable-display-domain-not-completely-mapped";return out;}
    out.reason=out.tracking_verified&&out.occupation_boundary_preserved?
        "frozen-stable-complete-display-domain-global-frontier-protected":"display-domain-or-global-frontier-not-certified";
    return out;
}

PiDisplayCalibration pi_channel_display_calibration(const NboPiCoupling& channel){
    PiDisplayCalibration calibration;
    if(channel.operator_kind!="canonical-same-operator"||!channel.localized_family_verified||
       channel.occupation_status!="available"||channel.ligand_family_nbo_ids.empty()||
       (channel.ligand_family!="internal-pi-bonding"&&channel.ligand_family!="internal-pi-antibonding")||
       (channel.centre_family!="valence-d"&&channel.centre_family!="additional-pd-acceptors")||
       channel.centre_projector_basis.values.empty()||channel.ligand_projector_basis.values.empty())return calibration;
    calibration.version="cov.pi-display.2026-10-03.v1";
    calibration.family="complete-internal-pi-pd/canonical-same-operator";
    calibration.energy_budget_ev=0.025;calibration.subspace_sin2_budget=0.02;
    // The predeclared 0.025 eV / 2% budget was checked on saved carbonyl
    // training scopes and a closed-shell cyanide held-out scope, with strong,
    // near-resonant, missing-operator and physical-spin controls retained.
    // It is a bounded display-error calibration for this structural family,
    // not a universal chemical-strength threshold or a deletion energy.
    calibration.validated=true;
    return calibration;
}

NboPiCouplingAnalysis analyse_nbo_pi_couplings(const Wavefunction& w,
        const NboIntegration& data,
        const std::vector<std::pair<std::size_t,std::size_t>>& strong_connectivity){
    NboPiCouplingAnalysis out;
    std::vector<std::string> projection_failures;
    std::optional<NboSalcModel> physical_operator_model;
    if(!data.dataset.association.compatible||data.canonical_fingerprint!=nbo_canonical_fingerprint(w)){
        out.status="rejected";out.reason="NBO association does not match immutable canonical identity";
        return out;
    }
    const auto n=static_cast<std::size_t>(w.basis_count);
    if(!n||!data.dataset.archive){out.reason="Complete same-source archive is unavailable";return out;}
    const auto* s=unique_matrix(data.dataset.archive->matrices,"OVERLAP",NboSpin::Total);
    if(!full(s,n)){out.reason="Complete archive AO overlap is unavailable";return out;}
    const auto ro=verify_nbo_restricted_open_shell(w,data);
    const bool shared_ro=ro.verified&&
        (w.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared||
         w.orbital_occupation_model==OrbitalOccupationModel::SharedIntegerDeterminant);
    for(auto spin:{NboSpin::Total,NboSpin::Alpha,NboSpin::Beta}){
        const auto* a=unique_matrix(data.dataset.matrices,"AONAO",spin);
        if(!a&&spin!=NboSpin::Total)
            a=unique_matrix(data.dataset.matrices,"AONAO",NboSpin::Total);
        const auto* u=unique_matrix(data.dataset.matrices,"NAOMO",spin);
        const auto* f=unique_matrix(data.dataset.archive->matrices,"FOCK",spin);
        if(!full(a,n)||!full(u,n)||!full(f,n))continue;
        const NboNaoValidation* validation=nullptr;
        for(const auto& candidate:data.dataset.nao_validation)
            if(candidate.spin==spin&&candidate.available&&
               (candidate.direct_fchk_coefficients||(shared_ro&&spin==NboSpin::Beta)))
                validation=&candidate;
        if(!validation||validation->effective_mo_columns!=n)continue;
        const M an=mat(*a),un=mat(*u),fn=an.transpose()*mat(*f)*an;
        const double orth=maxabs(an.transpose()*mat(*s)*an-M::Identity(n,n));
        if(orth>5e-5||maxabs(un.transpose()*un-M::Identity(n,n))>5e-5)continue;
        V energies(n),occupations(n);std::vector<bool> seen(n,false);
        bool occupation_provenance=true;
        for(const auto& mo:w.orbitals){
            if(!shared_ro&&(spin==NboSpin::Beta)!=(mo.spin==Spin::Beta))continue;
            if(shared_ro&&mo.spin!=Spin::Alpha)continue;
            const auto index=mo.source_orbital_index;
            if(index>=n||seen[index])continue;
            seen[index]=true;energies[index]=mo.energy_hartree;
            occupation_provenance=occupation_provenance&&
                (mo.occupation_provenance!=DataProvenance::Unavailable ||
                 w.electron_counts_provenance!=DataProvenance::Unavailable);
        }
        if(!std::all_of(seen.begin(),seen.end(),[](bool yes){return yes;}))continue;
        const double canonical_operator_residual=maxabs(fn-un*energies.asDiagonal()*un.transpose());
        double operator_error=canonical_operator_residual;
        std::string operator_kind="canonical-same-operator";
        double operator_validation_tolerance=2e-5;
        if(operator_error>2e-5) {
            if(!physical_operator_model)physical_operator_model=build_nbo_salc_model(w,data);
            const auto physical=std::find_if(physical_operator_model->energies.begin(),
                physical_operator_model->energies.end(),[&](const auto& e){
                    return e.spin==spin&&e.available&&e.printed_operator_verified;});
            const auto op=std::find_if(physical_operator_model->spin_operators.begin(),
                physical_operator_model->spin_operators.end(),[&](const auto& e){return e.spin==spin;});
            if(physical==physical_operator_model->energies.end()||
               op==physical_operator_model->spin_operators.end()||op->fock.size()!=n*n||
               op->basis.size()!=n) {
                projection_failures.push_back(std::string(nbo_spin_name(spin))+": canonical operator differs and independent printed physical-spin Fock validation is unavailable");
                continue;
            }
            bool same_basis=true;for(std::size_t i=0;i<n;++i)
                same_basis=same_basis&&op->basis[i].kind==NboOrbitalKind::NAO&&
                    op->basis[i].spin==spin&&op->basis[i].index==i;
            const M verified=Eigen::Map<const RM>(op->fock.data(),n,n);
            if(!same_basis||maxabs(fn-verified)>2e-5)continue;
            operator_error=std::max(maxabs(fn-verified),physical->hermiticity_error);
            operator_kind="independently-verified-physical-spin-fock";
            operator_validation_tolerance=5.1e-6;
        }
        bool density_validated=occupation_provenance &&
            validation->occupation_density_verified &&
            validation->canonical_occupations.size()==n &&
            std::all_of(validation->canonical_occupations.begin(),
                validation->canonical_occupations.end(),[](double x){return std::isfinite(x);});
        M density;
        if(density_validated){
            for(std::size_t col=0;col<n;++col)
                occupations[col]=validation->canonical_occupations[col];
            density=un*occupations.asDiagonal()*un.transpose();
        }
        if(!density_validated&&shared_ro&&physical_operator_model) {
            const auto op=std::find_if(physical_operator_model->spin_operators.begin(),
                physical_operator_model->spin_operators.end(),[&](const auto& e){return e.spin==spin;});
            if(op!=physical_operator_model->spin_operators.end()&&op->density.size()==n*n) {
                const std::size_t electrons=spin==NboSpin::Beta?w.beta_electrons:w.alpha_electrons;
                for(std::size_t i=0;i<n;++i)occupations[i]=i<electrons?1.0:0.0;
                const M reconstructed=un*occupations.asDiagonal()*un.transpose();
                const M independently_verified=Eigen::Map<const RM>(op->density.data(),n,n);
                if(maxabs(reconstructed-independently_verified)<=2e-5) {
                    density=independently_verified;density_validated=true;
                }
            }
        }
        std::optional<M> nbo_transform;
        const auto* naob=unique_matrix(data.dataset.matrices,"NAONBO",spin);
        const auto* aob=unique_matrix(data.dataset.matrices,"AONBO",spin);
        const auto* nbo_cap=nbo_capability(data,"nbo");
        if(nbo_cap&&nbo_cap->available()&&full(naob,n)&&full(aob,n)&&
           maxabs(an*mat(*naob)-mat(*aob))<=2e-5&&
           maxabs(mat(*naob).transpose()*mat(*naob)-M::Identity(n,n))<=5e-5)
            nbo_transform=mat(*naob);
        std::optional<M> nbo_fock;
        if(nbo_transform)nbo_fock=nbo_transform->transpose()*fn*(*nbo_transform);
        const E2Convention convention=nbo_fock?e2_convention(data,*nbo_fock,spin):E2Convention{};
        // Define resolvable source spectral clusters before any pi block is
        // removed. Near-degenerate members inside the source error bound form
        // one whole cluster; their internal spectral width is still tested.
        const M canonical_fock=un.transpose()*fn*un;
        const V canonical_expectations=canonical_fock.diagonal();
        const double source_matrix_bound=(fn-un*canonical_expectations.asDiagonal()*un.transpose()).norm()+
            canonical_expectations.cwiseAbs().maxCoeff()*(un.transpose()*un-M::Identity(n,n)).norm();
        const double source_cluster_tolerance=std::max(1e-6,4*std::max(operator_error,source_matrix_bound));
        const auto family=shells(data,spin);
        std::set<std::pair<std::size_t,std::size_t>> bonded;
        for(const auto& edge:strong_connectivity)bonded.emplace(std::minmax(edge.first,edge.second));
        struct Candidate{std::vector<std::size_t> centre,ligand,ci,li;M qc,ql;
            std::string ligand_family,centre_family;std::vector<std::size_t> family_nbo_ids;
            std::vector<NboPiAngularEvidence> angular;};
        std::vector<Candidate> candidates;
        const auto add_angular=[&](Candidate& c,std::size_t atom,const PiFrame& frame){
            NboPiAngularEvidence x;x.atom=atom;
            x.p_metric_min_eigenvalue=frame.minimum_p;
            x.max_sigma_leakage=frame.leakage;
            x.centre_projection=frame.centre;
            x.partition_residual=frame.residual;
            x.pi_generalized_eigenvalues=frame.spectrum;
            c.angular.push_back(x);
        };
        for(const auto& [atom,d_ids]:family.d){
            Candidate candidate;candidate.centre={atom};candidate.ci=d_ids;
            candidate.centre_family="valence-d";std::vector<std::size_t> extra_ids;
            // Empty valence acceptors can be labelled Ryd by NBO (including
            // transition-metal np shells). Search all noncore p/d radial shells;
            // the ligand projector still specifies the local pi symmetry.
            for(const auto& nao:data.dataset.naos)
                if(nao.spin==spin&&nao.atom==atom+1&&nao.id&&
                   nao.type.find("Cor")==std::string::npos&&
                   !nao.angular.empty()&&(nao.angular[0]=='p'||nao.angular[0]=='d')&&
                   std::find(candidate.ci.begin(),candidate.ci.end(),nao.id)==candidate.ci.end())
                    extra_ids.push_back(nao.id);
            candidate.qc=d_columns(n,candidate.ci);
            std::vector<M> parts;bool valid=true;
            for(const auto& [other,p_ids]:family.p){
                if(!bonded.count(std::minmax(atom,other))||atom==other)continue;
                const V3 axis=bond_axis(w,atom,other);if(axis.norm()<0.5)continue;
                std::string reason;
                auto frame=derived_pi_frame(w,data,other,p_ids,spin,axis,reason);
                if(!frame){projection_failures.push_back(std::string(nbo_spin_name(spin))+" atom "+std::to_string(other+1)+": "+reason);valid=false;break;}
                candidate.ligand.push_back(other);
                candidate.li.insert(candidate.li.end(),p_ids.begin(),p_ids.end());
                parts.push_back(frame->q);add_angular(candidate,other,*frame);
            }
            if(parts.empty()||!valid)continue;
            candidate.ql=M::Zero(n,2*parts.size());
            for(std::size_t k=0;k<parts.size();++k)
                candidate.ql.middleCols(2*k,2)=parts[k];
            candidates.push_back(candidate);
            if(!extra_ids.empty()){
                candidate.ci=std::move(extra_ids);candidate.qc=d_columns(n,candidate.ci);
                candidate.centre_family="additional-pd-acceptors";
                candidates.push_back(std::move(candidate));
            }
        }
        // Resolve the complete ligand's internal occupied pi and antibonding
        // pi families before reducing the d--ligand coupling. Atom-pi alone
        // combines these physically distinct spaces (e.g. both ends of any
        // heteronuclear multiple bond). This uses no ligand-name catalogue or
        // E2 printing threshold, and never assumes BD ordinal 1 means sigma.
        if(nbo_transform) {
            const std::size_t broad_count=candidates.size();
            for(std::size_t cindex=0;cindex<broad_count;++cindex) {
                const auto broad=candidates[cindex];
                std::set<std::size_t> fragment(broad.ligand.begin(),broad.ligand.end());
                bool changed=true;
                while(changed) {
                    changed=false;
                    for(const auto& [a,b]:bonded) {
                        if(family.d.count(a)||family.d.count(b))continue;
                        if(fragment.count(a)&&fragment.insert(b).second)changed=true;
                        if(fragment.count(b)&&fragment.insert(a).second)changed=true;
                    }
                }
                std::set<std::size_t> fragment_p_rows;
                for(const auto& row:data.dataset.naos)
                    if(row.spin==spin&&row.id&&row.atom&&fragment.count(row.atom-1))
                        fragment_p_rows.insert(row.id-1);
                for(const bool antibonding:{false,true}) {
                    std::vector<V> columns;std::vector<std::size_t> ids;
                    for(const auto& local:data.dataset.orbitals) {
                        if(local.spin!=spin||!local.id||local.id>static_cast<std::size_t>(nbo_transform->cols())||
                           local.kind!=(antibonding?"BD*":"BD"))continue;
                        const auto* desc=nbo_orbital(data,{NboOrbitalKind::NBO,spin,local.id-1});
                        if(!desc||desc->atoms.size()<2||
                           !std::all_of(desc->atoms.begin(),desc->atoms.end(),[&](auto atom){return fragment.count(atom);})||
                           !std::any_of(desc->channels.begin(),desc->channels.end(),[](const auto& evidence){
                               return evidence.status=="available"&&evidence.channel=="pi";}))continue;
                        const V original=nbo_transform->col(local.id-1);
                        double direct_p=0;for(auto id:broad.li)direct_p+=original[id-1]*original[id-1];
                        const double transverse=(broad.ql.transpose()*original).squaredNorm();
                        // Establish that this internal pi family belongs to the
                        // metal-facing pi geometry, not a sigma-facing orbital.
                        if(direct_p<1e-4||transverse/direct_p<0.9)continue;
                        V projected=V::Zero(n);for(auto row:fragment_p_rows)projected[row]=original[row];
                        if(projected.squaredNorm()<0.7)continue;
                        columns.push_back(std::move(projected));ids.push_back(local.id);
                    }
                    if(columns.empty())continue;
                    M raw(n,columns.size());for(std::size_t i=0;i<columns.size();++i)raw.col(i)=columns[i];
                    Eigen::JacobiSVD<M> basis(raw,Eigen::ComputeThinU);
                    const double cutoff=std::max(1e-8,1e-8*basis.singularValues()[0]);
                    std::size_t rank=0;for(double value:basis.singularValues())if(value>cutoff)++rank;
                    if(!rank||rank!=columns.size())continue; // incomplete dependent family is not silently truncated
                    Candidate split=broad;split.ligand.assign(fragment.begin(),fragment.end());
                    split.li.clear();for(auto row:fragment_p_rows)split.li.push_back(row+1);
                    split.ql=basis.matrixU().leftCols(rank);
                    split.ligand_family=antibonding?"internal-pi-antibonding":"internal-pi-bonding";
                    split.family_nbo_ids=std::move(ids);
                    candidates.push_back(std::move(split));
                }
            }
        }
        for(const auto& [atom,other]:bonded){
            const auto ca=family.p.find(atom),lb=family.p.find(other);
            if(ca==family.p.end()||lb==family.p.end())continue;
            const V3 axis=bond_axis(w,atom,other);if(axis.norm()<0.5)continue;
            std::string centre_reason,ligand_reason;
            auto centre_frame=derived_pi_frame(w,data,atom,ca->second,spin,axis,centre_reason);
            auto ligand_frame=derived_pi_frame(w,data,other,lb->second,spin,axis,ligand_reason);
            if(!centre_frame||!ligand_frame){
                projection_failures.push_back(std::string(nbo_spin_name(spin))+" atoms "+
                    std::to_string(atom+1)+"-"+std::to_string(other+1)+": "+centre_reason+" "+ligand_reason);
                continue;
            }
            Candidate candidate;candidate.centre={atom};candidate.ligand={other};
            candidate.ci.assign(ca->second.begin(),ca->second.end());
            candidate.li.assign(lb->second.begin(),lb->second.end());
            candidate.qc=centre_frame->q;
            candidate.ql=ligand_frame->q;
            add_angular(candidate,atom,*centre_frame);
            add_angular(candidate,other,*ligand_frame);
            candidates.push_back(std::move(candidate));
        }
        for(const auto& candidate:candidates){
            if(maxabs(candidate.qc.transpose()*candidate.ql)>1e-8)continue;
            const M block=candidate.qc.transpose()*fn*candidate.ql;
            Eigen::JacobiSVD<M> svd(block,Eigen::ComputeThinU|Eigen::ComputeThinV);
            std::size_t rank=0;for(auto value:svd.singularValues())if(value>1e-6)++rank;
            if(!rank)continue;
            NboPiCoupling record;
            record.spin=spin;record.centre_atoms=candidate.centre;
            record.ligand_family=candidate.ligand_family;record.centre_family=candidate.centre_family;
            if(!candidate.ligand_family.empty())record.ligand_space_kind=candidate.ligand_family;
            record.ligand_family_nbo_ids=candidate.family_nbo_ids;
            record.localized_family_verified=!candidate.ligand_family.empty();
            record.ligand_atoms=candidate.ligand;
            record.centre_nao_ids=candidate.ci;record.ligand_nao_ids=candidate.li;
            record.coupled_rank=rank;record.operator_max_error_hartree=operator_error;
            record.operator_kind=operator_kind;
            record.canonical_operator_residual_hartree=canonical_operator_residual;
            record.operator_validation_tolerance_hartree=operator_validation_tolerance;
            record.canonical_members_are_verified_shared_spatial=shared_ro;
            record.nao_orthogonality_error=orth;record.source=f->source;
            record.angular_projector_evidence=candidate.angular;
            record.minimum_centre_projection=1;
            for(const auto& evidence:candidate.angular){
                record.angular_leakage=std::max(record.angular_leakage,evidence.max_sigma_leakage);
                record.minimum_centre_projection=std::min(record.minimum_centre_projection,
                    evidence.centre_projection);
            }
            for(auto value:svd.singularValues())
                record.singular_values_hartree.push_back(value);
            std::ostringstream id;id<<"pi:"<<nbo_spin_name(spin)<<":";
            for(auto atom:candidate.centre)id<<atom<<',';id<<":";
            for(auto atom:candidate.ligand)id<<atom<<',';
            id<<":"<<candidate.ligand_family<<":"<<candidate.centre_family;record.id=id.str();
            const M qc=candidate.qc*svd.matrixU().leftCols(rank);
            const M ql=candidate.ql*svd.matrixV().leftCols(rank);
            const M pc=qc*qc.transpose(),pl=ql*ql.transpose();
            record.centre_onsite_hartree=(qc.transpose()*fn*qc).trace()/rank;
            record.ligand_onsite_hartree=(ql.transpose()*fn*ql).trace()/rank;
            const auto range=[&](const M& block){
                Eigen::SelfAdjointEigenSolver<M> solver((block+block.transpose())*0.5,
                    Eigen::EigenvaluesOnly);
                if(solver.info()!=Eigen::Success||!solver.eigenvalues().allFinite())
                    return std::optional<std::array<double,2>>{};
                return std::optional<std::array<double,2>>(
                    std::array<double,2>{solver.eigenvalues()[0],
                        solver.eigenvalues()[solver.eigenvalues().size()-1]});
            };
            const auto centre_fock_range=range(qc.transpose()*fn*qc),
                       ligand_fock_range=range(ql.transpose()*fn*ql);
            if(!centre_fock_range||!ligand_fock_range)continue;
            record.centre_onsite_range_hartree=*centre_fock_range;
            record.ligand_onsite_range_hartree=*ligand_fock_range;
            if(density_validated){
                const auto centre_density=qc.transpose()*density*qc,
                           ligand_density=ql.transpose()*density*ql;
                record.centre_occupation=centre_density.trace()/rank;
                record.ligand_occupation=ligand_density.trace()/rank;
                record.centre_occupation_range=range(centre_density);
                record.ligand_occupation_range=range(ligand_density);
                if(!record.centre_occupation_range||!record.ligand_occupation_range)
                    continue;
                record.occupation_status="available";
                record.occupation_reason="Validated same-spin canonical occupations reconstruct NAO density";
            }else record.occupation_reason=
                "Occupation provenance or independent canonical-to-NAO density validation unavailable; coupling retained without direction";
            for(const auto& group:canonical_groups(w,shared_ro?NboSpin::Alpha:spin))
                append_group(record,w,un,fn,candidate.qc,candidate.ql,
                    candidate.qc*candidate.qc.transpose(),candidate.ql*candidate.ql.transpose(),
                    density_validated?&occupations:nullptr,group);
            if(record.groups.empty())continue;
            E2Support forward,reverse;
            if(nbo_transform){
                forward=directed_e2(data,*nbo_transform,*nbo_fock,candidate.qc,candidate.ql,candidate.centre,
                                    candidate.ligand,spin,true,operator_error,convention);
                reverse=directed_e2(data,*nbo_transform,*nbo_fock,candidate.qc,candidate.ql,candidate.centre,
                                    candidate.ligand,spin,false,operator_error,convention);
                record.direction_projection_evidence=forward.rows;
                record.direction_projection_evidence.insert(record.direction_projection_evidence.end(),
                    reverse.rows.begin(),reverse.rows.end());
            }
            if(nbo_transform)record.direction_mapping_error_bound=std::min(1.0,4*n*std::max({orth,
                maxabs(nbo_transform->transpose()*(*nbo_transform)-M::Identity(n,n)),
                maxabs(an*(*nbo_transform)-mat(*aob))}));
            if(nbo_transform)for(auto& group:record.groups){
                const auto mapped=[&](const E2Support& support,bool donor){
                    std::set<std::size_t> ids;
                    for(const auto& row:support.rows)
                        if(row.energy_gap_hartree>0&&(row.status=="available"||row.reason=="large-Fock-to-gap-ratio-use-block-strength"))
                            ids.insert((donor?row.donor_nbo_id:row.acceptor_nbo_id)-1);
                    double trace=0;
                    for(auto id:ids)for(auto member:group.members){
                        const double overlap=nbo_transform->col(id).dot(un.col(w.orbitals[member].source_orbital_index));
                        trace+=overlap*overlap;
                    }
                    return trace/std::max<std::size_t>(1,group.members.size());
                };
                group.ligand_to_centre_donor_weight=mapped(forward,true);
                group.ligand_to_centre_acceptor_weight=mapped(forward,false);
                group.centre_to_ligand_donor_weight=mapped(reverse,true);
                group.centre_to_ligand_acceptor_weight=mapped(reverse,false);
            }
            record.direction_verified=forward.verified||reverse.verified;
            if(candidate.ligand_family.empty())for(const auto& edge:record.direction_projection_evidence){
                const bool forward_edge=edge.direction=="ligand_to_centre";
                const auto id=forward_edge?edge.donor_nbo_id:edge.acceptor_nbo_id;
                const auto local=std::find_if(data.dataset.orbitals.begin(),data.dataset.orbitals.end(),[&](const auto& o){return o.spin==spin&&o.id==id;});
                if(local!=data.dataset.orbitals.end()&&forward_edge&&local->kind=="LP")record.ligand_space_kind="localized-pi-lone-pair";
                else if(local!=data.dataset.orbitals.end()&&!forward_edge&&local->kind=="LV"&&record.ligand_space_kind=="transverse-p-basis")record.ligand_space_kind="localized-pi-vacancy";
            }
            if(forward.verified&&reverse.verified){
                record.direction="bidirectional";
                record.direction_evidence="Independent occupied-to-acceptor NBO pi channels in both directions; no net-charge inference";
            }else if(forward.verified||reverse.verified){
                const auto& evidence=forward.verified?forward:reverse;
                record.direction=forward.verified?"ligand_to_centre":"centre_to_ligand";
                record.direction_donor_weight=evidence.donor;record.direction_acceptor_weight=evidence.acceptor;
                record.direction_evidence="Complete ordered NBO occupied/acceptor search with same-spin Fock and pi-scope mapping; complex projected occupations are diagnostic only";
            }else if(record.centre_occupation_range&&record.ligand_occupation_range&&
                (*record.centre_occupation_range)[0]>(spin==NboSpin::Total?1.9:0.95)&&
                (*record.ligand_occupation_range)[0]>(spin==NboSpin::Total?1.9:0.95)){
                record.direction="occupied_space_mixing";
                record.direction_evidence="Both complete selected pi subspaces are occupied; no ordered occupied-to-acceptor channel in this scope; not a Pauli energy or whole-molecule absence claim";
            }else if(record.centre_occupation_range&&record.ligand_occupation_range&&
                std::abs(record.centre_onsite_range_hartree[0]-record.ligand_onsite_range_hartree[0])<1e-4&&
                std::abs(record.centre_onsite_range_hartree[1]-record.ligand_onsite_range_hartree[1])<1e-4&&
                std::abs((*record.centre_occupation_range)[0]-(*record.ligand_occupation_range)[0])<0.05&&
                std::abs((*record.centre_occupation_range)[1]-(*record.ligand_occupation_range)[1])<0.05){
                record.direction="symmetric_coupled";
                record.direction_evidence="Symmetric onsite and occupation spectra in the coupled pi scope; no donor/acceptor assignment";
            }else record.direction_evidence=nbo_transform?
                "No resolved ordered occupied-to-acceptor pi channel in the complete localized search; inspect recorded nonperturbative or scope limitations":
                "Complete same-source localized NBO transform unavailable; obtain occupied/acceptor evidence";
            const auto matrix=[&](const M& value,const char* kind){NboMatrix x;x.kind=kind;x.spin=spin;
                x.rows=value.rows();x.columns=value.cols();const RM rowmajor=value;
                x.values.assign(rowmajor.data(),rowmajor.data()+rowmajor.size());x.source=f->source;return x;};
            // Preserve complete declared scope before numerical SVD truncation.
            record.centre_projector_basis=matrix(candidate.qc,"pi-centre-projector-basis");
            record.ligand_projector_basis=matrix(candidate.ql,"pi-ligand-projector-basis");
            std::vector<std::size_t> canonical_source_columns;canonical_source_columns.reserve(w.orbitals.size());
            for(const auto& mo:w.orbitals)canonical_source_columns.push_back(mo.source_orbital_index);
            std::optional<NboMatrix> localized_columns;
            if(nbo_transform)localized_columns=matrix(*nbo_transform,"NBO-in-NAO");
            populate_pi_coupling_modes(record,matrix(fn,"Fock-NAO"),*u,
                localized_columns?&*localized_columns:nullptr,canonical_source_columns);
            for(auto& mode:record.modes){
                const M basis=mat(mode.centre_projector_basis);
                std::map<std::pair<std::string,std::string>,double> shell_weights;
                for(const auto& nao:data.dataset.naos)if(nao.spin==spin&&nao.id&&nao.id<=static_cast<std::size_t>(basis.rows())&&!nao.angular.empty())
                    shell_weights[{nao.type,nao.angular.substr(0,1)}]+=basis.row(nao.id-1).squaredNorm()/mode.rank;
                double p_weight=0,d_weight=0;
                for(const auto& [identity,weight]:shell_weights)if(weight>1e-12){
                    mode.centre_shells.push_back({identity.first,identity.second,weight});
                    if(identity.second=="p")p_weight+=weight;if(identity.second=="d")d_weight+=weight;
                }
                if(record.centre_family!="valence-d")mode.centre_space_kind=p_weight>1-1e-5?"auxiliary-p":d_weight>1-1e-5?"higher-radial-d":"auxiliary-pd";
            }
            std::vector<std::size_t> sorted(n);std::iota(sorted.begin(),sorted.end(),0);
            std::stable_sort(sorted.begin(),sorted.end(),[&](auto a,auto b){return energies[a]<energies[b];});
            std::vector<std::vector<std::size_t>> spectral_groups;
            for(auto i:sorted){if(spectral_groups.empty()||std::abs(energies[i]-energies[spectral_groups.back().back()])>source_cluster_tolerance)
                spectral_groups.push_back({});spectral_groups.back().push_back(i);}
            std::vector<double> saved_occupations;
            if(density_validated)saved_occupations.assign(occupations.data(),occupations.data()+occupations.size());
            record.frozen_operator=assess_pi_frozen_operator(matrix(fn,"Fock-NAO"),
                record.centre_projector_basis,record.ligand_projector_basis,*u,spectral_groups,
                saved_occupations,operator_error,&record.frozen_groups);
            for(auto& group:record.frozen_groups)for(auto& member:group.members){
                const auto mo=std::find_if(w.orbitals.begin(),w.orbitals.end(),[&](const auto& x){
                    return x.source_orbital_index==member&&(shared_ro?x.spin==Spin::Alpha:
                        ((spin==NboSpin::Beta)==(x.spin==Spin::Beta)));});
                if(mo==w.orbitals.end()){record.frozen_operator.available=false;
                    record.frozen_operator.reason="incomplete-canonical-source-index-map";continue;}
                member=static_cast<std::size_t>(mo-w.orbitals.begin());
            }
            out.couplings.push_back(std::move(record));
        }
        out.evidence.push_back(f->source);
    }
    if(out.couplings.empty()){
        out.status=out.evidence.empty()||!projection_failures.empty()?"insufficient_evidence":"not_applicable";
        out.reason=out.evidence.empty()?"No spin has a complete matched NAO transform and independently consistent same-source Fock operator":
            !projection_failures.empty()?"Same-source Fock verified; local pi projection evidence is incomplete":
            "No supported bonded valence p-pi or d-to-ligand-p-pi coupling block";
    }else{
        out.status="available";
        out.reason="Independently verified local-operator NAO projectors; each channel declares canonical or physical-spin operator identity; canonical groups and eigenvalue intervals remain separate";
    }
    if(!projection_failures.empty()){
        out.reason+="; candidate pi subspaces could not be verified";
        for(const auto& failure:projection_failures)out.reason+="; "+failure;
    }
    return out;
}
} // namespace cov
