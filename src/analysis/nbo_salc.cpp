#include "cov/molecular_point_group_frame.hpp"
#include "cov/nbo_salc.hpp"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace cov { namespace {
using M=Eigen::MatrixXd;
using V=Eigen::VectorXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
double err(const M& a){return a.size()?a.cwiseAbs().maxCoeff():0;}
M matrix(const NboMatrix& a){return Eigen::Map<const RM>(a.values.data(),a.rows,a.columns);}
const NboMatrix* find(const std::vector<NboMatrix>& a,const char* kind,NboSpin spin){
    const NboMatrix* p=nullptr;for(const auto& x:a)if(x.kind==kind&&x.spin==spin){if(p)return nullptr;p=&x;}return p;
}
bool complete(const NboMatrix* a,std::size_t n){return a&&a->rows==n&&a->columns==n&&a->values.size()==n*n&&std::all_of(a->values.begin(),a->values.end(),[](double v){return std::isfinite(v);});}
std::string compact(std::string s){s.erase(std::remove_if(s.begin(),s.end(),[](unsigned char c){return std::isspace(c);}),s.end());return s;}
std::string atoms_label(const Wavefunction& w,const std::vector<std::size_t>& a){std::ostringstream s;for(auto i:a){if(s.tellp()>0)s<<", ";s<<w.atoms[i].symbol<<i+1;}return s.str();}
std::string key_atoms(const std::vector<std::size_t>& a){std::string s;for(auto i:a)s+=":"+std::to_string(i);return s;}
std::vector<double> flat(const M& a){RM r=a;return {r.data(),r.data()+r.size()};}
std::array<double,9> product(const std::array<double,9>& a,const std::array<double,9>& b){std::array<double,9> c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
SymmetryOperation product(const SymmetryOperation& a,const SymmetryOperation& b){SymmetryOperation c;c.matrix=product(a.matrix,b.matrix);c.atom_permutation.resize(b.atom_permutation.size());for(std::size_t i=0;i<c.atom_permutation.size();++i)c.atom_permutation[i]=a.atom_permutation[b.atom_permutation[i]];return c;}
double op_error(const SymmetryOperation& a,const SymmetryOperation& b){if(a.atom_permutation!=b.atom_permutation)return 1;double e=0;for(int i=0;i<9;++i)e=std::max(e,std::abs(a.matrix[i]-b.matrix[i]));return e;}
struct Group {std::vector<SymmetryOperation> ops, generators;std::vector<std::size_t> parent, step, table;std::vector<std::vector<std::size_t>> classes;bool valid=false;double error=0;std::string detail;};
std::size_t lookup(const std::vector<SymmetryOperation>& ops,const SymmetryOperation& p){for(std::size_t i=0;i<ops.size();++i)if(op_error(ops[i],p)<2e-5)return i;return ops.size();}
Group make_group(const Wavefunction& w,const MolecularSymmetry& symmetry,const NboSalcOptions& options){
    Group out;SymmetryOperation identity;identity.atom_permutation.resize(w.atoms.size());std::iota(identity.atom_permutation.begin(),identity.atom_permutation.end(),0);
    out.ops={identity};out.parent={0};out.step={0};
    auto candidates=symmetry.operations;
    // Linear geometry supplies only E/i. Add a finite, explicitly recorded
    // sampling subgroup which resolves angular momenta through the actual lmax.
    if(symmetry.linear&&w.atoms.size()>1){
        Eigen::Vector3d axis(w.atoms.back().x-w.atoms.front().x,w.atoms.back().y-w.atoms.front().y,w.atoms.back().z-w.atoms.front().z);
        if(axis.norm()>1e-10){axis.normalize();unsigned lmax=0;for(const auto& s:w.shells)lmax=std::max(lmax,unsigned(s.angular_momentum));const int n=std::max(3,int(2*lmax+1));
            Eigen::Matrix3d rot=Eigen::AngleAxisd(2*3.14159265358979323846/n,axis).toRotationMatrix();
            Eigen::Vector3d seed=std::abs(axis.x())<.8?Eigen::Vector3d::UnitX():Eigen::Vector3d::UnitY();Eigen::Vector3d normal=axis.cross(seed).normalized();
            Eigen::Matrix3d reflect=Eigen::Matrix3d::Identity()-2*normal*normal.transpose();
            for(const auto& m:std::vector<Eigen::Matrix3d>{rot,reflect}){auto op=identity;for(int i=0;i<3;++i)for(int j=0;j<3;++j)op.matrix[3*i+j]=m(i,j);candidates.push_back(op);}
        }
    }
    for(const auto& candidate:candidates){
        if(candidate.atom_permutation.size()!=w.atoms.size())continue;
        bool legal=true;for(std::size_t i=0;i<w.atoms.size();++i){const auto j=candidate.atom_permutation[i];if(j>=w.atoms.size()||w.atoms[i].atomic_number!=w.atoms[j].atomic_number||std::abs(w.atoms[i].nuclear_charge-w.atoms[j].nuclear_charge)>1e-8){legal=false;break;}}
        if(!legal){out.detail="Symmetry exchanges inequivalent nuclei/ECP centres";return out;}
        if(lookup(out.ops,candidate)<out.ops.size())continue;
        out.generators.push_back(candidate);
        for(std::size_t i=0;i<out.ops.size();++i)for(std::size_t j=0;j<out.generators.size();++j){auto p=product(out.ops[i],out.generators[j]);if(lookup(out.ops,p)==out.ops.size()){
            if(out.ops.size()>=options.maximum_group_order){out.detail="Operation set does not close within finite-group capacity";return out;}
            double geom_error=0;for(std::size_t a=0;a<w.atoms.size();++a){const auto& from=w.atoms[a];const auto& to=w.atoms[p.atom_permutation[a]];double x[3]={from.x-symmetry.centre_bohr[0],from.y-symmetry.centre_bohr[1],from.z-symmetry.centre_bohr[2]};double y[3]={to.x-symmetry.centre_bohr[0],to.y-symmetry.centre_bohr[1],to.z-symmetry.centre_bohr[2]};double e=0;for(int r=0;r<3;++r){double d=-y[r];for(int k=0;k<3;++k)d+=p.matrix[3*r+k]*x[k];e+=d*d;}geom_error=std::max(geom_error,std::sqrt(e));}
            if(geom_error>symmetry.tolerance_bohr){out.detail="Products of approximate geometry operations fail atom mapping";return out;}
            p.max_mapping_error_bohr=geom_error;out.ops.push_back(p);out.parent.push_back(i);out.step.push_back(j);
        }}
    }
    const auto n=out.ops.size();out.table.resize(n*n);std::vector<std::size_t> inverse(n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){auto p=product(out.ops[i],out.ops[j]);auto k=lookup(out.ops,p);if(k==n){out.detail="Finite group multiplication did not close";return out;}out.table[i*n+j]=k;out.error=std::max(out.error,op_error(out.ops[k],p));if(k==0)inverse[i]=j;}
    std::vector<bool> used(n);for(std::size_t i=0;i<n;++i)if(!used[i]){std::set<std::size_t> c;for(std::size_t h=0;h<n;++h)c.insert(out.table[out.table[h*n+i]*n+inverse[h]]);out.classes.emplace_back(c.begin(),c.end());for(auto j:c)used[j]=true;}
    out.valid=true;out.detail=symmetry.linear?"Validated finite sampling subgroup of linear point group":"Validated complete finite group";return out;
}
NboElectronicSymmetryScope density_scope(const Wavefunction& w,const MolecularSymmetry& geometry,
    const Group& group,const NboSalcOptions& options){
    NboElectronicSymmetryScope out;out.geometry_group=geometry.point_group;out.naming_group=geometry.point_group;
    out.operations=group.ops;out.geometry_operations=group.ops;out.relative_tolerance=options.symmetry_tolerance;
    const auto n=std::size_t(w.basis_count);
    if(!group.valid||geometry.linear||w.ao_overlap.size()!=n*n){out.status="operation_frame_unavailable";return out;}
    const auto unpack=[&](const std::vector<double>& packed){
        if(packed.size()!=n*(n+1)/2)return M{};M p(n,n);std::size_t k=0;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j)p(i,j)=p(j,i)=packed[k++];return p.allFinite()?p:M{};};
    const M total=unpack(w.total_density_packed),spin=unpack(w.spin_density_packed);
    if(!total.size()){out.status="total_density_unavailable";return out;}
    const bool open=w.alpha_electrons!=w.beta_electrons||w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin;
    if(open&&!spin.size()){out.status="spin_density_unavailable";return out;}
    std::vector<const MolecularOrbital*> columns;for(const auto& mo:w.orbitals)if(mo.spin!=Spin::Beta&&mo.coefficients.size()==n)columns.push_back(&mo);
    if(columns.empty()){out.status="retained_metric_basis_unavailable";return out;}
    const M s=Eigen::Map<const RM>(w.ao_overlap.data(),n,n);M q(n,columns.size());
    for(std::size_t j=0;j<columns.size();++j)q.col(j)=Eigen::Map<const V>(columns[j]->coefficients.data(),n);
    const M gram=q.transpose()*s*q;
    if(!gram.allFinite()||err(gram-M::Identity(q.cols(),q.cols()))>options.metric_tolerance){out.status="retained_metric_basis_invalid";return out;}
    Eigen::SelfAdjointEigenSolver<M> eig((gram+gram.transpose())*.5);
    if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=1e-10){out.status="retained_metric_rank_invalid";return out;}
    q=(q*eig.eigenvectors()*eig.eigenvalues().array().rsqrt().matrix().asDiagonal()*eig.eigenvectors().transpose()).eval();
    const M sq=s*q,pt=sq.transpose()*total*sq,ps=spin.size()?M(sq.transpose()*spin*sq):M{};
    const auto density_complete=[&](const M& original,const M& represented){
        const M residual=original-q*represented*q.transpose(),weighted=s*residual;
        const double norm2=(weighted*weighted).trace();
        return std::sqrt(std::max(0.,norm2))/std::max(1e-12,represented.norm())<=options.metric_tolerance;};
    if(!density_complete(total,pt)||(spin.size()&&!density_complete(spin,ps))){out.status="density_outside_retained_basis";return out;}
    std::vector<bool> keep;std::vector<M> reps;
    for(const auto& operation:group.ops){
        const auto values=apply_orbital_symmetry_operation(w,operation,flat(q),q.cols());
        if(values.size()!=n*std::size_t(q.cols())){out.status="density_operation_unavailable";return out;}
        const M tq=Eigen::Map<const RM>(values.data(),n,q.cols()),d=sq.transpose()*tq;
        const M residual=tq-q*d,leakage=residual.transpose()*s*residual;
        if(err(d.transpose()*d-M::Identity(d.cols(),d.cols()))>options.metric_tolerance||leakage.diagonal().maxCoeff()>options.symmetry_tolerance*options.symmetry_tolerance){out.status="density_operation_not_closed";return out;}
        reps.push_back(d);
        const double rt=(d*pt*d.transpose()-pt).norm()/std::max(1e-12,pt.norm());
        const double rs=ps.size()?(d*ps*d.transpose()-ps).norm()/std::max(1e-12,ps.norm()):0;
        out.total_density_residuals.push_back(rt);if(ps.size())out.spin_density_residuals.push_back(rs);
        keep.push_back(rt<=out.relative_tolerance&&rs<=out.relative_tolerance);
    }
    for(std::size_t i=0;i<reps.size();++i)for(std::size_t j=0;j<reps.size();++j)
        if(err(reps[i]*reps[j]-reps[group.table[i*reps.size()+j]])>options.symmetry_tolerance){out.status="density_operation_relations_failed";return out;}
    out.density_checked=true;
    // Enumerate only subgroups contained in the accepted operation set. This
    // does not label a threshold-selected, non-closed set as a point group.
    const auto close=[&](std::vector<std::size_t> members){
        std::set<std::size_t> chosen(members.begin(),members.end());bool changed=true;
        while(changed){changed=false;const std::vector<std::size_t> old(chosen.begin(),chosen.end());
            for(auto i:old)for(auto j:old){const auto k=group.table[i*group.ops.size()+j];if(!keep[k])return std::vector<std::size_t>{};if(chosen.insert(k).second)changed=true;}}
        return std::vector<std::size_t>(chosen.begin(),chosen.end());};
    std::vector<std::size_t> accepted;for(std::size_t i=0;i<keep.size();++i)if(keep[i])accepted.push_back(i);
    std::vector<std::size_t> chosen=close(accepted);
    if(chosen.empty()){
        std::set<std::vector<std::size_t>> seen{{0}};std::vector<std::vector<std::size_t>> todo{{0}};
        for(std::size_t cursor=0;cursor<todo.size();++cursor){if(todo.size()>4096){out.status="density_subgroup_search_capacity";return out;}
            for(auto i:accepted){auto candidate=todo[cursor];candidate.push_back(i);candidate=close(candidate);if(!candidate.empty()&&seen.insert(candidate).second)todo.push_back(std::move(candidate));}}
        std::size_t largest=0,count=0;for(const auto& candidate:todo){if(candidate.size()>largest){largest=candidate.size();count=1;chosen=candidate;}else if(candidate.size()==largest)++count;}
        if(count!=1){out.status="ambiguous_maximal_density_subgroups";return out;}
    }
    if(chosen.size()==group.ops.size()){out.naming_scope_verified=true;out.status="geometry_group_preserves_density";out.retained_operation_indices=chosen;return out;}
    MolecularSymmetry subgroup=geometry;subgroup.operations.clear();for(auto i:chosen)subgroup.operations.push_back(group.ops[i]);
    std::vector<std::string> candidates{"C1","Cs","Ci","T","Td","Th","O","Oh","I","Ih"};
    for(std::size_t order=2;order<=group.ops.size();++order)for(const auto& family:{"C","D","S"})
        for(const auto& suffix:{"","v","h","d"})candidates.push_back(std::string(family)+std::to_string(order)+suffix);
    std::vector<std::string> matching;
    for(const auto& name:candidates){if(finite_point_group_order(name)!=chosen.size())continue;subgroup.point_group=name;
        if(molecular_point_group_irreps(w,subgroup).valid)matching.push_back(name);}
    // Standard aliases S2=Ci, C1v/C1h=Cs and D1=C2 are not separate scopes.
    matching.erase(std::remove(matching.begin(),matching.end(),"S2"),matching.end());
    if(matching.size()!=1){out.status="density_subgroup_name_ambiguous";return out;}
    out.naming_group=matching.front();out.operations=subgroup.operations;out.retained_operation_indices=chosen;
    out.reduced=true;out.naming_scope_verified=true;out.status="verified_approximate_density_subgroup";return out;
}

M deterministic_basis(const M& eigvectors,double tolerance){
    // Stable projector-column gauge, independent of eigenvector signs and
    // rotations within a repeated eigenvalue. Phase is a coordinate convention.
    M projector=eigvectors*eigvectors.transpose(),q=M::Zero(projector.rows(),eigvectors.cols());
    for(Eigen::Index k=0;k<q.cols();++k){
        V best;double largest=-1;
        for(Eigen::Index i=0;i<projector.cols();++i){V v=projector.col(i);for(int pass=0;pass<2;++pass)for(Eigen::Index j=0;j<k;++j)v-=q.col(j).dot(v)*q.col(j);const double norm=v.squaredNorm();if(norm>largest+1e-14){largest=norm;best=std::move(v);}}
        // A rank-d orthogonal projector always has a substantial remaining
        // column. Never normalize a tiny early column and amplify roundoff.
        if(largest<tolerance*tolerance)return {};
        best/=std::sqrt(largest);Eigen::Index pivot;best.cwiseAbs().maxCoeff(&pivot);if(best[pivot]<0)best=-best;q.col(k)=best;
    }
    if(err(q*q.transpose()-projector)>1e-10||err(q.transpose()*q-M::Identity(q.cols(),q.cols()))>1e-10)return {};
    return q;
}
std::vector<M> merge_group_connected_spaces(std::vector<M> spaces,const std::vector<M>& reps,double tolerance){
    // An approximate class sum can split a repeated eigenvalue. Its numerical
    // eigenvectors are not separate symmetry spaces when the actual group
    // mixes them. Repair the partition, never the underlying AO operations.
    // Frobenius block norms and normalized characters are invariant under
    // arbitrary orthogonal changes of gauge within each candidate block.
    const auto n=spaces.size();if(n<2)return spaces;
    std::vector<std::size_t> parent(n);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](std::size_t i){while(parent[i]!=i)i=parent[i];return i;};
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j){double coupling=0;for(const auto& r:reps)coupling=std::max(coupling,(spaces[i].transpose()*r*spaces[j]).norm());if(coupling>tolerance)parent[root(j)]=root(i);}
    auto collect=[&](){std::map<std::size_t,std::vector<std::size_t>> components;for(std::size_t i=0;i<n;++i)components[root(i)].push_back(i);std::vector<M> result;for(const auto& [key,indices]:components){Eigen::Index count=0;for(auto i:indices)count+=spaces[i].cols();M joined(spaces.front().rows(),count);Eigen::Index offset=0;for(auto i:indices){joined.middleCols(offset,spaces[i].cols())=spaces[i];offset+=spaces[i].cols();}auto q=deterministic_basis(joined,1e-9);if(q.size()==0)return std::vector<M>{};result.push_back(std::move(q));}return result;};
    spaces=collect();if(spaces.empty())return spaces;
    // Copies of the same one-dimensional irrep can remain individually
    // invariant despite numerical splitting. Restore the whole isotypic
    // projector using its full operation-character signature.
    std::vector<V> characters;std::vector<bool> invariant;
    for(const auto& q:spaces){V chars(reps.size());bool valid=true;for(std::size_t g=0;g<reps.size();++g){M local=q.transpose()*reps[g]*q;chars[g]=local.trace()/double(q.cols());if(err(reps[g]*q-q*local)>tolerance)valid=false;}characters.push_back(std::move(chars));invariant.push_back(valid);}
    const auto m=spaces.size();parent.resize(m);std::iota(parent.begin(),parent.end(),0);
    for(std::size_t i=0;i<m;++i)for(std::size_t j=i+1;j<m;++j)if(invariant[i]&&invariant[j]&&(characters[i]-characters[j]).cwiseAbs().maxCoeff()<=tolerance)parent[root(j)]=root(i);
    std::map<std::size_t,std::vector<std::size_t>> components;for(std::size_t i=0;i<m;++i)components[root(i)].push_back(i);std::vector<M> result;
    for(const auto& [key,indices]:components){Eigen::Index count=0;for(auto i:indices)count+=spaces[i].cols();M joined(spaces.front().rows(),count);Eigen::Index offset=0;for(auto i:indices){joined.middleCols(offset,spaces[i].cols())=spaces[i];offset+=spaces[i].cols();}auto q=deterministic_basis(joined,1e-9);if(q.size()==0)return {};result.push_back(std::move(q));}return result;
}
// A central class sum identifies isotypic spaces, but is scalar on their
// multiplicity coordinates. Split those coordinates with commuting operators.
// Every result is rechecked against the complete real character row; this also
// handles chemical E rows of real character norm 2 without splitting partners.
double copy_leakage(const M& q,const M& rep,const M& outside_metric){
    const M internal=rep*q-q*(q.transpose()*rep*q);
    const M residual_metric=internal.transpose()*internal+q.transpose()*outside_metric*q;
    return std::sqrt(std::max(0.0,residual_metric.diagonal().maxCoeff()));
}
std::vector<M> split_irrep_copies(const M& space,const std::vector<M>& reps,const std::vector<M>& outside_metrics,
    const PointGroupIrrepTable& table,const M& fock,const NboSalcOptions& options,
    std::string& gauge){
    if(!table.valid)return {space};
    std::vector<double> chars;for(const auto& d:reps)chars.push_back((space.transpose()*d*space).trace());
    const auto content=decompose_point_group_characters(table,chars,options.symmetry_tolerance);
    if(!content.valid)return {space};
    std::size_t row=table.rows.size(),copies=0;
    for(std::size_t i=0;i<content.multiplicities.size();++i)if(content.multiplicities[i]){
        if(row!=table.rows.size())return {space};row=i;copies=content.multiplicities[i];}
    if(row==table.rows.size()||copies<=1)return {space};
    const auto dimension=Eigen::Index(table.rows[row].dimension);
    std::vector<M> blocks{space};
    const auto refine=[&](const M& seed,double cluster_tolerance){
        std::vector<M> next;
        for(const auto& block:blocks){
            if(block.cols()==dimension){next.push_back(block);continue;}
            const M k=block.transpose()*seed*block;
            Eigen::SelfAdjointEigenSolver<M> eig((k+k.transpose())*.5);
            if(eig.info()!=Eigen::Success){next.push_back(block);continue;}
            std::vector<M> pieces;bool valid=true;
            for(Eigen::Index lo=0;lo<k.cols();){Eigen::Index hi=lo+1;
                while(hi<k.cols()&&std::abs(eig.eigenvalues()[hi]-eig.eigenvalues()[lo])<=cluster_tolerance)++hi;
                if((hi-lo)%dimension){valid=false;break;}
                const M projected=block*eig.eigenvectors().middleCols(lo,hi-lo);
                M q=deterministic_basis(projected,1e-9);if(!q.size()){valid=false;break;}
                std::vector<double> signature;for(std::size_t g=0;g<reps.size();++g){const auto& d=reps[g];
                    const M dq=d*q;
                    if(copy_leakage(q,d,outside_metrics[g])>options.symmetry_tolerance){valid=false;break;}
                    signature.push_back((q.transpose()*dq).trace());}
                if(!valid)break;
                const auto decomposition=decompose_point_group_characters(table,signature,options.symmetry_tolerance);
                if(!decomposition.valid||decomposition.multiplicities[row]!=std::size_t((hi-lo)/dimension)){valid=false;break;}
                pieces.push_back(std::move(q));lo=hi;
            }
            if(valid)next.insert(next.end(),pieces.begin(),pieces.end());else next.push_back(block);
        }
        blocks=std::move(next);
    };
    if(fock.size()){
        double commutator=0;for(const auto& d:reps)commutator=std::max(commutator,err(d.transpose()*fock*d-fock));
        if(commutator<=options.energy_tolerance_hartree){refine(fock,options.energy_tolerance_hartree);gauge="physical_fock_spectrum";}
        // A small absolute Fock asymmetry can rotate nearly coincident copies
        // appreciably. The Reynolds mean supplies a mathematical copy gauge;
        // every energy below still uses the original, unmodified operator.
        if(std::any_of(blocks.begin(),blocks.end(),[&](const M& q){return q.cols()!=dimension;})){
            M averaged=M::Zero(fock.rows(),fock.cols());for(const auto& d:reps)averaged+=d.transpose()*fock*d;
            averaged/=double(reps.size());refine(averaged,options.energy_tolerance_hartree);
            gauge="reynolds_fock_copy_gauge; original_operator_expectations";
        }
    }
    // Symmetric matrix units span all real symmetric seeds. Reynolds averaging
    // therefore spans the self-adjoint commutant, including equivalent copies.
    for(Eigen::Index i=0;i<space.rows();++i)for(Eigen::Index j=0;j<=i;++j){
        if(std::all_of(blocks.begin(),blocks.end(),[&](const M& q){return q.cols()==dimension;}))break;
        M seed=M::Zero(space.rows(),space.rows());
        for(const auto& d:reps){seed+=d.row(i).transpose()*d.row(j);if(i!=j)seed+=d.row(j).transpose()*d.row(i);}
        seed/=double(reps.size());refine(seed,options.eigenvalue_cluster_tolerance);
        gauge="reynolds_symmetric_matrix_units_v1; producer_NAO_order";
    }
    if(blocks.size()!=copies||std::any_of(blocks.begin(),blocks.end(),[&](const M& q){return q.cols()!=dimension;})){
        gauge="unresolved_isotypic_multiplicity";return {space};}
    M joined(space.rows(),space.cols());Eigen::Index offset=0;for(const auto& q:blocks){joined.middleCols(offset,q.cols())=q;offset+=q.cols();}
    if(err(joined.transpose()*joined-M::Identity(joined.cols(),joined.cols()))>2e-5||
       err(joined*joined.transpose()-space*space.transpose())>2e-5){gauge="copy_partition_rejected";return {space};}
    return blocks;
}
struct Energy {NboSalcEnergyEvidence evidence;M fock;};
M archive_coefficients(const Wavefunction& w,const NboAssociation& assoc,const M& x){
    const auto n=std::size_t(w.basis_count);
    if(x.rows()!=n||assoc.gaussian_row.size()!=n||assoc.coefficient_scale.size()!=n||
       w.gaussian_ao_transform.size()!=n)return {};
    std::vector<std::size_t> inv(n,n);
    for(std::size_t i=0;i<n;++i)if(w.gaussian_ao_transform[i].source_index<n)
        inv[w.gaussian_ao_transform[i].source_index]=i;
    M y(n,x.cols());
    for(std::size_t k=0;k<n;++k){
        if(assoc.gaussian_row[k]>=n||inv[assoc.gaussian_row[k]]>=n)return {};
        const auto i=inv[assoc.gaussian_row[k]];
        const double scale=assoc.coefficient_scale[k]*w.gaussian_ao_transform[i].coefficient_scale;
        if(!std::isfinite(scale)||std::abs(scale)<1e-20)return {};
        y.row(k)=x.row(i)/scale;
    }
    return y;
}
M verified_density(const Wavefunction& w,const NboIntegration& data,NboSpin spin,
                   const M& a,double tolerance,std::string& reason){
    const auto& d=data.dataset;const auto n=std::size_t(w.basis_count);
    const NboCanonicalEvidence* evidence=nullptr;
    for(const auto& item:d.association.canonical_evidence)if(item.spin==spin){
        if(evidence){reason="Ambiguous same-spin density association";return {};}
        evidence=&item;
    }
    if(!d.archive||!evidence||!evidence->density_verified){
        reason="No independently associated same-spin producer density";return {};
    }
    const auto* p=find(d.archive->matrices,"DENSITY",spin);
    const auto* s=find(d.archive->matrices,"OVERLAP",NboSpin::Total);
    if(!complete(p,n)||!complete(s,n)){
        reason="Missing complete same-spin producer density or overlap";return {};
    }
    const M density=matrix(*p);
    if(err(density-density.transpose())>tolerance){
        reason="Producer density is not Hermitian within tolerance";return {};
    }
    M b=archive_coefficients(w,d.association,a);
    if(!b.size()){reason="Missing or invalid AO convention mapping for density";return {};}
    if(d.archive->density_is_bond_order)b=(matrix(*s)*b).eval();
    reason="Independently associated same-spin archive density projected in the NAO metric";
    return b.transpose()*density*b;
}
Energy check_canonical_energy(const Wavefunction& w,const NboIntegration& data,NboSpin spin,const M& a,const NboSalcOptions& options){
    Energy e;e.evidence.spin=spin;const auto& d=data.dataset;const auto n=std::size_t(w.basis_count);if(!d.archive){e.evidence.detail="No same-source archive";return e;}
    e.evidence.input_units=d.archive->fock_input_units;const auto* f=find(d.archive->matrices,"FOCK",spin);const auto* s=find(d.archive->matrices,"OVERLAP",NboSpin::Total);
    if(!complete(f,n)){e.evidence.detail="Source does not provide a complete same-spin Fock matrix; side orbital wavefunctions remain available";return e;}
    if(!complete(s,n)){e.evidence.detail="Source does not provide a complete overlap matrix";return e;}e.evidence.source=f->source;
    const auto& assoc=d.association;const NboCanonicalEvidence* ev=nullptr;for(const auto& x:assoc.canonical_evidence)if(x.spin==spin)ev=&x;
    if(!ev||!ev->direct_fchk_coefficients){e.evidence.detail="No direct same-spin canonical coefficient identity";return e;}
    if(assoc.gaussian_row.size()!=n||assoc.coefficient_scale.size()!=n||w.gaussian_ao_transform.size()!=n){e.evidence.detail="Missing AO convention mapping";return e;}
    std::vector<std::size_t> indices;for(std::size_t i=0;i<w.orbitals.size();++i)if((spin==NboSpin::Beta)==(w.orbitals[i].spin==Spin::Beta))indices.push_back(i);
    if(indices.empty()){e.evidence.detail="No canonical columns for spin";return e;}
    M c(n,indices.size());V eps(indices.size());for(std::size_t j=0;j<indices.size();++j){const auto& mo=w.orbitals[indices[j]];if(mo.coefficients.size()!=n||!std::isfinite(mo.energy_hartree)){e.evidence.detail="Invalid canonical energy or coefficient dimensions";return e;}c.col(j)=Eigen::Map<const V>(mo.coefficients.data(),n);eps[j]=mo.energy_hartree;}
    c=archive_coefficients(w,assoc,c);const M an=archive_coefficients(w,assoc,a);if(c.size()==0||an.size()==0){e.evidence.detail="Invalid AO convention scale";return e;}
    const M fock=matrix(*f),overlap=matrix(*s);e.evidence.hermiticity_error=err(fock-fock.transpose());const M fc=fock*c,sc=overlap*c,residual=fc-sc*eps.asDiagonal();
    e.evidence.canonical_residual=err(residual);e.evidence.projected_residual=err(c.transpose()*residual);e.evidence.eigenvalue_error_hartree=err(c.transpose()*fc-eps.asDiagonal().toDenseMatrix());e.evidence.canonical_columns_checked=indices.size();
    Eigen::SelfAdjointEigenSolver<M> se(overlap);if(se.info()!=Eigen::Success){e.evidence.detail="Overlap spectral decomposition failed";return e;}
    // Numerical overlap rank is not the producer's retained canonical rank.
    // In particular, Gaussian can discard finite positive overlap directions.
    e.evidence.overlap_rank_threshold=std::max(1e-12,se.eigenvalues().maxCoeff()*1e-9);
    for(Eigen::Index i=0;i<se.eigenvalues().size();++i)if(se.eigenvalues()[i]>e.evidence.overlap_rank_threshold)++e.evidence.overlap_numerical_rank;
    const M gram=c.transpose()*sc;
    if(err(gram-M::Identity(c.cols(),c.cols()))>options.metric_tolerance){e.evidence.status="rejected_canonical_metric";e.evidence.detail="Canonical retained columns are not an orthonormal independent effective space";return e;}
    e.evidence.canonical_effective_rank=indices.size();e.evidence.effective_rank=indices.size();
    e.evidence.canonical_null_directions=n-indices.size();
    const auto solver=gram.ldlt();const M external=an-c*solver.solve(c.transpose()*overlap*an);
    const M external_metric=external.transpose()*overlap*external;
    e.evidence.outside_canonical_nao_norm=std::sqrt(std::max(0.0,external_metric.diagonal().maxCoeff()));
    e.evidence.outside_canonical_fock_coupling=err(c.transpose()*fock*external);
    e.evidence.nullspace_residual=err(residual-sc*solver.solve(c.transpose()*residual));
    const bool full_ok=e.evidence.canonical_residual<=options.energy_tolerance_hartree;
    if(e.evidence.canonical_null_directions&&e.evidence.outside_canonical_nao_norm>options.metric_tolerance){e.evidence.status="unverified_operator_outside_canonical";e.evidence.detail="NAO space extends beyond validated canonical effective space; exterior Fock action cannot be established from retained canonical eigenpairs; side energies withheld";return e;}
    if(e.evidence.hermiticity_error>options.energy_tolerance_hartree||e.evidence.projected_residual>options.energy_tolerance_hartree||e.evidence.eigenvalue_error_hartree>options.energy_tolerance_hartree||!full_ok){e.evidence.status="rejected_operator";e.evidence.detail=full_ok?"Fock/canonical projected operator mismatch":"Full FC-SC epsilon residual fails; projected/null-space diagnostics retained, side energies withheld";return e;}
    e.fock=an.transpose()*fock*an;e.evidence.available=true;e.evidence.canonical_same_operator=true;e.evidence.status="verified_same_operator";e.evidence.detail="Complete Fock, Hermiticity, full and projected canonical equations agree; not isolated-fragment energies";return e;
}

Energy check_energy(const Wavefunction& w,const NboIntegration& data,NboSpin spin,
                    const M& a,const NboSalcOptions& options) {
    auto e=check_canonical_energy(w,data,spin,a,options);
    e.evidence.canonical_operator_status=e.evidence.status;
    if(e.evidence.available)return e;
    const auto& d=data.dataset;
    // RO canonical eigenvalues describe an effective shared-spatial operator;
    // neither physical spin Fock matrix must diagonalize those orbitals.
    // Accept a distinct physical operator only with independently matched
    // density, complete local transforms and two printed diagonal checks.
    if(!d.archive || !d.archive->open_shell || spin==NboSpin::Total ||
       !d.association.compatible || e.evidence.status=="rejected_canonical_metric")return e;
    const auto n=static_cast<std::size_t>(w.basis_count);
    const auto* f=find(d.archive->matrices,"FOCK",spin);
    const auto* s=find(d.archive->matrices,"OVERLAP",NboSpin::Total);
    const auto* b=find(d.matrices,"AONBO",spin);
    if(!complete(f,n)||!complete(s,n)||!complete(b,n)||a.cols()!=Eigen::Index(n))return e;
    std::string density_reason;
    const M density=verified_density(w,data,spin,a,options.metric_tolerance,density_reason);
    if(density.rows()!=a.cols())return e;
    const M an=archive_coefficients(w,d.association,a);
    if(an.rows()!=Eigen::Index(n))return e;
    const M fock=matrix(*f),overlap=matrix(*s),bn=matrix(*b);
    e.evidence.source=f->source;
    e.evidence.hermiticity_error=err(fock-fock.transpose());
    if(e.evidence.hermiticity_error>options.energy_tolerance_hartree ||
       err(bn.transpose()*overlap*bn-M::Identity(n,n))>options.metric_tolerance)return e;
    const M fn=an.transpose()*fock*an;
    const M fb=bn.transpose()*fock*bn;
    // Producer NAO and NBO energy tables print five decimal places. This is
    // an independent print-precision gate, not a relaxed canonical residual.
    constexpr double printed_tolerance=5.1e-6;
    std::vector<bool> nao_seen(n,false),nbo_seen(n,false);
    bool valid=true;
    for(const auto& row:d.naos)if(row.spin==spin) {
        if(!row.id || row.id>n || nao_seen[row.id-1] || !row.energy_hartree ||
           !std::isfinite(*row.energy_hartree)){valid=false;continue;}
        const auto i=row.id-1;nao_seen[i]=true;++e.evidence.printed_nao_checked;
        e.evidence.printed_nao_error_hartree=std::max(e.evidence.printed_nao_error_hartree,
            std::abs(fn(i,i)-*row.energy_hartree));
        if(std::abs(density(i,i)-row.occupation)>4e-5)valid=false;
    }
    for(const auto& row:d.orbitals)if(row.spin==spin) {
        if(!row.id || row.id>n || nbo_seen[row.id-1] || !row.diagonal_fock_hartree ||
           !std::isfinite(*row.diagonal_fock_hartree)){valid=false;continue;}
        const auto i=row.id-1;nbo_seen[i]=true;++e.evidence.printed_nbo_checked;
        e.evidence.printed_nbo_error_hartree=std::max(e.evidence.printed_nbo_error_hartree,
            std::abs(fb(i,i)-*row.diagonal_fock_hartree));
    }
    // Printed off-diagonal magnitudes are additional rejection evidence, not
    // a claim that a thresholded table uniquely determines the full matrix.
    for(const auto& row:d.e2)if(row.spin==spin) {
        if(!row.donor || !row.acceptor || row.donor>n || row.acceptor>n ||
           !std::isfinite(row.fock_hartree) || !std::isfinite(row.energy_gap_hartree)){
            valid=false;continue;
        }
        const auto i=row.donor-1,j=row.acceptor-1;
        ++e.evidence.printed_couplings_checked;
        e.evidence.printed_coupling_error_hartree=std::max(e.evidence.printed_coupling_error_hartree,
            std::abs(std::abs(fb(i,j))-std::abs(row.fock_hartree)));
        e.evidence.printed_gap_error_hartree=std::max(e.evidence.printed_gap_error_hartree,
            std::abs(fb(j,j)-fb(i,i)-row.energy_gap_hartree));
    }
    if(!valid || !std::all_of(nao_seen.begin(),nao_seen.end(),[](bool v){return v;}) ||
       !std::all_of(nbo_seen.begin(),nbo_seen.end(),[](bool v){return v;})) {
        e.evidence.detail+="; independent physical spin-Fock verification lacks complete matched printed local energies/density";
        return e;
    }
    if(e.evidence.printed_nao_error_hartree>printed_tolerance ||
       e.evidence.printed_nbo_error_hartree>printed_tolerance ||
       e.evidence.printed_coupling_error_hartree>0.000501 ||
       e.evidence.printed_gap_error_hartree>0.00501) {
        e.evidence.status="rejected_printed_operator";
        e.evidence.detail+="; physical Fock disagrees with the matched printed NAO/NBO energy or coupling evidence";
        return e;
    }
    e.fock=fn;e.evidence.available=true;e.evidence.printed_operator_verified=true;
    e.evidence.status="verified_physical_spin_fock";
    e.evidence.detail="Associated physical spin-Fock expectations verified against complete NAO/NBO printed energies and available couplings; distinct from the effective canonical operator";
    return e;
}
struct Family {std::size_t fragment=0;std::string type,angular;std::vector<std::size_t> indices;};
std::string quoted(const std::string& s){std::ostringstream o;o<<'"';for(unsigned char c:s){if(c=='"'||c=='\\')o<<'\\'<<char(c);else if(c=='\n')o<<"\\n";else if(c=='\r')o<<"\\r";else if(c=='\t')o<<"\\t";else if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else o<<char(c);}o<<'"';return o.str();}
template<class T,class F>void json_array(std::ostream& o,const std::vector<T>& a,F fn){o<<'[';bool first=true;for(const auto& x:a){if(!first)o<<',';first=false;fn(x);}o<<']';}
void num(std::ostream& o,double x){if(std::isfinite(x))o<<x;else o<<"null";}
void optional(std::ostream& o,const std::optional<double>& x){if(x)num(o,*x);else o<<"null";}
void energy_provenance(std::ostream& o,const NboSalcEnergyEvidence& x) {
    o<<"\"canonical_same_operator\":"<<(x.canonical_same_operator?"true":"false")
     <<",\"printed_operator_verified\":"<<(x.printed_operator_verified?"true":"false")
     <<",\"canonical_operator_status\":"<<quoted(x.canonical_operator_status)
     <<",\"printed_nao_checked\":"<<x.printed_nao_checked
     <<",\"printed_nbo_checked\":"<<x.printed_nbo_checked
     <<",\"printed_couplings_checked\":"<<x.printed_couplings_checked
     <<",\"printed_nao_error_hartree\":";num(o,x.printed_nao_error_hartree);
    o<<",\"printed_nbo_error_hartree\":";num(o,x.printed_nbo_error_hartree);
    o<<",\"printed_coupling_error_hartree\":";num(o,x.printed_coupling_error_hartree);
    o<<",\"printed_gap_error_hartree\":";num(o,x.printed_gap_error_hartree);o<<',';
}
} // namespace

NboElectronicSymmetryScope analyse_nbo_electronic_symmetry_scope(const Wavefunction& w,const NboSalcOptions& options){
    const auto geometry=analyse_molecular_symmetry(w);const auto group=make_group(w,geometry,options);
    return density_scope(w,geometry,group,options);
}
NboSalcModel build_nbo_salc_model(const Wavefunction& w,const NboIntegration& data,const NboSalcOptions& options){
    NboSalcModel out;out.dataset_id=data.id;out.canonical_fingerprint=nbo_canonical_fingerprint(w);std::ostringstream key;key<<data.id<<':'<<out.canonical_fingerprint<<':'<<std::setprecision(17)<<options.metric_tolerance<<':'<<options.symmetry_tolerance<<':'<<options.energy_tolerance_hartree<<':'<<options.eigenvalue_cluster_tolerance<<':'<<options.maximum_group_order;out.cache_key=key.str()+"|salc-copy-v2";
    const auto n=std::size_t(w.basis_count);if(!n||w.ao_overlap.size()!=n*n||!data.dataset.association.compatible||data.canonical_fingerprint!=out.canonical_fingerprint){out.detail="Missing metric, rejected association, or changed immutable canonical identity";return out;}
    const M s=Eigen::Map<const RM>(w.ao_overlap.data(),n,n);const auto geometry=analyse_molecular_symmetry(w);out.point_group=geometry.point_group;auto group=make_group(w,geometry,options);out.group_verified=group.valid;out.group_closure_error=group.error;out.used_group=group.valid?(geometry.linear?"finite sampling subgroup of "+geometry.point_group:geometry.point_group):"unavailable";out.operations=group.valid?group.ops:std::vector<SymmetryOperation>{};if(!group.valid)out.diagnostics.push_back(group.detail);
    out.symmetry_scope=density_scope(w,geometry,group,options);
    MolecularSymmetry full_geometry=geometry;
    if(out.symmetry_scope.reduced){full_geometry.point_group=out.symmetry_scope.naming_group;full_geometry.operations=out.symmetry_scope.operations;
        group=make_group(w,full_geometry,options);out.group_verified=group.valid;out.group_closure_error=group.error;
        out.operations=group.ops;out.used_group=full_geometry.point_group;}
    full_geometry.operations=out.operations;
    const auto irrep_table=group.valid&&!geometry.linear?molecular_point_group_irreps(w,full_geometry):PointGroupIrrepTable{};
    std::vector<std::size_t> atom_group(w.atoms.size());std::iota(atom_group.begin(),atom_group.end(),0);auto root=[&](std::size_t a){while(atom_group[a]!=a)a=atom_group[a];return a;};
    if(group.valid&&w.atoms.size()!=2)for(const auto& op:group.ops)for(std::size_t a=0;a<w.atoms.size();++a){auto x=root(a),y=root(op.atom_permutation[a]);if(x!=y)atom_group[std::max(x,y)]=std::min(x,y);}
    std::map<std::size_t,std::vector<std::size_t>> atomsets;for(std::size_t a=0;a<w.atoms.size();++a)atomsets[root(a)].push_back(a);
    std::vector<std::vector<std::size_t>> ordered;for(const auto& p:atomsets)ordered.push_back(p.second);std::stable_sort(ordered.begin(),ordered.end(),[&](const auto& a,const auto& b){if(a.size()==1&&b.size()!=1)return true;if(a.size()!=1&&b.size()==1)return false;return w.atoms[a[0]].atomic_number>w.atoms[b[0]].atomic_number;});
    for(const auto& a:ordered){NboSalcFragment f;f.id="fragment"+key_atoms(a);f.atoms=a;f.label=atoms_label(w,a);f.side=out.fragments.empty()?0:1;for(auto atom:a)atom_group[atom]=out.fragments.size();out.fragments.push_back(f);}
    for(auto spin:{NboSpin::Total,NboSpin::Alpha,NboSpin::Beta}){
        std::vector<const NboOrbitalDescriptor*> descriptors;for(const auto& x:data.orbitals)if(x.ref.kind==NboOrbitalKind::NAO&&x.ref.spin==spin&&x.orthonormal_basis&&x.coefficients.size()==n&&x.atoms.size()==1)descriptors.push_back(&x);if(descriptors.empty())continue;
        std::sort(descriptors.begin(),descriptors.end(),[](auto a,auto b){return a->ref.index<b->ref.index;});const auto r=descriptors.size();M a(n,r);for(std::size_t j=0;j<r;++j)a.col(j)=Eigen::Map<const V>(descriptors[j]->coefficients.data(),n);const M sa=s*a;const double orth=err(a.transpose()*sa-M::Identity(r,r));out.orthogonality_error=std::max(out.orthogonality_error,orth);if(orth>options.metric_tolerance){out.diagnostics.push_back(std::string(nbo_spin_name(spin))+": NAO metric orthogonality rejected");continue;}
        auto energy=check_energy(w,data,spin,a,options);
        std::vector<M> representations,transformed_basis;bool transforms=group.valid;
        // Both metric isometry and family closure use the complete group, not
        // a coordinate-dependent choice of generators. Approximately invariant
        // fields can pass one generator while failing its powers.
        if(transforms)for(const auto& op:group.ops){const auto t=apply_orbital_symmetry_operation(w,op,flat(a),r);if(t.size()!=n*r){transforms=false;out.diagnostics.push_back("Full-group AO action unavailable; fixed NAO fallback");break;}const M ta=Eigen::Map<const RM>(t.data(),n,r);const double metric_error=err(ta.transpose()*s*ta-M::Identity(r,r));out.representation_error=std::max(out.representation_error,metric_error);if(metric_error>options.symmetry_tolerance){transforms=false;out.diagnostics.push_back("AO action fails metric isometry; fixed NAO fallback");break;}transformed_basis.push_back(ta);}
        if(transforms)for(const auto& action:transformed_basis)representations.push_back(sa.transpose()*action);
        if(energy.evidence.available&&transforms){double e=0;for(const auto& g:representations)e=std::max(e,err(g.transpose()*energy.fock*g-energy.fock));energy.evidence.fock_symmetry_error=e;energy.evidence.electronic_symmetry_verified=e<=options.energy_tolerance_hartree;}
        std::string density_reason;
        const M density=verified_density(w,data,spin,a,options.metric_tolerance,density_reason);
        NboSalcSpinOperator op;op.spin=spin;
        for(const auto* descriptor:descriptors)op.basis.push_back(descriptor->ref);
        if(energy.evidence.available)op.fock=flat(energy.fock);
        if(density.size())op.density=flat(density);
        op.energy_status=energy.evidence.detail;op.density_status=density_reason;
        out.spin_operators.push_back(std::move(op));
        out.diagnostics.push_back(std::string(nbo_spin_name(spin))+": SALC occupation: "+density_reason);
        // Canonical links and producer-density occupations have independent gates.
        // In particular, verified RO beta density needs no invented beta MO block.
        std::vector<std::size_t> canonical_indices;for(std::size_t j=0;j<w.orbitals.size();++j)if((spin==NboSpin::Beta)==(w.orbitals[j].spin==Spin::Beta))canonical_indices.push_back(j);
        M projections(r,canonical_indices.size());bool canonical_verified=false;for(const auto& ev:data.dataset.association.canonical_evidence)if(ev.spin==spin&&ev.direct_fchk_coefficients)canonical_verified=true;
        if(canonical_verified){M c(n,canonical_indices.size());for(std::size_t j=0;j<canonical_indices.size();++j)c.col(j)=Eigen::Map<const V>(w.orbitals[canonical_indices[j]].coefficients.data(),n);projections=sa.transpose()*c;}
        if(density.size()&&transforms){energy.evidence.density_symmetry_checked=true;for(const auto& g:representations)energy.evidence.density_symmetry_error=std::max(energy.evidence.density_symmetry_error,err(g.transpose()*density*g-density));energy.evidence.electronic_symmetry_verified=energy.evidence.electronic_symmetry_verified&&energy.evidence.density_symmetry_error<=options.symmetry_tolerance;}else energy.evidence.electronic_symmetry_verified=false;
        out.energies.push_back(energy.evidence);
        std::map<std::pair<std::size_t,std::string>,Family> familymap;for(std::size_t j=0;j<r;++j){const auto atom=descriptors[j]->atoms[0];if(atom>=atom_group.size())continue;const NboNao* row=nullptr;for(const auto& x:data.dataset.naos)if(x.spin==spin&&x.id==descriptors[j]->ref.index+1){row=&x;break;}const std::string type=row?compact(row->type):"unclassified";auto& f=familymap[{atom_group[atom],type}];f.fragment=atom_group[atom];f.type=type;f.angular=row?row->angular:"unclassified";f.indices.push_back(j);}
        const auto start=out.orbitals.size();
        for(const auto& entry:familymap){const auto& family=entry.second;const auto b=family.indices.size();const auto& fragment=out.fragments[family.fragment];std::vector<M> reps,outside_metrics;bool closed=transforms&&group.ops.size()>1;double closure=0;
            M family_basis(n,b);for(std::size_t j=0;j<b;++j)family_basis.col(j)=a.col(family.indices[j]);
            // Evaluate the actual field residual. Subtracting two nearly unit
            // norms assumes exactly orthonormal printed NAOs and turns their
            // small metric error into its square root, causing frame-dependent
            // false failures for ill-conditioned/rank-limited AO bases.
            const M original_family=family_basis;
            const M gram=family_basis.transpose()*s*family_basis;
            Eigen::SelfAdjointEigenSolver<M> metric_solver((gram+gram.transpose())*.5);
            M whitening=M::Identity(b,b);
            if(fragment.atoms.size()>1&&metric_solver.info()==Eigen::Success&&metric_solver.eigenvalues().minCoeff()>1e-10)
                whitening=metric_solver.eigenvectors()*metric_solver.eigenvalues().array().rsqrt().matrix().asDiagonal()*metric_solver.eigenvectors().transpose();
            else if(metric_solver.info()!=Eigen::Success)closed=false;
            family_basis=original_family*whitening;
            const M family_metric=s*family_basis;
            M family_fock,family_density;
            if(energy.evidence.available){M raw(b,b);for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)raw(i,j)=energy.fock(family.indices[i],family.indices[j]);family_fock=whitening*raw*whitening;}
            if(density.size()){M raw(b,b);for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)raw(i,j)=density(family.indices[i],family.indices[j]);family_density=whitening*raw*whitening;}
            if(closed)for(const auto& action:transformed_basis){M transformed(n,b);for(std::size_t j=0;j<b;++j)transformed.col(j)=action.col(family.indices[j]);transformed=(transformed*whitening).eval();const M local=family_metric.transpose()*transformed;reps.push_back(local);const M residual=transformed-family_basis*local;const M residual_metric=residual.transpose()*s*residual;outside_metrics.push_back(residual_metric);for(std::size_t j=0;j<b;++j)closure=std::max(closure,std::sqrt(std::max(0.0,residual_metric(j,j))));}
            if(closure>options.symmetry_tolerance)closed=false;
            if(closed){double relation=0;for(std::size_t i=0;i<reps.size();++i)for(std::size_t j=0;j<reps.size();++j)relation=std::max(relation,err(reps[i]*reps[j]-reps[group.table[i*group.ops.size()+j]]));closure=std::max(closure,relation);if(relation>options.symmetry_tolerance)closed=false;}
            std::vector<M> spaces;
            // Single-atom orbitals keep the literal producer NAO gauge. Never
            // mix O s/p merely because they carry the same one-dimensional irrep.
            if(closed&&fragment.atoms.size()>1){M central=M::Zero(b,b);std::size_t ci=0;for(const auto& cl:group.classes){M z=M::Zero(b,b);for(auto g:cl)z+=reps[g];z=(z+z.transpose()).eval()/(2*double(cl.size()));central+=std::sqrt(double(++ci)+1.6180339887498948)*z;}Eigen::SelfAdjointEigenSolver<M> eig(central);if(eig.info()==Eigen::Success){for(Eigen::Index lo=0;lo<eig.eigenvalues().size();){Eigen::Index hi=lo+1;while(hi<eig.eigenvalues().size()&&std::abs(eig.eigenvalues()[hi]-eig.eigenvalues()[lo])<=options.eigenvalue_cluster_tolerance*std::max(1.0,std::abs(eig.eigenvalues()[lo])))++hi;auto q=deterministic_basis(eig.eigenvectors().middleCols(lo,hi-lo),1e-9);if(!q.size()){spaces.clear();break;}spaces.push_back(std::move(q));lo=hi;}}}
            if(closed&&fragment.atoms.size()>1&&!spaces.empty())spaces=merge_group_connected_spaces(std::move(spaces),reps,options.symmetry_tolerance);
            if(spaces.empty()){spaces.push_back(M::Identity(b,b));if(fragment.atoms.size()>1)closed=false;}
            std::vector<std::string> copy_gauges;
            if(closed&&fragment.atoms.size()>1){std::vector<M> split;
                for(const auto& q:spaces){std::string gauge="stable_projector_columns";
                    auto copies=split_irrep_copies(q,reps,outside_metrics,irrep_table,family_fock,options,gauge);
                    for(auto& copy:copies){split.push_back(std::move(copy));copy_gauges.push_back(gauge);}}
                spaces=std::move(split);}
            // Acceptance is for the entire fixed basis, not just each block.
            // A bad gauge or noninvariant spectral split must never introduce
            // duplicate/nonorthogonal directions and inflate MO coverage.
            M joined(b,b);Eigen::Index offset=0;bool valid_partition=true;
            for(const auto& q:spaces){if(offset+q.cols()>Eigen::Index(b)){valid_partition=false;break;}joined.middleCols(offset,q.cols())=q;offset+=q.cols();if(closed)for(const auto& rep:reps)if(err(rep*q-q*(q.transpose()*rep*q))>options.symmetry_tolerance)valid_partition=false;}
            if(offset!=Eigen::Index(b)||err(joined.transpose()*joined-M::Identity(b,b))>options.metric_tolerance)valid_partition=false;
            if(!valid_partition){spaces={M::Identity(b,b)};closed=false;out.diagnostics.push_back(fragment.id+":"+family.type+": rejected spectral partition; preserving complete fixed NAO family");}
            std::size_t ordinal=0;for(const auto& q:spaces){NboSalcSubspace sub;sub.id="salc-copy-v2:"+out.used_group+":"+fragment.id+":"+nbo_spin_name(spin)+":"+family.type+":space"+std::to_string(++ordinal);sub.fragment_id=fragment.id;sub.spin=spin;sub.dimension=q.cols();sub.closure_error=closure;sub.orthogonality_error=err(q.transpose()*q-M::Identity(q.cols(),q.cols()));sub.symmetry_verified=closed;double retention=0;
                if(closed){for(std::size_t g=0;g<reps.size();++g){const auto& rep=reps[g];const M rq=rep*q;retention=std::max(retention,copy_leakage(q,rep,outside_metrics[g]));double character=(q.transpose()*rq).trace();sub.characters.push_back(character);sub.character_norm+=character*character/double(reps.size());}sub.closure_error=std::max(sub.closure_error,retention);if(retention>options.symmetry_tolerance)sub.symmetry_verified=false;const auto mult=std::size_t(std::llround(std::sqrt(sub.character_norm)));if(mult&&std::abs(sub.character_norm-double(mult*mult))<1e-3&&sub.dimension%mult==0){sub.multiplicity=mult;sub.irrep_dimension=sub.dimension/mult;}}
                if(sub.symmetry_verified&&irrep_table.valid){
                    const auto decomposition=decompose_point_group_characters(irrep_table,sub.characters,options.symmetry_tolerance);
                    sub.multiplicity=0;sub.irrep_dimension=0;
                    if(decomposition.valid){std::size_t row_count=0,row_index=0;for(std::size_t i=0;i<decomposition.multiplicities.size();++i)if(decomposition.multiplicities[i]){++row_count;row_index=i;}
                        if(row_count==1){sub.multiplicity=decomposition.multiplicities[row_index];sub.irrep_dimension=irrep_table.rows[row_index].dimension;}}
                }
                sub.copy_gauge=ordinal<=copy_gauges.size()?copy_gauges[ordinal-1]:"literal_producer_NAO";
                sub.basis_kind=fragment.atoms.size()>1?"derived_salc":"producer_nao";
                const M mapped=whitening*q;
                for(auto index:family.indices)sub.source_basis.push_back(descriptors[index]->ref);
                sub.source_to_derived=flat(mapped);
                if(family_density.size())sub.density=flat(q.transpose()*family_density*q);
                if(family_fock.size()){
                    const M k=q.transpose()*family_fock*q;sub.fock=flat(k);
                    Eigen::SelfAdjointEigenSolver<M> spectrum((k+k.transpose())*.5);
                    sub.energy_scalar_residual_hartree=err(k-M::Identity(k.rows(),k.cols())*(k.trace()/double(k.rows())));
                    sub.energy_spectral_width_hartree=spectrum.info()==Eigen::Success?spectrum.eigenvalues().maxCoeff()-spectrum.eigenvalues().minCoeff():std::numeric_limits<double>::infinity();
                    for(const auto& d:reps){const M local=q.transpose()*d*q;sub.energy_commutator_hartree=std::max(sub.energy_commutator_hartree,err(local.transpose()*k*local-k));}
                    sub.energy_copy_coupling_hartree=err(family_fock*q-q*k);
                    sub.energy_degeneracy_verified=sub.symmetry_verified&&sub.multiplicity==1&&
                        sub.energy_scalar_residual_hartree<=options.energy_tolerance_hartree&&sub.energy_spectral_width_hartree<=options.energy_tolerance_hartree&&
                        sub.energy_commutator_hartree<=options.energy_tolerance_hartree&&sub.energy_copy_coupling_hartree<=options.energy_tolerance_hartree;
                    sub.energy_degeneracy_status=sub.energy_degeneracy_verified?"verified_same_operator_degeneracy":"actual_operator_split_or_coupled";
                }
                sub.label=sub.symmetry_verified?"Symmetry channel "+std::to_string(ordinal):"Fixed NAO subspace";sub.detail=sub.symmetry_verified?"Geometry-verified invariant subspace; character/partner evidence retained. Phase gauge: stable projector columns in producer NAO order.":"No asserted SALC: trivial group, unsupported action, or literal family is not closed; fixed NAOs retained.";
                for(Eigen::Index col=0;col<q.cols();++col){NboSalcOrbital o;o.id=sub.id+":member"+std::to_string(col);o.fragment_id=fragment.id;o.subspace_id=sub.id;o.spin=spin;o.type=family.type;o.angular=family.angular;o.partner_index=col;o.partner_dimension=q.cols();o.atoms=fragment.atoms;o.symmetry_adapted=sub.symmetry_verified&&fragment.atoms.size()>1;o.detail=sub.detail;
                    V global=V::Zero(r);for(std::size_t row=0;row<b;++row){const double coeff=mapped(row,col);global[family.indices[row]]=coeff;if(coeff!=0)o.terms.push_back({descriptors[family.indices[row]]->ref,coeff});}
                    o.label=o.terms.size()==1?nbo_orbital(data,o.terms[0].orbital)->label:"SALC "+std::to_string(ordinal)+"."+std::to_string(col+1)+" "+family.type;
                    if(o.terms.size()==1){for(const auto& row:data.dataset.naos)if(row.spin==o.terms[0].orbital.spin&&row.id==o.terms[0].orbital.index+1){o.angular=row.angular;break;}}
                    else {o.angular="mixed directions";const auto begin=family.type.find('(');if(begin!=std::string::npos)for(std::size_t k=begin+1;k<family.type.size();++k)if(std::string("spdfghik").find(family.type[k])!=std::string::npos){o.angular=family.type.substr(k,1);break;}}
                    if(fragment.atoms.size()==1)o.partner_dimension=1;
                    if(o.symmetry_adapted&&b==2&&q.cols()==1&&family.angular=="s"&&o.terms.size()==2)o.label=(o.terms[0].coefficient*o.terms[1].coefficient>0?"In-phase ":"Out-of-phase ")+family.type;
                    if(!o.symmetry_adapted&&o.terms.size()>1)o.label="Fragment combination "+family.type;
                    if(energy.evidence.available){o.energy_hartree=(global.transpose()*energy.fock*global)(0,0);o.energy_semantics=energy.evidence.canonical_same_operator?"molecular-environment Fock/KS expectation; not a canonical or isolated-fragment eigenvalue":"physical spin-Fock expectation; distinct from the effective canonical MO operator";}
                    else o.detail+="; energy: "+energy.evidence.detail;
                    if(density.size())o.occupation=(global.transpose()*density*global)(0,0);else if(o.terms.size()==1){const auto* descriptor=nbo_orbital(data,o.terms[0].orbital);if(descriptor)o.occupation=descriptor->occupation;}
                    const auto index=out.orbitals.size();sub.orbital_indices.push_back(index);out.orbitals.push_back(std::move(o));if(canonical_verified)for(std::size_t j=0;j<canonical_indices.size();++j){const double coefficient=global.dot(projections.col(j));out.links.push_back({index,canonical_indices[j],coefficient,coefficient*coefficient});}
                }out.subspaces.push_back(std::move(sub));
            }
        }
        if(canonical_verified){std::map<std::size_t,double> shown;for(const auto& link:out.links)if(link.side_index>=start)shown[link.canonical_index]+=link.weight;for(std::size_t j=0;j<canonical_indices.size();++j){NboSalcCoverage c;c.canonical_index=canonical_indices[j];c.weight_sum=shown[c.canonical_index];const V residual=Eigen::Map<const V>(w.orbitals[c.canonical_index].coefficients.data(),n)-a*projections.col(j);c.residual_norm=std::sqrt(std::max(0.0,(residual.transpose()*s*residual)(0,0)));c.available=true;out.coverage.push_back(c);}}
    }
    out.restricted_open_shell=verify_nbo_restricted_open_shell(w,data);out.available=!out.orbitals.empty();out.status=out.available?"fixed_local_basis_available":"unavailable";out.detail="Fixed NAO/SALC basis; all canonical signed coefficients and raw weights retained without display renormalization. A fragment is a symmetry-related atom family, not an asserted chemical bond.";return out;
}
NboOrbitalSelection nbo_salc_selection(const NboSalcModel& model,std::size_t index){NboOrbitalSelection s;s.dataset_id=model.dataset_id;s.mode=NboSelectionMode::Combination;s.semantic_kind="salc";if(index<model.orbitals.size()){const auto& orbital=model.orbitals[index];s.terms=orbital.terms;s.label=orbital.label;s.group_id=orbital.subspace_id;s.source_id=orbital.id;s.spatial_spin=orbital.spatial_spin;if(s.spatial_spin)s.semantic_kind="spin_averaged_spatial_orbital";}return s;}
NboOrbitalSelection nbo_salc_component_selection(const NboSalcModel& model,const NboSalcLink& link){auto s=nbo_salc_selection(model,link.side_index);s.semantic_kind="salc_component";s.target_canonical_index=link.canonical_index;for(auto& term:s.terms)term.coefficient*=link.coefficient;s.label+=" component of canonical MO "+std::to_string(link.canonical_index+1);return s;}

std::string serialize_nbo_salc_json(const NboSalcModel& m){std::ostringstream o;o<<std::setprecision(17);o<<"{\"dataset_id\":"<<quoted(m.dataset_id)<<",\"canonical_fingerprint\":"<<quoted(m.canonical_fingerprint)<<",\"cache_key\":"<<quoted(m.cache_key)<<",\"available\":"<<(m.available?"true":"false")<<",\"status\":"<<quoted(m.status)<<",\"detail\":"<<quoted(m.detail)<<",\"point_group\":"<<quoted(m.point_group)<<",\"used_group\":"<<quoted(m.used_group)<<",\"group_verified\":"<<(m.group_verified?"true":"false")<<",\"group_closure_error\":";num(o,m.group_closure_error);o<<",\"representation_error\":";num(o,m.representation_error);o<<",\"orthogonality_error\":";num(o,m.orthogonality_error);
    const auto& scope=m.symmetry_scope;
    o<<",\"symmetry_scope\":{\"geometry_group\":"<<quoted(scope.geometry_group)<<",\"naming_group\":"<<quoted(scope.naming_group)<<",\"status\":"<<quoted(scope.status)
     <<",\"naming_scope_verified\":"<<(scope.naming_scope_verified?"true":"false")<<",\"density_checked\":"<<(scope.density_checked?"true":"false")<<",\"reduced\":"<<(scope.reduced?"true":"false")<<",\"relative_tolerance\":"<<scope.relative_tolerance;
    o<<",\"total_density_residuals\":";json_array(o,scope.total_density_residuals,[&](auto v){num(o,v);});
    o<<",\"spin_density_residuals\":";json_array(o,scope.spin_density_residuals,[&](auto v){num(o,v);});
    o<<",\"retained_operation_indices\":";json_array(o,scope.retained_operation_indices,[&](auto v){o<<v;});
    o<<",\"geometry_operations\":";json_array(o,scope.geometry_operations,[&](const auto& op){o<<"{\"matrix\":";std::vector<double> values(op.matrix.begin(),op.matrix.end());json_array(o,values,[&](auto v){num(o,v);});o<<",\"atom_permutation\":";json_array(o,op.atom_permutation,[&](auto v){o<<v;});o<<'}';});o<<'}';
    o<<",\"spin_averaged\":"<<(m.spin_averaged?"true":"false")<<",\"merged_spatial_count\":"<<m.merged_spatial_count<<",\"separate_spin_count\":"<<m.separate_spin_count;
    const auto& ro=m.restricted_open_shell;
    o<<",\"restricted_open_shell\":{\"verified\":"<<(ro.verified?"true":"false")<<",\"status\":"<<quoted(ro.status)<<",\"detail\":"<<quoted(ro.detail)<<",\"method\":"<<quoted(ro.method)<<",\"alpha_electrons\":"<<ro.alpha_electrons<<",\"beta_electrons\":"<<ro.beta_electrons<<",\"shared_columns\":"<<ro.shared_columns<<",\"coefficient_residual\":"<<ro.coefficient_residual<<",\"occupation_error\":"<<ro.occupation_error<<",\"density_error\":"<<ro.density_error<<",\"tolerance\":"<<ro.tolerance<<'}';
    o<<",\"phase_convention\":\"stable projector columns in producer NAO order; largest absolute coefficient positive; individual signs are gauge-dependent\",\"operations\":";json_array(o,m.operations,[&](const auto& x){o<<"{\"matrix\":[";for(int i=0;i<9;++i){if(i)o<<',';num(o,x.matrix[i]);}o<<"],\"atom_permutation\":";json_array(o,x.atom_permutation,[&](auto i){o<<i;});o<<",\"mapping_error_bohr\":";num(o,x.max_mapping_error_bohr);o<<'}';});
    o<<",\"fragments\":";json_array(o,m.fragments,[&](const auto& x){o<<"{\"id\":"<<quoted(x.id)<<",\"label\":"<<quoted(x.label)<<",\"side\":"<<x.side<<",\"atoms\":";json_array(o,x.atoms,[&](auto i){o<<i;});o<<'}';});
    o<<",\"subspaces\":";json_array(o,m.subspaces,[&](const auto& x){o<<"{\"id\":"<<quoted(x.id)<<",\"label\":"<<quoted(x.label)<<",\"fragment_id\":"<<quoted(x.fragment_id)<<",\"spin\":"<<quoted(nbo_spin_name(x.spin))<<",\"detail\":"<<quoted(x.detail)<<",\"dimension\":"<<x.dimension<<",\"irrep_dimension\":"<<x.irrep_dimension<<",\"multiplicity\":"<<x.multiplicity<<",\"symmetry_verified\":"<<(x.symmetry_verified?"true":"false")<<",\"closure_error\":";num(o,x.closure_error);o<<",\"orthogonality_error\":";num(o,x.orthogonality_error);o<<",\"character_norm\":";num(o,x.character_norm);o<<",\"characters\":";json_array(o,x.characters,[&](auto v){num(o,v);});o<<",\"copy_gauge\":"<<quoted(x.copy_gauge)<<",\"basis_kind\":"<<quoted(x.basis_kind)
        <<",\"energy_degeneracy_verified\":"<<(x.energy_degeneracy_verified?"true":"false")<<",\"energy_degeneracy_status\":"<<quoted(x.energy_degeneracy_status)
        <<",\"energy_scalar_residual_hartree\":";num(o,x.energy_scalar_residual_hartree);
        o<<",\"energy_spectral_width_hartree\":";num(o,x.energy_spectral_width_hartree);
        o<<",\"energy_commutator_hartree\":";num(o,x.energy_commutator_hartree);
        o<<",\"energy_copy_coupling_hartree\":";num(o,x.energy_copy_coupling_hartree);
        o<<",\"source_basis\":";json_array(o,x.source_basis,[&](const auto& ref){o<<"{\"kind\":"<<quoted(nbo_orbital_kind_name(ref.kind))<<",\"spin\":"<<quoted(nbo_spin_name(ref.spin))<<",\"index\":"<<ref.index<<'}';});
        o<<",\"source_to_derived\":";json_array(o,x.source_to_derived,[&](auto v){num(o,v);});
        o<<",\"fock\":";json_array(o,x.fock,[&](auto v){num(o,v);});
        o<<",\"density\":";json_array(o,x.density,[&](auto v){num(o,v);});
        o<<",\"orbital_indices\":";json_array(o,x.orbital_indices,[&](auto v){o<<v;});o<<'}';});
    o<<",\"orbitals\":";json_array(o,m.orbitals,[&](const auto& x){o<<"{\"id\":"<<quoted(x.id)<<",\"label\":"<<quoted(x.label)<<",\"fragment_id\":"<<quoted(x.fragment_id)<<",\"subspace_id\":"<<quoted(x.subspace_id)<<",\"type\":"<<quoted(x.type)<<",\"angular\":"<<quoted(x.angular)<<",\"detail\":"<<quoted(x.detail)<<",\"spin\":"<<quoted(nbo_spin_name(x.spin))<<",\"symmetry_adapted\":"<<(x.symmetry_adapted?"true":"false")<<",\"partner_index\":"<<x.partner_index<<",\"partner_dimension\":"<<x.partner_dimension<<",\"energy_hartree\":";optional(o,x.energy_hartree);o<<",\"occupation\":";optional(o,x.occupation);o<<",\"energy_semantics\":"<<quoted(x.energy_semantics)<<",\"spin_correspondence_status\":"<<quoted(x.spin_correspondence_status)<<",\"spatial_spin\":"<<(x.spatial_spin?serialize_nbo_spatial_spin_json(*x.spatial_spin):"null")<<",\"atoms\":";json_array(o,x.atoms,[&](auto i){o<<i;});o<<",\"terms\":";json_array(o,x.terms,[&](const auto& t){o<<"{\"kind\":"<<quoted(nbo_orbital_kind_name(t.orbital.kind))<<",\"spin\":"<<quoted(nbo_spin_name(t.orbital.spin))<<",\"index\":"<<t.orbital.index<<",\"coefficient\":";num(o,t.coefficient);o<<'}';});o<<'}';});
    o<<",\"links\":";json_array(o,m.links,[&](const auto& x){o<<"{\"side_index\":"<<x.side_index<<",\"canonical_index\":"<<x.canonical_index<<",\"coefficient\":";num(o,x.coefficient);o<<",\"weight\":";num(o,x.weight);o<<",\"source_spin\":"<<quoted(nbo_spin_name(x.source_spin))<<",\"source_side_indices\":";json_array(o,x.source_side_indices,[&](auto i){o<<i;});o<<",\"source_coefficients\":";json_array(o,x.source_coefficients,[&](double v){num(o,v);});o<<",\"basis_mapping\":";json_array(o,x.basis_mapping,[&](double v){num(o,v);});o<<'}';});
    o<<",\"coverage\":";json_array(o,m.coverage,[&](const auto& x){o<<"{\"canonical_index\":"<<x.canonical_index<<",\"available\":"<<(x.available?"true":"false")<<",\"weight_sum\":";num(o,x.weight_sum);o<<",\"residual_norm\":";num(o,x.residual_norm);o<<'}';});
    o<<",\"energies\":";json_array(o,m.energies,[&](const auto& x){o<<'{';energy_provenance(o,x);o<<"\"spin\":"<<quoted(nbo_spin_name(x.spin))<<",\"available\":"<<(x.available?"true":"false")<<",\"status\":"<<quoted(x.status)<<",\"detail\":"<<quoted(x.detail)<<",\"input_units\":"<<quoted(x.input_units)<<",\"canonical_columns_checked\":"<<x.canonical_columns_checked<<",\"effective_rank\":"<<x.effective_rank<<",\"effective_rank_semantics\":\"validated canonical column count; not overlap eigenvalue-floor rank\",\"overlap_numerical_rank\":"<<x.overlap_numerical_rank<<",\"canonical_effective_rank\":"<<x.canonical_effective_rank<<",\"canonical_null_directions\":"<<x.canonical_null_directions<<",\"overlap_rank_threshold\":";num(o,x.overlap_rank_threshold);o<<",\"outside_canonical_nao_norm\":";num(o,x.outside_canonical_nao_norm);o<<",\"outside_canonical_fock_coupling\":";num(o,x.outside_canonical_fock_coupling);o<<",\"nullspace_residual_semantics\":\"residual outside canonical dual metric projector\",\"hermiticity_error\":";num(o,x.hermiticity_error);o<<",\"canonical_residual\":";num(o,x.canonical_residual);o<<",\"projected_residual\":";num(o,x.projected_residual);o<<",\"nullspace_residual\":";num(o,x.nullspace_residual);o<<",\"eigenvalue_error_hartree\":";num(o,x.eigenvalue_error_hartree);o<<",\"fock_symmetry_error\":";num(o,x.fock_symmetry_error);o<<",\"density_symmetry_error\":";num(o,x.density_symmetry_error);o<<",\"density_symmetry_checked\":"<<(x.density_symmetry_checked?"true":"false");o<<",\"electronic_symmetry_verified\":"<<(x.electronic_symmetry_verified?"true":"false")<<",\"source_path\":"<<quoted(x.source.path)<<",\"source_line\":"<<x.source.line_begin<<'}';});
    o<<",\"diagnostics\":";json_array(o,m.diagnostics,[&](const auto& x){o<<quoted(x);});o<<'}';return o.str();}
} // namespace cov
