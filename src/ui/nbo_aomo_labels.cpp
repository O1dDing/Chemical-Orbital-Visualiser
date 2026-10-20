#include "cov/nbo_aomo_labels.hpp"
#include "cov/d2h_orbital_characters.hpp"
#include "cov/open_profile.hpp"
#include "cov/orbital_symmetry_scope.hpp"
#include "cov/orbital_symmetry.hpp"
#include "cov/molecular_point_group_frame.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/atomic_orbital_symmetry.hpp"
#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <iomanip>

namespace cov::ui {
namespace {
constexpr double character_tolerance=2e-4;
constexpr double metric_tolerance=2e-5;
constexpr double energy_order_tolerance=2e-5;
thread_local std::size_t canonical_name_revision=0;
using Vec=std::array<double,3>;
using Mat=std::array<double,9>;
using Dense=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
struct NamedCharacters {std::string label;std::size_t dimension;double norm=1;std::vector<double> values;};
double dot(const Vec& a,const Vec& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double det(const Mat& a){return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);}
double trace(const Mat& a){return a[0]+a[4]+a[8];}
double matrix_error(const Mat& a,const Mat& b){double e=0;for(std::size_t i=0;i<9;++i)e=std::max(e,std::abs(a[i]-b[i]));return e;}
Mat multiply(const Mat& a,const Mat& b){Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
const Mat identity{1,0,0,0,1,0,0,0,1};
const Mat inversion{-1,0,0,0,-1,0,0,0,-1};
std::string normalized(std::string s){std::string o;for(unsigned char c:s)if(!std::isspace(c)&&c!='_')o+=char(std::tolower(c));return o;}
void retain_dominant_name(NboAomoName& name){
    if(name.verified||!name.decomposition_verified||name.point_group.empty()||name.components.empty())return;
    double first=0,second=0;std::string irrep;
    for(const auto& c:name.components){if(c.weight>first){second=first;first=c.weight;irrep=c.irrep;}else second=std::max(second,c.weight);}
    if(first>=.95&&first-second>=.90){name.approximate_dominant_label=true;name.dominant_irrep=irrep;name.dominant_weight=first;}
}
struct CanonicalActionCache {
    const Wavefunction* wavefunction=nullptr;
    std::vector<double> packed,metric_packed;
    std::vector<std::vector<double>> transformed,metric_transformed;
};
struct Frame {
    std::string group, detail;
    std::vector<SymmetryOperation> ops;
    bool valid=false, linear=false;
    Vec axis{};
    unsigned lmax=0;
    PointGroupIrrepTable table;
    mutable std::shared_ptr<CanonicalActionCache> canonical_action;
};
Frame frame_for(const Wavefunction& w,const NboSalcModel* model){
    Frame f;
    if(!model||!model->group_verified)return f;
    f.group=model->symmetry_scope.reduced?model->used_group:model->point_group;f.linear=f.group=="Dinfh"||f.group=="Cinfv";
    if(model->used_group!=f.group&&(!f.linear||model->used_group!="finite sampling subgroup of "+f.group))return f;
    f.ops=model->operations;const auto n=f.ops.size();
    const auto expected_order=finite_point_group_order(f.group);
    if(!n||(!f.linear&&(!expected_order||expected_order!=n)))return f;
    if(f.linear){
        const auto geometry=analyse_molecular_symmetry(w);
        if(!geometry.linear||geometry.point_group!=f.group)return f;
        double best=0;for(const auto& a:w.atoms){Vec v{a.x-geometry.centre_bohr[0],a.y-geometry.centre_bohr[1],a.z-geometry.centre_bohr[2]};if(dot(v,v)>best){best=dot(v,v);f.axis=v;}}
        if(best<1e-20)return f;for(auto& v:f.axis)v/=std::sqrt(best);
        for(const auto& shell:w.shells)f.lmax=std::max(f.lmax,unsigned(shell.angular_momentum));
        // Sampling must resolve every supported angular frequency, without aliasing.
        if(f.lmax>4||n<(f.group=="Dinfh"?4:2)*std::max(3u,2*f.lmax+1))return f;
    }
    // SALC closure products intentionally carry only matrix/permutation data;
    // kind/order/axis metadata must not be used here.
    bool have_e=false;
    for(std::size_t i=0;i<n;++i){const auto& op=f.ops[i];
        if(op.atom_permutation.size()!=w.atoms.size())return f;
        std::set<std::size_t> seen;for(std::size_t a=0;a<w.atoms.size();++a){const auto b=op.atom_permutation[a];if(b>=w.atoms.size()||w.atoms[a].atomic_number!=w.atoms[b].atomic_number||!seen.insert(b).second)return f;}
        Mat trans{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)trans[3*a+b]=op.matrix[3*b+a];
        if(matrix_error(multiply(trans,op.matrix),identity)>metric_tolerance)return f;
        if(matrix_error(op.matrix,identity)<metric_tolerance)have_e=true;
    }
    if(!have_e)return f;
    // Recheck matrix closure, not just the advertised point-group string.
    for(const auto& a:f.ops)for(const auto& b:f.ops){const auto ab=multiply(a.matrix,b.matrix);if(std::none_of(f.ops.begin(),f.ops.end(),[&](const auto& c){return matrix_error(ab,c.matrix)<metric_tolerance;}))return f;}
    f.valid=true;
    if(!f.linear&&f.valid){
        auto group=analyse_molecular_symmetry(w);group.point_group=f.group;group.operations=f.ops;
        f.table=molecular_point_group_irreps(w,group);
        f.valid=f.table.valid;
        f.detail=f.table.axis_detail+(f.table.valid?"":"; "+f.table.reason);
    }
    if(f.detail.empty())f.detail=f.linear?"Validated angular-momentum-resolving sampling of "+f.group:"Validated full finite "+f.group+" operation matrices";
    if(model->symmetry_scope.reduced)f.detail+="; approximate total/spin-density symmetry scope; nuclear geometry="+model->point_group+"; physical Fock invariance is independently tested";
    return f;
}
Frame canonical_frame(const Wavefunction& w){
    // Use the same matrix, closure, character and metric gates as an attached
    // SALC frame. A geometry label alone does not certify orbital symmetry.
    const auto geometry=analyse_molecular_symmetry(w);
    NboSalcModel frame;
    frame.point_group=frame.used_group=geometry.point_group;
    frame.group_verified=!geometry.operations.empty();
    frame.operations=geometry.operations;
    const auto electronic=analyse_nbo_electronic_symmetry_scope(w);
    if(electronic.reduced){frame.symmetry_scope=electronic;frame.used_group=electronic.naming_group;frame.operations=electronic.operations;return frame_for(w,&frame);}
    if(!geometry.linear&&frame.group_verified){
        const auto complete=complete_molecular_point_group(w,geometry);
        frame.operations=complete.operations;frame.group_verified=!complete.operations.empty();
    }
    if(geometry.linear){
        Vec axis{};double length2=0;
        for(const auto& atom:w.atoms){
            Vec v{atom.x-geometry.centre_bohr[0],atom.y-geometry.centre_bohr[1],atom.z-geometry.centre_bohr[2]};
            if(dot(v,v)>length2){length2=dot(v,v);axis=v;}
        }
        if(length2<=1e-20)return frame_for(w,&frame);
        for(auto& x:axis)x/=std::sqrt(length2);
        const Vec seed=std::abs(axis[0])<.8?Vec{1,0,0}:Vec{0,1,0};
        Vec normal{axis[1]*seed[2]-axis[2]*seed[1],axis[2]*seed[0]-axis[0]*seed[2],axis[0]*seed[1]-axis[1]*seed[0]};
        const double normal_length=std::sqrt(dot(normal,normal));
        for(auto& x:normal)x/=normal_length;
        Mat reflection=identity;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)reflection[3*a+b]-=2*normal[a]*normal[b];
        unsigned lmax=0;for(const auto& shell:w.shells)lmax=std::max(lmax,unsigned(shell.angular_momentum));
        // Match the validated SALC sampling resolution; no angular order is
        // inferred from molecule identity or an orbital number.
        const unsigned order=std::max(3u,2*lmax+1);
        const auto inversion_op=std::find_if(geometry.operations.begin(),geometry.operations.end(),
            [](const auto& op){return op.kind==SymmetryOperationKind::Inversion;});
        const bool centrosymmetric=geometry.point_group=="Dinfh";
        if(centrosymmetric&&inversion_op==geometry.operations.end())return frame_for(w,&frame);
        frame.operations.clear();
        const Mat cross{0,-axis[2],axis[1],axis[2],0,-axis[0],-axis[1],axis[0],0};
        for(unsigned k=0;k<order;++k){
            const double angle=2*3.14159265358979323846*k/order,c=std::cos(angle),s=std::sin(angle);
            Mat rotation{};
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                rotation[3*a+b]=c*identity[3*a+b]+(1-c)*axis[a]*axis[b]+s*cross[3*a+b];
            for(int mirror=0;mirror<2;++mirror)for(int invert=0;invert<(centrosymmetric?2:1);++invert){
                SymmetryOperation op;op.matrix=mirror?multiply(rotation,reflection):rotation;
                if(invert){for(auto& x:op.matrix)x=-x;op.atom_permutation=inversion_op->atom_permutation;}
                else{op.atom_permutation.resize(w.atoms.size());std::iota(op.atom_permutation.begin(),op.atom_permutation.end(),0);}
                frame.operations.push_back(std::move(op));
            }
        }
        frame.used_group="finite sampling subgroup of "+geometry.point_group;
    }
    return frame_for(w,&frame);
}
// For a linear molecule O(2) characters are 2 cos(m theta), with an
// inversion parity in Dinfh and a reflection sign for Sigma. Compare EVERY
// sampled operation, not a visually inferred sigma/pi or one selected axis.
double metric(const Wavefunction& w,const std::vector<double>& a,const std::vector<double>& b){
    const auto n=std::size_t(w.basis_count);double v=0;for(std::size_t i=0;i<n;++i){double row=0;for(std::size_t j=0;j<n;++j)row+=w.ao_overlap[i*n+j]*b[j];v+=a[i]*row;}return v;
}
const CanonicalActionCache& canonical_action(const Wavefunction& w,const Frame& f,std::size_t operation){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size();
    const auto metric_columns=[&](const std::vector<double>& packed){
        std::vector<double> result(n*m,0);
        Eigen::Map<Dense>(result.data(),n,m).noalias()=
            Eigen::Map<const Dense>(w.ao_overlap.data(),n,n)*Eigen::Map<const Dense>(packed.data(),n,m);
        return result;
    };
    if(!f.canonical_action||f.canonical_action->wavefunction!=&w){
        f.canonical_action=std::make_shared<CanonicalActionCache>();auto& cache=*f.canonical_action;cache.wavefunction=&w;cache.packed.assign(n*m,0);
        for(std::size_t col=0;col<m;++col)if(w.orbitals[col].coefficients.size()==n)for(std::size_t row=0;row<n;++row)cache.packed[row*m+col]=w.orbitals[col].coefficients[row];
        cache.metric_packed=metric_columns(cache.packed);cache.transformed.resize(f.ops.size());cache.metric_transformed.resize(f.ops.size());
    }
    auto& cache=*f.canonical_action;
    if(cache.transformed[operation].empty()){
        cache.transformed[operation]=apply_orbital_symmetry_operation(w,f.ops[operation],cache.packed,m);
        if(cache.transformed[operation].size()==n*m)cache.metric_transformed[operation]=metric_columns(cache.transformed[operation]);
    }
    return cache;
}
// Complete S-orthonormal subspace character and leakage test. This examines
// the actual canonical coefficients, including core and all virtual rows.
struct CharacterFailure {
    std::string status,detail,quantity;
    std::optional<double> value,limit;
};
std::vector<double> characters(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& indices,CharacterFailure* failure=nullptr){
    const auto fail=[&](const char* status,const char* quantity,double value,double limit){
        if(failure){failure->status=status;failure->quantity=quantity;failure->value=value;failure->limit=limit;std::ostringstream text;text<<std::setprecision(12)<<quantity<<'='<<value<<"; limit="<<limit;failure->detail=text.str();}
        return std::vector<double>{};
    };
    const auto n=std::size_t(w.basis_count),k=indices.size(),m=w.orbitals.size();
    if(!f.valid||!n||!k||w.ao_overlap.size()!=n*n)return fail("symmetry_evidence_unavailable","AO metric dimension",double(w.ao_overlap.size()),double(n*n));
    const auto& action=canonical_action(w,f,0);Dense q(n,k),sq(n,k);
    for(std::size_t col=0;col<k;++col){const auto i=indices[col];if(i>=m||w.orbitals[i].coefficients.size()!=n)return {};
        for(std::size_t r=0;r<n;++r){q(r,col)=w.orbitals[i].coefficients[r];sq(r,col)=action.metric_packed[r*m+i];}}
    const Dense gram=q.transpose()*sq;
    if(!gram.allFinite()||(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff()>metric_tolerance)
        return fail("source_metric_not_orthonormal","source Gram error",(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff(),metric_tolerance);
    std::vector<double> result;
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);
        if(cache.transformed[g].size()!=n*m||cache.metric_transformed[g].size()!=n*m)return {};
        Dense tq(n,k),stq(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t r=0;r<n;++r){tq(r,c)=cache.transformed[g][r*m+indices[c]];stq(r,c)=cache.metric_transformed[g][r*m+indices[c]];}
        const Dense d=sq.transpose()*tq;
        const Dense residual=tq-q*d,residual_metric=stq-sq*d;
        for(std::size_t c=0;c<k;++c){const double norm=tq.col(c).dot(stq.col(c)),loss=residual.col(c).dot(residual_metric.col(c));
            if(!std::isfinite(norm)||std::abs(norm-1)>metric_tolerance)
                return fail("symmetry_action_not_isometric","transformed norm error",std::abs(norm-1),metric_tolerance);
            if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance)
                return fail("subspace_not_closed","S-metric leakage squared",loss,character_tolerance*character_tolerance);}
        result.push_back(d.trace());
    }return result;
}
// Residual of the real central idempotent on each unchanged source column.
// Both canonical MOs and reconstructed source SALCs use this same S-metric
// test after their complete containing span passes the closure/Gram gates.
std::vector<double> projection_residuals(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& members,std::size_t dimension,double norm,
    const std::vector<double>& row_characters){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size(),k=members.size();
    if(!f.valid||!n||!k||w.ao_overlap.size()!=n*n||row_characters.size()!=f.ops.size()||norm<=0)return {};
    Dense projected=Dense::Zero(n,k);
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);
        if(cache.transformed[g].size()!=n*m)return {};
        const double factor=double(dimension)*row_characters[g]/(double(f.ops.size())*norm);
        for(std::size_t c=0;c<k;++c){if(members[c]>=m||w.orbitals[members[c]].coefficients.size()!=n)return {};
            for(std::size_t a=0;a<n;++a)projected(a,c)+=factor*cache.transformed[g][a*m+members[c]];}}
    Dense residual(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a)residual(a,c)=w.orbitals[members[c]].coefficients[a]-projected(a,c);
    const Dense sr=Eigen::Map<const Dense>(w.ao_overlap.data(),n,n)*residual;
    std::vector<double> losses;for(std::size_t c=0;c<k;++c)losses.push_back(residual.col(c).dot(sr.col(c)));
    return losses;
}
bool pure_projection(double loss){return std::isfinite(loss)&&loss>=-1e-10&&loss<=character_tolerance*character_tolerance;}
void retain_projection_residual(NboAomoName& name,double loss){
    if(std::isfinite(loss)&&loss>=-1e-10&&(!name.projection_residual||loss<*name.projection_residual))
        name.projection_residual=std::max(0.,loss);
}
// Preserve each producer column. The derived central-projector components are
// expressed in its original containing basis, using the actual S Gram matrix.
// No orbital is rotated, replaced, renormalised or assigned a dominant irrep.
void retain_decomposition(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& members,const std::vector<NamedCharacters>& rows,
    std::vector<NboAomoName>& names,const std::string& source_kind="canonical",
    const Wavefunction* canonical_basis=nullptr,const std::vector<std::size_t>& basis_members={}){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size(),k=members.size();
    const auto fail=[&](const char* status){for(auto index:members){names[index].decomposition_status=status;names[index].components.clear();names[index].decomposition_verified=false;}};
    const auto& source=canonical_basis?*canonical_basis:w;
    const auto& source_members=canonical_basis?basis_members:members;
    if(!k||rows.empty()||source_members.empty()||characters(source,f,source_members).empty()||
       (canonical_basis&&source_members.size()!=n)){fail("containing_span_not_verified");return;}
    const auto& cache=canonical_action(w,f,0);Dense q(n,k),sq(n,k);
    for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a){q(a,c)=cache.packed[a*m+members[c]];sq(a,c)=cache.metric_packed[a*m+members[c]];}
    const Dense s=Eigen::Map<const Dense>(w.ao_overlap.data(),n,n);
    Dense basis(n,source_members.size());
    for(std::size_t c=0;c<source_members.size();++c)for(std::size_t a=0;a<n;++a)basis(a,c)=source.orbitals[source_members[c]].coefficients[a];
    const Dense sbasis=s*basis,gram=basis.transpose()*sbasis;
    const Dense target_gram=q.transpose()*sq;
    for(std::size_t c=0;c<k;++c)if(!std::isfinite(target_gram(c,c))||target_gram(c,c)<=0){fail("source_metric_invalid");return;}
    Eigen::LDLT<Dense> solve(gram);
    if(solve.info()!=Eigen::Success||!solve.isPositive()){fail("source_gram_solve_failed");return;}
    std::vector<Dense> projected,coordinates,metric_projected;
    Dense sum=Dense::Zero(n,k);double coordinate_loss=0;
    for(const auto& row:rows){
        if(row.values.size()!=f.ops.size()||row.norm<=0){fail("projector_table_unavailable");return;}
        Dense p=Dense::Zero(n,k);
        for(std::size_t g=0;g<f.ops.size();++g){const auto& action=canonical_action(w,f,g);
            if(action.transformed[g].size()!=n*m){fail("projector_action_unavailable");return;}
            const double factor=double(row.dimension)*row.values[g]/(double(f.ops.size())*row.norm);
            for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a)p(a,c)+=factor*action.transformed[g][a*m+members[c]];
        }
        Dense a=solve.solve(sbasis.transpose()*p),residual=p-basis*a,sres=s*residual;
        if(solve.info()!=Eigen::Success||!a.allFinite()||!p.allFinite()){fail("source_gram_solve_failed");return;}
        for(std::size_t c=0;c<k;++c){const double loss=residual.col(c).dot(sres.col(c))/target_gram(c,c);
            if(!std::isfinite(loss)||loss< -1e-10){fail("component_metric_invalid");return;}
            coordinate_loss=std::max(coordinate_loss,loss);}
        sum+=basis*a;metric_projected.push_back(s*p);projected.push_back(std::move(p));coordinates.push_back(std::move(a));
    }
    if(coordinate_loss>character_tolerance*character_tolerance){fail("component_outside_source_span");return;}
    const Dense residual=q-sum,sres=s*residual;
    for(std::size_t c=0;c<k;++c){
        auto& name=names[members[c]];const double norm=target_gram(c,c),loss=residual.col(c).dot(sres.col(c))/norm;
        double weights=0,cross=0;std::vector<NboIrrepComponent> components;
        for(std::size_t r=0;r<rows.size();++r){
            const double weight=projected[r].col(c).dot(metric_projected[r].col(c))/norm;
            if(!std::isfinite(weight)||weight< -1e-10){fail("component_metric_invalid");return;}
            NboIrrepComponent component;component.irrep=rows[r].label;component.dimension=rows[r].dimension;component.weight=std::max(0.,weight);
            for(std::size_t j=0;j<source_members.size();++j)component.source_coefficients.push_back(coordinates[r](j,c));
            weights+=component.weight;components.push_back(std::move(component));
            for(std::size_t t=0;t<r;++t)cross=std::max(cross,std::abs(projected[r].col(c).dot(metric_projected[t].col(c))/norm));
        }
        name.decomposition_reconstruction_residual=std::sqrt(std::max(0.,loss));
        name.decomposition_orthogonality_error=cross;name.decomposition_weight_sum_error=std::abs(weights-1);
        if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance||
           cross>metric_tolerance||std::abs(weights-1)>metric_tolerance){
            name.decomposition_status="projector_reconstruction_or_orthogonality_failed";continue;
        }
        name.components=std::move(components);name.decomposition_verified=true;
        name.component_source_kind=source_kind;name.component_source_members=source_members;
        name.point_group=f.group;
        name.decomposition_status="verified_source_symmetry_composition";
    }
}
// Every edge is measured from the actual group action in the AO metric.
// The scaled edge bound prevents many individually small couplings from
// evading the complete-column leakage test. Energy never excludes partners.
std::vector<std::vector<std::size_t>> connected_blocks(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& rows){
    const auto n=std::size_t(w.basis_count),k=rows.size(),m=w.orbitals.size();
    const auto singles=[&](){std::vector<std::vector<std::size_t>> v;for(auto i:rows)v.push_back({i});return v;};
    if(k<2||!f.valid||!n||w.ao_overlap.size()!=n*n)return singles();
    const auto& action=canonical_action(w,f,0);Dense q(n,k),sq(n,k);
    for(std::size_t c=0;c<k;++c){if(rows[c]>=m||w.orbitals[rows[c]].coefficients.size()!=n)return singles();
        for(std::size_t r=0;r<n;++r){q(r,c)=w.orbitals[rows[c]].coefficients[r];sq(r,c)=action.metric_packed[r*m+rows[c]];}}
    const Dense gram=q.transpose()*sq;
    if(!gram.allFinite()||(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff()>metric_tolerance)return singles();
    std::vector<std::size_t> parent(k);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](std::size_t a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;};
    const double edge_bound=character_tolerance/std::sqrt(double(k));
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);if(cache.transformed[g].size()!=n*m)return singles();
        Dense tq(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t r=0;r<n;++r)tq(r,c)=cache.transformed[g][r*m+rows[c]];
        const Dense d=sq.transpose()*tq;if(!d.allFinite())return singles();
        for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<a;++b)
            if(std::max(std::abs(d(a,b)),std::abs(d(b,a)))>edge_bound)parent[root(b)]=root(a);
    }
    std::map<std::size_t,std::vector<std::size_t>> groups;
    for(std::size_t a=0;a<k;++a)groups[root(a)].push_back(rows[a]);
    std::vector<std::vector<std::size_t>> result;for(auto& [key,members]:groups)result.push_back(std::move(members));return result;
}
std::vector<std::vector<std::size_t>> canonical_blocks(const Wavefunction& w,const Frame& f){
    std::vector<std::vector<std::size_t>> result;
    for(const auto spin:{Spin::Alpha,Spin::Beta}){std::vector<std::size_t> rows;
        for(std::size_t i=0;i<w.orbitals.size();++i)if(w.orbitals[i].spin==spin)rows.push_back(i);
        auto blocks=connected_blocks(w,f,rows);result.insert(result.end(),blocks.begin(),blocks.end());
    }return result;
}
// Restrict the measured action to unchanged, individually pure source columns.
// A graph component is only a proposal: its dimension, complete character row,
// S Gram, transformed norms and leakage must all certify one actual copy.
// Energy never removes a coupling or supplies missing partners.
std::vector<std::vector<std::size_t>> action_hierarchy_candidates(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& pure,std::size_t dimension){
    std::vector<std::vector<std::size_t>> candidates;
    const auto k=pure.size(),n=std::size_t(w.basis_count),m=w.orbitals.size();
    if(!dimension||k<dimension)return candidates;
    if(dimension==1){for(auto i:pure)candidates.push_back({i});return candidates;}
    const auto& initial=canonical_action(w,f,0);Dense sq(n,k),strength=Dense::Zero(k,k);
    for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a)sq(a,c)=initial.metric_packed[a*m+pure[c]];
    for(std::size_t g=0;g<f.ops.size();++g){const auto& action=canonical_action(w,f,g);if(action.transformed[g].size()!=n*m)return {};
        Dense tq(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a)tq(a,c)=action.transformed[g][a*m+pure[c]];
        const Dense d=sq.transpose()*tq;if(!d.allFinite())return {};
        for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<a;++b)strength(a,b)=std::max(strength(a,b),std::max(std::abs(d(a,b)),std::abs(d(b,a))));
    }
    struct Edge{double weight;std::size_t a,b;};std::vector<Edge> edges;
    for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<a;++b)if(strength(a,b)>0)edges.push_back({strength(a,b),a,b});
    std::sort(edges.begin(),edges.end(),[](const auto& a,const auto& b){if(a.weight!=b.weight)return a.weight>b.weight;return std::tie(a.a,a.b)<std::tie(b.a,b.b);});
    std::vector<std::size_t> parent(k);std::iota(parent.begin(),parent.end(),0);
    std::vector<std::vector<std::size_t>> clusters(k);for(std::size_t a=0;a<k;++a)clusters[a]={pure[a]};
    auto root=[&](std::size_t a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;};
    for(const auto& edge:edges){auto a=root(edge.a),b=root(edge.b);if(a==b)continue;
        parent[b]=a;clusters[a].insert(clusters[a].end(),clusters[b].begin(),clusters[b].end());clusters[b].clear();
        if(clusters[a].size()==dimension){auto candidate=clusters[a];std::sort(candidate.begin(),candidate.end());candidates.push_back(std::move(candidate));}}
    return candidates;
}
std::vector<std::vector<std::size_t>> verified_source_copies(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& pure,std::size_t dimension,double norm,
    const std::vector<double>& row_characters,std::size_t maximum_copies,
    std::set<std::size_t>* separable_members=nullptr,
    std::vector<NboAomoName>* names=nullptr){
    std::vector<std::vector<std::size_t>> result;
    if(names)for(auto i:pure)(*names)[i].partner_status="original_partner_subset_not_verified";
    auto candidates=connected_blocks(w,f,pure);
    auto separable=candidates;for(auto& candidate:separable)std::sort(candidate.begin(),candidate.end());
    // A conservative global graph can join nearly independent copies through
    // tiny entries even when each source copy passes the unchanged full-action
    // residual gate. Descending measured couplings propose additional subsets;
    // neither an energy window nor a relaxed acceptance threshold is involved.
    const auto hierarchy=action_hierarchy_candidates(w,f,pure,dimension);
    candidates.insert(candidates.end(),hierarchy.begin(),hierarchy.end());
    for(auto& candidate:candidates)std::sort(candidate.begin(),candidate.end());
    std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
    for(const auto& candidate:candidates){
        if(candidate.size()!=dimension)continue;
        CharacterFailure failure;
        const auto values=characters(w,f,candidate,&failure);
        if(values.size()!=row_characters.size()){
            if(names&&!failure.status.empty())for(auto i:candidate){auto& name=(*names)[i];
                name.partner_status=failure.status;name.partner_failure_quantity=failure.quantity;
                name.partner_failure_value=failure.value;name.partner_failure_limit=failure.limit;}
            continue;
        }
        bool matches=true;for(std::size_t g=0;g<values.size();++g)
            if(std::abs(values[g]-row_characters[g])>character_tolerance)matches=false;
        const auto losses=projection_residuals(w,f,candidate,dimension,norm,row_characters);
        if(losses.size()!=candidate.size()||!std::all_of(losses.begin(),losses.end(),pure_projection))matches=false;
        if(matches)result.push_back(candidate);
    }
    // Reject overlapping alternative source-copy assignments and inconsistent
    // counting evidence rather than letting a graph tie choose scientific IDs.
    std::set<std::size_t> used;
    for(const auto& candidate:result)for(auto i:candidate)if(!used.insert(i).second){
        if(names)for(auto j:pure)(*names)[j].partner_status="overlapping_source_copy_assignments";
        return {};
    }
    // Connected components are disjoint. Reject inconsistent counting evidence
    // instead of choosing a subset of equally plausible source copies.
    if(result.size()>maximum_copies){
        if(names)for(auto i:pure)(*names)[i].partner_status="source_copy_count_inconsistent";
        result.clear();
    }
    if(names)for(const auto& candidate:result)for(auto i:candidate){auto& name=(*names)[i];
        name.partner_status="verified_original_partner_block";name.partner_failure_quantity.clear();
        name.partner_failure_value.reset();name.partner_failure_limit.reset();}
    if(separable_members)for(const auto& candidate:result)if(std::find(separable.begin(),separable.end(),candidate)!=separable.end())
        separable_members->insert(candidate.begin(),candidate.end());
    return result;
}
struct Unit {
    std::vector<std::size_t> members;std::string irrep,scope,detail;
    double energy=0;bool energy_available=true;std::size_t copies=1;
    bool copies_resolved=true,counting_only=false;
    double lower=std::numeric_limits<double>::quiet_NaN(),upper=std::numeric_limits<double>::quiet_NaN();
};
void certify_salc_energy_bounds(Unit& unit,const Wavefunction& side,const Wavefunction& canonical,const NboSalcModel& model){
    if(unit.members.empty())return;const auto spin=model.orbitals[unit.members.front()].spin;
    if(spin==NboSpin::Total&&std::any_of(canonical.orbitals.begin(),canonical.orbitals.end(),[](const auto& mo){return mo.spin==Spin::Beta;}))return;
    const NboSalcEnergyEvidence* proof=nullptr;for(const auto& energy:model.energies)if(energy.spin==spin&&energy.available&&energy.electronic_symmetry_verified&&energy.status=="verified_same_operator")proof=&energy;
    if(!proof)return;
    const auto n=std::size_t(canonical.basis_count),k=unit.members.size();if(!n||canonical.ao_overlap.size()!=n*n)return;
    std::vector<std::size_t> columns;std::map<std::size_t,std::size_t> column_map;
    for(std::size_t i=0;i<canonical.orbitals.size();++i)if(canonical.orbitals[i].spin==(spin==NboSpin::Beta?Spin::Beta:Spin::Alpha)){if(canonical.orbitals[i].coefficients.size()!=n||!std::isfinite(canonical.orbitals[i].energy_hartree))return;column_map[i]=columns.size();columns.push_back(i);}
    if(columns.empty()||columns.size()!=proof->canonical_columns_checked)return;
    const auto m=columns.size();std::vector<double> projections(k*m,0);std::vector<bool> seen(k*m,false);
    for(const auto& link:model.links){auto member=std::find(unit.members.begin(),unit.members.end(),link.side_index);if(member==unit.members.end())continue;const auto col=column_map.find(link.canonical_index);if(col==column_map.end())continue;const auto at=std::size_t(member-unit.members.begin())*m+col->second;if(seen[at]||!std::isfinite(link.coefficient))return;seen[at]=true;projections[at]=link.coefficient;}
    if(std::find(seen.begin(),seen.end(),false)!=seen.end())return;
    // The producer-validated canonical eigenbasis must actually span each
    // retained SALC, in the same S metric and spin. Missing virtual columns
    // cannot silently turn a projected operator into a certified spectrum.
    for(std::size_t a=0;a<k;++a){const auto row=unit.members[a];if(row>=side.orbitals.size()||side.orbitals[row].coefficients.size()!=n)return;auto residual=side.orbitals[row].coefficients;for(std::size_t j=0;j<m;++j)for(std::size_t r=0;r<n;++r)residual[r]-=projections[a*m+j]*canonical.orbitals[columns[j]].coefficients[r];const double loss=metric(canonical,residual,residual);if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance)return;}
    std::vector<double> f(k*k,0);double gram_error=0;
    for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<k;++b){double gram=0;for(std::size_t j=0;j<m;++j){const double product=projections[a*m+j]*projections[b*m+j];gram+=product;f[a*k+b]+=product*canonical.orbitals[columns[j]].energy_hartree;}gram_error=std::max(gram_error,std::abs(gram-(a==b?1.:0.)));}
    if(!std::isfinite(gram_error)||gram_error>metric_tolerance||double(k)*gram_error>=.1)return;
    double lower=std::numeric_limits<double>::infinity(),upper=-lower;
    for(std::size_t a=0;a<k;++a){const auto energy=model.orbitals[unit.members[a]].energy_hartree;if(!energy||!std::isfinite(*energy)||std::abs(f[a*k+a]-*energy)>energy_order_tolerance)return;double radius=0;for(std::size_t b=0;b<k;++b)if(a!=b)radius+=std::abs(f[a*k+b]);lower=std::min(lower,f[a*k+a]-radius);upper=std::max(upper,f[a*k+a]+radius);}
    // Gershgorin bounds concern the complete projected Hermitian operator,
    // not the diagonal expectations or their mean. Enlarge for metric and
    // validated producer residuals; never rotate the displayed SALC basis.
    const double defect=double(k)*gram_error;const double padding=std::max(std::abs(lower),std::abs(upper))*defect/(1-defect)+energy_order_tolerance+std::abs(proof->eigenvalue_error_hartree);
    if(!std::isfinite(lower)||!std::isfinite(upper)||!std::isfinite(padding))return;unit.lower=lower-padding;unit.upper=upper+padding;
    std::ostringstream detail;detail<<"; same-spin complete-canonical projected Fock spectral bounds=["<<unit.lower<<','<<unit.upper<<"] Ha (Gershgorin with metric/error margin)";unit.detail+=detail.str();
}
std::string orbital_label(const std::string& irrep){
    auto label=format_symmetry_unicode(normalized(irrep));
    // These labels name one-electron orbitals. Preserve the machine irrep
    // identifier while using the spectroscopic lowercase orbital convention.
    const std::pair<const char*,const char*> symbols[]={{"Σ","σ"},{"Π","π"},{"Δ","δ"},{"Φ","φ"},{"Γ","γ"}};
    for(const auto& [upper,lower]:symbols)if(label.starts_with(upper)){label.replace(0,std::string(upper).size(),lower);break;}
    return label;
}
std::string numbered_orbital_label(const NboAomoName& name){
    const auto irrep=orbital_label(name.irrep);
    if(!name.ordinal)return irrep;
    // Atomic copy order is not the principal quantum number in n-l notation.
    if(name.point_group=="SO(3)")return irrep+"#"+std::to_string(name.ordinal);
    return std::to_string(name.ordinal)+irrep;
}
void assign_ordinals(std::vector<Unit> units,std::vector<NboAomoName>& names){
    for(auto& u:units)if(u.copies_resolved&&u.energy_available&&std::isfinite(u.energy)&&(!u.irrep.empty()||u.scope.starts_with("canonical ")))u.lower=u.upper=u.energy;
    for(std::size_t target=0;target<units.size();++target){const auto& u=units[target];if(u.irrep.empty()||u.counting_only)continue;
        std::size_t ordinal=1;bool resolved=u.copies_resolved&&u.energy_available&&std::isfinite(u.lower)&&std::isfinite(u.upper);
        std::string reason=resolved?"verified_complete_set_order":!u.copies_resolved?"source_copy_membership_unresolved":
            !u.energy_available?"source_energy_unavailable":"source_spectral_bounds_unavailable";
        std::vector<std::size_t> blockers;
        for(std::size_t other=0;resolved&&other<units.size();++other){if(other==target)continue;const auto& v=units[other];if(v.scope!=u.scope||(!v.irrep.empty()&&v.irrep!=u.irrep))continue;
            if(std::isfinite(v.lower)&&v.lower>u.upper+energy_order_tolerance)continue;
            if(!v.irrep.empty()&&std::isfinite(v.upper)&&v.upper<u.lower-energy_order_tolerance){ordinal+=v.copies;continue;}
            if(!v.irrep.empty()&&v.copies_resolved&&v.energy_available&&std::isfinite(v.energy)){
                // Two independent occurrences at unresolved equal energy
                // have no physical copy order. Source row is not evidence.
                if(std::abs(v.energy-u.energy)<=energy_order_tolerance){resolved=false;reason="same_irrep_energy_order_unresolved";blockers=v.members;}
                else if(v.energy<u.energy)ordinal+=v.copies;
            }else {
                resolved=false;blockers=v.members;
                reason=v.irrep.empty()?"earlier_irrep_count_unavailable":
                    !std::isfinite(v.lower)||!std::isfinite(v.upper)?"other_copy_spectral_bounds_unavailable":
                    "overlapping_irrep_spectral_bounds";
            }
        }
        // Ranking bounds use every relevant occurrence, not just the first
        // rejection. Close-energy relationships need not be transitive.
        bool bounded=u.copies_resolved&&u.energy_available&&std::isfinite(u.lower)&&std::isfinite(u.upper);
        std::size_t lower=1,upper=1;
        for(std::size_t other=0;bounded&&other<units.size();++other){if(other==target)continue;const auto& v=units[other];
            if(v.scope!=u.scope||(!v.irrep.empty()&&v.irrep!=u.irrep))continue;
            if(std::isfinite(v.lower)&&v.lower>u.upper+energy_order_tolerance)continue;
            if(v.irrep.empty()||!v.copies||!std::isfinite(v.lower)||!std::isfinite(v.upper)){bounded=false;break;}
            if(v.upper<u.lower-energy_order_tolerance){lower+=v.copies;upper+=v.copies;}
            else upper+=v.copies;
        }
        for(auto index:u.members){auto& name=names[index];name.irrep=u.irrep;name.verified=true;name.ordinal=resolved?ordinal:0;name.representation_multiplicity=u.copies;name.ordinal_status=reason;name.ordinal_blocking_members=blockers;name.label=(name.ordinal?std::to_string(name.ordinal):std::string{})+orbital_label(u.irrep);name.detail=u.detail+"; ordinal scope="+u.scope+"; complete-set counts with certified subspace energy ordering; ordinal status="+reason;
            if(bounded){name.complete_set_ordinal_lower=lower;name.complete_set_ordinal_upper=upper;}}
    }
}
bool retain_atomic_components(const Wavefunction& target,std::size_t index,
    const AtomicOrbitalAssignment& assignment,const Wavefunction& source,
    const std::vector<std::size_t>& source_members,const std::string& source_kind,NboAomoName& name){
    const auto n=std::size_t(target.basis_count),k=source_members.size();
    if(!assignment.decomposition_verified||!n||!k||target.orbitals[index].coefficients.size()!=n||
       source.basis_count!=n||source.ao_overlap.size()!=n*n)return false;
    const Dense s=Eigen::Map<const Dense>(source.ao_overlap.data(),n,n);
    Dense basis(n,k);for(std::size_t c=0;c<k;++c){if(source_members[c]>=source.orbitals.size()||source.orbitals[source_members[c]].coefficients.size()!=n)return false;
        for(std::size_t a=0;a<n;++a)basis(a,c)=source.orbitals[source_members[c]].coefficients[a];}
    const Dense sbasis=s*basis,gram=basis.transpose()*sbasis;
    if(!gram.allFinite()||(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff()>metric_tolerance)return false;
    Eigen::LDLT<Dense> solve(gram);if(solve.info()!=Eigen::Success||!solve.isPositive())return false;
    Eigen::VectorXd original=Eigen::Map<const Eigen::VectorXd>(target.orbitals[index].coefficients.data(),n),sum=Eigen::VectorXd::Zero(n);
    const double norm=original.dot(s*original);if(!std::isfinite(norm)||norm<=0)return false;
    std::vector<NboIrrepComponent> components;std::vector<Eigen::VectorXd> projected;double cross=0,weight_sum=0;
    for(const auto& component:assignment.components){if(component.ao_coefficients.size()!=n||component.angular_momentum<0)return false;
        const Eigen::VectorXd p=Eigen::Map<const Eigen::VectorXd>(component.ao_coefficients.data(),n),a=solve.solve(sbasis.transpose()*p),represented=basis*a,residual=p-represented;
        const double loss=residual.dot(s*residual)/norm;
        if(!a.allFinite()||!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance)return false;
        const double weight=p.dot(s*p)/norm;if(!std::isfinite(weight)||weight< -1e-10)return false;
        for(const auto& other:projected)cross=std::max(cross,std::abs(p.dot(s*other)/norm));
        NboIrrepComponent saved;saved.irrep=component.label;saved.dimension=2*component.angular_momentum+1;saved.weight=std::max(0.,weight);
        saved.source_coefficients.assign(a.data(),a.data()+a.size());components.push_back(std::move(saved));projected.push_back(p);sum+=represented;weight_sum+=std::max(0.,weight);
    }
    const Eigen::VectorXd residual=original-sum;const double loss=residual.dot(s*residual)/norm;
    name.decomposition_reconstruction_residual=std::sqrt(std::max(0.,loss));name.decomposition_orthogonality_error=cross;name.decomposition_weight_sum_error=std::abs(weight_sum-1);
    if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance||cross>metric_tolerance||std::abs(weight_sum-1)>metric_tolerance)return false;
    name.decomposition_verified=true;name.decomposition_status="verified_atomic_angular_composition";
    name.components=std::move(components);name.component_source_kind=source_kind;name.component_source_members=source_members;return true;
}
void atomic_name_rows(const Wavefunction& target,const AtomicOrbitalSymmetryResult& analysis,
    const Wavefunction& canonical,std::vector<NboAomoName>& names,const NboSalcModel* salc=nullptr){
    const bool explicit_spin=canonical.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin||std::any_of(canonical.orbitals.begin(),canonical.orbitals.end(),[](const auto& mo){return mo.spin==Spin::Beta;});
    const bool shared_spatial=canonical.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared&&!explicit_spin;
    for(std::size_t i=0;i<names.size();++i){auto& name=names[i];name.point_group="SO(3)";
        name.ordinal_scope="atomic radial-copy order; not principal n";name.status=analysis.status;name.detail=analysis.detail;
        if(!analysis.available||i>=analysis.orbitals.size()){name.label=salc?"SALC "+std::to_string(i+1):canonical_mo_source_label(target,i);continue;}
        const auto& assignment=analysis.orbitals[i];name.status=assignment.status;name.detail+="; "+assignment.detail+"; "+assignment.partner_status;
        name.verified=assignment.angular_momentum_verified;name.irrep=assignment.label;name.containing_members=assignment.containing_members;
        if(name.containing_members.empty())name.containing_members={i};
        name.representation_multiplicity=assignment.representation_multiplicity;
        if(name.verified&&std::isfinite(assignment.projection_residual_squared))name.projection_residual=assignment.projection_residual_squared;
        if(assignment.containing_span_closed&&name.verified)name.containing_irreps.push_back({assignment.label,std::size_t(2*assignment.angular_momentum+1),assignment.representation_multiplicity});
        if(assignment.partner_block_verified){name.partner_block_id=(salc?"salc:":"canonical:")+assignment.partner_block_id;name.partner_block_size=assignment.containing_members.size();}
        name.ordinal=assignment.radial_copy_ordinal;
        // Mixed columns may carry another unresolved radial occurrence of this
        // l. A sequence over only the pure columns would silently omit it.
        for(std::size_t j=0;name.ordinal&&j<analysis.orbitals.size();++j)if(target.orbitals[j].spin==target.orbitals[i].spin&&!analysis.orbitals[j].angular_momentum_verified){
            if(!analysis.orbitals[j].decomposition_verified)name.ordinal=0;
            for(const auto& component:analysis.orbitals[j].components)if(component.angular_momentum==assignment.angular_momentum&&component.weight>character_tolerance*character_tolerance)name.ordinal=0;
        }
        name.complete_set_ordinal=name.ordinal;
        name.partner_status=assignment.partner_status;
        name.ordinal_status=name.ordinal?"verified_atomic_radial_copy_order":"atomic_radial_copy_order_unresolved";
        name.complete_set_ordinal_status=name.ordinal_status;
        if(name.ordinal){name.complete_set_ordinal_lower=name.ordinal;name.complete_set_ordinal_upper=name.ordinal;}
        std::vector<std::size_t> source_members;
        for(std::size_t j=0;j<canonical.orbitals.size();++j)if(canonical.orbitals[j].spin==(shared_spatial?Spin::Alpha:target.orbitals[i].spin)&&canonical.orbitals[j].coefficients.size()==target.basis_count)source_members.push_back(j);
        const bool allow_canonical=!salc||salc->orbitals[i].spin!=NboSpin::Total||!explicit_spin;
        bool retained=allow_canonical&&retain_atomic_components(target,i,assignment,canonical,source_members,"canonical",name);
        if(!retained&&salc){const auto local=assignment.partner_block_verified?assignment.containing_members:std::vector<std::size_t>{i};retained=retain_atomic_components(target,i,assignment,target,local,"salc",name);}
        if(!retained)name.decomposition_status=assignment.decomposition_verified?"atomic_component_source_reconstruction_unavailable":assignment.status;
        if(!salc)name.label=canonical_mo_display_label(target,i,&name);
        else name.label=(name.verified?numbered_orbital_label(name)+" ":std::string{})+"[SALC "+std::to_string(i+1)+"]";
    }
}
NboAomoNames atomic_names(const Wavefunction& w,const NboIntegration& data,const NboSalcModel* salc){
    NboAomoNames out;out.canonical.resize(w.orbitals.size());atomic_name_rows(w,analyse_atomic_orbital_symmetry(w),w,out.canonical);
    if(!salc)return out;out.salc.resize(salc->orbitals.size());
    const bool associated=(salc->dataset_id.empty()||data.id.empty()||salc->dataset_id==data.id)&&
        (salc->canonical_fingerprint.empty()||data.canonical_fingerprint.empty()||salc->canonical_fingerprint==data.canonical_fingerprint);
    Wavefunction side;side.atoms=w.atoms;side.shells=w.shells;side.primitives=w.primitives;side.basis_count=w.basis_count;side.ao_overlap=w.ao_overlap;side.orbitals.resize(salc->orbitals.size());
    for(std::size_t i=0;i<salc->orbitals.size();++i){const auto& orbital=salc->orbitals[i];auto& mo=side.orbitals[i];mo.energy_hartree=orbital.energy_hartree.value_or(std::numeric_limits<double>::quiet_NaN());mo.spin=orbital.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        bool valid=associated&&!orbital.terms.empty();mo.coefficients.assign(w.basis_count,0);
        for(const auto& term:orbital.terms){const auto* descriptor=nbo_orbital(data,term.orbital);if(!descriptor||descriptor->coefficients.size()!=w.basis_count||!std::isfinite(term.coefficient)){valid=false;break;}for(std::size_t a=0;a<w.basis_count;++a)mo.coefficients[a]+=term.coefficient*descriptor->coefficients[a];}
        if(!valid)mo.coefficients.clear();}
    atomic_name_rows(side,analyse_atomic_orbital_symmetry(side),w,out.salc,salc);return out;
}
} // namespace

NboAomoNames build_nbo_aomo_names(const Wavefunction& w,const NboIntegration& data,const NboSalcModel* salc){
    if(w.atoms.size()==1)return atomic_names(w,data,salc);
    NboAomoNames out;out.canonical.resize(w.orbitals.size());if(salc)out.salc.resize(salc->orbitals.size());
    const bool open=std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& mo){return mo.spin==Spin::Beta;});
    const bool associated=!salc || ((salc->dataset_id.empty()||data.id.empty()||salc->dataset_id==data.id)&&
        (salc->canonical_fingerprint.empty()||data.canonical_fingerprint.empty()||salc->canonical_fingerprint==data.canonical_fingerprint));
    const auto frame=salc?frame_for(w,associated?salc:nullptr):canonical_frame(w);std::vector<Unit> units;std::vector<bool> taken(w.orbitals.size(),false);
    std::vector<NamedCharacters> known;
    for(const auto& row:frame.table.rows)known.push_back({row.label,row.dimension,row.character_norm,row.characters});
    if(frame.valid&&frame.linear){
        const char* symbols[]={"Sigma","Pi","Delta","Phi","Gamma"};
        for(unsigned l=0;l<=frame.lmax;++l)for(int parity:{1,-1})for(int sign:{1,-1}){
            if(frame.group=="Cinfv"&&parity<0)continue;if(l&&sign<0)continue;
            NamedCharacters row;row.label=symbols[l];row.dimension=l?2:1;
            if(frame.group=="Dinfh")row.label+=parity>0?"_g":"_u";
            if(!l)row.label+=sign>0?"+":"-";
            for(const auto& op:frame.ops){Vec z{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)z[a]+=op.matrix[3*a+b]*frame.axis[b];
                const bool reverse=dot(z,frame.axis)<0;Mat r=op.matrix;if(reverse)for(auto& x:r)x=-x;
                const bool reflect=det(r)<0;
                double chi=l?(reflect?0:2*std::cos(l*std::acos(std::clamp((trace(r)-1)/2,-1.,1.)))):(reflect?sign:1);
                if(reverse)chi*=parity;row.values.push_back(chi);
            }known.push_back(std::move(row));
        }
    }
    const auto decompose=[&](const std::vector<double>& values,std::size_t dimension){
        std::vector<std::pair<std::size_t,std::size_t>> content;std::vector<double> reconstructed(values.size(),0);std::size_t dimensions=0;
        if(values.empty()||known.empty())return content;
        for(std::size_t r=0;r<known.size();++r){const auto& row=known[r];if(row.values.size()!=values.size()||row.norm<=0)return decltype(content){};
            double inner=std::inner_product(values.begin(),values.end(),row.values.begin(),0.)/(double(values.size())*row.norm);
            if(!std::isfinite(inner))return decltype(content){};const auto copies=std::llround(inner);
            if(copies<0||std::abs(inner-double(copies))>character_tolerance)return decltype(content){};
            if(!copies)continue;content.emplace_back(r,std::size_t(copies));dimensions+=row.dimension*std::size_t(copies);
            for(std::size_t g=0;g<values.size();++g)reconstructed[g]+=double(copies)*row.values[g];
        }
        if(dimensions!=dimension)return decltype(content){};
        for(std::size_t g=0;g<values.size();++g)if(std::abs(values[g]-reconstructed[g])>character_tolerance)return decltype(content){};
        return content;
    };
    const auto blocks=canonical_blocks(w,frame);
    for(std::size_t b=0;b<blocks.size();++b){const auto& members=blocks[b];if(members.empty())continue;
        CharacterFailure failure;const auto values=characters(w,frame,members,&failure);const auto content=decompose(values,members.size());
        const auto scope=w.orbitals[members.front()].spin==Spin::Beta?"canonical beta":"canonical alpha";
        double low=std::numeric_limits<double>::infinity(),high=-low;bool have_energy=true;
        for(auto i:members){const double e=w.orbitals[i].energy_hartree;have_energy=have_energy&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);
            auto& name=out.canonical[i];name.label="?";name.containing_members=members;
            name.status=!frame.valid?"operation_frame_unavailable":values.empty()?(failure.status.empty()?"symmetry_evidence_unavailable":failure.status):content.empty()?"character_decomposition_unavailable":"mixed_source_member";
            name.detail="Actual full AO coefficients; "+frame.detail+"; "+name.status+(failure.detail.empty()?"":"; "+failure.detail);
            for(const auto& [r,copies]:content)name.containing_irreps.push_back({known[r].label,known[r].dimension,copies});
        }
        if(content.empty()){
            for(auto i:members){Unit u;u.members={i};u.scope=scope;u.energy=w.orbitals[i].energy_hartree;u.energy_available=std::isfinite(u.energy);units.push_back(std::move(u));}
            continue;
        }
        retain_decomposition(w,frame,members,known,out.canonical);
        std::vector<std::size_t> proved_counts(known.size(),0),separable_counts(known.size(),0);std::set<std::size_t> proved_members,separable_members;
        for(const auto& [r,copies]:content){const auto& row=known[r];std::vector<std::size_t> pure;
            // Central idempotent projector, with real conjugate-pair norm.
            // It tests each *unchanged source column*. No rotated orbital is
            // substituted for a mixed source MO merely to obtain a name.
            const auto losses=projection_residuals(w,frame,members,row.dimension,row.norm,row.values);
            for(std::size_t c=0;c<losses.size();++c){retain_projection_residual(out.canonical[members[c]],losses[c]);if(pure_projection(losses[c]))pure.push_back(members[c]);}
            for(auto i:pure){auto& name=out.canonical[i];name.verified=true;name.irrep=row.label;name.point_group=frame.group;name.representation_multiplicity=copies;
                name.status="verified_isotypic_member";name.detail="Source column lies in verified isotypic projection; independent occurrence membership unresolved; "+frame.detail;}
            const auto source_copies=verified_source_copies(w,frame,pure,row.dimension,row.norm,row.values,copies,&separable_members,&out.canonical);
            for(std::size_t copy=0;copy<source_copies.size();++copy){const auto& partners=source_copies[copy];
                if(std::any_of(partners.begin(),partners.end(),[&](auto i){return proved_members.contains(i);}))continue;
                Unit u;u.members=partners;u.irrep=row.label;u.scope=scope;u.energy=0;u.energy_available=true;
                u.detail="Verified complete irreducible source span under every operation; "+frame.detail;
                for(auto i:partners){u.energy+=w.orbitals[i].energy_hartree/double(partners.size());u.energy_available=u.energy_available&&std::isfinite(w.orbitals[i].energy_hartree);proved_members.insert(i);
                    auto& name=out.canonical[i];name.status="verified_irrep";name.partner_block_id="canonical-irrep:"+std::to_string(b)+":"+row.label+":"+std::to_string(copy);name.partner_block_size=partners.size();}
                ++proved_counts[r];
                if(std::all_of(partners.begin(),partners.end(),[&](auto i){return separable_members.contains(i);}))++separable_counts[r];
                units.push_back(std::move(u));
            }
        }
        // Remove only certified complete source copies. A narrower energy
        // range for the remaining counts requires independent closure and an
        // exact character decomposition equal to the original minus these
        // copies. Failed residual evidence retains the original wide bounds.
        std::vector<std::size_t> remaining;for(auto i:members)if(!proved_members.contains(i))remaining.push_back(i);
        // Keep a previously certified enclosing remainder when finer accepted
        // source copies make their joint complement fail the strict closure
        // gate. Additional evidence must not replace a valid narrow energy
        // certificate with the entire original span's much wider interval.
        std::vector<std::size_t> enclosing;for(auto i:members)if(!separable_members.contains(i))enclosing.push_back(i);
        std::vector<std::pair<std::size_t,std::size_t>> enclosing_content;
        for(const auto& [r,copies]:content)if(copies>separable_counts[r])enclosing_content.emplace_back(r,copies-separable_counts[r]);
        const bool enclosing_verified=!enclosing.empty()&&decompose(characters(w,frame,enclosing),enclosing.size())==enclosing_content;
        if(enclosing_verified){low=std::numeric_limits<double>::infinity();high=-low;have_energy=true;
            for(auto i:enclosing){const double e=w.orbitals[i].energy_hartree;have_energy=have_energy&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);}}
        std::vector<std::pair<std::size_t,std::size_t>> expected;
        for(const auto& [r,copies]:content)if(copies>proved_counts[r])expected.emplace_back(r,copies-proved_counts[r]);
        const bool residual_verified=!remaining.empty()&&decompose(characters(w,frame,remaining),remaining.size())==expected;
        if(residual_verified){low=std::numeric_limits<double>::infinity();high=-low;have_energy=true;
            for(auto i:remaining){const double e=w.orbitals[i].energy_hartree;have_energy=have_energy&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);}}
        for(const auto& [r,copies]:expected){Unit count;count.members=residual_verified?remaining:members;count.irrep=known[r].label;count.scope=scope;count.copies=copies;count.copies_resolved=false;count.counting_only=true;
            if(have_energy){count.lower=low;count.upper=high;}units.push_back(std::move(count));}
    }
    assign_ordinals(std::move(units),out.canonical);
    // A large action-connected span can have an overly broad residual energy
    // interval even when the actual earlier source space is independently
    // closed. Energy proposes an ordering prefix only; full group action,
    // metric and integer character counts must certify it. This never creates
    // an irrep, partner membership, an energy degeneracy, or a rotated source.
    std::map<std::string,std::vector<std::size_t>> unordered_copies;
    for(std::size_t i=0;i<out.canonical.size();++i){const auto& name=out.canonical[i];
        if(name.verified&&!name.ordinal&&!name.partner_block_id.empty())unordered_copies[name.partner_block_id].push_back(i);}
    std::map<std::vector<std::size_t>,std::vector<std::pair<std::size_t,std::size_t>>> prefix_counts;
    std::map<std::vector<std::size_t>,std::string> prefix_failures;
    for(const auto& [id,members]:unordered_copies){const auto& first=out.canonical[members.front()];
        if(members.size()!=first.partner_block_size||first.representation_multiplicity!=1)continue;
        double low=std::numeric_limits<double>::infinity(),high=-low;bool finite=true;
        for(auto i:members){const double e=w.orbitals[i].energy_hartree;finite=finite&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);}
        const auto ordering_failure=[&](const std::string& status,const std::vector<std::size_t>& blockers){
            for(auto i:members){out.canonical[i].ordinal_status=status;out.canonical[i].ordinal_blocking_members=blockers;}};
        if(!finite){ordering_failure("source_energy_unavailable",members);continue;}
        std::vector<std::size_t> prefix;bool separated=true;
        for(std::size_t i=0;i<w.orbitals.size();++i){if(w.orbitals[i].spin!=w.orbitals[members.front()].spin)continue;
            const double e=w.orbitals[i].energy_hartree;if(!std::isfinite(e)){separated=false;ordering_failure("ordering_source_energy_unavailable",{i});break;}
            if(e<low-energy_order_tolerance)prefix.push_back(i);
            else if(e<=high+energy_order_tolerance&&std::find(members.begin(),members.end(),i)==members.end()&&
                (!out.canonical[i].verified||out.canonical[i].irrep==first.irrep)){separated=false;
                    ordering_failure(out.canonical[i].verified?"same_irrep_energy_order_unresolved":"overlapping_unresolved_source_irrep",{i});break;}
        }
        if(!separated)continue;
        const auto count_prefix=[&](const std::vector<std::size_t>& indices)->const std::vector<std::pair<std::size_t,std::size_t>>&{
            if(!indices.empty()&&!prefix_counts.contains(indices)){
                CharacterFailure failure;const auto values=characters(w,frame,indices,&failure);
                prefix_counts[indices]=decompose(values,indices.size());
                prefix_failures[indices]=failure.status.empty()?"character_counts_unavailable":failure.status;
            }
            return prefix_counts[indices];};
        auto counts=count_prefix(prefix);bool omitted_zero_spans=false;
        if(!prefix.empty()&&counts.empty()){
            // A cut through an unrelated irreducible space cannot obstruct
            // this irrep's count. Omit only complete containing spans already
            // certified by their integer character content to have zero copies
            // of the target; never discard a small individual-column weight.
            std::vector<std::size_t> relevant;
            for(auto i:prefix){const auto& content=out.canonical[i].containing_irreps;
                if(content.empty()||std::any_of(content.begin(),content.end(),[&](const auto& irrep){return irrep.irrep==first.irrep;}))relevant.push_back(i);}
            omitted_zero_spans=relevant.size()!=prefix.size();prefix=std::move(relevant);counts=count_prefix(prefix);
        }
        if(!prefix.empty()&&counts.empty()){ordering_failure("earlier_"+prefix_failures[prefix],prefix);continue;}
        std::size_t ordinal=1;for(const auto& [r,copies]:counts)if(known[r].label==first.irrep)ordinal+=copies;
        for(auto i:members){auto& name=out.canonical[i];name.ordinal=ordinal;
            name.ordinal_status="verified_complete_set_prefix_order";name.ordinal_blocking_members.clear();
            name.complete_set_ordinal_lower=ordinal;name.complete_set_ordinal_upper=ordinal;
            name.detail+="; ordinal certified by complete-group characters of "+std::to_string(prefix.size())+" earlier same-spin source columns; "+std::to_string(ordinal-1)+" complete earlier copies of "+first.irrep+"; actual prefix Gram/isometry/closure gates passed; source energy separation tolerance=2e-5 Ha";
            if(omitted_zero_spans)name.detail+="; omitted only independently closed containing spans with zero target-irrep copies";}
    }
    for(auto& name:out.canonical){if(name.verified||name.decomposition_verified)name.point_group=frame.group; retain_dominant_name(name);name.complete_set_ordinal=name.ordinal;name.ordinal_scope="complete canonical set";
        if(name.ordinal_status.empty())name.ordinal_status=name.verified?"original_partner_block_unverified":"source_irrep_not_verified";
        if(name.partner_status.empty())name.partner_status=name.verified?"original_partner_subset_not_verified":"source_irrep_not_verified";
        name.complete_set_ordinal_status=name.ordinal_status;}
    for(std::size_t i=0;i<out.canonical.size();++i){auto& name=out.canonical[i];if(name.verified||name.decomposition_verified)name.label=canonical_mo_display_label(w,i,&name);else if(open)name.label+=w.orbitals[i].spin==Spin::Beta?" [beta]":" [alpha]";}
    if(!salc)return out;
    Wavefunction side;side.atoms=w.atoms;side.shells=w.shells;side.primitives=w.primitives;side.basis_count=w.basis_count;side.ao_overlap=w.ao_overlap;side.orbitals.resize(salc->orbitals.size());
    for(std::size_t i=0;i<salc->orbitals.size();++i){const auto& orbital=salc->orbitals[i];auto& mo=side.orbitals[i];mo.energy_hartree=orbital.energy_hartree.value_or(std::numeric_limits<double>::quiet_NaN());mo.spin=orbital.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        bool valid=associated&&!orbital.terms.empty();mo.coefficients.assign(w.basis_count,0);
        for(const auto& term:orbital.terms){const auto* descriptor=nbo_orbital(data,term.orbital);if(!descriptor||descriptor->coefficients.size()!=w.basis_count||!std::isfinite(term.coefficient)){valid=false;break;}for(std::size_t a=0;a<w.basis_count;++a)mo.coefficients[a]+=term.coefficient*descriptor->coefficients[a];}
        if(!valid)mo.coefficients.clear();
    }
    std::vector<std::vector<std::size_t>> side_blocks;std::vector<std::size_t> side_subspaces;
    std::vector<bool> side_actual(salc->subspaces.size(),false),side_complete(salc->subspaces.size(),false);
    std::vector<CharacterFailure> side_failures(salc->subspaces.size());
    for(std::size_t s=0;s<salc->subspaces.size();++s){const auto& sub=salc->subspaces[s];bool have_coefficients=!sub.orbital_indices.empty();for(auto i:sub.orbital_indices)if(i>=side.orbitals.size()||side.orbitals[i].coefficients.size()!=w.basis_count)have_coefficients=false;
        side_actual[s]=have_coefficients;
        if(have_coefficients){const auto values=characters(side,frame,sub.orbital_indices,&side_failures[s]);
            side_complete[s]=!decompose(values,sub.orbital_indices.size()).empty();}
        const auto pieces=side_complete[s]?connected_blocks(side,frame,sub.orbital_indices):std::vector<std::vector<std::size_t>>{sub.orbital_indices};
        for(const auto& piece:pieces){side_blocks.push_back(piece);side_subspaces.push_back(s);}}

    std::vector<bool> used(salc->orbitals.size(),false);units.clear();
    for(std::size_t b=0;b<side_blocks.size();++b){const auto sub_index=side_subspaces[b];const auto& sub=salc->subspaces[sub_index];Unit u;u.scope=sub.fragment_id+":"+nbo_spin_name(sub.spin);u.detail="Fixed fragment subspace; "+frame.detail;
        bool indices_ok=!side_blocks[b].empty();for(auto i:side_blocks[b])if(i>=out.salc.size()||used[i]||salc->orbitals[i].fragment_id!=sub.fragment_id||salc->orbitals[i].spin!=sub.spin)indices_ok=false;if(!indices_ok)continue;u.members=side_blocks[b];
        const bool verified=frame.valid&&sub.symmetry_verified&&sub.dimension==sub.orbital_indices.size()&&std::isfinite(sub.closure_error)&&sub.closure_error<=character_tolerance&&std::isfinite(sub.orthogonality_error)&&sub.orthogonality_error<=metric_tolerance;
        const bool actual=side_actual[sub_index];
        std::vector<double> values;
        if(verified&&actual&&side_complete[sub_index])values=characters(side,frame,u.members);
        // A stored signature can stand in for omitted AO columns only for its
        // complete recorded span. Failed actual closure/metric evidence never
        // falls back to a stored signature that happens to look irreducible.
        if(verified&&!actual&&u.members==sub.orbital_indices)values=sub.characters;
        const auto content=decompose(values,u.members.size());
        for(auto i:u.members){used[i]=true;const auto& o=salc->orbitals[i];if(o.energy_hartree&&std::isfinite(*o.energy_hartree))u.energy+=*o.energy_hartree/double(u.members.size());else u.energy_available=false;out.salc[i].label="?";out.salc[i].detail=u.detail;
            auto& name=out.salc[i];name.containing_members=u.members;
            const auto& failure=side_failures[sub_index];
            name.status=!frame.valid?"operation_frame_unavailable":actual&&!failure.status.empty()?failure.status:!verified?"subspace_not_verified":values.empty()?(actual?"subspace_not_closed":"stored_characters_unavailable"):content.empty()?"character_decomposition_unavailable":actual?"mixed_source_member":"stored_span_unclassified";
            name.detail+="; "+name.status;
            if(!failure.detail.empty())name.detail+="; "+failure.detail;
            for(const auto& [r,copies]:content)name.containing_irreps.push_back({known[r].label,known[r].dimension,copies});
        }
        if(content.empty()){if(verified)certify_salc_energy_bounds(u,side,w,*salc);units.push_back(std::move(u));continue;}
        if(!actual){
            if(content.size()==1){const auto& [r,copies]=content.front();u.irrep=known[r].label;u.copies=copies;u.copies_resolved=copies==1;
                u.detail+="; complete stored characters verified; irrep occurrence count="+std::to_string(copies)+(u.copies_resolved?"; one ordinal shared by stored partners":"; repeated irrep span; individual copy numbering unresolved");
                for(auto i:u.members){auto& name=out.salc[i];name.status=u.copies_resolved?"verified_stored_irrep_span":"verified_stored_isotypic_span";
                    if(u.copies_resolved){name.partner_block_id="salc-irrep-block:"+std::to_string(b);name.partner_block_size=u.members.size();}}
                if(!u.copies_resolved)certify_salc_energy_bounds(u,side,w,*salc);
                units.push_back(std::move(u));
            }else{
                // Full stored reducible characters provide counts but cannot
                // identify any particular source column without AO evidence.
                for(const auto& [r,copies]:content){Unit count=u;count.irrep=known[r].label;count.copies=copies;count.copies_resolved=false;count.counting_only=true;certify_salc_energy_bounds(count,side,w,*salc);units.push_back(std::move(count));}
            }
            continue;
        }
        retain_decomposition(side,frame,u.members,known,out.salc,"salc");
        std::vector<std::size_t> proved_counts(known.size(),0);std::set<std::size_t> proved_members;
        for(const auto& [r,copies]:content){const auto& row=known[r];std::vector<std::size_t> pure;
            const auto losses=projection_residuals(side,frame,u.members,row.dimension,row.norm,row.values);
            for(std::size_t c=0;c<losses.size();++c){auto& name=out.salc[u.members[c]];retain_projection_residual(name,losses[c]);if(pure_projection(losses[c]))pure.push_back(u.members[c]);}
            Unit evidence=u;certify_salc_energy_bounds(evidence,side,w,*salc);
            for(auto i:pure){auto& name=out.salc[i];name.verified=true;name.irrep=row.label;name.label=orbital_label(row.label);name.representation_multiplicity=copies;
                name.status="verified_isotypic_member";name.detail=(sub.basis_kind=="derived_salc"?"Mapped derived SALC column":"Unchanged producer column")+std::string(" lies in verified S-metric isotypic projection; independent occurrence membership unresolved; ")+evidence.detail;}
            const auto source_copies=verified_source_copies(side,frame,pure,row.dimension,row.norm,row.values,copies,nullptr,&out.salc);
            for(std::size_t copy=0;copy<source_copies.size();++copy){const auto& partners=source_copies[copy];
                if(std::any_of(partners.begin(),partners.end(),[&](auto i){return proved_members.contains(i);}))continue;
                Unit named;named.members=partners;named.irrep=row.label;named.scope=u.scope;
                named.detail="Verified irreducible source SALC span under every operation; "+frame.detail;
                for(auto i:partners){const auto& energy=salc->orbitals[i].energy_hartree;if(energy&&std::isfinite(*energy))named.energy+=*energy/double(partners.size());else named.energy_available=false;proved_members.insert(i);
                    auto& name=out.salc[i];name.status="verified_irrep";name.partner_block_id="salc-irrep:"+std::to_string(b)+":"+row.label+":"+std::to_string(copy);name.partner_block_size=partners.size();}
                ++proved_counts[r];
                units.push_back(std::move(named));
            }
        }
        std::vector<std::size_t> remaining;for(auto i:u.members)if(!proved_members.contains(i))remaining.push_back(i);
        std::vector<std::pair<std::size_t,std::size_t>> expected;
        for(const auto& [r,copies]:content)if(copies>proved_counts[r])expected.emplace_back(r,copies-proved_counts[r]);
        const bool residual_verified=!remaining.empty()&&decompose(characters(side,frame,remaining),remaining.size())==expected;
        for(const auto& [r,copies]:expected){Unit count=u;count.members=residual_verified?remaining:u.members;count.irrep=known[r].label;count.copies=copies;count.copies_resolved=false;count.counting_only=true;
            certify_salc_energy_bounds(count,side,w,*salc);units.push_back(std::move(count));}
    }
    for(std::size_t i=0;i<out.salc.size();++i)if(!used[i]){const auto& o=salc->orbitals[i];Unit u;u.members={i};u.scope=o.fragment_id+":"+nbo_spin_name(o.spin);u.energy=o.energy_hartree.value_or(0);u.energy_available=o.energy_hartree&&std::isfinite(*o.energy_hartree);units.push_back(u);out.salc[i].label="?";out.salc[i].status="subspace_unavailable";out.salc[i].detail="No validated containing SALC subspace";}
    // A fragment's span need not be closed under the whole molecular group.
    // Its actual vector still has a global composition if a complete same-spin
    // canonical basis independently closes and reconstructs every component.
    // This adds composition evidence only; local names/counts are untouched.
    if(frame.valid&&associated)for(const auto spin:{NboSpin::Total,NboSpin::Alpha,NboSpin::Beta}){
        if(spin==NboSpin::Total&&(open||w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin))continue;
        std::vector<std::size_t> targets,source_members;
        for(std::size_t i=0;i<out.salc.size();++i)if(!out.salc[i].decomposition_verified&&
            salc->orbitals[i].spin==spin&&side.orbitals[i].coefficients.size()==w.basis_count)targets.push_back(i);
        // A producer-declared shared canonical spatial set represents both
        // spin channels. This is only a coordinate expansion of the actual
        // SALC vector; its source spin and local identity remain unchanged.
        const bool shared_spatial=w.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared&&!open;
        for(std::size_t i=0;i<w.orbitals.size();++i)if(w.orbitals[i].spin==(spin==NboSpin::Beta&&!shared_spatial?Spin::Beta:Spin::Alpha))source_members.push_back(i);
        if(!targets.empty())retain_decomposition(side,frame,targets,known,out.salc,"canonical",&w,source_members);
    }
    assign_ordinals(std::move(units),out.salc);
    for(std::size_t i=0;i<out.salc.size();++i)if(out.salc[i].verified&&!out.salc[i].ordinal)out.salc[i].label+=" [SALC "+std::to_string(i+1)+"]";
    for(std::size_t i=0;i<out.salc.size();++i)if(open && !salc->orbitals[i].spatial_spin)out.salc[i].label+=salc->orbitals[i].spin==NboSpin::Beta?" [beta]":salc->orbitals[i].spin==NboSpin::Alpha?" [alpha]":" [total]";
    for(std::size_t i=0;i<out.salc.size();++i){auto& name=out.salc[i];if(name.verified||name.decomposition_verified)name.point_group=frame.group; retain_dominant_name(name);
        if(name.approximate_dominant_label){name.label="\xE2\x89\x88"+orbital_label(name.dominant_irrep)+" [SALC "+std::to_string(i+1)+"]";
            if(salc->orbitals[i].spin!=NboSpin::Total)name.label+=" ["+std::string(nbo_spin_name(salc->orbitals[i].spin))+"]";}
        name.complete_set_ordinal=name.ordinal;name.ordinal_scope="complete fragment spin set";
        if(name.ordinal_status.empty())name.ordinal_status=name.verified?"original_partner_block_unverified":"source_irrep_not_verified";
        if(name.partner_status.empty())name.partner_status=name.partner_block_id.empty()?"original_partner_subset_not_verified":"verified_stored_partner_span";
        name.complete_set_ordinal_status=name.ordinal_status;}
    (void)data; // Identity is immutable and belongs to the caller's attachment.
    return out;
}

NboAomoNames nbo_aomo_names_for_view(const Wavefunction& w,const NboAomoNames& source,
    const std::vector<std::size_t>& canonical_indices,const std::vector<std::size_t>& salc_indices,
    const NboSalcModel* model,const std::string& scope,
    const std::vector<double>* display_energies){
    auto out=source;
    const auto number=[&](std::vector<NboAomoName>& names,const std::vector<std::size_t>& indices,bool side){
        struct Occurrence {std::string family,block;std::vector<std::size_t> members,visible;double energy=0;bool available=true;};
        std::map<std::string,Occurrence> groups;std::set<std::size_t> selected(indices.begin(),indices.end());
        const auto valid=[&](std::size_t i){return i<names.size()&&(side?model&&i<model->orbitals.size():i<w.orbitals.size());};
        const auto energy=[&](std::size_t i)->std::optional<double> {
            if(side)return model->orbitals[i].energy_hartree;
            if(display_energies)return i<display_energies->size()?
                std::optional<double>((*display_energies)[i]):std::nullopt;
            return w.orbitals[i].energy_hartree;
        };
        const auto family=[&](std::size_t i){return side?model->orbitals[i].fragment_id+":"+nbo_spin_name(model->orbitals[i].spin):
            w.orbitals[i].spin==Spin::Beta?std::string("canonical beta"):std::string("canonical alpha");};
        std::vector<std::size_t> rows;for(auto i:selected)if(valid(i)){
            auto& name=names[i];name.ordinal=0;name.visible_partner_count=0;
            name.ordinal_scope=scope+(name.point_group=="SO(3)"?"; atomic radial-copy order; not principal n":"");
            name.ordinal_status=!name.verified?"source_irrep_not_verified":name.partner_block_id.empty()?"original_partner_block_unverified":"certified_copy_pending_view_order";
            rows.push_back(i);
        }
        std::stable_sort(rows.begin(),rows.end(),[&](auto i,auto j){const auto a=energy(i),b=energy(j);const bool av=a&&std::isfinite(*a),bv=b&&std::isfinite(*b);
            if(av!=bv)return av;if(av&&*a!=*b)return *a<*b;return i<j;});
        for(std::size_t j=0;j<rows.size();++j)names[rows[j]].view_row_ordinal=j+1;
        // Read every certified partner from the immutable complete model. Hidden
        // members retain their contribution to the copy's trace/size sort key.
        for(std::size_t i=0;i<names.size();++i){if(!valid(i))continue;const auto& name=names[i];
            if(!name.verified||name.irrep.empty()||name.partner_block_id.empty()||name.representation_multiplicity!=1)continue;
            auto& group=groups[family(i)+":"+name.irrep+":"+name.partner_block_id];group.family=family(i)+":"+name.irrep;group.block=name.partner_block_id;group.members.push_back(i);
            if(selected.contains(i))group.visible.push_back(i);
            const auto e=energy(i);group.available=group.available&&e&&std::isfinite(*e);if(e&&std::isfinite(*e))group.energy+=*e;
        }
        std::vector<Occurrence> ordered;
        for(auto& [key,g]:groups){if(g.visible.empty())continue;
            for(auto i:g.visible)names[i].visible_partner_count=g.visible.size();
            if(g.members.size()!=names[g.members.front()].partner_block_size){for(auto i:g.visible)names[i].ordinal_status="complete_partner_evidence_inconsistent";continue;}
            if(!g.available){for(auto i:g.visible)names[i].ordinal_status="nonquantitative_copy_entry";continue;}
            g.energy/=double(g.members.size());ordered.push_back(std::move(g));}
        std::stable_sort(ordered.begin(),ordered.end(),[](const auto& a,const auto& b){if(a.family!=b.family)return a.family<b.family;if(a.energy!=b.energy)return a.energy<b.energy;return a.block<b.block;});
        std::map<std::string,std::size_t> counters;
        for(const auto& g:ordered){const auto ordinal=++counters[g.family];for(auto i:g.visible){auto& name=names[i];name.ordinal=ordinal;name.ordinal_status="display_order_convention";
            name.detail+="; display ordinal in "+scope+": certified occurrences ordered by full-copy mean energy; exact ties use stable copy identity (display convention); visible members="+std::to_string(g.visible.size())+"/"+std::to_string(g.members.size());}}
        for(auto i:rows){auto& name=names[i];
            if(!side)name.label=canonical_mo_display_label(w,i,&name);
            else if(name.approximate_dominant_label){name.label="\xE2\x89\x88"+orbital_label(name.dominant_irrep)+" [SALC "+std::to_string(i+1)+"]";
                if(model->orbitals[i].spin!=NboSpin::Total)name.label+=" ["+std::string(nbo_spin_name(model->orbitals[i].spin))+"]";}
            else if(name.verified){name.label=numbered_orbital_label(name);if(!name.ordinal)name.label+=" [SALC "+std::to_string(i+1)+"]";
                if(model->orbitals[i].spin!=NboSpin::Total)name.label+=" ["+std::string(nbo_spin_name(model->orbitals[i].spin))+"]";}
        }
    };
    number(out.canonical,canonical_indices,false);number(out.salc,salc_indices,true);return out;
}

std::shared_ptr<const NboAomoNames> canonical_mo_names(const Wavefunction& w){
    struct Cache {
        const Wavefunction* wavefunction=nullptr;
        const MolecularOrbital* orbitals=nullptr;
        const Atom* atoms=nullptr;
        const DerivedOrbitalSymmetryAssignment* assignments=nullptr;
        std::size_t orbital_count=0,atom_count=0,assignment_count=0,revision=0;
        std::string group,enrichment;
        std::shared_ptr<const NboAomoNames> names;
    };
    static thread_local Cache cache;
    if(!cache.names||cache.revision!=canonical_name_revision||cache.wavefunction!=&w||cache.orbitals!=w.orbitals.data()||
       cache.atoms!=w.atoms.data()||cache.assignments!=w.derived_orbital_symmetry_assignments.data()||
       cache.orbital_count!=w.orbitals.size()||cache.atom_count!=w.atoms.size()||
       cache.assignment_count!=w.derived_orbital_symmetry_assignments.size()||
       cache.group!=w.point_group_detected||cache.enrichment!=w.enrichment_source){
        cache.wavefunction=&w;cache.orbitals=w.orbitals.data();cache.atoms=w.atoms.data();
        cache.assignments=w.derived_orbital_symmetry_assignments.data();
        cache.orbital_count=w.orbitals.size();cache.atom_count=w.atoms.size();
        cache.assignment_count=w.derived_orbital_symmetry_assignments.size();
        cache.group=w.point_group_detected;cache.enrichment=w.enrichment_source;
        cache.revision=canonical_name_revision;
        OpenProfile profile;
        cache.names=std::make_shared<const NboAomoNames>(build_nbo_aomo_names(w,NboIntegration{},nullptr));
        profile.stage("standalone-orbital-names");
    }
    return cache.names;
}

void invalidate_canonical_mo_names_cache(){++canonical_name_revision;}

std::string canonical_mo_source_label(const Wavefunction& w,std::size_t index){
    if(index>=w.orbitals.size())return "MO ?";
    const auto& mo=w.orbitals[index];
    const bool source=mo.source_orbital_index!=std::numeric_limits<std::size_t>::max();
    std::string label=std::string(source?"MO ":"MO [list] ")+std::to_string((source?mo.source_orbital_index:index)+1);
    const bool explicit_spin=w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin||
        std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& o){return o.spin==Spin::Beta;});
    if(explicit_spin)label+=mo.spin==Spin::Beta?" [beta]":" [alpha]";
    return label;
}

std::string canonical_mo_display_label(const Wavefunction& w,std::size_t index,const NboAomoName* name){
    if(index>=w.orbitals.size())return "MO ?";
    const auto standalone=name?std::shared_ptr<const NboAomoNames>{}:canonical_mo_names(w);
    if(!name&&standalone&&index<standalone->canonical.size())name=&standalone->canonical[index];
    if(name&&name->verified&&name->ordinal&&!name->irrep.empty()){
        auto label=numbered_orbital_label(*name);
        const bool explicit_spin=w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin||
            std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& o){return o.spin==Spin::Beta;});
        if(explicit_spin)
            label+=w.orbitals[index].spin==Spin::Beta?" [beta]":" [alpha]";
        return label;
    }
    if(name&&name->verified&&!name->irrep.empty()){
        auto source=canonical_mo_source_label(w,index);std::string suffix;
        for(const std::string spin:{" [alpha]"," [beta]"})if(source.ends_with(spin)){source.resize(source.size()-spin.size());suffix=spin;break;}
        return orbital_label(name->irrep)+" ["+source+"]"+suffix;
    }
    if(name&&name->approximate_dominant_label&&!name->dominant_irrep.empty()){
        return "\xE2\x89\x88"+orbital_label(name->dominant_irrep)+" ["+canonical_mo_source_label(w,index)+"]";}
    return canonical_mo_source_label(w,index);
}

std::string canonical_mo_current_irrep(const Wavefunction& w,std::size_t index,const NboAomoName* name){
    if(index>=w.orbitals.size())return "?";
    const auto standalone=name?std::shared_ptr<const NboAomoNames>{}:canonical_mo_names(w);
    if(!name&&standalone&&index<standalone->canonical.size())name=&standalone->canonical[index];
    return name?orbital_irrep_display_label(*name):"?";
}

std::string orbital_irrep_display_label(const NboAomoName& name){
    return name.verified&&!name.irrep.empty()?orbital_label(name.irrep):name.approximate_dominant_label?"\xE2\x89\x88"+orbital_label(name.dominant_irrep):"?";
}

std::string serialize_orbital_name_json(const NboAomoName& name){
    const auto quote=[](const std::string& text){std::ostringstream s;s<<'"';for(unsigned char c:text){
        if(c=='"'||c=='\\')s<<'\\'<<char(c);else if(c<32)s<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c);else s<<char(c);}s<<'"';return s.str();};
    std::ostringstream out;out<<std::setprecision(17)<<"{\"label\":"<<quote(name.label)
        <<",\"irrep\":"<<quote(name.irrep)<<",\"point_group\":"<<quote(name.point_group)
        <<",\"ordinal\":"<<name.ordinal<<",\"ordinal_scope\":"<<quote(name.ordinal_scope)
        <<",\"ordinal_energy_tolerance_hartree\":"<<energy_order_tolerance
        <<",\"view_row_ordinal\":"<<name.view_row_ordinal<<",\"visible_partner_count\":"<<name.visible_partner_count
        <<",\"approximate_dominant_label\":"<<(name.approximate_dominant_label?"true":"false")
        <<",\"dominant_irrep\":"<<quote(name.dominant_irrep)<<",\"dominant_weight\":"<<name.dominant_weight
        <<",\"complete_set_ordinal\":"<<name.complete_set_ordinal
        <<",\"ordinal_status\":"<<quote(name.ordinal_status)
        <<",\"complete_set_ordinal_status\":"<<quote(name.complete_set_ordinal_status)
        <<",\"verified\":"<<(name.verified?"true":"false")<<",\"status\":"<<quote(name.status)
        <<",\"representation_multiplicity\":"<<name.representation_multiplicity
        <<",\"partner_block_id\":"<<quote(name.partner_block_id)<<",\"partner_block_size\":"<<name.partner_block_size
        <<",\"partner_status\":"<<quote(name.partner_status)
        <<",\"partner_failure_quantity\":"<<quote(name.partner_failure_quantity)
        <<",\"partner_failure_value\":";
    const auto optional_number=[&](const std::optional<double>& value){if(value&&std::isfinite(*value))out<<*value;else out<<"null";};
    optional_number(name.partner_failure_value);out<<",\"partner_failure_limit\":";optional_number(name.partner_failure_limit);
    out<<",\"complete_set_ordinal_bounds\":";
    if(name.complete_set_ordinal_lower&&name.complete_set_ordinal_upper)
        out<<"{\"lower\":"<<*name.complete_set_ordinal_lower<<",\"upper\":"<<*name.complete_set_ordinal_upper
           <<",\"meaning\":\"conservative possible ranks; not a chosen ordinal or a claim that every interior rank is attainable\"}";
    else out<<"null";
    out<<",\"ordinal_blocking_members\":[";for(std::size_t i=0;i<name.ordinal_blocking_members.size();++i){if(i)out<<',';out<<name.ordinal_blocking_members[i];}
    out<<"],\"projection_residual_squared\":";
    if(name.projection_residual)out<<*name.projection_residual;else out<<"null";
    out<<",\"decomposition_verified\":"<<(name.decomposition_verified?"true":"false")
        <<",\"decomposition_status\":"<<quote(name.decomposition_status)
        <<",\"decomposition_reconstruction_residual\":";
    const auto number=[&](double value){if(std::isfinite(value))out<<value;else out<<"null";};
    number(name.decomposition_reconstruction_residual);
    out<<",\"decomposition_orthogonality_error\":";number(name.decomposition_orthogonality_error);
    out<<",\"decomposition_weight_sum_error\":";number(name.decomposition_weight_sum_error);
    out<<",\"decomposition_thresholds\":{\"source_residual_squared\":"<<character_tolerance*character_tolerance
        <<",\"metric_inner_product\":"<<metric_tolerance<<",\"metric_weight_sum\":"<<metric_tolerance<<'}';
    out<<",\"component_source_kind\":"<<quote(name.component_source_kind)<<",\"component_source_members\":[";
    for(std::size_t i=0;i<name.component_source_members.size();++i){if(i)out<<',';out<<name.component_source_members[i];}
    out<<"],\"components\":[";
    for(std::size_t r=0;r<name.components.size();++r){if(r)out<<',';const auto& component=name.components[r];
        out<<"{\"irrep\":"<<quote(component.irrep)<<",\"dimension\":"<<component.dimension<<",\"weight\":"<<component.weight<<",\"source_coefficients\":[";
        for(std::size_t i=0;i<component.source_coefficients.size();++i){if(i)out<<',';out<<component.source_coefficients[i];}out<<"]}";}
    out<<']';
    out<<",\"containing_members\":[";for(std::size_t i=0;i<name.containing_members.size();++i){if(i)out<<',';out<<name.containing_members[i];}
    out<<"],\"containing_irreps\":[";for(std::size_t i=0;i<name.containing_irreps.size();++i){if(i)out<<',';const auto& r=name.containing_irreps[i];
        out<<"{\"irrep\":"<<quote(r.irrep)<<",\"dimension\":"<<r.dimension<<",\"multiplicity\":"<<r.multiplicity<<'}';}
    out<<"],\"detail\":"<<quote(name.detail)<<'}';return out.str();
}
std::string serialize_orbital_names_json(const NboAomoNames& names){
    std::ostringstream out;out<<"{\"schema\":\"cov.orbital.display-names.v2\",\"canonical\":[";
    const auto rows=[&](const std::vector<NboAomoName>& entries){for(std::size_t i=0;i<entries.size();++i){if(i)out<<',';
        auto row=serialize_orbital_name_json(entries[i]);out<<"{\"index\":"<<i<<','<<row.substr(1);}};
    rows(names.canonical);out<<"],\"salc\":[";rows(names.salc);out<<"]}";return out.str();
}
} // namespace cov::ui
