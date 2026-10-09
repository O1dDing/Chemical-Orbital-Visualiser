#include "cov/nbo_salc.hpp"
#include "cov/ao_angular_basis.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
void near(double a,double b,const char* message,double tol=1e-8){require(std::abs(a-b)<tol,message);}
cov::NboMatrix matrix(const char* kind,std::vector<double> v,cov::NboSpin spin=cov::NboSpin::Total){cov::NboMatrix m;m.kind=kind;m.spin=spin;m.rows=m.columns=3;m.values=std::move(v);return m;}
struct Fixture {cov::Wavefunction w;cov::NboIntegration i;};
Fixture fixture(){Fixture f;auto& w=f.w;auto& i=f.i;const double h=std::sqrt(.5);
    w.basis_count=3;w.atoms={{"O",8,0,0,0,8},{"H",1,-1,0,1,1},{"H",1,1,0,1,1}};w.shells={{0,0,1,0,0,0},{1,1,1,1,0,0},{2,2,1,2,0,0}};w.primitives={{1,1},{1,1},{1,1}};w.gaussian_ao_transform={{0,1,1},{1,1,1},{2,1,1}};w.ao_overlap={1,0,0,0,1,0,0,0,1};w.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;
    const std::vector<std::vector<double>> coefficients={{1,0,0},{0,h,h},{0,h,-h}};
    for(std::size_t j=0;j<3;++j){cov::MolecularOrbital mo;mo.coefficients=mo.gaussian_source_coefficients=coefficients[j];mo.energy_hartree=j==0?-1:j==1?-.5:.3;mo.occupation=j<2?2:0;mo.occupation_provenance=cov::DataProvenance::Producer;mo.source_orbital_index=j;w.orbitals.push_back(mo);cov::NboOrbitalDescriptor d;d.ref={cov::NboOrbitalKind::NAO,cov::NboSpin::Total,j};d.id="NAO:"+std::to_string(j);d.label=j==0?"O 2s":"H 1s "+std::to_string(j);d.atoms={j};d.coefficients.assign(3,0);d.coefficients[j]=1;d.orthonormal_basis=true;d.metric_norm2=1;d.occupation=j==0?2:1;i.orbitals.push_back(d);cov::NboNao row;row.id=j+1;row.atom=j+1;row.type=j==0?"Val(2s)":"Val(1s)";row.angular="s";row.occupation=*d.occupation;i.dataset.naos.push_back(row);}
    cov::NboArchive ar;ar.atoms=w.atoms;ar.basis_count=3;ar.density_is_bond_order=true;ar.matrices={matrix("OVERLAP",w.ao_overlap),matrix("FOCK",{-1,0,0,0,-.1,-.4,0,-.4,-.1}),matrix("DENSITY",{2,0,0,0,1,1,0,1,1})};i.dataset.archive=ar;i.id="salc-synthetic";i.dataset.association.compatible=true;i.dataset.association.gaussian_row={0,1,2};i.dataset.association.coefficient_scale={1,1,1};cov::NboCanonicalEvidence ev;ev.direct_fchk_coefficients=true;ev.density_verified=true;i.dataset.association.canonical_evidence.push_back(ev);i.canonical_fingerprint=cov::nbo_canonical_fingerprint(w);return f;
}
void refresh(Fixture& f){f.i.canonical_fingerprint=cov::nbo_canonical_fingerprint(f.w);}
std::vector<double> salc_link_weights(const cov::NboSalcModel& m,std::size_t mo){std::vector<double> out;for(const auto& l:m.links)if(l.canonical_index==mo)out.push_back(l.weight);std::sort(out.begin(),out.end());return out;}
void same_weights(const cov::NboSalcModel& a,const cov::NboSalcModel& b){for(std::size_t mo=0;mo<3;++mo){const auto x=salc_link_weights(a,mo),y=salc_link_weights(b,mo);require(x.size()==y.size(),"weight count changed");for(std::size_t j=0;j<x.size();++j)near(x[j],y[j],"rotation/permutation changed projection weights");}}
double polynomial(const cov::ao_angular::Polynomial& p,const std::array<double,3>& r){double value=0;for(const auto& t:p)value+=t.coefficient*std::pow(r[0],t.powers[0])*std::pow(r[1],t.powers[1])*std::pow(r[2],t.powers[2]);return value;}
void angular_tests(){
    const double angle=.619,c=std::cos(angle),s=std::sin(angle),d=1-c;const std::array<double,3> u={1/std::sqrt(14.0),2/std::sqrt(14.0),3/std::sqrt(14.0)};
    cov::SymmetryOperation op;op.atom_permutation={0};op.matrix={c+u[0]*u[0]*d,u[0]*u[1]*d-u[2]*s,u[0]*u[2]*d+u[1]*s,u[1]*u[0]*d+u[2]*s,c+u[1]*u[1]*d,u[1]*u[2]*d-u[0]*s,u[2]*u[0]*d-u[1]*s,u[2]*u[1]*d+u[0]*s,c+u[2]*u[2]*d};
    for(int l=0;l<=4;++l)for(int pure=0;pure<=1;++pure){cov::Wavefunction w;w.atoms={{"C",6,0,0,0,6}};w.shells={{0,0,1,0,static_cast<std::uint8_t>(l),static_cast<std::uint8_t>(pure)}};w.primitives={{1,1}};const auto n=cov::shell_basis_count(w.shells[0]);w.basis_count=n;std::vector<double> eye(n*n);for(std::size_t i=0;i<n;++i)eye[i*n+i]=1;const auto t=cov::apply_orbital_symmetry_operation(w,op,eye,n);require(t.size()==n*n,"s through g angular operation unavailable");std::vector<cov::ao_angular::Polynomial> basis;for(std::size_t i=0;i<n;++i){if(pure)basis.push_back(cov::ao_angular::solid_harmonic(l,int(i)));else{const auto p=cov::ao_angular::cartesian_components(l)[i];const double norm=std::sqrt(cov::ao_angular::odd_factorial(2*p[0]-1)*cov::ao_angular::odd_factorial(2*p[1]-1)*cov::ao_angular::odd_factorial(2*p[2]-1));basis.push_back({{p,1/norm}});}}
        for(const auto& point:std::vector<std::array<double,3>>{{.27,-.61,.83},{-.79,.41,.52},{1.1,-.3,-.2}}){std::array<double,3> pulled{};for(int i=0;i<3;++i)for(int k=0;k<3;++k)pulled[i]+=op.matrix[3*k+i]*point[k];for(std::size_t j=0;j<n;++j){double result=0;for(std::size_t i=0;i<n;++i)result+=polynomial(basis[i],point)*t[i*n+j];near(result,polynomial(basis[j],pulled),"independent polynomial angular transformation mismatch",2e-8);}}
    }
}
}
int main(){try{
    auto scoped=fixture();scoped.w.total_density_packed={2,0,1,0,1,1};
    auto scope=cov::analyse_nbo_electronic_symmetry_scope(scoped.w);
    require(scope.naming_scope_verified&&!scope.reduced&&scope.geometry_group=="C2v"&&scope.naming_group=="C2v","Symmetric density lost the complete nuclear group");
    scoped.w.total_density_packed={2,0,1.1,0,0,0.9};
    scope=cov::analyse_nbo_electronic_symmetry_scope(scoped.w);
    require(scope.naming_scope_verified&&scope.reduced&&scope.geometry_group=="C2v"&&scope.naming_group=="Cs"&&scope.operations.size()==2&&scope.geometry_operations.size()==4,"Electronic density subgroup must preserve nuclear identity and carry actual closed operations");
    scoped.w.alpha_electrons=2;scoped.w.beta_electrons=1;scoped.w.spin_density_packed.clear();
    require(!cov::analyse_nbo_electronic_symmetry_scope(scoped.w).naming_scope_verified,"Missing open-shell spin density must not certify a shared electronic scope");
    auto independent_density=fixture();
    independent_density.i.dataset.archive->density_is_bond_order=true;
    independent_density.w.orbitals[0].occupation=0;
    refresh(independent_density);
    const auto producer_density_model=cov::build_nbo_salc_model(independent_density.w,independent_density.i);
    require(producer_density_model.orbitals.front().occupation.has_value(),"verified producer density occupation lost");
    near(*producer_density_model.orbitals.front().occupation,2,
         "SALC overwrote independent producer density with unverified canonical occupation");
    angular_tests();
    auto f=fixture();const auto hash=f.i.canonical_fingerprint;const auto m=cov::build_nbo_salc_model(f.w,f.i);require(m.available&&m.group_verified,"fixed SALC unavailable");require(m.orbitals.size()==3&&m.fragments.size()==2,"wrong full-space or fragment dimension");require(m.fragments[0].atoms==std::vector<std::size_t>{0}&&m.fragments[0].side==0,"central atom layout");std::size_t combinations=0;for(const auto& o:m.orbitals)if(o.symmetry_adapted){++combinations;require(o.terms.size()==2,"two equivalent s atoms must produce actual combinations");for(const auto& t:o.terms)near(t.coefficient*t.coefficient,.5,"SALC normalization");require(o.energy_hartree.has_value(),"missing verified SALC expectation");}require(combinations==2,"expected two fixed peripheral SALCs");require(m.energies.size()==1&&m.energies[0].available,"same-operator energy gate");for(const auto& c:m.coverage){near(c.weight_sum,1,"raw coverage");near(c.residual_norm,0,"signed reconstruction");}require(cov::nbo_canonical_fingerprint(f.w)==hash,"canonical mutation");
    const auto again=cov::build_nbo_salc_model(f.w,f.i);require(cov::serialize_nbo_salc_json(again)==cov::serialize_nbo_salc_json(m),"unstable cached model identity/gauge");
    for(std::size_t j=0;j<m.orbitals.size();++j){const auto selection=cov::nbo_salc_selection(m,j);require(selection.mode==cov::NboSelectionMode::Combination&&!selection.normalize,"SALC selection changed amplitudes");auto v=cov::make_nbo_selection_view(f.i,f.w,selection);require(v.available,"SALC terms not renderable in existing selector");near(v.metric_norm2[0],1,"selected SALC is not normalized");}
    auto phase=fixture();for(auto& x:phase.w.orbitals[1].coefficients)x=-x;for(auto& x:phase.w.orbitals[1].gaussian_source_coefficients)x=-x;refresh(phase);auto phased=cov::build_nbo_salc_model(phase.w,phase.i);same_weights(m,phased);for(std::size_t j=0;j<m.orbitals.size();++j){require(m.orbitals[j].id==phased.orbitals[j].id,"SALC depends on canonical phase");require(m.orbitals[j].terms.size()==phased.orbitals[j].terms.size(),"SALC gauge depends on MO");}for(std::size_t j=0;j<m.links.size();++j)if(m.links[j].canonical_index==1)near(m.links[j].coefficient,-phased.links[j].coefficient,"signed link phase lost");
    auto rotated=fixture();const double c=std::cos(.721),s=std::sin(.721);for(auto& at:rotated.w.atoms){const auto x=at.x,z=at.z;at.x=c*x+s*z;at.z=-s*x+c*z;}refresh(rotated);same_weights(m,cov::build_nbo_salc_model(rotated.w,rotated.i));
    auto permuted=fixture();std::swap(permuted.w.atoms[1],permuted.w.atoms[2]);refresh(permuted);same_weights(m,cov::build_nbo_salc_model(permuted.w,permuted.i));
    auto degenerate=fixture();degenerate.w.orbitals[1].energy_hartree=degenerate.w.orbitals[2].energy_hartree=-.5;degenerate.w.orbitals[1].occupation=degenerate.w.orbitals[2].occupation=1;degenerate.i.dataset.archive->matrices[2].values={2,0,0,0,1,0,0,0,1};degenerate.i.dataset.archive->matrices[1].values={-1,0,0,0,-.5,0,0,0,-.5};refresh(degenerate);const auto dm=cov::build_nbo_salc_model(degenerate.w,degenerate.i);auto mixed=degenerate;const double ca=std::cos(.417),sa=std::sin(.417);for(std::size_t j=0;j<3;++j){const double a=degenerate.w.orbitals[1].coefficients[j],b=degenerate.w.orbitals[2].coefficients[j];mixed.w.orbitals[1].coefficients[j]=ca*a+sa*b;mixed.w.orbitals[2].coefficients[j]=-sa*a+ca*b;}for(auto& mo:mixed.w.orbitals)mo.gaussian_source_coefficients=mo.coefficients;refresh(mixed);const auto mm=cov::build_nbo_salc_model(mixed.w,mixed.i);require(mm.energies[0].available,"degenerate canonical rotation rejected");for(std::size_t side=0;side<dm.orbitals.size();++side){require(dm.orbitals[side].id==mm.orbitals[side].id,"SALC identity depends on degenerate canonical gauge");double before=0,after=0;for(const auto& link:dm.links)if(link.side_index==side&&link.canonical_index>0)before+=link.weight;for(const auto& link:mm.links)if(link.side_index==side&&link.canonical_index>0)after+=link.weight;near(before,after,"degenerate-space summed weight changed");}
    auto missing=fixture();std::erase_if(missing.i.dataset.archive->matrices,[](const auto& x){return x.kind=="FOCK";});auto noenergy=cov::build_nbo_salc_model(missing.w,missing.i);require(noenergy.available&&!noenergy.energies[0].available,"missing Fock fabricated energy or destroyed basis");for(const auto& o:noenergy.orbitals)require(!o.energy_hartree,"missing Fock has side energy");
    auto wrong=fixture();wrong.i.dataset.archive->matrices[1].values[0]+=.2;auto rejected=cov::build_nbo_salc_model(wrong.w,wrong.i);require(rejected.available&&!rejected.energies[0].available&&rejected.energies[0].canonical_residual>.1,"wrong Fock accepted");
    // Printed NAOs can pass the declared metric gate without being exactly
    // orthonormal. Closure must measure the field residual, not sqrt(metric error).
    auto rounded=fixture();for(auto& d:rounded.i.orbitals)for(auto& x:d.coefficients)x*=1-1e-7;
    auto rm=cov::build_nbo_salc_model(rounded.w,rounded.i);require(rm.available&&rm.orbitals.size()==3,"rounded NAO basis lost");
    std::size_t rounded_salc=0;for(const auto& o:rm.orbitals)if(o.symmetry_adapted)++rounded_salc;require(rounded_salc==2,"NAO rounding falsely destroys closed SALC family");
    auto broken=fixture();for(auto& d:broken.i.orbitals)for(auto& x:d.coefficients)x*=.99;require(!cov::build_nbo_salc_model(broken.w,broken.i).available,"rank/metric-invalid NAO normalization accepted");
    auto rotate_h=[](Fixture& x){const auto a=x.w.orbitals[1].coefficients,b=x.w.orbitals[2].coefficients;const double c=std::cos(.31),s=std::sin(.31);for(std::size_t j=0;j<3;++j){x.w.orbitals[1].coefficients[j]=c*a[j]+s*b[j];x.w.orbitals[2].coefficients[j]=-s*a[j]+c*b[j];}for(auto& mo:x.w.orbitals)mo.gaussian_source_coefficients=mo.coefficients;refresh(x);};
    auto density_broken=fixture();density_broken.w.orbitals[1].energy_hartree=density_broken.w.orbitals[2].energy_hartree=-.5;density_broken.i.dataset.archive->matrices[1].values={-1,0,0,0,-.5,0,0,0,-.5};rotate_h(density_broken);
    auto& dp=density_broken.i.dataset.archive->matrices[2].values;std::fill(dp.begin(),dp.end(),0);
    for(const auto& mo:density_broken.w.orbitals)for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k)dp[j*3+k]+=mo.occupation*mo.coefficients[j]*mo.coefficients[k];
    const auto db=cov::build_nbo_salc_model(density_broken.w,density_broken.i);require(db.available&&db.group_verified&&db.energies[0].available,"broken density destroyed valid geometry/operator");require(db.energies[0].fock_symmetry_error<1e-8&&db.energies[0].density_symmetry_checked&&db.energies[0].density_symmetry_error>.1&&!db.energies[0].electronic_symmetry_verified,"broken density incorrectly inherits geometric symmetry");
    auto fock_broken=fixture();rotate_h(fock_broken);fock_broken.w.orbitals[1].occupation=fock_broken.w.orbitals[2].occupation=1;fock_broken.i.dataset.archive->matrices[2].values={2,0,0,0,1,0,0,0,1};auto& fv=fock_broken.i.dataset.archive->matrices[1].values;std::fill(fv.begin(),fv.end(),0);for(const auto& mo:fock_broken.w.orbitals)for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k)fv[j*3+k]+=mo.energy_hartree*mo.coefficients[j]*mo.coefficients[k];refresh(fock_broken);
    const auto fb=cov::build_nbo_salc_model(fock_broken.w,fock_broken.i);require(fb.available&&fb.group_verified&&fb.energies[0].available,"broken Fock destroyed valid same-operator identity");require(fb.energies[0].fock_symmetry_error>.1&&fb.energies[0].density_symmetry_error<1e-8&&!fb.energies[0].electronic_symmetry_verified,"broken Fock incorrectly inherits geometric symmetry");
    auto truncated=fixture();truncated.w.orbitals.pop_back();refresh(truncated);auto tm=cov::build_nbo_salc_model(truncated.w,truncated.i);require(tm.available&&!tm.energies[0].available,"unverified Fock action beyond canonical span supplied side energy");require(tm.energies[0].overlap_numerical_rank==3&&tm.energies[0].canonical_effective_rank==2&&tm.energies[0].canonical_null_directions==1,"overlap numerical rank confused with retained canonical rank");require(tm.energies[0].outside_canonical_nao_norm>.5,"missing canonical direction hidden by normalization");
    auto stale=fixture();stale.w.orbitals[0].coefficients[0]+=.1;require(!cov::build_nbo_salc_model(stale.w,stale.i).available,"stale canonical identity accepted");
    auto unclosed=fixture();unclosed.i.dataset.naos[2].type="Ryd(2s)";auto raw=cov::build_nbo_salc_model(unclosed.w,unclosed.i);require(raw.available&&raw.orbitals.size()==3,"nonclosed family erased directions");for(const auto& o:raw.orbitals)require(!o.symmetry_adapted,"nonclosed literal family asserted SALC");
    auto linear=fixture();linear.w.atoms[1].z=linear.w.atoms[2].z=0;refresh(linear);auto lm=cov::build_nbo_salc_model(linear.w,linear.i);require(lm.available&&lm.group_verified&&lm.operations.size()>2&&lm.used_group.find("finite sampling")!=std::string::npos,"linear symmetry silently limited to E/i");
    auto open=fixture();open.w.orbital_occupation_model=cov::OrbitalOccupationModel::ExplicitSpin;open.i.dataset.archive->open_shell=true;open.i.dataset.association.canonical_evidence.clear();open.i.orbitals.clear();open.i.dataset.naos.clear();const auto original=open.w.orbitals;open.w.orbitals.clear();const auto descriptors=f.i.orbitals;const auto naos=f.i.dataset.naos;auto fock=open.i.dataset.archive->matrices[1];open.i.dataset.archive->matrices.resize(1);
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){cov::NboCanonicalEvidence ev;ev.spin=spin;ev.direct_fchk_coefficients=true;ev.density_verified=true;open.i.dataset.association.canonical_evidence.push_back(ev);fock.spin=spin;open.i.dataset.archive->matrices.push_back(fock);open.i.dataset.archive->matrices.push_back(matrix("DENSITY",{1,0,0,0,.5,.5,0,.5,.5},spin));for(auto mo:original){mo.spin=spin==cov::NboSpin::Beta?cov::Spin::Beta:cov::Spin::Alpha;mo.occupation/=2;open.w.orbitals.push_back(mo);}for(auto d:descriptors){d.ref.spin=spin;if(d.occupation)*d.occupation/=2;open.i.orbitals.push_back(d);}for(auto row:naos){row.spin=spin;row.occupation/=2;open.i.dataset.naos.push_back(row);}}
    refresh(open);auto om=cov::build_nbo_salc_model(open.w,open.i);require(om.orbitals.size()==6&&om.energies.size()==2&&om.energies[0].available&&om.energies[1].available,"open-shell blocks lost");for(const auto& o:om.orbitals)for(const auto& term:o.terms)require(term.orbital.spin==o.spin,"SALC render spin lost");
    auto withheld=fixture();withheld.i.dataset.association.canonical_evidence.front().density_verified=false;
    const auto wm=cov::build_nbo_salc_model(withheld.w,withheld.i);
    for(const auto& o:wm.orbitals)if(o.terms.size()>1)require(!o.occupation,"unverified full density fabricated combination occupation");
    require(!wm.energies.front().density_symmetry_checked,"unverified density supplied electronic symmetry evidence");
    // RO beta has verified density but no independent canonical beta columns.
    for(auto& ev:open.i.dataset.association.canonical_evidence)if(ev.spin==cov::NboSpin::Beta)ev.direct_fchk_coefficients=false;
    std::erase_if(open.w.orbitals,[](const auto& mo){return mo.spin==cov::Spin::Beta;});refresh(open);
    const auto ro=cov::build_nbo_salc_model(open.w,open.i);
    std::size_t beta_count=0;for(std::size_t j=0;j<ro.orbitals.size();++j)if(ro.orbitals[j].spin==cov::NboSpin::Beta){
        ++beta_count;require(ro.orbitals[j].occupation.has_value(),"RO beta density occupation incorrectly requires canonical beta");
        require(!ro.orbitals[j].energy_hartree,"RO beta fabricated canonical-validated energy");
        for(const auto& link:ro.links)require(link.side_index!=j,"RO beta fabricated canonical link");
    }require(beta_count==3,"RO beta fixed space lost");
    // Physical spin operators need not diagonalize a shared RO canonical set.
    // Independently printed diagonals in two full local bases verify their
    // side expectations without manufacturing canonical beta columns.
    auto physical=open;
    physical.i.dataset.orbitals.clear();physical.i.dataset.matrices.clear();
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}) {
        auto fit=std::find_if(physical.i.dataset.archive->matrices.begin(),
            physical.i.dataset.archive->matrices.end(),[&](const auto& x){return x.kind=="FOCK"&&x.spin==spin;});
        fit->values[0]+=spin==cov::NboSpin::Alpha?.13:-.09;
        const double h=std::sqrt(.5);
        const auto transform=matrix("AONBO",{1,0,0,0,h,h,0,h,-h},spin);
        physical.i.dataset.matrices.push_back(transform);
        for(auto& row:physical.i.dataset.naos)if(row.spin==spin)
            row.energy_hartree=fit->values[(row.id-1)*3+row.id-1];
        for(std::size_t j=0;j<3;++j) {
            cov::NboOrbital row;row.id=j+1;row.spin=spin;row.diagonal_fock_hartree=0;
            for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b)
                *row.diagonal_fock_hartree+=transform.values[a*3+j]*fit->values[a*3+b]*transform.values[b*3+j];
            physical.i.dataset.orbitals.push_back(row);
        }
    }
    const auto pm=cov::build_nbo_salc_model(physical.w,physical.i);
    require(pm.energies.size()==2,"physical spin blocks lost");
    for(const auto& ev:pm.energies)require(ev.available&&ev.printed_operator_verified&&
        !ev.canonical_same_operator&&ev.electronic_symmetry_verified&&ev.fock_symmetry_error<=2e-5&&ev.density_symmetry_checked&&ev.density_symmetry_error<=2e-4&&
        ev.printed_nao_checked==3&&ev.printed_nbo_checked==3,
        "physical Fock symmetry must be checked independently without asserting canonical same-operator identity");
    for(std::size_t j=0;j<pm.orbitals.size();++j)if(pm.orbitals[j].spin==cov::NboSpin::Beta){
        require(pm.orbitals[j].energy_hartree.has_value(),"independently verified RO beta energy lost");
        for(const auto& link:pm.links)require(link.side_index!=j,"physical beta energy invented canonical beta links");
    }
    require(cov::nbo_canonical_fingerprint(physical.w)==physical.i.canonical_fingerprint,
        "physical operator route altered canonical data");
    auto physical_split=physical;
    for(auto& matrix:physical_split.i.dataset.archive->matrices)if(matrix.kind=="FOCK"){matrix.values[4]+=.05;matrix.values[8]-=.05;}
    for(auto& row:physical_split.i.dataset.naos){if(row.id==2)*row.energy_hartree+=.05;if(row.id==3)*row.energy_hartree-=.05;}
    const auto split_physical=cov::build_nbo_salc_model(physical_split.w,physical_split.i);
    for(const auto& ev:split_physical.energies)require(ev.available&&ev.printed_operator_verified&&!ev.canonical_same_operator&&ev.fock_symmetry_error>.09&&!ev.electronic_symmetry_verified,"Actual independently verified spin Fock breaking geometry must remain noninvariant");
    auto physical_density=physical;
    for(auto& matrix:physical_density.i.dataset.archive->matrices)if(matrix.kind=="DENSITY"){matrix.values[4]+=.1;matrix.values[8]-=.1;}
    for(auto& row:physical_density.i.dataset.naos){if(row.id==2)row.occupation+=.1;if(row.id==3)row.occupation-=.1;}
    const auto density_physical=cov::build_nbo_salc_model(physical_density.w,physical_density.i);
    for(const auto& ev:density_physical.energies)require(ev.available&&ev.printed_operator_verified&&ev.fock_symmetry_error<2e-5&&ev.density_symmetry_error>.19&&!ev.electronic_symmetry_verified,"Invariant physical Fock cannot conceal broken actual density symmetry");
    auto absent_print=physical;absent_print.i.dataset.naos.front().energy_hartree.reset();
    require(!cov::build_nbo_salc_model(absent_print.w,absent_print.i).energies[0].available,
        "partial printed evidence incorrectly enabled physical spin energy");
    auto wrong_local=physical;wrong_local.i.dataset.orbitals.front().diagonal_fock_hartree=5;
    require(!cov::build_nbo_salc_model(wrong_local.w,wrong_local.i).energies[0].available,
        "NBO second-basis mismatch accepted");
    auto wrong_spin=physical;
    std::swap(wrong_spin.i.dataset.archive->matrices[1].values,
              wrong_spin.i.dataset.archive->matrices[3].values);
    require(!cov::build_nbo_salc_model(wrong_spin.w,wrong_spin.i).energies[0].available,
        "swapped physical spin Fock accepted");
    auto wrong_offdiag=physical;wrong_offdiag.i.dataset.archive->matrices[1].values[5]+=.03;
    wrong_offdiag.i.dataset.archive->matrices[1].values[7]+=.03;
    require(!cov::build_nbo_salc_model(wrong_offdiag.w,wrong_offdiag.i).energies[0].available,
        "Hermitian off-diagonal perturbation escaped independent NBO check");
    auto bad_coupling=physical;cov::NboE2 e2;e2.spin=cov::NboSpin::Alpha;
    e2.donor=1;e2.acceptor=2;e2.fock_hartree=2;e2.energy_gap_hartree=.37;
    bad_coupling.i.dataset.e2.push_back(e2);
    require(!cov::build_nbo_salc_model(bad_coupling.w,bad_coupling.i).energies[0].available,
        "printed coupling disagreement accepted");
    auto wrong_units=physical;for(auto& value:wrong_units.i.dataset.archive->matrices[1].values)value*=27.211386245981;
    require(!cov::build_nbo_salc_model(wrong_units.w,wrong_units.i).energies[0].available,
        "incorrect physical Fock units accepted");
    const auto path=std::filesystem::temp_directory_path()/"cov_salc_ev_smoke.47";{std::ofstream o(path);o<<"$GENNBO NATOMS=1 NBAS=1 UPPER OPEN EV $END\n$COORD\nfixture\n1 1 0 0 0\n$END\n$BASIS CENTER=1 LABEL=1 $END\n$CONTRACT NCOMP=1 NPRIM=1 NPTR=1 EXP=1 CS=1 $END\n$FOCK -27.211386245981 -54.422772491962 $END\n";}const auto ar=cov::read_nbo_archive(path);std::filesystem::remove(path);require(ar.fock_input_units=="eV"&&ar.matrices.size()==2,"EV/OPEN archive metadata");near(ar.matrices[0].values[0],-1,"alpha eV conversion");near(ar.matrices[1].values[0],-2,"beta eV conversion");
    std::cout<<"NBO SALC core: fixed SALCs, metric/group closure, signed reconstruction, rotation/permutation/phase invariance, spin, energies and failure gates passed\n";return 0;
}catch(const std::exception& e){std::cerr<<"NBO SALC core failed: "<<e.what()<<'\n';return 1;}}
