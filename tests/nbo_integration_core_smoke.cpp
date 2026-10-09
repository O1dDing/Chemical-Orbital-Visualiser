#include "cov/nbo_integration.hpp"
#include "cov/wavefunction_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool x,const char* message){if(!x)throw std::runtime_error(message);}
void near(double x,double y,const char* message){require(std::abs(x-y)<1e-9,message);}
cov::NboMatrix mat(const char* kind,std::vector<double> values){cov::NboMatrix m;m.kind=kind;m.rows=m.columns=2;m.values=std::move(values);return m;}
struct Fixture{cov::Wavefunction w;cov::NboDataset d;};
Fixture fixture(){Fixture f;auto& w=f.w;auto& d=f.d;const double a=std::sqrt(.8),b=std::sqrt(.2);const std::vector<double> identity={1,0,0,1},rotation={a,b,b,-a};
    w.source=cov::WavefunctionSource::Fchk;w.basis_count=2;w.atoms={{"H",1,0,0,0,1},{"H",1,0,0,2,1}};w.shells={{0,0,1,0,0,0},{1,1,1,1,0,0}};w.primitives={{1,1},{1,1}};w.gaussian_ao_transform={{0,1,1},{1,1,1}};w.ao_overlap=identity;w.alpha_electrons=w.beta_electrons=1;w.electron_counts_provenance=cov::DataProvenance::Producer;w.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;
    for(std::size_t i=0;i<2;++i){cov::MolecularOrbital mo;mo.source_orbital_index=i;mo.occupation=i?0:2;mo.occupation_provenance=cov::DataProvenance::Derived;mo.coefficients=mo.gaussian_source_coefficients=i?std::vector<double>{b,-a}:std::vector<double>{a,b};w.orbitals.push_back(mo);}
    cov::NboArchive ar;ar.atoms=w.atoms;ar.basis_count=2;ar.density_is_bond_order=true;ar.centers={1,2};ar.labels={1,1};ar.ncomp={1,1};ar.nprim={1,1};ar.nptr={1,2};ar.exponents=ar.cs={1,1};ar.matrices={mat("OVERLAP",identity),mat("LCAOMO",rotation),mat("DENSITY",{1.6,.8,.8,.4})};d.archive=ar;
    for(const auto* key:{"AONAO","AONBO","NAONBO","AONHO","NAONHO","NHONBO","AOPNAO","NLMOMO"})d.matrices.push_back(mat(key,identity));for(const auto* key:{"NAOMO","NBOMO","AONLMO","NBONLMO"})d.matrices.push_back(mat(key,rotation));
    for(std::size_t i=0;i<2;++i){cov::NboNao n;n.id=n.atom=i+1;n.symbol="H";n.angular="s";n.type="Val(1s)";n.occupation=i?.4:1.6;d.naos.push_back(n);cov::NboOrbital o;o.id=i+1;o.atoms={i+1};o.label="LP H "+std::to_string(i+1);o.kind="LP";o.occupation=n.occupation;d.orbitals.push_back(o);}cov::NboNlmo nl;nl.id=1;nl.occupation=2;nl.parent_percent=80;nl.parent_label=d.orbitals[0].label;nl.parent_nbo=1;d.nlmos.push_back(nl);return f;
}
bool available(const cov::NboIntegration& d,const char* key){const auto* c=cov::nbo_capability(d,key);return c&&c->available();}
void erase_matrix(cov::NboDataset& d,const std::string& key){std::erase_if(d.matrices,[&](const auto& m){return m.kind==key;});}
void reconstruction(const cov::NboIntegration& d,const cov::Wavefunction& w,cov::NboOrbitalKind kind){for(std::size_t i=0;i<w.orbitals.size();++i){const auto links=cov::nbo_links_for_mo(d,i,kind);if(links.empty())continue;cov::NboOrbitalSelection s;s.dataset_id=d.id;s.mode=cov::NboSelectionMode::PartialSum;s.target_canonical_index=i;for(const auto& x:links)s.terms.push_back({x.orbital,x.coefficient});auto view=cov::make_nbo_selection_view(d,w,s);require(view.available&&view.reconstruction_error&&*view.reconstruction_error<2e-5,"full signed basis sum does not reconstruct MO");}}
void rejected_ordered_geometry_preserves_canonical(Fixture f) {
    const auto original=f.w;
    const auto original_id=cov::integrate_nbo(f.w,{}).canonical_fingerprint;
    const auto rejected=cov::integrate_nbo(f.w,f.d);
    const auto* association=cov::nbo_capability(rejected,"source_association");
    require(association&&association->state==cov::NboCapabilityState::Rejected,
            "Reordered/rotated source was silently registered");
    require(rejected.dataset.association.status=="incompatible_geometry",
            "Ordered geometry mismatch lost its exact rejection reason");
    require(available(rejected,"canonical")&&available(rejected,"gaussian_ao")&&
            !available(rejected,"nao")&&!available(rejected,"nbo"),
            "Rejected geometry corrupted canonical fallback or exposed NBO fields");
    require(rejected.canonical_fingerprint==original_id&&
            cov::integrate_nbo(f.w,{}).canonical_fingerprint==original_id,
            "Rejected source changed canonical identity");
    require(f.w.ao_overlap==original.ao_overlap&&f.w.orbitals.size()==original.orbitals.size(),
            "Rejected source changed canonical dimensions/metric");
    for(std::size_t i=0;i<original.orbitals.size();++i) {
        const auto& a=f.w.orbitals[i];const auto& b=original.orbitals[i];
        require(a.coefficients==b.coefficients&&a.gaussian_source_coefficients==b.gaussian_source_coefficients&&
                a.occupation==b.occupation&&a.source_orbital_index==b.source_orbital_index,
                "Rejected source changed canonical coefficients/occupation/order");
        const auto* retained=cov::nbo_orbital(rejected,{cov::NboOrbitalKind::Canonical,cov::NboSpin::Total,i});
        require(retained&&retained->coefficients==b.coefficients,
                "Fallback descriptor did not retain the canonical field");
    }
}
void ordered_geometry_rejection_controls() {
    // Relabel the two source atoms, retaining AO column identities by relabelling
    // their centres. This is the same molecule with a different ordered source.
    auto permuted=fixture();std::swap(permuted.d.archive->atoms[0],permuted.d.archive->atoms[1]);
    for(auto& centre:permuted.d.archive->centers)centre=3-centre;
    for(auto& row:permuted.d.naos)row.atom=3-row.atom;
    for(auto& orbital:permuted.d.orbitals)for(auto& atom:orbital.atoms)atom=3-atom;
    rejected_ordered_geometry_preserves_canonical(std::move(permuted));
    // A 45-degree rigid rotation about y. The fixture has only s AOs, so its
    // coefficient columns need no angular rotation; only the source frame differs.
    auto rotated=fixture();const double c=std::sqrt(.5);
    for(auto& atom:rotated.d.archive->atoms) {
        const double x=atom.x,z=atom.z;atom.x=c*(x+z);atom.z=c*(z-x);
    }
    rejected_ordered_geometry_preserves_canonical(std::move(rotated));
}
}
int main(int argc,char** argv){try{
    auto identity_fixture=fixture();
    const auto identity_before=cov::integrate_nbo(identity_fixture.w,identity_fixture.d);
    identity_fixture.d.naos.front().type="Ryd(2s)";
    require(cov::integrate_nbo(identity_fixture.w,identity_fixture.d).id!=identity_before.id,
            "NAO report classification change retained stale integration identity");
    cov::NboOrbitalSelection stale_report_selection;stale_report_selection.dataset_id=identity_before.id;
    stale_report_selection.terms.push_back({{cov::NboOrbitalKind::NAO,cov::NboSpin::Total,0},1});
    require(!cov::make_nbo_selection_view(cov::integrate_nbo(identity_fixture.w,identity_fixture.d),identity_fixture.w,stale_report_selection).available,
            "selection accepted an identity from the previous report contents");
    const auto changed_id=[&](auto mutate,const char* reason){auto copy=fixture();mutate(copy.d);
        require(cov::integrate_nbo(copy.w,copy.d).id!=identity_before.id,reason);};
    changed_id([](auto& d){d.naos[0].energy_hartree=-.3;},"NAO energy omitted from identity");
    changed_id([](auto& d){d.nlmos[0].parent_nbo=2;},"NLMO parent omitted from identity");
    changed_id([](auto& d){d.orbitals[0].components.push_back({1,50,.7,"s"});},"NHO component evidence omitted from identity");
    changed_id([](auto& d){cov::NboE2Section e;e.printing_threshold=1;d.e2_sections.push_back(e);},"E2 missing-value semantics omitted from identity");
    changed_id([](auto& d){d.archive->density_is_bond_order=false;},"density convention omitted from identity");
    changed_id([](auto& d){d.archive->atoms[0].x+=.1;},"archive association identity omitted from identity");
    changed_id([](auto& d){d.matrices[0].rows=1;},"matrix dimensions omitted from identity");
    auto density_identity=fixture();density_identity.w.total_density_packed={1.6,.8,.4};
    density_identity.w.total_density_provenance=cov::DataProvenance::Producer;
    require(cov::nbo_canonical_fingerprint(density_identity.w)!=identity_before.canonical_fingerprint,
            "producer density omitted from immutable canonical identity");
    identity_fixture=fixture();identity_fixture.d.nlmos.front().parent_percent=70;
    require(cov::integrate_nbo(identity_fixture.w,identity_fixture.d).id!=identity_before.id,
            "NLMO parent composition change retained stale integration identity");
    ordered_geometry_rejection_controls();
    auto f=fixture();const auto before=f.w.orbitals[0].coefficients;auto d=cov::integrate_nbo(f.w,f.d);
    for(const auto* k:{"canonical","gaussian_ao","source_association","nao","aomo","nbo","nho","nlmo","pnao"})require(available(d,k),k);
    for(auto kind:{cov::NboOrbitalKind::GaussianAO,cov::NboOrbitalKind::NAO,cov::NboOrbitalKind::NBO,cov::NboOrbitalKind::NHO,cov::NboOrbitalKind::NLMO})reconstruction(d,f.w,kind);
    auto links=cov::nbo_links_for_mo(d,0);require(links.size()==2,"full mapping lost rows");auto v=cov::make_nbo_selection_view(d,f.w,cov::nbo_component_selection(d,links[0]));require(v.available,"weighted component unavailable");near(v.metric_norm2[0],.8,"component silently normalized");
    auto left=cov::make_nbo_selection_view(d,f.w,cov::nbo_fragment_selection(d,0,{0}));auto right=cov::make_nbo_selection_view(d,f.w,cov::nbo_fragment_selection(d,0,{1}));require(left.available&&right.available,"fragment partial sums unavailable");for(std::size_t i=0;i<2;++i)near(left.wavefunction.orbitals[0].coefficients[i]+right.wavefunction.orbitals[0].coefficients[i],before[i],"fragment complement lost interference amplitude");
    const cov::NboOrbitalRef nl{cov::NboOrbitalKind::NLMO,cov::NboSpin::Total,0};auto terms=cov::nbo_nlmo_components(d,nl);require(terms.size()==2&&cov::nbo_nlmo_parent(d,nl).has_value(),"NLMO signed decomposition/printed parent missing");cov::NboOrbitalSelection s;s.dataset_id=d.id;s.mode=cov::NboSelectionMode::PartialSum;s.terms=terms;v=cov::make_nbo_selection_view(d,f.w,s);require(v.available,"NLMO sum unavailable");for(std::size_t i=0;i<2;++i)near(v.wavefunction.orbitals[0].coefficients[i],cov::nbo_orbital(d,nl)->coefficients[i],"NLMO main plus tail not full orbital");require(cov::nbo_nho_components(d,{cov::NboOrbitalKind::NHO,cov::NboSpin::Total,0}).size()==2,"NHO NAO hybrid composition missing");
    s.mode=cov::NboSelectionMode::Overlay;v=cov::make_nbo_selection_view(d,f.w,s);require(v.available&&v.wavefunction.orbitals.size()==2,"overlay collapsed orbitals");s.dataset_id="wrong";require(!cov::make_nbo_selection_view(d,f.w,s).available,"stale dataset accepted");s.dataset_id=d.id;auto changed=f.w;changed.orbitals[0].coefficients[0]+=.1;require(!cov::make_nbo_selection_view(d,changed,s).available,"changed canonical accepted");
    s.mode=cov::NboSelectionMode::Combination;s.terms={{terms[0].orbital,1},{terms[0].orbital,-1}};s.normalize=true;require(!cov::make_nbo_selection_view(d,f.w,s).available,"zero normalized sum accepted");
    auto separate=d;auto alpha=*cov::nbo_orbital(d,terms[0].orbital),beta=alpha;alpha.ref.spin=cov::NboSpin::Alpha;beta.ref.spin=cov::NboSpin::Beta;separate.orbitals.push_back(alpha);separate.orbitals.push_back(beta);s.normalize=false;s.terms={{alpha.ref,1},{beta.ref,1}};require(!cov::make_nbo_selection_view(separate,f.w,s).available,"coherent alpha beta sum accepted");s.mode=cov::NboSelectionMode::Overlay;require(cov::make_nbo_selection_view(separate,f.w,s).available,"separate alpha beta overlay rejected");
    auto bad=f;erase_matrix(bad.d,"NAOMO");auto degraded=cov::integrate_nbo(bad.w,bad.d);require(!available(degraded,"aomo")&&available(degraded,"nao")&&available(degraded,"nbo"),"missing NAOMO destroyed independent capabilities");bad=f;for(auto& m:bad.d.matrices)if(m.kind=="NAOMO")m.values[0]+=.1;degraded=cov::integrate_nbo(bad.w,bad.d);require(!available(degraded,"aomo")&&available(degraded,"nao")&&available(degraded,"nlmo"),"bad NAOMO isolation failed");bad=f;for(auto& m:bad.d.matrices)if(m.kind=="AONLMO")m.values[0]+=.1;degraded=cov::integrate_nbo(bad.w,bad.d);require(!available(degraded,"nlmo")&&available(degraded,"aomo")&&available(degraded,"nbo"),"bad NLMO isolation failed");
    auto fallback=cov::integrate_nbo(f.w,{});require(available(fallback,"canonical")&&available(fallback,"gaussian_ao")&&!available(fallback,"nbo"),"canonical fallback unavailable");require(f.w.orbitals[0].coefficients==before,"canonical mutated");
    auto reduced=f;reduced.d.naos.resize(1);reduced.d.orbitals.resize(1);reduced.d.nlmos[0].occupation=1.6;reduced.d.nlmos[0].parent_percent=100;reduced.d.matrices.clear();
    for(const auto* kind:{"AONAO","AONBO","AONHO","AONLMO","AOPNAO"}){auto m=mat(kind,{1,0});m.columns=1;reduced.d.matrices.push_back(m);}
    for(const auto* kind:{"NAONBO","NAONHO","NAONLMO","NHONBO","NBONLMO"}){auto m=mat(kind,{1});m.rows=m.columns=1;reduced.d.matrices.push_back(m);}
    for(const auto* kind:{"NAOMO","NBOMO","NLMOMO"}){auto m=mat(kind,{std::sqrt(.8),std::sqrt(.2)});m.rows=1;reduced.d.matrices.push_back(m);}
    const auto partial=cov::integrate_nbo(reduced.w,reduced.d);require(available(partial,"aomo")&&available(partial,"nho")&&available(partial,"nlmo")&&!available(partial,"aomo_full"),"rectangular local subspace capabilities incorrect");require(partial.dataset.mo_decompositions.size()==2,"reduced subspace lost canonical identities");const auto& pd=partial.dataset.mo_decompositions[0];require(pd.available&&pd.rows.size()==1&&pd.projection_residual_norm&&pd.status=="verified_partial_subspace_projection","partial projection falsely reported complete");near(*pd.weight_sum,.8,"retained subspace weight renormalized");near(*pd.projection_residual_norm,std::sqrt(.2),"missing canonical direction residual not explicit");
    if(argc==3){auto w=cov::parse_wavefunction(argv[1]);auto discovery=cov::discover_nbo_inputs({argv[2]});require(discovery.candidates.size()==1,"single analysis directory discovery is ambiguous");auto real=cov::read_nbo_integration(w,discovery.candidates.front());for(const auto* k:{"source_association","aomo","nho","nlmo"}){if(!available(real,k)){for(const auto& x:real.diagnostics)std::cerr<<x<<'\n';throw std::runtime_error(k);}}for(auto kind:{cov::NboOrbitalKind::GaussianAO,cov::NboOrbitalKind::NAO,cov::NboOrbitalKind::NHO,cov::NboOrbitalKind::NLMO})reconstruction(real,w,kind);}
    const auto rendered=cov::make_nbo_wavefunction(partial.dataset,reduced.w);
    require(rendered.orbitals.size()==1&&rendered.basis_count==2,"rectangular NBO renderer padded the local space");
    require(rendered.orbitals[0].coefficients==std::vector<double>({1,0}),"rectangular NBO renderer changed coefficients");
    near(rendered.orbitals[0].occupation,1.6,"rectangular NBO renderer changed occupation");
    require(std::isnan(rendered.orbitals[0].energy_hartree),"NBO renderer invented canonical energy");
    auto rejects_render=[&](cov::NboDataset data,cov::Wavefunction wave,const char* message){
        bool rejected=false;try{cov::make_nbo_wavefunction(data,wave);}catch(const std::exception&){rejected=true;}
        require(rejected,message);
    };
    auto damaged=partial.dataset;damaged.orbitals[0].occupation+=.1;
    rejects_render(damaged,reduced.w,"stale association accepted changed printed occupation");
    damaged=partial.dataset;damaged.orbitals.push_back(damaged.orbitals[0]);
    rejects_render(damaged,reduced.w,"duplicate NBO identity accepted");
    damaged=partial.dataset;damaged.orbitals.clear();
    rejects_render(damaged,reduced.w,"missing NBO identity accepted");
    damaged=partial.dataset;damaged.orbitals[0].id=2;
    rejects_render(damaged,reduced.w,"out of subspace NBO identity accepted");
    damaged=partial.dataset;damaged.orbitals[0].occupation=std::numeric_limits<double>::quiet_NaN();
    rejects_render(damaged,reduced.w,"nonfinite NBO occupation accepted");
    damaged=partial.dataset;damaged.orbitals[0].spin=cov::NboSpin::Beta;
    rejects_render(damaged,reduced.w,"wrong NBO report spin accepted");
    auto aonbo=[](cov::NboDataset& data)->cov::NboMatrix&{return *std::find_if(data.matrices.begin(),data.matrices.end(),[](const auto& m){return m.kind=="AONBO";});};
    damaged=partial.dataset;damaged.matrices.push_back(aonbo(damaged));
    rejects_render(damaged,reduced.w,"duplicate AONBO identity accepted");
    damaged=partial.dataset;aonbo(damaged).values[0]=.8;
    rejects_render(damaged,reduced.w,"nonorthonormal AONBO accepted");
    damaged=partial.dataset;aonbo(damaged).values[0]=std::numeric_limits<double>::quiet_NaN();
    rejects_render(damaged,reduced.w,"nonfinite AONBO accepted");
    damaged=partial.dataset;aonbo(damaged).values.pop_back();
    rejects_render(damaged,reduced.w,"truncated AONBO accepted");
    damaged=partial.dataset;aonbo(damaged).columns=0;
    rejects_render(damaged,reduced.w,"empty AONBO accepted");
    damaged=partial.dataset;aonbo(damaged).spin=cov::NboSpin::Beta;
    rejects_render(damaged,reduced.w,"wrong AONBO spin accepted");
    damaged=partial.dataset;for(auto& m:damaged.archive->matrices)if(m.kind=="DENSITY")m.values[0]+=.1;
    rejects_render(damaged,reduced.w,"stale association accepted changed source density");
    auto changed_geometry=reduced.w;changed_geometry.atoms[1].z+=.1;
    rejects_render(partial.dataset,changed_geometry,"stale association accepted changed canonical geometry");
    damaged=partial.dataset;for(auto& m:damaged.matrices)if(m.kind=="NAOMO")m.values[0]+=.1;
    require(cov::make_nbo_wavefunction(damaged,reduced.w).orbitals.size()==1,"unrelated transform broke NBO rendering");
    // The two archive density conventions must agree in a nonorthogonal AO metric.
    const double metric_offdiag=.25,k=1/std::sqrt(1-metric_offdiag*metric_offdiag);
    const std::vector<double> a={1,-metric_offdiag*k,0,k},at={1,0,-metric_offdiag*k,k};
    const std::vector<double> metric={1,metric_offdiag,metric_offdiag,1};
    auto product=[](const auto& x,const auto& y){std::vector<double> z(4);for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int l=0;l<2;++l)z[i*2+j]+=x[i*2+l]*y[l*2+j];return z;};
    auto metric_wave=reduced.w;metric_wave.ao_overlap=metric;
    for(auto& mo:metric_wave.orbitals){const auto c=mo.coefficients;mo.coefficients=mo.gaussian_source_coefficients={a[0]*c[0]+a[1]*c[1],a[2]*c[0]+a[3]*c[1]};}
    auto metric_data=partial.dataset;
    for(auto& m:metric_data.archive->matrices){if(m.kind=="OVERLAP")m.values=metric;else if(m.kind=="LCAOMO")m.values=product(a,m.values);else if(m.kind=="DENSITY")m.values=product(product(a,m.values),at);}
    near(cov::make_nbo_wavefunction(metric_data,metric_wave).orbitals[0].occupation,1.6,"BODM density convention rejected");
    metric_data.archive->density_is_bond_order=false;
    for(auto& m:metric_data.archive->matrices)if(m.kind=="DENSITY")m.values=product(product(metric,m.values),metric);
    near(cov::make_nbo_wavefunction(metric_data,metric_wave).orbitals[0].occupation,1.6,"covariant density convention rejected");
    require(reduced.w.orbitals.size()==2&&reduced.w.orbitals[0].coefficients==before,"NBO renderer changed canonical space");
    std::cout<<"NBO integration core passed: signed reconstruction, rectangular rendering, metric/density rejection controls and immutable canonical\n";return 0;
}catch(const std::exception& e){std::cerr<<"NBO integration core failed: "<<e.what()<<'\n';return 1;}}
