#include "cov/nbo.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
void near(double a,double b,const char* why){require(std::abs(a-b)<1e-10,why);}
cov::NboMatrix mat(const char* kind,const std::vector<double>& values,cov::NboSpin spin=cov::NboSpin::Total){cov::NboMatrix x;x.kind=kind;x.rows=x.columns=2;x.values=values;x.spin=spin;return x;}
struct Fixture {cov::Wavefunction w;cov::NboDataset d;};
Fixture fixture(){
    Fixture f;auto& w=f.w;auto& d=f.d;const double a=std::sqrt(.8),b=std::sqrt(.2);
    w.source=cov::WavefunctionSource::Fchk;w.basis_count=2;w.atoms={{"H",1,0,0,0,1},{"H",1,0,0,2,1}};
    w.shells={{0,0,1,0,0,0},{1,1,1,1,0,0}};w.primitives={{1,1},{1,1}};w.gaussian_ao_transform={{0,1,1},{1,1,1}};w.ao_overlap={1,0,0,1};
    w.alpha_electrons=w.beta_electrons=1;w.electron_counts_provenance=cov::DataProvenance::Producer;w.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;
    for(std::size_t i=0;i<2;++i){cov::MolecularOrbital mo;mo.source_orbital_index=i;mo.occupation=i?0:2;mo.occupation_provenance=cov::DataProvenance::Derived;mo.coefficients=mo.gaussian_source_coefficients=i?std::vector<double>{b,-a}:std::vector<double>{a,b};w.orbitals.push_back(mo);}
    cov::NboArchive ar;ar.atoms=w.atoms;ar.basis_count=2;ar.density_is_bond_order=true;ar.centers={1,2};ar.labels={1,1};ar.ncomp={1,1};ar.nprim={1,1};ar.nptr={1,2};ar.exponents={1,1};ar.cs={1,1};
    ar.matrices={mat("OVERLAP",{1,0,0,1}),mat("LCAOMO",{a,b,b,-a}),mat("DENSITY",{1.6,.8,.8,.4})};d.archive=ar;
    d.matrices={mat("AONBO",{1,0,0,1}),mat("NBOMO",{a,b,b,-a}),mat("AONAO",{1,0,0,1}),mat("NAOMO",{a,b,b,-a}),mat("NAONBO",{1,0,0,1})};
    for(std::size_t i=0;i<2;++i){cov::NboOrbital o;o.id=i+1;o.occupation=i?.4:1.6;d.orbitals.push_back(o);cov::NboNao n;n.id=n.atom=i+1;n.symbol="H";n.angular="s";n.type="Val( 1s)";n.occupation=o.occupation;d.naos.push_back(n);}return f;
}
}
int main(){try{
    auto f=fixture();auto original=f.w.orbitals;require(cov::associate_nbo(f.d,f.w).compatible,"complete NAO association");
    const auto* v=cov::nbo_mo_decomposition(f.d,0);require(v&&v->available&&v->rows.size()==2&&v->atoms.size()==2&&v->shells.size()==2,"complete decomposition identity");near(*v->weight_sum,1,"normalization");near(v->rows[0].weight,.8,"NAOMO direction");near(*v->rows[0].electron_contribution,1.6,"occupation-weighted contribution");require(v->rows[0].principal_n==1&&v->rows[0].angular_l==0,"literal shell classification");
    require(f.w.orbitals[0].gaussian_source_coefficients==original[0].gaussian_source_coefficients&&f.w.orbitals[0].occupation==original[0].occupation,"canonical source mutation");
    require(!cov::nbo_mo_decomposition(f.d,7),"out-of-range canonical lookup");
    auto phase=fixture();for(auto& x:phase.w.orbitals[0].gaussian_source_coefficients)x=-x;for(auto& x:phase.w.orbitals[0].coefficients)x=-x;const auto phased_source=phase.w.orbitals[0].gaussian_source_coefficients;require(cov::associate_nbo(phase.d,phase.w).compatible,"whole-orbital phase rejected");require(phase.w.orbitals[0].gaussian_source_coefficients==phased_source,"phase comparison modified source");v=cov::nbo_mo_decomposition(phase.d,0);near(v->rows[0].weight,.8,"phase changed weights");require(phase.d.association.provenance_status=="producer_step_unverified","numeric match claimed known producer step");
    auto bad=fixture();bad.d.matrices[3].values[0]+=.01;require(!cov::associate_nbo(bad.d,bad.w).compatible,"wrong transform accepted");require(!cov::nbo_mo_decomposition(bad.d,0)->available,"failed association leaves stale decomposition");
    bad=fixture();bad.d.naos[1].id=1;require(!cov::associate_nbo(bad.d,bad.w).compatible,"duplicate NAO identity accepted");
    bad=fixture();bad.d.naos[1].spin=cov::NboSpin::Beta;require(!cov::associate_nbo(bad.d,bad.w).compatible,"NAO labels merged across spins");
    bad=fixture();bad.d.naos[0].occupation+=.1;require(!cov::associate_nbo(bad.d,bad.w).compatible,"wrong NAO population accepted");
    bad=fixture();bad.d.matrices.pop_back();require(cov::associate_nbo(bad.d,bad.w).compatible,"optional NAONBO incorrectly required");
    bad=fixture();bad.d.matrices.erase(bad.d.matrices.begin()+2);require(!cov::associate_nbo(bad.d,bad.w).compatible,"NAOMO without AONAO accepted");
    bad=fixture();bad.d.matrices.resize(2);require(cov::associate_nbo(bad.d,bad.w).compatible,"old data compatibility");require(!cov::nbo_mo_decomposition(bad.d,0)->available&&!cov::nbo_mo_decomposition(bad.d,0)->weight_sum,"missing matrices became zero contribution");
    bad=fixture();bad.d.naos[0].type="Ryd( unknown)";require(cov::associate_nbo(bad.d,bad.w).compatible,"unclassified label must retain numerical data");v=cov::nbo_mo_decomposition(bad.d,0);require(!v->rows[0].principal_n&&!v->rows[0].angular_l&&v->shells[0].label.find("unclassified")!=std::string::npos,"invented shell label");
    bad=fixture();bad.d.matrices.push_back(bad.d.matrices[3]);require(!cov::associate_nbo(bad.d,bad.w).compatible,"duplicate matrix accepted");
    bad=fixture();bad.d.archive->matrices[1].values[1]=bad.d.archive->matrices[1].values[3]=0;bad.d.matrices[1].values[1]=bad.d.matrices[1].values[3]=0;bad.d.matrices[3].values[1]=bad.d.matrices[3].values[3]=0;bad.w.orbitals.pop_back();require(!cov::associate_nbo(bad.d,bad.w).compatible,"zero padding hid missing physical MO in full-rank NAO basis");require(!bad.d.nao_validation.empty()&&bad.d.nao_validation[0].column_status[1]=="null_padding","null padding provenance lost");
    // OPEN archive with a shared AONAO, directly checked alpha and density-only
    // beta. Alpha electronic contributions must use one electron, not the
    // original shared FCHK occupation of two.
    auto op=fixture();op.d.archive->open_shell=true;op.d.archive->matrices.resize(1);op.d.matrices.clear();op.d.naos.clear();op.d.orbitals.clear();const double a=std::sqrt(.8),b=std::sqrt(.2);
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){op.d.archive->matrices.push_back(mat("LCAOMO",{a,b,b,-a},spin));op.d.archive->matrices.push_back(mat("DENSITY",{.8,.4,.4,.2},spin));for(const auto& kind:{"AONBO","NAONBO"})op.d.matrices.push_back(mat(kind,{1,0,0,1},spin));for(const auto& kind:{"NBOMO","NAOMO"})op.d.matrices.push_back(mat(kind,{a,b,b,-a},spin));for(std::size_t i=0;i<2;++i){cov::NboNao n;n.id=n.atom=i+1;n.symbol="H";n.angular="s";n.type="Val(1s)";n.occupation=i?.2:.8;n.spin=spin;op.d.naos.push_back(n);cov::NboOrbital o;o.id=i+1;o.occupation=n.occupation;o.spin=spin;op.d.orbitals.push_back(o);}}
    op.d.matrices.push_back(mat("AONAO",{1,0,0,1}));op.w.total_density_packed={1.6,.8,.4};op.w.spin_density_packed={0,0,0};op.w.total_density_provenance=op.w.spin_density_provenance=cov::DataProvenance::Producer;
    require(cov::associate_nbo(op.d,op.w).compatible,"OPEN shared AONAO");v=cov::nbo_mo_decomposition(op.d,0);require(v&&v->available&&v->occupation==1,"alpha spin occupation confused with shared canonical occupation");near(*v->rows[0].electron_contribution,.8,"alpha electron decomposition");require(op.d.nao_validation.size()==2&&op.d.nao_validation[1].available&&!op.d.nao_validation[1].direct_fchk_coefficients,"beta archive-only boundary");require(op.d.mo_decompositions.size()==2,"guessed beta canonical columns created");
    auto unrelated_naomo=op.d;
    for(auto& matrix:unrelated_naomo.matrices)
        if(matrix.kind=="NAOMO"&&matrix.spin==cov::NboSpin::Alpha){matrix.values[0]+=.1;break;}
    auto independently_checked=unrelated_naomo;
    require(!cov::associate_nbo(independently_checked,op.w).compatible,
            "invalid NAOMO unexpectedly passed its independent full transform check");
    require(cov::make_nbo_wavefunction(unrelated_naomo,op.w).orbitals.size()==4,
            "rejected NAOMO incorrectly suppressed independently validated AONBO rendering");
    auto independent=op.d;
    independent.matrices.push_back(mat("AOPNAO",{1,0,0,1})); // shared RO sidecar unrelated to NBO rendering
    const auto nbo_view=cov::make_nbo_wavefunction(independent,op.w);
    require(nbo_view.orbitals.size()==4 &&
            std::count_if(nbo_view.orbitals.begin(),nbo_view.orbitals.end(),
                [](const auto& mo){return mo.spin==cov::Spin::Alpha;})==2 &&
            std::count_if(nbo_view.orbitals.begin(),nbo_view.orbitals.end(),
                [](const auto& mo){return mo.spin==cov::Spin::Beta;})==2,
            "independently rejected NAOMO/shared AOPNAO suppressed verified RO AONBO rendering");
    auto damaged=independent;
    for(auto& matrix:damaged.matrices)
        if(matrix.kind=="AONBO"&&matrix.spin==cov::NboSpin::Alpha){matrix.values[0]+=.1;break;}
    bool rejected_aonbo=false;
    try{(void)cov::make_nbo_wavefunction(damaged,op.w);}
    catch(const std::runtime_error& error){
        rejected_aonbo=std::string(error.what()).find("invalid_aonbo_metric")!=std::string::npos;
    }
    require(rejected_aonbo,"damaged actual AONBO coefficients passed render-time strict recheck");
    const auto json=cov::serialize_nbo_json(op.d);require(json.find("\"mo_decompositions\"")!=std::string::npos&&json.find("verified_archive_only")!=std::string::npos,"decomposition JSON omissions");
    std::cout<<"NBO focus core: passed full contributions, literal shells, spin/density closure, missing data and rejection contracts\n";return 0;
}catch(const std::exception& e){std::cerr<<"NBO focus core failed: "<<e.what()<<'\n';return 1;}}
