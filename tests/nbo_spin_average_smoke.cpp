#include "cov/nbo_spin_average.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using M=Eigen::MatrixXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
void near(double a,double b,const char* why){require(std::abs(a-b)<1e-10,why);}
std::vector<double> flat(const M& m){RM r=m;return {r.data(),r.data()+r.size()};}
struct Fixture {cov::Wavefunction w;cov::NboIntegration data;cov::NboSalcModel raw;M common,fa,fb,pa,pb;};
Fixture fixture(){Fixture f;auto& w=f.w;auto& data=f.data;auto& raw=f.raw;
    const M eye=M::Identity(4,4);w.basis_count=4;w.atoms={{"C",6,0,0,0,6}};w.ao_overlap=flat(eye);
    w.source=cov::WavefunctionSource::Fchk;w.source_route="SP ROPBE1PBE synthetic";
    w.alpha_electrons=3;w.beta_electrons=2;w.multiplicity=2;
    w.electron_counts_provenance=w.multiplicity_provenance=cov::DataProvenance::Producer;
    w.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;
    for(int i=0;i<4;++i){cov::MolecularOrbital mo;mo.coefficients.assign(4,0);mo.coefficients[i]=1;
        mo.gaussian_source_coefficients=mo.coefficients;mo.source_orbital_index=i;mo.energy_hartree=-1+i*.4;
        mo.occupation=i<2?2:i==2?1:0;mo.occupation_provenance=cov::DataProvenance::Producer;w.orbitals.push_back(mo);}
    data.id="ro-spatial-fixture";data.canonical_fingerprint=cov::nbo_canonical_fingerprint(w);data.dataset.association.compatible=true;
    cov::NboArchive ar;ar.basis_count=4;ar.open_shell=true;ar.density_is_bond_order=true;
    const auto add=[&](const char* name,cov::NboSpin spin,const M& x){cov::NboMatrix m;m.kind=name;m.spin=spin;m.rows=m.columns=4;m.values=flat(x);ar.matrices.push_back(m);};
    f.pa=eye;f.pa(3,3)=0;f.pb=f.pa;f.pb(2,2)=0;
    add("OVERLAP",cov::NboSpin::Total,eye);
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){add("LCAOMO",spin,eye);add("DENSITY",spin,spin==cov::NboSpin::Alpha?f.pa:f.pb);
        cov::NboCanonicalEvidence ev;ev.spin=spin;ev.direct_fchk_coefficients=spin==cov::NboSpin::Alpha;ev.density_verified=true;data.dataset.association.canonical_evidence.push_back(ev);}
    data.dataset.archive=ar;raw.available=true;raw.dataset_id=data.id;raw.canonical_fingerprint=data.canonical_fingerprint;raw.cache_key="raw";
    f.common=eye;const double a=.31;f.common(0,0)=f.common(3,3)=std::cos(a);f.common(3,0)=std::sin(a);f.common(0,3)=-std::sin(a);
    M rotation=eye;const double b=.63;rotation(0,0)=rotation(1,1)=std::cos(b);rotation(1,0)=std::sin(b);rotation(0,1)=-std::sin(b);
    f.fa=eye;f.fb=eye;for(int i=0;i<4;++i){f.fa(i,i)=1+i;f.fb(i,i)=5+2*i;}
    f.fa(0,1)=f.fa(1,0)=.17;f.fb(0,1)=f.fb(1,0)=.46;
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){
        const M basis=spin==cov::NboSpin::Alpha?f.common:M(f.common*rotation);
        const M ff=basis.transpose()*(spin==cov::NboSpin::Alpha?f.fa:f.fb)*basis;
        const M pp=basis.transpose()*(spin==cov::NboSpin::Alpha?f.pa:f.pb)*basis;
        cov::NboSalcSpinOperator op;op.spin=spin;op.fock=flat(ff);op.density=flat(pp);
        cov::NboSalcSubspace sub;sub.id=std::string(cov::nbo_spin_name(spin))+":space";sub.fragment_id="atom0";sub.spin=spin;sub.dimension=4;
        for(int i=0;i<4;++i){cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NAO,spin,std::size_t(i)};
            descriptor.id=std::string(cov::nbo_spin_name(spin))+":"+std::to_string(i);descriptor.label=descriptor.id;descriptor.atoms={0};
            descriptor.coefficients.resize(4);for(int k=0;k<4;++k)descriptor.coefficients[k]=basis(k,i);descriptor.orthonormal_basis=true;
            data.orbitals.push_back(descriptor);op.basis.push_back(descriptor.ref);
            cov::NboSalcOrbital o;o.id=descriptor.id;o.label=descriptor.label;o.type="Val(test)";o.spin=spin;o.atoms={0};o.fragment_id="atom0";o.subspace_id=sub.id;
            o.terms={{descriptor.ref,1}};o.energy_hartree=ff(i,i);o.occupation=pp(i,i);sub.orbital_indices.push_back(raw.orbitals.size());
            if(spin==cov::NboSpin::Alpha)for(int mo=0;mo<4;++mo){double c=basis(mo,i);raw.links.push_back({raw.orbitals.size(),std::size_t(mo),c,c*c});}
            raw.orbitals.push_back(o);}
        raw.subspaces.push_back(sub);raw.spin_operators.push_back(op);}
    for(std::size_t i=0;i<4;++i){w.gaussian_ao_transform.push_back({i,1,1});
        data.dataset.association.gaussian_row.push_back(i);data.dataset.association.coefficient_scale.push_back(1);}
    for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){
        cov::NboMatrix m;m.kind="FOCK";m.spin=spin;m.rows=m.columns=4;
        m.values=flat(spin==cov::NboSpin::Alpha?f.fa:f.fb);m.source.path="fixture/FILE.47";m.source.block="FOCK";
        data.dataset.archive->matrices.push_back(m);
        cov::NboSalcEnergyEvidence ev;ev.spin=spin;ev.available=true;ev.printed_operator_verified=true;
        ev.printed_nao_checked=ev.printed_nbo_checked=4;ev.source=m.source;raw.energies.push_back(ev);}
    data.canonical_fingerprint=cov::nbo_canonical_fingerprint(w);raw.canonical_fingerprint=data.canonical_fingerprint;
    return f;
}
}
int main(){try{
    auto f=fixture();const auto before=cov::serialize_nbo_salc_json(f.raw),fingerprint=f.data.canonical_fingerprint;
    const auto common=cov::build_nbo_ro_common_energy(f.w,f.data,f.raw);
    require(common.available&&common.orbitals.size()==4,"Common source energy qualification failed");
    require(!common.source_step_identity_verified,"Numerical equivalence incorrectly claimed original SCF identity");
    for(std::size_t i=0;i<4;++i){near(*common.orbitals[i].common_energy_hartree,(f.fa(i,i)+f.fb(i,i))/2,"Wrong source-MO common expectation");
        near(common.orbitals[i].source_energy_hartree,f.w.orbitals[i].energy_hartree,"Original source energy lost");}
    near(common.maximum_offdiagonal_hartree,.315,"Off-diagonal evidence lost; expectations mislabeled eigenvalues");
    require(cov::serialize_nbo_ro_common_energy_json(common).find("associated_archive_spin_average_expectation")!=std::string::npos,"Common operator semantics absent");
    auto bad_common=f;bad_common.raw.spin_operators[1].fock[0]+=.1;
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="operator_projection_mismatch","Stale same-fingerprint operator accepted");
    bad_common=f;bad_common.raw.energies[1].printed_nbo_checked=3;
    require(!cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).available,"Incomplete printed qualification accepted");
    bad_common=f;bad_common.data.dataset.archive->matrices.back().source.path="other/FILE.47";
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="operator_source_mismatch","Different analysis Focks averaged");
    bad_common=f;bad_common.data.dataset.association.gaussian_row[1]=0;
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="invalid_ao_mapping","Nonbijective mapping accepted");
    bad_common=f;bad_common.w.orbitals[0].coefficients[0]=.99;
    bad_common.data.canonical_fingerprint=cov::nbo_canonical_fingerprint(bad_common.w);bad_common.raw.canonical_fingerprint=bad_common.data.canonical_fingerprint;
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="canonical_linkage_mismatch","Compatible flag replaced source coefficient proof");
    bad_common=f;bad_common.raw.spin_operators[1].density[0]+=.1;
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="canonical_spin_density_mismatch","Projected density mismatch accepted");
    bad_common=f;bad_common.data.dataset.archive->matrices.back().values[1]+=.1;
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="nonhermitian_spin_fock","Non-Hermitian archive Fock accepted");
    bad_common=f;bad_common.data.dataset.archive->matrices.back().values.pop_back();
    require(cov::build_nbo_ro_common_energy(bad_common.w,bad_common.data,bad_common.raw).status=="complete_spin_fock_missing","Partial beta Fock accepted");
    auto averaged=cov::build_nbo_spin_averaged_model(f.w,f.data,f.raw);
    require(averaged.restricted_open_shell.verified&&averaged.spin_averaged&&averaged.orbitals.size()==4,"RO spatial merge missing members");
    require(averaged.merged_spatial_count==4&&averaged.separate_spin_count==0,"partial rotated space not completed");
    const M expected_f=f.common.transpose()*((f.fa+f.fb)/2)*f.common,expected_p=f.common.transpose()*(f.pa+f.pb)*f.common;
    require(averaged.subspaces.size()==1&&averaged.subspaces[0].fock.size()==16&&averaged.subspaces[0].density.size()==16,"Common subspace matrices missing");
    require((Eigen::Map<const RM>(averaged.subspaces[0].fock.data(),4,4)-expected_f).norm()<1e-10,"Subspace inherited alpha Fock instead of full common-basis mean");
    require((Eigen::Map<const RM>(averaged.subspaces[0].density.data(),4,4)-expected_p).norm()<1e-10,"Subspace inherited alpha density instead of common-basis total");
    require(!averaged.subspaces[0].energy_degeneracy_verified,"Split common spin operators falsely certified degenerate");
    bool fractional=false,rotated=false;
    for(std::size_t i=0;i<4;++i){const auto& o=averaged.orbitals[i];require(o.spatial_spin.has_value(),"Missing provenance");
        near(*o.energy_hartree,expected_f(i,i),"Averaged source diagonals instead of common Fock operator");
        near(*o.occupation,expected_p(i,i),"Occupation must sum after common-basis transform");
        fractional|=std::abs(*o.occupation-std::round(*o.occupation))>.01;rotated|=o.spatial_spin->dimension>1;
        auto selection=cov::nbo_salc_selection(averaged,i);auto view=cov::make_nbo_selection_view(f.data,f.w,selection);
        require(view.available&&selection.spatial_spin.has_value(),"Common basis selection unavailable");
        for(int row=0;row<4;++row)near(view.wavefunction.orbitals[0].coefficients[row],f.common(row,i),"Source field changed during merge");
        require(cov::serialize_nbo_selection_json(view).find("source_fock_matrix_hartree")!=std::string::npos,"Selection export lost mapping operators");
    }
    require(fractional&&rotated,"Fixture did not exercise fractional/rotated evidence");
    for(const auto& link:averaged.links){near(link.coefficient,f.common(link.canonical_index,link.side_index),"Common projection changed");near(link.weight,link.coefficient*link.coefficient,"Composition channels summed as percentages");}
    require(cov::serialize_nbo_salc_json(f.raw)==before&&cov::nbo_canonical_fingerprint(f.w)==fingerprint,"Source data mutated");
    auto energy_missing=f;energy_missing.raw.spin_operators[1].fock.clear();energy_missing.raw.spin_operators[1].energy_status="missing beta Fock";
    for(auto& o:energy_missing.raw.orbitals)if(o.spin==cov::NboSpin::Beta)o.energy_hartree.reset();
    auto unavailable=cov::build_nbo_spin_averaged_model(energy_missing.w,energy_missing.data,energy_missing.raw);
    require(unavailable.merged_spatial_count==4,"Missing energy incorrectly prevents proven spatial merge");
    for(const auto& sub:unavailable.subspaces)require(sub.fock.empty()&&!sub.energy_degeneracy_verified&&!sub.density.empty(),"Missing beta operator must clear quantitative common-Fock matrix without erasing density");
    for(const auto& o:unavailable.orbitals)require(!o.energy_hartree&&o.occupation&&o.spatial_spin->energy_status.find("missing beta Fock")!=std::string::npos,"Fabricated missing spin mean");
    for(const auto* method:{"UPBE1PBE","UHF","UKS","RHF","PBE0"}){auto wrong=f;wrong.w.source_route=std::string("SP ")+method+" synthetic";wrong.data.canonical_fingerprint=cov::nbo_canonical_fingerprint(wrong.w);
        require(!cov::build_nbo_ro_common_energy(wrong.w,wrong.data,wrong.raw).available,"Non-RO source given common RO energies");
        require(!cov::build_nbo_spin_averaged_model(wrong.w,wrong.data,wrong.raw).spin_averaged,"Unverified/U/closed-shell method merged");}
    auto invalid=f;invalid.data.dataset.association.compatible=false;require(!cov::verify_nbo_restricted_open_shell(invalid.w,invalid.data).verified,"Rejected association accepted");
    invalid=f;invalid.w.orbitals[2].occupation=.5;invalid.data.canonical_fingerprint=cov::nbo_canonical_fingerprint(invalid.w);require(!cov::verify_nbo_restricted_open_shell(invalid.w,invalid.data).verified,"Wrong integer occupation accepted");
    invalid=f;for(auto& m:invalid.data.dataset.archive->matrices)if(m.kind=="DENSITY"&&m.spin==cov::NboSpin::Beta)m.values[0]-=.1;
    require(!cov::verify_nbo_restricted_open_shell(invalid.w,invalid.data).verified,"Density evidence flag replaced actual shared-density proof");
    invalid=f;std::erase_if(invalid.data.dataset.archive->matrices,[](const auto& m){return m.kind=="LCAOMO"&&m.spin==cov::NboSpin::Beta;});
    require(!cov::verify_nbo_restricted_open_shell(invalid.w,invalid.data).verified,"Missing beta spatial evidence accepted");
    auto incomplete=f;incomplete.raw.subspaces[1].orbital_indices={4};
    const auto separated=cov::build_nbo_spin_averaged_model(incomplete.w,incomplete.data,incomplete.raw);
    require(separated.merged_spatial_count==2&&separated.separate_spin_count==4,"Different single local orbitals forcibly paired by label or energy");
    auto phase=f;for(auto& d:phase.data.orbitals)if(d.ref.spin==cov::NboSpin::Beta)for(auto& c:d.coefficients)c=-c;
    const auto phase_result=cov::build_nbo_spin_averaged_model(phase.w,phase.data,phase.raw);
    for(std::size_t i=0;i<4;++i)near(*phase_result.orbitals[i].energy_hartree,*averaged.orbitals[i].energy_hartree,"Local phase changes spin mean");
    // Total records are aliases, not a third population added to alpha+beta.
    auto total=f;auto alias=total.raw.orbitals[2];alias.spin=cov::NboSpin::Total;alias.id="total-alias";alias.occupation=99;
    total.raw.orbitals.push_back(alias);
    const auto total_result=cov::build_nbo_spin_averaged_model(total.w,total.data,total.raw);
    require(total_result.orbitals.size()==4&&total_result.orbitals[2].spatial_spin->total_aliases.size()==1,"Total record counted twice");
    near(*total_result.orbitals[2].occupation,*averaged.orbitals[2].occupation,"Total alias population added twice");
    total=f;cov::NboSalcSubspace totals;totals.spin=cov::NboSpin::Total;totals.id="total-space";
    for(int i=0;i<4;++i){auto o=total.raw.orbitals[i];o.id="total:"+std::to_string(i);o.spin=cov::NboSpin::Total;o.occupation=99;
        o.terms.clear();const double angle=.47;const int partner=(i+2)%4;
        o.terms.push_back({{cov::NboOrbitalKind::NAO,cov::NboSpin::Alpha,std::size_t(i)},std::cos(angle)});
        o.terms.push_back({{cov::NboOrbitalKind::NAO,cov::NboSpin::Alpha,std::size_t(partner)},i<2?std::sin(angle):-std::sin(angle)});
        totals.orbital_indices.push_back(total.raw.orbitals.size());total.raw.orbitals.push_back(o);}
    total.raw.subspaces.push_back(totals);
    const auto rotated_total=cov::build_nbo_spin_averaged_model(total.w,total.data,total.raw);
    require(rotated_total.orbitals.size()==4,"Rotated Total subspace duplicated");
    for(std::size_t i=0;i<4;++i){const auto& info=*rotated_total.orbitals[i].spatial_spin;
        require(info.total_aliases.size()==4&&info.total_alias_mapping.size()==4,"Total block provenance lost");
        Eigen::VectorXd reconstructed=Eigen::VectorXd::Zero(4);
        for(std::size_t j=0;j<4;++j)for(const auto& term:info.total_aliases[j].terms)
            reconstructed+=info.total_alias_mapping[j]*term.coefficient*f.common.col(term.orbital.index);
        require((reconstructed-f.common.col(i)).norm()<1e-10,"Total alias mapping not reproducible");
        near(*rotated_total.orbitals[i].occupation,*averaged.orbitals[i].occupation,"Rotated Total population added twice");}
    auto stale=f;stale.raw.dataset_id="previous-attachment";
    require(!cov::build_nbo_spin_averaged_model(stale.w,stale.data,stale.raw).spin_averaged,"Stale source model accepted");
    auto explicit_ro=f;explicit_ro.w.orbital_occupation_model=cov::OrbitalOccupationModel::ExplicitSpin;
    for(auto& o:explicit_ro.w.orbitals)o.occupation=o.source_orbital_index<3?1:0;
    for(int i=0;i<4;++i){auto o=explicit_ro.w.orbitals[i];o.spin=cov::Spin::Beta;o.occupation=i<2?1:0;explicit_ro.w.orbitals.push_back(o);}
    explicit_ro.data.canonical_fingerprint=cov::nbo_canonical_fingerprint(explicit_ro.w);explicit_ro.raw.canonical_fingerprint=explicit_ro.data.canonical_fingerprint;
    for(auto& ev:explicit_ro.data.dataset.association.canonical_evidence)ev.direct_fchk_coefficients=true;
    for(int i=0;i<4;++i){const auto& descriptor=explicit_ro.data.orbitals[4+i];for(int mo=0;mo<4;++mo){double c=descriptor.coefficients[mo];explicit_ro.raw.links.push_back({std::size_t(4+i),std::size_t(4+mo),c,c*c});}}
    const auto explicit_average=cov::build_nbo_spin_averaged_model(explicit_ro.w,explicit_ro.data,explicit_ro.raw);
    require(explicit_average.spin_averaged,"Explicit-spin RO not recognized");bool beta_component=false;
    const auto explicit_common=cov::build_nbo_ro_common_energy(explicit_ro.w,explicit_ro.data,explicit_ro.raw);
    require(explicit_common.available&&explicit_common.orbitals.size()==8,"Verified explicit-spin RO source energies unavailable");
    for(std::size_t i=0;i<4;++i)near(*explicit_common.orbitals[i].common_energy_hartree,*explicit_common.orbitals[i+4].common_energy_hartree,"Explicit RO channels do not share common operator");
    for(const auto& link:explicit_average.links)if(link.canonical_index>=4&&link.weight>.01){
        auto selection=cov::nbo_salc_component_selection(explicit_average,link);
        auto view=cov::make_nbo_selection_view(explicit_ro.data,explicit_ro.w,selection);
        require(view.available&&view.wavefunction.orbitals[0].spin==cov::Spin::Beta,"Common alpha rendering basis changed beta component identity");
        near(link.coefficient,f.common(link.canonical_index-4,link.side_index),"Signed beta projection not transformed to common basis");beta_component=true;}
    require(beta_component,"Beta component test not exercised");
    std::cout<<"nbo_spin_average_smoke passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
