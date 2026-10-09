#include "cov/pi_coupling.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace {
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
}
int main(){try{
    // Deliberate missing-descriptor fault after valid transform/operator gates.
    // This fixture tests failure semantics, not the chemistry of a molecule.
    cov::Wavefunction w;w.basis_count=6;
    w.atoms={{"C",6,0,0,0,6},{"C",6,0,0,2,6}};
    w.ao_overlap.assign(36,0);
    cov::NboIntegration d;d.dataset.association.compatible=true;
    cov::NboMatrix eye;eye.rows=eye.columns=6;eye.values.assign(36,0);
    for(std::size_t k=0;k<6;++k){
        eye.values[k*6+k]=1;w.ao_overlap[k*6+k]=1;
        cov::MolecularOrbital mo;mo.source_orbital_index=k;mo.energy_hartree=1;
        mo.coefficients.assign(6,0);mo.coefficients[k]=1;w.orbitals.push_back(mo);
        cov::NboNao row;row.id=k+1;row.atom=k/3+1;row.type="Val(2p)";
        row.angular=k%3==0?"px":k%3==1?"py":"pz";d.dataset.naos.push_back(row);
    }
    cov::NboArchive archive;
    for(const char* kind:{"OVERLAP","FOCK"}){eye.kind=kind;archive.matrices.push_back(eye);}
    d.dataset.archive=archive;
    for(const char* kind:{"AONAO","NAOMO"}){eye.kind=kind;d.dataset.matrices.push_back(eye);}
    cov::NboNaoValidation v;v.available=v.direct_fchk_coefficients=true;v.effective_mo_columns=6;
    d.dataset.nao_validation.push_back(v);d.canonical_fingerprint=cov::nbo_canonical_fingerprint(w);
    const auto absent=cov::analyse_nbo_pi_couplings(w,d,{});
    require(absent.status=="not_applicable"&&!absent.evidence.empty(),"no-candidate control lost operator evidence");
    const auto failed=cov::analyse_nbo_pi_couplings(w,d,{{0,1}});
    require(!failed.evidence.empty(),"fault did not reach the angular gate");
    require(failed.status=="insufficient_evidence","missing angular evidence misreported as not applicable");
    require(failed.reason.find("descriptor")!=std::string::npos,"specific angular failure reason lost");
    // Exact two-level fixtures test the diagnostic, not a fitted chemistry label.
    const auto frozen=[](double gap,double coupling,bool merged=false){
        cov::NboMatrix f;f.rows=f.columns=2;f.values={0,coupling,coupling,gap};
        cov::NboMatrix m;m.rows=2;m.columns=1;m.values={1,0};
        auto l=m;l.values={0,1};
        const double theta=0.5*std::atan2(2*coupling,gap),c=std::cos(theta),s=std::sin(theta);
        cov::NboMatrix eigen;eigen.rows=eigen.columns=2;eigen.values={c,s,-s,c};
        return cov::assess_pi_frozen_operator(f,m,l,eigen,
            merged?std::vector<std::vector<std::size_t>>{{0,1}}:std::vector<std::vector<std::size_t>>{{0},{1}},
            {2,0},1e-12);
    };
    const auto small=frozen(1,.001);
    require(small.available&&small.tracking_verified&&small.occupation_boundary_preserved,"isolated exact block tracking");
    require(std::abs(small.max_energy_shift_hartree-(std::sqrt(1.000004)-1)/2)<1e-10,"maximum frozen spectral shift");
    require(small.max_subspace_sin2<.000002,"small mixing subspace angle");
    const auto resonant=frozen(.000001,.001);
    require(resonant.max_subspace_sin2>.49,"small coupling near resonance retained");
    const auto merged=frozen(.001,.1,true);
    require(merged.max_subspace_sin2<1e-12&&merged.max_group_width_change_hartree>.19,"full-space invariance cannot hide internal splitting");
    cov::NboMatrix sf;sf.rows=sf.columns=4;
    sf.values={0,.001,0,0,.001,1,0,0,0,0,3,.001,0,0,.001,3.000001};
    cov::NboMatrix sm;sm.rows=4;sm.columns=2;sm.values={1,0,0,0,0,1,0,0};
    auto sl=sm;sl.values={0,0,1,0,0,0,0,1};
    const double t1=.5*std::atan2(.002,1.),t2=.5*std::atan2(.002,.000001);
    cov::NboMatrix sc;sc.rows=sc.columns=4;
    sc.values={std::cos(t1),std::sin(t1),0,0,-std::sin(t1),std::cos(t1),0,0,
        0,0,std::cos(t2),std::sin(t2),0,0,-std::sin(t2),std::cos(t2)};
    cov::NboPiCoupling scope;
    scope.frozen_operator=cov::assess_pi_frozen_operator(sf,sm,sl,sc,{{0},{1},{2},{3}},{2,0,0,0},1e-12,&scope.frozen_groups);
    const auto local=cov::assess_pi_frozen_display_scope(scope,{0,1});
    require(local.available&&local.tracking_verified&&local.max_subspace_sin2<.000002,
        "unrelated high virtual resonance cannot redefine the fixed low display domain");
    require(local.full_spectrum_max_energy_shift_hartree>local.max_energy_shift_hartree,
        "full-spectrum diagnostic retained beside scoped maximum");
    auto partial=scope;partial.frozen_groups[0].members={0,1};partial.frozen_groups.erase(partial.frozen_groups.begin()+1);
    require(!cov::assess_pi_frozen_display_scope(partial,{0}).available,"display domain cannot cut a complete degenerate cluster");
    auto overlap=eye;overlap.rows=overlap.columns=2;overlap.values={1,0,0,1};
    require(!cov::assess_pi_frozen_operator(overlap,overlap,overlap,overlap,{{0},{1}},{2,0},0).available,
        "overlapping projectors cannot be subtracted twice");
    cov::NboPiCoupling calibrated;
    calibrated.localized_family_verified=true;calibrated.occupation_status="available";
    calibrated.ligand_family="internal-pi-bonding";calibrated.centre_family="additional-pd-acceptors";
    calibrated.ligand_family_nbo_ids={1};calibrated.centre_projector_basis=overlap;
    calibrated.ligand_projector_basis=overlap;
    const auto budget=cov::pi_channel_display_calibration(calibrated);
    require(budget.validated&&budget.energy_budget_ev==.025&&budget.subspace_sin2_budget==.02,
        "saved-family display calibration budget");
    require(cov::pi_display_negligible(small,budget)&&!cov::pi_display_negligible(resonant,budget),
        "validated budget distinguishes small isolated mixing from resonance");
    require(!cov::pi_display_negligible(calibrated.frozen_operator,budget),
        "structural eligibility alone cannot replace missing frozen operator");
    calibrated.ligand_family.clear();
    require(!cov::pi_channel_display_calibration(calibrated).validated,"broad family outside calibration");
    calibrated.ligand_family="internal-pi-antibonding";calibrated.operator_kind="physical-spin-fock";
    require(!cov::pi_channel_display_calibration(calibrated).validated,"different operator outside calibration");
    std::cout<<"NBO pi evidence: no candidate and failed angular projection remain distinct\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
