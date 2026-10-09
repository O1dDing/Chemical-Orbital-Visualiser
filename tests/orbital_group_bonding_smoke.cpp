#include "cov/orbital_group_bonding.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
void check(bool x,const char* s){if(!x)throw std::runtime_error(s);}
void close(double a,double b,const char* s){check(std::abs(a-b)<2e-10,s);}
cov::Wavefunction fixture(){
    cov::Wavefunction w;w.basis_count=4;
    for(unsigned i=0;i<4;++i){cov::Atom a;a.symbol="H";a.atomic_number=1;a.x=i*2.;w.atoms.push_back(a);cov::Shell shell;shell.atom_index=i;shell.basis_offset=i;shell.angular_momentum=0;w.shells.push_back(shell);}
    w.ao_overlap={1,.05,0,0,.05,1,0,0,0,0,1,.3,0,0,.3,1};
    for(unsigned i=0;i<2;++i){cov::MolecularOrbital m;m.coefficients.assign(4,0);m.occupation=i?.4:2.;m.occupation_provenance=cov::DataProvenance::Producer;m.energy_hartree=-.3+i*1e-6;
        const auto a=i*2;const double s=i?.3:.05,sign=i?1.:-1.;
        m.coefficients[a]=1/std::sqrt(2*(1+sign*s));m.coefficients[a+1]=sign*m.coefficients[a];w.orbitals.push_back(m);}
    return w;
}
}
int main(){try{
    auto w=fixture();const auto original=w.orbitals;
    const auto scope=cov::make_fixed_bonding_scope(w,{{0,1},{2,3},{1,0}},"all_two_fragments");
    check(scope.edges.size()==2,"fixed pairs not deduplicated");
    cov::OrbitalGroupBondingOptions options;options.expected_dimension=2;
    const auto base=cov::analyse_orbital_group_bonding(w,scope,{0,1},options);
    check(base.status==cov::OrbitalGroupBondingStatus::Mixed && base.trace>0,"positive trace hid opposite-sign spectrum");
    close(base.eigenvalues[0],-.05/.95,"negative eigenvalue incorrect");close(base.eigenvalues[1],.3/1.3,"positive eigenvalue incorrect");
    check(base.occupation_available,"source occupation missing");
    const double c=std::cos(.731),s=std::sin(.731);
    for(unsigned mu=0;mu<4;++mu){w.orbitals[0].coefficients[mu]=c*original[0].coefficients[mu]+s*original[1].coefficients[mu];w.orbitals[1].coefficients[mu]=-s*original[0].coefficients[mu]+c*original[1].coefficients[mu];}
    const std::vector<double> occupation{2*c*c+.4*s*s,(-2+.4)*c*s,(-2+.4)*c*s,2*s*s+.4*c*c};
    const auto rotated=cov::analyse_orbital_group_bonding(w,scope,{0,1},options,occupation);
    for(unsigned i=0;i<2;++i)close(rotated.eigenvalues[i],base.eigenvalues[i],"group rotation changed spectrum");
    close(rotated.trace,base.trace,"rotation changed trace");close(rotated.electron_weighted_trace,base.electron_weighted_trace,"occupation matrix not rotated with unequal occupations");
    close(rotated.electron_count,base.electron_count,"rotation changed electron count");
    for(double& x:w.orbitals[0].coefficients)x=-x;
    auto phased_occupation=occupation;phased_occupation[1]*=-1;phased_occupation[2]*=-1;
    const auto phased=cov::analyse_orbital_group_bonding(w,scope,{0,1},options,phased_occupation);
    close(phased.eigenvalues[0],base.eigenvalues[0],"independent phase changed spectrum");close(phased.electron_weighted_trace,base.electron_weighted_trace,"phase changed physical density weighting");
    const auto exchanged=cov::analyse_orbital_group_bonding(w,scope,{1,0},options,{phased_occupation[3],phased_occupation[2],phased_occupation[1],phased_occupation[0]});
    close(exchanged.eigenvalues[1],base.eigenvalues[1],"member exchange changed spectrum");
    w.orbitals=original;for(double& x:w.orbitals[0].coefficients)x*=2;for(double& x:w.orbitals[1].coefficients)x*=.5;
    const auto rescaled=cov::analyse_orbital_group_bonding(w,scope,{0,1},options,{.5,0,0,1.6});
    close(rescaled.eigenvalues[0],base.eigenvalues[0],"metric normalization changed physical spectrum");close(rescaled.electron_weighted_trace,base.electron_weighted_trace,"metric occupation transform incorrect");
    w.orbitals=original;for(auto& mo:w.orbitals)mo.occupation=0;
    const auto virtual_group=cov::analyse_orbital_group_bonding(w,scope,{0,1},options);
    check(virtual_group.status==cov::OrbitalGroupBondingStatus::Mixed,"virtual group was called nonbonding");close(virtual_group.electron_weighted_trace,0,"virtual occupied trace not zero");
    const auto local=cov::make_fixed_bonding_scope(w,{{2,3}},"local_2_3");
    const auto local_result=cov::analyse_orbital_group_bonding(w,local,{0,1},options);
    check(local_result.scope_id!=base.scope_id && local_result.status==cov::OrbitalGroupBondingStatus::Unresolved,"local rank-zero direction inherited global mixed label");
    check(cov::analyse_orbital_group_bonding(w,scope,{0},options).status==cov::OrbitalGroupBondingStatus::IncompleteGroup,"incomplete group certified");
    w.orbitals[1].coefficients=w.orbitals[0].coefficients;
    check(cov::analyse_orbital_group_bonding(w,scope,{0,1},options).status==cov::OrbitalGroupBondingStatus::IncompleteGroup,"rank-deficient source group certified");
    w.orbitals=original;w.ao_overlap.clear();check(cov::analyse_orbital_group_bonding(w,scope,{0,1},options).status==cov::OrbitalGroupBondingStatus::Unavailable,"missing metric became zero indicator");
    w=fixture();w.atoms.resize(1);check(cov::analyse_orbital_group_bonding(w,scope,{0},{}).status==cov::OrbitalGroupBondingStatus::NotApplicable,"single atom acquired interatomic label");
    w=fixture();w.ao_overlap[1]=w.ao_overlap[4]=.3;w.orbitals[0].coefficients={1/std::sqrt(2.6),1/std::sqrt(2.6),0,0};w.orbitals[1].coefficients={0,0,1/std::sqrt(2.6),1/std::sqrt(2.6)};
    const auto repeated=cov::analyse_orbital_group_bonding(w,scope,{0,1},options);
    check(repeated.dimension==2 && repeated.status==cov::OrbitalGroupBondingStatus::Positive,"repeated independent copies collapsed");close(repeated.eigenvalues[0],repeated.eigenvalues[1],"equal repeated copy spectrum changed");
    cov::InteractionGraph graph;cov::InteractionEdge ordinary;ordinary.atom_a=0;ordinary.atom_b=1;ordinary.kind=cov::InteractionKind::CovalentConnectivity;ordinary.strength=cov::InteractionStrength::StrongConnectivity;graph.edges.push_back(ordinary);
    auto weak=ordinary;weak.atom_a=2;weak.atom_b=3;weak.kind=cov::InteractionKind::HydrogenBond;weak.strength=cov::InteractionStrength::WeakContact;graph.edges.push_back(weak);
    check(cov::make_fixed_bonding_scope(w,graph).edges.size()==1,"weak relation entered ordinary fixed skeleton");
    std::cout<<"orbital_group_bonding_smoke ok\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
