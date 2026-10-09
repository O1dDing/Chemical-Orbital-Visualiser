#include "cov/orbital_symmetry_components.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool condition,const char* why){if(!condition)throw std::runtime_error(why);}
int main(){try {
    cov::Wavefunction w;w.basis_count=2;w.ao_overlap={2,0,0,3};
    w.atoms={{"H",1,0,0,0}};
    const double a=std::sqrt(.7),b=std::sqrt(.3);
    cov::MolecularOrbital m;m.coefficients={a/std::sqrt(2.),b/std::sqrt(3.)};
    m.energy_hartree=-.4;m.occupation=2;m.occupation_provenance=cov::DataProvenance::Producer;
    w.orbitals.push_back(m);m.coefficients={-b/std::sqrt(2.),a/std::sqrt(3.)};
    m.energy_hartree=.3;m.occupation=0;w.orbitals.push_back(m);
    const auto before=w;
    cov::ui::NboAomoName name;name.decomposition_verified=true;name.point_group="Cs";
    name.component_source_kind="canonical";name.component_source_members={0,1};
    name.components={{"A'",1,.7,{a*a,-a*b}},{"A''",1,.3,{b*b,a*b}}};
    const auto data=cov::ui::canonical_component_dataset(w);
    require(data.capabilities.empty(),"canonical adapter invented NBO capabilities");
    std::vector<double> sum(2,0);
    for(std::size_t i=0;i<2;++i) {
        const auto selection=cov::ui::symmetry_component_selection(w,nullptr,nullptr,name,i,0);
        require(selection && !selection->normalize && selection->target_canonical_index==0,"component identity or norm changed");
        const auto view=cov::make_nbo_selection_view(data,w,*selection);
        require(view.available && view.wavefunction.orbitals.size()==1,"non-NBO component cannot render");
        const auto& c=view.wavefunction.orbitals.front().coefficients;
        require(std::abs(view.metric_norm2.front()-name.components[i].weight)<1e-12,"component weight not reproduced by field");
        for(std::size_t j=0;j<2;++j)sum[j]+=c[j];
        const auto attached=cov::ui::symmetry_component_selection(w,&data,nullptr,name,i,0);
        require(attached && attached->terms.size()==selection->terms.size(),"attached/nonattached source terms differ");
        require(cov::ui::is_symmetry_component(*attached),"component semantic missing");
        const auto json=cov::serialize_nbo_selection_json(view);
        require(json.find("canonical_symmetry_component")!=std::string::npos,"export lost component identity");
    }
    for(std::size_t i=0;i<2;++i)require(std::abs(sum[i]-w.orbitals[0].coefficients[i])<1e-12,"component sum changed source MO");
    auto stale=data;stale.canonical_fingerprint="wrong";
    const auto selected=cov::ui::symmetry_component_selection(w,&stale,nullptr,name,0,0);
    require(!cov::make_nbo_selection_view(stale,w,*selected).available,"stale dataset rendered");
    name.decomposition_verified=false;
    require(!cov::ui::symmetry_component_selection(w,&data,nullptr,name,0,0),"unverified decomposition rendered");
    name.decomposition_verified=true;name.component_source_members[1]=22;
    require(!cov::ui::symmetry_component_selection(w,&data,nullptr,name,0,0),"bad source member accepted");
    for(std::size_t i=0;i<2;++i)require(w.orbitals[i].coefficients==before.orbitals[i].coefficients &&
        w.orbitals[i].energy_hartree==before.orbitals[i].energy_hartree &&
        w.orbitals[i].occupation==before.orbitals[i].occupation,"source wavefunction modified");
    std::cout<<"symmetry component field, export, fallback, stale and rejected evidence passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
