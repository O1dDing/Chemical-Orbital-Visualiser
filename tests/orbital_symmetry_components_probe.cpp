#include "cov/orbital_symmetry_components.hpp"
#include "cov/wavefunction_io.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv){try {
    if(argc!=5)throw std::runtime_error("canonical.fchk attachment_directory_or_NONE canonical_or_salc zero_based_index");
    const auto w=cov::parse_wavefunction(argv[1]);const auto before=w;
    const bool attached=std::string(argv[2])!="NONE",canonical=std::string(argv[3])=="canonical";
    const std::size_t index=std::stoull(argv[4]);
    cov::NboIntegration data;cov::NboSalcModel salc;
    cov::ui::NboAomoNames names;
    if(attached) {
        const auto found=cov::discover_nbo_inputs({argv[2]});
        if(found.candidates.size()!=1)throw std::runtime_error("Expected one unambiguous attachment");
        data=cov::read_nbo_integration(w,found.candidates.front());salc=cov::build_nbo_salc_model(w,data);
        names=cov::ui::build_nbo_aomo_names(w,data,&salc);
    } else {data=cov::ui::canonical_component_dataset(w);names=*cov::ui::canonical_mo_names(w);}
    const auto& list=canonical?names.canonical:names.salc;
    if(index>=list.size())throw std::runtime_error("Source outside actual list");
    const auto& name=list[index];
    if(!name.decomposition_verified)throw std::runtime_error("No verified decomposition: "+name.decomposition_status);
    std::vector<double> target;
    if(canonical)target=w.orbitals[index].coefficients;
    else {
        const auto view=cov::make_nbo_selection_view(data,w,cov::nbo_salc_selection(salc,index));
        if(!view.available)throw std::runtime_error("Actual SALC cannot render");
        target=view.wavefunction.orbitals.front().coefficients;
    }
    const auto norm=[&](const std::vector<double>& x){double sum=0;const auto n=w.basis_count;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)sum+=x[i]*w.ao_overlap[i*n+j]*x[j];return sum;};
    std::vector<double> reconstructed(w.basis_count,0);double maximum_weight_error=0;
    std::cout<<std::setprecision(17)<<"{\"name\":"<<cov::ui::serialize_orbital_name_json(name)<<",\"components\":[";
    bool comma=false;
    for(std::size_t i=0;i<name.components.size();++i) {
        if(name.components[i].weight<=0)continue;
        const auto selection=cov::ui::symmetry_component_selection(w,attached?&data:nullptr,
            attached?&salc:nullptr,name,i,canonical?std::optional(index):std::nullopt,
            canonical?std::nullopt:std::optional(index));
        if(!selection)throw std::runtime_error("Cannot reconstruct component source references");
        const auto view=cov::make_nbo_selection_view(data,w,*selection);
        // Tiny components are still exported as source coefficients. The field
        // evaluator's positive-norm guard may legitimately decline numerical dust.
        if(!view.available) {
            if(name.components[i].weight>1e-10)throw std::runtime_error("Verified component render rejected: "+view.detail);
            continue;
        }
        const auto& coefficients=view.wavefunction.orbitals.front().coefficients;
        for(std::size_t j=0;j<coefficients.size();++j)reconstructed[j]+=coefficients[j];
        maximum_weight_error=std::max(maximum_weight_error,std::abs(norm(coefficients)/norm(target)-name.components[i].weight));
        if(comma)std::cout<<',';comma=true;std::cout<<cov::serialize_nbo_selection_json(view);
    }
    for(std::size_t j=0;j<target.size();++j)reconstructed[j]-=target[j];
    const double residual=std::sqrt(std::max(0.,norm(reconstructed)/norm(target)));
    bool unchanged=w.ao_overlap==before.ao_overlap;
    for(std::size_t i=0;i<w.orbitals.size();++i)unchanged=unchanged&&w.orbitals[i].coefficients==before.orbitals[i].coefficients&&
        w.orbitals[i].energy_hartree==before.orbitals[i].energy_hartree&&w.orbitals[i].occupation==before.orbitals[i].occupation;
    std::cout<<"],\"rendered_sum_residual\":"<<residual<<",\"maximum_weight_error\":"<<maximum_weight_error
        <<",\"canonical_immutable\":"<<(unchanged?"true":"false")<<'}';
    return unchanged&&residual<2e-5&&maximum_weight_error<2e-5?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
