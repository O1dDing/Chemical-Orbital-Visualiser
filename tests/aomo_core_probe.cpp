#include "cov/nbo_integration.hpp"
#include "cov/wavefunction_io.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv){try{
    if(argc!=4)throw std::runtime_error("Usage: cov_aomo_core_probe canonical.fchk result_directory output.json");
    const auto w=cov::parse_wavefunction(argv[1]);const auto before=w;const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1){std::cerr<<"candidate_count="<<found.candidates.size()<<'\n';for(const auto& c:found.candidates)std::cerr<<c.label<<'\n';throw std::runtime_error("Explicit unique analysis candidate required");}
    const auto integrated=cov::read_nbo_integration(w,found.candidates[0]);std::ofstream out(argv[3]);if(!out)throw std::runtime_error("Cannot write probe evidence");out<<std::setprecision(17)<<"{\"integration\":"<<cov::serialize_nbo_integration_json(integrated)<<",\"canonical\":{\"basis_count\":"<<w.basis_count<<",\"orbitals\":[";
    for(std::size_t i=0;i<w.orbitals.size();++i){const auto& m=w.orbitals[i];if(i)out<<',';out<<"{\"index\":"<<i<<",\"source_index\":"<<m.source_orbital_index<<",\"spin\":\""<<(m.spin==cov::Spin::Beta?"beta":"alpha")<<"\",\"occupation\":"<<m.occupation<<",\"coefficients\":[";for(std::size_t k=0;k<m.coefficients.size();++k){if(k)out<<',';out<<m.coefficients[k];}out<<"],\"gaussian_coefficients\":[";for(std::size_t k=0;k<m.gaussian_source_coefficients.size();++k){if(k)out<<',';out<<m.gaussian_source_coefficients[k];}out<<"]}";}out<<"],\"gaussian_transform\":[";for(std::size_t i=0;i<w.gaussian_ao_transform.size();++i){if(i)out<<',';const auto& t=w.gaussian_ao_transform[i];out<<"{\"source_index\":"<<t.source_index<<",\"scale\":"<<t.coefficient_scale<<'}';}out<<"]},\"reconstruction\":[";
    bool first=true;for(auto kind:{cov::NboOrbitalKind::GaussianAO,cov::NboOrbitalKind::NAO,cov::NboOrbitalKind::NHO,cov::NboOrbitalKind::NBO,cov::NboOrbitalKind::NLMO})for(std::size_t i=0;i<w.orbitals.size();++i){const auto links=cov::nbo_links_for_mo(integrated,i,kind);if(links.empty())continue;cov::NboOrbitalSelection s;s.dataset_id=integrated.id;s.mode=cov::NboSelectionMode::PartialSum;s.target_canonical_index=i;for(const auto& link:links)s.terms.push_back({link.orbital,link.coefficient});const auto view=cov::make_nbo_selection_view(integrated,w,s);if(!first)out<<',';first=false;out<<"{\"kind\":\""<<cov::nbo_orbital_kind_name(kind)<<"\",\"canonical_index\":"<<i<<",\"available\":"<<(view.available?"true":"false")<<",\"error\":";if(view.reconstruction_error)out<<*view.reconstruction_error;else out<<"null";out<<'}';}
    bool preserved=w.orbitals.size()==before.orbitals.size();for(std::size_t i=0;i<w.orbitals.size();++i)preserved=preserved&&w.orbitals[i].coefficients==before.orbitals[i].coefficients&&w.orbitals[i].gaussian_source_coefficients==before.orbitals[i].gaussian_source_coefficients&&w.orbitals[i].energy_hartree==before.orbitals[i].energy_hartree&&w.orbitals[i].occupation==before.orbitals[i].occupation;out<<"],\"canonical_preserved\":"<<(preserved?"true":"false")<<'}';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
