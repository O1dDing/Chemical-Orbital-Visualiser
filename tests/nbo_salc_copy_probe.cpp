#include "cov/nbo_salc.hpp"
#include "cov/wavefunction_io.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <map>

namespace {
using M=Eigen::MatrixXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
double error(const M& a){return a.size()?a.cwiseAbs().maxCoeff():0;}
}
int main(int argc,char** argv){try{
    require(argc==4,"Usage: cov_nbo_salc_copy_probe canonical.fchk analysis_directory output.json");
    const auto w=cov::parse_wavefunction(argv[1]);const auto fingerprint=cov::nbo_canonical_fingerprint(w);
    const auto inputs=cov::discover_nbo_inputs({argv[2]});require(inputs.candidates.size()==1,"Expected one associated analysis input");
    const auto data=cov::read_nbo_integration(w,inputs.candidates.front());
    const auto model=cov::build_nbo_salc_model(w,data);require(model.available,"SALC model unavailable");
    const M metric=Eigen::Map<const RM>(w.ao_overlap.data(),w.basis_count,w.basis_count);
    std::map<std::string,std::vector<M>> families;std::size_t copies=0,repeated=0;
    for(const auto& sub:model.subspaces){
        require(sub.source_to_derived.size()==sub.source_basis.size()*sub.dimension,"Missing explicit source transform");
        const M u=Eigen::Map<const RM>(sub.source_to_derived.data(),sub.source_basis.size(),sub.dimension);
        M source(w.basis_count,sub.source_basis.size());
        for(std::size_t j=0;j<sub.source_basis.size();++j){const auto* orbital=cov::nbo_orbital(data,sub.source_basis[j]);require(orbital!=nullptr,"Source reference lost");source.col(j)=Eigen::Map<const Eigen::VectorXd>(orbital->coefficients.data(),w.basis_count);}
        const M q=source*u;require(error(q.transpose()*metric*q-M::Identity(sub.dimension,sub.dimension))<=5e-5,"Derived Gram invalid");
        for(std::size_t j=0;j<sub.dimension;++j){const auto& orbital=model.orbitals[sub.orbital_indices[j]];
            Eigen::VectorXd actual=Eigen::VectorXd::Zero(w.basis_count);
            for(const auto& term:orbital.terms)actual+=term.coefficient*Eigen::Map<const Eigen::VectorXd>(cov::nbo_orbital(data,term.orbital)->coefficients.data(),w.basis_count);
            require((actual-q.col(j)).norm()<1e-10,"Selection terms do not implement declared mapping");
            if(!sub.fock.empty())require(orbital.energy_hartree&&std::abs(*orbital.energy_hartree-sub.fock[j*sub.dimension+j])<1e-9,"Energy does not equal transformed full Fock diagonal");
            if(!sub.density.empty())require(orbital.occupation&&std::abs(*orbital.occupation-sub.density[j*sub.dimension+j])<1e-9,"Occupation does not equal transformed density diagonal");
            for(const auto& link:model.links)if(link.side_index==sub.orbital_indices[j]){
                const auto& canonical=w.orbitals[link.canonical_index];const double amplitude=q.col(j).dot(metric*Eigen::Map<const Eigen::VectorXd>(canonical.coefficients.data(),w.basis_count));
                require(std::abs(amplitude-link.coefficient)<1e-9,"Signed link does not follow mapped derived field");}
        }
        if(sub.symmetry_verified&&sub.irrep_dimension>1){if(sub.multiplicity==1)++copies;else ++repeated;}
    }
    require(cov::nbo_canonical_fingerprint(w)==fingerprint,"Source canonical changed");
    std::ofstream output(argv[3]);require(bool(output),"Cannot write output");output<<cov::serialize_nbo_salc_json(model);output.close();
    std::cout<<"Verified multidimensional copies="<<copies<<", unresolved repeated spaces="<<repeated<<", source unchanged\n";
    require(repeated==0,"Repeated multidimensional spaces remain unresolved");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
