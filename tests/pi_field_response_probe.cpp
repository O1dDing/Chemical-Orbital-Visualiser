#include "cov/pi_field_response.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/wavefunction_io.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iomanip>
namespace {
void matrix(std::ostream& out,const cov::NboMatrix& m){out<<"{\"rows\":"<<m.rows<<",\"columns\":"<<m.columns<<",\"values\":[";for(std::size_t i=0;i<m.values.size();++i){if(i)out<<',';out<<m.values[i];}out<<"]}";}
void rotate(cov::NboMatrix& m,double angle){for(std::size_t j=1;j<m.columns;++j)for(std::size_t i=0;i<m.rows;++i){auto& a=m.values[i*m.columns+j-1];auto& b=m.values[i*m.columns+j];const double x=a,y=b;a=std::cos(angle)*x+std::sin(angle)*y;b=-std::sin(angle)*x+std::cos(angle)*y;}}
void sensitivity(std::ostream& out,const cov::PiFieldResponseAnalysis& result){out<<'[';bool first=true;for(double threshold:{.4,.5,.6})for(auto response:result.responses){if(!first)out<<',';first=false;response.minimum_group_spectral_fraction=threshold;out<<"{\"minimum_group_spectral_fraction\":"<<threshold<<",\"groups\":[";for(std::size_t i=0;i<response.groups.size();++i){if(i)out<<',';out<<cov::pi_field_group_assessment_json(cov::assess_pi_field_group(response,response.groups[i].members));}out<<"]}";}out<<']';}
}
int main(int argc,char** argv){try{
 if(argc<4||argc>5)throw std::runtime_error("Usage: cov_pi_field_response_probe canonical.fchk analysis_directory output.json [rotate-projectors|negative-fock]");
 const auto w=cov::parse_wavefunction(argv[1]);const auto identity=cov::nbo_canonical_fingerprint(w);
 const auto found=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});if(found.candidates.size()!=1)throw std::runtime_error("Exactly one source analysis required");
 auto data=cov::read_nbo_integration(w,found.candidates.front());cov::annotate_nbo_bond_channels(data,w);
 const auto route=cov::route_chemistry(w,&data);std::vector<const cov::NboPiCoupling*> channels;for(const auto& r:route.pi_couplings)if(r.available())channels.push_back(&*r.value);
 std::vector<cov::NboPiCoupling> rotated;if(argc==5&&std::string(argv[4])=="rotate-projectors"){for(const auto* channel:channels){rotated.push_back(*channel);rotate(rotated.back().centre_projector_basis,.731);rotate(rotated.back().ligand_projector_basis,-.417);}channels.clear();for(const auto& c:rotated)channels.push_back(&c);}
 if(argc==5&&std::string(argv[4])=="negative-fock"&&data.dataset.archive){for(auto& matrix:data.dataset.archive->matrices)if(matrix.kind=="FOCK"&&matrix.columns>1){matrix.values[1]+=.1;break;}}
 const auto result=cov::analyse_pi_field_responses(w,data,channels);std::ofstream output(argv[3]);if(!output)throw std::runtime_error("Cannot open output");output<<std::setprecision(17);
 output<<"{\"canonical_immutable\":"<<(identity==cov::nbo_canonical_fingerprint(w)?"true":"false")<<",\"analysis\":"<<cov::pi_field_response_analysis_json(result,true)<<",\"mapping_sensitivity\":";sensitivity(output,result);
 output<<",\"scope_geometry\":[";bool first=true;for(const auto* c:channels){if(c->centre_family!="valence-d")continue;if(!first)output<<',';first=false;output<<"{\"ligand_rank\":"<<c->ligand_projector_basis.columns<<",\"spin\":\""<<cov::nbo_spin_name(c->spin)<<"\",\"localized_family_verified\":"<<(c->localized_family_verified?"true":"false")<<",\"angular_leakage\":"<<c->angular_leakage<<",\"angular_evidence\":[";for(std::size_t i=0;i<c->angular_projector_evidence.size();++i){if(i)output<<',';const auto& e=c->angular_projector_evidence[i];output<<"{\"atom\":"<<e.atom<<",\"max_sigma_leakage\":"<<e.max_sigma_leakage<<",\"partition_residual\":"<<e.partition_residual<<"}";}output<<']';if(argc==5&&std::string(argv[4])=="diagnostic-matrices"){output<<",\"centre_projector_basis\":";matrix(output,c->centre_projector_basis);output<<",\"ligand_projector_basis\":";matrix(output,c->ligand_projector_basis);}output<<'}';}output<<"]}";
 return output?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
