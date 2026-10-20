#include "cov/pi_field_response.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
using namespace cov;
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
NboMatrix zeros(std::size_t r,std::size_t c){NboMatrix m;m.rows=r;m.columns=c;m.values.assign(r*c,0);return m;}
double& at(NboMatrix& m,std::size_t r,std::size_t c){return m.values[r*m.columns+c];}
PiFieldResponseInput pair(double ligand_energy,double metal_occ,double ligand_occ){PiFieldResponseInput x;x.fock=zeros(2,2);at(x.fock,1,1)=ligand_energy;at(x.fock,0,1)=at(x.fock,1,0)=.2;
 x.metal_d_basis=zeros(2,1);at(x.metal_d_basis,0,0)=1;x.ligand_pi_basis=zeros(2,1);at(x.ligand_pi_basis,1,0)=1;x.total_density=zeros(2,2);at(x.total_density,0,0)=metal_occ;at(x.total_density,1,1)=ligand_occ;
 x.canonical_columns=zeros(2,2);at(x.canonical_columns,0,0)=at(x.canonical_columns,1,1)=1;x.canonical_indices={0,1};x.density_verified=x.shared_spatial_verified=true;x.operator_error_hartree=1e-12;return x;}
void rotate(NboMatrix& m,std::size_t a,std::size_t b,double angle){for(std::size_t r=0;r<m.rows;++r){double x=at(m,r,a),y=at(m,r,b);at(m,r,a)=std::cos(angle)*x+std::sin(angle)*y;at(m,r,b)=-std::sin(angle)*x+std::cos(angle)*y;}}
}
int main(){try{
 auto donor=assess_pi_field_response(pair(-1,2,2));require(donor.available&&donor.role==PiFieldRole::DonorDominant,"filled-filled ligand field donation lost");require(std::abs(donor.mean_shift_hartree-(std::sqrt(1.16)-1)/2)<1e-10,"two-level donor eigenvalue response wrong");
 auto ai=pair(1,1,0);ai.occupied_metal_backbond_evidence=true;auto acceptor=assess_pi_field_response(ai);require(acceptor.available&&acceptor.role==PiFieldRole::AcceptorDominant&&acceptor.occupied_metal_backbond_supported,"occupied-d acceptor response wrong");
 // Uniform energy scaling leaves eigenvectors and electron populations
 // unchanged. Hartree errors must not become projection/population cutoffs.
 auto scale_base=ai;scale_base.operator_error_hartree=1e-5;const auto unscaled=assess_pi_field_response(scale_base);auto scaled=scale_base;for(double& v:scaled.fock.values)v*=1e4;scaled.operator_error_hartree*=1e4;const auto scaled_result=assess_pi_field_response(scaled);
 require(unscaled.available&&scaled_result.available&&std::abs(unscaled.tracking.minimum_d_projection-scaled_result.tracking.minimum_d_projection)<1e-10,"energy scaling changed dimensionless d-reference availability");
 require(unscaled.occupied_metal_backbond_supported&&scaled_result.occupied_metal_backbond_supported,"energy scaling changed occupied-backbond qualification");
 auto population_noise=ai;at(population_noise.total_density,0,0)=1e-6;require(!assess_pi_field_response(population_noise).occupied_metal_backbond_supported,"density precision noise counted as occupied metal donor");
 auto empty=pair(1,0,0);auto d0=assess_pi_field_response(empty);require(d0.available&&d0.role==PiFieldRole::AcceptorDominant&&!d0.occupied_metal_backbond_supported,"empty metal falsely donates or loses acceptor field");
 auto inverted=assess_pi_field_response(pair(1,2,2));require(inverted.available&&inverted.role==PiFieldRole::Mixed,"occupied high ligand level mislabeled acceptor");
 auto no=pair(1,1,0);at(no.fock,0,1)=at(no.fock,1,0)=0;auto zero=assess_pi_field_response(no);require(zero.available&&zero.role==PiFieldRole::Negligible,"zero coupling failed");
 auto bad=pair(1,1,0);at(bad.ligand_pi_basis,0,0)=.2;require(!assess_pi_field_response(bad).available,"overlapping partitions accepted");
 auto nonhermitian=pair(1,1,0);at(nonhermitian.fock,0,1)+=.01;require(!assess_pi_field_response(nonhermitian).available,"non-Hermitian physical operator silently symmetrized");
 PiFieldResponseInput cancel;cancel.fock=zeros(3,3);at(cancel.fock,1,1)=-1;at(cancel.fock,2,2)=1;at(cancel.fock,0,1)=at(cancel.fock,1,0)=at(cancel.fock,0,2)=at(cancel.fock,2,0)=.2;
 cancel.metal_d_basis=zeros(3,1);at(cancel.metal_d_basis,0,0)=1;cancel.ligand_pi_basis=zeros(3,2);at(cancel.ligand_pi_basis,1,0)=at(cancel.ligand_pi_basis,2,1)=1;
 cancel.total_density=zeros(3,3);at(cancel.total_density,0,0)=1;at(cancel.total_density,1,1)=2;cancel.density_verified=true;auto cancelled=assess_pi_field_response(cancel);
 require(cancelled.available&&cancelled.role==PiFieldRole::Mixed&&std::abs(cancelled.mean_shift_hartree)<1e-10&&cancelled.gross_response_hartree>.01,"cancelling large fields mislabeled weak");
 PiFieldResponseInput multi;multi.fock=zeros(4,4);at(multi.fock,1,1)=.3;at(multi.fock,2,2)=-1;at(multi.fock,3,3)=1;at(multi.fock,0,2)=at(multi.fock,2,0)=.2;
 multi.metal_d_basis=zeros(4,2);at(multi.metal_d_basis,0,0)=at(multi.metal_d_basis,1,1)=1;multi.ligand_pi_basis=zeros(4,2);at(multi.ligand_pi_basis,2,0)=at(multi.ligand_pi_basis,3,1)=1;
 multi.total_density=zeros(4,4);at(multi.total_density,0,0)=1;at(multi.total_density,2,2)=2;multi.density_verified=true;multi.canonical_columns=zeros(4,4);for(int i=0;i<4;++i)at(multi.canonical_columns,i,i)=1;multi.canonical_indices={0,1,2,3};auto m=assess_pi_field_response(multi);
 require(m.available,"multi-d response failed");require(assess_pi_field_group(m,{0}).applicable,"responsive d group missing");require(!assess_pi_field_group(m,{1}).applicable,"nonresponsive sigma/eg group inherited global role");require(!assess_pi_field_group(m,{3}).applicable,"ligand-only group inherited global role");
 auto tail=m;tail.canonical_d_coordinates.values[3]=.08;tail.canonical_tracked_d_coordinates.values[3]=.01;require(!assess_pi_field_group(tail,{3}).applicable,"small high-virtual d tail inherited tracked-level field role");
 auto rotated=multi;rotate(rotated.metal_d_basis,0,1,.71);rotate(rotated.ligand_pi_basis,0,1,-.42);auto mr=assess_pi_field_response(rotated);
 require(mr.available&&mr.role==m.role&&std::abs(mr.trace_shift_hartree-m.trace_shift_hartree)<1e-10,"complete-subspace basis rotation changed role");
 require(assess_pi_field_group(mr,{0}).role==assess_pi_field_group(m,{0}).role&&!assess_pi_field_group(mr,{1}).applicable,"group response not covariant");
 // Environment orbitals are kept rather than reduced to the d+ligand block.
 auto env=multi;at(env.fock,2,3)=at(env.fock,3,2)=.7;env.ligand_pi_basis=zeros(4,1);at(env.ligand_pi_basis,2,0)=1;auto full=assess_pi_field_response(env);
 require(full.available&&full.full_environment_retained&&std::abs(full.trace_shift_hartree-m.trace_shift_hartree)>.001,"external environment dropped");
 // One d direction is spread over four environment eigenstates, while the
 // other has two 50%-d partners. Top individual weights give a rank-1 pair;
 // the complete reference must instead include both independent directions.
 PiFieldResponseInput spread;spread.fock=zeros(7,7);at(spread.fock,0,2)=at(spread.fock,2,0)=1;for(int i:{1,3,4,5})at(spread.fock,i,i)=3;
 at(spread.fock,1,3)=at(spread.fock,3,1)=at(spread.fock,3,4)=at(spread.fock,4,3)=at(spread.fock,4,5)=at(spread.fock,5,4)=1;at(spread.fock,6,6)=8;
 spread.metal_d_basis=zeros(7,2);at(spread.metal_d_basis,0,0)=at(spread.metal_d_basis,1,1)=1;spread.ligand_pi_basis=zeros(7,1);at(spread.ligand_pi_basis,6,0)=1;
 const auto complete=assess_pi_field_response(spread);require(complete.available&&complete.tracking.minimum_d_projection>.3,"initial reference duplicated one d direction and omitted another");
 std::cout<<"pi_field_response_smoke: all scientific invariance, occupation, cancellation and applicability checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
