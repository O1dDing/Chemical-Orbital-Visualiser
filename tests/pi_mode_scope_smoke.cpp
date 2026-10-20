#include "cov/pi_coupling.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
using namespace cov;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
NboMatrix matrix(std::size_t rows,std::size_t cols){NboMatrix m;m.rows=rows;m.columns=cols;m.values.assign(rows*cols,0);return m;}
double& at(NboMatrix& m,std::size_t r,std::size_t c){return m.values[r*m.columns+c];}
struct Fixture {NboPiCoupling channel;NboMatrix fock,canonical,localized;};
Fixture fixture(double first,double second,bool grouped=false){
 Fixture x;x.fock=matrix(4,4);x.canonical=matrix(4,4);x.localized=matrix(4,4);
 auto& c=x.channel;c.id="synthetic-independent-fragment-blocks";c.spin=NboSpin::Total;
 c.centre_family="valence-d";c.ligand_space_kind="internal-pi-antibonding";c.occupation_status="available";
 c.direction_mapping_error_bound=1e-10;c.operator_max_error_hartree=1e-12;
 c.centre_projector_basis=matrix(4,2);c.ligand_projector_basis=matrix(4,2);
 for(std::size_t b=0;b<2;++b){const auto d=2*b,a=d+1;const double t=b?second:first;
  at(x.localized,d,d)=at(x.localized,a,a)=1;at(c.centre_projector_basis,d,b)=1;at(c.ligand_projector_basis,a,b)=1;
  at(x.fock,d,d)=-.2;at(x.fock,a,a)=.2;at(x.fock,d,a)=at(x.fock,a,d)=t;
  const double theta=.5*std::atan2(t,.2),co=std::cos(theta),si=std::sin(theta);
  at(x.canonical,d,d)=co;at(x.canonical,a,d)=-si;at(x.canonical,d,a)=si;at(x.canonical,a,a)=co;
  NboPiDirectionProjection edge;edge.donor_nbo_id=d+1;edge.acceptor_nbo_id=a+1;edge.direction="centre_to_ligand";
  edge.energy_gap_hartree=.4;edge.fock_hartree=t;edge.status="available";c.direction_projection_evidence.push_back(edge);
 }
 if(grouped){NboPiCanonicalGroup lo,hi;lo.members={0,2};hi.members={1,3};lo.occupation_per_mo=2;hi.occupation_per_mo=0;c.groups={lo,hi};}
 else for(std::size_t i=0;i<4;++i){NboPiCanonicalGroup g;g.members={i};g.occupation_per_mo=i%2?0:2;c.groups.push_back(g);}
 populate_pi_coupling_modes(c,x.fock,x.canonical,&x.localized);return x;
}
void rotate_columns(NboMatrix& m,std::size_t a,std::size_t b,double angle){
 const double c=std::cos(angle),s=std::sin(angle);
 for(std::size_t r=0;r<m.rows;++r){const double x=at(m,r,a),y=at(m,r,b);at(m,r,a)=c*x+s*y;at(m,r,b)=-s*x+c*y;}
}
}
int main(){using namespace cov;
 auto independent=fixture(.07,.13);
 require(independent.channel.modes.size()==2,"independent singular blocks lost");
 const auto own=assess_pi_mode_pair(independent.channel,{0},{1});
 require(own.verified&&own.direction_verified&&own.direction=="centre_to_ligand"&&own.ordinary_display_eligible,"own-mode actual ordered edge rejected");
 require(own.two_endpoint_relation,"isolated dominant two-endpoint relation lost");
 const auto crossed=assess_pi_mode_pair(independent.channel,{0},{3});
 require(!crossed.verified&&!crossed.direction_verified&&!crossed.ordinary_display_eligible,"independent blocks acquired a Cartesian counterpart");
 PiEndpointSymmetryEvidence ag{true,"same-source-same-operations-alpha","Ag"},au=ag;au.irrep="Au";
 require(!assess_pi_mode_pair(independent.channel,{0},{1},ag,au).verified,"same-domain incompatible irreps accepted");
 au.scope_id="different-operation-domain";
 require(assess_pi_mode_pair(independent.channel,{0},{1},ag,au).verified,"labels from different domains used as a false veto");
 auto degenerate=fixture(.1,.1,true);require(degenerate.channel.modes.size()==1&&degenerate.channel.modes[0].rank==2,"degenerate mode split into arbitrary SVD vectors");
 const auto before=assess_pi_mode_pair(degenerate.channel,{0,2},{1,3});
 require(before.verified&&before.ordinary_display_eligible&&before.matched_edge_ids.size()==2,"complete degenerate same-edge group rejected");
 rotate_columns(degenerate.canonical,0,2,.317);rotate_columns(degenerate.canonical,1,3,-.811);
 rotate_columns(degenerate.channel.centre_projector_basis,0,1,.293);rotate_columns(degenerate.channel.ligand_projector_basis,0,1,-.517);
 populate_pi_coupling_modes(degenerate.channel,degenerate.fock,degenerate.canonical,&degenerate.localized);
 const auto after=assess_pi_mode_pair(degenerate.channel,{2,0},{3,1});
 require(after.verified&&after.direction==before.direction&&after.ordinary_display_eligible==before.ordinary_display_eligible&&after.matched_edge_ids==before.matched_edge_ids,"degenerate rotations changed evidence classification");
 require(std::abs(after.lower_coverage-before.lower_coverage)<1e-12&&std::abs(after.lower_role_coverage-before.lower_role_coverage)<1e-12&&std::abs(after.lower_cross_fock_mean_hartree-before.lower_cross_fock_mean_hartree)<1e-12,"complete-group trace not invariant under rotation");
 // Degenerate singular values alone do not connect two distinct ordered edges.
 auto same_sigma=fixture(.1,.1);
 const auto crossed_edge=assess_pi_mode_pair(same_sigma.channel,{0},{3});
 require(!crossed_edge.verified&&!crossed_edge.direction_verified&&crossed_edge.matched_edge_ids.empty()&&!crossed_edge.ordinary_display_eligible,"independent equal-singular-value blocks acquired a false relation");
 // Same-mode support cannot cross-wire a donor from one ordered source edge
 // to an acceptor from another, even when each separate role union is nonzero.
 auto unlinked=independent.channel;
 for(auto& mode:unlinked.modes){if(mode.ordered_edges.empty())continue;auto first=mode.ordered_edges.front(),second=first;
  first.id="edge-a";second.id="edge-b";first.groups.erase(std::remove_if(first.groups.begin(),first.groups.end(),[](const auto& g){return g.group_index!=0;}),first.groups.end());
  second.groups.erase(std::remove_if(second.groups.begin(),second.groups.end(),[](const auto& g){return g.group_index!=1;}),second.groups.end());mode.ordered_edges={first,second};}
 const auto no_edge=assess_pi_mode_pair(unlinked,{0},{1});
 require(no_edge.verified&&!no_edge.direction_verified&&no_edge.matched_edge_ids.empty(),"independent ordered role unions manufactured a direction");
 auto filled=independent.channel;filled.groups[1].occupation_per_mo=2;
 const auto occupied=assess_pi_mode_pair(filled,{0},{1});
 require(occupied.direction=="occupied_space_mixing"&&!occupied.ordinary_display_eligible,"two occupied endpoints incorrectly called donation/backdonation");
 auto auxiliary=independent.channel;auxiliary.centre_family="additional-pd-acceptors";
 for(auto& mode:auxiliary.modes)mode.centre_space_kind="auxiliary-p";
 const auto aux=assess_pi_mode_pair(auxiliary,{0},{1});
 require(aux.verified&&aux.direction_verified&&!aux.ordinary_display_eligible&&aux.channel_family=="auxiliary-p-pi","auxiliary pi either erased or promoted to d-pi primary");
 for(auto& mode:auxiliary.modes)mode.centre_space_kind="higher-radial-d";
 require(assess_pi_mode_pair(auxiliary,{0},{1}).channel_family=="higher-radial-pi","higher-radial d incorrectly called p coupling");
 auto transverse=independent.channel;transverse.ligand_space_kind="transverse-p-basis";
 require(!assess_pi_mode_pair(transverse,{0},{1}).ordinary_display_eligible,"transverse mathematical p basis promoted to primary chemical pi");
 auto lonepair=independent.channel;lonepair.ligand_space_kind="localized-pi-lone-pair";
 require(assess_pi_mode_pair(lonepair,{0},{1}).ordinary_display_eligible,"valid transverse lone-pair family globally excluded");
 auto missing=fixture(.07,.13);populate_pi_coupling_modes(missing.channel,missing.fock,missing.canonical);
 const auto no_source=assess_pi_mode_pair(missing.channel,{0},{1});
 require(no_source.verified&&!no_source.direction_verified&&!no_source.ordinary_display_eligible,"missing localized source manufactured a direction");
 // A complete group can span more than one nondegenerate fragment mode.
 auto multi=fixture(.07,.13,true);const auto whole=assess_pi_mode_pair(multi.channel,{0,2},{1,3});
 require(whole.verified&&whole.shared_mode_ids.size()==2&&whole.ordinary_display_eligible,"per-mode null directions incorrectly veto a complete multi-mode group");
 // Three principal canonical groups form one mode network. Keep that actual
 // network rather than inventing all of its pairwise energy-gap descriptors.
 NboPiCoupling network;network.id="three-group-network";network.spin=NboSpin::Total;
 network.centre_family="valence-d";network.ligand_space_kind="internal-pi-antibonding";
 network.occupation_status="available";network.direction_mapping_error_bound=1e-10;
 network.centre_projector_basis=matrix(3,1);network.ligand_projector_basis=matrix(3,1);
 at(network.centre_projector_basis,0,0)=1;at(network.ligand_projector_basis,1,0)=1;
 auto c=matrix(3,3),f=matrix(3,3),nbo=matrix(3,3);
 at(c,0,0)=std::sqrt(2./3);at(c,0,1)=at(c,0,2)=-1/std::sqrt(6.);
 at(c,1,1)=1/std::sqrt(2.);at(c,1,2)=-1/std::sqrt(2.);
 for(std::size_t j=0;j<3;++j){at(c,2,j)=1/std::sqrt(3.);at(nbo,j,j)=1;
  NboPiCanonicalGroup g;g.members={j};g.occupation_per_mo=j==2?0:2;network.groups.push_back(g);}
 const double energies[]{-.3,-.1,.4};
 for(std::size_t r=0;r<3;++r)for(std::size_t s=0;s<3;++s)for(std::size_t j=0;j<3;++j)at(f,r,s)+=at(c,r,j)*energies[j]*at(c,s,j);
 NboPiDirectionProjection source;source.donor_nbo_id=1;source.acceptor_nbo_id=2;source.direction="centre_to_ligand";
 source.status="available";source.energy_gap_hartree=at(f,1,1)-at(f,0,0);network.direction_projection_evidence={source};
 populate_pi_coupling_modes(network,f,c,&nbo);require(network.modes.size()==1,"network mode unavailable");
 const auto reduction=assess_pi_mode_pair(network,{1},{2});
 require(reduction.verified&&!reduction.two_endpoint_relation&&!reduction.ordinary_display_eligible,"three principal nodes expanded into an ordinary two-endpoint relation");
 const auto graph=assess_pi_mode_network(network,network.modes.front());
 require(graph.verified&&graph.nodes.size()==3&&!graph.two_endpoint_relation&&graph.direction_verified&&graph.ordinary_display_eligible,"real multi-group direction discarded instead of retained as a network");
 const auto encoded=pi_mode_network_assessment_json(graph);
 require(encoded.find("three-group-network")!=std::string::npos&&encoded.find("\"members\":[2]")!=std::string::npos&&encoded.find("\"matched_edge_ids\":[\"total:1:2\"]")!=std::string::npos,"network source identities lost in compact serialization");
 std::cout<<"pi mode scope: independent-block, complete-group, rotation, ordered-edge, symmetry, family and occupied-space checks passed\n";
}
