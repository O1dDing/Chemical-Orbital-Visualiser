#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
void close(double actual,double expected,const char* message) {
    require(std::abs(actual-expected)<1e-12,message);
}
cov::NboNaoContribution row(std::size_t atom,const char* type,int n,int l,double weight) {
    cov::NboNaoContribution r;r.atom=atom;r.type=type;r.principal_n=n;r.angular_l=l;r.weight=weight;return r;
}
cov::RoutedAnalysis route(const cov::Wavefunction& wf,
    const std::vector<std::vector<cov::NboNaoContribution>>& rows) {
    cov::RoutedAnalysis r;r.canonical_fingerprint=cov::nbo_canonical_fingerprint(wf);
    for(std::size_t i=0;i<rows.size();++i) {
        cov::RoutedResult<cov::RoutedMoComposition> entry;
        entry.status=cov::RoutedStatus::Available;entry.provider=cov::RoutedProvider::Nbo;
        cov::RoutedMoComposition c;c.canonical_index=i;c.complete=true;c.rows=rows[i];entry.value=c;
        r.mo_composition.push_back(entry);
    }
    return r;
}
}

int main() {try {
    cov::Wavefunction wf;wf.atoms.resize(2);wf.atoms[0].atomic_number=24;wf.atoms[1].atomic_number=17;
    wf.orbitals.resize(2);wf.orbitals[0].coefficients={1,0};wf.orbitals[1].coefficients={0,1};
    const auto source_rows=[](double metal_p,double ligand_p) {
        return std::vector<cov::NboNaoContribution>{
            row(1,"Val(3d)",3,2,0),row(1,"Ryd(4s)",4,0,0),
            row(1,"Ryd(4p)",4,1,metal_p),row(1,"Ryd(4d)",4,2,0),
            row(1,"Ryd(5p)",5,1,0),row(2,"Val(3p)",3,1,ligand_p)};
    };
    auto routed=route(wf,{source_rows(.3,.7),source_rows(.7,.3)});
    const auto shells=cov::mo_current_radial_shells(wf,routed,{0});
    require(shells.size()==3,"Current Cr source space must be 3d,4s,4p only");
    require(std::any_of(shells.begin(),shells.end(),[](const auto& s){return s.n==4&&s.l==1;}),
        "Required empty np may be Ryd-labelled");
    require(std::none_of(shells.begin(),shells.end(),[](const auto& s){return (s.n==4&&s.l==2)||s.n==5;}),
        "Higher radial d/p must not expand the current space");
    const auto ledger=cov::mo_group_composition_ledger(wf,&routed,{0,1},{0},shells);
    require(ledger.complete && ledger.member_count==2,"Complete two-member ledger unavailable");
    close(ledger.centre_current_p,.5,"Current p must use complete NAO norm");
    close(ledger.ligand_valence_p,.5,"Ligand Val p must use complete norm");
    close(ledger.weight_sum,1,"Complete norm must close");
    close(ledger.centre_current_s+ledger.centre_current_p+ledger.centre_current_d+
        ledger.centre_current_f+ledger.centre_other+ledger.ligand_valence+
        ledger.ligand_other+ledger.other_atoms+ledger.core+ledger.unresolved,1,"Exclusive buckets do not close");

    // Chemical applicability is independent of composition percentages. A
    // disconnected ion is not a ligand, and ordinary organic centres are not metals.
    cov::Wavefunction scope_w;scope_w.atoms.resize(5);scope_w.orbitals.resize(1);
    scope_w.atoms[0].atomic_number=11;scope_w.atoms[1].atomic_number=8;
    scope_w.atoms[2].atomic_number=1;scope_w.atoms[3].atomic_number=17;
    scope_w.atoms[4].atomic_number=24;
    auto scope_r=route(scope_w,{{row(1,"Val(3s)",3,0,.3),row(2,"Val(2p)",2,1,.4),
                              row(4,"Val(3p)",3,1,.3)}});
    scope_r.interaction_graph.status=cov::RoutedStatus::Available;
    scope_r.interaction_graph.value=cov::InteractionGraph{};
    for(const auto pair:std::vector<std::pair<std::size_t,std::size_t>>{{0,1},{1,2},{1,4}}) {
        cov::InteractionEdge e;e.atom_a=pair.first;e.atom_b=pair.second;
        e.kind=pair.first==0?cov::InteractionKind::CoordinationContact:cov::InteractionKind::CovalentConnectivity;
        e.strength=cov::InteractionStrength::StrongConnectivity;
        scope_r.interaction_graph.value->edges.push_back(e);
    }
    const auto scope=cov::mo_composition_scope(scope_w,&scope_r,{0});
    require(scope.applicable && scope.ligand_atoms==std::vector<std::size_t>{1,2},
        "General-metal composition must include connected ligand fragment, stop at other metals");
    require(scope.other_atoms==std::vector<std::size_t>{3,4},"Other ions/metals must not become ligand space");
    const auto scoped=cov::mo_group_composition_ledger(scope_w,&scope_r,{0},{0},{{0,3,0}});
    close(scoped.other_atoms,.3,"Disconnected counterion must retain its own denominator share");
    close(scoped.ligand_valence,.4,"Ligand cannot absorb counterion contribution");
    for(const int z:{1,6,7,8,15,16}) {
        scope_w.atoms[0].atomic_number=z;
        scope_r.canonical_fingerprint=cov::nbo_canonical_fingerprint(scope_w);
        require(!cov::mo_composition_scope(scope_w,&scope_r,{0}).applicable,
            "Nonmetal centre must never activate a metal-ligand composition panel");
    }

    // Rotate an orthonormal two-member NAO subspace. Individual weights change,
    // but the complete group trace and its shared denominator remain invariant.
    const double c=std::cos(.63),s=std::sin(.63),a=std::sqrt(.3),b=std::sqrt(.7);
    auto rotated=route(wf,{source_rows(std::pow(c*a-s*b,2),std::pow(c*b+s*a,2)),
                          source_rows(std::pow(s*a+c*b,2),std::pow(s*b-c*a,2))});
    const auto rotated_ledger=cov::mo_group_composition_ledger(wf,&rotated,{0,1},{0},shells);
    close(rotated_ledger.centre_current_p,ledger.centre_current_p,"Group rotation changes metal p trace");
    close(rotated_ledger.ligand_valence_p,ledger.ligand_valence_p,"Group rotation changes ligand p trace");

    routed.mo_composition[0].value->rows={row(1,"Val(3d)",3,2,.0014),
        row(1,"Ryd(4p)",4,1,.00047),row(2,"Val(3s)",3,0,.99313),row(1,"Cor(3p)",3,1,.005)};
    const auto deep=cov::mo_group_composition_ledger(wf,&routed,{0},{0},shells);
    require(deep.complete,"Deep ligand source must remain measurable");
    close(deep.ligand_valence_s,.99313,"Deep ligand s must not disappear from the denominator");
    close(deep.centre_current_p,.00047,"Tiny actual np must not be renormalized into dominance");
    routed.mo_composition[0].value->rows={row(1,"Ryd(4d)",4,2,.074),
        row(1,"Val(3d)",3,2,.0001),row(2,"Ryd(4p)",4,1,.9259)};
    const auto high=cov::mo_group_composition_ledger(wf,&routed,{0},{0},shells);
    close(high.centre_other,.074,"All non-core metal must not become current valence");
    close(high.centre_current_d,.0001,"Current d cannot absorb higher radial d");

    routed.mo_composition[1].status=cov::RoutedStatus::Insufficient;
    const auto partial=cov::mo_group_composition_ledger(wf,&routed,{0,1},{0},shells);
    require(partial.available && !partial.complete && partial.unresolved>=.5,
        "A missing member must remain explicit, not masquerade as measured zero");
    const auto missing=cov::mo_group_composition_ledger(wf,nullptr,{0},{0},shells);
    require(!missing.available && !missing.complete && missing.unresolved==1,"Missing complete source must be unavailable");
    routed.canonical_fingerprint="different-source";
    require(!cov::mo_group_composition_ledger(wf,&routed,{0},{0},shells).available,
        "Foreign canonical identity must be rejected");
    require(!cov::analyse_mo_sigma_framework(wf,nullptr,&routed,{0}).available,
        "Missing sigma source must be unavailable");
    // Complete occupied LP and centre-ligand BD form a fixed sigma projector.
    // A low-centre ligand donor remains visible to this independent space.
    cov::Wavefunction sw;sw.atoms.resize(2);sw.atoms[0].atomic_number=24;
    sw.atoms[1].atomic_number=7;sw.atoms[1].z=3;sw.basis_count=5;
    sw.ao_overlap.assign(25,0);sw.orbitals.resize(5);
    for(std::size_t i=0;i<5;++i) {
        sw.ao_overlap[i*5+i]=1;sw.orbitals[i].coefficients.assign(5,0);
        sw.orbitals[i].coefficients[i]=1;sw.orbitals[i].occupation=2;
    }
    cov::Shell ligand_s;ligand_s.atom_index=1;ligand_s.basis_offset=0;ligand_s.angular_momentum=0;
    cov::Shell ligand_p=ligand_s;ligand_p.basis_offset=1;ligand_p.angular_momentum=1;
    cov::Shell metal_s=ligand_s;metal_s.atom_index=0;metal_s.basis_offset=4;
    sw.shells={ligand_s,ligand_p,metal_s};
    cov::NboIntegration si;si.canonical_fingerprint=cov::nbo_canonical_fingerprint(sw);
    si.dataset.association.compatible=true;
    cov::NboMatrix nm;nm.kind="NAONBO";nm.rows=nm.columns=5;nm.values.assign(25,0);
    for(std::size_t i=0;i<3;++i)nm.values[i*5+i]=1;
    nm.values[3*5+3]=.8;nm.values[4*5+3]=.6;
    nm.values[3*5+4]=-.6;nm.values[4*5+4]=.8;si.dataset.matrices.push_back(nm);
    for(std::size_t i=0;i<5;++i) {
        cov::NboNao nao;nao.id=i+1;nao.atom=i<4?2:1;
        nao.type="Val";nao.angular=i==0||i==4?"s":"p";si.dataset.naos.push_back(nao);
    }
    for(const auto id:{std::size_t{0},std::size_t{3}}) {
        cov::NboOrbital literal;literal.id=id+1;literal.kind=id==0?"LP":"BD";
        literal.occupation=2;si.dataset.orbitals.push_back(literal);
        cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NBO,cov::NboSpin::Total,id};
        descriptor.orthonormal_basis=true;descriptor.atoms=id==0?std::vector<std::size_t>{1}:
            std::vector<std::size_t>{0,1};descriptor.coefficients.resize(5);
        for(std::size_t i=0;i<5;++i) {
            descriptor.coefficients[i]=nm.values[i*5+id];cov::NboMoLink link;
            link.orbital=descriptor.ref;link.canonical_index=i;link.coefficient=descriptor.coefficients[i];
            link.weight=link.coefficient*link.coefficient;si.links.push_back(link);
        }
        si.orbitals.push_back(descriptor);
    }
    cov::RoutedAnalysis sr;sr.canonical_fingerprint=si.canonical_fingerprint;
    sr.interaction_graph.status=cov::RoutedStatus::Available;
    sr.interaction_graph.value=cov::InteractionGraph{};
    cov::InteractionEdge edge;edge.atom_a=0;edge.atom_b=1;
    edge.strength=cov::InteractionStrength::StrongConnectivity;
    sr.interaction_graph.value->edges.push_back(edge);
    const auto sigma=cov::analyse_mo_sigma_framework(sw,&si,&sr,{0});
    require(sigma.available && sigma.rank==2,"Complete LP and centre-ligand BD sigma source unavailable");
    close(sigma.canonical_weights[0],1,"Ligand LP cannot require a centre population floor");
    close(sigma.canonical_weights[3],.64,"Sigma projection must retain original canonical norm");
    close(sigma.occupied_trace,2,"Complete donor occupied trace does not close");
    // The complete spin source must project onto a declared shared spatial
    // canonical basis even when the legacy physical-spin link list is empty.
    sw.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;
    const auto first_sources=si.orbitals;
    const auto first_literals=si.dataset.orbitals;
    auto beta_matrix=nm;beta_matrix.spin=cov::NboSpin::Beta;si.dataset.matrices.push_back(beta_matrix);
    for(auto descriptor:first_sources) {
        descriptor.ref.spin=cov::NboSpin::Beta;si.orbitals.push_back(descriptor);
    }
    for(auto literal:first_literals){literal.spin=cov::NboSpin::Beta;literal.occupation=1;si.dataset.orbitals.push_back(literal);}
    const auto first_naos=si.dataset.naos;
    for(auto nao:first_naos){nao.spin=cov::NboSpin::Beta;si.dataset.naos.push_back(nao);}
    // Both explicit spin blocks, with Total aliases excluded.
    for(auto& d:si.orbitals)if(d.ref.spin==cov::NboSpin::Total)d.ref.spin=cov::NboSpin::Alpha;
    for(auto& o:si.dataset.orbitals)if(o.spin==cov::NboSpin::Total)o.spin=cov::NboSpin::Alpha;
    for(auto& nao:si.dataset.naos)if(nao.spin==cov::NboSpin::Total)nao.spin=cov::NboSpin::Alpha;
    si.dataset.matrices[0].spin=cov::NboSpin::Alpha;si.links.clear();
    sr.canonical_fingerprint=si.canonical_fingerprint=cov::nbo_canonical_fingerprint(sw);
    const auto shared_sigma=cov::analyse_mo_sigma_framework(sw,&si,&sr,{0});
    require(shared_sigma.available && shared_sigma.rank==2 && shared_sigma.sources.size()==4,
        "Verified shared spatial canonical basis must support both source spin blocks");
    close(shared_sigma.canonical_weights[3],.64,"Shared alpha/beta projections must not double spatial norm");
    sw.orbitals.pop_back();
    sr.canonical_fingerprint=si.canonical_fingerprint=cov::nbo_canonical_fingerprint(sw);
    require(!cov::analyse_mo_sigma_framework(sw,&si,&sr,{0}).available,
        "An incomplete sigma source column must not become a low contribution");
    std::cout<<"complete-NAO closure, radial scope, degenerate rotation, deep s, high radial and missing-source checks passed\n";
    return 0;
} catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}}
