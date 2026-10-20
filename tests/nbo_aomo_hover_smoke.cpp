#include "cov/nbo_aomo_hover.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::string joined(const std::vector<std::string>& lines) {
    require(lines.size()>=2 && lines.size()<=4,"hover must have two to four lines");
    std::string result;
    for(const auto& line:lines){require(!line.empty(),"empty hover line");result+=line+'\n';}
    require(result.find("Ha")==std::string::npos,"energy leaked into concise hover");
    require(result.find("occupation")==std::string::npos,"occupation leaked into concise hover");
    return result;
}
void require_glyphs(const std::string& value,const std::string& seed) {
    for(std::size_t i=0;i<value.size();) {
        const auto c=static_cast<unsigned char>(value[i]);
        const std::size_t count=c<128?1:c<224?2:c<240?3:4;
        if(c>=128)require(seed.find(value.substr(i,count))!=std::string::npos,"localized hover glyph missing from font seed");
        i+=count;
    }
}
}
int main(){try {
    using namespace cov;
    using namespace cov::ui;
    Wavefunction wf;wf.atoms={{"O",8},{"H",1},{"H",1}};
    wf.orbitals.resize(1);
    NboIntegration data;data.id="fixture";
    NboAomoNode mo;mo.label="1b2";mo.canonical_index=0;
    mo.symmetry_irrep="b2";mo.symmetry_name_verified=true;
    mo.symmetry_ordinal=1;
    auto render=[&](const NboAomoNode& n,const NboSalcModel* s=nullptr,Language l=Language::English){return joined(nbo_aomo_hover_lines(n,wf,data,s,l));};
    require(render(mo).find("unresolved")!=std::string::npos,"missing chemistry must stay unresolved");
    auto& chemistry=wf.orbitals[0].chemistry;
    chemistry.available=true;
    chemistry.channel.status=ChemistryStatus::Determined;
    chemistry.channel.dominant=OrbitalAngularFamily::Pi;
    chemistry.bonding.status=ChemistryStatus::Determined;
    chemistry.bonding.dominant=OrbitalBondingRole::Antibonding;
    OrbitalPairInteraction pair;pair.atom_a=0;pair.atom_b=1;pair.total_mayer_index=2.0;
    chemistry.interactions.push_back(pair);wf.orbitals[0].occupation=2;
    const auto definite=render(mo);
    require(definite.find("\xCF\x80")!=std::string::npos,"existing pi evidence lost");
    require(definite.find("antibonding")!=std::string::npos,"existing role evidence lost");
    require(definite.find("atom-pair axes")!=std::string::npos && definite.find("atom-pair overlap")!=std::string::npos,"pair context lost");
    require(definite.find("double")==std::string::npos && definite.find("single")==std::string::npos,"one orbital or Mayer index became a whole bond order");
    chemistry.channel.status=ChemistryStatus::Percentages;
    chemistry.channel.sigma=.2;chemistry.channel.pi=.8;chemistry.channel.undetermined=0;
    chemistry.bonding.status=ChemistryStatus::Percentages;
    chemistry.bonding.bonding=.3;chemistry.bonding.antibonding=.7;
    chemistry.bonding.undetermined=0;
    const auto mixed=render(mo);
    require(mixed.find("Mixed angular")!=std::string::npos && mixed.find("Mixed bonding")!=std::string::npos,"percentage evidence became a definite assignment");
    require(mixed.find("σ 20% · π 80%")!=std::string::npos &&
            mixed.find("bonding 30% · antibonding 70%")!=std::string::npos,
            "hover lost the orbital percentages displayed in the detail panel");
    require(mixed.find("Mainly")==std::string::npos,"mixed evidence became a pure channel/role");
    wf.orbitals[0].spin=Spin::Beta;
    require(render(mo).find("[\xCE\xB2]")!=std::string::npos,"spin identity lost");
    mo.label="1b2 [beta]";
    const auto spin_label=render(mo);
    require(spin_label.find("beta")==std::string::npos && spin_label.find("[\xCE\xB2]")==spin_label.rfind("[\xCE\xB2]"),"spin suffix not localized or duplicated");
    mo.label="1b2 [\xCE\xB2]";
    require(render(mo)==spin_label,"symbolic spin suffix duplicated");
    chemistry.available=false;
    require(render(mo).find("Mixed")==std::string::npos,"unavailable chemistry exposed stale classification");

    NboNao row;row.id=1;row.atom=1;row.symbol="O";row.type="Val(2p)";row.angular="px";
    data.dataset.naos.push_back(row);
    NboAomoNode ao;ao.label="O1 2p";ao.orbital=NboOrbitalRef{NboOrbitalKind::NAO,NboSpin::Total,0};
    const auto atomic=render(ao);
    require(atomic.find("O1 2px")!=std::string::npos && atomic.find("Val(")==std::string::npos,"explicit shell/direction not compacted");
    require(atomic.find("Atomic basis")!=std::string::npos && atomic.find("bonding")==std::string::npos && atomic.find("\xCF\x83")==std::string::npos,"atomic p direction became sigma bonding");
    data.dataset.naos[0].type="Ryd";
    require(render(ao).find("Ryd px")!=std::string::npos && render(ao).find("2px")==std::string::npos,"principal shell fabricated from unnumbered direction");
    data.dataset.naos[0].type="Ryd(10p)";
    require(render(ao).find("10px")!=std::string::npos,"multi-digit producer principal shell lost");
    data.dataset.naos[0].type="Val(2p)";

    NboSalcModel model;model.available=true;model.dataset_id=data.id;
    NboSalcOrbital salc;salc.atoms={1,2};salc.type="Val(1s)";salc.angular="s";
    salc.terms={{{NboOrbitalKind::NAO,NboSpin::Total,1},.7071067811865475},{{NboOrbitalKind::NAO,NboSpin::Total,2},.7071067811865475}};
    model.orbitals.push_back(salc);
    NboAomoNode sn;sn.label="1a1";sn.salc_index=0;sn.symmetry_irrep="a1";sn.symmetry_name_verified=true;
    const auto phase=render(sn,&model);
    require(phase.find("H2, H3 1s")!=std::string::npos,"SALC atom/shell identity missing");
    require(phase.find("In-phase basis")!=std::string::npos &&
        phase.find("NAO basis phases alone do not determine bonding")!=std::string::npos,
        "basis phase confused with bonding");
    for(auto& term:model.orbitals[0].terms)term.coefficient=-term.coefficient;
    require(render(sn,&model)==phase,"global SALC phase changed its chemical description");
    model.orbitals[0].terms[1].coefficient*=-1;
    require(render(sn,&model).find("Out-of-phase basis")!=std::string::npos,"opposite basis signs lost");
    model.orbitals[0].terms.push_back({{NboOrbitalKind::NAO,NboSpin::Total,3},.2});
    require(render(sn,&model).find("Mixed basis phases")!=std::string::npos,"multicentre sign pattern overclassified");
    sn.symmetry_name_verified=false;
    require(render(sn,&model).find("Verified symmetry")==std::string::npos,"unverified symmetry presented as verified");
    model.dataset_id="another attachment";
    require(render(sn,&model).find("unresolved")!=std::string::npos,"stale SALC attachment was accepted");
    model.dataset_id=data.id;
    model.orbitals[0].atoms={0};model.orbitals[0].type="Val(2p)";model.orbitals[0].angular="px";
    model.orbitals[0].terms={{{NboOrbitalKind::NAO,NboSpin::Total,0},1}};
    require(render(sn,&model).find("Atomic basis")!=std::string::npos && render(sn,&model).find("phase")==std::string::npos,"single-atom NAO inherited molecular SALC classification");
    auto unordered=mo;unordered.symmetry_ordinal=0;
    auto unknown=unordered;unknown.symmetry_name_verified=false;unknown.symmetry_irrep.clear();
    const auto source_index=wf.orbitals[0].source_orbital_index;
    wf.orbitals[0].source_orbital_index=7;
    for(auto language:{Language::English,Language::ChineseSimplified,Language::Japanese,Language::French}) {
        const auto localized=render(ao,nullptr,language);
        require((language==Language::English||localized!=atomic) && localized.find("px")!=std::string::npos,"localized atomic direction identity lost");
        const auto seed=nbo_aomo_hover_glyph_seed(language);
        require_glyphs(localized,seed);
        require_glyphs(render(mo,nullptr,language),seed);
        require_glyphs(render(sn,&model,language),seed);
        const auto without_order=render(unordered,nullptr,language);
        require(without_order.find("MO 8")!=std::string::npos &&
                without_order.find("b₂")!=std::string::npos && without_order.find("1b₂")==std::string::npos,
                "known symmetry without order must remain auxiliary to the true source identity");
        require_glyphs(without_order,seed);
        const auto without_symmetry=render(unknown,nullptr,language);
        require(without_symmetry.find("MO 8")!=std::string::npos && without_symmetry.find("b₂")==std::string::npos,
                "missing symmetry must not inherit the old node title's scientific name");
        require_glyphs(without_symmetry,seed);
    }
    wf.orbitals[0].source_orbital_index=source_index;
    NboAomoNode group;group.label="Set (3)";group.group_header=true;group.canonical_index=0;
    require(render(group).find("subspace")!=std::string::npos && render(group).find("[\xCE\xB2]")==std::string::npos,"group header misrepresented as an individual orbital");
    require(wf.orbitals[0].occupation==2 && data.dataset.naos[0].angular=="px","hover changed source data");
    std::cout<<"nbo_aomo_hover_smoke passed\n";return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
