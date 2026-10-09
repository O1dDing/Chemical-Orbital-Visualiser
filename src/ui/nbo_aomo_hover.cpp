#include "cov/nbo_aomo_hover.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/nbo_aomo_text.hpp"
#include "cov/orbital_chemistry_summary.hpp"
#include "cov/orbital_ui_text.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <regex>
#include <utility>
#include <sstream>
#include <iomanip>

namespace cov::ui {
namespace {
enum class Word : std::size_t {
    Unknown, Atomic, AtomicCombination, Direction, Spherical, Symmetry,
    Subspace, Mainly, MixedAngular, MixedRole, PairAxes, OrbitalAnalysis,
    PairOverlap, Bonding, Antibonding, Nonbonding, InPhase, OutPhase,
    MixedPhase, PhaseUnknown, PhaseCaution, Count
};
using Row=std::array<const char*,4>;
constexpr std::array<Row,static_cast<std::size_t>(Word::Count)> words{{
    Row{"Character unresolved", "性质未确定", "性質未確定", "Caractère indéterminé"},
    Row{"Atomic basis function", "原子基函数", "原子基底関数", "Fonction de base atomique"},
    Row{"Atomic NAO combination", "原子 NAO 组合", "原子 NAO の組合せ", "Combinaison de NAO atomiques"},
    Row{"Direction: calculation frame", "方向：计算坐标系", "方向：計算の座標系", "Direction : repère du calcul"},
    Row{"Spherical s component", "球对称 s 分量", "球対称の s 成分", "Composante s sphérique"},
    Row{"Symmetry: ", "对称类型：", "対称性：", "Symétrie : "},
    Row{"Orbital subspace; expand to inspect members", "轨道子空间；展开查看成员", "軌道部分空間：展開して各軌道を確認", "Sous-espace orbital ; développer pour voir les membres"},
    Row{"Mainly ", "主要为 ", "主に ", "Principalement "},
    Row{"Mixed angular character", "混合角向性质", "混合した角度特性", "Caractère angulaire mixte"},
    Row{"Mixed bonding character", "混合成键性质", "混合した結合性", "Caractère de liaison mixte"},
    Row{" (atom-pair axes)", "（原子对轴）", "（原子対の軸）", " (axes des paires atomiques)"},
    Row{" (orbital analysis)", "（轨道分析）", "（軌道解析）", " (analyse orbitale)"},
    Row{" (atom-pair overlap)", "（原子对重叠分析）", "（原子対の重なり解析）", " (recouvrement des paires atomiques)"},
    Row{"bonding", "成键", "結合性", "liant"},
    Row{"antibonding", "反键", "反結合性", "antiliant"},
    Row{"nonbonding", "非键", "非結合性", "non liant"},
    Row{"In-phase basis combination", "基函数同相组合", "基底関数の同位相の組合せ", "Combinaison de base en phase"},
    Row{"Out-of-phase basis combination", "基函数反相组合", "基底関数の逆位相の組合せ", "Combinaison de base en opposition de phase"},
    Row{"Mixed basis phases", "基函数相位混合", "基底関数の位相が混在", "Phases de base mixtes"},
    Row{"Basis phase unresolved", "基函数相位未确定", "基底関数の位相は未確定", "Phase de base indéterminée"},
    Row{"NAO basis phases alone do not determine bonding", "仅凭 NAO 基底相位不能判断成键", "NAO 基底の位相だけでは結合性を判断できません", "Les phases des NAO seules ne déterminent pas le caractère liant"}
}};
std::size_t language_index(Language language) {
    const auto i=static_cast<std::size_t>(language);
    return i<4?i:0;
}
std::string text(Word word,Language language) {
    return words[static_cast<std::size_t>(word)][language_index(language)];
}
std::string atoms_text(const Wavefunction& wf,const std::vector<std::size_t>& atoms) {
    std::string result;
    for(std::size_t i=0;i<atoms.size();++i) {
        if(i==4) {result+=" +"+std::to_string(atoms.size()-i);break;}
        if(i)result+=", ";
        const auto a=atoms[i];
        result+=(a<wf.atoms.size()?wf.atoms[a].symbol:"?")+std::to_string(a+1);
    }
    return result;
}
std::string spin_text(NboSpin spin) {
    return spin==NboSpin::Total?"":spin==NboSpin::Alpha?" [α]":" [β]";
}
bool has_spin(const std::string& title) {
    return title.find("[alpha]")!=std::string::npos || title.find("[beta]")!=std::string::npos ||
        title.find("[α]")!=std::string::npos || title.find("[β]")!=std::string::npos;
}
std::string concise_spin(std::string title) {
    for(const auto& replacement:{std::pair{"[alpha]","[α]"},std::pair{"[beta]","[β]"}}) {
        std::size_t at=0;
        while((at=title.find(replacement.first,at))!=std::string::npos) {
            title.replace(at,std::char_traits<char>::length(replacement.first),replacement.second);
            at+=std::char_traits<char>::length(replacement.second);
        }
    }
    return title;
}
std::string shell_direction(const std::string& type,const std::string& angular) {
    // Only the explicitly printed principal shell is usable. Gaussian AO
    // contraction order and an unnumbered p direction cannot supply n.
    static const std::regex shell(R"(\(\s*(\d+)([spdfghik])\s*\))");
    std::smatch match;
    if(std::regex_search(type,match,shell)) {
        const auto n=match[1].str(),family=match[2].str();
        if(angular.empty())return n+family;
        if(angular.front()==family.front())return n+angular;
        if(angular.rfind(n+family,0)==0)return angular;
    }
    return type+(type.empty() || angular.empty()?"":" ")+angular;
}
const NboNao* nao_row(const NboIntegration& data,const NboOrbitalRef& ref) {
    if(ref.kind!=NboOrbitalKind::NAO && ref.kind!=NboOrbitalKind::PNAO)return nullptr;
    for(const auto& row:data.dataset.naos)
        if(row.id && row.id-1==ref.index && row.spin==ref.spin)return &row;
    return nullptr;
}
std::string basis_title(const Wavefunction& wf,const NboNao& row,NboSpin spin) {
    const auto atom=row.atom?atoms_text(wf,{row.atom-1}):std::string("?");
    return atom+" "+shell_direction(row.type,row.angular)+spin_text(spin);
}
void symmetry_line(std::vector<std::string>& lines,const NboAomoNode& node,Language language) {
    if(node.symmetry_name_verified && !node.symmetry_irrep.empty())
        lines.push_back(text(Word::Symmetry,language)+node.symmetry_irrep);
}
const char* angular_symbol(OrbitalAngularFamily family) {
    switch(family) {
    case OrbitalAngularFamily::Sigma:return "σ";
    case OrbitalAngularFamily::Pi:return "π";
    case OrbitalAngularFamily::Delta:return "δ";
    case OrbitalAngularFamily::Phi:return "φ";
    default:return nullptr;
    }
}
std::string role_text(OrbitalBondingRole role,Language language) {
    switch(role) {
    case OrbitalBondingRole::Bonding:return text(Word::Bonding,language);
    case OrbitalBondingRole::Antibonding:return text(Word::Antibonding,language);
    case OrbitalBondingRole::Nonbonding:return text(Word::Nonbonding,language);
    default:return {};
    }
}
Word phase_word(const NboSalcOrbital& orbital) {
    // Signs are relative to the actual stored basis, not spatial lobe phases.
    // Multiplying the complete SALC by -1 must not change this description.
    bool positive=false,negative=false;
    std::size_t nonzero=0;
    for(const auto& term:orbital.terms) {
        if(!std::isfinite(term.coefficient))return Word::PhaseUnknown;
        if(term.coefficient==0)continue;
        ++nonzero;
        positive=positive || term.coefficient>0;
        negative=negative || term.coefficient<0;
    }
    if(nonzero<2)return Word::PhaseUnknown;
    if(!(positive&&negative))return Word::InPhase;
    return nonzero==2?Word::OutPhase:Word::MixedPhase;
}
} // namespace

std::vector<std::string> nbo_aomo_hover_lines(const NboAomoNode& node,
    const Wavefunction& wf,const NboIntegration& data,const NboSalcModel* model,
    Language language) {
    std::vector<std::string> lines;
    const auto title=concise_spin(node.individual_label.empty()?node.label:node.individual_label);
    if(node.weak_display_container) {
        lines={title,aomo_text(language,"Click to expand real members")};
        for(auto member:node.member_canonical_indices)if(member<wf.orbitals.size()) {
            std::ostringstream value;value<<std::setprecision(9)<<canonical_mo_source_label(wf,member)
                <<": E="<<wf.orbitals[member].energy_hartree<<" Ha; n="<<wf.orbitals[member].occupation;
            lines.push_back(value.str());
        }
        return lines;
    }
    if(node.group_header) return {title,text(Word::Subspace,language)};
    if(node.spatial_spin){
        const auto value=[&](const std::optional<double>& v){if(!v)return std::string(aomo_text(language,"Unavailable"));
            std::ostringstream out;out<<std::setprecision(9)<<*v;return out.str();};
        lines={title};
        for(const auto& channel:node.spatial_spin->channels){
            lines.push_back(std::string(nbo_spin_name(channel.spin))+": E="+value(channel.energy_hartree)+" Ha; n="+value(channel.occupation));
            // Rotated source diagonals differ from the common-basis expectation.
            if(channel.members.size()==1){const auto& source=channel.members.front();
                lines.push_back(source.id+": E="+value(source.energy_hartree)+" Ha; n="+value(source.occupation));}
        }
        lines.push_back("E=(Eα+Eβ)/2: "+value(node.spatial_spin->energy_hartree)+" Ha; n=nα+nβ: "+value(node.spatial_spin->occupation));
        if(!node.spatial_spin->energy_hartree)lines.push_back(node.spatial_spin->energy_status);
        return lines;
    }
    if(node.canonical_index) {
        const auto index=*node.canonical_index;
        if(index>=wf.orbitals.size())return {title,text(Word::Unknown,language)};
        const auto& mo=wf.orbitals[index];
        NboAomoName display_name;display_name.label=node.individual_label.empty()?node.label:node.individual_label;
        display_name.irrep=node.symmetry_irrep;display_name.ordinal=node.symmetry_ordinal;
        display_name.verified=node.symmetry_name_verified;
        std::string identity=node.individual_label.empty()?
            concise_spin(canonical_mo_display_label(wf,index,&display_name)):title;
        const bool spin_resolved=wf.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin ||
            std::any_of(wf.orbitals.begin(),wf.orbitals.end(),[](const auto& o){return o.spin==Spin::Beta;});
        if(spin_resolved && !has_spin(identity))
            identity+=spin_text(mo.spin==Spin::Beta?NboSpin::Beta:NboSpin::Alpha);
        lines.push_back(identity);
        if(display_name.verified&&!display_name.ordinal)
            lines.push_back(text(Word::Symmetry,language)+orbital_irrep_display_label(display_name));
        const auto& chemistry=mo.chemistry;
        if(chemistry.available) {
            std::string pair_identity;
            if(chemistry.interactions.size()==1) {
                const auto& pair=chemistry.interactions.front();
                pair_identity=" ["+atoms_text(wf,{pair.atom_a,pair.atom_b})+"]";
            }
            const auto context=text(chemistry.interactions.empty()?Word::OrbitalAnalysis:Word::PairAxes,language)+pair_identity;
            if(chemistry.channel.status==ChemistryStatus::Determined) {
                if(const auto* symbol=angular_symbol(chemistry.channel.dominant))
                    lines.push_back(text(Word::Mainly,language)+symbol+context);
            } else if(chemistry.channel.status==ChemistryStatus::Percentages)
                lines.push_back(text(Word::MixedAngular,language)+
                    orbital_channel_fraction_summary(chemistry.channel)+context);
            // Existing bonding evidence is based on per-orbital atom-pair
            // overlap, not on a whole-bond order or an occupation count.
            if(!node.bonding_scope_status.empty()) {
                if(node.bonding_class!=BondingClass::Unclassified)
                    lines.push_back(std::string(aomo_text(language,"Skeleton group bonding character"))+": "+
                        localised_bonding_class(node.bonding_class,language));
            } else if(chemistry.bonding.status==ChemistryStatus::Determined) {
                const auto role=role_text(chemistry.bonding.dominant,language);
                if(!role.empty())lines.push_back(text(Word::Mainly,language)+role+text(Word::PairOverlap,language)+pair_identity);
            } else if(chemistry.bonding.status==ChemistryStatus::Percentages)
                lines.push_back(text(Word::MixedRole,language)+
                    orbital_bonding_fraction_summary(chemistry.bonding,
                        role_text(OrbitalBondingRole::Bonding,language).c_str(),
                        role_text(OrbitalBondingRole::Antibonding,language).c_str(),
                        role_text(OrbitalBondingRole::Nonbonding,language).c_str())+
                    text(Word::PairOverlap,language)+pair_identity);
        }
        if(lines.size()==1)lines.push_back(text(Word::Unknown,language));
        symmetry_line(lines,node,language);
        return lines;
    }
    if(node.salc_index) {
        if(!model || !model->available || model->dataset_id!=data.id ||
           *node.salc_index>=model->orbitals.size())
            return {title,text(Word::Unknown,language)};
        const auto& orbital=model->orbitals[*node.salc_index];
        if(orbital.atoms.size()==1) {
            const auto* row=orbital.terms.size()==1?nao_row(data,orbital.terms.front().orbital):nullptr;
            lines.push_back(row?basis_title(wf,*row,orbital.spin):
                atoms_text(wf,orbital.atoms)+" "+shell_direction(orbital.type,orbital.angular)+spin_text(orbital.spin));
            lines.push_back(text(orbital.terms.size()==1?Word::Atomic:Word::AtomicCombination,language)+
                (orbital.terms.size()==1?" (NAO)":""));
            lines.push_back(text(orbital.angular=="s"?Word::Spherical:Word::Direction,language));
        } else {
            const auto spin=has_spin(title)?std::string{}:spin_text(orbital.spin);
            lines.push_back(title+" · SALC "+atoms_text(wf,orbital.atoms)+" "+shell_direction(orbital.type,orbital.angular)+spin);
            lines.push_back(text(phase_word(orbital),language));
            lines.push_back(text(Word::PhaseCaution,language));
            symmetry_line(lines,node,language);
        }
        return lines;
    }
    if(node.orbital) {
        const auto& ref=*node.orbital;
        const auto* row=nao_row(data,ref);
        const auto* descriptor=nbo_orbital(data,ref);
        lines.push_back(row?basis_title(wf,*row,ref.spin):
            descriptor?concise_spin(descriptor->label):title);
        if(ref.kind==NboOrbitalKind::GaussianAO || ref.kind==NboOrbitalKind::NAO || ref.kind==NboOrbitalKind::PNAO) {
            lines.push_back(text(Word::Atomic,language)+" ("+nbo_orbital_kind_name(ref.kind)+")");
            lines.push_back(text(row && row->angular=="s"?Word::Spherical:Word::Direction,language));
        } else lines.push_back(text(Word::Unknown,language));
        return lines;
    }
    return {title,text(Word::Unknown,language)};
}

std::string nbo_aomo_hover_glyph_seed(Language language) {
    // Common scientific names use Unicode subscripts and spectroscopic signs;
    // source-number fallback and known-irrep auxiliary lines share this font.
    std::string result="MO NAO SALC AO PNAO α β σ π δ φ γ Σ Π Δ Φ Γ ∞ · + − ? [] () 0123456789 ₀₁₂₃₄₅₆₇₈₉ ′ ″";
    for(const auto& row:words){result+=' ';result+=row[language_index(language)];}
    return result;
}
} // namespace cov::ui
