#include "cov/mo_diagram.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/ligand_field.hpp"
#include "cov/local_geometry.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/pi_pair_evidence.hpp"
#include "cov/point_group_catalog.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace cov {

namespace {
double diagram_orbital_energy(const Wavefunction& w,const MODiagramOptions& options,std::size_t index) {
    return options.canonical_display_energies.size()==w.orbitals.size()
        ?options.canonical_display_energies[index]:w.orbitals[index].energy_hartree;
}
bool composition_metal(int z) noexcept {
    if((z>=21 && z<=30)||(z>=39 && z<=48)||(z>=57 && z<=80)||
       (z>=89 && z<=112))return true;
    switch(z) {
        case 3:case 4:case 11:case 12:case 13:case 19:case 20:
        case 31:case 37:case 38:case 49:case 50:case 55:case 56:
        case 81:case 82:case 83:case 87:case 88:case 113:case 114:
        case 115:case 116:return true;
        default:return false;
    }
}
}

MOCompositionScope mo_composition_scope(const Wavefunction& w,const RoutedAnalysis* routed,
                                        const std::vector<std::size_t>& centres) {
    MOCompositionScope result;
    if(centres.empty()){result.detail="no-explicit-composition-centre";return result;}
    for(auto a:centres) {
        if(a>=w.atoms.size() || !composition_metal(w.atoms[a].atomic_number)) {
            result.detail="composition-centre-is-not-a-metal";return result;
        }
        if(std::find(result.centre_atoms.begin(),result.centre_atoms.end(),a)==result.centre_atoms.end())
            result.centre_atoms.push_back(a);
    }
    const bool linked=routed && routed->canonical_fingerprint==nbo_canonical_fingerprint(w) &&
        routed->interaction_graph.available();
    const auto graph=linked?make_fixed_bonding_scope(w,*routed->interaction_graph.value):
        make_fixed_bonding_scope(w);
    std::vector<std::vector<std::size_t>> adjacency(w.atoms.size());
    for(const auto& [a,b]:graph.edges)if(a<w.atoms.size() && b<w.atoms.size()) {
        adjacency[a].push_back(b);adjacency[b].push_back(a);
    }
    std::vector<bool> visited(w.atoms.size(),false);
    std::queue<std::size_t> pending;
    for(auto a:result.centre_atoms){visited[a]=true;pending.push(a);}
    // Walk only verified skeletal/coordination connectivity. Other metals are
    // boundaries, not ligand atoms; connected ligand fragments retain all atoms.
    while(!pending.empty()) {
        const auto a=pending.front();pending.pop();
        for(auto b:adjacency[a])if(!visited[b] && !composition_metal(w.atoms[b].atomic_number)) {
            visited[b]=true;pending.push(b);result.ligand_atoms.push_back(b);
        }
    }
    for(std::size_t a=0;a<w.atoms.size();++a)if(!visited[a])result.other_atoms.push_back(a);
    std::sort(result.ligand_atoms.begin(),result.ligand_atoms.end());
    result.applicable=!result.ligand_atoms.empty();
    result.detail=result.applicable?"actual-metal-and-connected-ligand-fragments":
        "no-supported-associated-ligand-fragment";
    return result;
}

std::vector<MOCurrentRadialShell> mo_current_radial_shells(
    const Wavefunction& wavefunction,const RoutedAnalysis& routed,
    const std::vector<std::size_t>& centres) {
    if(routed.canonical_fingerprint!=nbo_canonical_fingerprint(wavefunction))return {};
    std::map<std::pair<std::size_t,int>,int> valence;
    std::set<std::tuple<std::size_t,int,int>> present;
    for(const auto& entry:routed.mo_composition) {
        if(!entry.available() || entry.provider!=RoutedProvider::Nbo)continue;
        for(const auto& row:entry.value->rows) {
            if(!row.atom || !row.principal_n || !row.angular_l || *row.principal_n<=0 ||
               *row.angular_l<0 || *row.angular_l>3 ||
               std::find(centres.begin(),centres.end(),row.atom-1)==centres.end())continue;
            const auto atom=row.atom-1;
            present.emplace(atom,*row.principal_n,*row.angular_l);
            if(row.type.rfind("Val",0)!=0)continue;
            const auto key=std::pair{atom,*row.angular_l};
            auto [it,inserted]=valence.emplace(key,*row.principal_n);
            if(!inserted)it->second=std::min(it->second,*row.principal_n);
        }
    }
    std::vector<MOCurrentRadialShell> result;
    for(const auto atom:centres) {
        if(atom>=wavefunction.atoms.size())continue;
        int outer_n=0;
        for(int l=0;l<=3;++l)if(const auto it=valence.find({atom,l});it!=valence.end()) {
            result.push_back({atom,it->second,l,"producer-Val-current-radial-shell"});
            // The source valence period supplies necessary ns/np, including
            // empty Ryd-labelled shells. Never follow a retained MO to higher nd/nf.
            outer_n=std::max(outer_n,it->second+(l==2?1:l==3?2:0));
        }
        for(int l=0;l<=1;++l)if(outer_n>0 && !valence.contains({atom,l}) &&
            present.contains({atom,outer_n,l}))
            result.push_back({atom,outer_n,l,"source-valence-period-required-ns-np-framework"});
    }
    return result;
}

MOGroupCompositionLedger mo_group_composition_ledger(
    const Wavefunction& wavefunction,const RoutedAnalysis* routed,
    const std::vector<std::size_t>& members,const std::vector<std::size_t>& centres,
    const std::vector<MOCurrentRadialShell>& shells) {
    MOGroupCompositionLedger result;result.member_count=members.size();
    if(!routed || members.empty() ||
       routed->canonical_fingerprint!=nbo_canonical_fingerprint(wavefunction)) {
        result.detail="complete-NAO-source-unavailable-or-identity-mismatch";
        result.unresolved=1;return result;
    }
    result.source="complete-NAO-original-norm-member-average";
    const auto scope=mo_composition_scope(wavefunction,routed,centres);
    bool complete=true;
    const double divisor=static_cast<double>(members.size());
    for(const auto member:members) {
        if(member>=routed->mo_composition.size()) {
            complete=false;result.unresolved+=1/divisor;continue;
        }
        const auto& entry=routed->mo_composition[member];
        if(!entry.available() || entry.provider!=RoutedProvider::Nbo) {
            complete=false;result.unresolved+=1/divisor;continue;
        }
        result.available=true;complete=complete && entry.value->complete;
        double sum=0;
        for(const auto& row:entry.value->rows) {
            if(!std::isfinite(row.weight) || row.weight<0) {complete=false;continue;}
            sum+=row.weight;
            const double weight=row.weight/divisor;
            if(!row.atom || row.atom>wavefunction.atoms.size()) {
                result.unresolved+=weight;continue;
            }
            if(row.type.rfind("Cor",0)==0) {result.core+=weight;continue;}
            const auto atom=row.atom-1;
            const bool centre=std::find(centres.begin(),centres.end(),atom)!=centres.end();
            if(centre) {
                const bool current=row.principal_n && row.angular_l &&
                    std::any_of(shells.begin(),shells.end(),[&](const auto& shell) {
                        return shell.atom==atom && shell.n==*row.principal_n && shell.l==*row.angular_l;
                    });
                if(!current)result.centre_other+=weight;
                else switch(*row.angular_l) {
                    case 0:result.centre_current_s+=weight;break;
                    case 1:result.centre_current_p+=weight;break;
                    case 2:result.centre_current_d+=weight;break;
                    case 3:result.centre_current_f+=weight;break;
                    default:result.centre_other+=weight;
                }
            } else if(scope.applicable &&
                      std::find(scope.ligand_atoms.begin(),scope.ligand_atoms.end(),atom)==scope.ligand_atoms.end()) {
                result.other_atoms+=weight;
            } else if(row.type.rfind("Val",0)==0) {
                result.ligand_valence+=weight;
                if(row.angular_l==0)result.ligand_valence_s+=weight;
                if(row.angular_l==1)result.ligand_valence_p+=weight;
            } else result.ligand_other+=weight;
        }
        result.weight_sum+=sum/divisor;
        result.normalization_error=std::max(result.normalization_error,std::abs(sum-1));
        if(sum<1)result.unresolved+=(1-sum)/divisor;
        // Numerical closure validation, never a chemical display threshold.
        if(std::abs(sum-1)>1e-5)complete=false;
    }
    result.complete=result.available && complete;
    result.status=result.complete?"complete":result.available?"partial":"unavailable";
    result.detail=result.complete?"exclusive-complete-NAO-buckets; original-member-norm":
        "missing-or-incomplete-members; do-not-infer-low-coverage";
    return result;
}

OrbitalAnnotation annotate_orbital_legacy(const MolecularOrbital& orbital);
DiagramSelectionPlan build_valence_selection_plan_legacy(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata);
MODiagramData build_mo_diagram_data_legacy(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options);

namespace {

bool chemistry_available(const Wavefunction& wavefunction) {
    return std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const MolecularOrbital& orbital) {
            return orbital.chemistry.available;
        });
}

void attach_local_geometry_metadata(const Wavefunction& wavefunction,
                                    MODiagramData& data) {
    data.local_geometries.clear();
    for (const auto& environment:analyse_local_molecular_geometries(wavefunction)) {
        const auto* descriptor=coordination_geometry_descriptor(
            environment.geometry_id);
        if (descriptor==nullptr) continue;
        LocalGeometryDiagramDescriptor item;
        item.centre_atom=environment.centre_atom;
        item.geometry_id=std::string(descriptor->machine_id);
        item.geometry_name=std::string(descriptor->name);
        item.point_group=std::string(descriptor->point_group);
        item.neighbour_atoms=environment.neighbour_atoms;
        item.confidence=environment.confidence;
        item.angular_rms=environment.angular_rms;
        item.shape_measure=environment.shape_measure;
        item.radial_cv=environment.radial_cv;
        data.local_geometries.push_back(std::move(item));
    }
}

void attach_electronic_state_metadata(const Wavefunction& wavefunction,
                                      MODiagramData& data) {
    auto& state=data.electronic_state;
    state.source=wavefunction.source;
    state.charge=wavefunction.charge;
    state.multiplicity=wavefunction.multiplicity;
    state.alpha_electrons=wavefunction.alpha_electrons;
    state.beta_electrons=wavefunction.beta_electrons;
    state.charge_provenance=wavefunction.charge_provenance;
    state.multiplicity_provenance=wavefunction.multiplicity_provenance;
    state.electron_counts_provenance=wavefunction.electron_counts_provenance;
    state.scf_convergence=wavefunction.scf_convergence;
    state.scf_convergence_provenance=wavefunction.scf_convergence_provenance;
    state.stability=wavefunction.stability;
    state.stability_provenance=wavefunction.stability_provenance;
    state.stability_detail=wavefunction.stability_detail;
    state.spin_squared_before=wavefunction.spin_squared_before_annihilation;
    state.spin_squared_after=wavefunction.spin_squared_after_annihilation;
    state.spin_squared_provenance=wavefunction.spin_squared_provenance;
    state.atomic_partial_charges=wavefunction.atomic_partial_charges;
    state.atomic_partial_charge_scheme=wavefunction.atomic_partial_charge_scheme;
    state.atomic_partial_charge_provenance=
        wavefunction.atomic_partial_charge_provenance;
    state.point_group_detected=wavefunction.point_group_detected;
    state.point_group_used=wavefunction.point_group_used;
    state.point_group_provenance=wavefunction.point_group_provenance;
    state.point_group_detected_provenance=wavefunction.point_group_detected_provenance;
    state.point_group_used_provenance=wavefunction.point_group_used_provenance;
    state.source_title=wavefunction.source_title;
    state.source_route=wavefunction.source_route;
    state.enrichment_source=wavefunction.enrichment_source;
}

bool occupied(const MolecularOrbital& orbital,const double threshold) {
    return static_cast<double>(orbital.occupation)>threshold;
}

std::string family_name(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        default: return "unavailable";
    }
}

std::string script_number(const std::size_t value,const bool superscript) {
    static constexpr const char* superscript_digits[]={
        "⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
    static constexpr const char* subscript_digits[]={
        "₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};
    const std::string plain=std::to_string(value);
    std::string result;
    for (const char digit:plain) {
        const auto index=static_cast<std::size_t>(digit-'0');
        result+=superscript?superscript_digits[index]:subscript_digits[index];
    }
    return result;
}

std::string delocalised_pi_label(const std::size_t atoms,
                                 const double electrons) {
    const auto rounded=static_cast<std::size_t>(
        std::max(0LL,std::llround(electrons)));
    return std::string("Π")+script_number(atoms,true)+
           script_number(rounded,false);
}

OrbitalAnnotation chemistry_annotation(const MolecularOrbital& orbital) {
    const auto& chemistry=orbital.chemistry;
    OrbitalAnnotation result;
    if (!chemistry.available) return result;

    if (chemistry.channel.status==ChemistryStatus::Determined ||
        chemistry.channel.status==ChemistryStatus::Percentages) {
        result.family=family_name(chemistry.channel.dominant);
        result.family_source=AnnotationSource::Derived;
        result.family_confidence=chemistry.confidence;
    }

    if (chemistry.bonding.status==ChemistryStatus::Determined ||
        chemistry.bonding.status==ChemistryStatus::Percentages) {
        switch (chemistry.bonding.dominant) {
            case OrbitalBondingRole::Bonding:
                result.bonding_class=BondingClass::Bonding; break;
            case OrbitalBondingRole::Antibonding:
                result.bonding_class=BondingClass::Antibonding; break;
            case OrbitalBondingRole::Nonbonding:
                result.bonding_class=BondingClass::Nonbonding; break;
            default: break;
        }
        result.bonding_source=AnnotationSource::Derived;
        result.bonding_confidence=chemistry.confidence;
    }

    if (chemistry.multicentre_assignment_available && !chemistry.multicentre_label.empty()) {
        result.multicentre.available=true;
        result.multicentre.centres=
            chemistry.multicentre_participating_atoms;
        result.multicentre.electrons=
            chemistry.multicentre_participating_electrons;
        result.multicentre.atom_indices.assign(
            chemistry.multicentre_participating_atom_indices.begin(),
            chemistry.multicentre_participating_atom_indices.end());
        result.multicentre.label=chemistry.multicentre_label;
        result.multicentre.channel_count=chemistry.multicentre_channel_count;
        result.multicentre.source_subspace_id=
            chemistry.multicentre_source_subspace_id;
        result.multicentre.source_subspace_electron_count=
            chemistry.multicentre_source_electron_count;
        result.multicentre.source=AnnotationSource::Derived;
        result.multicentre.confidence=chemistry.multicentre_confidence;
        result.multicentre.heuristic=false;
    }

    if (std::isfinite(chemistry.delocalised_pi_weight) &&
        chemistry.delocalised_pi_weight>0.0 &&
        chemistry.delocalised_participating_atoms>2u &&
        !chemistry.delocalised_family_orbitals.empty() &&
        !chemistry.delocalised_family_id.empty()) {
        result.delocalised_pi.available=true;
        result.delocalised_pi.participating_atoms=
            chemistry.delocalised_participating_atoms;
        result.delocalised_pi.participating_electrons=
            chemistry.delocalised_participating_electrons;
        result.delocalised_pi.atom_indices.assign(
            chemistry.delocalised_participating_atom_indices.begin(),
            chemistry.delocalised_participating_atom_indices.end());
        result.delocalised_pi.orbital_indices.assign(
            chemistry.delocalised_family_orbitals.begin(),
            chemistry.delocalised_family_orbitals.end());
        result.delocalised_pi.family_id=
            chemistry.delocalised_family_id;
        // The generic orbital-chemistry record only establishes pi-family
        // membership.  A topology is authoritative only when a concrete
        // DelocalisedPiAssignment is matched below.
        result.delocalised_pi.topology_available=false;
        result.delocalised_pi.orientation_channels=
            chemistry.delocalised_orientation_channels;
        result.delocalised_pi.cyclic_topology=
            chemistry.delocalised_cyclic_topology;
        result.delocalised_pi.label=delocalised_pi_label(
            chemistry.delocalised_participating_atoms,
            chemistry.delocalised_participating_electrons);
        result.delocalised_pi.source=AnnotationSource::Derived;
        result.delocalised_pi.confidence=
            chemistry.delocalised_pi_confidence;
        result.delocalised_pi.heuristic=false;
    }
    result.heuristic=false;
    return result;
}

std::vector<std::size_t> expand_degenerate(
    std::vector<std::size_t> selected,
    const std::vector<OrbitalLabel>& labels) {
    std::set<std::size_t> expanded(selected.begin(),selected.end());
    for (const auto index:selected) {
        if (index>=labels.size()) continue;
        const auto& label=labels[index];
        if (label.group_size<=1u) continue;
        const std::size_t begin=label.group_base_number>0u
            ?label.group_base_number-1u:index;
        for (std::size_t member=0;member<label.group_size;++member) {
            if (begin+member<labels.size()) expanded.insert(begin+member);
        }
    }
    return {expanded.begin(),expanded.end()};
}

double group_layout_energy(
    const std::vector<OrbitalMetadata>& metadata,
    const std::vector<OrbitalLabel>& labels,
    const std::size_t index) {
    if (index>=metadata.size() || index>=labels.size()) return 0.0;
    const auto& label=labels[index];
    if (label.group_size<=1u) return metadata[index].energy_hartree;
    const std::size_t begin=label.group_base_number>0u
        ?label.group_base_number-1u:index;
    double sum=0.0;
    std::size_t count=0;
    for (std::size_t member=0;member<label.group_size;++member) {
        if (begin+member>=metadata.size()) break;
        sum+=metadata[begin+member].energy_hartree;
        ++count;
    }
    return count>0u?sum/static_cast<double>(count)
                   :metadata[index].energy_hartree;
}

bool transition_metal(const int z) noexcept {
    return (z>=21 && z<=30) || (z>=39 && z<=48) ||
           (z>=72 && z<=80) || (z>=104 && z<=112);
}

enum class LocalLigandPiRole {
    Unresolved,
    SigmaOnly,
    Donor,
    Acceptor,
    Ambiguous,
};

struct LigandScope {
    bool available=false;
    std::size_t metal=0;
    std::set<std::size_t> direct_donors;
    std::set<std::size_t> fragment_atoms;
    LocalLigandPiRole pi_role=LocalLigandPiRole::Unresolved;
};

bool any_metal_ligand_pair(const Wavefunction& wavefunction,
                           const OrbitalPairInteraction& interaction) {
    if (interaction.atom_a>=wavefunction.atoms.size() ||
        interaction.atom_b>=wavefunction.atoms.size()) return false;
    const bool a=transition_metal(
        wavefunction.atoms[interaction.atom_a].atomic_number);
    const bool b=transition_metal(
        wavefunction.atoms[interaction.atom_b].atomic_number);
    return a!=b;
}

LigandScope make_ligand_scope(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment) {
    LigandScope scope;
    if (!environment.available() ||
        environment.metal_atom>=wavefunction.atoms.size()) return scope;
    scope.available=true;
    scope.metal=environment.metal_atom;
    scope.direct_donors.insert(
        environment.ligand_atoms.begin(),environment.ligand_atoms.end());

    // Grow each ligand fragment without crossing a transition-metal centre.
    // This includes O in CO and N in CN for ligand-character analysis while
    // keeping the direct M-donor interaction restricted to the first shell.
    std::vector<std::size_t> pending(
        scope.direct_donors.begin(),scope.direct_donors.end());
    scope.fragment_atoms=scope.direct_donors;
    while (!pending.empty()) {
        const std::size_t atom=pending.back();
        pending.pop_back();
        for (const auto& bond:wavefunction.bond_orders) {
            if (std::abs(bond.mayer_order)<0.05) continue;
            std::size_t other=wavefunction.atoms.size();
            if (bond.atom_a==atom) other=bond.atom_b;
            else if (bond.atom_b==atom) other=bond.atom_a;
            if (other>=wavefunction.atoms.size() || other==scope.metal ||
                transition_metal(wavefunction.atoms[other].atomic_number)) {
                continue;
            }
            if (scope.fragment_atoms.insert(other).second) {
                pending.push_back(other);
            }
        }
    }

    std::size_t sigma_only=0u;
    std::size_t donor=0u;
    std::size_t acceptor=0u;
    std::size_t ambiguous=0u;
    for (const auto atom:scope.direct_donors) {
        if (atom>=wavefunction.atoms.size()) continue;
        const int z=wavefunction.atoms[atom].atomic_number;
        std::size_t external_neighbour_count=0u;
        bool multiple_heavy_bond=false;
        for (const auto& bond:wavefunction.bond_orders) {
            std::size_t other=wavefunction.atoms.size();
            if (bond.atom_a==atom) other=bond.atom_b;
            else if (bond.atom_b==atom) other=bond.atom_a;
            if (other>=wavefunction.atoms.size() || other==scope.metal ||
                transition_metal(wavefunction.atoms[other].atomic_number)) {
                continue;
            }
            const int other_z=wavefunction.atoms[other].atomic_number;
            if (std::abs(bond.mayer_order)>=0.10) {
                ++external_neighbour_count;
            }
            if (other_z>1 && std::abs(bond.mayer_order)>=0.10) {
                multiple_heavy_bond=multiple_heavy_bond ||
                    std::abs(bond.mayer_order)>=1.15;
            }
        }
        if (z==6 && multiple_heavy_bond) {
            ++acceptor; // CO, CN and related unsaturated carbon donors.
        } else if ((z==7 || z==15) && !multiple_heavy_bond &&
                   external_neighbour_count>=3u) {
            ++sigma_only; // Saturated NH3/amine/phosphine-like donors.
        } else if ((z==7 || z==15) && multiple_heavy_bond) {
            ++acceptor;
        } else if (z==9 || z==17 || z==35 || z==53 ||
                   z==8 || z==16 || z==34) {
            ++donor;
        } else {
            ++ambiguous;
        }
    }
    const std::size_t classified=sigma_only+donor+acceptor+ambiguous;
    if (classified==0u) scope.pi_role=LocalLigandPiRole::Unresolved;
    else if (sigma_only==classified) scope.pi_role=LocalLigandPiRole::SigmaOnly;
    else if (acceptor>0u && donor==0u && ambiguous==0u) {
        scope.pi_role=LocalLigandPiRole::Acceptor;
    } else if (donor>0u && acceptor==0u && ambiguous==0u) {
        scope.pi_role=LocalLigandPiRole::Donor;
    } else {
        scope.pi_role=LocalLigandPiRole::Ambiguous;
    }
    return scope;
}

bool scoped_metal_ligand_pair(
    const Wavefunction& wavefunction,
    const OrbitalPairInteraction& interaction,
    const LigandScope* scope) {
    if (scope==nullptr || !scope->available) {
        return any_metal_ligand_pair(wavefunction,interaction);
    }
    const std::size_t a=interaction.atom_a;
    const std::size_t b=interaction.atom_b;
    return (a==scope->metal && scope->direct_donors.count(b)>0u) ||
           (b==scope->metal && scope->direct_donors.count(a)>0u);
}

double bounded_overlap_character(const double value) noexcept {
    // Individual canonical-MO Mulliken terms are basis dependent and may be
    // arbitrarily large for diffuse virtuals.  tanh preserves sign and the
    // linear weak-coupling regime while preventing one outlier from taking
    // over a group score.
    return std::tanh(value);
}

std::string normalised_symmetry(std::string value) {
    value.erase(std::remove_if(value.begin(),value.end(),[](unsigned char c) {
        return c==' ' || c=='\t' || c=='\r' || c=='\n';
    }),value.end());
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char c) {
        if (c>='A' && c<='Z') return static_cast<char>(c-'A'+'a');
        return static_cast<char>(c);
    });
    if (value=="?" || value=="n/a" || value=="na" ||
        value=="none" || value=="-") return {};
    return value;
}

struct GroupCandidate {
    MODiagramLevel level;
    std::size_t base_index=0;
    bool selected_by_reference=false;
    bool selected_by_raw=false;
    bool include=false;
    bool locally_grouped=false;
    bool suppressed_spin_counterpart=false;
    std::uint8_t local_irrep_copy=0;
    bool locally_classified=false;
    std::optional<LocalAngularProjectionWorkspace> local_projection;
    std::optional<LigandFieldEnvironment> local_environment;
};

Spin group_spin(const Wavefunction& wavefunction,
                const GroupCandidate& group) {
    if (group.level.member_indices.empty() ||
        group.level.member_indices.front()>=wavefunction.orbitals.size()) {
        return Spin::Alpha;
    }
    return wavefunction.orbitals[group.level.member_indices.front()].spin;
}

int dominant_metal_family(const MODiagramLevel& level) {
    const std::array<double,3> weights{
        level.metal_s_weight,level.metal_p_weight,level.metal_d_weight};
    const auto found=std::max_element(weights.begin(),weights.end());
    if (found==weights.end()) return -1;
    const int family=static_cast<int>(std::distance(weights.begin(),found));
    const double floor=family==2?kLocalDIrrepWeightFloor:0.08;
    return *found>=floor?family:-1;
}

GroupCandidate merge_local_groups(const Wavefunction& wavefunction,
                                  const GroupCandidate& left,
                                  const GroupCandidate& right) {
    GroupCandidate result=left;
    auto& level=result.level;
    const auto left_count=std::max<std::size_t>(
        1u,left.level.member_indices.size());
    const auto right_count=std::max<std::size_t>(
        1u,right.level.member_indices.size());
    const double total=static_cast<double>(left_count+right_count);
    const auto weighted=[&](const double a,const double b) {
        return (a*static_cast<double>(left_count)+
                b*static_cast<double>(right_count))/total;
    };
    level.layout_energy_hartree=weighted(
        left.level.layout_energy_hartree,right.level.layout_energy_hartree);
    level.metadata.energy_hartree=level.layout_energy_hartree;
    level.total_occupation=left.level.total_occupation+
                           right.level.total_occupation;
    level.metadata.occupation=static_cast<float>(
        level.total_occupation/total);
    level.metal_s_weight=weighted(
        left.level.metal_s_weight,right.level.metal_s_weight);
    level.metal_p_weight=weighted(
        left.level.metal_p_weight,right.level.metal_p_weight);
    level.metal_d_weight=weighted(
        left.level.metal_d_weight,right.level.metal_d_weight);
    level.direct_ligand_p_weight=weighted(
        left.level.direct_ligand_p_weight,
        right.level.direct_ligand_p_weight);
    level.ligand_p_weight=weighted(
        left.level.ligand_p_weight,right.level.ligand_p_weight);
    level.sigma_fraction=weighted(
        left.level.sigma_fraction,right.level.sigma_fraction);
    level.pi_fraction=weighted(
        left.level.pi_fraction,right.level.pi_fraction);
    level.metal_ligand_overlap=weighted(
        left.level.metal_ligand_overlap,right.level.metal_ligand_overlap);
    level.homo=left.level.homo || right.level.homo;
    level.lumo=left.level.lumo || right.level.lumo;
    level.metadata.selected=left.level.metadata.selected ||
                            right.level.metadata.selected;
    level.raw_data_fallback=left.level.raw_data_fallback ||
                            right.level.raw_data_fallback;
    level.member_indices.insert(level.member_indices.end(),
        right.level.member_indices.begin(),right.level.member_indices.end());
    level.member_electrons.insert(level.member_electrons.end(),
        right.level.member_electrons.begin(),right.level.member_electrons.end());
    level.metadata.degeneracy_size=level.member_indices.size();
    if (!level.member_electrons.empty()) level.electrons=level.member_electrons.front();

    double minimum=std::numeric_limits<double>::infinity();
    double maximum=-std::numeric_limits<double>::infinity();
    for (const auto index:level.member_indices) {
        if (index>=wavefunction.orbitals.size()) continue;
        minimum=std::min(minimum,wavefunction.orbitals[index].energy_hartree);
        maximum=std::max(maximum,wavefunction.orbitals[index].energy_hartree);
    }
    level.energy_spread_hartree=
        std::isfinite(minimum) && std::isfinite(maximum)
            ?std::max(0.0,maximum-minimum):0.0;
    if (level.annotation.family=="unavailable" &&
        right.level.annotation.family!="unavailable") {
        level.annotation=right.level.annotation;
    }
    result.selected_by_reference=left.selected_by_reference ||
                                 right.selected_by_reference;
    result.selected_by_raw=left.selected_by_raw || right.selected_by_raw;
    result.include=left.include || right.include;
    result.locally_grouped=true;
    if (left.local_irrep_copy!=right.local_irrep_copy) {
        result.local_irrep_copy=0;
    }
    result.locally_classified=left.locally_classified &&
                              right.locally_classified;
    // Merging changes the actual target span. Never keep the left member's
    // projection proof attached to the enlarged group.
    level.metadata.symmetry_view=aggregate_molecular_symmetry(wavefunction,level.member_indices);
    result.local_irrep_copy=0; result.locally_classified=false;
    if (result.local_projection && result.local_environment) {
        const auto& environment=*result.local_environment;
        const int family=dominant_metal_family(level);
        if(family>=0) {
            level.metadata.symmetry_view=evaluate_local_orbital_symmetry(
                *result.local_projection,environment,level.member_indices,family);
            const auto& assignment=level.metadata.symmetry_view.local_assignment;
            result.local_irrep_copy=assignment?assignment->copy_index:0u;
            result.locally_classified=static_cast<bool>(assignment);
        }
    }
    return result;
}

bool merge_resolved_five_d_run(const Wavefunction& wavefunction,
                               std::vector<GroupCandidate>& groups) {
    constexpr double d_tolerance_hartree=0.005;
    const auto try_run=[&](const std::vector<std::size_t>& run) {
        if (run.size()!=5u) return false;
        std::array<double,4> gaps{};
        for (std::size_t i=0;i<gaps.size();++i) {
            gaps[i]=groups[run[i+1u]].level.layout_energy_hartree-
                    groups[run[i]].level.layout_energy_hartree;
        }
        const auto largest=std::max_element(gaps.begin(),gaps.end());
        const std::size_t split=static_cast<std::size_t>(
            std::distance(gaps.begin(),largest))+1u;
        if (split!=2u && split!=3u) return false;
        double second=0.0;
        for (std::size_t i=0;i<gaps.size();++i) {
            if (i+1u==split) continue;
            second=std::max(second,gaps[i]);
        }
        if (*largest<1.0e-6 || *largest<1.25*second) return false;

        GroupCandidate lower=groups[run.front()];
        for (std::size_t i=1u;i<split;++i) {
            lower=merge_local_groups(wavefunction,lower,groups[run[i]]);
        }
        GroupCandidate upper=groups[run[split]];
        for (std::size_t i=split+1u;i<run.size();++i) {
            upper=merge_local_groups(wavefunction,upper,groups[run[i]]);
        }
        const std::set<std::size_t> consumed(run.begin(),run.end());
        std::vector<GroupCandidate> rebuilt;
        rebuilt.reserve(groups.size()-3u);
        for (std::size_t i=0;i<groups.size();++i) {
            if (i==run.front()) rebuilt.push_back(lower);
            else if (i==run[split]) rebuilt.push_back(upper);
            else if (consumed.count(i)==0u) rebuilt.push_back(groups[i]);
        }
        groups=std::move(rebuilt);
        std::stable_sort(groups.begin(),groups.end(),[&](const auto& a,const auto& b) {
            const Spin a_spin=group_spin(wavefunction,a);
            const Spin b_spin=group_spin(wavefunction,b);
            if (a_spin!=b_spin) return a_spin==Spin::Alpha;
            return a.level.layout_energy_hartree<b.level.layout_energy_hartree;
        });
        return true;
    };

    for (const auto spin:{Spin::Alpha,Spin::Beta}) {
        std::vector<std::size_t> run;
        auto flush=[&]() {
            const bool merged=try_run(run);
            run.clear();
            return merged;
        };
        for (std::size_t i=0;i<groups.size();++i) {
            if (group_spin(wavefunction,groups[i])!=spin) continue;
            const int family=dominant_metal_family(groups[i].level);
            const auto& scope=groups[i].level.metadata.symmetry_view;
            const bool candidate=family==2 &&
                groups[i].level.member_indices.size()==1u &&
                !scope.local_assignment && orbital_symmetry_missing_local_input(scope);
            if (!candidate) {
                if (family==2 && flush()) return true;
                continue;
            }
            if (!run.empty() &&
                groups[i].level.layout_energy_hartree-
                    groups[run.back()].level.layout_energy_hartree>
                        d_tolerance_hartree) {
                if (flush()) return true;
            }
            run.push_back(i);
        }
        if (flush()) return true;
    }
    return false;
}

void merge_local_pseudodegenerate_groups(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment,
    std::vector<GroupCandidate>& groups) {
    if (!environment.available() ||
        !environment.equivalent_ligand_elements) return;
    const std::string point_group=environment.local_point_group();
    const bool unique_two_three_d=
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,2u).has_value() &&
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,3u).has_value();
    if (unique_two_three_d) {
        while (merge_resolved_five_d_run(wavefunction,groups)) {
            // Re-evaluate indices after forming a non-adjacent
            // five-dimensional d run into two resolved subspaces.
        }
    }
    for (std::size_t i=0;i+1u<groups.size();) {
        auto& left=groups[i];
        auto& right=groups[i+1u];
        const std::string left_symmetry=normalised_symmetry(
            left.level.metadata.symmetry_view.label);
        const std::string right_symmetry=normalised_symmetry(
            right.level.metadata.symmetry_view.label);
        const bool local_copies_compatible=
            (!left.locally_classified && !right.locally_classified) ||
            (left.locally_classified && right.locally_classified &&
             left.local_irrep_copy!=0u &&
             left.local_irrep_copy==right.local_irrep_copy);
        const bool labels_compatible=
            (left_symmetry.empty() && right_symmetry.empty()) ||
            (!left_symmetry.empty() && left_symmetry==right_symmetry &&
             compatible_symmetry_scopes(left.level.metadata.symmetry_view,right.level.metadata.symmetry_view) &&
             local_copies_compatible);
        if (group_spin(wavefunction,left)!=group_spin(wavefunction,right) ||
            !labels_compatible) {
            ++i;
            continue;
        }
        const int family=dominant_metal_family(left.level);
        if (family<1 || family!=dominant_metal_family(right.level)) {
            ++i;
            continue;
        }
        const std::size_t left_size=left.level.member_indices.size();
        const std::size_t right_size=right.level.member_indices.size();
        const std::size_t combined=left_size+right_size;
        const bool unlabeled=left_symmetry.empty() && right_symmetry.empty();
        const bool dimension_unambiguous=!unlabeled ||
            classify_local_irrep_by_dimension(
                point_group,static_cast<MetalAOShell>(family),combined)
                .has_value();
        const bool expected_p=family==1 && combined<=3u &&
            (left_size<3u || right_size<3u);
        const bool expected_d=family==2 &&
            ((combined==3u && (left_size==1u || right_size==1u)) ||
             (combined==2u && left_size==1u && right_size==1u));
        const double gap=std::abs(
            right.level.layout_energy_hartree-left.level.layout_energy_hartree);
        const double local_tolerance_hartree=family==1?0.015:0.005;
        if (!dimension_unambiguous || (!expected_p && !expected_d) ||
            gap>local_tolerance_hartree) {
            ++i;
            continue;
        }
        groups[i]=merge_local_groups(wavefunction,left,right);
        groups.erase(groups.begin()+static_cast<std::ptrdiff_t>(i+1u));
    }
}

void recover_local_ligand_field_symmetry(
    const LocalAngularProjectionWorkspace& workspace,
    const LigandFieldEnvironment& environment,
    std::vector<GroupCandidate>& groups) {
    if (!environment.available()) return;
    for(auto& group:groups) {
        auto& level=group.level;
        group.local_projection=workspace;group.local_environment=environment;
        const int family=dominant_metal_family(level);
        if(family<0)continue;
        level.metadata.symmetry_view=evaluate_local_orbital_symmetry(
            workspace,environment,level.member_indices,family);
        const auto& assignment=level.metadata.symmetry_view.local_assignment;
        group.local_irrep_copy=assignment?assignment->copy_index:0u;
        group.locally_classified=static_cast<bool>(assignment);
    }
}

class SpinOverlapWorkspace {
public:
    explicit SpinOverlapWorkspace(const Wavefunction& wavefunction)
        : wavefunction_(wavefunction),
          basis_count_(wavefunction.basis_count),
          alpha_position_(wavefunction.orbitals.size(),invalid_index()),
          beta_position_(wavefunction.orbitals.size(),invalid_index()) {
        if (basis_count_==0u ||
            wavefunction_.ao_overlap.size()!=basis_count_*basis_count_) {
            return;
        }
        for (std::size_t orbital=0;orbital<wavefunction_.orbitals.size();++orbital) {
            const auto& item=wavefunction_.orbitals[orbital];
            if (item.coefficients.size()!=basis_count_) continue;
            if (item.spin==Spin::Beta) {
                beta_position_[orbital]=beta_orbitals_.size();
                beta_orbitals_.push_back(orbital);
            } else {
                alpha_position_[orbital]=alpha_orbitals_.size();
                alpha_orbitals_.push_back(orbital);
            }
        }
        if (alpha_orbitals_.empty() || beta_orbitals_.empty()) return;

        // The old pairwise implementation recomputed S*C_beta for every
        // alpha/beta group candidate.  Open-shell FCHK files contain a full
        // alpha block followed by a full beta block, so that repeated matrix
        // vector product dominated every UI frame.  Transform each beta MO
        // once, then reuse ordinary dot products for every subspace score and
        // member match.  The mathematical S-metric overlap is unchanged.
        transformed_beta_.assign(
            beta_orbitals_.size()*basis_count_,0.0);
        for (std::size_t beta=0;beta<beta_orbitals_.size();++beta) {
            const auto& coefficients=
                wavefunction_.orbitals[beta_orbitals_[beta]].coefficients;
            double* transformed=transformed_beta_.data()+beta*basis_count_;
            for (std::size_t mu=0;mu<basis_count_;++mu) {
                double value=0.0;
                for (std::size_t nu=0;nu<basis_count_;++nu) {
                    value+=wavefunction_.ao_overlap[mu*basis_count_+nu]*
                           static_cast<double>(coefficients[nu]);
                }
                transformed[mu]=value;
            }
        }

        overlaps_.assign(alpha_orbitals_.size()*beta_orbitals_.size(),0.0);
        for (std::size_t alpha=0;alpha<alpha_orbitals_.size();++alpha) {
            const auto& coefficients=
                wavefunction_.orbitals[alpha_orbitals_[alpha]].coefficients;
            for (std::size_t beta=0;beta<beta_orbitals_.size();++beta) {
                const double* transformed=
                    transformed_beta_.data()+beta*basis_count_;
                double value=0.0;
                for (std::size_t mu=0;mu<basis_count_;++mu) {
                    value+=static_cast<double>(coefficients[mu])*
                           transformed[mu];
                }
                overlaps_[alpha*beta_orbitals_.size()+beta]=value;
            }
        }
    }

    [[nodiscard]] double overlap(const std::size_t left,
                                 const std::size_t right) const noexcept {
        if (left>=alpha_position_.size() || right>=beta_position_.size()) {
            return 0.0;
        }
        std::size_t alpha=alpha_position_[left];
        std::size_t beta=beta_position_[right];
        if (alpha==invalid_index() || beta==invalid_index()) {
            if (right>=alpha_position_.size() || left>=beta_position_.size()) {
                return 0.0;
            }
            alpha=alpha_position_[right];
            beta=beta_position_[left];
        }
        if (alpha==invalid_index() || beta==invalid_index() ||
            beta_orbitals_.empty()) return 0.0;
        return overlaps_[alpha*beta_orbitals_.size()+beta];
    }

private:
    [[nodiscard]] static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }

    const Wavefunction& wavefunction_;
    std::size_t basis_count_=0u;
    std::vector<std::size_t> alpha_orbitals_;
    std::vector<std::size_t> beta_orbitals_;
    std::vector<std::size_t> alpha_position_;
    std::vector<std::size_t> beta_position_;
    std::vector<double> transformed_beta_;
    std::vector<double> overlaps_;
};

double subspace_overlap(const SpinOverlapWorkspace& workspace,
                        const std::vector<std::size_t>& alpha_members,
                        const std::vector<std::size_t>& beta_members) {
    if (alpha_members.empty() ||
        alpha_members.size()!=beta_members.size()) {
        return 0.0;
    }
    double squared=0.0;
    for (const auto a:alpha_members) {
        for (const auto b:beta_members) {
            const double overlap=workspace.overlap(a,b);
            squared+=overlap*overlap;
        }
    }
    return squared/static_cast<double>(alpha_members.size());
}

double subspace_overlap(const SpinOverlapWorkspace& workspace,
                        const GroupCandidate& alpha,
                        const GroupCandidate& beta) {
    return subspace_overlap(workspace,alpha.level.member_indices,
                            beta.level.member_indices);
}

void combine_spin_occupations(const Wavefunction& wavefunction,
                              const SpinOverlapWorkspace& workspace,
                              GroupCandidate& alpha,
                              const GroupCandidate& beta,
                              const double occupation_threshold) {
    const std::string alpha_symmetry=normalised_symmetry(
        alpha.level.metadata.symmetry_view.label);
    const std::string beta_symmetry=normalised_symmetry(
        beta.level.metadata.symmetry_view.label);
    if ((alpha_symmetry.empty() || alpha_symmetry=="?" ||
         alpha_symmetry=="n/a") &&
        !beta_symmetry.empty() && beta_symmetry!="?" &&
        beta_symmetry!="n/a") {
        alpha.level.metadata.symmetry_view=candidate_orbital_symmetry(
            beta.level.metadata.symmetry_view,OrbitalSymmetryOrigin::SpinCounterpartCandidate,
            alpha.level.member_indices,subspace_overlap(workspace,alpha,beta));
    }
    struct Match { std::size_t a=0; std::size_t b=0; double score=0.0; };
    std::vector<Match> matches;
    for (std::size_t a=0;a<alpha.level.member_indices.size();++a) {
        for (std::size_t b=0;b<beta.level.member_indices.size();++b) {
            matches.push_back({a,b,std::abs(workspace.overlap(
                alpha.level.member_indices[a],
                beta.level.member_indices[b]))});
        }
    }
    std::sort(matches.begin(),matches.end(),[](const auto& a,const auto& b) {
        return a.score>b.score;
    });
    std::vector<std::size_t> beta_for_alpha(
        alpha.level.member_indices.size(),std::numeric_limits<std::size_t>::max());
    std::set<std::size_t> used_beta;
    for (const auto& match:matches) {
        if (beta_for_alpha[match.a]!=std::numeric_limits<std::size_t>::max() ||
            used_beta.count(match.b)) continue;
        beta_for_alpha[match.a]=match.b;
        used_beta.insert(match.b);
    }
    alpha.level.member_electrons.assign(
        alpha.level.member_indices.size(),ElectronGlyphs{});
    alpha.level.member_spin_counterparts.assign(
        alpha.level.member_indices.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::size_t a=0;a<alpha.level.member_indices.size();++a) {
        const auto alpha_index=alpha.level.member_indices[a];
        if (alpha_index<wavefunction.orbitals.size() &&
            wavefunction.orbitals[alpha_index].occupation>
                occupation_threshold) {
            alpha.level.member_electrons[a].alpha=1;
        }
        const auto matched=beta_for_alpha[a];
        if (matched<beta.level.member_indices.size()) {
            const auto beta_index=beta.level.member_indices[matched];
            alpha.level.member_spin_counterparts[a]=beta_index;
            if (beta_index<wavefunction.orbitals.size() &&
                wavefunction.orbitals[beta_index].occupation>
                    occupation_threshold) {
                alpha.level.member_electrons[a].beta=1;
            }
        }
    }
    if (!alpha.level.member_electrons.empty()) {
        alpha.level.electrons=alpha.level.member_electrons.front();
    }
    alpha.level.total_occupation+=beta.level.total_occupation;
    alpha.level.metadata.occupation=static_cast<float>(
        alpha.level.total_occupation/
        static_cast<double>(std::max<std::size_t>(
            1u,alpha.level.member_indices.size())));
    alpha.level.metadata.selected=alpha.level.metadata.selected ||
                                  beta.level.metadata.selected;
    alpha.selected_by_reference=alpha.selected_by_reference ||
                                beta.selected_by_reference;
    alpha.selected_by_raw=alpha.selected_by_raw || beta.selected_by_raw;
    alpha.include=alpha.include || beta.include;
}

struct SpinCollapseResult {
    std::size_t paired_groups=0u;
    std::size_t paired_members=0u;
};

struct SpinPairCandidate {
    std::size_t alpha=0u;
    std::size_t beta=0u;
    double score=0.0;
};

std::vector<SpinPairCandidate> maximum_cardinality_spin_matching(
    const std::vector<SpinPairCandidate>& candidates,
    const std::size_t group_count) {
    if (candidates.empty()) return {};

    std::vector<std::size_t> alpha_groups;
    std::vector<std::size_t> beta_groups;
    for (const auto& candidate:candidates) {
        alpha_groups.push_back(candidate.alpha);
        beta_groups.push_back(candidate.beta);
    }
    std::sort(alpha_groups.begin(),alpha_groups.end());
    alpha_groups.erase(
        std::unique(alpha_groups.begin(),alpha_groups.end()),alpha_groups.end());
    std::sort(beta_groups.begin(),beta_groups.end());
    beta_groups.erase(
        std::unique(beta_groups.begin(),beta_groups.end()),beta_groups.end());

    const int source=0;
    const int alpha_begin=1;
    const int beta_begin=alpha_begin+static_cast<int>(alpha_groups.size());
    const int sink=beta_begin+static_cast<int>(beta_groups.size());
    const int node_count=sink+1;
    std::vector<int> alpha_node(group_count,-1);
    std::vector<int> beta_node(group_count,-1);
    for (std::size_t i=0;i<alpha_groups.size();++i) {
        if (alpha_groups[i]<group_count) {
            alpha_node[alpha_groups[i]]=alpha_begin+static_cast<int>(i);
        }
    }
    for (std::size_t i=0;i<beta_groups.size();++i) {
        if (beta_groups[i]<group_count) {
            beta_node[beta_groups[i]]=beta_begin+static_cast<int>(i);
        }
    }

    struct FlowEdge {
        int to=0;
        std::size_t reverse=0u;
        int capacity=0;
        double cost=0.0;
    };
    std::vector<std::vector<FlowEdge>> graph(
        static_cast<std::size_t>(node_count));
    const auto add_edge=[&](const int from,const int to,const double cost) {
        const std::size_t forward=graph[static_cast<std::size_t>(from)].size();
        const std::size_t reverse=graph[static_cast<std::size_t>(to)].size();
        graph[static_cast<std::size_t>(from)].push_back(
            {to,reverse,1,cost});
        graph[static_cast<std::size_t>(to)].push_back(
            {from,forward,0,-cost});
        return forward;
    };
    for (const auto group:alpha_groups) {
        add_edge(source,alpha_node[group],0.0);
    }
    for (const auto group:beta_groups) {
        add_edge(beta_node[group],sink,0.0);
    }
    struct PairArc {
        std::size_t candidate=0u;
        int from=0;
        std::size_t edge=0u;
    };
    std::vector<PairArc> pair_arcs;
    pair_arcs.reserve(candidates.size());
    for (std::size_t i=0;i<candidates.size();++i) {
        const int from=alpha_node[candidates[i].alpha];
        const int to=beta_node[candidates[i].beta];
        pair_arcs.push_back({i,from,add_edge(from,to,-candidates[i].score)});
    }

    // Successive shortest augmenting paths continue until no source/sink
    // path remains.  Consequently cardinality is the primary objective; the
    // negative overlap cost makes total S-metric overlap the secondary one.
    // Residual reverse arcs permit an earlier choice to be replaced, which is
    // precisely the reassignment the former sorted-edge greedy pass lacked.
    constexpr double infinity=std::numeric_limits<double>::infinity();
    for (;;) {
        std::vector<double> distance(
            static_cast<std::size_t>(node_count),infinity);
        std::vector<int> previous_node(
            static_cast<std::size_t>(node_count),-1);
        std::vector<std::size_t> previous_edge(
            static_cast<std::size_t>(node_count),0u);
        std::vector<bool> queued(static_cast<std::size_t>(node_count),false);
        std::queue<int> pending;
        distance[static_cast<std::size_t>(source)]=0.0;
        pending.push(source);
        queued[static_cast<std::size_t>(source)]=true;
        while (!pending.empty()) {
            const int from=pending.front();
            pending.pop();
            queued[static_cast<std::size_t>(from)]=false;
            const auto& edges=graph[static_cast<std::size_t>(from)];
            for (std::size_t edge_index=0;edge_index<edges.size();++edge_index) {
                const auto& edge=edges[edge_index];
                if (edge.capacity<=0) continue;
                const double proposed=distance[static_cast<std::size_t>(from)]+
                                      edge.cost;
                if (proposed>=distance[static_cast<std::size_t>(edge.to)]-
                                 1.0e-12) continue;
                distance[static_cast<std::size_t>(edge.to)]=proposed;
                previous_node[static_cast<std::size_t>(edge.to)]=from;
                previous_edge[static_cast<std::size_t>(edge.to)]=edge_index;
                if (!queued[static_cast<std::size_t>(edge.to)]) {
                    pending.push(edge.to);
                    queued[static_cast<std::size_t>(edge.to)]=true;
                }
            }
        }
        if (previous_node[static_cast<std::size_t>(sink)]<0) break;
        for (int node=sink;node!=source;
             node=previous_node[static_cast<std::size_t>(node)]) {
            const int from=previous_node[static_cast<std::size_t>(node)];
            const std::size_t edge_index=
                previous_edge[static_cast<std::size_t>(node)];
            auto& edge=graph[static_cast<std::size_t>(from)][edge_index];
            --edge.capacity;
            ++graph[static_cast<std::size_t>(node)][edge.reverse].capacity;
        }
    }

    std::vector<SpinPairCandidate> result;
    for (const auto& arc:pair_arcs) {
        if (graph[static_cast<std::size_t>(arc.from)][arc.edge].capacity==0) {
            result.push_back(candidates[arc.candidate]);
        }
    }
    std::sort(result.begin(),result.end(),[](const auto& left,const auto& right) {
        if (left.alpha!=right.alpha) return left.alpha<right.alpha;
        return left.beta<right.beta;
    });
    return result;
}

std::optional<std::size_t> formal_irrep_dimension(
    const std::string& point_group,const std::string& symmetry) {
    const auto* definition=find_point_group(point_group);
    if (definition==nullptr || symmetry.empty()) return std::nullopt;
    const auto found=std::find_if(
        definition->irreps.begin(),definition->irreps.end(),[&](const auto& irrep) {
            return normalised_symmetry(std::string(irrep.label))==symmetry;
        });
    if (found==definition->irreps.end() || found->dimension==0u) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(found->dimension);
}

std::vector<std::vector<std::size_t>> group_subsets_with_member_count(
    const std::vector<GroupCandidate>& groups,
    const std::vector<std::size_t>& candidates,
    const std::size_t target) {
    std::vector<std::vector<std::size_t>> result;
    std::vector<std::size_t> current;
    const auto visit=[&](const auto& self,const std::size_t begin,
                         const std::size_t count)->void {
        if (count==target) {
            if (current.size()>=2u) result.push_back(current);
            return;
        }
        for (std::size_t i=begin;i<candidates.size();++i) {
            const std::size_t group=candidates[i];
            const std::size_t size=groups[group].level.member_indices.size();
            if (size==0u || count+size>target) continue;
            current.push_back(group);
            self(self,i+1u,count+size);
            current.pop_back();
        }
    };
    visit(visit,0u,0u);
    return result;
}

GroupCandidate merge_spin_partition(
    const Wavefunction& wavefunction,
    const std::vector<GroupCandidate>& groups,
    const std::vector<std::size_t>& partition) {
    GroupCandidate merged=groups[partition.front()];
    for (std::size_t i=1u;i<partition.size();++i) {
        merged=merge_local_groups(wavefunction,merged,groups[partition[i]]);
    }
    return merged;
}

void collapse_split_spin_partitions(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const SpinOverlapWorkspace& overlap_workspace,
    std::vector<GroupCandidate>& groups,
    std::set<std::size_t>& used_alpha,
    std::set<std::size_t>& used_beta,
    SpinCollapseResult& result) {
    struct ResidualCandidate {
        std::vector<std::size_t> alpha;
        std::vector<std::size_t> beta;
        double score=0.0;
    };
    std::vector<ResidualCandidate> candidates;
    const auto add_full_against_split=[&](const Spin full_spin) {
        const Spin split_spin=full_spin==Spin::Alpha?Spin::Beta:Spin::Alpha;
        const auto& used_full=full_spin==Spin::Alpha?used_alpha:used_beta;
        const auto& used_split=split_spin==Spin::Alpha?used_alpha:used_beta;
        for (std::size_t full=0;full<groups.size();++full) {
            if (group_spin(wavefunction,groups[full])!=full_spin ||
                used_full.count(full)>0u) continue;
            const std::string symmetry=normalised_symmetry(
                groups[full].level.metadata.symmetry_view.label);
            const auto& full_explanation=groups[full].level.metadata.symmetry_view;
            const auto dimension=formal_irrep_dimension(full_explanation.point_group,symmetry);
            if (!dimension || *dimension<2u ||
                groups[full].level.member_indices.size()!=*dimension) continue;
            std::vector<std::size_t> pieces;
            for (std::size_t split=0;split<groups.size();++split) {
                if (group_spin(wavefunction,groups[split])!=split_spin ||
                    used_split.count(split)>0u ||
                    groups[split].level.member_indices.empty() ||
                    groups[split].level.member_indices.size()>=*dimension ||
                    normalised_symmetry(
                        groups[split].level.metadata.symmetry_view.label)!=symmetry) continue;
                if (!compatible_symmetry_scopes(full_explanation,groups[split].level.metadata.symmetry_view))
                    continue;
                pieces.push_back(split);
            }
            for (auto partition:group_subsets_with_member_count(
                     groups,pieces,*dimension)) {
                std::vector<std::size_t> alpha;
                std::vector<std::size_t> beta;
                if (full_spin==Spin::Alpha) {
                    alpha={full};
                    beta=std::move(partition);
                } else {
                    alpha=std::move(partition);
                    beta={full};
                }
                std::vector<std::size_t> alpha_members;
                std::vector<std::size_t> beta_members;
                for (const auto group:alpha) {
                    alpha_members.insert(alpha_members.end(),
                        groups[group].level.member_indices.begin(),
                        groups[group].level.member_indices.end());
                }
                for (const auto group:beta) {
                    beta_members.insert(beta_members.end(),
                        groups[group].level.member_indices.begin(),
                        groups[group].level.member_indices.end());
                }
                const double score=subspace_overlap(
                    overlap_workspace,alpha_members,beta_members);
                if (score>=0.80) {
                    candidates.push_back(
                        {std::move(alpha),std::move(beta),score});
                }
            }
        }
    };
    add_full_against_split(Spin::Alpha);
    add_full_against_split(Spin::Beta);
    std::sort(candidates.begin(),candidates.end(),[](const auto& left,
                                                     const auto& right) {
        if (left.score!=right.score) return left.score>right.score;
        if (left.alpha!=right.alpha) return left.alpha<right.alpha;
        return left.beta<right.beta;
    });

    for (const auto& candidate:candidates) {
        if (std::any_of(candidate.alpha.begin(),candidate.alpha.end(),
                [&](const auto group) {return used_alpha.count(group)>0u;}) ||
            std::any_of(candidate.beta.begin(),candidate.beta.end(),
                [&](const auto group) {return used_beta.count(group)>0u;})) {
            continue;
        }
        GroupCandidate alpha=merge_spin_partition(
            wavefunction,groups,candidate.alpha);
        const GroupCandidate beta=merge_spin_partition(
            wavefunction,groups,candidate.beta);
        combine_spin_occupations(
            wavefunction,overlap_workspace,alpha,beta,
            options.filter.occupation_threshold);

        const std::size_t alpha_primary=candidate.alpha.front();
        groups[alpha_primary]=std::move(alpha);
        groups[alpha_primary].suppressed_spin_counterpart=false;
        for (std::size_t i=1u;i<candidate.alpha.size();++i) {
            auto& piece=groups[candidate.alpha[i]];
            piece.include=false;
            piece.suppressed_spin_counterpart=true;
        }
        for (const auto group:candidate.beta) {
            groups[group].include=false;
            groups[group].suppressed_spin_counterpart=true;
        }
        used_alpha.insert(candidate.alpha.begin(),candidate.alpha.end());
        used_beta.insert(candidate.beta.begin(),candidate.beta.end());
        ++result.paired_groups;
        result.paired_members+=
            groups[alpha_primary].level.member_indices.size();
    }
}

SpinCollapseResult collapse_spin_counterparts(
                                const Wavefunction& wavefunction,
                                const MODiagramOptions& options,
    std::vector<GroupCandidate>& groups) {
    SpinCollapseResult result;
    const bool has_alpha_orbitals=std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const auto& orbital) {return orbital.spin==Spin::Alpha;});
    const bool has_beta_orbitals=std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const auto& orbital) {return orbital.spin==Spin::Beta;});
    if (!has_alpha_orbitals || !has_beta_orbitals ||
        wavefunction.ao_overlap.empty()) return result;
    const SpinOverlapWorkspace overlap_workspace(wavefunction);
    std::vector<SpinPairCandidate> candidates;
    for (std::size_t a=0;a<groups.size();++a) {
        if (group_spin(wavefunction,groups[a])!=Spin::Alpha) continue;
        for (std::size_t b=0;b<groups.size();++b) {
            if (group_spin(wavefunction,groups[b])!=Spin::Beta ||
                groups[a].level.member_indices.size()!=
                    groups[b].level.member_indices.size()) continue;
            const std::string sa=normalised_symmetry(
                groups[a].level.metadata.symmetry_view.label);
            const std::string sb=normalised_symmetry(
                groups[b].level.metadata.symmetry_view.label);
            const bool same_scope=compatible_symmetry_scopes(
                groups[a].level.metadata.symmetry_view,groups[b].level.metadata.symmetry_view);
            if (same_scope && !sa.empty() && !sb.empty() && sa!=sb) continue;
            const int alpha_family=dominant_metal_family(groups[a].level);
            const int beta_family=dominant_metal_family(groups[b].level);
            const bool same_valid_irrep=same_scope && !sa.empty() && !sb.empty() && sa==sb;
            if (!same_valid_irrep && alpha_family>=0 && beta_family>=0 &&
                alpha_family!=beta_family) continue;
            const double score=subspace_overlap(
                overlap_workspace,groups[a],groups[b]);
            // The squared S-metric subspace score is invariant to rotations
            // inside a degenerate block.  UDFT spin polarisation can lower it
            // well below the old 0.50 cutoff even for the same spatial d
            // block; metal-family and symmetry compatibility guard the more
            // permissive floor against unrelated matches.
            if (score>=0.30) candidates.push_back({a,b,score});
        }
    }
    const auto pairs=maximum_cardinality_spin_matching(
        candidates,groups.size());
    std::set<std::size_t> used_alpha;
    std::set<std::size_t> used_beta;
    for (const auto& pair:pairs) {
        used_alpha.insert(pair.alpha);
        used_beta.insert(pair.beta);
        combine_spin_occupations(
            wavefunction,overlap_workspace,
            groups[pair.alpha],groups[pair.beta],
            options.filter.occupation_threshold);
        groups[pair.beta].include=false;
        groups[pair.beta].suppressed_spin_counterpart=true;
        ++result.paired_groups;
        result.paired_members+=groups[pair.alpha].level.member_indices.size();
    }
    collapse_split_spin_partitions(
        wavefunction,options,overlap_workspace,groups,
        used_alpha,used_beta,result);
    return result;
}

double member_metal_weight(const Wavefunction& wavefunction,
                           const MolecularOrbital& orbital,
                           const int angular_momentum,
                           const LigandScope* scope) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index>=wavefunction.atoms.size()) continue;
        if (scope!=nullptr && scope->available) {
            if (contribution.atom_index!=scope->metal) continue;
        } else if (!transition_metal(
                       wavefunction.atoms[contribution.atom_index].atomic_number)) {
            continue;
        }
        if (contribution.angular_momentum==angular_momentum) {
            result+=contribution.weight;
        }
    }
    return result;
}

double member_ligand_p_weight(const Wavefunction& wavefunction,
                              const MolecularOrbital& orbital,
                              const LigandScope* scope,
                              const bool direct_only) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index>=wavefunction.atoms.size()) continue;
        if (scope!=nullptr && scope->available) {
            const auto& allowed=direct_only?scope->direct_donors
                                           :scope->fragment_atoms;
            if (allowed.count(contribution.atom_index)==0u) continue;
        } else if (transition_metal(
                       wavefunction.atoms[contribution.atom_index].atomic_number)) {
            continue;
        }
        if (contribution.angular_momentum==1) result+=contribution.weight;
    }
    return result;
}

GroupCandidate make_group_candidate(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const FrontierOrbitals& frontier,
    const std::vector<OrbitalMetadata>& metadata,
    const std::vector<OrbitalAnnotation>& annotations,
    const std::vector<OrbitalLabel>& labels,
    const std::size_t base,
    const LigandScope* scope,
    const std::size_t forced_group_size=0u,
    const std::vector<MOCurrentRadialShell>* current_shells=nullptr,
    const std::vector<std::size_t>* current_centres=nullptr) {
    GroupCandidate candidate;
    candidate.base_index=base;
    if (base>=wavefunction.orbitals.size() || base>=metadata.size() ||
        base>=labels.size()) return candidate;

    const std::size_t requested=forced_group_size>0u
        ?forced_group_size
        :std::max<std::size_t>(1u,labels[base].group_size);
    const std::size_t end=std::min(wavefunction.orbitals.size(),base+requested);
    auto& level=candidate.level;
    level.metadata=metadata[base];
    level.annotation=annotations[base];
    level.chemistry=wavefunction.orbitals[base].chemistry;
    level.metadata.display_label=std::to_string(base+1u);
    level.metadata.degeneracy_size=end-base;
    level.metadata.selected=false;

    double energy_sum=0.0;
    double energy_min=std::numeric_limits<double>::infinity();
    double energy_max=-std::numeric_limits<double>::infinity();
    double sigma_weighted=0.0;
    double pi_weighted=0.0;
    double channel_weight=0.0;
    double overlap=0.0;
    double valence_weight=0.0;
    bool all_valence=true;
    bool protected_occupied=false;
    bool cation_frontier_group=false;

    for (std::size_t index=base;index<end;++index) {
        const auto& orbital=wavefunction.orbitals[index];
        level.member_indices.push_back(index);
        level.member_electrons.push_back(electron_glyphs_for_orbital(
            orbital,frontier.separate_spin_sets));
        level.total_occupation+=static_cast<double>(orbital.occupation);
        const double display_energy=diagram_orbital_energy(wavefunction,options,index);
        energy_sum+=display_energy;
        energy_min=std::min(energy_min,display_energy);
        energy_max=std::max(energy_max,display_energy);
        level.homo=level.homo || (frontier.homo && *frontier.homo==index);
        level.lumo=level.lumo || (frontier.lumo && *frontier.lumo==index);
        level.metadata.selected=level.metadata.selected ||
                                options.selected_index==index;
        candidate.selected_by_reference=candidate.selected_by_reference ||
            (orbital.chemistry.available &&
             orbital.chemistry.valence_manifold &&
             !confidently_deep_core_orbital(orbital,options.filter));
        protected_occupied=protected_occupied ||
            (occupied(orbital,options.filter.occupation_threshold) &&
             !confidently_deep_core_orbital(orbital,options.filter));
        if (wavefunction.charge>0 &&
            wavefunction.charge_provenance!=DataProvenance::Unavailable &&
            !occupied(orbital,options.filter.occupation_threshold)) {
            cation_frontier_group=cation_frontier_group ||
                (frontier.alpha_lumo && *frontier.alpha_lumo==index) ||
                (frontier.separate_spin_sets && frontier.beta_lumo &&
                 *frontier.beta_lumo==index);
        }
        all_valence=all_valence && orbital.chemistry.available &&
                    orbital.chemistry.valence_manifold;
        valence_weight+=orbital.chemistry.valence_weight;

        level.metal_s_weight+=member_metal_weight(wavefunction,orbital,0,scope);
        level.metal_p_weight+=member_metal_weight(wavefunction,orbital,1,scope);
        level.metal_d_weight+=member_metal_weight(wavefunction,orbital,2,scope);
        level.direct_ligand_p_weight+=member_ligand_p_weight(
            wavefunction,orbital,scope,true);
        level.ligand_p_weight+=member_ligand_p_weight(
            wavefunction,orbital,scope,false);

        for (const auto& interaction:orbital.chemistry.interactions) {
            if (!scoped_metal_ligand_pair(wavefunction,interaction,scope)) continue;
            const double bounded=bounded_overlap_character(
                interaction.overlap_character);
            const double magnitude=std::abs(bounded);
            const double determined=interaction.channel.sigma+
                                    interaction.channel.pi+
                                    interaction.channel.delta+
                                    interaction.channel.phi;
            sigma_weighted+=magnitude*interaction.channel.sigma;
            pi_weighted+=magnitude*interaction.channel.pi;
            channel_weight+=magnitude*determined;
            overlap+=bounded;
        }
    }

    const double count=static_cast<double>(std::max<std::size_t>(1u,end-base));
    level.layout_energy_hartree=energy_sum/count;
    level.metadata.energy_hartree=level.layout_energy_hartree;
    level.energy_spread_hartree=std::max(0.0,energy_max-energy_min);
    level.metadata.occupation=static_cast<float>(level.total_occupation/count);
    level.electrons=level.member_electrons.empty()?ElectronGlyphs{}
                                                   :level.member_electrons.front();
    level.metal_s_weight/=count;
    level.metal_p_weight/=count;
    level.metal_d_weight/=count;
    level.direct_ligand_p_weight/=count;
    level.ligand_p_weight/=count;
    const double donor_count=scope!=nullptr && scope->available
        ?static_cast<double>(std::max<std::size_t>(1u,scope->direct_donors.size()))
        :1.0;
    level.metal_ligand_overlap=std::clamp(
        overlap/(count*donor_count),-1.0,1.0);
    valence_weight/=count;
    if (channel_weight>1.0e-12) {
        level.sigma_fraction=std::clamp(sigma_weighted/channel_weight,0.0,1.0);
        level.pi_fraction=std::clamp(pi_weighted/channel_weight,0.0,1.0);
    }

    level.metadata.symmetry_view=aggregate_molecular_symmetry(
        wavefunction,level.member_indices);

    if(current_shells && current_centres) {
        level.composition=mo_group_composition_ledger(wavefunction,options.routed,
            level.member_indices,*current_centres,*current_shells);
        if(level.composition.complete) {
            level.metal_s_weight=level.composition.centre_current_s;
            level.metal_p_weight=level.composition.centre_current_p;
            level.metal_d_weight=level.composition.centre_current_d;
            level.ligand_p_weight=level.composition.ligand_valence_p;
            level.direct_ligand_p_weight=0;
            if(scope && scope->available)for(const auto member:level.member_indices)
                for(const auto& row:options.routed->mo_composition[member].value->rows)
                    if(row.atom && scope->direct_donors.contains(row.atom-1) &&
                       row.type.rfind("Val",0)==0 && row.angular_l==1)
                        level.direct_ligand_p_weight+=row.weight/count;
        }
    }
    const double metal_valence=level.metal_s_weight+
                               level.metal_p_weight+
                               level.metal_d_weight;
    const bool occupied_group=level.total_occupation>options.filter.occupation_threshold;
    const double virtual_ceiling=frontier.lumo && *frontier.lumo<wavefunction.orbitals.size()
        ?diagram_orbital_energy(wavefunction,options,*frontier.lumo)+
             std::min(0.75,options.filter.virtual_window_hartree)
        :std::numeric_limits<double>::infinity();
    const bool energy_relevant=occupied_group ||
        level.layout_energy_hartree<=virtual_ceiling;
    const bool raw_relevant=energy_relevant &&
        level.layout_energy_hartree>=
            options.filter.core_energy_cutoff_hartree &&
        (metal_valence>=0.08 ||
         (valence_weight>=0.20 &&
          std::abs(level.metal_ligand_overlap)>=
              options.weak_metal_ligand_overlap));
    candidate.selected_by_raw=!candidate.selected_by_reference && raw_relevant;
    level.raw_data_fallback=candidate.selected_by_raw;
    const bool selected_inspection_only=
        options.hide_ligand_centred_intermediates &&
        level.metadata.selected;
    candidate.include=(candidate.selected_by_reference &&
                       level.layout_energy_hartree>=
                           options.filter.core_energy_cutoff_hartree &&
                       energy_relevant) ||
                      candidate.selected_by_raw ||
                      protected_occupied ||
                      cation_frontier_group ||
                      (level.metadata.selected && !selected_inspection_only);

    // A sigma-framework group that the minimal canonical manifold missed is
    // recovered from the complete MO block rather than silently disappearing.
    if (!candidate.include && energy_relevant &&
        level.layout_energy_hartree>=
            options.filter.core_energy_cutoff_hartree &&
        level.sigma_fraction>=0.55 &&
        std::abs(level.metal_ligand_overlap)>=options.weak_metal_ligand_overlap) {
        candidate.include=true;
        candidate.selected_by_raw=true;
        level.raw_data_fallback=true;
    }

    if (all_valence && level.annotation.family=="unavailable") {
        if (level.pi_fraction>=0.55) level.annotation.family="pi";
        else if (level.sigma_fraction>=0.55) level.annotation.family="sigma";
        if (level.annotation.family!="unavailable") {
            level.annotation.family_source=AnnotationSource::Derived;
            level.annotation.family_confidence=std::max(
                level.pi_fraction,level.sigma_fraction);
        }
    }

    // This aggregate is specifically a central-metal/first-shell ligand
    // descriptor.  Outside a resolved ligand scope the accumulated overlap is
    // identically zero, which is absence of evidence rather than evidence for
    // a nonbonding MO.  Preserve the per-orbital chemistry annotation in that
    // case instead of manufacturing a 0%-confidence "nonbonding" label for
    // ordinary main-group and deep-core orbitals.
    const bool measured_metal_ligand_channel=
        scope!=nullptr && scope->available && channel_weight>1.0e-10 &&
        metal_valence>1.0e-6 &&
        (level.direct_ligand_p_weight+level.ligand_p_weight)>1.0e-6;
    if (measured_metal_ligand_channel) {
        const double magnitude=std::abs(level.metal_ligand_overlap);
        const double evidence=std::clamp(
            channel_weight/(count*donor_count*0.08),0.0,1.0);
        BondingClass proposed=BondingClass::Unclassified;
        double proposed_confidence=0.0;
        if (magnitude<options.weak_metal_ligand_overlap) {
            proposed=BondingClass::Nonbonding;
            proposed_confidence=evidence*std::clamp(
                1.0-magnitude/std::max(
                    1.0e-12,options.weak_metal_ligand_overlap),0.0,1.0);
        } else if (level.metal_ligand_overlap>0.0) {
            proposed=BondingClass::Bonding;
            proposed_confidence=evidence*std::clamp(
                magnitude/0.10,0.0,1.0);
        } else {
            proposed=BondingClass::Antibonding;
            proposed_confidence=evidence*std::clamp(
                magnitude/0.10,0.0,1.0);
        }
        // A categorical word at 0% confidence is worse than an honest UND.
        // Only promote this aggregate descriptor once the measured evidence
        // clears the same human-readable confidence floor used elsewhere.
        if (proposed_confidence>=0.20 &&
            proposed_confidence>=level.annotation.bonding_confidence) {
            level.annotation.bonding_class=proposed;
            level.annotation.bonding_confidence=proposed_confidence;
            level.annotation.bonding_source=AnnotationSource::Derived;
        }
    }
    return candidate;
}

struct RawPiPair {
    OrbitalEnergyGapKind gap_kind=OrbitalEnergyGapKind::PiPartner;
    std::shared_ptr<const WeakCrystalFieldAssessment> crystal_field_evidence;
    std::shared_ptr<const PiPartnerAssessment> evidence;
    std::vector<std::string> equivalent_channel_ids;
    std::size_t lower=0;
    std::size_t upper=0;
    PiInteractionKind kind=PiInteractionKind::Coupled;
    double split=0.0;
    double confidence=0.0;
    double display_strength_hartree=std::numeric_limits<double>::quiet_NaN();
    std::size_t retained=0;
    std::string symmetry;
};

bool local_pi_irrep(const std::string& point_group,
                    const std::string& symmetry) {
    if (symmetry.empty() || symmetry=="?" || symmetry=="n/a") return false;
    const auto decomposition=decompose_metal_ao_shell(
        point_group,MetalAOShell::D);
    if (!decomposition) return false;
    return std::any_of(
        decomposition->begin(),decomposition->end(),[&](const auto& copy) {
            return normalised_symmetry(std::string(copy.label))==symmetry;
        });
}

std::vector<RawPiPair> find_pi_pairs(
    const Wavefunction& wavefunction,
    std::vector<GroupCandidate>& groups,
    const MODiagramOptions& options,
    const std::string& point_group,
    const LigandScope* scope,
    std::vector<PiPartnerAssessment>& candidate_evidence) {
    struct ScoredPair { RawPiPair pair; double score=0.0; };
    std::vector<ScoredPair> scored;
    std::set<std::pair<std::string,std::string>> rejected_reasons;
    const LocalLigandPiRole role=scope!=nullptr
        ?scope->pi_role:LocalLigandPiRole::Unresolved;
    const LigandPiPrior prior=role==LocalLigandPiRole::SigmaOnly?LigandPiPrior::SigmaOnly:
        role==LocalLigandPiRole::Donor?LigandPiPrior::Donor:
        role==LocalLigandPiRole::Acceptor?LigandPiPrior::Acceptor:
        role==LocalLigandPiRole::Ambiguous?LigandPiPrior::Ambiguous:LigandPiPrior::Unresolved;
    const bool matched_route=options.routed &&
        options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction);
    const auto same_members=[](std::vector<std::size_t> a,std::vector<std::size_t> b) {
        std::sort(a.begin(),a.end());std::sort(b.begin(),b.end());return a==b;
    };
    const auto routed_members=[&](const auto& a,const auto* b) {
        if(!matched_route||scope==nullptr||!scope->available)return false;
        for(const auto& routed:options.routed->pi_couplings) {
            if(!routed.available()||routed.provider!=RoutedProvider::Nbo)continue;
            const auto& channel=*routed.value;
            if(channel.centre_atoms.size()!=1||channel.centre_atoms.front()!=scope->metal)continue;
            const auto contains=[&](const auto& level){return std::any_of(channel.groups.begin(),channel.groups.end(),
                [&](const auto& group){return same_members(group.members,level.member_indices);});};
            if(contains(a)&&(b==nullptr||contains(*b)))return true;
        }
        return false;
    };
    for (std::size_t lower=0;lower<groups.size();++lower) {
        if (groups[lower].suppressed_spin_counterpart) continue;
        const auto& a=groups[lower].level;
        if (!routed_members(a,static_cast<const MODiagramLevel*>(nullptr)) &&
            (a.pi_fraction<0.60 || a.metal_d_weight+a.ligand_p_weight<0.18)) continue;
        const std::string symmetry_a=normalised_symmetry(a.metadata.symmetry_view.label);
        for (std::size_t upper=lower+1u;upper<groups.size();++upper) {
            if (groups[upper].suppressed_spin_counterpart) continue;
            const auto& b=groups[upper].level;
            const bool routed_pair=routed_members(a,&b);
            if (b.layout_energy_hartree<=a.layout_energy_hartree ||
                group_spin(wavefunction,groups[upper])!=
                    group_spin(wavefunction,groups[lower]) ||
                (!routed_pair&&(b.metadata.degeneracy_size!=a.metadata.degeneracy_size ||
                 b.pi_fraction<0.60 || b.metal_d_weight+b.ligand_p_weight<0.18))) continue;
            const std::string symmetry_b=normalised_symmetry(
                b.metadata.symmetry_view.label);
            const bool a_known=orbital_symmetry_is_local(a.metadata.symmetry_view) &&
                a.metadata.symmetry_view.point_group==point_group && local_pi_irrep(point_group,symmetry_a);
            const bool b_known=orbital_symmetry_is_local(b.metadata.symmetry_view) &&
                b.metadata.symmetry_view.point_group==point_group && local_pi_irrep(point_group,symmetry_b);
            if (a_known && b_known && !compatible_symmetry_scopes(
                    a.metadata.symmetry_view,b.metadata.symmetry_view)) continue;
            std::string symmetry;
            if (a_known && b_known && symmetry_a==symmetry_b) {
                symmetry=symmetry_a;
            } else if (a_known &&
                       (symmetry_b.empty() || symmetry_b=="?" ||
                        symmetry_b=="n/a")) {
                symmetry=symmetry_a;
            } else if (b_known &&
                       (symmetry_a.empty() || symmetry_a=="?" ||
                        symmetry_a=="n/a")) {
                symmetry=symmetry_b;
            } else if(!routed_pair) {
                continue;
            }
            if(routed_pair) {
                const bool a_local=orbital_symmetry_is_local(a.metadata.symmetry_view);
                const bool b_local=orbital_symmetry_is_local(b.metadata.symmetry_view);
                if(a_local&&b_local&&(!compatible_symmetry_scopes(a.metadata.symmetry_view,b.metadata.symmetry_view)||
                   symmetry_a!=symmetry_b))continue;
                // A common multi-rank channel is not a proof of an unknown
                // row's irrep. Keep the relation independent of label recovery.
                if(!a_known||!b_known||symmetry_a!=symmetry_b)symmetry.clear();
            } else if (!local_pi_irrep(point_group,symmetry)) continue;
            std::vector<PiPartnerAssessment> assessments;
            if(matched_route)for(const auto& routed:options.routed->pi_couplings) {
                if(!routed.available()||routed.provider!=RoutedProvider::Nbo)continue;
                const auto& channel=*routed.value;
                if(scope==nullptr||!scope->available||channel.centre_atoms.size()!=1||
                   channel.centre_atoms.front()!=scope->metal)continue;
                const auto first=std::find_if(channel.groups.begin(),channel.groups.end(),
                    [&](const auto& g){return same_members(g.members,a.member_indices);});
                const auto second=std::find_if(channel.groups.begin(),channel.groups.end(),
                    [&](const auto& g){return same_members(g.members,b.member_indices);});
                if(first==channel.groups.end()||second==channel.groups.end()||first->spin!=second->spin) {
                    PiPartnerAssessment missing;
                    missing.input_valid=true;missing.prior=prior;
                    missing.lower={a.layout_energy_hartree,a.metal_d_weight,a.ligand_p_weight,a.pi_fraction,a.metal_ligand_overlap};
                    missing.upper={b.layout_energy_hartree,b.metal_d_weight,b.ligand_p_weight,b.pi_fraction,b.metal_ligand_overlap};
                    missing.splitting_hartree=b.layout_energy_hartree-a.layout_energy_hartree;
                    missing.detail="complete-canonical-member-set-not-matched-in-local-channel";
                    missing.channel.channel_id=channel.id;
                    missing.channel.canonical_fingerprint=options.routed->canonical_fingerprint;
                    missing.channel.spin=nbo_spin_name(channel.spin);
                    missing.channel.operator_kind=channel.operator_kind;
                    // Unmatched scope is already present in routed data. It
                    // is not a candidate for every unrelated pair of levels.
                    continue;
                }
                PiPartnerChannelEvidence provenance;
                provenance.channel_id=channel.id;
                provenance.canonical_fingerprint=options.routed->canonical_fingerprint;
                provenance.spin=nbo_spin_name(channel.spin);
                provenance.lower_members=a.member_indices;provenance.upper_members=b.member_indices;
                provenance.same_operator_verified=true;
                provenance.operator_kind=channel.operator_kind;
                // Use full-molecule source irreps, never recovered local row
                // labels, to veto different strict representations.
                const auto molecular_a=aggregate_molecular_symmetry(wavefunction,a.member_indices);
                const auto molecular_b=aggregate_molecular_symmetry(wavefunction,b.member_indices);
                PiEndpointSymmetryEvidence strict_a,strict_b;
                const bool common_molecular_scope=compatible_symmetry_scopes(molecular_a,molecular_b) &&
                    !orbital_symmetry_is_local(molecular_a) &&
                    ((molecular_a.origin==OrbitalSymmetryOrigin::MolecularOperations &&
                      molecular_a.molecular_assignment && molecular_b.molecular_assignment) ||
                     (molecular_a.origin==OrbitalSymmetryOrigin::Producer &&
                      !molecular_a.source_path.empty() && molecular_a.source_line_begin>0));
                if(common_molecular_scope) {
                    strict_a={true,options.routed->canonical_fingerprint+":full-molecular-operation-domain",
                        normalised_symmetry(molecular_a.label)};
                    strict_b={true,strict_a.scope_id,normalised_symmetry(molecular_b.label)};
                }
                provenance.mode_assessment=assess_pi_mode_pair(channel,a.member_indices,b.member_indices,
                    strict_a,strict_b);
                if(!provenance.mode_assessment.two_endpoint_relation)continue;
                const auto valid_members=[&](const auto& level) {
                    return level.member_indices.size()==level.metadata.degeneracy_size &&
                        std::all_of(level.member_indices.begin(),level.member_indices.end(),[&](auto i){
                            return i<wavefunction.orbitals.size() &&
                                (channel.canonical_members_are_verified_shared_spatial?
                                    wavefunction.orbitals[i].spin==Spin::Alpha:
                                    ((channel.spin==NboSpin::Beta)==(wavefunction.orbitals[i].spin==Spin::Beta)));});
                };
                provenance.complete_membership_verified=valid_members(a)&&valid_members(b);
                provenance.occupations_verified=channel.occupation_status=="available";
                provenance.direction_reference=channel.direction_reference;
                provenance.frozen_operator=channel.frozen_operator;
                provenance.direction=provenance.mode_assessment.direction;
                provenance.direction_verified=provenance.mode_assessment.direction_verified;
                provenance.direction_reference+="; linked-source-edge/common-mode: "+
                    provenance.mode_assessment.reason;
                provenance.lower_character=first->character;provenance.upper_character=second->character;
                provenance.lower_cross_fock_max_hartree=first->cross_fock_max_hartree;
                provenance.upper_cross_fock_min_hartree=second->cross_fock_min_hartree;
                provenance.operator_error_hartree=channel.operator_max_error_hartree;
                provenance.operator_tolerance_hartree=channel.operator_validation_tolerance_hartree;
                assessments.push_back(assess_pi_channel_partner(
                    {a.layout_energy_hartree,a.metal_d_weight,a.ligand_p_weight,a.pi_fraction,a.metal_ligand_overlap},
                    {b.layout_energy_hartree,b.metal_d_weight,b.ligand_p_weight,b.pi_fraction,b.metal_ligand_overlap},
                    prior,provenance));
            }
            if(assessments.empty() && matched_route)continue;
            if(assessments.empty())assessments.push_back(assess_pi_partner(
                {a.layout_energy_hartree,a.metal_d_weight,a.ligand_p_weight,a.pi_fraction,a.metal_ligand_overlap},
                {b.layout_energy_hartree,b.metal_d_weight,b.ligand_p_weight,b.pi_fraction,b.metal_ligand_overlap},
                prior,options.weak_pi_split_hartree,options.weak_metal_ligand_overlap));
            for(auto& assessment:assessments) {
            assessment.channel.lower_members=a.member_indices;
            assessment.channel.upper_members=b.member_indices;
            if(assessment.channel.canonical_fingerprint.empty())
                assessment.channel.canonical_fingerprint=matched_route?
                    options.routed->canonical_fingerprint:nbo_canonical_fingerprint(wavefunction);
            if(assessment.channel.spin.empty())assessment.channel.spin=
                group_spin(wavefunction,groups[lower])==Spin::Beta?"beta":"alpha";
            if(assessment.accepted || rejected_reasons.insert(
                {assessment.channel.channel_id,assessment.detail}).second)
                candidate_evidence.push_back(assessment);
            auto evidence=std::make_shared<const PiPartnerAssessment>(std::move(assessment));
            if(!evidence->accepted)continue;
            const double score=evidence->ranking_score;
            RawPiPair pair;pair.evidence=evidence;
            pair.lower=lower;pair.upper=upper;pair.split=evidence->splitting_hartree;
            pair.symmetry=symmetry.empty()?"":(a_known?a.metadata.symmetry_view.label:b.metadata.symmetry_view.label);
            switch(evidence->direction) {
                case PiPairDirection::Donor:pair.kind=PiInteractionKind::Donor;break;
                case PiPairDirection::Acceptor:pair.kind=PiInteractionKind::Acceptor;break;
                case PiPairDirection::WeakNearNonbonding:pair.kind=PiInteractionKind::WeakNearNonbonding;break;
                default:pair.kind=PiInteractionKind::Coupled;break;
            }
            pair.confidence=evidence->support_score;
            const double a_metal=a.metal_s_weight+a.metal_p_weight+a.metal_d_weight;
            const double b_metal=b.metal_s_weight+b.metal_p_weight+b.metal_d_weight;
            pair.retained=b_metal>a_metal?upper:lower;
            scored.push_back({pair,score});
            }
        }
    }

    // A weak separation between different local d irreps is a crystal-field
    // gap. It is independent of same-irrep ligand pi partner matching and
    // never changes the pi pair count or overwrites the whole-MO annotation.
    {
    for (std::size_t first=0;first<groups.size();++first) {
        if (groups[first].suppressed_spin_counterpart) continue;
        const auto& a=groups[first].level;
        if (a.metal_d_weight<0.60) continue;
        for (std::size_t second=first+1u;second<groups.size();++second) {
            if (groups[second].suppressed_spin_counterpart) continue;
            const auto& b=groups[second].level;
            if (b.metal_d_weight<0.60 ||
                group_spin(wavefunction,groups[second])!=
                    group_spin(wavefunction,groups[first])) continue;
            const std::string sa=normalised_symmetry(a.metadata.symmetry_view.label);
            const std::string sb=normalised_symmetry(b.metadata.symmetry_view.label);
            if (!orbital_symmetry_is_local(a.metadata.symmetry_view) ||
                !orbital_symmetry_is_local(b.metadata.symmetry_view) ||
                !compatible_symmetry_scopes(a.metadata.symmetry_view,b.metadata.symmetry_view)) continue;
            const bool tetrahedral=(sa=="e" && sb=="t2") ||
                                   (sa=="t2" && sb=="e");
            const bool octahedral=(sa=="eg" && sb=="t2g") ||
                                  (sa=="t2g" && sb=="eg");
            const auto local_group=normalised_symmetry(a.metadata.symmetry_view.point_group);
            if (!(tetrahedral && local_group=="td") && !(octahedral && local_group=="oh")) continue;
            auto assessment=assess_weak_crystal_field(
                {a.layout_energy_hartree,a.metal_d_weight,a.ligand_p_weight,a.pi_fraction,a.metal_ligand_overlap},
                {b.layout_energy_hartree,b.metal_d_weight,b.ligand_p_weight,b.pi_fraction,b.metal_ligand_overlap},
                options.weak_crystal_field_split_hartree,options.weak_crystal_field_overlap);
            if (assessment.first.energy_hartree>assessment.second.energy_hartree)
                std::swap(assessment.first,assessment.second);
            auto evidence=std::make_shared<const WeakCrystalFieldAssessment>(std::move(assessment));
            if(!evidence->accepted)continue;
            RawPiPair pair;pair.gap_kind=OrbitalEnergyGapKind::CrystalField;
            pair.crystal_field_evidence=evidence;
            if (a.layout_energy_hartree<=b.layout_energy_hartree) {
                pair.lower=first;
                pair.upper=second;
            } else {
                pair.lower=second;
                pair.upper=first;
            }
            pair.kind=PiInteractionKind::WeakNearNonbonding;
            pair.split=evidence->splitting_hartree;
            pair.confidence=evidence->support_score;
            pair.retained=a.metal_d_weight>=b.metal_d_weight?first:second;
            pair.symmetry=tetrahedral?"E/T2":"Eg/T2g";
            scored.push_back({pair,10.0+pair.confidence});
        }
    }
    }

    std::sort(scored.begin(),scored.end(),[](const auto& a,const auto& b) {
        if(a.score!=b.score)return a.score>b.score;
        const std::string a_id=a.pair.evidence?a.pair.evidence->channel.channel_id:"";
        const std::string b_id=b.pair.evidence?b.pair.evidence->channel.channel_id:"";
        if(a_id!=b_id)return a_id<b_id;
        if(a.pair.lower!=b.pair.lower)return a.pair.lower<b.pair.lower;
        return a.pair.upper<b.pair.upper;
    });
    const auto source_channel=[&](const RawPiPair& pair)->const NboPiCoupling* {
        if(!matched_route || !pair.evidence)return nullptr;
        for(const auto& row:options.routed->pi_couplings)
            if(row.available() && row.value->id==pair.evidence->channel.channel_id)return &*row.value;
        return nullptr;
    };
    const auto same_projector=[](const NboMatrix& a,const NboMatrix& b) {
        if(!a.rows || !a.columns || a.rows!=b.rows || a.columns!=b.columns ||
           a.values.size()!=a.rows*a.columns || b.values.size()!=b.rows*b.columns)return false;
        // Both bases were orthonormalised in the same verified NAO metric.
        // Compare spaces, permitting phase, member permutation and rotation.
        double residual=0;
        for(std::size_t j=0;j<b.columns;++j) {
            std::vector<double> coordinates(a.columns,0);
            for(std::size_t k=0;k<a.columns;++k)for(std::size_t r=0;r<a.rows;++r)
                coordinates[k]+=a.values[r*a.columns+k]*b.values[r*b.columns+j];
            for(std::size_t r=0;r<a.rows;++r) {
                double value=b.values[r*b.columns+j];
                for(std::size_t k=0;k<a.columns;++k)value-=a.values[r*a.columns+k]*coordinates[k];
                residual+=value*value;
            }
        }
        return residual<4e-10*double(b.columns);
    };
    const auto equivalent_relation=[&](const RawPiPair& a,const RawPiPair& b) {
        if(a.gap_kind!=OrbitalEnergyGapKind::PiPartner || b.gap_kind!=a.gap_kind ||
           a.lower!=b.lower || a.upper!=b.upper || a.kind!=b.kind)return false;
        const auto* ca=source_channel(a);const auto* cb=source_channel(b);
        return ca && cb && ca->spin==cb->spin && ca->operator_kind==cb->operator_kind &&
            ca->direction_verified==cb->direction_verified &&
            same_projector(ca->centre_projector_basis,cb->centre_projector_basis) &&
            same_projector(ca->ligand_projector_basis,cb->ligand_projector_basis);
    };
    std::set<std::size_t> used_crystal_field;
    std::vector<RawPiPair> result;
    for (const auto& item:scored) {
        const auto duplicate=std::find_if(result.begin(),result.end(),[&](const auto& prior) {
            return equivalent_relation(prior,item.pair);
        });
        if(duplicate!=result.end()) {
            duplicate->equivalent_channel_ids.push_back(item.pair.evidence->channel.channel_id);
            continue;
        }
        // A canonical group can participate in more than one physical channel.
        // Only the independent weak crystal-field screen is a disjoint match.
        if(item.pair.gap_kind==OrbitalEnergyGapKind::CrystalField) {
            if(used_crystal_field.count(item.pair.lower)||used_crystal_field.count(item.pair.upper))continue;
            used_crystal_field.insert(item.pair.lower);used_crystal_field.insert(item.pair.upper);
        }
        const bool ordinary_eligible=!item.pair.evidence ||
            item.pair.evidence->channel.channel_id.empty() ||
            item.pair.evidence->channel.mode_assessment.ordinary_display_eligible;
        if(options.aomo_scope!=1 || ordinary_eligible) {
            groups[item.pair.lower].include=true;
            groups[item.pair.upper].include=true;
        }
        for (const auto group_index:{item.pair.lower,item.pair.upper}) {
            // Shared participation in a multi-rank Fock channel cannot
            // independently establish a missing row symmetry label.
            if(item.pair.evidence&&!item.pair.evidence->channel.channel_id.empty())continue;
            auto& symmetry=groups[group_index].level.metadata.symmetry_view.label;
            const std::string current=normalised_symmetry(symmetry);
            if ((current.empty() || current=="?" || current=="n/a") &&
                !item.pair.symmetry.empty()) {
                const auto other=group_index==item.pair.lower?item.pair.upper:item.pair.lower;
                const auto& source=groups[other].level.metadata.symmetry_view;
                // A compound E/T2 gap name describes two different levels; it
                // is never an individual row's orbital irrep.
                if (item.pair.symmetry.find('/')==std::string::npos && !source.label.empty()) {
                    groups[group_index].level.metadata.symmetry_view=candidate_orbital_symmetry(
                        source,OrbitalSymmetryOrigin::PiPartnerCandidate,
                        groups[group_index].level.member_indices,item.score);
                    groups[group_index].level.metadata.symmetry_view.pi_partner_evidence=item.pair.evidence;
                }
            }
        }
        result.push_back(item.pair);
        if(item.pair.evidence)result.back().equivalent_channel_ids.push_back(item.pair.evidence->channel.channel_id);
    }
    return result;
}

bool compact_pi_family_eligible(
    const Wavefunction& wavefunction,
    const DelocalisedPiAssignment& assignment) {
    if (assignment.orbitals.empty()) return false;
    if (assignment.atoms.size()<3u) return false;

    // A fully occupied transverse ligand-p manifold in XeF2/I3-/XeF4 is not
    // the active 3c sigma space and must not take over the compact diagram.
    const bool has_virtual=std::any_of(
        assignment.orbitals.begin(),assignment.orbitals.end(),
        [&](const auto orbital) {
            return orbital<wavefunction.orbitals.size() &&
                !occupied(wavefunction.orbitals[orbital],0.05);
        });
    if (!has_virtual) return false;
    switch (assignment.topology) {
        case DelocalisedPiTopology::Path:
        case DelocalisedPiTopology::Cycle:
        case DelocalisedPiTopology::Spiro:
        case DelocalisedPiTopology::HapticMetal:
        case DelocalisedPiTopology::SymmetryDirectSum:
        case DelocalisedPiTopology::MultiChannel:
            return true;
        case DelocalisedPiTopology::BranchedResonance:
            // A three-arm donor star and a four-centre resonance projector are
            // graph- and electron-count-isomorphic.  Use the measured
            // occupation of the centre's selected oriented-p column: weak
            // BF3-like donation remains in the full valence view, whereas a
            // genuinely occupied carbonate/nitrate/guanidinium centre is a
            // useful pi focus.  The interval below 0.90 is intentionally the
            // conservative side of the cross-family calibration rather than
            // an element or molecule-name exception.
            return assignment.branch_centre_projected_occupation>=0.90;
        default:
            return assignment.cyclic_topology;
    }
}

std::set<std::size_t> compact_pi_family_indices(
    const Wavefunction& wavefunction) {
    std::set<std::size_t> result;
    for (const auto& assignment:wavefunction.delocalised_pi_assignments) {
        if (!compact_pi_family_eligible(wavefunction,assignment)) continue;
        for (const auto orbital:assignment.orbitals) {
            if (orbital<wavefunction.orbitals.size()) result.insert(orbital);
        }
    }
    return result;
}

std::set<std::size_t> compact_multicentre_indices(
    const Wavefunction& wavefunction) {
    std::set<std::size_t> result;
    std::set<std::string> seen_subspaces;
    const auto conjugated_pi=compact_pi_family_indices(wavefunction);
    for (const auto& assignment:wavefunction.multicentre_assignments) {
        if (assignment.orbitals.empty() ||
            assignment.provenance==DataProvenance::Unavailable) continue;
        const bool shared_bridge_subspace=
            !assignment.source_subspace_id.empty();
        if (!assignment.geometry_qualified_framework &&
            !shared_bridge_subspace) {
            continue;
        }
        // Generic three-centre fallback can overlap a complete conjugated-pi
        // projector (acrolein is a real example).  Such an overlap is one
        // active space described twice, so pi wins.  Genuine axial 3c4e sets
        // in XeF2/I3- are disjoint from their fully occupied transverse p
        // manifold, and equivalent bridge channels carry an explicit shared
        // source id, so both retain multicentre priority without name rules.
        const bool overlaps_conjugated_pi=std::any_of(
            assignment.orbitals.begin(),assignment.orbitals.end(),
            [&](const auto orbital){return conjugated_pi.count(orbital)!=0u;});
        if (overlaps_conjugated_pi && !shared_bridge_subspace) {
            continue;
        }
        if (shared_bridge_subspace &&
            !seen_subspaces.insert(assignment.source_subspace_id).second) {
            continue;
        }
        for (const auto orbital:assignment.orbitals) {
            if (orbital<wavefunction.orbitals.size()) result.insert(orbital);
        }
    }
    return result;
}

std::set<std::size_t> compact_active_indices(
    const Wavefunction& wavefunction,
    const MODiagramMode mode) {
    if (mode==MODiagramMode::DelocalisedPiFamilyOnly) {
        return compact_pi_family_indices(wavefunction);
    }
    if (mode==MODiagramMode::MulticentreActiveSpaceOnly) {
        return compact_multicentre_indices(wavefunction);
    }
    return {};
}

} // namespace

const char* mo_diagram_mode_name(const MODiagramMode mode) noexcept {
    switch (mode) {
        case MODiagramMode::DelocalisedPiFamilyOnly:
            return "delocalised-pi-family-only";
        case MODiagramMode::MulticentreActiveSpaceOnly:
            return "multicentre-active-space-only";
        default:
            return "valence-central";
    }
}

const char* mo_diagram_mode_title(const MODiagramMode mode) noexcept {
    switch (mode) {
        case MODiagramMode::DelocalisedPiFamilyOnly:
            return "Delocalised pi MO diagram";
        case MODiagramMode::MulticentreActiveSpaceOnly:
            return "Multicentre active-space MO diagram";
        default:
            return "Valence MO diagram";
    }
}

MODiagramMode preferred_compact_mo_diagram_mode(
    const Wavefunction& wavefunction,
    const bool compact) noexcept {
    if (!compact) return MODiagramMode::ValenceCentral;
    try {
        const bool haptic_pi_focus=std::any_of(
            wavefunction.delocalised_pi_assignments.begin(),
            wavefunction.delocalised_pi_assignments.end(),
            [&](const auto& assignment) {
                return assignment.topology==
                           DelocalisedPiTopology::HapticMetal &&
                       compact_pi_family_eligible(wavefunction,assignment);
            });
        if (haptic_pi_focus) {
            return MODiagramMode::DelocalisedPiFamilyOnly;
        }
        if (analyse_ligand_field_environment(wavefunction).available()) {
            return MODiagramMode::ValenceCentral;
        }
        if (!compact_multicentre_indices(wavefunction).empty()) {
            return MODiagramMode::MulticentreActiveSpaceOnly;
        }
        return compact_pi_family_indices(wavefunction).empty()
            ?MODiagramMode::ValenceCentral
            :MODiagramMode::DelocalisedPiFamilyOnly;
    } catch (...) {
        return MODiagramMode::ValenceCentral;
    }
}

const char* pi_interaction_kind_name(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return "pi donation";
        case PiInteractionKind::Acceptor: return "pi back-donation";
        case PiInteractionKind::WeakNearNonbonding:
            return "weak pi interaction; approximately nonbonding";
        default: return "pi coupling";
    }
}

OrbitalAnnotation annotate_orbital(const MolecularOrbital& orbital) {
    if (orbital.chemistry.available) {
        return chemistry_annotation(orbital);
    }
    return annotate_orbital_legacy(orbital);
}

DiagramSelectionPlan build_valence_selection_plan(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata) {
    if (!chemistry_available(wavefunction)) {
        return build_valence_selection_plan_legacy(
            wavefunction,options,metadata);
    }

    DiagramSelectionPlan plan;
    if (options.mode==MODiagramMode::DelocalisedPiFamilyOnly ||
        options.mode==MODiagramMode::MulticentreActiveSpaceOnly) {
        const auto family=compact_active_indices(wavefunction,options.mode);
        plan.included_indices.assign(family.begin(),family.end());
        plan.hidden_count=metadata.size()>plan.included_indices.size()
            ?metadata.size()-plan.included_indices.size():0u;
        for (const auto index:plan.included_indices) {
            if (index>=wavefunction.orbitals.size()) continue;
            if (occupied(wavefunction.orbitals[index],
                         options.filter.occupation_threshold)) {
                ++plan.valence_occupied_count;
            } else {
                ++plan.frontier_virtual_count;
            }
        }
        plan.summary=plan.included_indices.empty()
            ?"compact active space unavailable; valence fallback required"
            :(options.mode==MODiagramMode::DelocalisedPiFamilyOnly
                ?"complete delocalised pi family only"
                :"complete multicentre active space only");
        return plan;
    }
    if(options.routed &&
       options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction) &&
       std::any_of(options.routed->mo_composition.begin(),
                   options.routed->mo_composition.end(),[](const auto& result){
                       return result.available()&&result.provider==RoutedProvider::Nbo;
                   })){
        const auto frontier=find_frontier_orbitals(
            wavefunction.orbitals,options.filter.occupation_threshold);
        const auto labels=build_orbital_labels(wavefunction.orbitals,
            point_group_limited_degeneracy(wavefunction,options.degeneracy));
        std::size_t core_hidden=0,unresolved_occupied=0;
        for(std::size_t i=0;i<wavefunction.orbitals.size();++i){
            const auto& mo=wavefunction.orbitals[i];
            const bool occ=occupied(mo,options.filter.occupation_threshold);
            const auto balance=routed_mo_shell_balance(*options.routed,i);
            const bool deep_core=balance&&balance->complete&&occ&&
                balance->core>=0.70&&balance->core>=balance->valence+0.10;
            if(deep_core){++core_hidden;continue;}
            const bool frontier_virtual=(frontier.alpha_lumo&&*frontier.alpha_lumo==i)||
                (frontier.beta_lumo&&*frontier.beta_lumo==i);
            const bool valence_virtual=balance&&balance->valence>=0.05;
            if(occ||frontier_virtual||valence_virtual){
                plan.included_indices.push_back(i);
                if(occ){++plan.valence_occupied_count;if(!balance)++unresolved_occupied;}
                else ++plan.frontier_virtual_count;
            }
        }
        plan.included_indices=expand_degenerate(std::move(plan.included_indices),labels);
        plan.hidden_count=metadata.size()>plan.included_indices.size()?metadata.size()-
            plan.included_indices.size():0;
        std::ostringstream summary;
        summary<<"validated NAO Cor/Val/Ryd valence scope; "
               <<plan.included_indices.size()<<'/'<<metadata.size()
               <<" canonical MOs; complete deep-core hidden="<<core_hidden
               <<"; unresolved occupied retained="<<unresolved_occupied
               <<"; complete degeneracy groups retained";
        plan.summary=summary.str();return plan;
    }
    const LigandFieldEnvironment ligand_field=
        analyse_ligand_field_environment(wavefunction);
    const LigandScope ligand_scope=make_ligand_scope(
        wavefunction,ligand_field);
    const LigandScope* scope=ligand_scope.available?&ligand_scope:nullptr;
    const auto labels=build_orbital_labels(
        wavefunction.orbitals,
        point_group_limited_degeneracy(wavefunction,options.degeneracy));
    const auto frontier=find_frontier_orbitals(
        wavefunction.orbitals,options.filter.occupation_threshold);
    std::vector<std::size_t> cation_frontier;
    if (wavefunction.charge>0 &&
        wavefunction.charge_provenance!=DataProvenance::Unavailable) {
        if (frontier.alpha_lumo) cation_frontier.push_back(*frontier.alpha_lumo);
        if (frontier.separate_spin_sets && frontier.beta_lumo) {
            cation_frontier.push_back(*frontier.beta_lumo);
        }
        cation_frontier=expand_degenerate(
            std::move(cation_frontier),labels);
    }
    std::size_t raw_supplements=0u;
    std::size_t occupied_safeguards=0u;
    std::size_t cation_vacancies=0u;
    for (std::size_t i=0;i<wavefunction.orbitals.size();++i) {
        const auto& orbital=wavefunction.orbitals[i];
        if (orbital.chemistry.available &&
            orbital.chemistry.valence_manifold &&
            !confidently_deep_core_orbital(orbital,options.filter)) {
            plan.included_indices.push_back(i);
            continue;
        }

        // The minimal S-metric reference chooses a compact chemical subspace;
        // it is not permission to discard occupied canonical MOs.  Only a
        // confidently assigned deep-core level may be folded automatically.
        if (occupied(orbital,options.filter.occupation_threshold) &&
            !confidently_deep_core_orbital(orbital,options.filter)) {
            plan.included_indices.push_back(i);
            ++occupied_safeguards;
            continue;
        }
        if (std::binary_search(cation_frontier.begin(),
                               cation_frontier.end(),i)) {
            plan.included_indices.push_back(i);
            ++cation_vacancies;
            continue;
        }

        double metal=0.0;
        for (const auto& contribution:orbital.chemistry.ao_contributions) {
            const bool selected_metal=scope!=nullptr
                ?contribution.atom_index==scope->metal
                :(contribution.atom_index<wavefunction.atoms.size() &&
                  transition_metal(wavefunction.atoms[
                      contribution.atom_index].atomic_number));
            if (selected_metal && contribution.angular_momentum<=2) {
                metal+=contribution.weight;
            }
        }
        double interaction=0.0;
        for (const auto& pair:orbital.chemistry.interactions) {
            if (scoped_metal_ligand_pair(wavefunction,pair,scope)) {
                interaction+=std::abs(bounded_overlap_character(
                    pair.overlap_character));
            }
        }
        if (metal>=0.08 || interaction>=options.weak_metal_ligand_overlap ||
            (i==options.selected_index &&
             !options.hide_ligand_centred_intermediates)) {
            plan.included_indices.push_back(i);
            ++raw_supplements;
        }
    }
    plan.included_indices=expand_degenerate(
        std::move(plan.included_indices),labels);
    plan.hidden_count=metadata.size()>plan.included_indices.size()
        ?metadata.size()-plan.included_indices.size():0u;

    for (const auto index:plan.included_indices) {
        if (index>=wavefunction.orbitals.size()) continue;
        if (occupied(wavefunction.orbitals[index],
                     options.filter.occupation_threshold)) {
            ++plan.valence_occupied_count;
        } else {
            ++plan.frontier_virtual_count;
        }
    }

    std::ostringstream summary;
    summary<<"chemical-valence reference: "
           <<plan.included_indices.size()<<'/'<<metadata.size()
           <<" canonical MOs; occupied="<<plan.valence_occupied_count
           <<"; virtual="<<plan.frontier_virtual_count
           <<"; raw-MO supplements="<<raw_supplements
           <<"; occupied safeguards="<<occupied_safeguards
           <<"; cation vacancies="<<cation_vacancies
           <<"; confident deep-core/polarisation/Rydberg hidden";
    plan.summary=summary.str();
    return plan;
}

static MODiagramData build_mo_diagram_data_impl(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options) {
    if ((options.mode==MODiagramMode::DelocalisedPiFamilyOnly ||
         options.mode==MODiagramMode::MulticentreActiveSpaceOnly) &&
        compact_active_indices(wavefunction,options.mode).empty()) {
        MODiagramOptions fallback=options;
        fallback.mode=MODiagramMode::ValenceCentral;
        auto data=build_mo_diagram_data_impl(wavefunction,fallback);
        data.selection.summary=
            "compact active space unavailable; "+data.selection.summary;
        return data;
    }
    if (!chemistry_available(wavefunction)) {
        auto data=build_mo_diagram_data_legacy(wavefunction,options);
        attach_local_geometry_metadata(wavefunction,data);
        attach_electronic_state_metadata(wavefunction,data);
        return data;
    }

    MODiagramData data;
    data.mode=options.mode;
    data.plan=choose_diagram_plan(wavefunction);
    data.plan.machine_reason=
        "COV S-metric minimal atomic chemical-valence reference";
    data.frontier=find_frontier_orbitals(
        wavefunction.orbitals,options.filter.occupation_threshold);
    if(options.canonical_display_energies.size()==wavefunction.orbitals.size()) {
        const bool separate=data.frontier.separate_spin_sets;
        data.frontier={};data.frontier.separate_spin_sets=separate;
        for(std::size_t i=0;i<wavefunction.orbitals.size();++i) {
            const auto& orbital=wavefunction.orbitals[i];
            const bool filled=occupied(orbital,options.filter.occupation_threshold);
            const auto choose=[&](std::optional<std::size_t>& entry) {
                if(!entry || (filled?options.canonical_display_energies[i]>options.canonical_display_energies[*entry]:
                                    options.canonical_display_energies[i]<options.canonical_display_energies[*entry]))entry=i;
            };
            choose(filled?data.frontier.homo:data.frontier.lumo);
            if(orbital.spin==Spin::Beta)choose(filled?data.frontier.beta_homo:data.frontier.beta_lumo);
            else choose(filled?data.frontier.alpha_homo:data.frontier.alpha_lumo);
        }
    }
    data.metadata=build_orbital_metadata(
        wavefunction,options.selected_index,
        options.degeneracy,options.filter);
    if(options.canonical_display_energies.size()==wavefunction.orbitals.size())
        for(std::size_t i=0;i<data.metadata.size();++i)
            data.metadata[i].energy_hartree=options.canonical_display_energies[i];
    if(options.routed &&
       options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction))
        for(std::size_t i=0;i<data.metadata.size();++i)
            if(const auto balance=routed_mo_shell_balance(*options.routed,i)){
                const bool occ=occupied(wavefunction.orbitals[i],
                                        options.filter.occupation_threshold);
                if(balance->complete&&occ&&balance->core>=0.70&&
                   balance->core>=balance->valence+0.10)
                    data.metadata[i].region=OrbitalRegion::Core;
                else if(occ||balance->valence>=0.05)
                    data.metadata[i].region=OrbitalRegion::Valence;
            }

    const LigandFieldEnvironment ligand_field=
        analyse_ligand_field_environment(wavefunction);
    const LigandScope ligand_scope=make_ligand_scope(
        wavefunction,ligand_field);
    auto current_centres=options.display_centre_atoms;
    if(current_centres.empty() && ligand_scope.available)current_centres.push_back(ligand_scope.metal);
    data.composition_scope=mo_composition_scope(wavefunction,options.routed,current_centres);
    if(options.routed)data.current_radial_shells=mo_current_radial_shells(
        wavefunction,*options.routed,current_centres);
    data.sigma_framework=analyse_mo_sigma_framework(wavefunction,options.nbo_source,
        options.routed,current_centres);
    if(options.routed && options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction))
        for(const auto& entry:options.routed->pi_couplings)if(entry.available() && entry.provider==RoutedProvider::Nbo) {
            const auto& channel=*entry.value;
            if(!current_centres.empty() && std::none_of(channel.centre_atoms.begin(),channel.centre_atoms.end(),
                [&](auto atom){return std::find(current_centres.begin(),current_centres.end(),atom)!=current_centres.end();}))continue;
            for(const auto& mode:channel.modes) {
                auto network=assess_pi_mode_network(channel,mode);
                if(network.verified)data.pi_mode_networks.push_back(std::move(network));
            }
        }
    if (ligand_field.available()) {
        data.ligand_field_point_group=ligand_field.local_point_group();
        data.ligand_field_geometry_id=ligand_field.geometry_machine_id();
        data.ligand_field_geometry_name=ligand_field.geometry_name();
        data.ligand_field_coordination_number=
            ligand_field.coordination_number();
        data.ligand_field_metal_atom=ligand_field.metal_atom;
        data.ligand_field_ligand_atoms=ligand_field.ligand_atoms;
        data.ligand_field_confidence=ligand_field.confidence;
        data.ligand_field_angular_rms=ligand_field.angular_rms;
        data.ligand_field_shape_measure=ligand_field.shape_measure;
        data.ligand_field_radial_cv=ligand_field.radial_cv;
    }
    attach_local_geometry_metadata(wavefunction,data);
    attach_electronic_state_metadata(wavefunction,data);

    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& orbital:wavefunction.orbitals) {
        data.annotations.push_back(annotate_orbital(orbital));
    }
    // Preserve the full family-level orientation evidence in every canonical
    // member annotation.  Canonical members inside a degenerate subspace are
    // deliberately not assigned arbitrary individual directions.
    for (const auto& assignment:wavefunction.delocalised_pi_assignments) {
        for (const auto orbital:assignment.orbitals) {
            if (orbital>=data.annotations.size()) continue;
            auto& pi=data.annotations[orbital].delocalised_pi;
            if (!pi.available || pi.family_id!=assignment.family_id) continue;
            pi.topology_available=
                assignment.provenance!=DataProvenance::Unavailable &&
                assignment.topology!=DelocalisedPiTopology::Unknown;
            pi.topology=assignment.topology;
            pi.orientation_channels=assignment.orientation_channels.size();
            pi.cyclic_topology=assignment.cyclic_topology;
            pi.orientation_channel_details=assignment.orientation_channels;
            pi.topology_graph=assignment.topology_graph;
        }
    }

    data.selection=build_valence_selection_plan(
        wavefunction,options,data.metadata);
    const auto labels=build_orbital_labels(
        wavefunction.orbitals,
        point_group_limited_degeneracy(wavefunction,options.degeneracy));
    std::vector<GroupCandidate> groups;
    groups.reserve(wavefunction.orbitals.size());
    for (std::size_t base=0;base<wavefunction.orbitals.size();) {
        groups.push_back(make_group_candidate(
            wavefunction,options,data.frontier,data.metadata,
            data.annotations,labels,base,
            ligand_scope.available?&ligand_scope:nullptr,0u,
            &data.current_radial_shells,&current_centres));
        const std::size_t step=base<labels.size()
            ?std::max<std::size_t>(1u,labels[base].group_size):1u;
        base+=step;
    }

    const bool active_space_mode=
        options.mode==MODiagramMode::DelocalisedPiFamilyOnly ||
        options.mode==MODiagramMode::MulticentreActiveSpaceOnly;
    const auto compact_family=active_space_mode
        ?compact_active_indices(wavefunction,options.mode)
        :std::set<std::size_t>{};
    if (active_space_mode) {
        for (auto& group:groups) {
            group.include=std::any_of(
                group.level.member_indices.begin(),
                group.level.member_indices.end(),
                [&](const auto member){
                    return compact_family.count(member)!=0u;
                });
            group.selected_by_reference=group.include;
            group.selected_by_raw=false;
            group.level.raw_data_fallback=false;
        }
    }

    std::optional<LocalAngularProjectionWorkspace> local_projection;
    if (ligand_field.available()) {
        local_projection.emplace(wavefunction,ligand_field.metal_atom,
            ligand_field.rotation_reference_to_input);
        recover_local_ligand_field_symmetry(
            *local_projection,ligand_field,groups);
        merge_local_pseudodegenerate_groups(
            wavefunction,ligand_field,groups);
        recover_local_ligand_field_symmetry(
            *local_projection,ligand_field,groups);
    }

    // Some FCHK producers omit member-level irreps even though the complete
    // exactly-degenerate subspace has an unambiguous ligand-field label.  A
    // recovered Eg/T2g/E/T2 label belongs to every canonical member of that
    // subspace; copy it back to the member metadata so per-line selection and
    // tooltips do not regress to N/A.
    for (const auto& group:groups) {
        const std::string recovered=group.level.metadata.symmetry_view.label;
        const std::string normalised=normalised_symmetry(recovered);
        if (normalised.empty() || normalised=="?" || normalised=="n/a") continue;
        for (const auto member:group.level.member_indices) {
            if (member>=data.metadata.size()) continue;
            const std::string current=normalised_symmetry(
                data.metadata[member].symmetry_view.label);
            if (current.empty() || current=="?" || current=="n/a") {
                data.metadata[member].symmetry_view=group.level.metadata.symmetry_view;
            }
        }
    }

    const SpinCollapseResult spin_collapse=collapse_spin_counterparts(
        wavefunction,options,groups);
    data.spin_counterpart_pair_count=spin_collapse.paired_groups;

    // Matching may recover a label from either spin channel.  Apply it to
    // both canonical members so the browser, hover text and exported member
    // metadata agree with the spatial-row diagram.
    for (const auto& group:groups) {
        const std::string recovered=group.level.metadata.symmetry_view.label;
        const std::string normalised=normalised_symmetry(recovered);
        if (normalised.empty() || normalised=="?" || normalised=="n/a") continue;
        for (const auto member:group.level.member_indices) {
            if (member<data.metadata.size()) {
                const std::string current=normalised_symmetry(
                    data.metadata[member].symmetry_view.label);
                if (current.empty() || current=="?" || current=="n/a") {
                    data.metadata[member].symmetry_view=group.level.metadata.symmetry_view;
                }
            }
        }
        for (const auto counterpart:group.level.member_spin_counterparts) {
            if (counterpart>=data.metadata.size()) continue;
            const std::string current=normalised_symmetry(
                data.metadata[counterpart].symmetry_view.label);
            if (current.empty() || current=="?" || current=="n/a") {
                const std::array<std::size_t,1> target{counterpart};
                data.metadata[counterpart].symmetry_view=candidate_orbital_symmetry(
                    group.level.metadata.symmetry_view,OrbitalSymmetryOrigin::SpinCounterpartCandidate,
                    target,std::numeric_limits<double>::quiet_NaN());
            }
        }
    }

    if (active_space_mode) {
        // Degeneracy grouping is a display convenience, never permission to
        // enlarge a chemically selected active space. Rebuild every retained
        // row from exact family members after local-irrep recovery and spin
        // collapse. This preserves valid grouping inside the family while an
        // accidentally near-degenerate outsider cannot leak to screen or
        // PNG/SVG/JSON/CSV exports.
        const auto no_counterpart=std::numeric_limits<std::size_t>::max();
        for (auto& group:groups) {
            const GroupCandidate original=group;
            std::vector<std::size_t> members;
            for (const auto member:original.level.member_indices) {
                if (compact_family.count(member)!=0u) members.push_back(member);
            }
            if (members.empty()) {
                group.include=false;
                continue;
            }

            GroupCandidate exact;
            bool initialised=false;
            for (const auto member:members) {
                auto singleton=make_group_candidate(
                    wavefunction,options,data.frontier,data.metadata,
                    data.annotations,labels,member,
                    ligand_scope.available?&ligand_scope:nullptr,1u,
                    &data.current_radial_shells,&current_centres);
                if (!initialised) {
                    exact=std::move(singleton);
                    initialised=true;
                } else {
                    exact=merge_local_groups(wavefunction,exact,singleton);
                }
            }

            exact.level.member_electrons.clear();
            exact.level.member_spin_counterparts.clear();
            exact.level.total_occupation=0.0;
            for (const auto member:members) {
                const auto found=std::find(
                    original.level.member_indices.begin(),
                    original.level.member_indices.end(),member);
                const auto position=static_cast<std::size_t>(std::distance(
                    original.level.member_indices.begin(),found));
                const auto counterpart=
                    found!=original.level.member_indices.end() &&
                    position<original.level.member_spin_counterparts.size()
                        ?original.level.member_spin_counterparts[position]
                        :no_counterpart;
                if (counterpart!=no_counterpart &&
                    compact_family.count(counterpart)!=0u &&
                    position<original.level.member_electrons.size()) {
                    exact.level.member_electrons.push_back(
                        original.level.member_electrons[position]);
                    exact.level.member_spin_counterparts.push_back(counterpart);
                    exact.level.total_occupation+=
                        static_cast<double>(wavefunction.orbitals[member].occupation)+
                        static_cast<double>(wavefunction.orbitals[counterpart].occupation);
                } else {
                    exact.level.member_electrons.push_back(
                        electron_glyphs_for_orbital(
                            wavefunction.orbitals[member],
                            data.frontier.separate_spin_sets));
                    exact.level.member_spin_counterparts.push_back(
                        no_counterpart);
                    exact.level.total_occupation+=static_cast<double>(
                        wavefunction.orbitals[member].occupation);
                }
            }
            exact.level.metadata.occupation=static_cast<float>(
                exact.level.total_occupation/static_cast<double>(
                    std::max<std::size_t>(1u,members.size())));
            if (!exact.level.member_electrons.empty()) {
                exact.level.electrons=exact.level.member_electrons.front();
            }
            // Rebuilding a subset invalidates the original subspace proof.
            // A complete identical target can retain it; otherwise project the
            // actual retained target again using the shared immutable workspace.
            exact.level.metadata.symmetry_view=aggregate_molecular_symmetry(wavefunction,members);
            if (members==original.level.member_indices) {
                exact.level.metadata.symmetry_view=original.level.metadata.symmetry_view;
            } else if (local_projection) {
                const int family=dominant_metal_family(exact.level);
                if(family>=0) exact.level.metadata.symmetry_view=evaluate_local_orbital_symmetry(
                    *local_projection,ligand_field,members,family);
            }
            exact.include=original.include;
            exact.selected_by_reference=original.selected_by_reference;
            exact.selected_by_raw=false;
            exact.suppressed_spin_counterpart=
                original.suppressed_spin_counterpart;
            exact.locally_grouped=members.size()>1u &&
                original.locally_grouped;
            exact.local_irrep_copy=exact.level.metadata.symmetry_view.local_assignment
                ?exact.level.metadata.symmetry_view.local_assignment->copy_index:0u;
            exact.locally_classified=orbital_symmetry_is_local(exact.level.metadata.symmetry_view);
            group=std::move(exact);
        }
    }

    // Local regrouping/spin collapse can change member sets. Recompute the
    // same ledger over the final immutable canonical source membership.
    for(auto& group:groups) {
        auto composition_members=group.level.member_indices;
        for(const auto counterpart:group.level.member_spin_counterparts)
            if(counterpart<wavefunction.orbitals.size() &&
               std::find(composition_members.begin(),composition_members.end(),counterpart)==composition_members.end())
                composition_members.push_back(counterpart);
        group.level.composition=mo_group_composition_ledger(wavefunction,options.routed,
            composition_members,current_centres,data.current_radial_shells);
        if(group.level.composition.complete) {
            group.level.metal_s_weight=group.level.composition.centre_current_s;
            group.level.metal_p_weight=group.level.composition.centre_current_p;
            group.level.metal_d_weight=group.level.composition.centre_current_d;
            group.level.ligand_p_weight=group.level.composition.ligand_valence_p;
        }
    }
    std::vector<RawPiPair> raw_pairs;
    if(!options.canonical_display_energies.empty())
        std::stable_sort(groups.begin(),groups.end(),[&](const auto& a,const auto& b) {
            const auto as=group_spin(wavefunction,a),bs=group_spin(wavefunction,b);
            if(as!=bs)return as==Spin::Alpha;
            return a.level.layout_energy_hartree<b.level.layout_energy_hartree;
        });
    if (!active_space_mode) {
        raw_pairs=find_pi_pairs(
            wavefunction,groups,options,data.ligand_field_point_group,
            ligand_scope.available?&ligand_scope:nullptr,data.pi_partner_candidates);
    } else {
        // Pair detection may enrich labels and, in the ordinary valence view,
        // deliberately brings both sides of a ligand-field interaction into
        // view.  A requested family/active-space diagram has a stricter
        // membership contract: run pairing on a copy, then retain only pairs
        // whose two rows were already members of that exact active space.
        // Thus neither screen nor export can leak an unrelated MO into a
        // pi-only or multicentre-only diagram.
        auto paired_groups=groups;
        const auto candidates=find_pi_pairs(
            wavefunction,paired_groups,options,data.ligand_field_point_group,
            ligand_scope.available?&ligand_scope:nullptr,data.pi_partner_candidates);
        for (const auto& pair:candidates) {
            if (pair.lower>=groups.size() || pair.upper>=groups.size() ||
                !groups[pair.lower].include || !groups[pair.upper].include) {
                continue;
            }
            groups[pair.lower]=std::move(paired_groups[pair.lower]);
            groups[pair.upper]=std::move(paired_groups[pair.upper]);
            groups[pair.lower].include=true;
            groups[pair.upper].include=true;
            raw_pairs.push_back(pair);
        }
    }
    for (const auto& pair:raw_pairs) {
        for (const auto group_index:{pair.lower,pair.upper}) {
            if (group_index>=groups.size()) continue;
            const auto& group=groups[group_index].level;
            const std::string recovered=group.metadata.symmetry_view.label;
            const std::string normalised=normalised_symmetry(recovered);
            if (normalised.empty() || normalised=="?" ||
                normalised=="n/a") continue;
            for (const auto member:group.member_indices) {
                if (member>=data.metadata.size()) continue;
                const std::string current=normalised_symmetry(
                    data.metadata[member].symmetry_view.label);
                if (current.empty() || current=="?" || current=="n/a") {
                    data.metadata[member].symmetry_view=group.metadata.symmetry_view;
                }
            }
        }
    }

    // Keep the reduced diagram readable.  The budget applies to degenerate
    // rows rather than canonical MOs.  In the compact ligand-field view the
    // canonical row set is purely chemical and therefore independent of the
    // MO currently inspected in the browser or 3-D viewport.  Both members
    // of a resolved donor/acceptor pair are never removed; the selected row
    // is pinned only in the expanded view.
    std::size_t row_budget=options.max_levels>0u
        ?std::max<std::size_t>(4u,options.max_levels)
        :std::clamp<std::size_t>(2u*options.neighbourhood+1u,10u,48u);
    std::vector<bool> essential(groups.size(),false);
    for (std::size_t group=0;group<groups.size();++group) {
        essential[group]=
            (options.mode==MODiagramMode::DelocalisedPiFamilyOnly ||
             options.mode==MODiagramMode::MulticentreActiveSpaceOnly)
                ?groups[group].include
                :(!options.hide_ligand_centred_intermediates &&
                  groups[group].level.metadata.selected);
        for (const auto member:groups[group].level.member_indices) {
            if (member>=wavefunction.orbitals.size()) continue;
            const auto& orbital=wavefunction.orbitals[member];
            if (!options.hide_ligand_centred_intermediates &&
                occupied(orbital,options.filter.occupation_threshold) &&
                !confidently_deep_core_orbital(orbital,options.filter)) {
                essential[group]=true;
            }
            if (wavefunction.charge>0 &&
                wavefunction.charge_provenance!=DataProvenance::Unavailable &&
                ((!data.frontier.separate_spin_sets && data.frontier.alpha_lumo &&
                  *data.frontier.alpha_lumo==member) ||
                 (data.frontier.separate_spin_sets &&
                  ((data.frontier.alpha_lumo &&
                    *data.frontier.alpha_lumo==member) ||
                   (data.frontier.beta_lumo &&
                    *data.frontier.beta_lumo==member))))) {
                essential[group]=true;
            }
        }
    }
    for (const auto& pair:raw_pairs) {
        if(options.aomo_scope==1 && pair.evidence &&
           !pair.evidence->channel.channel_id.empty() &&
           !pair.evidence->channel.mode_assessment.ordinary_display_eligible)continue;
        essential[pair.retained]=true;
        if (pair.gap_kind==OrbitalEnergyGapKind::CrystalField ||
            pair.kind!=PiInteractionKind::WeakNearNonbonding) {
            essential[pair.lower]=true;
            essential[pair.upper]=true;
        }
    }
    std::size_t hidden_intermediate_count=0u;
    if (options.hide_ligand_centred_intermediates) {
        const std::string& local_group=data.ligand_field_point_group;
        std::map<std::string,std::size_t> d_irrep_multiplicity;
        if (const auto decomposition=decompose_metal_ao_shell(
                local_group,MetalAOShell::D)) {
            for (const auto& copy:*decomposition) {
                ++d_irrep_multiplicity[normalised_symmetry(
                    std::string(copy.label))];
            }
        }
        const bool supported_ligand_field=!d_irrep_multiplicity.empty();
        std::map<std::string,std::vector<std::size_t>> candidates;
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            const auto& explanation=groups[group].level.metadata.symmetry_view;
            if (!orbital_symmetry_is_local(explanation) || explanation.point_group!=local_group)
                continue;
            const std::string symmetry=normalised_symmetry(
                groups[group].level.metadata.symmetry_view.label);
            if (d_irrep_multiplicity.count(symmetry)) {
                candidates[symmetry].push_back(group);
            }
        }
        std::set<std::size_t> anchors;
        for (auto& [symmetry,indices]:candidates) {
            std::sort(indices.begin(),indices.end(),[&](const auto a,const auto b) {
                const auto& left=groups[a].level;
                const auto& right=groups[b].level;
                const double left_score=left.metal_d_weight+
                    0.25*left.direct_ligand_p_weight;
                const double right_score=right.metal_d_weight+
                    0.25*right.direct_ligand_p_weight;
                if (left_score!=right_score) return left_score>right_score;
                return a<b;
            });
            const std::size_t needed=std::max<std::size_t>(
                1u,d_irrep_multiplicity[symmetry]);
            for (std::size_t i=0;i<std::min(needed,indices.size());++i) {
                anchors.insert(indices[i]);
            }
        }
        for (const auto group:anchors) essential[group]=true;

        // Preserve one bonding and one antibonding representative for each
        // central-metal s/p/d sigma-framework symmetry.  The measured sigma
        // channel is the gate; the point-group catalogue supplies the allowed
        // local labels without hard-coding Td/Oh names.
        using SigmaKey=std::pair<std::string,bool>;
        std::map<SigmaKey,std::pair<std::size_t,double>> sigma_representatives;
        std::set<std::string> local_spd_labels;
        if (const auto decomposition=decompose_metal_spd(local_group)) {
            for (const auto* block:{&decomposition->s,&decomposition->p,
                                   &decomposition->d}) {
                for (const auto& copy:*block) {
                    local_spd_labels.insert(normalised_symmetry(
                        std::string(copy.label)));
                }
            }
        }
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            const auto& level=groups[group].level;
            const double metal_sp=level.metal_s_weight+level.metal_p_weight;
            const double metal_spd=metal_sp+level.metal_d_weight;
            const std::string symmetry=normalised_symmetry(
                level.metadata.symmetry_view.label);
            if (symmetry.empty() || symmetry=="?" || symmetry=="n/a") continue;
            const bool local_sigma_label=orbital_symmetry_is_local(level.metadata.symmetry_view) &&
                level.metadata.symmetry_view.point_group==local_group && local_spd_labels.count(symmetry)>0u;
            const double sigma_floor=local_sigma_label?0.50:0.55;
            if (level.sigma_fraction<sigma_floor || metal_spd<0.025) continue;
            const SigmaKey key{
                symmetry,
                level.metal_ligand_overlap<0.0};
            const double balanced_mixing=2.0*std::min(
                metal_spd,level.direct_ligand_p_weight);
            const double representative_score=balanced_mixing+
                0.20*std::min(1.0,std::abs(level.metal_ligand_overlap))+
                0.05*metal_spd;
            const auto current=sigma_representatives.find(key);
            if (current==sigma_representatives.end() ||
                representative_score>current->second.second) {
                sigma_representatives[key]={group,representative_score};
            }
        }
        for (const auto& [key,representative]:sigma_representatives) {
            (void)key;
            essential[representative.first]=true;
        }

        // One recovered row for every copy in the formal d decomposition is a
        // complete minimal ligand-field manifold.  Repeated irreps (for
        // example 2A1 in C2v or 2E in C3v) are counted explicitly.
        const bool complete_d_manifold=supported_ligand_field &&
            std::all_of(d_irrep_multiplicity.begin(),
                d_irrep_multiplicity.end(),
                [&](const auto& expected) {
                    const auto found=candidates.find(expected.first);
                    return found!=candidates.end() &&
                           found->second.size()>=expected.second;
                });
        if (complete_d_manifold) {
            for (std::size_t group=0;group<groups.size();++group) {
                if (!groups[group].include) continue;
                if (groups[group].level.metal_d_weight>=0.15) {
                    essential[group]=true;
                }
                if (!essential[group]) {
                    groups[group].include=false;
                    ++hidden_intermediate_count;
                }
            }
        }
    }
    std::size_t included_count=static_cast<std::size_t>(std::count_if(
        groups.begin(),groups.end(),[](const auto& group) {
            return group.include;
        }));
    if (included_count>row_budget) {
        const double frontier_energy=data.frontier.homo && data.frontier.lumo
            ?0.5*(diagram_orbital_energy(wavefunction,options,*data.frontier.homo)+
                  diagram_orbital_energy(wavefunction,options,*data.frontier.lumo))
            :0.0;
        struct RankedGroup { std::size_t index=0; double score=0.0; };
        std::vector<RankedGroup> ranked;
        std::size_t essential_count=0u;
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            if (essential[group]) {
                ++essential_count;
                continue;
            }
            const auto& level=groups[group].level;
            const double metal=level.metal_s_weight+level.metal_p_weight+
                               level.metal_d_weight;
            const double score=4.0*metal+1.5*level.ligand_p_weight+
                0.35*std::max(level.sigma_fraction,level.pi_fraction)+
                1.5*std::min(0.20,std::abs(level.metal_ligand_overlap))-
                0.20*std::abs(level.layout_energy_hartree-frontier_energy)+
                (level.total_occupation>
                     options.filter.occupation_threshold?0.35:0.0);
            ranked.push_back({group,score});
            groups[group].include=false;
        }
        std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
            if (a.score!=b.score) return a.score>b.score;
            return a.index<b.index;
        });
        const std::size_t available=row_budget>essential_count
            ?row_budget-essential_count:0u;
        data.selection.protected_overflow_count=
            essential_count>row_budget?essential_count-row_budget:0u;
        for (std::size_t i=0;i<std::min(available,ranked.size());++i) {
            groups[ranked[i].index].include=true;
        }
    }

    // AOMO presets own the central population as well as the side population.
    // This final membership pass does not change source coefficients, grouping,
    // occupations, or the scientific classification of a folded background.
    if(options.aomo_scope && !active_space_mode) {
        const bool full=options.aomo_scope==3;
        const bool research=options.aomo_scope==2;
        std::set<std::size_t> protected_pi;
        for(const auto& pair:raw_pairs)
            if(pair.evidence && pair.evidence->accepted &&
               pair.evidence->channel.mode_assessment.ordinary_display_eligible) {
                protected_pi.insert(pair.lower);protected_pi.insert(pair.upper);
            }
        for(const auto& network:data.pi_mode_networks)if(network.ordinary_display_eligible)
            for(const auto& node:network.nodes)if(node.primary)
                for(std::size_t index=0;index<groups.size();++index)
                    if(std::all_of(node.members.begin(),node.members.end(),[&](auto member) {
                        return mo_diagram_level_covers_orbital(groups[index].level,member);
                    }))protected_pi.insert(index);
        const auto centre_coverage=[](const MOGroupCompositionLedger& c) {
            return c.centre_current_s+c.centre_current_p+c.centre_current_d+c.centre_current_f;
        };
        struct SigmaSupport {std::size_t index=0,members=0;double trace=0;bool deep=false;};
        std::vector<SigmaSupport> sigma_support;
        if(data.sigma_framework.available)for(std::size_t index=0;index<groups.size();++index) {
            auto& group=groups[index];if(group.suppressed_spin_counterpart)continue;
            std::set<std::size_t> members(group.level.member_indices.begin(),group.level.member_indices.end());
            for(const auto counterpart:group.level.member_spin_counterparts)
                if(counterpart<wavefunction.orbitals.size())members.insert(counterpart);
            double trace=0;
            for(const auto member:members)if(member<data.sigma_framework.canonical_weights.size())
                trace+=data.sigma_framework.canonical_weights[member];
            group.level.display_decision.sigma_coverage=members.empty()?0:trace/members.size();
            if(group.level.total_occupation<=options.filter.occupation_threshold ||
               !group.level.composition.complete || group.level.composition.core>=.70)continue;
            const auto& c=group.level.composition;
            const double frontier=data.frontier.homo?diagram_orbital_energy(wavefunction,options,*data.frontier.homo):
                group.level.layout_energy_hartree;
            const bool deep=c.ligand_valence_s>=.80 && centre_coverage(c)<.20 &&
                group.level.layout_energy_hartree<frontier-.20;
            sigma_support.push_back({index,members.size(),trace,deep});
        }
        // The energy reference is the lowest empty state with measured
        // current-shell support. An external diffuse global LUMO is still
        // reachable in the browser/research view, but cannot set this window.
        double effective_lumo=std::numeric_limits<double>::infinity();
        for(const auto& group:groups)if(!group.suppressed_spin_counterpart &&
            group.level.total_occupation<=options.filter.occupation_threshold &&
            group.level.composition.complete && centre_coverage(group.level.composition)>=0.05)
            effective_lumo=std::min(effective_lumo,group.level.layout_energy_hartree);
        const double virtual_window=std::min(0.75,options.filter.virtual_window_hartree);
        std::set<std::size_t> shell_recovery;
        if(!full && !research && !current_centres.empty()) {
            // Recover the low-energy current-shell trace, not the entire
            // arbitrarily large virtual basis. The energy kernel defines a
            // reading budget, not a bond energy or numerical-zero threshold.
            // Each angular family has at most two complete radial manifolds'
            // worth of canonical members, plus separately verified major pairs.
            for(int l=0;l<=3;++l) {
                struct Support {std::size_t index=0,members=0;double score=0;};
                std::vector<Support> ranked;double target=0;
                for(std::size_t i=0;i<groups.size();++i) {
                    const auto& group=groups[i];const auto& level=group.level;
                    if(group.suppressed_spin_counterpart || !level.composition.complete ||
                       level.total_occupation>options.filter.occupation_threshold ||
                       level.composition.core>=0.70 || !std::isfinite(effective_lumo))continue;
                    const auto& c=level.composition;
                    const double weight=l==0?c.centre_current_s:l==1?c.centre_current_p:
                        l==2?c.centre_current_d:c.centre_current_f;
                    const double gap=std::max(0.0,level.layout_energy_hartree-effective_lumo);
                    const double score=weight*level.member_indices.size()/
                        (1+std::pow(gap/std::max(0.05,virtual_window),2));
                    target+=score;
                    if(centre_coverage(c)>=0.05 && weight>0)ranked.push_back(
                        {i,level.member_indices.size(),score});
                }
                std::stable_sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
                    return a.score>b.score;
                });
                double recovered=0;std::size_t members=0;
                const std::size_t member_budget=2u*(2u*static_cast<std::size_t>(l)+1u)*
                    current_centres.size();
                for(const auto& support:ranked) {
                    if(recovered>=0.90*target)break;
                    if(members && members+support.members>member_budget)continue;
                    shell_recovery.insert(support.index);
                    recovered+=support.score;members+=support.members;
                }
            }
        }
        for(std::size_t index=0;index<groups.size();++index) {
            auto& candidate=groups[index];
            if(candidate.suppressed_spin_counterpart)continue;
            const auto& members=candidate.level.member_indices;
            const auto& composition=candidate.level.composition;
            auto& decision=candidate.level.display_decision;
            decision.coverage=centre_coverage(composition);
            decision.major_relation=protected_pi.contains(index);
            const bool is_occupied=candidate.level.total_occupation>options.filter.occupation_threshold;
            decision.energy_window=is_occupied || candidate.level.layout_energy_hartree<=effective_lumo+virtual_window;
            decision.frontier=candidate.level.homo || (candidate.level.lumo &&
                (!composition.complete || decision.coverage>=0.05));
            if(!is_occupied && composition.complete && decision.coverage>=0.05 &&
               candidate.level.layout_energy_hartree==effective_lumo) {
                decision.frontier=true;
                decision.reason_codes.push_back("effective-current-space-lowest-empty-state");
            }
            const bool core_group=composition.complete?composition.core>=0.70:
                std::all_of(members.begin(),members.end(),[&](auto member) {
                    return member<data.metadata.size() && data.metadata[member].region==OrbitalRegion::Core;
                });
            if(full || research) {
                candidate.include=full || options.show_core_background || !core_group;
                decision.reason_codes.push_back(full?"all-source-groups":"research-source-groups");
            } else if(core_group) {
                candidate.include=options.show_core_background;
                decision.reason_codes.push_back(candidate.include?"requested-core-background":"folded-core-background");
            } else if(options.show_fragment_background) {
                candidate.include=true;decision.reason_codes.push_back("requested-fragment-background");
            } else if(composition.complete && !current_centres.empty()) {
                const double frontier=data.frontier.homo?diagram_orbital_energy(wavefunction,options,*data.frontier.homo):
                    candidate.level.layout_energy_hartree;
                const bool deep_background=composition.ligand_valence_s>=0.80 && decision.coverage<0.20 &&
                    candidate.level.layout_energy_hartree<frontier-0.20;
                candidate.include=decision.frontier || decision.major_relation ||
                    (!deep_background && (is_occupied?decision.coverage>=0.05:shell_recovery.contains(index)));
                if(decision.frontier)decision.reason_codes.push_back("retained-supported-frontier");
                if(decision.major_relation)decision.reason_codes.push_back("retained-verified-major-mode-relation");
                if(shell_recovery.contains(index) && candidate.include)
                    decision.reason_codes.push_back("retained-current-shell-low-energy-trace-budget");
                if(is_occupied && candidate.include && !deep_background)
                    decision.reason_codes.push_back("retained-occupied-current-shell-support");
                if(!candidate.include) {
                    ++hidden_intermediate_count;
                    decision.reason_codes.push_back(deep_background?"folded-deep-ligand-valence-s-background":
                        decision.coverage<0.05?"folded-outside-current-shell-space":"folded-current-shell-display-budget");
                }
                if(candidate.level.lumo && !decision.frontier)
                    decision.reason_codes.push_back("global-LUMO-outside-current-space-browser-reachable");
            } else decision.reason_codes.push_back(composition.complete?
                "legacy-membership-no-centre-scope":"legacy-membership-complete-NAO-unavailable");
            decision.included=candidate.include;
        }
        // Recover a fixed occupied donor projector, whose NBO columns and
        // canonical closure were validated before any view pruning. This is
        // an overlapping channel reading, never a redefined MO population.
        // Existing displayed groups count towards the same finite budget.
        auto& sigma=data.sigma_framework;
        sigma.display_member_budget=2*sigma.rank;
        for(const auto& support:sigma_support)
            if(groups[support.index].include && support.trace>0) {
                sigma.retained_occupied_trace+=support.trace;
                if(groups[support.index].level.display_decision.sigma_coverage>=.05)
                    sigma.retained_members+=support.members;
            }
        if(!full && !research && sigma.available) {
            std::stable_sort(sigma_support.begin(),sigma_support.end(),[](const auto& a,const auto& b) {
                return a.trace>b.trace;
            });
            for(const auto& support:sigma_support) {
                if(sigma.retained_occupied_trace>=sigma.trace_target*sigma.occupied_trace)break;
                auto& group=groups[support.index];
                if(group.include || support.deep || group.level.display_decision.sigma_coverage<.05 ||
                   sigma.retained_members+support.members>sigma.display_member_budget)continue;
                group.include=true;group.level.display_decision.included=true;
                auto& reasons=group.level.display_decision.reason_codes;
                std::erase_if(reasons,[](const auto& reason){return reason.rfind("folded-",0)==0;});
                reasons.push_back("retained-verified-occupied-sigma-donor-space-budget");
                sigma.retained_occupied_trace+=support.trace;sigma.retained_members+=support.members;
                if(hidden_intermediate_count)--hidden_intermediate_count;
            }
        }
    }

    // Freeze the scientific source domain before weak-display decisions. A
    // folded relation cannot shrink its own validation domain and pass itself.
    std::vector<std::size_t> stable_display_members;
    for(const auto& group:groups)if(group.include && !group.suppressed_spin_counterpart)
        stable_display_members.insert(stable_display_members.end(),group.level.member_indices.begin(),group.level.member_indices.end());
    std::map<std::string,PiFrozenOperatorAssessment> frozen_display_scopes;
    if(options.routed)for(auto& pair:raw_pairs)if(pair.evidence && !pair.evidence->channel.channel_id.empty()) {
        const auto& id=pair.evidence->channel.channel_id;
        const NboPiCoupling* channel=nullptr;
        for(const auto& entry:options.routed->pi_couplings)if(entry.available() && entry.value->id==id){channel=&*entry.value;break;}
        if(!channel)continue;
        const auto& mode=pair.evidence->channel.mode_assessment;
        const double lower_strength=std::abs(mode.lower_cross_fock_mean_hartree)*
            pair.evidence->channel.lower_members.size();
        const double upper_strength=std::abs(mode.upper_cross_fock_mean_hartree)*
            pair.evidence->channel.upper_members.size();
        pair.display_strength_hartree=std::min(lower_strength,upper_strength);
        if(!frozen_display_scopes.contains(id)) {
            std::vector<std::size_t> same_spin_members;
            for(auto member:stable_display_members)if(member<wavefunction.orbitals.size() &&
                (channel->spin==NboSpin::Beta)==(wavefunction.orbitals[member].spin==Spin::Beta))same_spin_members.push_back(member);
            frozen_display_scopes[id]=assess_pi_frozen_display_scope(*channel,same_spin_members);
        }
        auto result=*pair.evidence;
        result.channel.frozen_operator=frozen_display_scopes.at(id);
        result.channel.display_calibration=pi_channel_display_calibration(*channel);
        result.weak=pi_display_negligible(result.channel.frozen_operator,result.channel.display_calibration);
        pair.evidence=std::make_shared<const PiPartnerAssessment>(result);
        for(auto& candidate:data.pi_partner_candidates)if(candidate.accepted && candidate.channel.channel_id==id &&
            candidate.channel.lower_members==result.channel.lower_members && candidate.channel.upper_members==result.channel.upper_members)
            candidate=result;
    }

    std::vector<std::size_t> group_to_level(
        groups.size(),std::numeric_limits<std::size_t>::max());
    std::vector<double> axis_energies;
    for (std::size_t group=0;group<groups.size();++group) {
        auto& level=groups[group].level;
        level.display_decision.included=groups[group].include && !groups[group].suppressed_spin_counterpart;
        if(groups[group].suppressed_spin_counterpart)
            level.display_decision.reason_codes.push_back("paired-opposite-spin-source-group");
        if(level.display_decision.reason_codes.empty())
            level.display_decision.reason_codes.push_back(level.display_decision.included?
                "retained-legacy-or-active-space-membership":"folded-legacy-or-active-space-membership");
        data.group_audit.push_back({level.member_indices,level.member_spin_counterparts,
            level.layout_energy_hartree,level.total_occupation,level.composition,level.display_decision});
        if (!groups[group].include) continue;
        group_to_level[group]=data.levels.size();
        axis_energies.push_back(groups[group].level.layout_energy_hartree);
        data.levels.push_back(groups[group].level);
    }

    for (const auto& pair:raw_pairs) {
        const std::size_t lower=group_to_level[pair.lower];
        const std::size_t upper=group_to_level[pair.upper];
        const std::size_t retained=group_to_level[pair.retained];
        if (retained==std::numeric_limits<std::size_t>::max()) continue;
        PiInteractionDescriptor descriptor;
        descriptor.lower_level=lower==std::numeric_limits<std::size_t>::max()
            ?retained:lower;
        descriptor.upper_level=upper==std::numeric_limits<std::size_t>::max()
            ?retained:upper;
        descriptor.lower_orbitals=groups[pair.lower].level.member_indices;
        descriptor.upper_orbitals=groups[pair.upper].level.member_indices;
        descriptor.symmetry=pair.symmetry.empty()&&(!pair.evidence||pair.evidence->channel.channel_id.empty())
            ?groups[pair.lower].level.metadata.symmetry_view.label:pair.symmetry;
        descriptor.kind=pair.kind;
        descriptor.orbital_evidence=pair.evidence;
        descriptor.equivalent_channel_ids=pair.equivalent_channel_ids;
        descriptor.crystal_field_evidence=pair.crystal_field_evidence;
        descriptor.gap_kind=pair.gap_kind;
        descriptor.lower_symmetry_scope=groups[pair.lower].level.metadata.symmetry_view;
        descriptor.upper_symmetry_scope=groups[pair.upper].level.metadata.symmetry_view;
        descriptor.lower_energy_hartree=groups[pair.lower].level.layout_energy_hartree;
        descriptor.upper_energy_hartree=groups[pair.upper].level.layout_energy_hartree;
        descriptor.lower_energy_spread_hartree=groups[pair.lower].level.energy_spread_hartree;
        descriptor.upper_energy_spread_hartree=groups[pair.upper].level.energy_spread_hartree;
        descriptor.weak_split_threshold_hartree=pair.gap_kind==OrbitalEnergyGapKind::CrystalField?
            options.weak_crystal_field_split_hartree:options.weak_pi_split_hartree;
        descriptor.weak_overlap_threshold=pair.gap_kind==OrbitalEnergyGapKind::CrystalField?
            options.weak_crystal_field_overlap:options.weak_metal_ligand_overlap;
        descriptor.splitting_hartree=pair.split;
        descriptor.confidence=pair.confidence;
        descriptor.display_strength_hartree=pair.display_strength_hartree;
        descriptor.lower_visible=lower!=std::numeric_limits<std::size_t>::max();
        descriptor.upper_visible=upper!=std::numeric_limits<std::size_t>::max();
        descriptor.retained_level=retained;
        if (descriptor.gap_kind==OrbitalEnergyGapKind::CrystalField)
            data.crystal_field_gaps.push_back(std::move(descriptor));
        else data.pi_interactions.push_back(std::move(descriptor));
    }

    std::size_t raw_groups=0u;
    for (const auto& level:data.levels) {
        if (level.raw_data_fallback) ++raw_groups;
    }
    data.spin_counterpart_unmatched_visible=static_cast<std::size_t>(
        std::count_if(groups.begin(),groups.end(),[&](const auto& group) {
            return group.include && !group.suppressed_spin_counterpart &&
                   group_spin(wavefunction,group)==Spin::Beta;
        }));
    data.spin_counterparts_collapsed=
        data.spin_counterpart_pair_count>0u &&
        data.spin_counterpart_unmatched_visible==0u;
    data.spin_counterparts_partial=
        data.spin_counterpart_pair_count>0u &&
        data.spin_counterpart_unmatched_visible>0u;
    std::ostringstream summary;
    summary<<(options.mode==MODiagramMode::DelocalisedPiFamilyOnly
              ?"delocalised pi family groups: "
              :(options.mode==MODiagramMode::MulticentreActiveSpaceOnly
                  ?"multicentre active-space groups: "
                  :"ligand-field valence groups: "))<<data.levels.size()
           <<"; local field="<<(data.ligand_field_point_group.empty()
                 ?"unresolved":data.ligand_field_point_group)
           <<"; geometry="<<(data.ligand_field_geometry_id.empty()
                ?"unresolved":data.ligand_field_geometry_id)
           <<"; CN="<<data.ligand_field_coordination_number
           <<"; spin counterparts="<<(data.spin_counterparts_collapsed
                ?"collapsed":(data.spin_counterparts_partial
                    ?"partial":"separate"))
           <<" (pairs="<<data.spin_counterpart_pair_count
           <<", unmatched-visible="<<data.spin_counterpart_unmatched_visible
           <<')'
           <<"; pi pairs="<<data.pi_interactions.size()
           <<"; crystal-field gaps="<<data.crystal_field_gaps.size()
           <<"; protected row overflow="
           <<data.selection.protected_overflow_count
           <<"; raw-MO recovered groups="<<raw_groups
           <<"; intermediate groups hidden="<<hidden_intermediate_count
           <<"; degeneracies collapsed";
    data.selection.summary=summary.str();

    data.energy_transform=build_energy_transform(
        axis_energies,options.energy_axis_mode,
        options.nonlinear_minimum_gap_weight);
    for (auto& level:data.levels) {
        level.metadata.molecular_member_symmetries.clear();
        auto members=level.member_indices;
        for (const auto counterpart:level.member_spin_counterparts) {
            if (counterpart<wavefunction.orbitals.size() &&
                std::find(members.begin(),members.end(),counterpart)==members.end())
                members.push_back(counterpart);
        }
        for (const auto index:members)
            level.metadata.molecular_member_symmetries.push_back(molecular_orbital_symmetry(wavefunction,index));
    }
    return data;
}

MODiagramData build_mo_diagram_data(const Wavefunction& wavefunction,
                                   const MODiagramOptions& options) {
    auto prepared=options;
    NboRoCommonEnergyModel common;
    if(options.ro_common_energy)common=*options.ro_common_energy;
    else if(options.nbo_source && wavefunction.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared) {
        const auto raw=options.source_salc_model?*options.source_salc_model:
            build_nbo_salc_model(wavefunction,*options.nbo_source);
        common=build_nbo_ro_common_energy(wavefunction,*options.nbo_source,raw);
    }
    const bool common_view=options.use_ro_common_energy && common.available &&
        common.canonical_fingerprint==nbo_canonical_fingerprint(wavefunction) &&
        common.orbitals.size()==wavefunction.orbitals.size() && chemistry_available(wavefunction);
    prepared.canonical_display_energies.clear();
    if(common_view) {
        prepared.canonical_display_energies.resize(wavefunction.orbitals.size());
        for(const auto& entry:common.orbitals)
            prepared.canonical_display_energies.at(entry.canonical_index)=entry.common_energy_hartree.value();
    }
    auto data=build_mo_diagram_data_impl(wavefunction,prepared);
    data.ro_common_energy=std::move(common);data.using_ro_common_energy=common_view;
    if(options.pi_field_response)data.pi_field_response=*options.pi_field_response;
    else if(options.nbo_source && options.routed &&
            options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction)) {
        std::vector<const NboPiCoupling*> channels;
        for(const auto& entry:options.routed->pi_couplings)
            if(entry.available() && entry.provider==RoutedProvider::Nbo)channels.push_back(&*entry.value);
        if(!channels.empty())data.pi_field_response=analyse_pi_field_responses(
            wavefunction,*options.nbo_source,channels);
    }
    // Counts describe the complete final scrollable central graph, not the
    // earlier valence selector or the currently clipped canvas pixels.
    std::set<std::size_t> final_members;
    for(const auto& level:data.levels) {
        for(const auto member:level.member_indices)if(member<wavefunction.orbitals.size())final_members.insert(member);
        if(level.member_indices.empty() && level.metadata.orbital_index<wavefunction.orbitals.size())
            final_members.insert(level.metadata.orbital_index);
        for(const auto member:level.member_spin_counterparts)
            if(member<wavefunction.orbitals.size())final_members.insert(member);
    }
    data.selection.included_indices.assign(final_members.begin(),final_members.end());
    data.selection.final_group_count=data.levels.size();
    data.selection.final_member_count=final_members.size();
    data.selection.counts_are_final=true;
    data.selection.valence_occupied_count=0;
    data.selection.frontier_virtual_count=0;
    for(const auto member:final_members)
        if(occupied(wavefunction.orbitals[member],options.filter.occupation_threshold))
            ++data.selection.valence_occupied_count;
        else ++data.selection.frontier_virtual_count;
    data.selection.hidden_count=wavefunction.orbitals.size()-final_members.size();
    const bool routed=options.routed &&
        options.routed->canonical_fingerprint==nbo_canonical_fingerprint(wavefunction) &&
        options.routed->interaction_graph.available();
    const auto scope=routed?make_fixed_bonding_scope(wavefunction,*options.routed->interaction_graph.value):
        make_fixed_bonding_scope(wavefunction);
    for(auto& level:data.levels) {
        auto members=level.member_indices;
        if(members.empty())members.push_back(level.metadata.orbital_index);
        for(auto index:level.member_spin_counterparts)
            if(index<wavefunction.orbitals.size() && std::find(members.begin(),members.end(),index)==members.end())members.push_back(index);
        std::map<Spin,std::vector<std::size_t>> blocks;
        for(auto index:members)if(index<wavefunction.orbitals.size())blocks[wavefunction.orbitals[index].spin].push_back(index);
        bool positive=!blocks.empty(),negative=!blocks.empty(),mixed=false;
        for(const auto& [spin,indices]:blocks) {
            OrbitalGroupBondingOptions group_options;group_options.expected_dimension=indices.size();
            auto result=analyse_orbital_group_bonding(wavefunction,scope,indices,group_options);
            positive&=result.status==OrbitalGroupBondingStatus::Positive;
            negative&=result.status==OrbitalGroupBondingStatus::Negative;
            mixed|=result.status==OrbitalGroupBondingStatus::Mixed;
            level.bonding_scopes.push_back(std::move(result));
        }
        // Opposite spin blocks can differ; never average away a sign conflict.
        bool any_positive=false,any_negative=false;
        for(const auto& result:level.bonding_scopes) {
            any_positive|=result.status==OrbitalGroupBondingStatus::Positive;
            any_negative|=result.status==OrbitalGroupBondingStatus::Negative;
        }
        mixed|=any_positive&&any_negative;
        const auto role=mixed?BondingClass::Mixed:positive?BondingClass::Bonding:
            negative?BondingClass::Antibonding:BondingClass::Unclassified;
        level.annotation.bonding_class=role;
        level.annotation.bonding_source=role==BondingClass::Unclassified?
            AnnotationSource::Unavailable:AnnotationSource::Derived;
        level.annotation.bonding_confidence=0; // numerical bound lives with the actual group operator
        for(auto index:members)if(index<data.annotations.size()) {
            data.annotations[index].bonding_class=role;
            data.annotations[index].bonding_source=level.annotation.bonding_source;
            data.annotations[index].bonding_confidence=0;
        }
    }
    return data;
}

MetalLigandDetailAvailability metal_ligand_detail_availability(
    const Wavefunction& wavefunction, const MODiagramData& data,
    const MODiagramLevel& level) {
    MetalLigandDetailAvailability result;
    result.composition=data.composition_scope.applicable &&
        level.composition.available && level.composition.complete;
    if(!data.composition_scope.detail.empty() && !data.composition_scope.applicable) {
        result.scope=ChemistryStatus::NotApplicable;
        return result;
    }
    if (wavefunction.atoms.empty()) return result;
    bool metal=false, ligand=false;
    for (const auto& atom:wavefunction.atoms) {
        if (atom.atomic_number<=0) return result;
        if (transition_metal(atom.atomic_number)) metal=true;
        else ligand=true;
    }
    if (!metal || !ligand) {
        result.scope=ChemistryStatus::NotApplicable;
        return result;
    }
    if (data.ligand_field_point_group.empty() ||
        data.ligand_field_metal_atom>=wavefunction.atoms.size() ||
        !transition_metal(wavefunction.atoms[data.ligand_field_metal_atom].atomic_number) ||
        data.ligand_field_ligand_atoms.empty() || level.member_indices.empty()) return result;
    LigandScope scope;
    scope.available=true;
    scope.metal=data.ligand_field_metal_atom;
    for (const auto atom:data.ligand_field_ligand_atoms) {
        if (atom>=wavefunction.atoms.size() || atom==scope.metal ||
            transition_metal(wavefunction.atoms[atom].atomic_number)) return result;
        scope.direct_donors.insert(atom);
    }
    result.scope=ChemistryStatus::Determined;
    result.populations=true;
    result.overlap=true;
    bool valid_channels=true;
    double channel_weight=0.0;
    for (const auto index:level.member_indices) {
        if (index>=wavefunction.orbitals.size() ||
            !wavefunction.orbitals[index].chemistry.available) {
            result.populations=false;
            result.overlap=false;
            valid_channels=false;
            continue;
        }
        const auto& chemistry=wavefunction.orbitals[index].chemistry;
        if (chemistry.ao_contributions.empty()) result.populations=false;
        for (const auto& contribution:chemistry.ao_contributions) {
            if (contribution.atom_index>=wavefunction.atoms.size() ||
                !std::isfinite(contribution.weight)) result.populations=false;
        }
        bool found_pair=false;
        for (const auto& pair:chemistry.interactions) {
            if (!scoped_metal_ligand_pair(wavefunction,pair,&scope)) continue;
            found_pair=true;
            if (!std::isfinite(pair.overlap_character)) {
                result.overlap=false;
                valid_channels=false;
                continue;
            }
            const auto& channel=pair.channel;
            if (!std::isfinite(channel.sigma) || !std::isfinite(channel.pi) ||
                !std::isfinite(channel.delta) || !std::isfinite(channel.phi)) {
                valid_channels=false;
                continue;
            }
            const double resolved=channel.sigma+channel.pi+channel.delta+channel.phi;
            const double weight=std::abs(bounded_overlap_character(pair.overlap_character))*resolved;
            if (channel.status==ChemistryStatus::Determined ||
                channel.status==ChemistryStatus::Percentages) channel_weight+=weight;
            else if (weight>1.0e-12) valid_channels=false;
        }
        if (!found_pair) {
            result.overlap=false;
            valid_channels=false;
        }
    }
    result.populations=result.populations &&
        std::isfinite(level.metal_s_weight) && std::isfinite(level.metal_p_weight) &&
        std::isfinite(level.metal_d_weight) && std::isfinite(level.ligand_p_weight);
    result.overlap=result.overlap && std::isfinite(level.metal_ligand_overlap);
    result.channels=valid_channels && channel_weight>1.0e-12 &&
        std::isfinite(level.sigma_fraction) && std::isfinite(level.pi_fraction);
    return result;
}

} // namespace cov
