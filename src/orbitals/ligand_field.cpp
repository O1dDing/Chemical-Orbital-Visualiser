#include "cov/ligand_field.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/pi_pair_evidence.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace cov {
namespace {

bool transition_metal(const int z) noexcept {
    return (z>=21 && z<=30) || (z>=39 && z<=48) ||
           (z>=72 && z<=80) || (z>=104 && z<=112);
}

struct Neighbour {
    std::size_t atom = 0;
    double mayer = 0.0;
    double distance = 0.0;
    std::array<double,3> direction{};
};

double dot(const std::array<double,3>& a,
           const std::array<double,3>& b) noexcept {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

Neighbour make_neighbour(const Wavefunction& wavefunction,
                         const std::size_t metal,
                         const std::size_t atom,
                         const double mayer) {
    const auto& m=wavefunction.atoms[metal];
    const auto& a=wavefunction.atoms[atom];
    const std::array<double,3> displacement{a.x-m.x,a.y-m.y,a.z-m.z};
    const double distance=std::sqrt(dot(displacement,displacement));
    Neighbour result;
    result.atom=atom;
    result.mayer=std::abs(mayer);
    result.distance=distance;
    if (distance>1.0e-12) {
        result.direction={displacement[0]/distance,
                          displacement[1]/distance,
                          displacement[2]/distance};
    }
    return result;
}

double tetrahedral_confidence(const std::vector<Neighbour>& shell);
double octahedral_confidence(const std::vector<Neighbour>& shell);

std::vector<Neighbour> first_shell(const Wavefunction& wavefunction,
                                   const std::size_t metal) {
    std::vector<Neighbour> candidates;
    for (const auto& bond:wavefunction.bond_orders) {
        std::size_t other=wavefunction.atoms.size();
        if (bond.atom_a==metal) other=bond.atom_b;
        else if (bond.atom_b==metal) other=bond.atom_a;
        if (other>=wavefunction.atoms.size() ||
            transition_metal(wavefunction.atoms[other].atomic_number)) continue;
        const double magnitude=std::abs(bond.mayer_order);
        if (magnitude<0.02) continue;
        candidates.push_back(make_neighbour(
            wavefunction,metal,other,bond.mayer_order));
    }
    if (candidates.empty()) return {};
    candidates.erase(std::remove_if(
        candidates.begin(),candidates.end(),[](const auto& item) {
            return item.distance<=1.0e-12;
        }),candidates.end());
    std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
        if (a.distance!=b.distance) return a.distance<b.distance;
        return a.mayer>b.mayer;
    });

    // A carbonyl oxygen, for example, can carry a small density-level Mayer
    // coupling to the metal but lies behind the directly coordinated carbon.
    // Retain only the nearest atom along each ligand ray.
    std::vector<Neighbour> radial;
    for (const auto& candidate:candidates) {
        const bool shadowed=std::any_of(
            radial.begin(),radial.end(),[&](const auto& retained) {
                return dot(candidate.direction,retained.direction)>0.94;
            });
        if (!shadowed) radial.push_back(candidate);
    }
    if (radial.size()>6u) {
        std::stable_sort(radial.begin(),radial.end(),[](const auto& a,const auto& b) {
            const double sa=a.mayer/std::max(1.0e-9,a.distance);
            const double sb=b.mayer/std::max(1.0e-9,b.distance);
            return sa>sb;
        });
        radial.resize(6u);
    }

    // Set the electronic cutoff only after collinear through-ligand atoms
    // have been removed.  Otherwise a remote M...O/N density coupling can
    // raise the floor and delete the nearer direct donor before radial
    // shadowing has a chance to act.
    const double strongest=std::accumulate(
        radial.begin(),radial.end(),0.0,[](const double value,const auto& item) {
            return std::max(value,item.mayer);
        });
    const double floor=std::max(0.035,0.08*strongest);
    auto electronically_strong=radial;
    electronically_strong.erase(std::remove_if(
        electronically_strong.begin(),electronically_strong.end(),
        [&](const auto& item) {return item.mayer<floor;}),
        electronically_strong.end());
    const auto supported=[](const auto& shell) {
        return (shell.size()==4u && tetrahedral_confidence(shell)>=0.55) ||
               (shell.size()==6u && octahedral_confidence(shell)>=0.55);
    };
    if (supported(electronically_strong)) return electronically_strong;
    // Ionic first-shell bonds can all lie below 0.035.  The low-Mayer retry
    // is accepted only when the complete radial-nearest set independently
    // forms a supported Td/Oh geometry; this does not lower the global floor.
    if (supported(radial)) return radial;
    return electronically_strong;
}

double tetrahedral_confidence(const std::vector<Neighbour>& shell) {
    if (shell.size()!=4u) return 0.0;
    double squared=0.0;
    std::size_t count=0u;
    for (std::size_t i=0;i<shell.size();++i) {
        for (std::size_t j=i+1u;j<shell.size();++j) {
            const double error=dot(shell[i].direction,shell[j].direction)+1.0/3.0;
            squared+=error*error;
            ++count;
        }
    }
    const double rms=std::sqrt(squared/static_cast<double>(count));
    return std::clamp(1.0-rms/0.30,0.0,1.0);
}

double octahedral_confidence(const std::vector<Neighbour>& shell) {
    if (shell.size()!=6u) return 0.0;
    struct Pair { std::size_t a=0; std::size_t b=0; double value=0.0; };
    std::vector<Pair> pairs;
    for (std::size_t i=0;i<shell.size();++i) {
        for (std::size_t j=i+1u;j<shell.size();++j) {
            pairs.push_back({i,j,dot(shell[i].direction,shell[j].direction)});
        }
    }
    std::sort(pairs.begin(),pairs.end(),[](const auto& a,const auto& b) {
        return a.value<b.value;
    });
    std::set<std::size_t> used;
    std::vector<Pair> opposite;
    for (const auto& pair:pairs) {
        if (used.count(pair.a) || used.count(pair.b)) continue;
        used.insert(pair.a);
        used.insert(pair.b);
        opposite.push_back(pair);
        if (opposite.size()==3u) break;
    }
    if (opposite.size()!=3u ||
        std::any_of(opposite.begin(),opposite.end(),[](const auto& pair) {
            return pair.value>-0.72;
        })) return 0.0;

    double opposite_squared=0.0;
    for (const auto& pair:opposite) {
        const double error=pair.value+1.0;
        opposite_squared+=error*error;
    }
    double right_squared=0.0;
    std::size_t right_count=0u;
    for (const auto& pair:pairs) {
        const bool is_opposite=std::any_of(
            opposite.begin(),opposite.end(),[&](const auto& selected) {
                return selected.a==pair.a && selected.b==pair.b;
            });
        if (is_opposite) continue;
        right_squared+=pair.value*pair.value;
        ++right_count;
    }
    const double opposite_rms=std::sqrt(opposite_squared/3.0);
    const double right_rms=std::sqrt(
        right_squared/static_cast<double>(std::max<std::size_t>(1u,right_count)));
    return std::clamp(1.0-0.45*opposite_rms/0.25-
                          0.55*right_rms/0.35,0.0,1.0);
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

double central_metal_weight(const MolecularOrbital& orbital,
                            const std::size_t metal_atom,
                            const int angular_momentum) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index==metal_atom &&
            contribution.angular_momentum==angular_momentum) {
            result+=contribution.weight;
        }
    }
    return result;
}

struct LocalOrbitalGroup {
    OrbitalSymmetryExplanation explanation;
    std::vector<std::size_t> members;
    Spin spin=Spin::Alpha;
    double energy=0.0;
    std::array<double,3> metal{};
    double ligand_p=0.0;
    double pi_fraction=0.0;
    double metal_ligand_overlap=0.0;
    std::string symmetry;
    std::uint8_t local_irrep_copy=0;
    bool locally_classified=false;
};

enum class LocalPiRole {
    Unresolved,
    SigmaOnly,
    Donor,
    Acceptor,
    Ambiguous,
};

LocalPiRole local_pi_role(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment) {
    std::size_t sigma_only=0u;
    std::size_t donor=0u;
    std::size_t acceptor=0u;
    std::size_t ambiguous=0u;
    for (const auto atom:environment.ligand_atoms) {
        if (atom>=wavefunction.atoms.size()) continue;
        const int z=wavefunction.atoms[atom].atomic_number;
        std::size_t external_neighbour_count=0u;
        bool multiple_heavy_bond=false;
        for (const auto& bond:wavefunction.bond_orders) {
            std::size_t other=wavefunction.atoms.size();
            if (bond.atom_a==atom) other=bond.atom_b;
            else if (bond.atom_b==atom) other=bond.atom_a;
            if (other>=wavefunction.atoms.size() ||
                other==environment.metal_atom ||
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
            ++acceptor;
        } else if ((z==7 || z==15) && !multiple_heavy_bond &&
                   external_neighbour_count>=3u) {
            ++sigma_only;
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
    if (classified==0u) return LocalPiRole::Unresolved;
    if (sigma_only==classified) return LocalPiRole::SigmaOnly;
    if (acceptor>0u && donor==0u && ambiguous==0u) {
        return LocalPiRole::Acceptor;
    }
    if (donor>0u && acceptor==0u && ambiguous==0u) {
        return LocalPiRole::Donor;
    }
    return LocalPiRole::Ambiguous;
}

bool direct_metal_ligand_pair(
    const OrbitalPairInteraction& interaction,
    const LigandFieldEnvironment& environment) {
    const auto direct=[&](const std::size_t atom) {
        return std::find(environment.ligand_atoms.begin(),
                         environment.ligand_atoms.end(),atom)!=
               environment.ligand_atoms.end();
    };
    return (interaction.atom_a==environment.metal_atom &&
            direct(interaction.atom_b)) ||
           (interaction.atom_b==environment.metal_atom &&
            direct(interaction.atom_a));
}

class SpinLabelOverlapWorkspace {
public:
    explicit SpinLabelOverlapWorkspace(const Wavefunction& wavefunction)
        : basis_count_(wavefunction.basis_count),
          beta_position_(wavefunction.orbitals.size(),invalid_index()) {
        if (basis_count_==0u ||
            wavefunction.ao_overlap.size()!=basis_count_*basis_count_) return;
        for (std::size_t orbital=0;orbital<wavefunction.orbitals.size();++orbital) {
            const auto& item=wavefunction.orbitals[orbital];
            if (item.spin!=Spin::Beta || item.coefficients.size()!=basis_count_) {
                continue;
            }
            beta_position_[orbital]=beta_orbitals_.size();
            beta_orbitals_.push_back(orbital);
        }
        transformed_beta_.assign(beta_orbitals_.size()*basis_count_,0.0);
        for (std::size_t beta=0;beta<beta_orbitals_.size();++beta) {
            const auto& coefficients=
                wavefunction.orbitals[beta_orbitals_[beta]].coefficients;
            double* transformed=transformed_beta_.data()+beta*basis_count_;
            for (std::size_t mu=0;mu<basis_count_;++mu) {
                double value=0.0;
                for (std::size_t nu=0;nu<basis_count_;++nu) {
                    value+=wavefunction.ao_overlap[mu*basis_count_+nu]*
                           static_cast<double>(coefficients[nu]);
                }
                transformed[mu]=value;
            }
        }
    }

    [[nodiscard]] double overlap(const Wavefunction& wavefunction,
                                 const std::size_t alpha,
                                 const std::size_t beta) const noexcept {
        if (alpha>=wavefunction.orbitals.size() ||
            beta>=beta_position_.size() ||
            wavefunction.orbitals[alpha].coefficients.size()!=basis_count_) {
            return 0.0;
        }
        const std::size_t position=beta_position_[beta];
        if (position==invalid_index()) return 0.0;
        const double* transformed=
            transformed_beta_.data()+position*basis_count_;
        double result=0.0;
        for (std::size_t mu=0;mu<basis_count_;++mu) {
            result+=static_cast<double>(
                wavefunction.orbitals[alpha].coefficients[mu])*transformed[mu];
        }
        return result;
    }

private:
    [[nodiscard]] static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }

    std::size_t basis_count_=0u;
    std::vector<std::size_t> beta_orbitals_;
    std::vector<std::size_t> beta_position_;
    std::vector<double> transformed_beta_;
};

void propagate_spin_counterpart_labels(
    const Wavefunction& wavefunction,
    std::vector<OrbitalMetadata>& metadata) {
    const bool unrestricted=std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const auto& orbital) {return orbital.spin==Spin::Beta;});
    if (!unrestricted || wavefunction.ao_overlap.empty()) return;

    const SpinLabelOverlapWorkspace workspace(wavefunction);
    struct BestMatch {
        std::size_t index=std::numeric_limits<std::size_t>::max();
        double score=0.0;
        double runner_up=0.0;
    };
    std::vector<BestMatch> alpha_best(wavefunction.orbitals.size());
    std::vector<BestMatch> beta_best(wavefunction.orbitals.size());
    const auto consider=[](BestMatch& best,const std::size_t index,
                           const double score) {
        if (score>best.score) {
            best.runner_up=best.score;
            best.score=score;
            best.index=index;
        } else if (score>best.runner_up) {
            best.runner_up=score;
        }
    };
    for (std::size_t alpha=0;alpha<wavefunction.orbitals.size();++alpha) {
        if (wavefunction.orbitals[alpha].spin!=Spin::Alpha) continue;
        for (std::size_t beta=0;beta<wavefunction.orbitals.size();++beta) {
            if (wavefunction.orbitals[beta].spin!=Spin::Beta) continue;
            const double overlap=workspace.overlap(wavefunction,alpha,beta);
            const double score=overlap*overlap;
            consider(alpha_best[alpha],beta,score);
            consider(beta_best[beta],alpha,score);
        }
    }

    constexpr double minimum_score=0.90;
    constexpr double minimum_margin=0.10;
    for (std::size_t alpha=0;alpha<wavefunction.orbitals.size();++alpha) {
        if (wavefunction.orbitals[alpha].spin!=Spin::Alpha) continue;
        const auto& match=alpha_best[alpha];
        if (match.index>=wavefunction.orbitals.size() ||
            match.score<minimum_score ||
            match.score-match.runner_up<minimum_margin) continue;
        const std::size_t beta=match.index;
        const auto& reverse=beta_best[beta];
        if (reverse.index!=alpha || reverse.score<minimum_score ||
            reverse.score-reverse.runner_up<minimum_margin) continue;
        const bool alpha_labelled=!normalised_symmetry(
            metadata[alpha].symmetry_view.label).empty();
        const bool beta_labelled=!normalised_symmetry(
            metadata[beta].symmetry_view.label).empty();
        if (alpha_labelled==beta_labelled) continue;
        if (alpha_labelled) {
            const std::array<std::size_t,1> target{beta};
            metadata[beta].symmetry_view=candidate_orbital_symmetry(
                metadata[alpha].symmetry_view,OrbitalSymmetryOrigin::SpinCounterpartCandidate,target,match.score);
        } else {
            const std::array<std::size_t,1> target{alpha};
            metadata[alpha].symmetry_view=candidate_orbital_symmetry(
                metadata[beta].symmetry_view,OrbitalSymmetryOrigin::SpinCounterpartCandidate,target,match.score);
        }
    }
}

int dominant_family(const LocalOrbitalGroup& group) {
    const auto found=std::max_element(group.metal.begin(),group.metal.end());
    if (found==group.metal.end()) return -1;
    const int family=static_cast<int>(std::distance(group.metal.begin(),found));
    const double floor=family==2?kLocalDIrrepWeightFloor:0.08;
    return *found>=floor?family:-1;
}

LocalOrbitalGroup merge_groups(const LocalOrbitalGroup& left,
                               const LocalOrbitalGroup& right) {
    LocalOrbitalGroup result=left;
    const double left_count=static_cast<double>(left.members.size());
    const double right_count=static_cast<double>(right.members.size());
    const double total=std::max(1.0,left_count+right_count);
    result.energy=(left.energy*left_count+right.energy*right_count)/total;
    for (std::size_t family=0;family<result.metal.size();++family) {
        result.metal[family]=(
            left.metal[family]*left_count+right.metal[family]*right_count)/total;
    }
    result.ligand_p=(left.ligand_p*left_count+
                     right.ligand_p*right_count)/total;
    result.pi_fraction=(left.pi_fraction*left_count+
                        right.pi_fraction*right_count)/total;
    result.metal_ligand_overlap=(
        left.metal_ligand_overlap*left_count+
        right.metal_ligand_overlap*right_count)/total;
    result.members.insert(
        result.members.end(),right.members.begin(),right.members.end());
    if (result.symmetry.empty()) result.symmetry=right.symmetry;
    if (left.local_irrep_copy!=right.local_irrep_copy) {
        result.local_irrep_copy=0;
    }
    result.locally_classified=left.locally_classified &&
                              right.locally_classified;
    return result;
}

bool merge_resolved_five_d_run(std::vector<LocalOrbitalGroup>& groups) {
    constexpr double d_tolerance_hartree=0.005;
    const auto try_run=[&](const std::vector<std::size_t>& run) {
        if (run.size()!=5u) return false;
        std::array<double,4> gaps{};
        for (std::size_t i=0;i<gaps.size();++i) {
            gaps[i]=groups[run[i+1u]].energy-groups[run[i]].energy;
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

        LocalOrbitalGroup lower=groups[run.front()];
        for (std::size_t i=1u;i<split;++i) {
            lower=merge_groups(lower,groups[run[i]]);
        }
        LocalOrbitalGroup upper=groups[run[split]];
        for (std::size_t i=split+1u;i<run.size();++i) {
            upper=merge_groups(upper,groups[run[i]]);
        }
        const std::set<std::size_t> consumed(run.begin(),run.end());
        std::vector<LocalOrbitalGroup> rebuilt;
        rebuilt.reserve(groups.size()-3u);
        for (std::size_t i=0;i<groups.size();++i) {
            if (i==run.front()) rebuilt.push_back(lower);
            else if (i==run[split]) rebuilt.push_back(upper);
            else if (consumed.count(i)==0u) rebuilt.push_back(groups[i]);
        }
        groups=std::move(rebuilt);
        std::stable_sort(groups.begin(),groups.end(),[](const auto& a,const auto& b) {
            if (a.spin!=b.spin) return a.spin==Spin::Alpha;
            return a.energy<b.energy;
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
            if (groups[i].spin!=spin) continue;
            const int family=dominant_family(groups[i]);
            const auto& scope=groups[i].explanation;
            const bool candidate=family==2 && groups[i].members.size()==1u &&
                !scope.local_assignment && orbital_symmetry_missing_local_input(scope);
            if (!candidate) {
                if (family==2 && flush()) return true;
                continue;
            }
            if (!run.empty() && groups[i].energy-groups[run.back()].energy>
                                    d_tolerance_hartree) {
                if (flush()) return true;
            }
            run.push_back(i);
        }
        if (flush()) return true;
    }
    return false;
}

std::vector<LocalOrbitalGroup> initial_orbital_groups(
    const Wavefunction& wavefunction,
    const std::vector<OrbitalMetadata>& metadata,
    const LigandFieldEnvironment& environment) {
    std::vector<LocalOrbitalGroup> result;
    const std::size_t metal_atom=environment.metal_atom;
    const double donor_count=static_cast<double>(
        std::max<std::size_t>(1u,environment.ligand_atoms.size()));
    for (std::size_t base=0;base<wavefunction.orbitals.size();) {
        const Spin spin=wavefunction.orbitals[base].spin;
        const std::size_t requested=base<metadata.size()
            ?std::max<std::size_t>(1u,metadata[base].degeneracy_size):1u;
        std::size_t end=std::min(wavefunction.orbitals.size(),base+requested);
        for (std::size_t index=base+1u;index<end;++index) {
            if (wavefunction.orbitals[index].spin!=spin) {
                end=index;
                break;
            }
        }
        LocalOrbitalGroup group;
        double pi_weighted=0.0;
        double channel_weight=0.0;
        double overlap=0.0;
        group.spin=spin;
        // Local classification starts without borrowing a molecular label.
        // Identical text in a producer group need not mean the same local axes.
        for (std::size_t index=base;index<end;++index) {
            const auto& orbital=wavefunction.orbitals[index];
            group.members.push_back(index);
            group.energy+=orbital.energy_hartree;
            for (int family=0;family<3;++family) {
                group.metal[static_cast<std::size_t>(family)]+=
                    central_metal_weight(orbital,metal_atom,family);
            }
            for (const auto& contribution:orbital.chemistry.ao_contributions) {
                if (contribution.atom_index!=metal_atom &&
                    contribution.angular_momentum==1) {
                    group.ligand_p+=contribution.weight;
                }
            }
            for (const auto& interaction:orbital.chemistry.interactions) {
                if (!direct_metal_ligand_pair(interaction,environment)) {
                    continue;
                }
                const double bounded=std::tanh(interaction.overlap_character);
                const double magnitude=std::abs(bounded);
                const double determined=interaction.channel.sigma+
                    interaction.channel.pi+interaction.channel.delta+
                    interaction.channel.phi;
                pi_weighted+=magnitude*interaction.channel.pi;
                channel_weight+=magnitude*determined;
                overlap+=bounded;
            }
        }
        const double count=static_cast<double>(
            std::max<std::size_t>(1u,group.members.size()));
        group.energy/=count;
        for (auto& weight:group.metal) weight/=count;
        group.ligand_p/=count;
        if (channel_weight>1.0e-12) {
            group.pi_fraction=std::clamp(
                pi_weighted/channel_weight,0.0,1.0);
        }
        group.metal_ligand_overlap=std::clamp(
            overlap/(count*donor_count),-1.0,1.0);
        result.push_back(std::move(group));
        base=std::max(base+1u,end);
    }
    return result;
}

bool local_d_irrep(const std::string& point_group,
                   const std::string& symmetry) {
    const std::string wanted=normalised_symmetry(symmetry);
    if (wanted.empty()) return false;
    const auto decomposition=decompose_metal_ao_shell(
        point_group,MetalAOShell::D);
    return decomposition && std::any_of(
        decomposition->begin(),decomposition->end(),[&](const auto& copy) {
            return normalised_symmetry(std::string(copy.label))==wanted;
        });
}

void propagate_pi_partner_labels(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment,
    std::vector<LocalOrbitalGroup>& groups) {
    const std::string point_group=environment.local_point_group();
    const LocalPiRole role=local_pi_role(wavefunction,environment);
    const LigandPiPrior prior=role==LocalPiRole::SigmaOnly?LigandPiPrior::SigmaOnly:
        role==LocalPiRole::Donor?LigandPiPrior::Donor:
        role==LocalPiRole::Acceptor?LigandPiPrior::Acceptor:
        role==LocalPiRole::Ambiguous?LigandPiPrior::Ambiguous:LigandPiPrior::Unresolved;
    struct Candidate {
        std::size_t known=0u;
        std::size_t missing=0u;
        double score=0.0;
        std::shared_ptr<const PiPartnerAssessment> evidence;
    };
    std::vector<Candidate> candidates;
    for (std::size_t first=0;first<groups.size();++first) {
        for (std::size_t second=first+1u;second<groups.size();++second) {
            auto& a=groups[first];
            auto& b=groups[second];
            if (a.spin!=b.spin || a.members.size()!=b.members.size()) continue;
            const bool a_known=local_d_irrep(point_group,a.symmetry);
            const bool b_known=local_d_irrep(point_group,b.symmetry);
            if (a_known==b_known) continue;
            const auto& missing=a_known?b:a;
            if (!normalised_symmetry(missing.symmetry).empty()) continue;
            const auto& lower=a.energy<=b.energy?a:b;
            const auto& upper=a.energy<=b.energy?b:a;
            auto evidence=std::make_shared<const PiPartnerAssessment>(assess_pi_partner(
                {lower.energy,lower.metal[2],lower.ligand_p,lower.pi_fraction,lower.metal_ligand_overlap},
                {upper.energy,upper.metal[2],upper.ligand_p,upper.pi_fraction,upper.metal_ligand_overlap},prior));
            if(!evidence->accepted)continue;
            // Weak scalar similarity is not a whole-subspace irrep proof.
            // Only the independent metric/symmetry analyses may fill a label.
            if(evidence->channel.channel_id.empty())continue;
            candidates.push_back({a_known?first:second,a_known?second:first,evidence->ranking_score,evidence});
        }
    }
    std::sort(candidates.begin(),candidates.end(),[](const auto& left,
                                                     const auto& right) {
        return left.score>right.score;
    });
    std::set<std::size_t> used;
    for (const auto& candidate:candidates) {
        if (used.count(candidate.known) || used.count(candidate.missing)) {
            continue;
        }
        const auto& source=groups[candidate.known];
        auto& target=groups[candidate.missing];
        target.symmetry=source.symmetry;
        target.explanation=candidate_orbital_symmetry(source.explanation,
            OrbitalSymmetryOrigin::PiPartnerCandidate,target.members,candidate.score);
        target.explanation.pi_partner_evidence=candidate.evidence;
        used.insert(candidate.known);
        used.insert(candidate.missing);
    }
}

void merge_local_pseudodegenerate_groups(
    const std::string& point_group,
    std::vector<LocalOrbitalGroup>& groups) {
    const bool unique_two_three_d=
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,2u).has_value() &&
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,3u).has_value();
    if (unique_two_three_d) {
        while (merge_resolved_five_d_run(groups)) {
            // Restart because non-adjacent ligand-centred rows may have been
            // crossed while forming the two resolved d subspaces.
        }
    }
    for (std::size_t index=0;index+1u<groups.size();) {
        const auto& left=groups[index];
        const auto& right=groups[index+1u];
        const std::string left_symmetry=normalised_symmetry(left.symmetry);
        const std::string right_symmetry=normalised_symmetry(right.symmetry);
        const bool local_copies_compatible=
            (!left.locally_classified && !right.locally_classified) ||
            (left.locally_classified && right.locally_classified &&
             left.local_irrep_copy!=0u &&
             left.local_irrep_copy==right.local_irrep_copy);
        const bool labels_compatible=
            (left_symmetry.empty() && right_symmetry.empty()) ||
            (!left_symmetry.empty() && left_symmetry==right_symmetry &&
             local_copies_compatible);
        const int family=dominant_family(left);
        if (left.spin!=right.spin || !labels_compatible || family<1 ||
            family!=dominant_family(right)) {
            ++index;
            continue;
        }
        const std::size_t left_size=left.members.size();
        const std::size_t right_size=right.members.size();
        const std::size_t combined=left_size+right_size;
        const bool unlabeled=left_symmetry.empty() && right_symmetry.empty();
        const bool dimension_unambiguous=!unlabeled ||
            classify_local_irrep_by_dimension(
                point_group,static_cast<MetalAOShell>(family),combined)
                .has_value();
        // A C1/Cs producer commonly emits a local metal-p T2/T1u block as
        // three single rows.  Permit the intermediate 1+1 -> 2 merge so the
        // following pass can complete 2+1 -> 3; requiring combined==3 here
        // made that recovery impossible.
        const bool expected_p=family==1 && combined<=3u &&
            (left_size<3u || right_size<3u);
        const bool expected_d=family==2 &&
            ((combined==3u && (left_size==1u || right_size==1u)) ||
             (combined==2u && left_size==1u && right_size==1u));
        const double tolerance=family==1?0.015:0.005;
        if (!dimension_unambiguous || (!expected_p && !expected_d) ||
            std::abs(right.energy-left.energy)>tolerance) {
            ++index;
            continue;
        }
        groups[index]=merge_groups(left,right);
        groups.erase(groups.begin()+static_cast<std::ptrdiff_t>(index+1u));
    }
}

void recover_local_group_symmetry(
    const LocalAngularProjectionWorkspace& workspace,
    const LigandFieldEnvironment& environment,
    std::vector<LocalOrbitalGroup>& groups) {
    if (!environment.available()) return;
    for(auto& group:groups) {
        const int family=dominant_family(group);
        if(family<0) {
            group.symmetry.clear();group.explanation={};
            group.local_irrep_copy=0;group.locally_classified=false;
            continue;
        }
        group.explanation=evaluate_local_orbital_symmetry(workspace,environment,group.members,family);
        group.symmetry=group.explanation.label;
        group.local_irrep_copy=group.explanation.local_assignment
            ?group.explanation.local_assignment->copy_index:0u;
        group.locally_classified=static_cast<bool>(group.explanation.local_assignment);
    }
}

} // namespace

std::string LigandFieldEnvironment::local_point_group() const {
    const auto* descriptor=coordination_geometry_descriptor(geometry_id);
    return descriptor==nullptr?std::string{}:std::string(descriptor->point_group);
}

std::string LigandFieldEnvironment::geometry_machine_id() const {
    const auto* descriptor=coordination_geometry_descriptor(geometry_id);
    return descriptor==nullptr?std::string{}:std::string(descriptor->machine_id);
}

std::string LigandFieldEnvironment::geometry_name() const {
    const auto* descriptor=coordination_geometry_descriptor(geometry_id);
    return descriptor==nullptr?std::string{}:std::string(descriptor->name);
}

LigandFieldEnvironment analyse_ligand_field_environment(
    const Wavefunction& wavefunction) {
    LigandFieldEnvironment best;
    for (std::size_t metal=0;metal<wavefunction.atoms.size();++metal) {
        if (!transition_metal(wavefunction.atoms[metal].atomic_number)) continue;
        const CoordinationShell shell=extract_coordination_shell(
            wavefunction,metal);
        const GeometryMatch match=analyse_coordination_shell(shell);
        if (!match.accepted || match.ambiguous || !match.best) continue;

        LigandFieldEnvironment candidate;
        candidate.metal_atom=metal;
        candidate.geometry_id=match.best->id;
        candidate.angular_rms=match.best->angular_rms;
        candidate.shape_measure=match.best->shape_measure;
        candidate.radial_cv=match.best->radial_cv;
        candidate.rotation_reference_to_input=
            match.best->rotation_reference_to_input;
        candidate.ambiguous=match.ambiguous;
        switch (candidate.geometry_id) {
            case GeometryId::Tetrahedral4:
                candidate.geometry=LigandFieldGeometry::Tetrahedral;
                break;
            case GeometryId::Octahedral6:
                candidate.geometry=LigandFieldGeometry::Octahedral;
                break;
            default:
                candidate.geometry=LigandFieldGeometry::General;
                break;
        }
        for (const auto& item:shell.contacts) {
            candidate.ligand_atoms.push_back(item.atom_index);
        }

        constexpr double accepted_rms=0.40;
        const double fit=std::clamp(
            1.0-candidate.angular_rms/accepted_rms,0.0,1.0);
        double separation=1.0;
        if (match.runner_up) {
            separation=std::clamp(
                (match.runner_up->angular_rms-candidate.angular_rms)/0.12,
                0.0,1.0);
        }
        const double radial_penalty=0.15*std::clamp(
            candidate.radial_cv/0.80,0.0,1.0);
        candidate.confidence=std::clamp(
            fit*(0.75+0.25*separation)*(1.0-radial_penalty),0.0,1.0);
        if (!candidate.ligand_atoms.empty()) {
            const int first_z=wavefunction.atoms[
                candidate.ligand_atoms.front()].atomic_number;
            candidate.equivalent_ligand_elements=std::all_of(
                candidate.ligand_atoms.begin(),candidate.ligand_atoms.end(),
                [&](const auto atom) {
                    return wavefunction.atoms[atom].atomic_number==first_z;
                });
        }
        if (!best.available() || candidate.confidence>best.confidence) {
            best=std::move(candidate);
        }
    }
    return best;
}

void apply_local_ligand_field_symmetry(
    const Wavefunction& wavefunction,
    std::vector<OrbitalMetadata>& metadata) {
    if (metadata.size()!=wavefunction.orbitals.size()) return;
    const LigandFieldEnvironment environment=
        analyse_ligand_field_environment(wavefunction);
    if (!environment.available()) return;

    const LocalAngularProjectionWorkspace workspace(wavefunction,environment.metal_atom,
        environment.rotation_reference_to_input);
    auto groups=initial_orbital_groups(wavefunction,metadata,environment);
    recover_local_group_symmetry(workspace,environment,groups);
    if (environment.equivalent_ligand_elements) {
        merge_local_pseudodegenerate_groups(
            environment.local_point_group(),groups);
        recover_local_group_symmetry(workspace,environment,groups);
    }
    propagate_pi_partner_labels(wavefunction,environment,groups);
    for (const auto& group:groups) {
        const std::string& symmetry=group.symmetry;
        if (symmetry.empty()) continue;
        for (const auto member:group.members) {
            if (member>=metadata.size()) continue;
            metadata[member].symmetry_view=group.explanation;
        }
    }
    propagate_spin_counterpart_labels(wavefunction,metadata);
}

} // namespace cov
