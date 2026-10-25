/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#include "cov/bond_analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

std::size_t packed_index(std::size_t i, std::size_t j) {
    if (j > i) std::swap(i, j);
    return i * (i + 1u) / 2u + j;
}

std::vector<double> unpack_symmetric(const std::vector<double>& packed,
                                     const std::size_t n,
                                     const char* label) {
    const std::size_t expected = n * (n + 1u) / 2u;
    if (packed.size() != expected) {
        throw std::runtime_error(std::string(label) +
                                 " packed dimension does not match AO basis");
    }
    std::vector<double> full(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double value = packed[packed_index(i, j)];
            full[i * n + j] = value;
            full[j * n + i] = value;
        }
    }
    return full;
}

std::vector<std::uint32_t> basis_atom_map(const Wavefunction& wf) {
    const std::size_t n = wf.basis_count;
    const auto missing = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> map(n, missing);
    for (const auto& shell : wf.shells) {
        const std::size_t begin = shell.basis_offset;
        const std::size_t count = shell_basis_count(shell);
        if (begin + count > n || shell.atom_index >= wf.atoms.size()) {
            throw std::runtime_error("Shell/basis mapping is inconsistent during bond analysis");
        }
        for (std::size_t local = 0; local < count; ++local) {
            const std::size_t basis = begin + local;
            if (map[basis] != missing) {
                throw std::runtime_error("AO basis function is assigned to more than one shell");
            }
            map[basis] = shell.atom_index;
        }
    }
    for (const auto atom : map) {
        if (atom == missing) {
            throw std::runtime_error("AO basis function is not assigned to an atom");
        }
    }
    return map;
}

std::vector<double> multiply(const std::vector<double>& a,
                             const std::vector<double>& b,
                             const std::size_t n) {
    std::vector<double> out(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
            const double aik = a[i * n + k];
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < n; ++j) {
                out[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    return out;
}

std::vector<double> atomic_participation(const Wavefunction& wf,
                                         const std::vector<std::uint32_t>& basis_atom,
                                         const MolecularOrbital& mo) {
    const std::size_t n = wf.basis_count;
    if (mo.coefficients.size() != n || wf.ao_overlap.size() != n * n) return {};

    std::vector<double> sc(n, 0.0);
    for (std::size_t mu = 0; mu < n; ++mu) {
        double value = 0.0;
        for (std::size_t nu = 0; nu < n; ++nu) {
            value += wf.ao_overlap[mu * n + nu] *
                     static_cast<double>(mo.coefficients[nu]);
        }
        sc[mu] = value;
    }

    std::vector<double> atom_weight(wf.atoms.size(), 0.0);
    for (std::size_t mu = 0; mu < n; ++mu) {
        atom_weight[basis_atom[mu]] +=
            static_cast<double>(mo.coefficients[mu]) * sc[mu];
    }

    double abs_sum = 0.0;
    for (double& value : atom_weight) {
        value = std::abs(value);
        abs_sum += value;
    }
    if (!(abs_sum > 1.0e-12) || !std::isfinite(abs_sum)) return {};
    for (double& value : atom_weight) value /= abs_sum;
    return atom_weight;
}

using AtomTriple = std::array<std::uint32_t, 3>;

AtomTriple sorted_triple(std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    AtomTriple out{a,b,c};
    std::sort(out.begin(),out.end());
    return out;
}

int pair_support_count(const Wavefunction& wf,
                       const AtomTriple& atoms,
                       const double minimum) {
    int supported = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = i + 1; j < 3; ++j) {
            for (const auto& bond : wf.bond_orders) {
                const bool same =
                    (bond.atom_a == atoms[i] && bond.atom_b == atoms[j]) ||
                    (bond.atom_a == atoms[j] && bond.atom_b == atoms[i]);
                if (same && std::abs(bond.mayer_order) >= minimum) {
                    ++supported;
                    break;
                }
            }
        }
    }
    return supported;
}

bool near_occupation(const double value,
                     const double target,
                     const double tolerance) {
    return std::abs(value-target) <= tolerance;
}

struct ActiveMO {
    std::uint32_t index = 0;
    double occupation = 0.0;
    double energy = 0.0;
    double coverage = 0.0;
    double score = 0.0;
    bool full_three = false;
    bool pair_like = false;
};

std::vector<ActiveMO> active_mos_for_triple(
    const Wavefunction& wf,
    const std::vector<std::vector<double>>& weights,
    const AtomTriple& atoms,
    const BondAnalysisOptions& options) {
    std::vector<ActiveMO> active;
    for (std::size_t i = 0; i < wf.orbitals.size(); ++i) {
        const auto& mo = wf.orbitals[i];
        if (mo.spin != Spin::Alpha || i >= weights.size() || weights[i].empty()) continue;
        if (mo.energy_hartree < options.assignment_min_energy_hartree ||
            mo.energy_hartree > options.assignment_max_energy_hartree) continue;

        std::array<double,3> local{
            weights[i][atoms[0]], weights[i][atoms[1]], weights[i][atoms[2]]
        };
        const double coverage=local[0]+local[1]+local[2];
        if (coverage < options.assignment_triple_coverage) continue;
        std::sort(local.begin(),local.end(),std::greater<double>());
        const bool full_three = local[2] >= options.assignment_three_atom_floor;
        const bool pair_like =
            local[0] >= options.assignment_pair_atom_floor &&
            local[1] >= options.assignment_pair_atom_floor &&
            local[2] <= options.assignment_pair_third_ceiling;
        if (!full_three && !pair_like) continue;

        ActiveMO entry;
        entry.index=static_cast<std::uint32_t>(i);
        entry.occupation=static_cast<double>(mo.occupation);
        entry.energy=mo.energy_hartree;
        entry.coverage=coverage;
        entry.full_three=full_three;
        entry.pair_like=pair_like;
        entry.score=coverage + 0.35*(1.0-local[0]) + (full_three?0.15:0.05);
        active.push_back(entry);
    }

    std::sort(active.begin(),active.end(),[](const ActiveMO& a,const ActiveMO& b) {
        if (a.score != b.score) return a.score > b.score;
        return a.index < b.index;
    });
    if (active.size() > 12u) active.resize(12u);
    return active;
}

struct AssignmentChoice {
    MulticentreKind kind = MulticentreKind::Unclassified;
    std::array<std::uint32_t,3> orbitals{};
    double electron_count = 0.0;
    double score = -std::numeric_limits<double>::infinity();
    double mean_coverage = 0.0;
};

AssignmentChoice best_assignment_for_triple(
    const Wavefunction& wf,
    const AtomTriple& atoms,
    const std::vector<ActiveMO>& active,
    const std::set<std::uint32_t>& seed_orbitals,
    const BondAnalysisOptions& options,
    double& second_best_score) {
    AssignmentChoice best;
    second_best_score=-std::numeric_limits<double>::infinity();
    if (active.size()<3u || pair_support_count(wf,atoms,options.assignment_min_mayer_support)<2) {
        return best;
    }

    for (std::size_t ia=0; ia+2<active.size(); ++ia) {
        for (std::size_t ib=ia+1; ib+1<active.size(); ++ib) {
            for (std::size_t ic=ib+1; ic<active.size(); ++ic) {
                const std::array<const ActiveMO*,3> chosen{
                    &active[ia],&active[ib],&active[ic]
                };

                int full_count=0;
                int occupied_count=0;
                bool occupation_clean=true;
                bool contains_seed=false;
                double electrons=0.0;
                double score=0.0;
                double coverage=0.0;
                for (const ActiveMO* entry:chosen) {
                    full_count += entry->full_three?1:0;
                    contains_seed = contains_seed || seed_orbitals.count(entry->index)!=0u;
                    score += entry->score;
                    coverage += entry->coverage;
                    electrons += entry->occupation;
                    if (near_occupation(entry->occupation,2.0,options.assignment_occupation_tolerance)) {
                        ++occupied_count;
                    } else if (!near_occupation(entry->occupation,0.0,options.assignment_occupation_tolerance)) {
                        occupation_clean=false;
                    }
                }
                if (!occupation_clean || !contains_seed || full_count<2) continue;

                MulticentreKind kind=MulticentreKind::Unclassified;
                if (occupied_count==1 && near_occupation(electrons,2.0,options.assignment_occupation_tolerance)) {
                    kind=MulticentreKind::ThreeCentreTwoElectron;
                } else if (occupied_count==2 && near_occupation(electrons,4.0,options.assignment_occupation_tolerance)) {
                    kind=MulticentreKind::ThreeCentreFourElectron;
                } else {
                    continue;
                }

                // Prefer compact three-orbital active spaces with high atomic
                // coverage.  A tiny deterministic index term only breaks exact
                // numerical ties; it is far below the uniqueness threshold.
                score/=3.0;
                score += 0.01*static_cast<double>(full_count);
                score -= 1.0e-10*static_cast<double>(
                    chosen[0]->index + chosen[1]->index + chosen[2]->index);

                AssignmentChoice candidate;
                candidate.kind=kind;
                candidate.orbitals={chosen[0]->index,chosen[1]->index,chosen[2]->index};
                candidate.electron_count=electrons;
                candidate.score=score;
                candidate.mean_coverage=coverage/3.0;

                if (candidate.score>best.score) {
                    second_best_score=best.score;
                    best=candidate;
                } else if (candidate.score>second_best_score) {
                    second_best_score=candidate.score;
                }
            }
        }
    }
    return best;
}

void derive_strict_multicentre_assignments(
    Wavefunction& wf,
    const std::vector<std::uint32_t>& basis_atom,
    const BondAnalysisOptions& options) {
    wf.multicentre_assignments.clear();

    // Separate alpha/beta orbital blocks need spin-resolved active-space
    // electron counting. Until that path is implemented, leave them explicitly
    // unclassified instead of collapsing two spin sets into a false 3c2e/3c4e.
    for (const auto& mo:wf.orbitals) {
        if (mo.spin==Spin::Beta) return;
    }

    std::vector<std::vector<double>> weights;
    weights.reserve(wf.orbitals.size());
    for (const auto& mo:wf.orbitals) {
        weights.push_back(atomic_participation(wf,basis_atom,mo));
    }

    std::set<AtomTriple> triples;
    std::map<AtomTriple,std::set<std::uint32_t>> seeds;
    for (const auto& candidate:wf.multicentre_candidates) {
        if (candidate.atoms.size()<3u || candidate.participation.size()<3u) continue;
        const double top3=candidate.participation[0]+candidate.participation[1]+candidate.participation[2];
        if (top3<options.assignment_triple_coverage) continue;
        const AtomTriple triple=sorted_triple(candidate.atoms[0],candidate.atoms[1],candidate.atoms[2]);
        triples.insert(triple);
        seeds[triple].insert(candidate.orbital_index);
    }

    for (const AtomTriple& triple:triples) {
        const auto active=active_mos_for_triple(wf,weights,triple,options);
        double second=-std::numeric_limits<double>::infinity();
        const AssignmentChoice best=best_assignment_for_triple(
            wf,triple,active,seeds[triple],options,second);
        if (best.kind==MulticentreKind::Unclassified || !std::isfinite(best.score)) continue;

        const double margin=std::isfinite(second)?best.score-second:1.0;
        if (margin<options.assignment_min_uniqueness_margin) continue;

        MulticentreAssignment assignment;
        assignment.kind=best.kind;
        assignment.atoms.assign(triple.begin(),triple.end());
        assignment.orbitals.assign(best.orbitals.begin(),best.orbitals.end());
        std::sort(assignment.orbitals.begin(),assignment.orbitals.end(),[&](std::uint32_t a,std::uint32_t b) {
            return wf.orbitals[a].energy_hartree < wf.orbitals[b].energy_hartree;
        });
        assignment.electron_count=best.electron_count;
        const double uniqueness=std::clamp(margin/0.20,0.0,1.0);
        assignment.confidence=std::clamp(
            0.55 + 0.25*best.mean_coverage + 0.20*uniqueness,0.0,1.0);
        assignment.rationale=
            "unique three-orbital active subspace; >=2 three-atom MOs; "
            ">=2 Mayer-supported pairs; restricted occupations match electron count";
        assignment.provenance=DataProvenance::Derived;
        wf.multicentre_assignments.push_back(std::move(assignment));
    }
}

} // namespace

void derive_bond_and_multicentre_analysis(Wavefunction& wavefunction,
                                          const BondAnalysisOptions& options) {
    wavefunction.bond_orders.clear();
    wavefunction.bond_order_provenance = DataProvenance::Unavailable;
    wavefunction.multicentre_candidates.clear();
    wavefunction.multicentre_assignments.clear();

    const std::size_t n = wavefunction.basis_count;
    if (n == 0 || wavefunction.ao_overlap.size() != n * n ||
        wavefunction.total_density_packed.empty() || wavefunction.spin_density_packed.empty() ||
        wavefunction.total_density_diagnostics.status != NumericalStatus::Available ||
        wavefunction.spin_density_diagnostics.status != NumericalStatus::Available) {
        return;
    }

    const auto basis_atom = basis_atom_map(wavefunction);
    const auto total_density = unpack_symmetric(
        wavefunction.total_density_packed, n, "Total density");
    const auto total_ps = multiply(total_density, wavefunction.ao_overlap, n);

    std::vector<double> spin_ps;
    if (!wavefunction.spin_density_packed.empty()) {
        const auto spin_density = unpack_symmetric(
            wavefunction.spin_density_packed, n, "Spin density");
        spin_ps = multiply(spin_density, wavefunction.ao_overlap, n);
    }

    // Mayer pair index. For an unrestricted wavefunction the total/spin form
    // B_AB = sum[(PS)_mn(PS)_nm + (QS)_mn(QS)_nm] is equivalent to the
    // conventional 2*(alpha + beta) expression and reduces continuously to the
    // ordinary closed-shell formula when Q=0.
    for (std::size_t atom_a = 0; atom_a < wavefunction.atoms.size(); ++atom_a) {
        for (std::size_t atom_b = atom_a + 1; atom_b < wavefunction.atoms.size(); ++atom_b) {
            double order = 0.0;
            for (std::size_t mu = 0; mu < n; ++mu) {
                if (basis_atom[mu] != atom_a) continue;
                for (std::size_t nu = 0; nu < n; ++nu) {
                    if (basis_atom[nu] != atom_b) continue;
                    order += total_ps[mu * n + nu] * total_ps[nu * n + mu];
                    if (!spin_ps.empty()) {
                        order += spin_ps[mu * n + nu] * spin_ps[nu * n + mu];
                    }
                }
            }
            if (std::abs(order) >= options.record_abs_mayer_floor) {
                BondOrderRecord record;
                record.atom_a = static_cast<std::uint32_t>(atom_a);
                record.atom_b = static_cast<std::uint32_t>(atom_b);
                record.mayer_order = order;
                record.provenance = DataProvenance::Derived;
                wavefunction.bond_orders.push_back(record);
            }
        }
    }
    wavefunction.bond_order_provenance = DataProvenance::Derived;

    // Occupied-MO atomic participation. q_A = sum_{mu in A} C_mu (S C)_mu
    // is a Mulliken-like gross MO population and sums to one for an exactly
    // normalized orbital. Absolute q_A values are normalized only for the
    // purpose of flagging spatially shared candidates; the signed values are
    // not reinterpreted as charges or bond orders.
    for (std::size_t orbital_index = 0;
         orbital_index < wavefunction.orbitals.size(); ++orbital_index) {
        const auto& mo = wavefunction.orbitals[orbital_index];
        if (mo.occupation <= 1.0e-4f || mo.coefficients.size() != n) continue;

        const auto atom_weight=atomic_participation(wavefunction,basis_atom,mo);
        if (atom_weight.empty()) continue;

        std::vector<std::pair<double, std::uint32_t>> significant;
        double covered = 0.0;
        for (std::size_t atom = 0; atom < atom_weight.size(); ++atom) {
            if (atom_weight[atom] >= options.multicentre_atom_participation) {
                significant.emplace_back(atom_weight[atom],
                                         static_cast<std::uint32_t>(atom));
                covered += atom_weight[atom];
            }
        }
        if (significant.size() < 3 ||
            covered < options.multicentre_min_covered_weight) {
            continue;
        }

        std::sort(significant.begin(), significant.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        MulticentreCandidate candidate;
        candidate.orbital_index = static_cast<std::uint32_t>(orbital_index);
        candidate.occupation = static_cast<double>(mo.occupation);
        candidate.provenance = DataProvenance::Derived;
        candidate.atoms.reserve(significant.size());
        candidate.participation.reserve(significant.size());
        for (const auto& [weight, atom] : significant) {
            candidate.atoms.push_back(atom);
            candidate.participation.push_back(weight);
        }
        wavefunction.multicentre_candidates.push_back(std::move(candidate));
    }

    derive_strict_multicentre_assignments(wavefunction,basis_atom,options);
}

} // namespace cov
