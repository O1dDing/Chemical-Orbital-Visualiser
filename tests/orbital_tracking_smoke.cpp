#include "cov/orbital_tracking.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace {

cov::MolecularOrbital orbital(double energy, double occupation,
                              std::size_t atom, double weight) {
    cov::MolecularOrbital result;
    result.energy_hartree = energy;
    result.occupation = static_cast<float>(occupation);
    result.chemistry.available = true;
    cov::OrbitalAOContribution contribution;
    contribution.atom_index = static_cast<std::uint32_t>(atom);
    contribution.principal_n = 2;
    contribution.angular_momentum = 1;
    contribution.weight = weight;
    result.chemistry.ao_contributions.push_back(contribution);
    return result;
}

cov::MolecularOrbital mixed_orbital(
    const double energy,
    const std::vector<std::pair<std::size_t, double>>& contributions) {
    cov::MolecularOrbital result;
    result.energy_hartree = energy;
    result.chemistry.available = true;
    for (const auto& [atom, weight] : contributions) {
        cov::OrbitalAOContribution contribution;
        contribution.atom_index = static_cast<std::uint32_t>(atom);
        contribution.principal_n = 2;
        contribution.angular_momentum = 1;
        contribution.weight = weight;
        result.chemistry.ao_contributions.push_back(contribution);
    }
    return result;
}

cov::Wavefunction crossing_frame(bool crossed) {
    cov::Wavefunction wf;
    wf.atoms.push_back({"C", 6, 0.0, 0.0, 0.0, 6.0});
    wf.atoms.push_back({"O", 8, 2.3, 0.0, 0.0, 8.0});
    if (!crossed) {
        wf.orbitals.push_back(orbital(-0.40, 2.0, 0u, 1.0));
        wf.orbitals.push_back(orbital(-0.20, 0.0, 1u, 1.0));
    } else {
        // Energy ordering crosses, while the atom/angular chemical identity is
        // retained. An energy-index tracker would swap the two identities.
        wf.orbitals.push_back(orbital(-0.45, 0.0, 1u, 1.0));
        wf.orbitals.push_back(orbital(-0.15, 2.0, 0u, 1.0));
    }
    return wf;
}

cov::Wavefunction split_frame(const std::size_t dimension,
                              const bool degenerate,
                              const double centre_energy,
                              const double occupation) {
    cov::Wavefunction wf;
    const std::vector<std::string> symbols = {"C", "O", "N"};
    const std::vector<int> atomic_numbers = {6, 8, 7};
    for (std::size_t atom = 0u; atom < dimension; ++atom) {
        wf.atoms.push_back({symbols[atom], atomic_numbers[atom],
                            2.5 * static_cast<double>(atom), 0.0, 0.0,
                            static_cast<double>(atomic_numbers[atom])});
        const double offset = degenerate
            ? 0.0
            : 0.035 * (static_cast<double>(atom) -
                       0.5 * static_cast<double>(dimension - 1u));
        wf.orbitals.push_back(
            orbital(centre_energy + offset, occupation, atom, 1.0));
    }
    return wf;
}

bool single_union_match(const cov::OrbitalTrackingResult& result,
                        const std::size_t dimension) {
    if (!result.atom_mapping_compatible || result.matches.size() != 1u ||
        !result.unmatched_from.empty() || !result.unmatched_to.empty()) {
        return false;
    }
    std::vector<std::size_t> expected(dimension);
    for (std::size_t index = 0u; index < dimension; ++index) expected[index] = index;
    return result.matches.front().from_members == expected &&
           result.matches.front().to_members == expected;
}

bool same_correspondence(const cov::OrbitalTrackingResult& left,
                         const cov::OrbitalTrackingResult& right) {
    if (left.matches.size() != right.matches.size() ||
        left.unmatched_from != right.unmatched_from ||
        left.unmatched_to != right.unmatched_to ||
        left.atom_mapping_compatible != right.atom_mapping_compatible ||
        left.composite_optimisation_truncated != right.composite_optimisation_truncated)
        return false;
    for (std::size_t i = 0u; i < left.matches.size(); ++i) {
        const auto& a = left.matches[i];
        const auto& b = right.matches[i];
        if (a.from_members != b.from_members || a.to_members != b.to_members ||
            a.similarity != b.similarity || a.score != b.score ||
            a.ambiguous != b.ambiguous || a.source != b.source) return false;
    }
    return !left.tracking_budget_exhausted && !right.tracking_budget_exhausted;
}

bool preserves_unresolved(const cov::OrbitalTrackingResult& result,
                          const std::size_t from_count, const std::size_t to_count) {
    const auto complete = [](const auto& groups, const std::size_t count) {
        std::vector<std::size_t> seen(count, 0u);
        for (const auto& group : groups) {
            if (group.empty()) return false;
            for (const auto member : group) {
                if (member >= count || ++seen[member] != 1u) return false;
            }
        }
        for (const auto visits : seen) if (visits != 1u) return false;
        return true;
    };
    return result.tracking_budget_exhausted && result.matches.empty() &&
        result.unmatched_from.empty() && result.unmatched_to.empty() &&
        result.composite_matches_selected == 0u &&
        complete(result.unresolved_from, from_count) &&
        complete(result.unresolved_to, to_count);
}

bool stage_boundary_checks(const cov::Wavefunction& from,
                           const cov::Wavefunction& to,
                           cov::OrbitalTrackingOptions options = {}) {
    options.maximum_tracking_milliseconds = std::numeric_limits<double>::infinity();
    options.maximum_tracking_work_units = std::numeric_limits<std::size_t>::max();
    const auto reference = cov::track_orbital_subspaces(from, to, options);
    if (reference.tracking_budget_exhausted) return false;
    auto exact_options = options;
    exact_options.maximum_tracking_work_units = reference.tracking_work_units;
    if (!same_correspondence(reference,
            cov::track_orbital_subspaces(from, to, exact_options))) return false;
    std::size_t cumulative = 0u;
    for (std::size_t stage = 1u; stage < reference.tracking_work_units_by_stage.size(); ++stage) {
        const auto work = reference.tracking_work_units_by_stage[stage];
        cumulative += work;
        if (work == 0u) continue;
        auto limited = options;
        limited.maximum_tracking_work_units = cumulative - 1u;
        const auto stopped = cov::track_orbital_subspaces(from, to, limited);
        if (!preserves_unresolved(stopped, from.orbitals.size(), to.orbitals.size()) ||
            stopped.budget_exhausted_stage != static_cast<cov::OrbitalTrackingStage>(stage) ||
            stopped.budget_exhaustion_reason != cov::OrbitalTrackingBudgetExhaustion::WorkLimit ||
            stopped.tracking_work_units > limited.maximum_tracking_work_units) return false;
    }
    return true;
}

} // namespace

int main() {
    const auto first = crossing_frame(false);
    const auto second = crossing_frame(true);
    const auto result = cov::track_orbital_subspaces(first, second);
    if (!result.atom_mapping_compatible || result.matches.size() != 2u ||
        !result.unmatched_from.empty() || !result.unmatched_to.empty()) {
        std::cerr << "crossing trajectory match count regression\n";
        return EXIT_FAILURE;
    }
    bool occupied_identity = false;
    bool virtual_identity = false;
    for (const auto& match : result.matches) {
        if (match.from_members == std::vector<std::size_t>{0u} &&
            match.to_members == std::vector<std::size_t>{1u}) {
            occupied_identity = true;
        }
        if (match.from_members == std::vector<std::size_t>{1u} &&
            match.to_members == std::vector<std::size_t>{0u}) {
            virtual_identity = true;
        }
    }
    if (!occupied_identity || !virtual_identity) {
        std::cerr << "energy crossing incorrectly changed orbital identity\n";
        return EXIT_FAILURE;
    }


    // Two source groups may prefer the same target. The globally weaker source
    // must still receive its chemically valid second choice rather than being
    // dropped by a greedy first-choice assignment.
    cov::Wavefunction competition_from;
    competition_from.atoms = first.atoms;
    competition_from.orbitals = {
        orbital(-0.40, 0.0, 0u, 1.0),
        mixed_orbital(-0.20, {{0u, 0.80}, {1u, 0.60}}),
    };
    cov::Wavefunction competition_to;
    competition_to.atoms = first.atoms;
    competition_to.orbitals = {
        orbital(-0.45, 0.0, 0u, 1.0),
        orbital(-0.15, 0.0, 1u, 1.0),
    };
    const auto competition = cov::track_orbital_subspaces(
        competition_from, competition_to);
    bool second_choice = false;
    for (const auto& match : competition.matches) {
        if (match.from_members == std::vector<std::size_t>{1u} &&
            match.to_members == std::vector<std::size_t>{1u}) {
            second_choice = true;
        }
    }
    if (competition.matches.size() != 2u || !second_choice ||
        !competition.unmatched_from.empty() || !competition.unmatched_to.empty()) {
        std::cerr << "global orbital assignment lost a valid second choice\n";
        return EXIT_FAILURE;
    }

    // Check observable scores and ambiguity, not merely assignment counts.
    // The second source's locally preferred edge is consumed by the stronger
    // source, so its valid second choice must retain the ambiguity warning.
    for (const auto& match : competition.matches) {
        const bool first_source = match.from_members == std::vector<std::size_t>{0u};
        const double expected_score = first_source ? 0.99 : 0.59;
        const double expected_similarity = first_source ? 1.0 : 0.60;
        if (std::abs(match.score - expected_score) > 1e-12 ||
            std::abs(match.similarity - expected_similarity) > 1e-12 ||
            match.ambiguous != !first_source ||
            match.source != cov::OrbitalTrackingSource::SMetricChemicalDescriptor) {
            std::cerr << "ordinary score or ambiguity regression\n";
            return EXIT_FAILURE;
        }
    }

    // Incompatible atom order must still retain sorted, degenerate membership
    // groups on both sides, even when their chemical descriptors are unusable.
    auto incompatible_from = split_frame(2u, true, -0.30, 2.0);
    incompatible_from.orbitals.push_back(orbital(-0.80, 2.0, 0u, 1.0));
    auto incompatible_to = incompatible_from;
    std::swap(incompatible_to.atoms[0], incompatible_to.atoms[1]);
    incompatible_to.orbitals[0].chemistry.available = false;
    const auto incompatible = cov::track_orbital_subspaces(
        incompatible_from, incompatible_to);
    const std::vector<std::vector<std::size_t>> incompatible_groups{{2u}, {0u, 1u}};
    if (incompatible.atom_mapping_compatible || !incompatible.matches.empty() ||
        incompatible.unmatched_from != incompatible_groups ||
        incompatible.unmatched_to != incompatible_groups ||
        incompatible.composite_optimisation_truncated) {
        std::cerr << "incompatible frame membership regression\n";
        return EXIT_FAILURE;
    }

    // Missing chemistry supports no assignment; every original group must be
    // reported unmatched, with no optimizer truncation invented for this case.
    auto unavailable_to = second;
    for (auto& mo : unavailable_to.orbitals) mo.chemistry.available = false;
    const auto unavailable = cov::track_orbital_subspaces(first, unavailable_to);
    const std::vector<std::vector<std::size_t>> singleton_groups{{0u}, {1u}};
    if (!unavailable.atom_mapping_compatible || !unavailable.matches.empty() ||
        unavailable.unmatched_from != singleton_groups ||
        unavailable.unmatched_to != singleton_groups ||
        unavailable.composite_optimisation_truncated) {
        std::cerr << "unsupported ordinary assignment regression\n";
        return EXIT_FAILURE;
    }

    // Symmetry lowering may split an e subspace into two adjacent singlets.
    // Tracking must preserve the two-dimensional identity as a whole instead
    // of silently dropping it because the temporary group sizes differ.
    const auto e_symmetric = split_frame(2u, true, -0.30, 2.0);
    const auto e_split = split_frame(2u, false, -0.29, 2.0);
    const auto e_result = cov::track_orbital_subspaces(e_symmetric, e_split);
    if (!single_union_match(e_result, 2u)) {
        std::cerr << "e to singlet+singlet subspace tracking regression\n";
        return EXIT_FAILURE;
    }

    // A crossing singlet may temporarily interleave the two components in
    // energy order. They remain neighbours in the bounded local energy window
    // and must still be joined without consuming the crossing identity.
    cov::Wavefunction interleaved_from;
    interleaved_from.atoms = split_frame(3u, false, 0.0, 2.0).atoms;
    interleaved_from.orbitals = {
        orbital(-0.30, 2.0, 0u, 1.0),
        orbital(-0.30, 2.0, 1u, 1.0),
        orbital(-0.05, 2.0, 2u, 1.0),
    };
    cov::Wavefunction interleaved_to;
    interleaved_to.atoms = interleaved_from.atoms;
    interleaved_to.orbitals = {
        orbital(-0.32, 2.0, 0u, 1.0),
        orbital(-0.30, 2.0, 2u, 1.0),
        orbital(-0.28, 2.0, 1u, 1.0),
    };
    const auto interleaved = cov::track_orbital_subspaces(
        interleaved_from, interleaved_to);
    bool interleaved_union = false;
    bool crossing_singlet = false;
    for (const auto& match : interleaved.matches) {
        if (match.from_members == std::vector<std::size_t>{0u, 1u} &&
            match.to_members == std::vector<std::size_t>{0u, 2u}) {
            interleaved_union = true;
        }
        if (match.from_members == std::vector<std::size_t>{2u} &&
            match.to_members == std::vector<std::size_t>{1u}) {
            crossing_singlet = true;
        }
    }
    if (!interleaved_union || !crossing_singlet ||
        !interleaved.unmatched_from.empty() ||
        !interleaved.unmatched_to.empty()) {
        std::cerr << "interleaved local split subspace regression\n";
        return EXIT_FAILURE;
    }

    // A complete split/recombine may consume one side while the other has a
    // chemically distinct extra singlet. Preserve that unmatched member in
    // both transition directions without requesting an ordinary assignment.
    auto shorter_split = interleaved_to;
    shorter_split.orbitals.erase(shorter_split.orbitals.begin() + 1);
    const auto partial_lowering = cov::track_orbital_subspaces(
        interleaved_from, shorter_split);
    const auto partial_restoring = cov::track_orbital_subspaces(
        shorter_split, interleaved_from);
    const std::vector<std::vector<std::size_t>> extra_singlet{{2u}};
    const auto partial_union = [](const cov::OrbitalTrackingResult& value) {
        return value.atom_mapping_compatible && value.matches.size() == 1u &&
            value.composite_matches_selected == 1u &&
            value.matches.front().from_members == std::vector<std::size_t>{0u, 1u} &&
            value.matches.front().to_members == std::vector<std::size_t>{0u, 1u};
    };
    if (!partial_union(partial_lowering) || !partial_union(partial_restoring) ||
        partial_lowering.unmatched_from != extra_singlet ||
        !partial_lowering.unmatched_to.empty() ||
        !partial_restoring.unmatched_from.empty() ||
        partial_restoring.unmatched_to != extra_singlet) {
        std::cerr << "partial split/recombine unmatched membership regression\n";
        return EXIT_FAILURE;
    }

    // Missing chemistry in any split member makes the whole temporary union
    // unavailable. A later non-empty member must not resurrect a partial
    // descriptor and make it look like a complete e subspace.
    cov::Wavefunction incomplete_split = e_split;
    incomplete_split.orbitals[0].chemistry.available = false;
    incomplete_split.orbitals[0].chemistry.ao_contributions.clear();
    const auto incomplete = cov::track_orbital_subspaces(
        e_symmetric, incomplete_split);
    if (!incomplete.matches.empty() || incomplete.unmatched_from.size() != 1u ||
        incomplete.unmatched_to.size() != 2u) {
        std::cerr << "partial descriptor incorrectly formed a composite match\n";
        return EXIT_FAILURE;
    }

    // Composite selection is a weighted set-packing problem. The locally best
    // edge L0->{R1,R2} conflicts with both L0->{R0,R1} and L1->{R2,R3}; the
    // latter compatible pair has greater total utility and must win regardless
    // of candidate enumeration order.
    cov::Wavefunction packing_from;
    packing_from.atoms = split_frame(3u, false, 0.0, 2.0).atoms;
    packing_from.orbitals = {
        orbital(-0.30, 2.0, 0u, 1.0),
        orbital(-0.30, 2.0, 0u, 1.0),
        orbital(-0.10, 2.0, 2u, 1.0),
        orbital(-0.10, 2.0, 2u, 1.0),
    };
    cov::Wavefunction packing_to;
    packing_to.atoms = packing_from.atoms;
    packing_to.orbitals = {
        orbital(-0.34, 2.0, 1u, 1.0),
        orbital(-0.29, 2.0, 0u, 1.0),
        orbital(-0.11, 2.0, 0u, 1.0),
        orbital(-0.06, 2.0, 2u, 1.0),
    };
    const auto packing = cov::track_orbital_subspaces(
        packing_from, packing_to);
    bool first_global = false;
    bool second_global = false;
    for (const auto& match : packing.matches) {
        first_global = first_global ||
            (match.from_members == std::vector<std::size_t>{0u, 1u} &&
             match.to_members == std::vector<std::size_t>{0u, 1u});
        second_global = second_global ||
            (match.from_members == std::vector<std::size_t>{2u, 3u} &&
             match.to_members == std::vector<std::size_t>{2u, 3u});
    }
    if (!first_global || !second_global || packing.matches.size() != 2u ||
        packing.composite_matches_selected != 2u ||
        !packing.unmatched_from.empty() || !packing.unmatched_to.empty()) {
        std::cerr << "composite conflict component was not globally optimised\n";
        return EXIT_FAILURE;
    }

    // Exercise both directions of a t2 -> 1+1+1 -> t2 path with singly
    // occupied members. This also guards against losing SOMOs at either the
    // split or recombination step.
    const auto t2_symmetric_a = split_frame(3u, true, -0.22, 1.0);
    const auto t2_split = split_frame(3u, false, -0.20, 1.0);
    const auto t2_symmetric_b = split_frame(3u, true, -0.18, 1.0);
    const auto t2_lowering = cov::track_orbital_subspaces(
        t2_symmetric_a, t2_split);
    const auto t2_restoring = cov::track_orbital_subspaces(
        t2_split, t2_symmetric_b);
    if (!single_union_match(t2_lowering, 3u) ||
        !single_union_match(t2_restoring, 3u)) {
        std::cerr << "t2 split/recombine subspace tracking regression\n";
        return EXIT_FAILURE;
    }

    // The tracker is descriptive only: it must not mutate either canonical MO
    // list while forming a temporary composite subspace.
    if (t2_symmetric_a.orbitals[0].energy_hartree != -0.22 ||
        t2_split.orbitals[0].energy_hartree ==
            t2_split.orbitals[1].energy_hartree) {
        std::cerr << "orbital tracking mutated canonical MO data\n";
        return EXIT_FAILURE;
    }

    // Dense near-degenerate spectra previously enumerated every subset inside
    // the energy window. The bounded local pool and per-anchor candidate cap
    // provide a hard construction limit while preserving the t2 union.
    cov::Wavefunction dense_from;
    for (std::size_t atom = 0u; atom < 12u; ++atom) {
        dense_from.atoms.push_back({"C", 6, static_cast<double>(atom),
                                    0.0, 0.0, 6.0});
    }
    dense_from.orbitals = {
        orbital(-0.20, 0.0, 0u, 1.0),
        orbital(-0.20, 0.0, 1u, 1.0),
        orbital(-0.20, 0.0, 2u, 1.0),
    };
    for (std::size_t index = 0u; index < 36u; ++index) {
        dense_from.orbitals.push_back(orbital(
            -0.18 + 0.0002 * static_cast<double>(index), 0.0,
            (index + 3u) % dense_from.atoms.size(), 1.0));
    }
    cov::Wavefunction dense_to;
    dense_to.atoms = dense_from.atoms;
    for (std::size_t index = 0u; index < dense_from.orbitals.size(); ++index) {
        dense_to.orbitals.push_back(orbital(
            -0.205 + 0.0002 * static_cast<double>(index), 0.0,
            index % dense_to.atoms.size(), 1.0));
    }
    cov::OrbitalTrackingOptions dense_options;
    const auto dense_start = std::chrono::steady_clock::now();
    const auto dense = cov::track_orbital_subspaces(
        dense_from, dense_to, dense_options);
    const double dense_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - dense_start).count();
    if (dense.composite_candidates_considered >
            dense_options.maximum_composite_candidates_per_anchor ||
        dense.composite_matches_selected != 1u || dense_seconds > 2.0) {
        std::cerr << "dense composite candidate bound/performance regression\n";
        return EXIT_FAILURE;
    }

    // A chain of many overlapping split anchors can make one large conflict
    // component even though every anchor is capped. Force the generic
    // component guard and verify that it returns quickly, marks the result as
    // truncated, and uses the same conservative deterministic fallback twice.
    cov::Wavefunction chain_from;
    chain_from.atoms.push_back({"C", 6, 0.0, 0.0, 0.0, 6.0});
    for (std::size_t anchor = 0u; anchor < 24u; ++anchor) {
        const double energy = -0.24 + 0.01 * static_cast<double>(anchor);
        chain_from.orbitals.push_back(orbital(energy, 0.0, 0u, 1.0));
        chain_from.orbitals.push_back(orbital(energy, 0.0, 0u, 1.0));
    }
    cov::Wavefunction chain_to;
    chain_to.atoms = chain_from.atoms;
    for (std::size_t index = 0u; index < 48u; ++index) {
        chain_to.orbitals.push_back(orbital(
            -0.245 + 0.005 * static_cast<double>(index),
            0.0, 0u, 1.0));
    }
    cov::OrbitalTrackingOptions chain_options;
    chain_options.maximum_conflict_component_candidates = 20u;
    const auto chain_start = std::chrono::steady_clock::now();
    const auto chain_first = cov::track_orbital_subspaces(
        chain_from, chain_to, chain_options);
    const auto chain_second = cov::track_orbital_subspaces(
        chain_from, chain_to, chain_options);
    const double chain_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - chain_start).count();
    bool deterministic = chain_first.matches.size() == chain_second.matches.size();
    if (deterministic) {
        for (std::size_t index = 0u; index < chain_first.matches.size(); ++index) {
            deterministic = deterministic &&
                chain_first.matches[index].from_members ==
                    chain_second.matches[index].from_members &&
                chain_first.matches[index].to_members ==
                    chain_second.matches[index].to_members;
        }
    }
    if (!chain_first.composite_optimisation_truncated ||
        chain_first.composite_fallback_components == 0u || !deterministic ||
        chain_seconds > 2.0) {
        std::cerr << "large conflict fallback bound/determinism regression\n";
        return EXIT_FAILURE;
    }
    // Boundaries cover incomplete ordinary/composite candidate sets, matrix
    // preallocation, partially augmented Hungarian solutions and reporting.
    // Even a previously emitted match must roll back if a later phase stops.
    if (!stage_boundary_checks(first, second) ||
        !stage_boundary_checks(competition_from, competition_to) ||
        !stage_boundary_checks(interleaved_from, interleaved_to) ||
        !stage_boundary_checks(packing_from, packing_to)) {
        std::cerr << "end-to-end budget boundary/correspondence regression\n";
        return EXIT_FAILURE;
    }

    auto tied_from = competition_from;
    auto tied_to = competition_to;
    tied_from.orbitals[1] = orbital(-0.20, 0.0, 0u, 1.0);
    tied_to.orbitals[1] = orbital(-0.15, 0.0, 0u, 1.0);
    cov::OrbitalTrackingOptions tie_options;
    tie_options.energy_weight = 0.0;
    const auto tied = cov::track_orbital_subspaces(tied_from, tied_to, tie_options);
    if (!stage_boundary_checks(tied_from, tied_to, tie_options) ||
        tied.matches.size() != 2u || !tied.matches[0].ambiguous ||
        !tied.matches[1].ambiguous ||
        tied.matches[0].from_members != std::vector<std::size_t>{0u} ||
        tied.matches[0].to_members != std::vector<std::size_t>{0u} ||
        tied.matches[1].from_members != std::vector<std::size_t>{1u} ||
        tied.matches[1].to_members != std::vector<std::size_t>{1u}) {
        std::cerr << "equal-score Hungarian order changed under tracking budget\n";
        return EXIT_FAILURE;
    }

    // Isolated equal-energy anchors yield two independent conflict components.
    // Their local DP budgets both suffice, but they share one tracking budget.
    cov::Wavefunction components_from, components_to;
    components_from.atoms = first.atoms;
    components_to.atoms = first.atoms;
    for (std::size_t group = 0u; group < 2u; ++group) {
        const double energy = -5.0 + 3.0 * static_cast<double>(group);
        components_from.orbitals.push_back(orbital(energy, 2.0, group, 1.0));
        components_from.orbitals.push_back(orbital(energy, 2.0, group, 1.0));
        components_to.orbitals.push_back(orbital(energy - 0.02, 2.0, group, 1.0));
        components_to.orbitals.push_back(orbital(energy + 0.02, 2.0, group, 1.0));
    }
    cov::OrbitalTrackingOptions component_options;
    component_options.maximum_tracking_milliseconds = std::numeric_limits<double>::infinity();
    const auto components = cov::track_orbital_subspaces(
        components_from, components_to, component_options);
    const auto optimizer_stage = static_cast<std::size_t>(
        cov::OrbitalTrackingStage::CompositeOptimization);
    std::size_t before_optimization = 0u;
    for (std::size_t stage = 0u; stage < optimizer_stage; ++stage)
        before_optimization += components.tracking_work_units_by_stage[stage];
    component_options.maximum_tracking_work_units = before_optimization +
        components.tracking_work_units_by_stage[optimizer_stage] * 3u / 4u;
    const auto components_stopped = cov::track_orbital_subspaces(
        components_from, components_to, component_options);
    if (components.matches.size() != 2u || components.composite_components_completed != 2u ||
        components_stopped.composite_components_completed != 1u ||
        components_stopped.budget_exhausted_stage != cov::OrbitalTrackingStage::CompositeOptimization ||
        !preserves_unresolved(components_stopped, 4u, 4u)) {
        std::cerr << "independent conflict components did not share total budget\n";
        return EXIT_FAILURE;
    }

    cov::OrbitalTrackingOptions no_time;
    no_time.maximum_tracking_milliseconds = 0.0;
    const auto timed_out = cov::track_orbital_subspaces(first, second, no_time);
    if (!preserves_unresolved(timed_out, 2u, 2u) ||
        timed_out.budget_exhaustion_reason != cov::OrbitalTrackingBudgetExhaustion::TimeLimit ||
        timed_out.budget_exhausted_stage != cov::OrbitalTrackingStage::Descriptors) {
        std::cerr << "tracking deadline did not preserve unresolved memberships\n";
        return EXIT_FAILURE;
    }
    std::cout << "tracking work units: crossing=" << result.tracking_work_units
              << ", packing=" << packing.tracking_work_units
              << ", dense=" << dense.tracking_work_units
              << ", chain=" << chain_first.tracking_work_units << '\n';
    std::cout << "orbital subspace tracking smoke test passed\n";
    return EXIT_SUCCESS;
}
