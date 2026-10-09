#include "cov/orbital_tracking.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include <utility>

namespace cov {
namespace {

struct TrackingBudgetExceeded {};

class TrackingBudget {
public:
    TrackingBudget(const OrbitalTrackingOptions& options,
                   OrbitalTrackingResult& result)
        : options_(options), result_(result), started_(std::chrono::steady_clock::now()) {}

    void stage(const OrbitalTrackingStage stage) {
        stage_ = stage;
        check_clock();
    }

    void charge(const std::size_t units = 1u) {
        if (units > options_.maximum_tracking_work_units - result_.tracking_work_units)
            exhaust(OrbitalTrackingBudgetExhaustion::WorkLimit);
        result_.tracking_work_units += units;
        result_.tracking_work_units_by_stage[static_cast<std::size_t>(stage_)] += units;
        if (result_.tracking_work_units - last_clock_work_ >= 256u) {
            last_clock_work_ = result_.tracking_work_units;
            check_clock();
        }
    }

    void product(const std::size_t rows, const std::size_t columns) {
        const auto remaining = options_.maximum_tracking_work_units - result_.tracking_work_units;
        if (columns != 0u && rows > remaining / columns)
            exhaust(OrbitalTrackingBudgetExhaustion::WorkLimit);
        charge(rows * columns);
    }

    void component_completed() { ++result_.composite_components_completed; }

private:
    [[noreturn]] void exhaust(const OrbitalTrackingBudgetExhaustion reason) {
        result_.tracking_budget_exhausted = true;
        result_.budget_exhaustion_reason = reason;
        result_.budget_exhausted_stage = stage_;
        throw TrackingBudgetExceeded{};
    }
    void check_clock() {
        if (!(options_.maximum_tracking_milliseconds > 0.0) ||
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started_).count() >=
                options_.maximum_tracking_milliseconds)
            exhaust(OrbitalTrackingBudgetExhaustion::TimeLimit);
    }
    const OrbitalTrackingOptions& options_;
    OrbitalTrackingResult& result_;
    OrbitalTrackingStage stage_ = OrbitalTrackingStage::None;
    std::chrono::steady_clock::time_point started_;
    std::size_t last_clock_work_ = 0u;
};

struct Subspace {
    std::vector<std::size_t> members;
    Spin spin = Spin::Alpha;
    double energy = 0.0;
    double occupation = 0.0;
    std::vector<double> descriptor;
};

struct CompositeCandidate {
    std::vector<std::size_t> left_groups;
    std::vector<std::size_t> right_groups;
    Subspace left;
    Subspace right;
    double similarity = 0.0;
    double score = 0.0;

    [[nodiscard]] bool composite() const noexcept {
        return left_groups.size() > 1u || right_groups.size() > 1u;
    }
};

constexpr std::size_t angular_bins = 5u; // s,p,d,f,g+
constexpr std::size_t element_bins = 119u; // Z=0..118

std::size_t descriptor_dimension(const Wavefunction& wf) {
    return (wf.atoms.size() + element_bins) * angular_bins;
}

bool compatible_atoms(const Wavefunction& left, const Wavefunction& right) {
    if (left.atoms.size() != right.atoms.size()) return false;
    for (std::size_t i = 0; i < left.atoms.size(); ++i) {
        if (left.atoms[i].atomic_number != right.atoms[i].atomic_number) {
            return false;
        }
    }
    return true;
}

std::vector<double> orbital_descriptor(const Wavefunction& wf,
                                       const MolecularOrbital& orbital,
                                       TrackingBudget& budget) {
    budget.charge(descriptor_dimension(wf));
    std::vector<double> result(descriptor_dimension(wf), 0.0);
    if (!orbital.chemistry.available) return {};
    for (const auto& contribution : orbital.chemistry.ao_contributions) {
        budget.charge();
        if (contribution.atom_index >= wf.atoms.size() ||
            !std::isfinite(contribution.weight)) continue;
        const std::size_t family = static_cast<std::size_t>(std::clamp(
            contribution.angular_momentum, 0, 4));
        result[static_cast<std::size_t>(contribution.atom_index) * angular_bins +
               family] += std::max(0.0, contribution.weight);
        const int atomic_number = wf.atoms[contribution.atom_index].atomic_number;
        if (atomic_number >= 0 &&
            static_cast<std::size_t>(atomic_number) < element_bins) {
            const std::size_t element_offset = wf.atoms.size() * angular_bins;
            result[element_offset +
                   static_cast<std::size_t>(atomic_number) * angular_bins +
                   family] += std::max(0.0, contribution.weight);
        }
    }
    budget.product(2u, result.size());
    const double norm2 = std::inner_product(
        result.begin(), result.end(), result.begin(), 0.0);
    if (!(norm2 > 1.0e-14)) return {};
    const double inverse = 1.0 / std::sqrt(norm2);
    for (double& value : result) value *= inverse;
    return result;
}

std::vector<Subspace> make_subspaces(const Wavefunction& wf,
                                     const double tolerance) {
    std::vector<std::size_t> order(wf.orbitals.size());
    std::iota(order.begin(), order.end(), 0u);
    std::stable_sort(order.begin(), order.end(), [&](const auto a, const auto b) {
        if (wf.orbitals[a].spin != wf.orbitals[b].spin) {
            return wf.orbitals[a].spin == Spin::Alpha;
        }
        if (wf.orbitals[a].energy_hartree != wf.orbitals[b].energy_hartree) {
            return wf.orbitals[a].energy_hartree < wf.orbitals[b].energy_hartree;
        }
        return a < b;
    });

    std::vector<Subspace> result;
    for (const auto index : order) {
        const auto& orbital = wf.orbitals[index];
        const bool join = !result.empty() &&
            result.back().spin == orbital.spin &&
            std::abs(result.back().energy - orbital.energy_hartree) <= tolerance;
        if (!join) {
            Subspace group;
            group.spin = orbital.spin;
            group.energy = orbital.energy_hartree;
            result.push_back(std::move(group));
        }
        auto& group = result.back();
        const std::size_t previous = group.members.size();
        group.members.push_back(index);
        group.energy = (group.energy * static_cast<double>(previous) +
                        orbital.energy_hartree) /
                       static_cast<double>(previous + 1u);
        group.occupation += static_cast<double>(orbital.occupation);
    }

    return result;
}

void build_descriptors(const Wavefunction& wf, std::vector<Subspace>& groups,
                       TrackingBudget& budget) {
    for (auto& group : groups) {
        budget.charge(descriptor_dimension(wf));
        group.descriptor.assign(descriptor_dimension(wf), 0.0);
        bool available = true;
        for (const auto index : group.members) {
            const auto descriptor = orbital_descriptor(wf, wf.orbitals[index], budget);
            if (descriptor.size() != group.descriptor.size()) {
                available = false;
                break;
            }
            budget.charge(descriptor.size());
            for (std::size_t i = 0; i < descriptor.size(); ++i) {
                group.descriptor[i] += descriptor[i];
            }
        }
        if (!available) {
            group.descriptor.clear();
            continue;
        }
        budget.product(2u, group.descriptor.size());
        const double norm2 = std::inner_product(
            group.descriptor.begin(), group.descriptor.end(),
            group.descriptor.begin(), 0.0);
        if (!(norm2 > 1.0e-14)) {
            group.descriptor.clear();
            continue;
        }
        const double inverse = 1.0 / std::sqrt(norm2);
        for (double& value : group.descriptor) value *= inverse;
    }
}

double cosine_similarity(const std::vector<double>& left,
                         const std::vector<double>& right, TrackingBudget& budget) {
    if (left.empty() || left.size() != right.size()) return 0.0;
    budget.charge(left.size());
    return std::clamp(std::inner_product(
        left.begin(), left.end(), right.begin(), 0.0), 0.0, 1.0);
}

double match_score(const Subspace& left, const Subspace& right,
                   const OrbitalTrackingOptions& options,
                   double& similarity, TrackingBudget& budget) {
    budget.charge();
    if (left.spin != right.spin || left.members.size() != right.members.size()) {
        similarity = 0.0;
        return -std::numeric_limits<double>::infinity();
    }
    const double capacity = 2.0 * static_cast<double>(left.members.size());
    const double dimension = static_cast<double>(left.members.size());
    if (dimension > 0.0 &&
        std::abs(left.occupation - right.occupation) / dimension >
            options.maximum_average_occupation_change) {
        similarity = 0.0;
        return -std::numeric_limits<double>::infinity();
    }
    similarity = cosine_similarity(left.descriptor, right.descriptor, budget);
    const double occupation_difference = capacity > 0.0
        ? std::abs(left.occupation - right.occupation) / capacity : 0.0;
    const double energy_difference = std::abs(left.energy - right.energy) /
        std::max(1.0e-8, options.energy_scale_hartree);
    return similarity -
        options.occupation_weight * occupation_difference -
        options.energy_weight * std::min(1.0, energy_difference);
}

Subspace combine_subspaces(const std::vector<Subspace>& source,
                           const std::vector<std::size_t>& groups, TrackingBudget& budget) {
    Subspace result;
    if (groups.empty()) return result;
    result.spin = source[groups.front()].spin;
    std::size_t dimension = 0u;
    bool descriptor_available = true;
    for (const auto group_index : groups) {
        budget.charge();
        const auto& group = source[group_index];
        if (group.spin != result.spin) return {};
        const std::size_t group_dimension = group.members.size();
        budget.charge(group_dimension);
        result.members.insert(result.members.end(), group.members.begin(),
                              group.members.end());
        result.energy += group.energy * static_cast<double>(group_dimension);
        result.occupation += group.occupation;
        dimension += group_dimension;
        if (!descriptor_available) continue;
        if (group.descriptor.empty()) {
            descriptor_available = false;
            result.descriptor.clear();
            continue;
        }
        if (result.descriptor.empty()) {
            budget.charge(group.descriptor.size());
            result.descriptor.assign(group.descriptor.size(), 0.0);
        }
        if (result.descriptor.size() != group.descriptor.size()) {
            descriptor_available = false;
            result.descriptor.clear();
            continue;
        }
        // A degenerate subspace descriptor is the normalised sum of its
        // canonical members. Weighting by dimension reconstructs the same
        // additive chemical support when several split subspaces are joined.
        budget.charge(group.descriptor.size());
        for (std::size_t i = 0u; i < group.descriptor.size(); ++i) {
            result.descriptor[i] +=
                static_cast<double>(group_dimension) * group.descriptor[i];
        }
    }
    if (dimension == 0u) return {};
    result.energy /= static_cast<double>(dimension);
    if (descriptor_available && !result.descriptor.empty()) {
        budget.product(2u, result.descriptor.size());
        const double norm2 = std::inner_product(
            result.descriptor.begin(), result.descriptor.end(),
            result.descriptor.begin(), 0.0);
        if (norm2 > 1.0e-14) {
            const double inverse = 1.0 / std::sqrt(norm2);
            for (double& value : result.descriptor) value *= inverse;
        } else {
            result.descriptor.clear();
        }
    } else {
        result.descriptor.clear();
    }
    return result;
}

std::vector<std::vector<std::size_t>> local_unions(
    const std::vector<Subspace>& source, const Subspace& anchor,
    const OrbitalTrackingOptions& options, TrackingBudget& budget) {
    budget.charge();
    const std::size_t target_dimension = anchor.members.size();
    if (target_dimension <= 1u || anchor.descriptor.empty() ||
        options.maximum_split_components < 2u ||
        options.maximum_local_component_pool < 2u ||
        options.maximum_composite_candidates_per_anchor == 0u) {
        return {};
    }

    struct PoolMember {
        std::size_t group = 0u;
        double relevance = 0.0;
    };
    std::vector<PoolMember> pool;
    for (std::size_t index = 0u; index < source.size(); ++index) {
        budget.charge();
        const auto& group = source[index];
        if (group.spin != anchor.spin || group.descriptor.empty() ||
            group.members.empty() ||
            group.members.size() >= target_dimension) {
            continue;
        }
        const double energy_difference = std::abs(group.energy - anchor.energy);
        if (energy_difference > options.composite_energy_window_hartree +
                                    options.composite_component_span_hartree) {
            continue;
        }
        const double similarity = cosine_similarity(anchor.descriptor,
                                                    group.descriptor, budget);
        const double anchor_occupation = anchor.occupation /
            static_cast<double>(target_dimension);
        const double group_occupation = group.occupation /
            static_cast<double>(group.members.size());
        const double energy_penalty = std::min(
            1.0, energy_difference /
                     std::max(1.0e-8, options.composite_energy_window_hartree));
        const double occupation_penalty =
            0.5 * std::abs(anchor_occupation - group_occupation);
        pool.push_back({index, similarity - 0.10 * energy_penalty -
                                   0.10 * occupation_penalty});
    }
    std::stable_sort(pool.begin(), pool.end(), [&](const auto& left,
                                                   const auto& right) {
        budget.charge();
        if (left.relevance != right.relevance) {
            return left.relevance > right.relevance;
        }
        return left.group < right.group;
    });
    if (pool.size() > options.maximum_local_component_pool) {
        pool.resize(options.maximum_local_component_pool);
    }
    std::sort(pool.begin(), pool.end(), [&](const auto& left, const auto& right) {
        budget.charge();
        return left.group < right.group;
    });

    // The pool bound makes this exact subset enumeration finite:
    // at most 2^maximum_local_component_pool states per anchor. Dimension,
    // component-count, suffix-capacity and energy-span tests are admissible
    // pruning conditions, so no valid union inside the retained pool is lost.
    std::vector<std::size_t> suffix_dimension(pool.size() + 1u, 0u);
    for (std::size_t index = pool.size(); index-- > 0u;) {
        suffix_dimension[index] = suffix_dimension[index + 1u] +
            source[pool[index].group].members.size();
    }
    struct RankedUnion {
        std::vector<std::size_t> groups;
        double score = 0.0;
    };
    std::vector<RankedUnion> ranked;
    const auto ranked_before = [](const auto& left, const auto& right) {
        if (left.score != right.score) return left.score > right.score;
        return left.groups < right.groups;
    };
    // Every increasing-index subset is visited exactly once. Score it when
    // complete and retain only the exact best K under the existing score and
    // lexicographic tie order; no full subset set or full ranking is needed.
    const auto retain = [&](const std::vector<std::size_t>& groups) {
        const auto combined = combine_subspaces(source, groups, budget);
        if (combined.descriptor.empty() || combined.members.size() !=
                                               target_dimension) return;
        if (std::abs(anchor.energy - combined.energy) >
            options.composite_energy_window_hartree) return;
        const double dimension = static_cast<double>(target_dimension);
        if (std::abs(anchor.occupation - combined.occupation) / dimension >
            options.composite_occupation_window) return;
        double similarity = 0.0;
        const double score = match_score(anchor, combined, options, similarity, budget);
        if (!std::isfinite(score) || score <= 0.0 ||
            similarity < options.minimum_composite_similarity) return;
        RankedUnion item{groups, score};
        const auto position = std::lower_bound(ranked.begin(), ranked.end(),
                                               item, ranked_before);
        if (ranked.size() == options.maximum_composite_candidates_per_anchor &&
            position == ranked.end()) return;
        ranked.insert(position, std::move(item));
        if (ranked.size() > options.maximum_composite_candidates_per_anchor)
            ranked.pop_back();
    };
    std::vector<std::size_t> chosen;
    const auto extend = [&](const auto& self, const std::size_t position,
                            const std::size_t accumulated) -> void {
        budget.charge();
        if (accumulated == target_dimension) {
            if (chosen.size() >= 2u) retain(chosen);
            return;
        }
        if (position >= pool.size() ||
            chosen.size() >= options.maximum_split_components ||
            accumulated + suffix_dimension[position] < target_dimension) {
            return;
        }
        for (std::size_t item = position; item < pool.size(); ++item) {
            budget.charge();
            const auto group_index = pool[item].group;
            const auto& group = source[group_index];
            if (!chosen.empty()) {
                const double first_energy = source[chosen.front()].energy;
                if (group.energy - first_energy >
                    options.composite_component_span_hartree) {
                    break;
                }
            }
            const std::size_t next_dimension =
                accumulated + group.members.size();
            if (next_dimension > target_dimension) continue;
            chosen.push_back(group_index);
            self(self, item + 1u, next_dimension);
            chosen.pop_back();
        }
    };
    extend(extend, 0u, 0u);

    std::vector<std::vector<std::size_t>> result;
    result.reserve(ranked.size());
    for (auto& item : ranked) result.push_back(std::move(item.groups));
    return result;
}

struct OrdinaryPairScore {
    double score = -std::numeric_limits<double>::infinity();
    double similarity = 0.0;
    bool ready = false;
};

using OrdinaryPairScores = std::vector<std::vector<OrdinaryPairScore>>;

std::vector<CompositeCandidate> make_candidates(
    const std::vector<Subspace>& left,
    const std::vector<Subspace>& right,
    const OrbitalTrackingOptions& options,
    OrdinaryPairScores& ordinary_scores, TrackingBudget& budget) {
    std::vector<CompositeCandidate> result;
    // Preserve combine_subspaces' singleton arithmetic (including its second
    // normalization), while constructing each singleton only once per frame.
    std::vector<Subspace> left_singletons, right_singletons;
    std::vector<bool> left_raw_equivalent, right_raw_equivalent;
    left_singletons.reserve(left.size());
    right_singletons.reserve(right.size());
    left_raw_equivalent.reserve(left.size());
    right_raw_equivalent.reserve(right.size());
    const auto same_score_inputs = [](const Subspace& original,
                                      const Subspace& singleton) {
        return singleton.energy == original.energy &&
            singleton.occupation == original.occupation &&
            singleton.descriptor == original.descriptor;
    };
    for (std::size_t i = 0u; i < left.size(); ++i) {
        left_singletons.push_back(combine_subspaces(left, {i}, budget));
        left_raw_equivalent.push_back(same_score_inputs(left[i], left_singletons.back()));
    }
    for (std::size_t j = 0u; j < right.size(); ++j) {
        right_singletons.push_back(combine_subspaces(right, {j}, budget));
        right_raw_equivalent.push_back(same_score_inputs(right[j], right_singletons.back()));
    }
    std::set<std::pair<std::vector<std::size_t>,
                       std::vector<std::size_t>>> unique;
    const auto append = [&](std::vector<std::size_t> left_groups,
                            std::vector<std::size_t> right_groups) {
        budget.charge();
        if (!unique.emplace(left_groups, right_groups).second) return;
        CompositeCandidate candidate;
        candidate.left_groups = std::move(left_groups);
        candidate.right_groups = std::move(right_groups);
        if (candidate.left_groups.size() == 1u)
            budget.charge(left_singletons[candidate.left_groups.front()].descriptor.size());
        if (candidate.right_groups.size() == 1u)
            budget.charge(right_singletons[candidate.right_groups.front()].descriptor.size());
        candidate.left = candidate.left_groups.size() == 1u
            ? left_singletons[candidate.left_groups.front()]
            : combine_subspaces(left, candidate.left_groups, budget);
        candidate.right = candidate.right_groups.size() == 1u
            ? right_singletons[candidate.right_groups.front()]
            : combine_subspaces(right, candidate.right_groups, budget);
        if ((candidate.left_groups.size() > 1u ||
             candidate.right_groups.size() > 1u) &&
            std::abs(candidate.left.energy - candidate.right.energy) >
                options.composite_energy_window_hartree) {
            return;
        }
        if (candidate.composite()) {
            const double dimension =
                static_cast<double>(candidate.left.members.size());
            if (std::abs(candidate.left.occupation -
                         candidate.right.occupation) / dimension >
                options.composite_occupation_window) {
                return;
            }
        }
        candidate.score = match_score(candidate.left, candidate.right, options,
                                      candidate.similarity, budget);
        if (!candidate.composite()) {
            const auto i = candidate.left_groups.front();
            const auto j = candidate.right_groups.front();
            // The final ordinary assignment historically scores the original
            // groups, not their normalized singleton copies. Reuse a score only
            // when all arithmetic inputs agree exactly; otherwise score those
            // original groups lazily below to preserve rounding and tie order.
            if (left_raw_equivalent[i] && right_raw_equivalent[j]) {
                ordinary_scores[i][j] =
                    {candidate.score, candidate.similarity, true};
            }
        }
        const double minimum_similarity = candidate.composite()
            ? options.minimum_composite_similarity
            : options.minimum_similarity;
        if (std::isfinite(candidate.score) && candidate.score > 0.0 &&
            candidate.similarity >= minimum_similarity) {
            result.push_back(std::move(candidate));
        }
    };

    // Keep all ordinary one-to-one candidates available for ambiguity
    // assessment. The final ordinary assignment remains Hungarian.
    for (std::size_t i = 0u; i < left.size(); ++i) {
        for (std::size_t j = 0u; j < right.size(); ++j) {
            budget.charge();
            if (left[i].spin == right[j].spin &&
                left[i].members.size() == right[j].members.size()) {
                // This is the same hard occupation rejection as match_score.
                // Apply it before copying dense descriptors; rejected edges
                // never entered the candidate set or its tie ordering.
                const auto dimension = static_cast<double>(left[i].members.size());
                if (dimension > 0.0 &&
                    std::abs(left[i].occupation - right[j].occupation) / dimension >
                        options.maximum_average_occupation_change) continue;
                append({i}, {j});
            }
        }
    }

    // One intact subspace may correspond to lower-symmetry components with the
    // same spin and total dimension in one bounded energy neighbourhood. A
    // crossing level is allowed to interleave the components. The reverse
    // construction handles recombination on restoring symmetry.
    for (std::size_t i = 0u; i < left.size(); ++i) {
        for (auto groups : local_unions(right, left[i], options, budget)) {
            append({i}, std::move(groups));
        }
    }
    for (std::size_t j = 0u; j < right.size(); ++j) {
        for (auto groups : local_unions(left, right[j], options, budget)) {
            append(std::move(groups), {j});
        }
    }
    return result;
}

bool contains_group(const std::vector<std::size_t>& groups,
                    const std::size_t group) {
    return std::find(groups.begin(), groups.end(), group) != groups.end();
}

double candidate_utility(const CompositeCandidate& candidate) {
    return candidate.score * static_cast<double>(candidate.left.members.size());
}

class DisjointSet {
public:
    explicit DisjointSet(const std::size_t size) : parent_(size), rank_(size, 0u) {
        std::iota(parent_.begin(), parent_.end(), 0u);
    }

    std::size_t find(const std::size_t item) {
        if (parent_[item] != item) parent_[item] = find(parent_[item]);
        return parent_[item];
    }

    void unite(const std::size_t left, const std::size_t right) {
        std::size_t root_left = find(left);
        std::size_t root_right = find(right);
        if (root_left == root_right) return;
        if (rank_[root_left] < rank_[root_right]) std::swap(root_left, root_right);
        parent_[root_right] = root_left;
        if (rank_[root_left] == rank_[root_right]) ++rank_[root_left];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> rank_;
};

struct CompositeSelectionOutcome {
    std::vector<std::size_t> selected;
    std::size_t optimizer_states = 0u;
    std::size_t fallback_components = 0u;
    bool truncated = false;
};

std::vector<std::size_t> conservative_component_selection(
    const std::vector<std::size_t>& component,
    const std::vector<CompositeCandidate>& candidates,
    const OrbitalTrackingOptions& options, TrackingBudget& budget) {
    struct Ranking {
        std::size_t best = std::numeric_limits<std::size_t>::max();
        double best_score = -std::numeric_limits<double>::infinity();
        double second_score = -std::numeric_limits<double>::infinity();
    };
    std::map<std::pair<int, std::size_t>, Ranking> rankings;
    for (const auto candidate_index : component) {
        budget.charge();
        const auto& candidate = candidates[candidate_index];
        const auto update = [&](const int side, const std::size_t group) {
            budget.charge();
            auto& ranking = rankings[{side, group}];
            if (candidate.score > ranking.best_score + 1.0e-12) {
                ranking.second_score = ranking.best_score;
                ranking.best_score = candidate.score;
                ranking.best = candidate_index;
            } else {
                ranking.second_score = std::max(ranking.second_score,
                                                candidate.score);
            }
        };
        for (const auto group : candidate.left_groups) update(0, group);
        for (const auto group : candidate.right_groups) update(1, group);
    }

    std::vector<std::size_t> selected;
    for (const auto candidate_index : component) {
        budget.charge();
        const auto& candidate = candidates[candidate_index];
        bool uniquely_best = true;
        const auto check = [&](const int side, const std::size_t group) {
            budget.charge();
            const auto& ranking = rankings.at({side, group});
            return ranking.best == candidate_index &&
                (!std::isfinite(ranking.second_score) ||
                 ranking.best_score - ranking.second_score >=
                     options.ambiguity_margin);
        };
        for (const auto group : candidate.left_groups) {
            uniquely_best = uniquely_best && check(0, group);
        }
        for (const auto group : candidate.right_groups) {
            uniquely_best = uniquely_best && check(1, group);
        }
        if (uniquely_best) selected.push_back(candidate_index);
    }
    return selected;
}

CompositeSelectionOutcome optimal_composite_selection(
    const std::vector<CompositeCandidate>& candidates,
    const OrbitalTrackingOptions& options, TrackingBudget& budget) {
    CompositeSelectionOutcome outcome;
    std::vector<std::size_t> composites;
    for (std::size_t index = 0u; index < candidates.size(); ++index) {
        budget.charge();
        if (candidates[index].composite()) composites.push_back(index);
    }
    if (composites.empty()) return outcome;

    // Candidates in different conflict components share no source or target
    // group, so their optima are independent. Build those components without
    // an O(candidate^2) pair scan.
    DisjointSet disjoint(composites.size());
    std::map<std::size_t, std::size_t> left_owner;
    std::map<std::size_t, std::size_t> right_owner;
    for (std::size_t position = 0u; position < composites.size(); ++position) {
        budget.charge();
        const auto& candidate = candidates[composites[position]];
        for (const auto group : candidate.left_groups) {
            budget.charge();
            const auto [iterator, inserted] = left_owner.emplace(group, position);
            if (!inserted) disjoint.unite(position, iterator->second);
        }
        for (const auto group : candidate.right_groups) {
            budget.charge();
            const auto [iterator, inserted] = right_owner.emplace(group, position);
            if (!inserted) disjoint.unite(position, iterator->second);
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> components;
    for (std::size_t position = 0u; position < composites.size(); ++position) {
        budget.charge();
        components[disjoint.find(position)].push_back(composites[position]);
    }

    struct DpValue {
        double utility = 0.0;
        std::size_t covered_dimension = 0u;
        std::ptrdiff_t chosen_candidate = -2; // -1 excludes pivot; >=0 selects
        std::size_t pivot_token = 0u;
    };
    const auto better = [](const DpValue& left, const DpValue& right) {
        if (left.utility > right.utility + 1.0e-12) return true;
        if (right.utility > left.utility + 1.0e-12) return false;
        return left.covered_dimension > right.covered_dimension;
    };

    for (auto& [root, component] : components) {
        budget.charge();
        (void)root;
        std::stable_sort(component.begin(), component.end(),
                         [&](const auto left, const auto right) {
            budget.charge();
            return std::tie(candidates[left].left_groups,
                            candidates[left].right_groups) <
                   std::tie(candidates[right].left_groups,
                            candidates[right].right_groups);
        });
        if (component.size() > options.maximum_conflict_component_candidates ||
            options.maximum_optimizer_states == 0u ||
            options.maximum_optimizer_milliseconds <= 0.0) {
            const auto conservative = conservative_component_selection(
                component, candidates, options, budget);
            outcome.selected.insert(outcome.selected.end(), conservative.begin(),
                                    conservative.end());
            ++outcome.fallback_components;
            outcome.truncated = true;
            budget.component_completed();
            continue;
        }
        std::map<std::pair<int, std::size_t>, std::size_t> token_index;
        for (const auto candidate_index : component) {
            budget.charge();
            for (const auto group : candidates[candidate_index].left_groups) {
                token_index.emplace(std::make_pair(0, group), token_index.size());
            }
            for (const auto group : candidates[candidate_index].right_groups) {
                token_index.emplace(std::make_pair(1, group), token_index.size());
            }
        }
        const std::size_t token_words = (token_index.size() + 63u) / 64u;
        budget.product(component.size(), token_words);
        std::vector<std::vector<std::uint64_t>> masks(
            component.size(), std::vector<std::uint64_t>(token_words, 0u));
        for (std::size_t position = 0u; position < component.size(); ++position) {
            budget.charge();
            const auto& candidate = candidates[component[position]];
            const auto set_token = [&](const int side, const std::size_t group) {
                budget.charge();
                const std::size_t bit = token_index.at({side, group});
                masks[position][bit / 64u] |= std::uint64_t{1} << (bit % 64u);
            };
            for (const auto group : candidate.left_groups) set_token(0, group);
            for (const auto group : candidate.right_groups) set_token(1, group);
        }

        const std::size_t candidate_words = (component.size() + 63u) / 64u;
        budget.product(token_index.size(), candidate_words);
        std::vector<std::vector<std::uint64_t>> incidence(
            token_index.size(),
            std::vector<std::uint64_t>(candidate_words, 0u));
        for (std::size_t candidate = 0u; candidate < component.size(); ++candidate) {
            for (std::size_t token = 0u; token < token_index.size(); ++token) {
                budget.charge();
                if ((masks[candidate][token / 64u] &
                     (std::uint64_t{1} << (token % 64u))) != 0u) {
                    incidence[token][candidate / 64u] |=
                        std::uint64_t{1} << (candidate % 64u);
                }
            }
        }
        budget.product(component.size(), candidate_words);
        std::vector<std::vector<std::uint64_t>> conflicts(
            component.size(),
            std::vector<std::uint64_t>(candidate_words, 0u));
        for (std::size_t candidate = 0u; candidate < component.size(); ++candidate) {
            for (std::size_t token = 0u; token < token_index.size(); ++token) {
                budget.charge();
                if ((masks[candidate][token / 64u] &
                     (std::uint64_t{1} << (token % 64u))) == 0u) {
                    continue;
                }
                for (std::size_t word = 0u; word < candidate_words; ++word) {
                    budget.charge();
                    conflicts[candidate][word] |= incidence[token][word];
                }
            }
        }

        // Exact weighted set-packing DP. Branching on the most-contended group
        // removes a whole candidate clique at every step; this is dramatically
        // smaller than include/exclude-by-candidate for split anchors while
        // preserving the global optimum of the conflict component.
        std::map<std::vector<std::uint64_t>, DpValue> memo;
        std::size_t explored_states = 0u;
        bool aborted = false;
        const auto deadline = std::chrono::steady_clock::now() +
            std::chrono::duration<double, std::milli>(
                options.maximum_optimizer_milliseconds);
        const auto solve = [&](const auto& self,
                               const std::vector<std::uint64_t>& available)
                               -> DpValue {
            budget.charge();
            if (aborted) return {};
            budget.charge(available.size());
            bool any = false;
            for (const auto word : available) any = any || word != 0u;
            if (!any) return {};
            if (const auto iterator = memo.find(available);
                iterator != memo.end()) {
                return iterator->second;
            }
            if (explored_states >= options.maximum_optimizer_states ||
                (explored_states % 64u == 0u &&
                 std::chrono::steady_clock::now() >= deadline)) {
                aborted = true;
                return {};
            }
            ++explored_states;

            std::size_t pivot = 0u;
            std::size_t pivot_count = 0u;
            for (std::size_t token = 0u; token < incidence.size(); ++token) {
                budget.charge(candidate_words);
                std::size_t count = 0u;
                for (std::size_t word = 0u; word < candidate_words; ++word) {
                    count += static_cast<std::size_t>(std::popcount(
                        available[word] & incidence[token][word]));
                }
                if (count > pivot_count) {
                    pivot = token;
                    pivot_count = count;
                }
            }
            std::vector<std::uint64_t> next = available;
            budget.charge(candidate_words);
            for (std::size_t word = 0u; word < candidate_words; ++word) {
                next[word] &= ~incidence[pivot][word];
            }
            DpValue best = self(self, next);
            if (aborted) return {};
            best.chosen_candidate = -1;
            best.pivot_token = pivot;

            for (std::size_t candidate = 0u; candidate < component.size();
                 ++candidate) {
                budget.charge();
                const bool incident =
                    (available[candidate / 64u] &
                     incidence[pivot][candidate / 64u] &
                     (std::uint64_t{1} << (candidate % 64u))) != 0u;
                if (!incident) continue;
                next = available;
                budget.product(2u, candidate_words);
                for (std::size_t word = 0u; word < candidate_words; ++word) {
                    next[word] &= ~conflicts[candidate][word];
                }
                DpValue take = self(self, next);
                if (aborted) return {};
                const auto& item = candidates[component[candidate]];
                take.utility += candidate_utility(item);
                take.covered_dimension += item.left.members.size();
                take.chosen_candidate = static_cast<std::ptrdiff_t>(candidate);
                take.pivot_token = pivot;
                if (better(take, best)) best = take;
            }
            memo.emplace(available, best);
            return best;
        };

        std::vector<std::uint64_t> available(candidate_words,
                                             ~std::uint64_t{0});
        if (component.size() % 64u != 0u) {
            available.back() = (std::uint64_t{1} << (component.size() % 64u)) - 1u;
        }
        (void)solve(solve, available);
        outcome.optimizer_states += explored_states;
        if (aborted) {
            const auto conservative = conservative_component_selection(
                component, candidates, options, budget);
            outcome.selected.insert(outcome.selected.end(), conservative.begin(),
                                    conservative.end());
            ++outcome.fallback_components;
            outcome.truncated = true;
        } else {
            while (true) {
                budget.charge(candidate_words);
                const auto iterator = memo.find(available);
                if (iterator == memo.end()) break;
                const auto& decision = iterator->second;
                if (decision.chosen_candidate >= 0) {
                    const std::size_t candidate = static_cast<std::size_t>(
                        decision.chosen_candidate);
                    outcome.selected.push_back(component[candidate]);
                    for (std::size_t word = 0u; word < candidate_words; ++word) {
                        available[word] &= ~conflicts[candidate][word];
                    }
                } else {
                    for (std::size_t word = 0u; word < candidate_words; ++word) {
                        available[word] &=
                            ~incidence[decision.pivot_token][word];
                    }
                }
            }
        }
        budget.component_completed();
    }
    std::sort(outcome.selected.begin(), outcome.selected.end(),
              [&](const auto left, const auto right) {
                  budget.charge();
                  return left < right;
              });
    return outcome;
}

// Maximum-weight one-to-one assignment. Invalid and dummy cells have zero
// weight; every chemically admissible match is positive, so the square padded
// assignment naturally leaves unsupported groups unmatched without a greedy
// first-choice collision. The implementation is the rectangular Hungarian
// primal-dual algorithm applied to a square cost matrix.
std::vector<std::size_t> maximum_weight_assignment(
    const std::vector<std::vector<double>>& weights, TrackingBudget& budget) {
    const std::size_t n = weights.size();
    if (n == 0u) return {};
    std::vector<double> u(n + 1u, 0.0), v(n + 1u, 0.0);
    std::vector<std::size_t> p(n + 1u, 0u), way(n + 1u, 0u);
    for (std::size_t i = 1u; i <= n; ++i) {
        budget.charge(n);
        p[0] = i;
        std::size_t j0 = 0u;
        std::vector<double> minimum(n + 1u,
                                    std::numeric_limits<double>::infinity());
        std::vector<bool> used(n + 1u, false);
        do {
            used[j0] = true;
            const std::size_t i0 = p[j0];
            double delta = std::numeric_limits<double>::infinity();
            std::size_t j1 = 0u;
            for (std::size_t j = 1u; j <= n; ++j) {
                budget.charge();
                if (used[j]) continue;
                const double cost = -weights[i0 - 1u][j - 1u];
                const double current = cost - u[i0] - v[j];
                if (current < minimum[j]) {
                    minimum[j] = current;
                    way[j] = j0;
                }
                if (minimum[j] < delta) {
                    delta = minimum[j];
                    j1 = j;
                }
            }
            for (std::size_t j = 0u; j <= n; ++j) {
                budget.charge();
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minimum[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0u);
        do {
            budget.charge();
            const std::size_t j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0u);
    }
    std::vector<std::size_t> column_for_row(n, n);
    for (std::size_t j = 1u; j <= n; ++j) {
        if (p[j] != 0u) column_for_row[p[j] - 1u] = j - 1u;
    }
    return column_for_row;
}

} // namespace

OrbitalTrackingResult track_orbital_subspaces(
    const Wavefunction& from,
    const Wavefunction& to,
    const OrbitalTrackingOptions& options) {
    OrbitalTrackingResult result;
    TrackingBudget budget(options, result);
    result.atom_mapping_compatible = compatible_atoms(from, to);
    auto left = make_subspaces(from, options.degeneracy_tolerance_hartree);
    auto right = make_subspaces(to, options.degeneracy_tolerance_hartree);
    if (!result.atom_mapping_compatible) {
        for (const auto& item : left) result.unmatched_from.push_back(item.members);
        for (const auto& item : right) result.unmatched_to.push_back(item.members);
        return result;
    }

    try {
        budget.stage(OrbitalTrackingStage::Descriptors);
        build_descriptors(from, left, budget);
        build_descriptors(to, right, budget);
        budget.stage(OrbitalTrackingStage::Candidates);
        // TODO(perf) ✳: Profile candidate generation and final assignment as well as
        // composite search. Per-component DP limits do not bound total tracking time.
        // All search stages share this transaction's budget. Per-component DP
        // limits remain separate; exhausting the shared budget discards partial
        // matching decisions and reports every original membership as unresolved.
        budget.product(left.size(), right.size());
        OrdinaryPairScores ordinary_scores(
            left.size(), std::vector<OrdinaryPairScore>(right.size()));
        const auto candidates = make_candidates(left, right, options, ordinary_scores, budget);
        for (const auto& candidate : candidates) {
            budget.charge();
            if (candidate.composite()) ++result.composite_candidates_considered;
        }

        std::vector<bool> used_left(left.size(), false);
        std::vector<bool> used_right(right.size(), false);
        const auto append_match = [&](const CompositeCandidate& candidate) {
            double alternative = -std::numeric_limits<double>::infinity();
            for (const auto& other : candidates) {
                budget.charge();
                if (&other == &candidate) continue;
                bool conflict = false;
                for (const auto group : candidate.left_groups) {
                    conflict = conflict || contains_group(other.left_groups, group);
                }
                for (const auto group : candidate.right_groups) {
                    conflict = conflict || contains_group(other.right_groups, group);
                }
                if (conflict) alternative = std::max(alternative, other.score);
            }
            OrbitalSubspaceMatch match;
            match.from_members = candidate.left.members;
            match.to_members = candidate.right.members;
            match.similarity = candidate.similarity;
            match.score = candidate.score;
            match.ambiguous = std::isfinite(alternative) &&
                match.score - alternative < options.ambiguity_margin;
            // A split/recombined identity is intentionally reported as one
            // subspace match. Canonical MO indices, coefficients, energies and
            // symmetry labels in both Wavefunctions remain untouched.
            match.source = OrbitalTrackingSource::SMetricChemicalDescriptor;
            result.matches.push_back(std::move(match));
        };

        // Solve each independent composite-candidate conflict component exactly.
        // This removes greedy ordering dependence: a single high-scoring union is
        // rejected when two mutually compatible unions have greater total utility.
        budget.stage(OrbitalTrackingStage::CompositeOptimization);
        const auto composite_selection = optimal_composite_selection(candidates,
                                                                     options, budget);
        result.composite_optimizer_states = composite_selection.optimizer_states;
        result.composite_fallback_components =
            composite_selection.fallback_components;
        result.composite_optimisation_truncated = composite_selection.truncated;
        budget.stage(OrbitalTrackingStage::CompositeReporting);
        for (const auto candidate_index : composite_selection.selected) {
            budget.charge();
            const auto& candidate = candidates[candidate_index];
            append_match(candidate);
            for (const auto group : candidate.left_groups) used_left[group] = true;
            for (const auto group : candidate.right_groups) used_right[group] = true;
            ++result.composite_matches_selected;
        }

        const bool remaining_left =
            std::find(used_left.begin(), used_left.end(), false) != used_left.end();
        const bool remaining_right =
            std::find(used_right.begin(), used_right.end(), false) != used_right.end();
        // No ordinary assignment can remain when composites consumed either side.
        // Keep the original padded dimensions/order otherwise: deleting zero rows
        // or columns can change the historical choice between equal optima.
        if (remaining_left && remaining_right) {
            budget.stage(OrbitalTrackingStage::OrdinaryMatrix);
            const std::size_t padded = std::max(left.size(), right.size());
            budget.product(padded, padded);
            std::vector<std::vector<double>> weights(
                padded, std::vector<double>(padded, 0.0));
            bool any_weight = false;
            for (std::size_t i = 0; i < left.size(); ++i) {
                for (std::size_t j = 0; j < right.size(); ++j) {
                    budget.charge();
                    auto& cached = ordinary_scores[i][j];
                    if (used_left[i] || used_right[j]) {
                        // Consumed groups did not enter the old ordinary score
                        // matrix and must not become ambiguity alternatives now.
                        cached = {};
                        continue;
                    }
                    if (!cached.ready) {
                        cached.score = match_score(left[i], right[j], options,
                                                   cached.similarity, budget);
                        cached.ready = true;
                    }
                    const double similarity = cached.similarity;
                    const double score = cached.score;
                    if (std::isfinite(score) &&
                        similarity >= options.minimum_similarity && score > 0.0) {
                        weights[i][j] = score;
                        any_weight = true;
                    }
                }
            }
            budget.stage(OrbitalTrackingStage::OrdinaryAssignment);
            const auto assignment = any_weight ? maximum_weight_assignment(weights, budget)
                : std::vector<std::size_t>(padded, padded);
            budget.stage(OrbitalTrackingStage::OrdinaryReporting);
            for (std::size_t i = 0u; i < left.size(); ++i) {
                budget.charge();
                if (used_left[i]) continue;
                const std::size_t j = assignment[i];
                if (j >= right.size() || used_right[j] || weights[i][j] <= 0.0) continue;
                double alternative = -std::numeric_limits<double>::infinity();
                for (std::size_t other = 0u; other < right.size(); ++other) {
                    budget.charge();
                    if (other != j) alternative = std::max(alternative, ordinary_scores[i][other].score);
                }
                for (std::size_t other = 0u; other < left.size(); ++other) {
                    budget.charge();
                    if (other != i) alternative = std::max(alternative, ordinary_scores[other][j].score);
                }
                OrbitalSubspaceMatch match;
                match.from_members = left[i].members;
                match.to_members = right[j].members;
                match.similarity = ordinary_scores[i][j].similarity;
                match.score = ordinary_scores[i][j].score;
                match.ambiguous = std::isfinite(alternative) &&
                    match.score - alternative < options.ambiguity_margin;
                match.source = OrbitalTrackingSource::SMetricChemicalDescriptor;
                result.matches.push_back(std::move(match));
                used_left[i] = true;
                used_right[j] = true;
            }
        }
        // A last stage checkpoint also catches elapsed time in the last unpolled
        // batch. No partially reported assignment escapes this transaction.
        budget.stage(OrbitalTrackingStage::OrdinaryReporting);
        for (std::size_t i = 0; i < left.size(); ++i) {
            if (!used_left[i]) result.unmatched_from.push_back(left[i].members);
        }
        for (std::size_t i = 0; i < right.size(); ++i) {
            if (!used_right[i]) result.unmatched_to.push_back(right[i].members);
        }
    } catch (const TrackingBudgetExceeded&) {
        result.matches.clear();
        result.unmatched_from.clear();
        result.unmatched_to.clear();
        result.composite_matches_selected = 0u;
        for (const auto& item : left) result.unresolved_from.push_back(item.members);
        for (const auto& item : right) result.unresolved_to.push_back(item.members);
    }
    return result;
}

const char* orbital_tracking_source_name(
    const OrbitalTrackingSource source) noexcept {
    switch (source) {
        case OrbitalTrackingSource::SMetricChemicalDescriptor:
            return "S-metric chemical descriptor";
        default:
            return "unavailable";
    }
}

} // namespace cov
