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

#include "cov/pi_topology.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

namespace cov {
namespace {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;
constexpr double kTiny = 1.0e-12;

Vec3 subtract(const Atom& a, const Atom& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vec3& value) {
    return std::sqrt(std::max(0.0, dot(value, value)));
}

Vec3 normalised(Vec3 value) {
    const double length = norm(value);
    if (length <= kTiny) return {0.0, 0.0, 0.0};
    for (double& component : value) component /= length;
    return value;
}

Mat3 zero_matrix() {
    return {{{0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0}}};
}

Mat3 identity_matrix() {
    Mat3 result = zero_matrix();
    for (std::size_t i = 0; i < 3; ++i) result[i][i] = 1.0;
    return result;
}

void add_outer(Mat3& matrix, const Vec3& direction,
               const double scale = 1.0) {
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            matrix[row][column] +=
                scale * direction[row] * direction[column];
        }
    }
}

struct Eigen3 {
    std::array<double, 3> values{0.0, 0.0, 0.0};
    std::array<Vec3, 3> vectors{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}}};
};

Eigen3 symmetric_eigen(Mat3 matrix) {
    Mat3 vectors = identity_matrix();
    for (int iteration = 0; iteration < 32; ++iteration) {
        std::size_t p = 0;
        std::size_t q = 1;
        double largest = std::abs(matrix[p][q]);
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = row + 1; column < 3; ++column) {
                const double candidate = std::abs(matrix[row][column]);
                if (candidate > largest) {
                    largest = candidate;
                    p = row;
                    q = column;
                }
            }
        }
        if (largest < 1.0e-13) break;
        const double angle = 0.5 * std::atan2(
            2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);

        // Apply a Jacobi similarity rotation while preserving symmetry.
        const double app = matrix[p][p];
        const double aqq = matrix[q][q];
        const double apq = matrix[p][q];
        for (std::size_t k = 0; k < 3; ++k) {
            if (k == p || k == q) continue;
            const double akp = matrix[k][p];
            const double akq = matrix[k][q];
            matrix[k][p] = matrix[p][k] = cosine * akp - sine * akq;
            matrix[k][q] = matrix[q][k] = sine * akp + cosine * akq;
        }
        matrix[p][p] = cosine * cosine * app -
                       2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app +
                       2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = matrix[q][p] = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const double vkp = vectors[k][p];
            const double vkq = vectors[k][q];
            vectors[k][p] = cosine * vkp - sine * vkq;
            vectors[k][q] = sine * vkp + cosine * vkq;
        }
    }

    std::array<std::size_t, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](const auto a, const auto b) {
        return matrix[a][a] < matrix[b][b];
    });
    Eigen3 result;
    for (std::size_t out = 0; out < 3; ++out) {
        const auto source = order[out];
        result.values[out] = matrix[source][source];
        result.vectors[out] = normalised(
            {vectors[0][source], vectors[1][source], vectors[2][source]});
    }
    return result;
}

bool possible_main_group_p_centre(const Atom& atom) {
    const auto z = atom.atomic_number;
    if (z <= 2) return false;
    if ((z >= 21 && z <= 30) || (z >= 39 && z <= 48) ||
        (z >= 57 && z <= 80) || (z >= 89 && z <= 112)) {
        return false;
    }
    return true;
}

bool transition_or_f_block(const Atom& atom) {
    const auto z=atom.atomic_number;
    return (z>=21 && z<=30) || (z>=39 && z<=48) ||
           (z>=57 && z<=80) || (z>=89 && z<=112);
}

double absolute_pair_mayer(const Wavefunction& wavefunction,
                           std::uint32_t atom_a,
                           std::uint32_t atom_b) {
    if (atom_b<atom_a) std::swap(atom_a,atom_b);
    double strongest=0.0;
    for (const auto& record:wavefunction.bond_orders) {
        auto first=record.atom_a;
        auto second=record.atom_b;
        if (second<first) std::swap(first,second);
        if (first==atom_a && second==atom_b) {
            strongest=std::max(strongest,std::abs(record.mayer_order));
        }
    }
    return strongest;
}

struct LocalPSpace {
    std::uint32_t atom = 0;
    Mat3 projector = zero_matrix();
    std::size_t rank = 0;
    double confidence = 0.0;
};

LocalPSpace local_p_space(
    const Wavefunction& wavefunction,
    const std::uint32_t atom,
    const std::vector<std::vector<std::uint32_t>>& neighbours,
    const PiTopologyOptions& options) {
    LocalPSpace result;
    result.atom = atom;
    if (atom >= wavefunction.atoms.size() ||
        !possible_main_group_p_centre(wavefunction.atoms[atom]) ||
        neighbours[atom].empty()) {
        return result;
    }

    Mat3 direction_covariance = zero_matrix();
    std::size_t valid_directions = 0;
    bool has_hydrogen_neighbour=false;
    for (const auto neighbour : neighbours[atom]) {
        if (neighbour >= wavefunction.atoms.size()) continue;
        // A metal above a ligand p system is an interaction partner, not a
        // substituent defining the ligand atom's local hybridisation plane.
        // Including Fe in every Cp-carbon covariance makes a planar eta5 ring
        // look tetrahedral and destroys both ring channels before electronic
        // validation can run.
        if (transition_or_f_block(wavefunction.atoms[neighbour])) continue;
        has_hydrogen_neighbour=has_hydrogen_neighbour ||
            wavefunction.atoms[neighbour].atomic_number==1;
        const Vec3 direction = normalised(subtract(
            wavefunction.atoms[neighbour], wavefunction.atoms[atom]));
        if (norm(direction) <= kTiny) continue;
        add_outer(direction_covariance, direction);
        ++valid_directions;
    }
    if (valid_directions == 0) return result;

    const auto eigen = symmetric_eigen(direction_covariance);
    const double scale = std::max(kTiny, eigen.values[2]);
    const double first = eigen.values[0] / scale;
    const double second = eigen.values[1] / scale;
    if (second <= options.linear_rank_tolerance) {
        result.rank = 2;
        add_outer(result.projector, eigen.vectors[0]);
        add_outer(result.projector, eigen.vectors[1]);
        result.confidence = std::clamp(1.0 - second /
            std::max(kTiny, options.linear_rank_tolerance), 0.0, 1.0);
    } else {
        const double planar_limit=has_hydrogen_neighbour
            ?options.hydrogen_pyramidal_rank_tolerance
            :options.planar_rank_tolerance;
        if (first>planar_limit) return result;
        result.rank = 1;
        add_outer(result.projector, eigen.vectors[0]);
        result.confidence = std::clamp(1.0 - first /
            std::max(kTiny, planar_limit), 0.0, 1.0);
    }
    return result;
}

struct EdgeDirection {
    std::pair<std::uint32_t, std::uint32_t> edge;
    Vec3 direction{0.0, 0.0, 1.0};
    double confidence = 0.0;
};

std::vector<Vec3> common_directions(const LocalPSpace& a,
                                    const LocalPSpace& b,
                                    const PiTopologyOptions& options) {
    Mat3 sum = zero_matrix();
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            sum[row][column] = a.projector[row][column] +
                              b.projector[row][column];
        }
    }
    const auto eigen = symmetric_eigen(sum);
    std::vector<Vec3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        if (eigen.values[i] >= 2.0 - options.common_direction_tolerance) {
            result.push_back(eigen.vectors[i]);
        }
    }
    return result;
}

bool shares_atom(const EdgeDirection& a, const EdgeDirection& b) {
    return a.edge.first == b.edge.first || a.edge.first == b.edge.second ||
           a.edge.second == b.edge.first || a.edge.second == b.edge.second;
}

std::pair<Vec3, double> representative_axis(
    const std::vector<EdgeDirection>& directions,
    const std::vector<std::size_t>& members) {
    Mat3 covariance = zero_matrix();
    for (const auto member : members) {
        add_outer(covariance, directions[member].direction);
    }
    const auto eigen = symmetric_eigen(covariance);
    const double total = std::accumulate(
        eigen.values.begin(), eigen.values.end(), 0.0);
    const double coherence = total > kTiny ? eigen.values[2] / total : 0.0;
    return {eigen.vectors[2], coherence};
}

} // namespace

std::vector<OrientedPiNetwork> infer_oriented_pi_networks(
    const Wavefunction& wavefunction,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& bonded_pairs,
    const PiTopologyOptions& options) {
    std::vector<OrientedPiNetwork> result;
    if (wavefunction.atoms.size() < 2 || bonded_pairs.empty()) return result;

    std::vector<std::vector<std::uint32_t>> neighbours(
        wavefunction.atoms.size());
    std::set<std::pair<std::uint32_t, std::uint32_t>> unique_edges;
    for (auto edge : bonded_pairs) {
        if (edge.first == edge.second || edge.first >= wavefunction.atoms.size() ||
            edge.second >= wavefunction.atoms.size()) {
            continue;
        }
        if (edge.second < edge.first) std::swap(edge.first, edge.second);
        if (!unique_edges.insert(edge).second) continue;
        neighbours[edge.first].push_back(edge.second);
        neighbours[edge.second].push_back(edge.first);
    }

    std::vector<LocalPSpace> local(wavefunction.atoms.size());
    for (std::uint32_t atom = 0; atom < wavefunction.atoms.size(); ++atom) {
        local[atom] = local_p_space(wavefunction, atom, neighbours, options);
    }

    std::vector<EdgeDirection> directions;
    for (const auto& edge : unique_edges) {
        const auto& a = local[edge.first];
        const auto& b = local[edge.second];
        if (a.rank == 0 || b.rank == 0) continue;
        for (const auto& direction : common_directions(a, b, options)) {
            directions.push_back({edge, direction,
                std::min(a.confidence, b.confidence)});
        }
    }
    if (directions.empty()) return result;

    std::vector<bool> visited(directions.size(), false);
    for (std::size_t seed = 0; seed < directions.size(); ++seed) {
        if (visited[seed]) continue;
        std::vector<std::size_t> members;
        std::vector<std::size_t> pending{seed};
        visited[seed] = true;
        while (!pending.empty()) {
            const auto current = pending.back();
            pending.pop_back();
            members.push_back(current);
            for (std::size_t candidate = 0; candidate < directions.size();
                 ++candidate) {
                if (visited[candidate] ||
                    !shares_atom(directions[current], directions[candidate])) {
                    continue;
                }
                if (std::abs(dot(directions[current].direction,
                                 directions[candidate].direction)) <
                    options.adjacent_direction_cosine) {
                    continue;
                }
                visited[candidate] = true;
                pending.push_back(candidate);
            }
        }

        OrientedPiNetwork network;
        std::set<std::uint32_t> atoms;
        std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
        double confidence_sum = 0.0;
        for (const auto member : members) {
            atoms.insert(directions[member].edge.first);
            atoms.insert(directions[member].edge.second);
            edges.insert(directions[member].edge);
            confidence_sum += directions[member].confidence;
        }
        network.atoms.assign(atoms.begin(), atoms.end());
        network.edges.assign(edges.begin(), edges.end());
        const auto [axis, coherence] = representative_axis(directions, members);
        network.representative_direction = axis;
        network.orientation_coherence = coherence;
        network.globally_oriented = coherence >= options.cyclic_global_coherence;
        network.cyclic = network.globally_oriented &&
                         network.edges.size() >= network.atoms.size();
        if (network.cyclic) {
            // Pendant atoms whose p lone pairs happen to align with a cyclic
            // core (for example C6F5) are not members of the ring pi active
            // space.  Extract the graph 2-core only after a cycle has been
            // established; paths, allylic chains and resonance stars retain
            // their terminal atoms unchanged.
            const auto original_edges=network.edges;
            std::map<std::uint32_t,std::set<std::uint32_t>> adjacency;
            for (const auto& edge:network.edges) {
                adjacency[edge.first].insert(edge.second);
                adjacency[edge.second].insert(edge.first);
            }
            std::vector<std::uint32_t> leaf_queue;
            for (const auto& [atom,edge_neighbours]:adjacency) {
                if (edge_neighbours.size()<2u) leaf_queue.push_back(atom);
            }
            while (!leaf_queue.empty()) {
                const auto atom=leaf_queue.back();
                leaf_queue.pop_back();
                const auto found=adjacency.find(atom);
                if (found==adjacency.end() || found->second.size()>=2u) {
                    continue;
                }
                const auto edge_neighbours=found->second;
                adjacency.erase(found);
                for (const auto neighbour:edge_neighbours) {
                    const auto other=adjacency.find(neighbour);
                    if (other==adjacency.end()) continue;
                    other->second.erase(atom);
                    if (other->second.size()<2u) leaf_queue.push_back(neighbour);
                }
            }
            if (adjacency.size()>=3u) {
                std::set<std::uint32_t> core_atoms;
                for (const auto& [atom,edge_neighbours]:adjacency) {
                    (void)edge_neighbours;
                    core_atoms.insert(atom);
                }
                // The graph 2-core deliberately removes aligned pendant lone
                // pairs, but a genuinely conjugated exocyclic multiple bond
                // (fulvene, quinones, exocyclic imines, ...) belongs to the
                // same pi active space. Reattach only tails supported by a
                // strong direct electronic bond. Weak C--F p/lone-pair
                // coupling therefore remains outside aryl pi cores.
                bool extended=true;
                while (extended) {
                    extended=false;
                    for (const auto& edge:original_edges) {
                        const bool first_in=core_atoms.count(edge.first)!=0u;
                        const bool second_in=core_atoms.count(edge.second)!=0u;
                        if (first_in==second_in ||
                            absolute_pair_mayer(wavefunction,edge.first,
                                                edge.second)<0.80) {
                            continue;
                        }
                        core_atoms.insert(first_in?edge.second:edge.first);
                        extended=true;
                    }
                }
                std::vector<std::pair<std::uint32_t,std::uint32_t>> core_edges;
                for (const auto& edge:original_edges) {
                    if (core_atoms.count(edge.first)!=0u &&
                        core_atoms.count(edge.second)!=0u) {
                        core_edges.push_back(edge);
                    }
                }
                if (core_edges.size()>=core_atoms.size()) {
                    network.atoms.assign(core_atoms.begin(),core_atoms.end());
                    network.edges=std::move(core_edges);
                }
            }
        }
        network.confidence = std::clamp(
            0.55 * coherence + 0.45 * confidence_sum /
                static_cast<double>(members.size()), 0.0, 1.0);
        result.push_back(std::move(network));
    }

    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) {
        if (a.atoms != b.atoms) return a.atoms < b.atoms;
        return a.representative_direction < b.representative_direction;
    });
    return result;
}

} // namespace cov
