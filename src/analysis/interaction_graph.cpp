#include "cov/interaction_graph.hpp"

#include "cov/molecule_style.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace cov {
namespace {

constexpr double kRadiansToDegrees = 57.2957795130823208768;

using AtomPair = std::pair<std::uint32_t, std::uint32_t>;

AtomPair ordered_pair(std::uint32_t a, std::uint32_t b) noexcept {
    if (b < a) std::swap(a, b);
    return {a, b};
}

double distance_bohr(const Atom& a, const Atom& b) noexcept {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double radius_sum_bohr(const Atom& a, const Atom& b) noexcept {
    return (covalent_radius_angstrom(a.atomic_number) +
            covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
}

bool is_d_block(const int z) noexcept {
    return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) ||
           (z >= 72 && z <= 80) || (z >= 104 && z <= 112);
}

bool is_f_block(const int z) noexcept {
    return (z >= 57 && z <= 71) || (z >= 89 && z <= 103);
}

bool is_electropositive_main_group_metal(const int z) noexcept {
    switch (z) {
        case 3: case 4: case 11: case 12: case 13:
        case 19: case 20: case 31: case 37: case 38:
        case 49: case 50: case 55: case 56: case 81:
        case 82: case 83: case 87: case 88:
            return true;
        default:
            return false;
    }
}

bool is_coordination_centre_element(const int z) noexcept {
    return is_d_block(z) || is_f_block(z) ||
           is_electropositive_main_group_metal(z);
}

bool is_coordination_donor_element(const int z) noexcept {
    switch (z) {
        // Hydrogen-bearing structural links remain ordinary connectivity;
        // multicentre H semantics are carried by their independent hyperedge.
        case 5:  case 6:  case 7:  case 8:  case 9:
        case 14: case 15: case 16: case 17:
        case 32: case 33: case 34: case 35:
        case 50: case 51: case 52: case 53:
            return true;
        default:
            return false;
    }
}

bool is_hydrogen_bond_donor_element(const int z) noexcept {
    return z == 7 || z == 8 || z == 9 || z == 16 || z == 17;
}

bool is_hydrogen_bond_acceptor_element(const int z) noexcept {
    return z == 7 || z == 8 || z == 9 || z == 15 || z == 16 ||
           z == 17 || z == 33 || z == 34 || z == 35 || z == 53;
}

double mayer_for_pair(const std::map<AtomPair, const BondOrderRecord*>& records,
                      const std::uint32_t a,
                      const std::uint32_t b,
                      DataProvenance& provenance) noexcept {
    const auto found = records.find(ordered_pair(a, b));
    if (found == records.end()) {
        provenance = DataProvenance::Unavailable;
        return 0.0;
    }
    provenance = found->second->provenance;
    return found->second->mayer_order;
}

bool opposite_atomic_charge_pair(const InteractionEvidence& evidence,
                                 const std::size_t atom_count,
                                 const std::uint32_t a,
                                 const std::uint32_t b,
                                 const InteractionGraphOptions& options) noexcept {
    if (!evidence.has_atomic_charges(atom_count) ||
        a >= evidence.atomic_partial_charges.size() ||
        b >= evidence.atomic_partial_charges.size()) {
        return false;
    }
    const double qa = evidence.atomic_partial_charges[a];
    const double qb = evidence.atomic_partial_charges[b];
    return qa * qb < 0.0 && std::abs(qa) >= options.ionic_atomic_charge_floor &&
           std::abs(qb) >= options.ionic_atomic_charge_floor;
}

InteractionEdge make_edge(const Wavefunction& wf,
                          const std::map<AtomPair, const BondOrderRecord*>& records,
                          const std::uint32_t a,
                          const std::uint32_t b,
                          const InteractionKind kind,
                          const InteractionStrength strength,
                          const double confidence) {
    InteractionEdge edge;
    const AtomPair pair = ordered_pair(a, b);
    edge.atom_a = pair.first;
    edge.atom_b = pair.second;
    edge.kind = kind;
    edge.strength = strength;
    edge.distance_bohr = distance_bohr(wf.atoms[edge.atom_a], wf.atoms[edge.atom_b]);
    const double radii = radius_sum_bohr(wf.atoms[edge.atom_a], wf.atoms[edge.atom_b]);
    edge.covalent_radius_ratio = edge.distance_bohr / std::max(1.0e-12, radii);
    edge.mayer_order = mayer_for_pair(records, edge.atom_a, edge.atom_b,
                                      edge.electronic_provenance);
    edge.confidence = std::clamp(confidence, 0.0, 1.0);
    return edge;
}

using EdgeIndicesByPair = std::map<AtomPair, std::vector<std::uint32_t>>;

bool has_pair(const EdgeIndicesByPair& edge_by_pair,
              const std::uint32_t a,
              const std::uint32_t b) {
    return edge_by_pair.find(ordered_pair(a, b)) != edge_by_pair.end();
}

void add_or_replace_strong_edge(
    InteractionGraph& graph,
    EdgeIndicesByPair& edge_by_pair,
    InteractionEdge edge) {
    const AtomPair pair = ordered_pair(edge.atom_a, edge.atom_b);
    const auto found = edge_by_pair.find(pair);
    if (found == edge_by_pair.end()) {
        edge_by_pair[pair].push_back(
            static_cast<std::uint32_t>(graph.edges.size()));
        graph.edges.push_back(std::move(edge));
        return;
    }

    // Repeated evidence for the same semantic layer updates that layer rather
    // than drawing it twice.
    for (const auto index : found->second) {
        if (index < graph.edges.size() && graph.edges[index].kind == edge.kind) {
            graph.edges[index] = std::move(edge);
            return;
        }
    }

    // Strong structural connectivity and higher-level support are independent
    // semantic layers. A 3c/4c or cage interpretation must not erase the
    // ordinary molecular skeleton; keeping both also lets the UI filter the
    // support overlay without making the molecule fall apart.
    for (const auto index : found->second) {
        if (index < graph.edges.size() &&
            graph.edges[index].strength == InteractionStrength::WeakContact) {
            graph.edges[index] = std::move(edge);
            return;
        }
    }
    found->second.push_back(static_cast<std::uint32_t>(graph.edges.size()));
    graph.edges.push_back(std::move(edge));
}

bool donor_is_shielded(
    const Wavefunction& wf,
    const std::uint32_t centre,
    const std::uint32_t donor,
    const std::vector<std::vector<std::uint32_t>>& covalent_adjacency,
    const double margin_angstrom) {
    const double centre_donor = distance_bohr(wf.atoms[centre], wf.atoms[donor]);
    const double margin_bohr = margin_angstrom * kAngstromToBohr;
    for (const std::uint32_t neighbour : covalent_adjacency[donor]) {
        if (neighbour == centre) continue;
        const double centre_neighbour = distance_bohr(wf.atoms[centre], wf.atoms[neighbour]);
        if (centre_neighbour + margin_bohr < centre_donor) return true;
    }
    return false;
}

FragmentAnalysis derive_fragments(const std::size_t atom_count,
                                  const std::vector<InteractionEdge>& edges,
                                  const InteractionEvidence& evidence) {
    FragmentAnalysis analysis;
    analysis.atom_to_fragment.assign(atom_count, FragmentAnalysis::no_fragment);

    std::vector<std::vector<std::uint32_t>> adjacency(atom_count);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        if (edge.strength != InteractionStrength::StrongConnectivity ||
            edge.atom_a >= atom_count || edge.atom_b >= atom_count) {
            continue;
        }
        adjacency[edge.atom_a].push_back(edge.atom_b);
        adjacency[edge.atom_b].push_back(edge.atom_a);
    }

    for (std::uint32_t seed = 0; seed < atom_count; ++seed) {
        if (analysis.atom_to_fragment[seed] != FragmentAnalysis::no_fragment) continue;
        MolecularFragment fragment;
        fragment.index = static_cast<std::uint32_t>(analysis.fragments.size());

        std::queue<std::uint32_t> queue;
        queue.push(seed);
        analysis.atom_to_fragment[seed] = fragment.index;
        while (!queue.empty()) {
            const std::uint32_t atom = queue.front();
            queue.pop();
            fragment.atoms.push_back(atom);
            for (const std::uint32_t next : adjacency[atom]) {
                if (analysis.atom_to_fragment[next] != FragmentAnalysis::no_fragment) continue;
                analysis.atom_to_fragment[next] = fragment.index;
                queue.push(next);
            }
        }
        std::sort(fragment.atoms.begin(), fragment.atoms.end());
        analysis.fragments.push_back(std::move(fragment));
    }

    for (std::uint32_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
        const auto& edge = edges[edge_index];
        if (edge.atom_a >= atom_count || edge.atom_b >= atom_count) continue;
        const std::uint32_t fa = analysis.atom_to_fragment[edge.atom_a];
        const std::uint32_t fb = analysis.atom_to_fragment[edge.atom_b];
        if (edge.strength == InteractionStrength::StrongConnectivity && fa == fb &&
            fa < analysis.fragments.size()) {
            analysis.fragments[fa].strong_edge_indices.push_back(edge_index);
        } else if (edge.strength == InteractionStrength::WeakContact && fa != fb) {
            analysis.weak_interfragment_edge_indices.push_back(edge_index);
        }
    }

    if (evidence.has_atomic_charges(atom_count)) {
        for (auto& fragment : analysis.fragments) {
            fragment.evidenced_charge = 0.0;
            for (const std::uint32_t atom : fragment.atoms) {
                fragment.evidenced_charge += evidence.atomic_partial_charges[atom];
            }
            fragment.charge_provenance = evidence.atomic_charge_provenance;
        }
    }
    return analysis;
}

double angle_degrees(const Atom& a, const Atom& vertex, const Atom& b) noexcept {
    const std::array<double, 3> va{a.x - vertex.x, a.y - vertex.y, a.z - vertex.z};
    const std::array<double, 3> vb{b.x - vertex.x, b.y - vertex.y, b.z - vertex.z};
    const double na = std::sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
    const double nb = std::sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);
    if (!(na > 1.0e-12) || !(nb > 1.0e-12)) return 0.0;
    const double cosine = std::clamp(
        (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) / (na * nb),
        -1.0, 1.0);
    return std::acos(cosine) * kRadiansToDegrees;
}

bool hydrogen_bond_geometry(
    const Wavefunction& wf,
    const std::uint32_t a,
    const std::uint32_t b,
    const std::vector<std::vector<std::uint32_t>>& strong_adjacency,
    const InteractionGraphOptions& options) {
    const std::uint32_t hydrogen = wf.atoms[a].atomic_number == 1 ? a : b;
    const std::uint32_t acceptor = hydrogen == a ? b : a;
    if (wf.atoms[hydrogen].atomic_number != 1 ||
        !is_hydrogen_bond_acceptor_element(wf.atoms[acceptor].atomic_number)) {
        return false;
    }
    const double h_acceptor_angstrom =
        distance_bohr(wf.atoms[hydrogen], wf.atoms[acceptor]) / kAngstromToBohr;
    if (h_acceptor_angstrom > options.hydrogen_bond_max_h_acceptor_angstrom) return false;

    for (const std::uint32_t donor : strong_adjacency[hydrogen]) {
        if (donor == acceptor ||
            !is_hydrogen_bond_donor_element(wf.atoms[donor].atomic_number)) {
            continue;
        }
        if (angle_degrees(wf.atoms[donor], wf.atoms[hydrogen],
                          wf.atoms[acceptor]) >= options.hydrogen_bond_min_angle_degrees) {
            return true;
        }
    }
    return false;
}

struct WeakCandidate {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    std::uint32_t fragment_a = 0;
    std::uint32_t fragment_b = 0;
    double surface_gap_angstrom = 0.0;
    bool hydrogen_bond = false;
};

struct ValidatedMulticentre {
    const MulticentreAssignment* assignment = nullptr;
    std::set<AtomPair> support_pairs;
    std::set<AtomPair> member_pairs;
};

ValidatedMulticentre validate_multicentre_assignment(
    const Wavefunction& wf,
    const MulticentreAssignment& assignment,
    const InteractionGraphOptions& options) {
    ValidatedMulticentre validated;
    if (assignment.kind == MulticentreKind::Unclassified ||
        assignment.atoms.size() < 3u) {
        return validated;
    }
    std::vector<std::uint32_t> atoms = assignment.atoms;
    std::sort(atoms.begin(), atoms.end());
    if (atoms.back() >= wf.atoms.size() ||
        std::adjacent_find(atoms.begin(), atoms.end()) != atoms.end()) {
        return validated;
    }

    struct Candidate {
        AtomPair pair;
        double distance = 0.0;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < atoms.size(); ++j) {
            const AtomPair pair = ordered_pair(atoms[i], atoms[j]);
            validated.member_pairs.insert(pair);
            const double distance = distance_bohr(
                wf.atoms[pair.first], wf.atoms[pair.second]);
            if (distance <= options.multicentre_distance_factor *
                                radius_sum_bohr(wf.atoms[pair.first],
                                                wf.atoms[pair.second])) {
                candidates.push_back({pair, distance});
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a,
                                                       const Candidate& b) {
        if (a.distance != b.distance) return a.distance < b.distance;
        return a.pair < b.pair;
    });

    // Use a distance minimum-spanning tree to identify adjacent members, then
    // restore all symmetry-equivalent edges within a small distance tolerance.
    // Thus linear F-H-F and B-H-B retain two supports, while equilateral H3+
    // retains all three rather than acquiring an arbitrary broken-symmetry tree.
    std::map<std::uint32_t, std::uint32_t> parent;
    for (const std::uint32_t atom : atoms) parent[atom] = atom;
    const auto find_root = [&](std::uint32_t atom, const auto& self) -> std::uint32_t {
        const std::uint32_t p = parent[atom];
        if (p == atom) return atom;
        parent[atom] = self(p, self);
        return parent[atom];
    };

    double longest_tree_edge = 0.0;
    std::size_t tree_edges = 0;
    for (const auto& candidate : candidates) {
        const std::uint32_t ra = find_root(candidate.pair.first, find_root);
        const std::uint32_t rb = find_root(candidate.pair.second, find_root);
        if (ra == rb) continue;
        parent[rb] = ra;
        longest_tree_edge = std::max(longest_tree_edge, candidate.distance);
        ++tree_edges;
        if (tree_edges + 1u == atoms.size()) break;
    }
    if (tree_edges + 1u != atoms.size()) {
        validated.member_pairs.clear();
        return validated;
    }

    constexpr double symmetry_distance_tolerance = 1.05;
    for (const auto& candidate : candidates) {
        if (candidate.distance <= symmetry_distance_tolerance * longest_tree_edge) {
            validated.support_pairs.insert(candidate.pair);
        }
    }
    validated.assignment = &assignment;
    return validated;
}

} // namespace

bool InteractionEvidence::has_atomic_charges(const std::size_t atom_count) const noexcept {
    if (atomic_charge_provenance == DataProvenance::Unavailable ||
        atomic_partial_charges.size() != atom_count) {
        return false;
    }
    return std::all_of(atomic_partial_charges.begin(), atomic_partial_charges.end(),
                       [](const double charge) { return std::isfinite(charge); });
}

bool interaction_merges_fragments(const InteractionKind kind) noexcept {
    return kind == InteractionKind::CovalentConnectivity ||
           kind == InteractionKind::CoordinationContact ||
           kind == InteractionKind::MulticentreSupport ||
           kind == InteractionKind::PolyhedralCageSupport;
}

bool interaction_is_ordinary_bond(const InteractionKind kind) noexcept {
    return kind == InteractionKind::CovalentConnectivity;
}

const char* interaction_kind_name(const InteractionKind kind) noexcept {
    switch (kind) {
        case InteractionKind::CovalentConnectivity: return "covalent-connectivity";
        case InteractionKind::CoordinationContact: return "coordination-contact";
        case InteractionKind::MulticentreSupport: return "multicentre-support";
        case InteractionKind::PolyhedralCageSupport: return "polyhedral-cage-support";
        case InteractionKind::HydrogenBond: return "hydrogen-bond";
        case InteractionKind::NoncovalentContact: return "noncovalent-contact";
        case InteractionKind::IonicContact: return "ionic-contact";
        default: return "ambiguous-contact";
    }
}

std::vector<PolyhedralCage> analyse_polyhedral_cages(const Wavefunction& wf) {
    std::vector<PolyhedralCage> result;
    const std::size_t n = wf.atoms.size();
    if (n < 12u) return result;

    std::vector<std::vector<std::uint32_t>> neighbours(n);
    std::vector<bool> heavy(n, false);
    for (std::size_t i = 0u; i < n; ++i) {
        heavy[i] = wf.atoms[i].atomic_number > 1;
    }
    constexpr double nearest_neighbour_radius_factor = 1.26;
    for (std::uint32_t a = 0u; a < n; ++a) {
        if (!heavy[a]) continue;
        for (std::uint32_t b = a + 1u; b < n; ++b) {
            if (!heavy[b]) continue;
            const double radii = radius_sum_bohr(wf.atoms[a], wf.atoms[b]);
            if (distance_bohr(wf.atoms[a], wf.atoms[b]) <=
                nearest_neighbour_radius_factor * radii) {
                neighbours[a].push_back(b);
                neighbours[b].push_back(a);
            }
        }
    }

    // Terminal/exo atoms cannot belong to an icosahedral cage. Repeatedly
    // remove vertices with fewer than four neighbours before component fits.
    std::vector<bool> active = heavy;
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<std::uint32_t> remove;
        for (std::uint32_t atom = 0u; atom < n; ++atom) {
            if (!active[atom]) continue;
            const std::size_t degree = static_cast<std::size_t>(std::count_if(
                neighbours[atom].begin(), neighbours[atom].end(),
                [&](const auto other) { return active[other]; }));
            if (degree < 4u) remove.push_back(atom);
        }
        for (const auto atom : remove) {
            active[atom] = false;
            changed = true;
        }
    }

    std::vector<bool> visited(n, false);
    for (std::uint32_t seed = 0u; seed < n; ++seed) {
        if (!active[seed] || visited[seed]) continue;
        std::vector<std::uint32_t> component;
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        visited[seed] = true;
        while (!queue.empty()) {
            const auto atom = queue.front();
            queue.pop();
            component.push_back(atom);
            for (const auto other : neighbours[atom]) {
                if (active[other] && !visited[other]) {
                    visited[other] = true;
                    queue.push(other);
                }
            }
        }
        if (component.size() != 12u) continue;
        std::sort(component.begin(), component.end());
        const std::set<std::uint32_t> members(component.begin(), component.end());

        PolyhedralCage cage;
        cage.atoms = component;
        bool topology_ok = true;
        std::vector<double> edge_lengths;
        for (const auto atom : component) {
            std::vector<std::uint32_t> local;
            for (const auto other : neighbours[atom]) {
                if (members.count(other) != 0u) local.push_back(other);
            }
            if (local.size() != 5u) {
                topology_ok = false;
                break;
            }
            for (const auto other : local) {
                if (atom >= other) continue;
                std::size_t common = 0u;
                for (const auto candidate : local) {
                    if (std::find(neighbours[other].begin(), neighbours[other].end(),
                                  candidate) != neighbours[other].end()) {
                        ++common;
                    }
                }
                if (common != 2u) {
                    topology_ok = false;
                    break;
                }
                cage.support_edges.push_back({atom, other});
                edge_lengths.push_back(distance_bohr(wf.atoms[atom], wf.atoms[other]));
            }
            if (!topology_ok) break;
        }
        if (!topology_ok || cage.support_edges.size() != 30u) continue;

        std::size_t electronically_supported = 0u;
        double electronic_sum = 0.0;
        if (wf.bond_order_provenance != DataProvenance::Unavailable) {
            for (const auto& pair : cage.support_edges) {
                const auto found = std::find_if(
                    wf.bond_orders.begin(), wf.bond_orders.end(),
                    [&](const auto& record) {
                        return ordered_pair(record.atom_a, record.atom_b) ==
                               ordered_pair(pair[0], pair[1]);
                    });
                if (found != wf.bond_orders.end() &&
                    found->provenance != DataProvenance::Unavailable &&
                    found->mayer_order >= 0.015) {
                    ++electronically_supported;
                    electronic_sum += found->mayer_order;
                }
            }
            if (electronically_supported < 24u ||
                electronic_sum / 30.0 < 0.025) {
                continue;
            }
        }

        std::array<double, 3> centroid{0.0, 0.0, 0.0};
        for (const auto atom : component) {
            centroid[0] += wf.atoms[atom].x;
            centroid[1] += wf.atoms[atom].y;
            centroid[2] += wf.atoms[atom].z;
        }
        for (auto& value : centroid) value /= 12.0;
        std::vector<double> radial;
        radial.reserve(12u);
        for (const auto atom : component) {
            const double dx = wf.atoms[atom].x - centroid[0];
            const double dy = wf.atoms[atom].y - centroid[1];
            const double dz = wf.atoms[atom].z - centroid[2];
            radial.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
        }
        const auto coefficient_of_variation = [](const std::vector<double>& values) {
            const double mean = std::accumulate(values.begin(), values.end(), 0.0) /
                                static_cast<double>(values.size());
            double variance = 0.0;
            for (const auto value : values) {
                const double delta = value - mean;
                variance += delta * delta;
            }
            variance /= static_cast<double>(values.size());
            return std::sqrt(variance) / std::max(1.0e-12, mean);
        };
        const double edge_cv = coefficient_of_variation(edge_lengths);
        const double radial_cv = coefficient_of_variation(radial);
        if (edge_cv > 0.16 || radial_cv > 0.16) continue;
        cage.confidence = std::clamp(1.0 - 2.2 * edge_cv - 1.8 * radial_cv,
                                     0.0, 1.0);
        if (wf.bond_order_provenance == DataProvenance::Unavailable) {
            cage.confidence = std::min(cage.confidence, 0.62);
        }
        cage.provenance = DataProvenance::Derived;
        result.push_back(std::move(cage));
    }
    return result;
}

InteractionGraph build_interaction_graph(const Wavefunction& wf,
                                         const InteractionEvidence& evidence,
                                         const InteractionGraphOptions& options) {
    InteractionGraph graph;
    graph.atom_count = wf.atoms.size();
    if (wf.atoms.empty()) return graph;

    InteractionEvidence resolved_evidence = evidence;
    if (!resolved_evidence.has_atomic_charges(wf.atoms.size()) &&
        wf.atomic_partial_charge_provenance != DataProvenance::Unavailable &&
        wf.atomic_partial_charges.size() == wf.atoms.size()) {
        resolved_evidence.atomic_partial_charges = wf.atomic_partial_charges;
        resolved_evidence.atomic_charge_provenance =
            wf.atomic_partial_charge_provenance;
    }

    std::map<AtomPair, const BondOrderRecord*> records;
    for (const auto& record : wf.bond_orders) {
        if (record.atom_a >= wf.atoms.size() || record.atom_b >= wf.atoms.size() ||
            record.atom_a == record.atom_b ||
            record.provenance==DataProvenance::Unavailable ||
            !std::isfinite(record.mayer_order)) {
            continue;
        }
        const AtomPair pair = ordered_pair(record.atom_a, record.atom_b);
        const auto found = records.find(pair);
        if (found == records.end() ||
            std::abs(record.mayer_order) > std::abs(found->second->mayer_order)) {
            records[pair] = &record;
        }
    }

    EdgeIndicesByPair edge_by_pair;
    std::map<AtomPair,const InteractionBondEvidence*> nbo_pairs;
    for(const auto& row:resolved_evidence.bonds)
        if(row.atom_a<wf.atoms.size() && row.atom_b<wf.atoms.size() && row.atom_a!=row.atom_b &&
           std::isfinite(row.wiberg_index) && row.wiberg_index>=0)
            nbo_pairs[ordered_pair(static_cast<std::uint32_t>(row.atom_a),
                                   static_cast<std::uint32_t>(row.atom_b))]=&row;
    const auto support_value=[&](const AtomPair& pair) {
        if(const auto nbo=nbo_pairs.find(pair);nbo!=nbo_pairs.end())return nbo->second->wiberg_index;
        if(const auto old=records.find(pair);old!=records.end())return old->second->mayer_order;
        return 0.0;
    };
    std::vector<InteractionEdge> coordination_candidates;
    std::vector<InteractionEdge> deferred_covalent_candidates;
    std::vector<std::vector<std::uint32_t>> covalent_adjacency(wf.atoms.size());

    std::vector<ValidatedMulticentre> multicentre;
    std::set<AtomPair> multicentre_support_pairs;
    std::set<AtomPair> multicentre_member_pairs;
    for (const auto& assignment : wf.multicentre_assignments) {
        auto validated = validate_multicentre_assignment(wf, assignment, options);
        if (validated.assignment == nullptr) continue;
        multicentre_support_pairs.insert(validated.support_pairs.begin(),
                                         validated.support_pairs.end());
        multicentre_member_pairs.insert(validated.member_pairs.begin(),
                                        validated.member_pairs.end());
        multicentre.push_back(std::move(validated));
    }
    std::set<AtomPair> suppressed_multicentre_pairs;
    std::set_difference(multicentre_member_pairs.begin(), multicentre_member_pairs.end(),
                        multicentre_support_pairs.begin(), multicentre_support_pairs.end(),
                        std::inserter(suppressed_multicentre_pairs,
                                      suppressed_multicentre_pairs.end()));

    // Reuse the renderer's already regression-tested structural-neighbour
    // envelope, while assigning richer semantics here. Determine X-H ligand
    // membership from the complete structural skeleton before classifying any
    // M...H edge, so the result cannot depend on bond-record iteration order.
    const auto structural_bonds=analyse_bonds(wf,resolved_evidence.bonds);
    std::set<std::uint32_t> ligand_bound_hydrogens;
    for (const auto& bond:structural_bonds) {
        const auto a=static_cast<std::uint32_t>(bond.atom_a);
        const auto b=static_cast<std::uint32_t>(bond.atom_b);
        if (wf.atoms[a].atomic_number==1 &&
            !is_coordination_centre_element(wf.atoms[b].atomic_number)) {
            ligand_bound_hydrogens.insert(a);
        }
        if (wf.atoms[b].atomic_number==1 &&
            !is_coordination_centre_element(wf.atoms[a].atomic_number)) {
            ligand_bound_hydrogens.insert(b);
        }
    }
    for (const auto& bond : structural_bonds) {
        const std::uint32_t a = static_cast<std::uint32_t>(bond.atom_a);
        const std::uint32_t b = static_cast<std::uint32_t>(bond.atom_b);
        if (suppressed_multicentre_pairs.count(ordered_pair(a, b)) != 0u) continue;
        const bool a_centre = is_coordination_centre_element(wf.atoms[a].atomic_number);
        const bool b_centre = is_coordination_centre_element(wf.atoms[b].atomic_number);
        const bool ligand_bound_metal_hydrogen=
            a_centre!=b_centre &&
            ((wf.atoms[a].atomic_number==1 &&
              ligand_bound_hydrogens.count(a)!=0u) ||
             (wf.atoms[b].atomic_number==1 &&
              ligand_bound_hydrogens.count(b)!=0u));
        if (ligand_bound_metal_hydrogen) continue;
        const bool metal_donor =
            (a_centre && !b_centre && is_coordination_donor_element(wf.atoms[b].atomic_number)) ||
            (b_centre && !a_centre && is_coordination_donor_element(wf.atoms[a].atomic_number));
        const bool electropositive_main_group_pair =
            (is_electropositive_main_group_metal(wf.atoms[a].atomic_number) &&
             is_coordination_donor_element(wf.atoms[b].atomic_number)) ||
            (is_electropositive_main_group_metal(wf.atoms[b].atomic_number) &&
             is_coordination_donor_element(wf.atoms[a].atomic_number));

        DataProvenance provenance=DataProvenance::Unavailable;
        (void)mayer_for_pair(records,a,b,provenance);
        const double support=support_value(ordered_pair(a,b));
        const bool charge_separated_weak_pair =
            electropositive_main_group_pair &&
            opposite_atomic_charge_pair(resolved_evidence, wf.atoms.size(), a, b, options) &&
            support < options.ionic_strong_bond_override_mayer;
        if (charge_separated_weak_pair) continue;

        InteractionEdge edge = make_edge(
            wf, records, a, b,
            metal_donor ? InteractionKind::CoordinationContact
                        : InteractionKind::CovalentConnectivity,
            InteractionStrength::StrongConnectivity,
            provenance == DataProvenance::Unavailable ? 0.58 : 0.86);
        if (metal_donor) {
            coordination_candidates.push_back(std::move(edge));
        } else if (resolved_evidence.has_atomic_charges(wf.atoms.size()) &&
                   provenance != DataProvenance::Unavailable &&
                   support < options.ionic_strong_bond_override_mayer) {
            // Borderline electronic links are evaluated only after a seed
            // graph of high-confidence molecular fragments exists. This
            // prevents a single H...F salt contact from erasing the very
            // fragment charges needed to recognise it.
            deferred_covalent_candidates.push_back(std::move(edge));
        } else {
            add_or_replace_strong_edge(graph, edge_by_pair, std::move(edge));
            covalent_adjacency[a].push_back(b);
            covalent_adjacency[b].push_back(a);
        }
    }

    // Supplement electronically supported low-order coordination edges that
    // deliberately fall below analyse_bonds()' ordinary renderer threshold.
    std::set<AtomPair> electronic_pairs;
    for(const auto& [pair,_]:records)electronic_pairs.insert(pair);
    for(const auto& [pair,_]:nbo_pairs)electronic_pairs.insert(pair);
    for (const auto& pair : electronic_pairs) {
        const double support=support_value(pair);
        if (has_pair(edge_by_pair, pair.first, pair.second) ||
            suppressed_multicentre_pairs.count(pair) != 0u ||
            support < options.coordination_mayer_floor) {
            continue;
        }
        const bool a_centre = is_coordination_centre_element(
            wf.atoms[pair.first].atomic_number);
        const bool b_centre = is_coordination_centre_element(
            wf.atoms[pair.second].atomic_number);
        if (a_centre == b_centre) continue;
        const std::uint32_t centre = a_centre ? pair.first : pair.second;
        const std::uint32_t donor = a_centre ? pair.second : pair.first;
        const double distance = distance_bohr(wf.atoms[centre], wf.atoms[donor]);
        if (distance > options.coordination_distance_factor *
                           radius_sum_bohr(wf.atoms[centre], wf.atoms[donor])) {
            continue;
        }
        if (wf.atoms[donor].atomic_number==1) {
            // A low-order but electronically supported direct M-H link is a
            // structural bond, not a dashed coordination contact.  Hydrogen
            // is deliberately excluded from the generic donor list, so this
            // dedicated path preserves weak terminal hydrides.  A hydrogen
            // already attached to another ordinary covalent neighbour is a
            // ligand C-H/O-H/H-H atom, however; its small through-space
            // metal population is agostic/dihydrogen interaction evidence,
            // not a second structural bond.
            const bool already_ligand_bound=
                ligand_bound_hydrogens.count(donor)!=0u || std::any_of(
                covalent_adjacency[donor].begin(),
                covalent_adjacency[donor].end(),
                [&](const std::uint32_t neighbour) {
                    return neighbour!=centre;
                });
            if (already_ligand_bound) continue;
            auto edge=make_edge(
                wf,records,centre,donor,
                InteractionKind::CovalentConnectivity,
                InteractionStrength::StrongConnectivity,
                std::clamp(0.55+0.8*support,0.55,0.82));
            add_or_replace_strong_edge(graph,edge_by_pair,std::move(edge));
            covalent_adjacency[centre].push_back(donor);
            covalent_adjacency[donor].push_back(centre);
            continue;
        }
        if (!is_coordination_donor_element(wf.atoms[donor].atomic_number)) continue;
        if (is_electropositive_main_group_metal(
                wf.atoms[centre].atomic_number) &&
            opposite_atomic_charge_pair(resolved_evidence, wf.atoms.size(), centre, donor, options) &&
            support < options.ionic_strong_bond_override_mayer) {
            continue;
        }
        coordination_candidates.push_back(make_edge(
            wf, records, centre, donor, InteractionKind::CoordinationContact,
            InteractionStrength::StrongConnectivity,
            std::clamp(0.55 + 0.8 * support, 0.55, 0.82)));
    }

    // A more distant atom hidden behind its ordinary covalent neighbour is not
    // a first-shell donor (e.g. O behind C in a metal carbonyl, or H behind O
    // in an aqua ligand).
    for (auto& edge : coordination_candidates) {
        const bool a_centre = is_coordination_centre_element(
            wf.atoms[edge.atom_a].atomic_number);
        const std::uint32_t centre = a_centre ? edge.atom_a : edge.atom_b;
        const std::uint32_t donor = a_centre ? edge.atom_b : edge.atom_a;
        if (donor_is_shielded(wf, centre, donor, covalent_adjacency,
                              options.coordination_shield_margin_angstrom)) {
            continue;
        }
        add_or_replace_strong_edge(graph, edge_by_pair, std::move(edge));
    }

    // Trusted multicentre assignments are retained as hyperedges. Spatially
    // adjacent member pairs form support connectivity, but distant terminal
    // atoms (F...F in linear HF2-, for example) do not become ordinary bonds.
    for (const auto& validated : multicentre) {
        const auto& assignment = *validated.assignment;
        MulticentreSupportGroup group;
        group.kind = assignment.kind;
        group.atoms = assignment.atoms;
        group.orbitals = assignment.orbitals;
        group.electron_count = assignment.electron_count;
        group.source_subspace_id = assignment.source_subspace_id;
        group.source_subspace_electron_count =
            assignment.source_subspace_electron_count;
        group.source_subspace_fraction = assignment.source_subspace_fraction;
        group.confidence = assignment.confidence;
        group.provenance = assignment.provenance;
        graph.multicentre_groups.push_back(std::move(group));

        for (const AtomPair& pair : validated.support_pairs) {
            add_or_replace_strong_edge(graph, edge_by_pair, make_edge(
                wf, records, pair.first, pair.second,
                InteractionKind::MulticentreSupport,
                InteractionStrength::StrongConnectivity,
                std::max(0.70, assignment.confidence)));
        }
    }

    graph.polyhedral_cages = analyse_polyhedral_cages(wf);
    for (const auto& cage : graph.polyhedral_cages) {
        for (const auto& pair : cage.support_edges) {
            add_or_replace_strong_edge(graph, edge_by_pair, make_edge(
                wf, records, pair[0], pair[1],
                InteractionKind::PolyhedralCageSupport,
                InteractionStrength::StrongConnectivity,
                std::max(0.70, cage.confidence)));
        }
    }

    const auto seed_fragments = derive_fragments(
        wf.atoms.size(), graph.edges, resolved_evidence);
    for (auto& edge : deferred_covalent_candidates) {
        const std::uint32_t fa = seed_fragments.atom_to_fragment[edge.atom_a];
        const std::uint32_t fb = seed_fragments.atom_to_fragment[edge.atom_b];
        bool ionic_bridge = false;
        if (fa != fb && fa < seed_fragments.fragments.size() &&
            fb < seed_fragments.fragments.size()) {
            const auto& left = seed_fragments.fragments[fa];
            const auto& right = seed_fragments.fragments[fb];
            const bool opposite_fragments =
                left.charge_provenance != DataProvenance::Unavailable &&
                right.charge_provenance != DataProvenance::Unavailable &&
                left.evidenced_charge * right.evidenced_charge < 0.0 &&
                std::abs(left.evidenced_charge) >= options.ionic_fragment_charge_floor &&
                std::abs(right.evidenced_charge) >= options.ionic_fragment_charge_floor;
            double strongest_a = 0.0;
            double strongest_b = 0.0;
            for (const auto& strong : graph.edges) {
                if (strong.strength != InteractionStrength::StrongConnectivity) continue;
                if (strong.atom_a == edge.atom_a || strong.atom_b == edge.atom_a) {
                    strongest_a = std::max(strongest_a, std::abs(support_value(ordered_pair(strong.atom_a,strong.atom_b))));
                }
                if (strong.atom_a == edge.atom_b || strong.atom_b == edge.atom_b) {
                    strongest_b = std::max(strongest_b, std::abs(support_value(ordered_pair(strong.atom_a,strong.atom_b))));
                }
            }
            const bool relatively_weak_at_both_ends = strongest_a > 0.0 &&
                strongest_b > 0.0 &&
                std::abs(support_value(ordered_pair(edge.atom_a,edge.atom_b))) < 0.45 * strongest_a &&
                std::abs(support_value(ordered_pair(edge.atom_a,edge.atom_b))) < 0.45 * strongest_b;
            ionic_bridge = opposite_fragments &&
                (edge.covalent_radius_ratio > 1.30 || relatively_weak_at_both_ends);
        }
        if (ionic_bridge) continue;
        covalent_adjacency[edge.atom_a].push_back(edge.atom_b);
        covalent_adjacency[edge.atom_b].push_back(edge.atom_a);
        add_or_replace_strong_edge(graph, edge_by_pair, std::move(edge));
    }

    graph.fragment_analysis = derive_fragments(
        wf.atoms.size(), graph.edges, resolved_evidence);

    std::vector<std::vector<std::uint32_t>> strong_adjacency(wf.atoms.size());
    for (const auto& edge : graph.edges) {
        if (edge.strength != InteractionStrength::StrongConnectivity) continue;
        strong_adjacency[edge.atom_a].push_back(edge.atom_b);
        strong_adjacency[edge.atom_b].push_back(edge.atom_a);
    }

    std::vector<WeakCandidate> weak_candidates;
    std::map<AtomPair, double> minimum_fragment_pair_gap;
    for (std::uint32_t a = 0; a < wf.atoms.size(); ++a) {
        for (std::uint32_t b = a + 1; b < wf.atoms.size(); ++b) {
            if (has_pair(edge_by_pair, a, b)) continue;
            const double distance = distance_bohr(wf.atoms[a], wf.atoms[b]);
            const double radii = radius_sum_bohr(wf.atoms[a], wf.atoms[b]);
            const double distance_angstrom = distance / kAngstromToBohr;
            if (distance_angstrom > options.weak_contact_max_angstrom ||
                distance > options.weak_contact_distance_factor * radii) {
                continue;
            }

            const std::uint32_t fa = graph.fragment_analysis.atom_to_fragment[a];
            const std::uint32_t fb = graph.fragment_analysis.atom_to_fragment[b];
            const bool hbond = hydrogen_bond_geometry(
                wf, a, b, strong_adjacency, options);
            if (fa == fb && !hbond) continue;

            WeakCandidate candidate;
            candidate.atom_a = a;
            candidate.atom_b = b;
            candidate.fragment_a = std::min(fa, fb);
            candidate.fragment_b = std::max(fa, fb);
            candidate.surface_gap_angstrom =
                distance_angstrom - radii / kAngstromToBohr;
            candidate.hydrogen_bond = hbond;
            weak_candidates.push_back(candidate);
            if (fa != fb) {
                const AtomPair fragments{candidate.fragment_a, candidate.fragment_b};
                const auto found = minimum_fragment_pair_gap.find(fragments);
                if (found == minimum_fragment_pair_gap.end()) {
                    minimum_fragment_pair_gap[fragments] = candidate.surface_gap_angstrom;
                } else {
                    found->second = std::min(found->second, candidate.surface_gap_angstrom);
                }
            }
        }
    }

    for (const auto& candidate : weak_candidates) {
        const bool interfragment = candidate.fragment_a != candidate.fragment_b;
        if (interfragment && !candidate.hydrogen_bond) {
            const AtomPair fragments{candidate.fragment_a, candidate.fragment_b};
            const double minimum_gap = minimum_fragment_pair_gap[fragments];
            if (candidate.surface_gap_angstrom >
                minimum_gap + options.weak_contact_shell_angstrom) {
                continue;
            }
        }

        InteractionKind kind = InteractionKind::NoncovalentContact;
        double confidence = 0.55;
        if (candidate.hydrogen_bond) {
            kind = InteractionKind::HydrogenBond;
            confidence = 0.78;
        } else if (interfragment &&
                   resolved_evidence.has_atomic_charges(wf.atoms.size())) {
            const auto& fa = graph.fragment_analysis.fragments[candidate.fragment_a];
            const auto& fb = graph.fragment_analysis.fragments[candidate.fragment_b];
            if (fa.evidenced_charge * fb.evidenced_charge < 0.0 &&
                std::abs(fa.evidenced_charge) >= options.ionic_fragment_charge_floor &&
                std::abs(fb.evidenced_charge) >= options.ionic_fragment_charge_floor) {
                kind = InteractionKind::IonicContact;
                confidence = 0.82;
            }
        }

        DataProvenance provenance = DataProvenance::Unavailable;
        const double mayer = mayer_for_pair(records, candidate.atom_a,
                                            candidate.atom_b, provenance);
        if (kind == InteractionKind::NoncovalentContact &&
            std::abs(mayer) >= options.weak_record_abs_mayer_floor) {
            kind = InteractionKind::AmbiguousContact;
            confidence = 0.45;
        }
        const std::uint32_t edge_index = static_cast<std::uint32_t>(graph.edges.size());
        graph.edges.push_back(make_edge(
            wf, records, candidate.atom_a, candidate.atom_b, kind,
            InteractionStrength::WeakContact, confidence));
        edge_by_pair[ordered_pair(candidate.atom_a, candidate.atom_b)]
            .push_back(edge_index);
    }

    for(auto& edge:graph.edges) {
        const auto pair=ordered_pair(edge.atom_a,edge.atom_b);
        if(const auto nbo=nbo_pairs.find(pair);nbo!=nbo_pairs.end()) {
            edge.wiberg_index=nbo->second->wiberg_index;
            edge.wiberg_source_path=nbo->second->source_path;
            edge.connectivity_method="Wiberg (NAO) and structural-neighbour geometry";
            edge.connectivity_source_path=nbo->second->source_path;
        } else {
            edge.connectivity_method=edge.electronic_provenance!=DataProvenance::Unavailable?
                "Mayer and structural-neighbour geometry":"covalent-radius geometry";
        }
    }
    // Re-derive once weak edges exist so interfragment edge indices are stable
    // and fragment charges remain attached to the final graph.
    graph.fragment_analysis = derive_fragments(
        wf.atoms.size(), graph.edges, resolved_evidence);
    return graph;
}

} // namespace cov
