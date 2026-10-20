#include "cov/molecule_style.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <queue>
#include <utility>
#include <vector>

namespace cov {
namespace {

constexpr double kMayerRenderFloor = 0.05;
// Mayer indices also contain non-local through-bond and pi-system coupling.
// Keep the electronic value as bond evidence, but do not turn those remote
// couplings into structural chords. The base 1.50 separates the next-neighbour
// shell in ordinary carbon rings. Strong direct electronic evidence can relax
// the envelope smoothly (needed by weakly bound one-electron bonds such as
// equilibrium H2+) without granting the same range to small remote pi Mayer
// terms.
constexpr double kElectronicAdjacencyFactor = 1.50;

// A positive pair index is not sufficient to make two atoms structural
// neighbours: delocalised density commonly gives a small Mayer value to the
// diagonal of a four- or five-membered ring.  Such a diagonal lies in the
// two-hop "lune" of two much shorter, already electronically supported local
// edges.  The constants below deliberately describe that topology rather than
// a particular molecule.  Near-equilateral three-membered rings do not meet
// the distance ratio, while a genuinely strong direct edge can override the
// weak/non-local test.
constexpr double kTwoHopLuneDistanceRatio = 0.82;
constexpr double kLocalLegRadiusFactor = 1.25;
constexpr double kLocalLegMayerFloor = 0.20;
constexpr double kDirectMayerOverride = 0.30;
constexpr double kRelativeDirectMayerOverride = 0.55;
// A weak edge can still close a genuine, distorted three-membered ring.  Such
// an edge has exactly one local two-hop witness, remains geometrically compact,
// and carries electronic support comparable with both legs.  By contrast a
// square diagonal has two witnesses, ordinary ring chords are too open, and a
// linear XY2 terminal--terminal contact is nearly the sum of its two legs.
constexpr double kCompactTriangleLongestEdgeRatio = 1.55;
constexpr double kCompactTriangleMaximumAngleDegrees = 105.0;
constexpr double kCompactTriangleRelativeMayerFloor = 0.75;

double electronic_adjacency_factor(const double mayer_order) noexcept {
    const double magnitude=std::abs(mayer_order);
    const double strong_fraction=std::clamp(
        (magnitude-0.20)/0.30,0.0,1.0);
    return kElectronicAdjacencyFactor+0.35*strong_fraction;
}

void prune_two_hop_electronic_chords(const Wavefunction& wavefunction,
                                     std::vector<BondVisual>& bonds) {
    const std::size_t atom_count = wavefunction.atoms.size();
    if (atom_count < 3u || bonds.size() < 3u) return;

    using AtomPair=std::pair<std::size_t,std::size_t>;
    const auto canonical_pair=[](const std::size_t a,const std::size_t b) {
        return a<b?AtomPair{a,b}:AtomPair{b,a};
    };
    std::map<AtomPair,std::size_t> pair_index;
    for (std::size_t index = 0u; index < bonds.size(); ++index) {
        const auto& bond = bonds[index];
        if (bond.atom_a >= atom_count || bond.atom_b >= atom_count) continue;
        pair_index[canonical_pair(bond.atom_a,bond.atom_b)]=index;
    }
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> adjacency(
        atom_count);
    for (const auto& [pair,index]:pair_index) {
        adjacency[pair.first].push_back({pair.second,index});
        adjacency[pair.second].push_back({pair.first,index});
    }

    const auto locally_credible = [&](const BondVisual& bond) {
        const Atom& a = wavefunction.atoms[bond.atom_a];
        const Atom& b = wavefunction.atoms[bond.atom_b];
        const double radius_sum_bohr =
            (covalent_radius_angstrom(a.atomic_number) +
             covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
        const double radius_ratio =
            bond.distance_bohr / std::max(1.0e-12, radius_sum_bohr);
        return radius_ratio <= kLocalLegRadiusFactor ||
               bond.bond_order >= kLocalLegMayerFloor;
    };

    std::vector<bool> remove(bonds.size(), false);
    for (std::size_t edge_index = 0u; edge_index < bonds.size(); ++edge_index) {
        const auto& edge = bonds[edge_index];
        if (edge.atom_a>=adjacency.size() || edge.atom_b>=adjacency.size()) {
            continue;
        }
        std::vector<std::pair<std::size_t,std::size_t>> two_hop_witnesses;
        for (const auto& [middle,left_index]:adjacency[edge.atom_a]) {
            if (middle == edge.atom_a || middle == edge.atom_b) continue;
            const auto right_found=pair_index.find(
                canonical_pair(middle,edge.atom_b));
            if (right_found==pair_index.end()) continue;
            const std::size_t right_index=right_found->second;
            const auto& left = bonds[left_index];
            const auto& right = bonds[right_index];
            if (!locally_credible(left) || !locally_credible(right)) continue;
            if (std::max(left.distance_bohr, right.distance_bohr) >
                kTwoHopLuneDistanceRatio * edge.distance_bohr) {
                continue;
            }
            two_hop_witnesses.push_back({left_index,right_index});
        }

        if (two_hop_witnesses.empty()) continue;

        const auto& [left_index,right_index]=two_hop_witnesses.front();
        const auto& left=bonds[left_index];
        const auto& right=bonds[right_index];
        const double weaker_leg=std::min(left.bond_order,right.bond_order);
        const bool independently_strong=
            edge.bond_order>=kDirectMayerOverride &&
            edge.bond_order>=kRelativeDirectMayerOverride*weaker_leg;
        if (independently_strong) continue;

        bool compact_three_cycle=false;
        if (two_hop_witnesses.size()==1u &&
            edge.bond_order>=kCompactTriangleRelativeMayerFloor*weaker_leg) {
            const std::size_t middle=
                left.atom_a==edge.atom_a?left.atom_b:left.atom_a;
            const Atom& a=wavefunction.atoms[edge.atom_a];
            const Atom& b=wavefunction.atoms[edge.atom_b];
            const Atom& m=wavefunction.atoms[middle];
            const double am=left.distance_bohr;
            const double bm=right.distance_bohr;
            const double longest_leg=std::max(am,bm);
            const double dot=(a.x-m.x)*(b.x-m.x)+
                             (a.y-m.y)*(b.y-m.y)+
                             (a.z-m.z)*(b.z-m.z);
            const double cosine=std::clamp(
                dot/std::max(1.0e-12,am*bm),-1.0,1.0);
            const double maximum_angle_cosine=std::cos(
                kCompactTriangleMaximumAngleDegrees*
                3.14159265358979323846/180.0);
            compact_three_cycle=
                edge.distance_bohr<=
                    kCompactTriangleLongestEdgeRatio*longest_leg &&
                cosine>=maximum_angle_cosine;
        }
        remove[edge_index]=!compact_three_cycle;
    }

    std::size_t output = 0u;
    for (std::size_t index = 0u; index < bonds.size(); ++index) {
        if (remove[index]) continue;
        if (output != index) bonds[output] = std::move(bonds[index]);
        ++output;
    }
    bonds.resize(output);
}

bool likely_pi_element(const int z) noexcept {
    switch (z) {
        case 5:  // B
        case 6:  // C
        case 7:  // N
        case 8:  // O
        case 14: // Si
        case 15: // P
        case 16: // S
            return true;
        default:
            return false;
    }
}

double atom_distance_bohr(const Atom& a, const Atom& b) noexcept {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<std::size_t> alternate_path(const std::vector<std::vector<std::size_t>>& adjacency,
                                        const std::size_t source,
                                        const std::size_t target,
                                        const std::size_t blocked_a,
                                        const std::size_t blocked_b,
                                        const std::size_t max_edges) {
    const std::size_t n = adjacency.size();
    std::vector<int> parent(n, -1);
    std::vector<std::size_t> depth(n, max_edges + 1);
    std::queue<std::size_t> queue;
    queue.push(source);
    depth[source] = 0;

    while (!queue.empty()) {
        const std::size_t current = queue.front();
        queue.pop();
        if (depth[current] >= max_edges) continue;

        for (const std::size_t next : adjacency[current]) {
            if ((current == blocked_a && next == blocked_b) ||
                (current == blocked_b && next == blocked_a)) {
                continue;
            }
            if (depth[next] <= depth[current] + 1) continue;
            parent[next] = static_cast<int>(current);
            depth[next] = depth[current] + 1;
            if (next == target) {
                std::vector<std::size_t> path;
                for (std::size_t p = target;; p = static_cast<std::size_t>(parent[p])) {
                    path.push_back(p);
                    if (p == source) break;
                }
                std::reverse(path.begin(), path.end());
                return path;
            }
            queue.push(next);
        }
    }
    return {};
}

std::size_t bond_index_between(const std::vector<BondVisual>& bonds,
                               const std::size_t a,
                               const std::size_t b) {
    for (std::size_t i = 0; i < bonds.size(); ++i) {
        if ((bonds[i].atom_a == a && bonds[i].atom_b == b) ||
            (bonds[i].atom_a == b && bonds[i].atom_b == a)) {
            return i;
        }
    }
    return bonds.size();
}

void mark_conservative_ring_delocalisation(const Wavefunction& wavefunction,
                                           std::vector<BondVisual>& bonds) {
    if (bonds.empty()) return;
    const std::size_t atom_count = wavefunction.atoms.size();

    std::vector<std::vector<std::size_t>> adjacency(atom_count);
    for (const auto& bond : bonds) {
        if (bond.atom_a >= atom_count || bond.atom_b >= atom_count) continue;
        adjacency[bond.atom_a].push_back(bond.atom_b);
        adjacency[bond.atom_b].push_back(bond.atom_a);
    }

    // Delocalisation remains a conservative rendering annotation. Electronic
    // pairwise bond order can improve connectivity, but it is not itself enough
    // to claim aromaticity. Compact 5-7 member pi-capable rings must also have
    // unusually short and near-equal edge lengths.
    for (std::size_t bi = 0; bi < bonds.size(); ++bi) {
        const BondVisual seed = bonds[bi];
        const auto path = alternate_path(adjacency, seed.atom_a, seed.atom_b,
                                         seed.atom_a, seed.atom_b, 6);
        if (path.size() < 5 || path.size() > 7) continue;

        std::vector<std::size_t> cycle_bonds;
        cycle_bonds.push_back(bi);
        bool valid = true;
        double min_ratio = 1.0e9;
        double max_ratio = 0.0;
        double ratio_sum = 0.0;
        std::size_t ratio_count = 0;

        for (const std::size_t atom_index : path) {
            if (atom_index >= atom_count ||
                !likely_pi_element(wavefunction.atoms[atom_index].atomic_number)) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        for (std::size_t pi = 0; pi + 1 < path.size(); ++pi) {
            const std::size_t bidx = bond_index_between(bonds, path[pi], path[pi + 1]);
            if (bidx >= bonds.size()) {
                valid = false;
                break;
            }
            cycle_bonds.push_back(bidx);
        }
        if (!valid) continue;

        for (const std::size_t bidx : cycle_bonds) {
            const auto& bond = bonds[bidx];
            const Atom& a = wavefunction.atoms[bond.atom_a];
            const Atom& b = wavefunction.atoms[bond.atom_b];
            const double single_length_bohr =
                (covalent_radius_angstrom(a.atomic_number) +
                 covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
            const double ratio = bond.distance_bohr / std::max(1.0e-12, single_length_bohr);
            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);
            ratio_sum += ratio;
            ++ratio_count;
        }

        if (!valid || ratio_count == 0) continue;
        const double mean_ratio = ratio_sum / static_cast<double>(ratio_count);
        const bool short_enough = mean_ratio >= 0.82 && mean_ratio <= 0.97;
        const bool near_equal = (max_ratio - min_ratio) <= 0.10;
        if (!short_enough || !near_equal) continue;

        for (const std::size_t bidx : cycle_bonds) bonds[bidx].delocalised = true;
    }
}

} // namespace

std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction,
                                    const std::vector<InteractionBondEvidence>& nbo_bonds) {
    std::vector<BondVisual> bonds;
    const std::size_t atom_count = wavefunction.atoms.size();

    using Pair=std::pair<std::size_t,std::size_t>;
    struct Electronic {double value=0;DataProvenance provenance=DataProvenance::Unavailable;
        std::string method,path;};
    std::map<Pair,Electronic> electronic;
    if(wavefunction.bond_order_provenance!=DataProvenance::Unavailable)
        for(const auto& row:wavefunction.bond_orders)
            if(row.atom_a<atom_count && row.atom_b<atom_count && row.atom_a!=row.atom_b &&
               row.provenance!=DataProvenance::Unavailable && std::isfinite(row.mayer_order))
                electronic[std::minmax(row.atom_a,row.atom_b)]={row.mayer_order,row.provenance,"Mayer",{}};
    for(const auto& row:nbo_bonds)
        if(row.atom_a<atom_count && row.atom_b<atom_count && row.atom_a!=row.atom_b &&
           std::isfinite(row.wiberg_index) && row.wiberg_index>=0)
            electronic[std::minmax(row.atom_a,row.atom_b)]={row.wiberg_index,DataProvenance::Derived,
                "Wiberg (NAO)",row.source_path};
    if (!electronic.empty() || wavefunction.bond_order_provenance!=DataProvenance::Unavailable) {
        // Electronic evidence is authoritative when available, while the
        // distance envelope determines whether that pair is a drawable
        // structural neighbour. This prevents delocalised Mayer coupling from
        // becoming chords across aromatic rings.
        for (const auto& [pair,record] : electronic) {
            const auto [i,j]=pair;
            const double electronic_strength=std::abs(record.value);
            if (i >= atom_count || j >= atom_count || i == j ||
                electronic_strength < kMayerRenderFloor) {
                continue;
            }
            const Atom& a = wavefunction.atoms[i];
            const Atom& b = wavefunction.atoms[j];
            const double distance_bohr = atom_distance_bohr(a, b);
            const double radius_sum_bohr =
                (covalent_radius_angstrom(a.atomic_number) +
                 covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
            // A negative pair contribution can rescue an obvious first-shell
            // skeleton edge, but it is not positive bonding evidence and must
            // never receive the relaxed electronic distance envelope.
            const double adjacency_factor=record.value<0.0
                ?kLocalLegRadiusFactor
                :electronic_adjacency_factor(electronic_strength);
            const double sanity_bohr =
                adjacency_factor*radius_sum_bohr;
            if (distance_bohr > sanity_bohr) continue;

            BondVisual bond;
            bond.atom_a = i;
            bond.atom_b = j;
            bond.distance_bohr = distance_bohr;
            // A signed Mayer contribution can be negative in diffuse or
            // heavily delocalised bases even for a short first-neighbour
            // contact. Connectivity uses its magnitude only inside the strict
            // geometry gate above; positive evidence alone gets the wider
            // adaptive envelope.
            bond.bond_order = electronic_strength;
            bond.provenance = record.provenance;
            bond.connectivity_method=record.method;
            bond.connectivity_source_path=record.path;
            bonds.push_back(bond);
        }
    }
    {
        // Capability is per pair: a partial electronic result cannot turn all
        // omitted pairs into exact zeros. Only genuinely missing pair values
        // use conservative geometry; an explicit zero stays authoritative.
        for (std::size_t i = 0; i < atom_count; ++i) {
            for (std::size_t j = i + 1; j < atom_count; ++j) {
                // Per-pair capability fallback: an actual WBI=0 vetoes the
                // geometry fallback; a missing NBO pair does not invent zero.
                if(electronic.contains({i,j}))continue;
                const Atom& a = wavefunction.atoms[i];
                const Atom& b = wavefunction.atoms[j];
                const double distance_bohr = atom_distance_bohr(a, b);
                const double cutoff_bohr =
                    1.25 * (covalent_radius_angstrom(a.atomic_number) +
                            covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
                if (distance_bohr <= cutoff_bohr) {
                    BondVisual bond;
                    bond.atom_a = i;
                    bond.atom_b = j;
                    bond.distance_bohr = distance_bohr;
                    bond.provenance = DataProvenance::Unavailable;
                    bond.connectivity_method="covalent-radius geometry";
                    bonds.push_back(bond);
                }
            }
        }
    }

    if (!electronic.empty() || wavefunction.bond_order_provenance != DataProvenance::Unavailable) {
        prune_two_hop_electronic_chords(wavefunction, bonds);
    }
    mark_conservative_ring_delocalisation(wavefunction, bonds);
    return bonds;
}

} // namespace cov
