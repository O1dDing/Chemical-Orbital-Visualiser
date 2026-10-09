#include "cov/orbital_chemistry.hpp"
#include "cov/molecule_style.hpp"
#include "cov/pi_topology.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace cov {
namespace {

enum class ReferenceClass : std::uint8_t {
    DeepCore = 0,
    Semicore,
    ChemicalValence,
};

struct ReferenceShell {
    int n = 0;
    int l = 0;
    int electrons = 0;
    ReferenceClass space = ReferenceClass::DeepCore;
    bool removed_by_ecp = false;
};

struct ReferenceColumn {
    std::vector<double> coefficients;
    std::uint32_t atom = 0;
    int n = 0;
    int l = 0;
    int component = 0;
    ReferenceClass space = ReferenceClass::DeepCore;
    std::string label;
};

struct AufbauEntry {
    int n;
    int l;
    int capacity;
};

constexpr std::array<AufbauEntry, 19> kAufbau{{
    {1,0,2}, {2,0,2}, {2,1,6}, {3,0,2}, {3,1,6},
    {4,0,2}, {3,2,10}, {4,1,6}, {5,0,2}, {4,2,10},
    {5,1,6}, {6,0,2}, {4,3,14}, {5,2,10}, {6,1,6},
    {7,0,2}, {5,3,14}, {6,2,10}, {7,1,6},
}};

int period_for_z(const int z) {
    if (z <= 0) return 0;
    if (z <= 2) return 1;
    if (z <= 10) return 2;
    if (z <= 18) return 3;
    if (z <= 36) return 4;
    if (z <= 54) return 5;
    if (z <= 86) return 6;
    return 7;
}

bool transition_metal(const int z) {
    return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) ||
           (z >= 72 && z <= 80) || (z >= 104 && z <= 112);
}

bool f_block(const int z) {
    return (z >= 57 && z <= 71) || (z >= 89 && z <= 103);
}

using ShellKey = std::pair<int,int>;

std::map<ShellKey,int> neutral_configuration(const int z) {
    std::map<ShellKey,int> result;
    int remaining = std::max(0, z);
    for (const auto& shell : kAufbau) {
        if (remaining <= 0) break;
        const int fill = std::min(shell.capacity, remaining);
        result[{shell.n,shell.l}] = fill;
        remaining -= fill;
    }
    return result;
}

std::vector<ReferenceShell> atomic_reference_shells(const Atom& atom) {
    const int z = atom.atomic_number;
    const int period = period_for_z(z);
    auto occupation = neutral_configuration(z);

    int explicit_charge = z;
    if (std::isfinite(atom.nuclear_charge) && atom.nuclear_charge > 0.0) {
        const int rounded = static_cast<int>(std::llround(atom.nuclear_charge));
        if (rounded > 0 && rounded <= z) explicit_charge = rounded;
    }

    int remove = std::max(0, z - explicit_charge);
    std::vector<ShellKey> occupied;
    for (const auto& [key, electrons] : occupation) {
        if (electrons > 0) occupied.push_back(key);
    }
    std::sort(occupied.begin(), occupied.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });

    std::set<ShellKey> removed;
    for (const auto& key : occupied) {
        if (remove <= 0) break;
        int& electrons = occupation[key];
        const int take = std::min(remove, electrons);
        electrons -= take;
        remove -= take;
        if (electrons == 0) removed.insert(key);
    }

    std::set<ShellKey> valence;
    std::set<ShellKey> semicore;
    if (period == 1) {
        valence.insert({1,0});
    } else if (period > 1) {
        valence.insert({period,0});
        valence.insert({period,1});
        if (transition_metal(z)) {
            valence.insert({period-1,2});
            semicore.insert({period-1,0});
            semicore.insert({period-1,1});
        } else if (f_block(z)) {
            valence.insert({period-2,3});
            valence.insert({period-1,2});
            semicore.insert({period-1,0});
            semicore.insert({period-1,1});
        } else if (period >= 4) {
            semicore.insert({period-1,2});
        }
    }

    std::set<ShellKey> keys;
    for (const auto& [key, electrons] : occupation) {
        if (electrons > 0 || removed.count(key) != 0u) keys.insert(key);
    }
    keys.insert(valence.begin(), valence.end());
    keys.insert(semicore.begin(), semicore.end());

    std::vector<ReferenceShell> result;
    result.reserve(keys.size());
    for (const auto& key : keys) {
        ReferenceShell shell;
        shell.n = key.first;
        shell.l = key.second;
        shell.electrons = occupation[key];
        shell.removed_by_ecp = removed.count(key) != 0u;
        if (valence.count(key) != 0u) shell.space = ReferenceClass::ChemicalValence;
        else if (semicore.count(key) != 0u) shell.space = ReferenceClass::Semicore;
        else shell.space = ReferenceClass::DeepCore;
        result.push_back(shell);
    }
    return result;
}

double shell_effective_exponent(const Wavefunction& wf, const Shell& shell) {
    const std::size_t begin = shell.primitive_offset;
    const std::size_t count = shell.primitive_count;
    if (count == 0u || begin + count > wf.primitives.size()) return 0.0;
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t p = 0; p < count; ++p) {
        const auto& primitive = wf.primitives[begin+p];
        const double weight = std::max(1.0e-12,
            std::abs(static_cast<double>(primitive.coefficient)));
        numerator += weight * std::max(1.0e-20,
            static_cast<double>(primitive.exponent));
        denominator += weight;
    }
    return denominator > 0.0 ? numerator / denominator : 0.0;
}

const char* angular_letter(const int l) {
    static constexpr const char* letters[] = {"s","p","d","f","g","h"};
    return l >= 0 && l < 6 ? letters[l] : "?";
}

std::string reference_label(const Wavefunction& wf,
                            const std::uint32_t atom,
                            const int n,
                            const int l) {
    std::ostringstream out;
    if (atom < wf.atoms.size()) out << wf.atoms[atom].symbol;
    else out << 'X';
    out << (atom + 1u) << ' ' << n << angular_letter(l);
    return out.str();
}

double s_dot(const std::vector<double>& a,
             const std::vector<double>& b,
             const std::vector<double>& overlap,
             const std::size_t n) {
    double result = 0.0;
    for (std::size_t mu = 0; mu < n; ++mu) {
        if (a[mu] == 0.0) continue;
        double sb = 0.0;
        for (std::size_t nu = 0; nu < n; ++nu) {
            sb += overlap[mu*n+nu] * b[nu];
        }
        result += a[mu] * sb;
    }
    return result;
}

double mo_s_dot(const std::vector<double>& a,
                const MolecularOrbital& mo,
                const std::vector<double>& overlap,
                const std::size_t n) {
    double result = 0.0;
    for (std::size_t mu = 0; mu < n; ++mu) {
        if (a[mu] == 0.0) continue;
        double sc = 0.0;
        for (std::size_t nu = 0; nu < n; ++nu) {
            sc += overlap[mu*n+nu] *
                  static_cast<double>(mo.coefficients[nu]);
        }
        result += a[mu] * sc;
    }
    return result;
}

std::vector<std::size_t> largest_gap_boundaries(
    const std::vector<std::pair<std::size_t,double>>& shells,
    const std::size_t group_count) {
    std::vector<std::pair<double,std::size_t>> gaps;
    for (std::size_t i = 0; i + 1u < shells.size(); ++i) {
        const double a = std::log(std::max(1.0e-30, shells[i].second));
        const double b = std::log(std::max(1.0e-30, shells[i+1u].second));
        gaps.push_back({a-b,i+1u});
    }
    std::sort(gaps.begin(), gaps.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });
    std::vector<std::size_t> boundaries;
    const std::size_t wanted = group_count > 0u ? group_count-1u : 0u;
    for (std::size_t i=0;i<std::min(wanted,gaps.size());++i) {
        boundaries.push_back(gaps[i].second);
    }
    std::sort(boundaries.begin(),boundaries.end());
    return boundaries;
}

std::vector<double> radial_weights(
    const Wavefunction& wf,
    const std::vector<std::size_t>& shell_indices) {
    std::vector<double> result(shell_indices.size(),0.0);
    if (shell_indices.empty()) return result;
    double max_log = -std::numeric_limits<double>::infinity();
    for (const auto s : shell_indices) {
        max_log = std::max(max_log,
            std::log(std::max(1.0e-30,
                shell_effective_exponent(wf,wf.shells[s]))));
    }
    double norm = 0.0;
    for (std::size_t i=0;i<shell_indices.size();++i) {
        const double loga = std::log(std::max(1.0e-30,
            shell_effective_exponent(wf,wf.shells[shell_indices[i]])));
        result[i] = std::exp(0.65 * (loga-max_log));
        norm += result[i]*result[i];
    }
    if (norm > 0.0) {
        const double inv = 1.0/std::sqrt(norm);
        for (double& v : result) v *= inv;
    }
    return result;
}

void append_standard_components(
    const Wavefunction& wf,
    const std::vector<std::size_t>& shell_indices,
    const std::vector<double>& radial,
    const ReferenceShell& ref,
    const std::uint32_t atom,
    std::vector<ReferenceColumn>& raw) {
    if (shell_indices.empty() || radial.size()!=shell_indices.size()) return;
    const Shell& first = wf.shells[shell_indices.front()];
    const int l = static_cast<int>(first.angular_momentum);
    for (const auto s : shell_indices) {
        if (wf.shells[s].angular_momentum != first.angular_momentum ||
            wf.shells[s].pure != first.pure ||
            shell_basis_count(wf.shells[s]) != shell_basis_count(first)) {
            return;
        }
    }

    const std::size_t nbf = wf.basis_count;
    const std::size_t local_count = shell_basis_count(first);

    auto make_direct = [&](const int standard_component,
                           const std::vector<std::pair<int,double>>& local) {
        ReferenceColumn column;
        column.coefficients.assign(nbf,0.0);
        column.atom=atom;
        column.n=ref.n;
        column.l=l;
        column.component=standard_component;
        column.space=ref.space;
        column.label=reference_label(wf,atom,ref.n,l);
        for (std::size_t r=0;r<shell_indices.size();++r) {
            const Shell& shell=wf.shells[shell_indices[r]];
            for (const auto& [component,coefficient] : local) {
                if (component < 0 ||
                    static_cast<std::size_t>(component)>=local_count) continue;
                column.coefficients[shell.basis_offset+
                    static_cast<std::size_t>(component)] += radial[r]*coefficient;
            }
        }
        raw.push_back(std::move(column));
    };

    if (l <= 1 || first.pure != 0u) {
        const std::size_t chemical_count =
            l >= 0 ? static_cast<std::size_t>(2*l+1) : 0u;
        const std::size_t count = std::min(local_count,chemical_count);
        for (std::size_t c=0;c<count;++c) {
            make_direct(static_cast<int>(c),{{static_cast<int>(c),1.0}});
        }
        return;
    }

    if (l == 2 && local_count == 6u) {
        make_direct(0,{{0,-1.0/std::sqrt(6.0)},
                       {1,-1.0/std::sqrt(6.0)},
                       {2, 2.0/std::sqrt(6.0)}});
        make_direct(1,{{4,1.0}});
        make_direct(2,{{5,1.0}});
        make_direct(3,{{0, 1.0/std::sqrt(2.0)},
                       {1,-1.0/std::sqrt(2.0)}});
        make_direct(4,{{3,1.0}});
    }
}

std::vector<ReferenceColumn> build_raw_reference(const Wavefunction& wf) {
    std::vector<ReferenceColumn> raw;
    for (std::uint32_t atom=0;atom<wf.atoms.size();++atom) {
        const auto refs_all=atomic_reference_shells(wf.atoms[atom]);
        for (int l=0;l<=4;++l) {
            std::vector<ReferenceShell> refs;
            for (const auto& ref:refs_all) {
                if (!ref.removed_by_ecp && ref.l==l) refs.push_back(ref);
            }
            if (refs.empty()) continue;
            std::sort(refs.begin(),refs.end(),[](const auto& a,const auto& b) {
                return a.n<b.n;
            });

            std::vector<std::pair<std::size_t,double>> working;
            for (std::size_t s=0;s<wf.shells.size();++s) {
                const auto& shell=wf.shells[s];
                if (shell.atom_index==atom &&
                    static_cast<int>(shell.angular_momentum)==l) {
                    working.push_back({s,shell_effective_exponent(wf,shell)});
                }
            }
            if (working.empty()) continue;
            std::sort(working.begin(),working.end(),[](const auto& a,const auto& b) {
                if (a.second!=b.second) return a.second>b.second;
                return a.first<b.first;
            });

            const std::size_t groups=std::min(refs.size(),working.size());
            if (groups==0u) continue;
            const auto boundaries=largest_gap_boundaries(working,groups);
            std::size_t begin=0;
            for (std::size_t group=0;group<groups;++group) {
                const std::size_t end =
                    group<boundaries.size()?boundaries[group]:working.size();
                std::vector<std::size_t> shell_indices;
                for (std::size_t i=begin;i<end;++i) {
                    shell_indices.push_back(working[i].first);
                }
                if (shell_indices.empty() && begin<working.size()) {
                    shell_indices.push_back(working[begin].first);
                }
                const auto radial=radial_weights(wf,shell_indices);
                const std::size_t ref_offset=refs.size()-groups;
                append_standard_components(
                    wf,shell_indices,radial,refs[ref_offset+group],atom,raw);
                begin=end;
            }
        }
    }
    return raw;
}

int reference_class_order(const ReferenceClass value) {
    switch (value) {
        case ReferenceClass::DeepCore: return 0;
        case ReferenceClass::Semicore: return 1;
        default: return 2;
    }
}

std::vector<ReferenceColumn> orthonormalise_reference(
    std::vector<ReferenceColumn> raw,
    const Wavefunction& wf) {
    const std::size_t n=wf.basis_count;
    std::stable_sort(raw.begin(),raw.end(),[](const auto& a,const auto& b) {
        const int ca=reference_class_order(a.space);
        const int cb=reference_class_order(b.space);
        if (ca!=cb) return ca<cb;
        if (a.atom!=b.atom) return a.atom<b.atom;
        if (a.n!=b.n) return a.n<b.n;
        if (a.l!=b.l) return a.l<b.l;
        return a.component<b.component;
    });

    std::vector<ReferenceColumn> q;
    for (auto& candidate:raw) {
        if (candidate.coefficients.size()!=n) continue;
        for (const auto& previous:q) {
            const double projection=s_dot(
                previous.coefficients,candidate.coefficients,wf.ao_overlap,n);
            for (std::size_t mu=0;mu<n;++mu) {
                candidate.coefficients[mu]-=
                    projection*previous.coefficients[mu];
            }
        }
        const double norm2=s_dot(
            candidate.coefficients,candidate.coefficients,wf.ao_overlap,n);
        if (!(norm2>1.0e-10) || !std::isfinite(norm2)) continue;
        const double inv=1.0/std::sqrt(norm2);
        for (double& v:candidate.coefficients) v*=inv;
        q.push_back(std::move(candidate));
    }
    return q;
}

std::vector<ReferenceColumn> core_orthogonalised_valence_p_reference(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw) {
    std::vector<ReferenceColumn> result=raw;
    for (auto& candidate:result) {
        if (candidate.space!=ReferenceClass::ChemicalValence ||
            candidate.l!=1 || candidate.coefficients.size()!=wf.basis_count) {
            continue;
        }
        std::vector<ReferenceColumn> lower_p;
        for (const auto& column:raw) {
            if (column.atom==candidate.atom && column.l==1 &&
                column.n<candidate.n &&
                column.coefficients.size()==wf.basis_count) {
                lower_p.push_back(column);
            }
        }
        const auto core=orthonormalise_reference(std::move(lower_p),wf);
        for (const auto& previous:core) {
            const double projection=s_dot(
                previous.coefficients,candidate.coefficients,
                wf.ao_overlap,wf.basis_count);
            for (std::size_t mu=0u;mu<wf.basis_count;++mu) {
                candidate.coefficients[mu]-=
                    projection*previous.coefficients[mu];
            }
        }
        const double norm2=s_dot(
            candidate.coefficients,candidate.coefficients,
            wf.ao_overlap,wf.basis_count);
        if (!(norm2>1.0e-10) || !std::isfinite(norm2)) {
            candidate.coefficients.clear();
            continue;
        }
        const double inverse=1.0/std::sqrt(norm2);
        for (double& value:candidate.coefficients) value*=inverse;
    }
    return result;
}

std::vector<std::uint32_t> basis_atom_map(const Wavefunction& wf) {
    const auto missing=std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> result(wf.basis_count,missing);
    for (const auto& shell:wf.shells) {
        const std::size_t count=shell_basis_count(shell);
        if (shell.basis_offset+count>wf.basis_count) return {};
        for (std::size_t c=0;c<count;++c) {
            result[shell.basis_offset+c]=shell.atom_index;
        }
    }
    if (std::find(result.begin(),result.end(),missing)!=result.end()) return {};
    return result;
}

std::vector<double> atomic_weights(
    const Wavefunction& wf,
    const MolecularOrbital& mo,
    const std::vector<std::uint32_t>& basis_atom) {
    const std::size_t n=wf.basis_count;
    std::vector<double> result(wf.atoms.size(),0.0);
    if (mo.coefficients.size()!=n ||
        wf.ao_overlap.size()!=n*n ||
        basis_atom.size()!=n) return result;
    std::vector<double> sc(n,0.0);
    for (std::size_t mu=0;mu<n;++mu) {
        for (std::size_t nu=0;nu<n;++nu) {
            sc[mu]+=wf.ao_overlap[mu*n+nu]*
                static_cast<double>(mo.coefficients[nu]);
        }
    }
    double total=0.0;
    for (std::size_t mu=0;mu<n;++mu) {
        result[basis_atom[mu]]+=
            static_cast<double>(mo.coefficients[mu])*sc[mu];
    }
    for (double& v:result) {
        v=std::abs(v);
        total+=v;
    }
    if (total>1.0e-12) for (double& v:result) v/=total;
    return result;
}

struct DegenerateGroup {
    Spin spin=Spin::Alpha;
    std::vector<std::size_t> orbitals;
    double score=0.0;
};

std::vector<DegenerateGroup> degeneracy_groups(
    const Wavefunction& wf,
    const double tolerance) {
    std::vector<DegenerateGroup> groups;
    std::size_t begin=0;
    while (begin<wf.orbitals.size()) {
        DegenerateGroup group;
        group.spin=wf.orbitals[begin].spin;
        group.orbitals.push_back(begin);
        std::size_t end=begin+1u;
        while (end<wf.orbitals.size()) {
            const auto& first=wf.orbitals[begin];
            const auto& current=wf.orbitals[end];
            if (current.spin!=first.spin) break;
            if (std::abs(current.energy_hartree-first.energy_hartree)>tolerance) break;
            if (!first.symmetry.empty() && !current.symmetry.empty() &&
                first.symmetry!=current.symmetry) break;
            group.orbitals.push_back(end);
            ++end;
        }
        groups.push_back(std::move(group));
        begin=end;
    }
    return groups;
}

void select_valence_manifold(
    Wavefunction& wf,
    const std::size_t reference_rank,
    const double tolerance) {
    if (reference_rank==0u || wf.orbitals.empty()) return;
    auto groups=degeneracy_groups(wf,tolerance);
    for (auto& group:groups) {
        for (const auto i:group.orbitals) {
            group.score+=wf.orbitals[i].chemistry.valence_weight;
        }
    }

    for (const Spin spin : {Spin::Alpha,Spin::Beta}) {
        std::vector<std::size_t> candidates;
        std::size_t orbital_count=0;
        for (std::size_t g=0;g<groups.size();++g) {
            if (groups[g].spin==spin) {
                candidates.push_back(g);
                orbital_count+=groups[g].orbitals.size();
            }
        }
        if (candidates.empty()) continue;
        const std::size_t target=std::min(reference_rank,orbital_count);
        const double neg=-std::numeric_limits<double>::infinity();
        std::vector<double> dp(orbital_count+1u,neg);
        std::vector<std::vector<std::size_t>> chosen(orbital_count+1u);
        dp[0]=0.0;
        for (const auto gi:candidates) {
            const std::size_t size=groups[gi].orbitals.size();
            const double value=groups[gi].score;
            for (std::size_t used=orbital_count+1u;used-- > size;) {
                if (!std::isfinite(dp[used-size])) continue;
                const double trial=dp[used-size]+value;
                if (trial>dp[used]+1.0e-12) {
                    dp[used]=trial;
                    chosen[used]=chosen[used-size];
                    chosen[used].push_back(gi);
                }
            }
        }

        std::size_t best=target;
        if (!std::isfinite(dp[best])) {
            std::size_t distance=orbital_count+1u;
            double best_score=neg;
            for (std::size_t count=0;count<=orbital_count;++count) {
                if (!std::isfinite(dp[count])) continue;
                const std::size_t d=count>target?count-target:target-count;
                if (d<distance || (d==distance && dp[count]>best_score)) {
                    distance=d;
                    best_score=dp[count];
                    best=count;
                }
            }
        }
        for (const auto gi:chosen[best]) {
            for (const auto i:groups[gi].orbitals) {
                wf.orbitals[i].chemistry.valence_manifold=true;
            }
        }
    }
}

double pair_mayer(const Wavefunction& wf,
                  const std::uint32_t a,
                  const std::uint32_t b) {
    for (const auto& record:wf.bond_orders) {
        if ((record.atom_a==a && record.atom_b==b) ||
            (record.atom_a==b && record.atom_b==a)) {
            return record.mayer_order;
        }
    }
    return 0.0;
}

double distance_bohr(const Atom& a,const Atom& b) {
    const double dx=a.x-b.x;
    const double dy=a.y-b.y;
    const double dz=a.z-b.z;
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}

std::array<double,3> unit_axis(const Atom& a,const Atom& b) {
    std::array<double,3> z{b.x-a.x,b.y-a.y,b.z-a.z};
    const double n=std::sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
    if (n>1.0e-14) for (double& v:z) v/=n;
    return z;
}

std::array<double,3> cross3(
    const std::array<double,3>& a,
    const std::array<double,3>& b) {
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

double dot3(const std::array<double,3>& a,
            const std::array<double,3>& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

struct RawChannels {
    double sigma=0.0;
    double pi=0.0;
    double delta=0.0;
    double phi=0.0;
    double undetermined=0.0;
};

using Mat3=std::array<double,9>;

Mat3 mat_add_scaled(Mat3 a,const Mat3& b,const double scale) {
    for (std::size_t i=0;i<9u;++i) a[i]+=scale*b[i];
    return a;
}

double frobenius(const Mat3& a,const Mat3& b) {
    double result=0.0;
    for (std::size_t i=0;i<9u;++i) result+=a[i]*b[i];
    return result;
}

Mat3 transform_tensor(const Mat3& t,
                      const std::array<double,3>& ex,
                      const std::array<double,3>& ey,
                      const std::array<double,3>& ez) {
    const std::array<std::array<double,3>,3> e{ex,ey,ez};
    Mat3 out{};
    for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
            double value=0.0;
            for (int a=0;a<3;++a) {
                for (int b=0;b<3;++b) {
                    value+=e[i][a]*t[3*a+b]*e[j][b];
                }
            }
            out[3*i+j]=value;
        }
    }
    return out;
}

void add_atom_channels(
    const std::uint32_t atom,
    const std::array<double,3>& axis,
    const std::vector<ReferenceColumn>& reference,
    const std::vector<double>& amplitudes,
    RawChannels& out) {
    std::map<std::pair<int,int>,std::vector<double>> by_shell;
    for (std::size_t r=0;r<reference.size();++r) {
        const auto& ref=reference[r];
        if (ref.atom!=atom || ref.space!=ReferenceClass::ChemicalValence) continue;
        auto& values=by_shell[{ref.n,ref.l}];
        const std::size_t needed=static_cast<std::size_t>(
            std::max(0,2*ref.l+1));
        if (values.size()<needed) values.resize(needed,0.0);
        if (ref.component>=0 &&
            static_cast<std::size_t>(ref.component)<values.size()) {
            values[static_cast<std::size_t>(ref.component)]=amplitudes[r];
        }
    }

    std::array<double,3> ez=axis;
    std::array<double,3> seed=
        std::abs(ez[2])<0.85?std::array<double,3>{0.0,0.0,1.0}
                            :std::array<double,3>{1.0,0.0,0.0};
    auto ex=cross3(seed,ez);
    double exn=std::sqrt(dot3(ex,ex));
    if (exn<1.0e-12) {
        out.undetermined+=1.0;
        return;
    }
    for (double& v:ex) v/=exn;
    auto ey=cross3(ez,ex);

    for (const auto& [key,values]:by_shell) {
        const int l=key.second;
        double norm=0.0;
        for (const double v:values) norm+=v*v;
        if (norm<=1.0e-14) continue;
        if (l==0) {
            out.sigma+=norm;
        } else if (l==1 && values.size()>=3u) {
            const std::array<double,3> p{values[0],values[1],values[2]};
            const double along=dot3(p,ez);
            const double sigma=along*along;
            out.sigma+=sigma;
            out.pi+=std::max(0.0,norm-sigma);
        } else if (l==2 && values.size()>=5u) {
            const double inv_sqrt6=1.0/std::sqrt(6.0);
            const double inv_sqrt2=1.0/std::sqrt(2.0);
            const Mat3 b0{-inv_sqrt6,0,0, 0,-inv_sqrt6,0, 0,0,2*inv_sqrt6};
            const Mat3 b1{0,0,inv_sqrt2, 0,0,0, inv_sqrt2,0,0};
            const Mat3 b2{0,0,0, 0,0,inv_sqrt2, 0,inv_sqrt2,0};
            const Mat3 b3{inv_sqrt2,0,0, 0,-inv_sqrt2,0, 0,0,0};
            const Mat3 b4{0,inv_sqrt2,0, inv_sqrt2,0,0, 0,0,0};
            Mat3 tensor{};
            tensor=mat_add_scaled(tensor,b0,values[0]);
            tensor=mat_add_scaled(tensor,b1,values[1]);
            tensor=mat_add_scaled(tensor,b2,values[2]);
            tensor=mat_add_scaled(tensor,b3,values[3]);
            tensor=mat_add_scaled(tensor,b4,values[4]);
            const Mat3 local=transform_tensor(tensor,ex,ey,ez);
            const double c0=frobenius(local,b0);
            const double c1=frobenius(local,b1);
            const double c2=frobenius(local,b2);
            const double c3=frobenius(local,b3);
            const double c4=frobenius(local,b4);
            out.sigma+=c0*c0;
            out.pi+=c1*c1+c2*c2;
            out.delta+=c3*c3+c4*c4;
        } else {
            if (std::abs(std::abs(ez[2])-1.0)<1.0e-6) {
                if (!values.empty()) out.sigma+=values[0]*values[0];
                if (values.size()>2u) out.pi+=values[1]*values[1]+values[2]*values[2];
                if (values.size()>4u) out.delta+=values[3]*values[3]+values[4]*values[4];
                if (values.size()>6u) out.phi+=values[5]*values[5]+values[6]*values[6];
            } else {
                out.undetermined+=norm;
            }
        }
    }
}

OrbitalChannelDistribution finish_channels(
    const RawChannels& raw,
    const OrbitalChemistryOptions& options) {
    OrbitalChannelDistribution result;
    const double total=raw.sigma+raw.pi+raw.delta+raw.phi+raw.undetermined;
    if (!(total>1.0e-14) || !std::isfinite(total)) {
        result.status=ChemistryStatus::Undetermined;
        return result;
    }
    result.sigma=raw.sigma/total;
    result.pi=raw.pi/total;
    result.delta=raw.delta/total;
    result.phi=raw.phi/total;
    result.undetermined=raw.undetermined/total;

    const std::array<std::pair<double,OrbitalAngularFamily>,4> values{{
        {result.sigma,OrbitalAngularFamily::Sigma},
        {result.pi,OrbitalAngularFamily::Pi},
        {result.delta,OrbitalAngularFamily::Delta},
        {result.phi,OrbitalAngularFamily::Phi},
    }};
    const auto best=std::max_element(values.begin(),values.end(),
        [](const auto& a,const auto& b){return a.first<b.first;});
    result.dominant=best->second;
    if (best->first>=options.determined_fraction &&
        result.undetermined<=0.10) {
        result.status=ChemistryStatus::Determined;
    } else if (result.undetermined<=
               options.maximum_undetermined_for_percentages) {
        result.status=ChemistryStatus::Percentages;
    } else {
        result.status=ChemistryStatus::Undetermined;
        result.dominant=OrbitalAngularFamily::Mixed;
    }
    return result;
}

OrbitalBondingDistribution pair_bonding(
    const double overlap_character,
    const bool applicable,
    const OrbitalChemistryOptions& options) {
    OrbitalBondingDistribution result;
    if (!applicable) {
        result.status=ChemistryStatus::NotApplicable;
        result.dominant=OrbitalBondingRole::NotApplicable;
        result.undetermined=0.0;
        return result;
    }
    if (overlap_character>=options.overlap_bonding_threshold) {
        result.bonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Bonding;
        result.status=ChemistryStatus::Determined;
    } else if (overlap_character<=-options.overlap_bonding_threshold) {
        result.antibonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Antibonding;
        result.status=ChemistryStatus::Determined;
    } else if (std::abs(overlap_character)<=
               options.nonbonding_overlap_threshold) {
        result.nonbonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Nonbonding;
        result.status=ChemistryStatus::Determined;
    } else {
        const double magnitude=std::abs(overlap_character);
        const double fraction=std::min(0.80,
            magnitude/options.overlap_bonding_threshold);
        if (overlap_character>0.0) result.bonding=fraction;
        else result.antibonding=fraction;
        result.nonbonding=1.0-fraction;
        result.undetermined=0.0;
        result.dominant=overlap_character>0.0
            ?OrbitalBondingRole::Bonding
            :OrbitalBondingRole::Antibonding;
        result.status=ChemistryStatus::Percentages;
    }
    return result;
}

OrbitalBondingDistribution aggregate_bonding(
    const std::vector<OrbitalPairInteraction>& interactions) {
    OrbitalBondingDistribution result;
    // The public default describes unavailable data (UND=1). An accumulator
    // starts at zero; otherwise each aggregate receives one fictitious UND
    // contribution that depends on the total interaction weight.
    result.undetermined=0.0;
    double total=0.0;
    for (const auto& interaction:interactions) {
        if (interaction.bonding.status==ChemistryStatus::NotApplicable) continue;
        const double weight=std::max(0.01,
            std::abs(interaction.overlap_character));
        result.bonding+=weight*interaction.bonding.bonding;
        result.antibonding+=weight*interaction.bonding.antibonding;
        result.nonbonding+=weight*interaction.bonding.nonbonding;
        result.undetermined+=weight*interaction.bonding.undetermined;
        total+=weight;
    }
    if (!(total>0.0)) {
        result.status=ChemistryStatus::Undetermined;
        result.undetermined=1.0;
        return result;
    }
    result.bonding/=total;
    result.antibonding/=total;
    result.nonbonding/=total;
    result.undetermined/=total;
    const std::array<std::pair<double,OrbitalBondingRole>,3> values{{
        {result.bonding,OrbitalBondingRole::Bonding},
        {result.antibonding,OrbitalBondingRole::Antibonding},
        {result.nonbonding,OrbitalBondingRole::Nonbonding},
    }};
    const auto best=std::max_element(values.begin(),values.end(),
        [](const auto& a,const auto& b){return a.first<b.first;});
    result.dominant=best->second;
    result.status=best->first>=0.82
        ?ChemistryStatus::Determined
        :ChemistryStatus::Percentages;
    return result;
}

OrbitalChannelDistribution aggregate_channels(
    const std::vector<OrbitalPairInteraction>& interactions,
    const OrbitalChemistryOptions& options) {
    RawChannels raw;
    for (const auto& interaction:interactions) {
        const double weight=std::max(0.01,
            std::abs(interaction.overlap_character));
        raw.sigma+=weight*interaction.channel.sigma;
        raw.pi+=weight*interaction.channel.pi;
        raw.delta+=weight*interaction.channel.delta;
        raw.phi+=weight*interaction.channel.phi;
        raw.undetermined+=weight*interaction.channel.undetermined;
    }
    return finish_channels(raw,options);
}

std::string family_symbol(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        default: return "mixed";
    }
}

bool near(const double value,const double target,const double tolerance) {
    return std::abs(value-target)<=tolerance;
}

bool pi_connectivity_evidence(const Wavefunction& wf,
                              const std::uint32_t a,
                              const std::uint32_t b) {
    if (a>=wf.atoms.size() || b>=wf.atoms.size() || a==b) return false;
    const double radii_bohr=
        (covalent_radius_angstrom(wf.atoms[a].atomic_number)+
         covalent_radius_angstrom(wf.atoms[b].atomic_number))*
        kAngstromToBohr;
    const double distance=distance_bohr(wf.atoms[a],wf.atoms[b]);

    if (wf.bond_order_provenance!=DataProvenance::Unavailable) {
        // Pi topology must follow the structural first-neighbour graph.  A
        // long or antibonding pair can carry a non-zero signed Mayer coupling
        // (for example the two terminal O atoms in CO2), but it must not close
        // a fictitious aromatic cycle.
        if (distance>1.45*radii_bohr) return false;
        for (const auto& record:wf.bond_orders) {
            if ((record.atom_a==a && record.atom_b==b) ||
                (record.atom_a==b && record.atom_b==a)) {
                return record.mayer_order>=0.05;
            }
        }
        return false;
    }

    // Geometry is used only when no electronic connectivity is available at
    // all.  The cutoff is deliberately tighter than the rendering fallback.
    return distance<=1.22*radii_bohr;
}

struct AtomicPReference {
    std::uint32_t atom=0;
    int principal_n=-1;
    std::array<const ReferenceColumn*,3> component{nullptr,nullptr,nullptr};
};

std::vector<AtomicPReference> atomic_valence_p_references(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw) {
    std::map<std::uint32_t,AtomicPReference> by_atom;
    for (const auto& column:raw) {
        if (column.space!=ReferenceClass::ChemicalValence ||
            column.l!=1 || column.component<0 || column.component>=3 ||
            column.atom>=wf.atoms.size() ||
            wf.atoms[column.atom].atomic_number<=2 ||
            transition_metal(wf.atoms[column.atom].atomic_number) ||
            f_block(wf.atoms[column.atom].atomic_number)) {
            continue;
        }
        auto& candidate=by_atom[column.atom];
        candidate.atom=column.atom;
        if (column.n>candidate.principal_n) {
            candidate.principal_n=column.n;
            candidate.component={nullptr,nullptr,nullptr};
        }
        if (column.n==candidate.principal_n) {
            candidate.component[static_cast<std::size_t>(column.component)]=
                &column;
        }
    }

    std::vector<AtomicPReference> result;
    for (const auto& [atom,candidate]:by_atom) {
        (void)atom;
        if (std::all_of(candidate.component.begin(),candidate.component.end(),
                        [](const auto* value){return value!=nullptr;})) {
            result.push_back(candidate);
        }
    }
    return result;
}

struct PiPlane {
    bool valid=false;
    std::array<double,3> origin{0.0,0.0,0.0};
    std::array<double,3> normal{0.0,0.0,1.0};
    double span_bohr=0.0;
    double rms_bohr=0.0;
    double max_deviation_bohr=0.0;
};

std::array<double,3> atom_position(const Atom& atom) {
    return {atom.x,atom.y,atom.z};
}

std::array<double,3> subtract3(const std::array<double,3>& a,
                               const std::array<double,3>& b) {
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

PiPlane fit_pi_plane(const Wavefunction& wf,
                     const std::vector<std::uint32_t>& atoms) {
    PiPlane result;
    if (atoms.size()<3u) return result;
    for (const auto atom:atoms) {
        if (atom>=wf.atoms.size()) return result;
        const auto position=atom_position(wf.atoms[atom]);
        for (std::size_t c=0;c<3u;++c) result.origin[c]+=position[c];
    }
    for (double& value:result.origin) {
        value/=static_cast<double>(atoms.size());
    }

    double best_area=0.0;
    std::array<double,3> best_normal{0.0,0.0,0.0};
    for (std::size_t i=0;i<atoms.size();++i) {
        for (std::size_t j=i+1u;j<atoms.size();++j) {
            result.span_bohr=std::max(result.span_bohr,
                distance_bohr(wf.atoms[atoms[i]],wf.atoms[atoms[j]]));
            for (std::size_t k=j+1u;k<atoms.size();++k) {
                const auto a=atom_position(wf.atoms[atoms[i]]);
                const auto ab=subtract3(atom_position(wf.atoms[atoms[j]]),a);
                const auto ac=subtract3(atom_position(wf.atoms[atoms[k]]),a);
                const auto normal=cross3(ab,ac);
                const double area=std::sqrt(dot3(normal,normal));
                if (area>best_area) {
                    best_area=area;
                    best_normal=normal;
                }
            }
        }
    }
    if (!(result.span_bohr>1.0e-6) ||
        best_area<0.015*result.span_bohr*result.span_bohr) {
        return result;
    }
    for (double& value:best_normal) value/=best_area;
    const auto sign_component=std::max_element(
        best_normal.begin(),best_normal.end(),
        [](const double a,const double b){return std::abs(a)<std::abs(b);});
    if (sign_component!=best_normal.end() && *sign_component<0.0) {
        for (double& value:best_normal) value=-value;
    }
    result.normal=best_normal;

    double sum2=0.0;
    for (const auto atom:atoms) {
        const double deviation=std::abs(dot3(
            subtract3(atom_position(wf.atoms[atom]),result.origin),
            result.normal));
        result.max_deviation_bohr=std::max(
            result.max_deviation_bohr,deviation);
        sum2+=deviation*deviation;
    }
    result.rms_bohr=std::sqrt(sum2/static_cast<double>(atoms.size()));
    const double rms_limit=std::max(0.12,0.035*result.span_bohr);
    const double max_limit=std::max(0.25,0.070*result.span_bohr);
    result.valid=result.rms_bohr<=rms_limit &&
                 result.max_deviation_bohr<=max_limit;
    return result;
}

bool substituents_follow_plane(const Wavefunction& wf,
                               const std::vector<std::uint32_t>& atoms,
                               const PiPlane& plane) {
    std::set<std::uint32_t> checked(atoms.begin(),atoms.end());
    for (const auto a:atoms) {
        for (std::uint32_t b=0;b<wf.atoms.size();++b) {
            if (pi_connectivity_evidence(wf,a,b)) checked.insert(b);
        }
    }
    const double limit=std::max(0.30,0.075*plane.span_bohr);
    for (const auto atom:checked) {
        const double deviation=std::abs(dot3(
            subtract3(atom_position(wf.atoms[atom]),plane.origin),
            plane.normal));
        if (deviation>limit) return false;
    }
    return true;
}

std::vector<std::vector<std::uint32_t>> planar_p_components(
    const Wavefunction& wf,
    const std::vector<AtomicPReference>& p_reference) {
    std::vector<std::vector<std::uint32_t>> result;
    std::set<std::uint32_t> visited;
    for (const auto& seed:p_reference) {
        if (visited.count(seed.atom)!=0u) continue;
        std::vector<std::uint32_t> component;
        std::vector<std::uint32_t> pending{seed.atom};
        visited.insert(seed.atom);
        while (!pending.empty()) {
            const auto atom=pending.back();
            pending.pop_back();
            component.push_back(atom);
            for (const auto& neighbour:p_reference) {
                if (visited.count(neighbour.atom)!=0u ||
                    !pi_connectivity_evidence(wf,atom,neighbour.atom)) {
                    continue;
                }
                visited.insert(neighbour.atom);
                pending.push_back(neighbour.atom);
            }
        }
        if (component.size()<3u) continue;
        std::sort(component.begin(),component.end());
        result.push_back(std::move(component));
    }
    return result;
}

std::vector<ReferenceColumn> perpendicular_p_subspace(
    const Wavefunction& wf,
    const std::vector<std::uint32_t>& atoms,
    const PiPlane& plane,
    const std::vector<AtomicPReference>& p_reference) {
    std::vector<ReferenceColumn> raw_perpendicular;
    raw_perpendicular.reserve(atoms.size());
    for (const auto atom:atoms) {
        const auto found=std::find_if(
            p_reference.begin(),p_reference.end(),
            [atom](const auto& candidate){return candidate.atom==atom;});
        if (found==p_reference.end()) return {};
        ReferenceColumn column;
        column.coefficients.assign(wf.basis_count,0.0);
        column.atom=atom;
        column.n=found->principal_n;
        column.l=1;
        column.component=2;
        column.space=ReferenceClass::ChemicalValence;
        column.label=reference_label(wf,atom,column.n,column.l)+" p_perp";
        for (std::size_t axis=0;axis<3u;++axis) {
            const auto* source=found->component[axis];
            if (source==nullptr ||
                source->coefficients.size()!=column.coefficients.size()) {
                return {};
            }
            for (std::size_t mu=0;mu<column.coefficients.size();++mu) {
                column.coefficients[mu]+=
                    plane.normal[axis]*source->coefficients[mu];
            }
        }
        raw_perpendicular.push_back(std::move(column));
    }
    auto orthonormal=orthonormalise_reference(
        std::move(raw_perpendicular),wf);
    if (orthonormal.size()!=atoms.size()) return {};
    return orthonormal;
}

struct OrientedPColumn {
    std::uint32_t atom=0;
    std::array<double,3> direction{0.0,0.0,1.0};
};

std::vector<ReferenceColumn> oriented_p_subspace(
    const Wavefunction& wf,
    const std::vector<OrientedPColumn>& requested,
    const std::vector<AtomicPReference>& p_reference) {
    std::vector<ReferenceColumn> raw_oriented;
    raw_oriented.reserve(requested.size());
    for (const auto& requested_column:requested) {
        const auto found=std::find_if(
            p_reference.begin(),p_reference.end(),
            [&](const auto& candidate){
                return candidate.atom==requested_column.atom;
            });
        if (found==p_reference.end()) return {};
        ReferenceColumn column;
        column.coefficients.assign(wf.basis_count,0.0);
        column.atom=requested_column.atom;
        column.n=found->principal_n;
        column.l=1;
        column.component=2;
        column.space=ReferenceClass::ChemicalValence;
        column.label=reference_label(
            wf,column.atom,column.n,column.l)+" p_oriented";
        for (std::size_t axis=0;axis<3u;++axis) {
            const auto* source=found->component[axis];
            if (source==nullptr ||
                source->coefficients.size()!=column.coefficients.size()) {
                return {};
            }
            for (std::size_t mu=0;mu<column.coefficients.size();++mu) {
                column.coefficients[mu]+=
                    requested_column.direction[axis]*source->coefficients[mu];
            }
        }
        raw_oriented.push_back(std::move(column));
    }
    auto orthonormal=orthonormalise_reference(std::move(raw_oriented),wf);
    if (orthonormal.size()!=requested.size()) return {};
    return orthonormal;
}

std::vector<std::pair<std::uint32_t,std::uint32_t>>
pi_topology_bonded_pairs(const Wavefunction& wf) {
    std::vector<std::pair<std::uint32_t,std::uint32_t>> result;
    if (wf.bond_order_provenance!=DataProvenance::Unavailable) {
        std::set<std::pair<std::uint32_t,std::uint32_t>> unique;
        for (const auto& record:wf.bond_orders) {
            auto a=record.atom_a;
            auto b=record.atom_b;
            if (a==b || a>=wf.atoms.size() || b>=wf.atoms.size() ||
                record.mayer_order<0.05) {
                continue;
            }
            if (b<a) std::swap(a,b);
            const double radii=(
                covalent_radius_angstrom(wf.atoms[a].atomic_number)+
                covalent_radius_angstrom(wf.atoms[b].atomic_number))*
                kAngstromToBohr;
            if (distance_bohr(wf.atoms[a],wf.atoms[b])>1.45*radii) continue;
            unique.insert({a,b});
        }
        result.assign(unique.begin(),unique.end());
        return result;
    }
    for (std::uint32_t a=0;a<wf.atoms.size();++a) {
        for (std::uint32_t b=a+1u;b<wf.atoms.size();++b) {
            if (pi_connectivity_evidence(wf,a,b)) result.emplace_back(a,b);
        }
    }
    return result;
}

PiTopologyGraphEvidence pi_ring_graph_evidence(
    const Wavefunction& wf,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& bonded_pairs,
    const std::vector<PiOrientationChannel>& channels) {
    PiTopologyGraphEvidence evidence;
    evidence.atom_count=wf.atoms.size();
    evidence.edges=bonded_pairs;
    evidence.bond_order_provenance=wf.bond_order_provenance;
    const bool electronic=wf.bond_order_provenance!=DataProvenance::Unavailable;
    evidence.source=electronic?PiTopologyGraphSource::MayerDistanceModel:
        PiTopologyGraphSource::CovalentDistanceModel;
    evidence.maximum_covalent_radius_factor=electronic?1.45:1.22;
    if(electronic)evidence.minimum_mayer_order=0.05;
    if(channels.size()<2u) {
        evidence.channel_association_search_complete=true;
        return evidence;
    }
    std::vector<std::size_t> degree(wf.atoms.size());
    for(const auto& [a,b]:bonded_pairs) {++degree[a];++degree[b];}
    for(std::uint32_t hub=0;hub<wf.atoms.size();++hub) {
        // Two simple rings sharing only the hub require four distinct
        // incident edges. The graph test makes no element-specific exception.
        if(degree[hub]<4u)continue;
        evidence.examined_hubs.push_back(hub);
        if(find_cycle_union_at_hub(static_cast<std::uint32_t>(wf.atoms.size()),bonded_pairs,hub))
            evidence.unscoped_cycle_union_hubs.push_back(hub);
    }
    for(std::size_t a=0;a<channels.size();++a)
    for(std::size_t b=a+1;b<channels.size();++b)
    for(const auto hub:evidence.unscoped_cycle_union_hubs) {
        const auto found=find_channel_cycle_union_at_hub(
            static_cast<std::uint32_t>(wf.atoms.size()),bonded_pairs,hub,
            {{channels[a].atoms,channels[b].atoms}});
        if(found.witness)evidence.channel_ring_witnesses.push_back({{{a,b}},*found.witness});
    }
    evidence.channel_association_search_complete=true;
    return evidence;
}

bool networks_share_atom(const OrientedPiNetwork& a,
                         const OrientedPiNetwork& b) {
    return std::any_of(a.atoms.begin(),a.atoms.end(),[&](const auto atom){
        return std::find(b.atoms.begin(),b.atoms.end(),atom)!=b.atoms.end();
    });
}

struct PiNetworkBundle {
    std::vector<std::size_t> networks;
    bool haptic_metal=false;
    std::set<std::uint32_t> haptic_metals;
    bool symmetry_direct_sum=false;
    std::uint32_t symmetry_hub=std::numeric_limits<std::uint32_t>::max();
    bool has_non_direct_sum_join=false;
};

std::vector<PiNetworkBundle> pi_network_bundles(
    const Wavefunction& wf,
    const std::vector<OrientedPiNetwork>& networks,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& bonded_pairs) {
    std::vector<std::vector<std::uint32_t>> adjacency(wf.atoms.size());
    for (const auto& pair:bonded_pairs) {
        if (pair.first>=adjacency.size() || pair.second>=adjacency.size()) {
            continue;
        }
        adjacency[pair.first].push_back(pair.second);
        adjacency[pair.second].push_back(pair.first);
    }
    std::vector<std::map<std::uint32_t,std::size_t>> bridge_counts(
        networks.size());
    for (std::size_t network_index=0u;network_index<networks.size();
         ++network_index) {
        const auto& network=networks[network_index];
        const std::set<std::uint32_t> members(
            network.atoms.begin(),network.atoms.end());
        for (const auto atom:network.atoms) {
            if (atom>=adjacency.size()) continue;
            for (const auto neighbour:adjacency[atom]) {
                if (members.count(neighbour)==0u) {
                    ++bridge_counts[network_index][neighbour];
                }
            }
        }
    }
    auto orthogonal_bridge_coupled=[&](const std::size_t a_index,
                            const std::size_t b_index) {
        const auto& a=networks[a_index];
        const auto& b=networks[b_index];
        if (std::abs(dot3(a.representative_direction,
                          b.representative_direction))>0.35) return false;
        const auto& a_bridges=bridge_counts[a_index];
        const auto& b_bridges=bridge_counts[b_index];
        for (const auto& [bridge,a_neighbours]:a_bridges) {
            const auto b_bridge=b_bridges.find(bridge);
            if (a_neighbours<2u || b_bridge==b_bridges.end() ||
                b_bridge->second<2u) {
                continue;
            }
            if (bridge>=wf.atoms.size()) continue;
            const int z=wf.atoms[bridge].atomic_number;
            if (z<=2 || transition_metal(z) || f_block(z)) continue;
            // This proposes an electronic bundle, not a ring classification.
            // Two neighbours in each component do not prove closed rings.
            // Structural ring paths and their channel correspondence are
            // established separately after the active space is selected.
            return true;
        }
        return false;
    };
    std::vector<std::set<std::uint32_t>> network_haptic_metals(
        networks.size());
    for (std::size_t network_index=0u;network_index<networks.size();
         ++network_index) {
        const auto& network=networks[network_index];
        if (!network.cyclic) continue;
        for (std::uint32_t metal=0u;metal<wf.atoms.size();++metal) {
            const int z=wf.atoms[metal].atomic_number;
            if (!transition_metal(z) && !f_block(z)) continue;
            std::size_t contacts=0u;
            for (const auto atom:network.atoms) {
                if (atom>=wf.atoms.size()) continue;
                const double mayer=std::abs(pair_mayer(wf,metal,atom));
                const double radius=(
                    covalent_radius_angstrom(z)+
                    covalent_radius_angstrom(
                        wf.atoms[atom].atomic_number))*kAngstromToBohr;
                const bool within_haptic_shell=
                    distance_bohr(wf.atoms[metal],wf.atoms[atom])<=
                        1.60*radius;
                const bool electronic_contact=
                    within_haptic_shell && mayer>=0.02;
                const bool geometry_fallback=
                    wf.bond_order_provenance==DataProvenance::Unavailable &&
                    within_haptic_shell;
                if (electronic_contact || geometry_fallback) ++contacts;
            }
            const std::size_t required=std::max<std::size_t>(
                3u,(3u*network.atoms.size()+4u)/5u);
            if (contacts>=required) {
                network_haptic_metals[network_index].insert(metal);
            }
        }
    }
    auto haptic_metal_coupled=[&](const std::size_t a_index,
                                  const std::size_t b_index) {
        if (a_index>=network_haptic_metals.size() ||
            b_index>=network_haptic_metals.size() ||
            networks_share_atom(networks[a_index],networks[b_index])) {
            return false;
        }
        const auto& a=network_haptic_metals[a_index];
        const auto& b=network_haptic_metals[b_index];
        return std::any_of(a.begin(),a.end(),[&](const auto metal) {
            return b.count(metal)!=0u;
        });
    };

    struct RootedCycleTraversal {
        std::uint32_t root=std::numeric_limits<std::uint32_t>::max();
        std::vector<std::vector<int>> labels;
        std::vector<std::uint32_t> mapped_atoms;
    };
    auto rooted_cycle_traversals=[&](const OrientedPiNetwork& network,
                                     const std::uint32_t hub) {
        std::vector<RootedCycleTraversal> traversals;
        const std::set<std::uint32_t> members(
            network.atoms.begin(),network.atoms.end());
        if (members.size()<3u || network.edges.size()!=members.size()) {
            return traversals;
        }
        std::map<std::uint32_t,std::vector<std::uint32_t>> internal;
        for (const auto& edge:network.edges) {
            if (members.count(edge.first)==0u ||
                members.count(edge.second)==0u) {
                return traversals;
            }
            internal[edge.first].push_back(edge.second);
            internal[edge.second].push_back(edge.first);
        }
        for (const auto atom:network.atoms) {
            if (internal[atom].size()!=2u) return traversals;
        }
        std::vector<std::uint32_t> roots;
        for (const auto atom:network.atoms) {
            if (atom<adjacency.size() &&
                std::find(adjacency[atom].begin(),adjacency[atom].end(),hub)!=
                    adjacency[atom].end()) {
                roots.push_back(atom);
            }
        }
        // A saturated common hub supplies one root per independent ring.  A
        // multi-contact root is haptic/spiro chemistry and is handled by the
        // corresponding predicates instead of a symmetry direct-sum claim.
        if (roots.size()!=1u) return traversals;
        const auto root=roots.front();
        for (const auto first:internal[root]) {
            std::vector<std::uint32_t> order{root};
            std::set<std::uint32_t> visited{root};
            auto previous=root;
            auto current=first;
            while (current!=root && visited.insert(current).second) {
                order.push_back(current);
                const auto& neighbours=internal[current];
                const auto next=neighbours[0]==previous
                    ?neighbours[1]:neighbours[0];
                previous=current;
                current=next;
            }
            if (current!=root || order.size()!=members.size()) continue;

            RootedCycleTraversal traversal;
            traversal.root=root;
            bool supported=true;
            for (const auto atom:order) {
                if (atom>=wf.atoms.size() || atom>=adjacency.size()) {
                    supported=false;
                    break;
                }
                std::vector<std::uint32_t> external;
                for (const auto neighbour:adjacency[atom]) {
                    if (members.count(neighbour)==0u && neighbour!=hub) {
                        external.push_back(neighbour);
                    }
                }
                // Stay conservative for branched substituents: without a
                // full molecular graph isomorphism they are separate pi
                // families, not asserted symmetry copies.
                if (external.size()>1u ||
                    (!external.empty() &&
                     adjacency[external.front()].size()!=1u)) {
                    supported=false;
                    break;
                }
                std::vector<int> label{wf.atoms[atom].atomic_number};
                traversal.mapped_atoms.push_back(atom);
                if (!external.empty()) {
                    if (external.front()>=wf.atoms.size()) {
                        supported=false;
                        break;
                    }
                    label.push_back(wf.atoms[external.front()].atomic_number);
                    traversal.mapped_atoms.push_back(external.front());
                } else {
                    label.push_back(0);
                }
                traversal.labels.push_back(std::move(label));
            }
            if (supported) traversals.push_back(std::move(traversal));
        }
        return traversals;
    };
    auto symmetry_related_geometry=[&](const RootedCycleTraversal& a,
                                       const RootedCycleTraversal& b,
                                       const std::uint32_t hub) {
        if (hub>=wf.atoms.size() || a.labels!=b.labels ||
            a.mapped_atoms.size()!=b.mapped_atoms.size() ||
            a.mapped_atoms.empty()) {
            return false;
        }
        const auto vector_from_hub=[&](const std::uint32_t atom) {
            return std::array<double,3>{
                wf.atoms[atom].x-wf.atoms[hub].x,
                wf.atoms[atom].y-wf.atoms[hub].y,
                wf.atoms[atom].z-wf.atoms[hub].z};
        };
        const auto gram_equivalent=[&](
            const std::vector<std::uint32_t>& source,
            const std::vector<std::uint32_t>& target) {
            if (source.size()!=target.size() || source.empty()) return false;
            double squared_error=0.0;
            double squared_scale=0.0;
            for (std::size_t i=0u;i<source.size();++i) {
                const auto ai=vector_from_hub(source[i]);
                const auto bi=vector_from_hub(target[i]);
                for (std::size_t j=0u;j<source.size();++j) {
                    const auto aj=vector_from_hub(source[j]);
                    const auto bj=vector_from_hub(target[j]);
                    const double ga=dot3(ai,aj);
                    const double gb=dot3(bi,bj);
                    const double scale=0.5*(std::abs(ga)+std::abs(gb));
                    squared_error+=(ga-gb)*(ga-gb);
                    squared_scale+=scale*scale;
                }
            }
            return squared_scale>1.0e-16 &&
                std::sqrt(squared_error/squared_scale)<=0.06;
        };

        // A local ring-to-ring fit is not enough to claim a molecular
        // symmetry direct sum: the same orthogonal map must also permute the
        // complete first shell of the common hub.  This rejects, for example,
        // C(Ph)2ClF even when the two phenyl fragments are individually
        // congruent.  Non-terminal non-ring hub substituents are deliberately
        // unsupported rather than over-claimed.
        const auto hub_neighbour_label=[&](const std::uint32_t neighbour) {
            std::pair<bool,std::vector<int>> result{
                false,{}};
            if (neighbour>=wf.atoms.size() ||
                neighbour>=adjacency.size()) return result;
            result.second={wf.atoms[neighbour].atomic_number,
                static_cast<int>(adjacency[neighbour].size())};
            std::size_t rooted_networks=0u;
            std::vector<std::vector<int>> best;
            for (std::size_t index=0u;index<networks.size();++index) {
                const auto traversals=rooted_cycle_traversals(
                    networks[index],hub);
                for (const auto& traversal:traversals) {
                    if (traversal.root!=neighbour) continue;
                    ++rooted_networks;
                    if (best.empty() || traversal.labels<best) {
                        best=traversal.labels;
                    }
                }
            }
            if (rooted_networks>0u) {
                result.second.push_back(-3);
                for (const auto& row:best) {
                    result.second.push_back(static_cast<int>(row.size()));
                    result.second.insert(
                        result.second.end(),row.begin(),row.end());
                }
                result.first=true;
            } else if (adjacency[neighbour].size()==1u) {
                result.first=true;
            }
            return result;
        };

        std::vector<std::uint32_t> source=a.mapped_atoms;
        std::vector<std::uint32_t> target=b.mapped_atoms;
        // A symmetry operation that exchanges the two chosen rooted rings
        // also exchanges their hub bonds.
        source.push_back(b.root);
        target.push_back(a.root);

        std::vector<std::uint32_t> remaining;
        for (const auto neighbour:adjacency[hub]) {
            if (neighbour!=a.root && neighbour!=b.root) {
                remaining.push_back(neighbour);
            }
        }
        std::vector<std::vector<int>> labels;
        labels.reserve(remaining.size());
        for (const auto neighbour:remaining) {
            auto label=hub_neighbour_label(neighbour);
            if (!label.first) return false;
            labels.push_back(std::move(label.second));
        }
        std::vector<bool> used(remaining.size(),false);
        const auto search=[&](const auto& self,const std::size_t index)->bool {
            if (index==remaining.size()) {
                return gram_equivalent(source,target);
            }
            source.push_back(remaining[index]);
            for (std::size_t candidate=0u;candidate<remaining.size();
                 ++candidate) {
                if (used[candidate] || labels[candidate]!=labels[index]) {
                    continue;
                }
                used[candidate]=true;
                target.push_back(remaining[candidate]);
                if (self(self,index+1u)) return true;
                target.pop_back();
                used[candidate]=false;
            }
            source.pop_back();
            return false;
        };
        return search(search,0u);
    };
    const auto no_hub=std::numeric_limits<std::uint32_t>::max();
    auto common_hub_cycle_coupled=[&](const std::size_t a_index,
                                      const std::size_t b_index) {
        const auto& a=networks[a_index];
        const auto& b=networks[b_index];
        if (!a.cyclic || !b.cyclic || a.atoms.size()!=b.atoms.size() ||
            networks_share_atom(a,b)) {
            return no_hub;
        }
        // Equivalent ligand-local rings attached to one saturated hub can be
        // mixed by the canonical symmetry of the whole molecule even though
        // they are not one conjugated ring.  Select their direct-sum
        // projector jointly (for example the four C6F5 rings of B(C6F5)4-),
        // while retaining separate cyclic orientation channels.  Requiring a
        // common external hub and cyclic networks excludes independent
        // fragments and every CO ligand-field case.
        for (const auto& [hub,a_neighbours]:bridge_counts[a_index]) {
            const auto found=bridge_counts[b_index].find(hub);
            if (a_neighbours<1u || found==bridge_counts[b_index].end() ||
                found->second<1u || hub>=wf.atoms.size()) {
                continue;
            }
            const int z=wf.atoms[hub].atomic_number;
            if (z<=2 || transition_metal(z) || f_block(z)) continue;
            const auto a_traversals=rooted_cycle_traversals(a,hub);
            const auto b_traversals=rooted_cycle_traversals(b,hub);
            for (const auto& a_traversal:a_traversals) {
                for (const auto& b_traversal:b_traversals) {
                    if (symmetry_related_geometry(
                            a_traversal,b_traversal,hub)) {
                        return hub;
                    }
                }
            }
        }
        return no_hub;
    };
    std::vector<PiNetworkBundle> result;
    std::vector<bool> visited(networks.size(),false);
    for (std::size_t seed=0;seed<networks.size();++seed) {
        if (visited[seed]) continue;
        PiNetworkBundle bundle;
        bundle.haptic_metals=network_haptic_metals[seed];
        bundle.haptic_metal=!bundle.haptic_metals.empty();
        std::vector<std::size_t> pending{seed};
        visited[seed]=true;
        while (!pending.empty()) {
            const auto current=pending.back();
            pending.pop_back();
            bundle.networks.push_back(current);
            for (std::size_t candidate=0;candidate<networks.size();++candidate) {
                if (visited[candidate]) continue;
                const bool shared=networks_share_atom(
                    networks[current],networks[candidate]);
                const bool orthogonal_bridge=orthogonal_bridge_coupled(current,candidate);
                const bool haptic=haptic_metal_coupled(current,candidate);
                const auto hub=common_hub_cycle_coupled(current,candidate);
                const bool hub_coupled=hub!=no_hub &&
                    (bundle.symmetry_hub==no_hub ||
                     bundle.symmetry_hub==hub);
                if (!shared && !orthogonal_bridge && !haptic && !hub_coupled) continue;
                bundle.haptic_metal=bundle.haptic_metal || haptic;
                bundle.haptic_metals.insert(
                    network_haptic_metals[candidate].begin(),
                    network_haptic_metals[candidate].end());
                bundle.haptic_metal=!bundle.haptic_metals.empty();
                if (hub_coupled && bundle.symmetry_hub==no_hub) {
                    bundle.symmetry_hub=hub;
                }
                bundle.has_non_direct_sum_join=
                    bundle.has_non_direct_sum_join || shared || orthogonal_bridge || haptic;
                visited[candidate]=true;
                pending.push_back(candidate);
            }
        }
        bundle.symmetry_direct_sum=bundle.symmetry_hub!=no_hub &&
            !bundle.has_non_direct_sum_join;
        result.push_back(std::move(bundle));
    }
    return result;
}

std::vector<OrientedPColumn> oriented_columns_for_bundle(
    const std::vector<OrientedPiNetwork>& networks,
    const PiNetworkBundle& bundle) {
    std::vector<OrientedPColumn> result;
    for (const auto network_index:bundle.networks) {
        if (network_index>=networks.size()) continue;
        const auto& network=networks[network_index];
        for (const auto atom:network.atoms) {
            const bool duplicate=std::any_of(
                result.begin(),result.end(),[&](const auto& existing){
                    return existing.atom==atom &&
                        std::abs(dot3(existing.direction,
                                      network.representative_direction))>0.985;
                });
            if (!duplicate) {
                result.push_back({atom,network.representative_direction});
            }
        }
    }
    return result;
}

struct PiOrbitalSelection {
    std::vector<std::size_t> orbitals;
    std::vector<double> weights;
    std::vector<std::vector<double>> projection;
    double coverage=0.0;
    double minimum_atom_coverage=0.0;
    double minimum_pivot=0.0;
};

PiOrbitalSelection select_perpendicular_p_orbitals(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& perpendicular,
    const double degeneracy_tolerance,
    const std::set<std::size_t>& excluded) {
    PiOrbitalSelection result;
    const std::size_t target=perpendicular.size();
    if (target==0u || wf.orbitals.empty()) return result;

    const bool has_beta=std::any_of(
        wf.orbitals.begin(),wf.orbitals.end(),
        [](const auto& mo){return mo.spin==Spin::Beta;});
    std::vector<Spin> spins{Spin::Alpha};
    if (has_beta) spins.push_back(Spin::Beta);

    result.weights.assign(wf.orbitals.size(),0.0);
    result.projection.assign(
        wf.orbitals.size(),std::vector<double>(target,0.0));
    std::map<Spin,double> complete_weight;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& mo=wf.orbitals[i];
        if (mo.coefficients.size()!=wf.basis_count) continue;
        for (std::size_t p=0;p<target;++p) {
            const double amplitude=mo_s_dot(
                perpendicular[p].coefficients,mo,
                wf.ao_overlap,wf.basis_count);
            result.projection[i][p]=amplitude;
            result.weights[i]+=amplitude*amplitude;
        }
        complete_weight[mo.spin]+=result.weights[i];
    }
    for (const auto spin:spins) {
        const double weight=complete_weight[spin];
        if (weight<0.80*static_cast<double>(target) ||
            weight>1.20*static_cast<double>(target)) {
            return {};
        }
    }

    auto groups=degeneracy_groups(wf,degeneracy_tolerance);
    std::map<Spin,std::vector<std::size_t>> selected_by_spin;
    for (const auto spin:spins) {
        std::vector<std::size_t> candidates;
        for (std::size_t g=0;g<groups.size();++g) {
            auto& group=groups[g];
            bool eligible=group.spin==spin;
            group.score=0.0;
            for (const auto index:group.orbitals) {
                if (index>=wf.orbitals.size() ||
                    excluded.count(index)!=0u ||
                    (!wf.orbitals[index].chemistry.valence_manifold &&
                     wf.orbitals[index].chemistry.valence_weight<0.10)) {
                    eligible=false;
                    break;
                }
                group.score+=result.weights[index];
            }
            if (eligible && group.orbitals.size()<=target) {
                candidates.push_back(g);
            }
        }

        const double negative=-std::numeric_limits<double>::infinity();
        std::vector<double> dp(target+1u,negative);
        std::vector<std::vector<std::size_t>> chosen(target+1u);
        dp[0]=0.0;
        for (const auto group_index:candidates) {
            const auto size=groups[group_index].orbitals.size();
            for (std::size_t used=target+1u;used-- > size;) {
                if (!std::isfinite(dp[used-size])) continue;
                const double trial=dp[used-size]+groups[group_index].score;
                if (trial>dp[used]+1.0e-12) {
                    dp[used]=trial;
                    chosen[used]=chosen[used-size];
                    chosen[used].push_back(group_index);
                }
            }
        }
        if (!std::isfinite(dp[target])) return {};
        auto& spin_selection=selected_by_spin[spin];
        for (const auto group_index:chosen[target]) {
            spin_selection.insert(
                spin_selection.end(),
                groups[group_index].orbitals.begin(),
                groups[group_index].orbitals.end());
        }

        auto full_rank_selection=[&](const std::vector<std::size_t>& selection) {
            if (selection.size()!=target) return false;
            double spin_weight=0.0;
            for (const auto index:selection) {
                if (result.weights[index]<0.12) return false;
                spin_weight+=result.weights[index];
            }
            if (spin_weight/complete_weight[spin]<0.62) return false;
            for (std::size_t p=0;p<target;++p) {
                double atom_coverage=0.0;
                for (const auto index:selection) {
                    const double amplitude=result.projection[index][p];
                    atom_coverage+=amplitude*amplitude;
                }
                if (atom_coverage<0.30) return false;
            }
            std::vector<std::vector<double>> residuals;
            for (const auto index:selection) {
                residuals.push_back(result.projection[index]);
            }
            for (std::size_t pivot=0;pivot<target;++pivot) {
                std::size_t best=pivot;
                double best_norm2=-1.0;
                for (std::size_t candidate=pivot;candidate<target;++candidate) {
                    const double norm2=std::inner_product(
                        residuals[candidate].begin(),residuals[candidate].end(),
                        residuals[candidate].begin(),0.0);
                    if (norm2>best_norm2) {best_norm2=norm2;best=candidate;}
                }
                if (best_norm2<1.0e-3) return false;
                std::swap(residuals[pivot],residuals[best]);
                const double inverse=1.0/std::sqrt(best_norm2);
                for (double& value:residuals[pivot]) value*=inverse;
                for (std::size_t candidate=pivot+1u;candidate<target;++candidate) {
                    const double projection=std::inner_product(
                        residuals[pivot].begin(),residuals[pivot].end(),
                        residuals[candidate].begin(),0.0);
                    for (std::size_t p=0;p<target;++p) {
                        residuals[candidate][p]-=projection*residuals[pivot][p];
                    }
                }
            }
            return true;
        };

        // A maximum-trace N-set can be rank deficient when bonding and
        // antibonding character is distributed across nearby canonical MOs.
        // Retry with a block-pivoted selection.  Exact degenerate groups remain
        // indivisible, while the residual score depends only on their spanned
        // subspace and is therefore invariant to rotations inside the block.
        if (!full_rank_selection(spin_selection)) {
            spin_selection.clear();
            std::vector<std::vector<double>> basis;
            std::set<std::size_t> used_groups;
            while (spin_selection.size()<target) {
                std::size_t best_group=groups.size();
                double best_score=-1.0;
                std::vector<std::vector<double>> best_vectors;
                for (const auto group_index:candidates) {
                    if (used_groups.count(group_index)!=0u) continue;
                    const auto& orbitals=groups[group_index].orbitals;
                    if (spin_selection.size()+orbitals.size()>target) continue;
                    std::vector<std::vector<double>> trial;
                    double score=0.0;
                    bool independent=true;
                    for (const auto index:orbitals) {
                        auto residual=result.projection[index];
                        for (const auto& vector:basis) {
                            const double overlap=std::inner_product(
                                vector.begin(),vector.end(),residual.begin(),0.0);
                            for (std::size_t p=0;p<target;++p) {
                                residual[p]-=overlap*vector[p];
                            }
                        }
                        for (const auto& vector:trial) {
                            const double overlap=std::inner_product(
                                vector.begin(),vector.end(),residual.begin(),0.0);
                            for (std::size_t p=0;p<target;++p) {
                                residual[p]-=overlap*vector[p];
                            }
                        }
                        const double norm2=std::inner_product(
                            residual.begin(),residual.end(),residual.begin(),0.0);
                        if (norm2<1.0e-8) {independent=false;break;}
                        score+=norm2;
                        const double inverse=1.0/std::sqrt(norm2);
                        for (double& value:residual) value*=inverse;
                        trial.push_back(std::move(residual));
                    }
                    if (!independent) continue;
                    score+=0.01*groups[group_index].score;
                    if (score>best_score+1.0e-12) {
                        best_score=score;
                        best_group=group_index;
                        best_vectors=std::move(trial);
                    }
                }
                if (best_group==groups.size()) break;
                used_groups.insert(best_group);
                spin_selection.insert(
                    spin_selection.end(),groups[best_group].orbitals.begin(),
                    groups[best_group].orbitals.end());
                basis.insert(
                    basis.end(),best_vectors.begin(),best_vectors.end());
            }
            if (!full_rank_selection(spin_selection)) return {};
        }
        std::sort(spin_selection.begin(),spin_selection.end());
        if (spin_selection.size()!=target) return {};
        result.orbitals.insert(
            result.orbitals.end(),spin_selection.begin(),spin_selection.end());
    }
    std::sort(result.orbitals.begin(),result.orbitals.end());
    if (result.orbitals.size()!=target*spins.size()) return {};

    double selected_weight=0.0;
    for (const auto spin:spins) {
        double spin_weight=0.0;
        for (const auto index:selected_by_spin[spin]) {
            if (result.weights[index]<0.12) return {};
            spin_weight+=result.weights[index];
        }
        if (spin_weight/complete_weight[spin]<0.62) return {};
        selected_weight+=spin_weight;
    }
    result.coverage=selected_weight/
        static_cast<double>(target*spins.size());
    if (result.coverage<0.62) return {};

    result.minimum_atom_coverage=std::numeric_limits<double>::infinity();
    for (const auto spin:spins) {
        for (std::size_t p=0;p<target;++p) {
            double atom_coverage=0.0;
            for (const auto index:selected_by_spin[spin]) {
                const double amplitude=result.projection[index][p];
                atom_coverage+=amplitude*amplitude;
            }
            result.minimum_atom_coverage=std::min(
                result.minimum_atom_coverage,atom_coverage);
        }
    }
    if (result.minimum_atom_coverage<0.30) return {};

    result.minimum_pivot=std::numeric_limits<double>::infinity();
    for (const auto spin:spins) {
        // Per-spin pivoted Gram-Schmidt rejects a high-trace but rank-deficient
        // selection. Alpha and beta subspaces are independently complete.
        std::vector<std::vector<double>> residuals;
        residuals.reserve(target);
        for (const auto index:selected_by_spin[spin]) {
            residuals.push_back(result.projection[index]);
        }
        for (std::size_t pivot=0;pivot<target;++pivot) {
            std::size_t best=pivot;
            double best_norm2=-1.0;
            for (std::size_t candidate=pivot;candidate<target;++candidate) {
                const double norm2=std::inner_product(
                    residuals[candidate].begin(),residuals[candidate].end(),
                    residuals[candidate].begin(),0.0);
                if (norm2>best_norm2) {
                    best_norm2=norm2;
                    best=candidate;
                }
            }
            if (best_norm2<1.0e-3) return {};
            std::swap(residuals[pivot],residuals[best]);
            result.minimum_pivot=std::min(result.minimum_pivot,best_norm2);
            const double inverse=1.0/std::sqrt(best_norm2);
            for (double& value:residuals[pivot]) value*=inverse;
            for (std::size_t candidate=pivot+1u;candidate<target;++candidate) {
                const double projection=std::inner_product(
                    residuals[pivot].begin(),residuals[pivot].end(),
                    residuals[candidate].begin(),0.0);
                for (std::size_t p=0;p<target;++p) {
                    residuals[candidate][p]-=
                        projection*residuals[pivot][p];
                }
            }
        }
    }
    return result;
}

std::string delocalised_pi_family_id(
    const std::vector<std::uint32_t>& atoms,
    const std::size_t bundle_index) {
    std::ostringstream out;
    out<<"delocalised-pi";
    for (const auto atom:atoms) out<<':'<<(atom+1u);
    out<<":channel-bundle-"<<(bundle_index+1u);
    return out.str();
}

void attach_planar_p_delocalised_families(
    Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw,
    const OrbitalChemistryOptions& options) {
    const auto p_reference=atomic_valence_p_references(wf,raw);
    if (p_reference.size()<2u) return;

    const auto bonded_pairs=pi_topology_bonded_pairs(wf);
    const auto networks=infer_oriented_pi_networks(wf,bonded_pairs);
    if (networks.empty()) return;
    const auto bundles=pi_network_bundles(wf,networks,bonded_pairs);
    std::set<std::size_t> excluded;

    for (std::size_t bundle_index=0;bundle_index<bundles.size();
         ++bundle_index) {
        const auto& bundle=bundles[bundle_index];
        const auto requested=oriented_columns_for_bundle(networks,bundle);
        if (requested.size()<2u) continue;
        const auto oriented=oriented_p_subspace(wf,requested,p_reference);
        if (oriented.size()!=requested.size()) continue;
        const auto selected=select_perpendicular_p_orbitals(
            wf,oriented,options.degeneracy_tolerance_hartree,excluded);
        const bool has_beta=std::any_of(
            wf.orbitals.begin(),wf.orbitals.end(),
            [](const auto& mo){return mo.spin==Spin::Beta;});
        const std::size_t expected_orbitals=
            oriented.size()*(has_beta?2u:1u);
        if (selected.orbitals.size()!=expected_orbitals) continue;

        bool conflict=false;
        for (const auto index:selected.orbitals) {
            if (index>=wf.orbitals.size()) {
                conflict=true;
                break;
            }
            if (!wf.orbitals[index].chemistry.delocalised_family_id.empty()) {
                conflict=true;
                break;
            }
        }
        if (conflict) continue;

        double electrons=0.0;
        for (const auto index:selected.orbitals) {
            electrons+=static_cast<double>(wf.orbitals[index].occupation);
        }
        const double maximum_electrons=
            2.0*static_cast<double>(oriented.size());
        if (electrons<0.5 || electrons>maximum_electrons+0.20) continue;

        std::set<std::uint32_t> atom_set;
        for (const auto& column:requested) atom_set.insert(column.atom);
        const std::vector<std::uint32_t> atoms(atom_set.begin(),atom_set.end());
        double topology_confidence=0.0;
        double topology_coherence=0.0;
        bool cyclic=false;
        for (const auto network_index:bundle.networks) {
            topology_confidence+=networks[network_index].confidence;
            topology_coherence+=networks[network_index].orientation_coherence;
            cyclic=cyclic || networks[network_index].cyclic;
        }
        topology_confidence/=
            static_cast<double>(bundle.networks.size());
        topology_coherence/=
            static_cast<double>(bundle.networks.size());

        DelocalisedPiAssignment assignment;
        assignment.family_id=delocalised_pi_family_id(atoms,bundle_index);
        assignment.atoms=atoms;
        for (const auto index:selected.orbitals) {
            assignment.orbitals.push_back(static_cast<std::uint32_t>(index));
        }
        assignment.electron_count=electrons;
        for (const auto network_index:bundle.networks) {
            PiOrientationChannel channel;
            channel.atoms=networks[network_index].atoms;
            channel.direction=networks[network_index].representative_direction;
            channel.coherence=networks[network_index].orientation_coherence;
            channel.cyclic=networks[network_index].cyclic;
            assignment.orientation_channels.push_back(std::move(channel));
        }
        assignment.topology_graph=pi_ring_graph_evidence(
            wf,bonded_pairs,assignment.orientation_channels);
        if (bundle.haptic_metal) {
            assignment.topology=DelocalisedPiTopology::HapticMetal;
        } else if (bundle.symmetry_direct_sum) {
            assignment.topology=DelocalisedPiTopology::SymmetryDirectSum;
        } else if (!assignment.topology_graph.channel_ring_witnesses.empty()) {
            assignment.topology=DelocalisedPiTopology::Spiro;
        } else if (bundle.networks.size()>1u) {
            assignment.topology=DelocalisedPiTopology::MultiChannel;
        } else if (cyclic) {
            assignment.topology=DelocalisedPiTopology::Cycle;
        } else {
            std::map<std::uint32_t,std::size_t> degree;
            for (const auto& edge:
                 networks[bundle.networks.front()].edges) {
                ++degree[edge.first];
                ++degree[edge.second];
            }
            std::vector<std::uint32_t> branch_centres;
            for (const auto& [atom,value]:degree) {
                if (value>2u) branch_centres.push_back(atom);
            }
            assignment.topology=branch_centres.size()==1u
                ?DelocalisedPiTopology::BranchedResonance
                :DelocalisedPiTopology::Path;
            if (branch_centres.size()==1u) {
                for (std::size_t p=0u;p<requested.size();++p) {
                    if (requested[p].atom!=branch_centres.front()) continue;
                    struct SpinProjectedOccupation {
                        double projected_electrons=0.0;
                        double coverage=0.0;
                    };
                    std::map<Spin,SpinProjectedOccupation> by_spin;
                    for (const auto orbital:selected.orbitals) {
                        if (orbital>=selected.projection.size() ||
                            p>=selected.projection[orbital].size()) continue;
                        const double amplitude=
                            selected.projection[orbital][p];
                        const double weight=amplitude*amplitude;
                        auto& projected=by_spin[wf.orbitals[orbital].spin];
                        projected.coverage+=weight;
                        projected.projected_electrons+=
                            weight*static_cast<double>(
                                wf.orbitals[orbital].occupation);
                    }
                    double projected_occupation=0.0;
                    bool available=false;
                    for (const auto& [spin,projected]:by_spin) {
                        (void)spin;
                        if (projected.coverage<=1.0e-10) continue;
                        projected_occupation+=
                            projected.projected_electrons/
                            projected.coverage;
                        available=true;
                    }
                    if (available) {
                        assignment.branch_centre_projected_occupation=
                            std::clamp(projected_occupation,0.0,2.0);
                    }
                    break;
                }
            }
        }
        assignment.cyclic_topology=cyclic;
        assignment.plane_normal=
            networks[bundle.networks.front()].representative_direction;
        assignment.plane_rms_bohr=0.0;
        assignment.subspace_coverage=selected.coverage;
        assignment.confidence=std::clamp(
            0.35+0.30*selected.coverage+
            0.15*std::min(1.0,selected.minimum_atom_coverage)+
            0.10*topology_confidence+0.10*topology_coherence,0.0,1.0);
        std::ostringstream rationale;
        rationale<<"connected locally oriented main-group valence-p network; "
                 <<bundle.networks.size()
                 <<" orientation channel(s); full-rank "
                 <<"S-metric active subspace; exact spin- and "
                 <<"degeneracy-preserving selection including virtual members";
        if (cyclic) rationale<<"; globally coherent cyclic p topology";
        else rationale<<"; no aromaticity claim";
        if (!assignment.topology_graph.channel_ring_witnesses.empty()) {
            rationale<<"; explicit skeletal ring paths with channel-specific correspondence"
                     <<" in the stated first-neighbour graph model";
        } else if (bundle.networks.size()>1u) {
            rationale<<"; multiple orientation channels do not by themselves establish a spiro ring union";
        }
        assignment.rationale=rationale.str();
        assignment.provenance=DataProvenance::Derived;

        for (const auto index:selected.orbitals) {
            auto& chemistry=wf.orbitals[index].chemistry;
            chemistry.valence_manifold=true;
            chemistry.family_symbol="pi";
            chemistry.delocalised_participating_atoms=atoms.size();
            chemistry.delocalised_participating_electrons=electrons;
            chemistry.delocalised_participating_atom_indices=atoms;
            chemistry.delocalised_family_orbitals=assignment.orbitals;
            chemistry.delocalised_family_id=assignment.family_id;
            chemistry.delocalised_orientation_channels=
                assignment.orientation_channels.size();
            chemistry.delocalised_cyclic_topology=
                assignment.cyclic_topology;
            chemistry.delocalised_pi_confidence=std::max(
                chemistry.delocalised_pi_confidence,
                assignment.confidence);
            chemistry.delocalised_pi_weight=std::clamp(
                selected.weights[index],0.0,1.0);
            // Preserve the separately computed atom-pair angular character.
            // Family membership must not relabel every member as a pure pi MO.
            chemistry.confidence=std::max(
                chemistry.confidence,assignment.confidence);
            if (chemistry.method.find("oriented p active-subspace")==
                std::string::npos) {
                chemistry.method+=
                    "; oriented p active-subspace family";
            }
            excluded.insert(index);
        }
        wf.delocalised_pi_assignments.push_back(std::move(assignment));
    }
}

const ReferenceColumn* outer_valence_s_reference(
    const std::vector<ReferenceColumn>& raw,
    const std::uint32_t atom) {
    const ReferenceColumn* result=nullptr;
    for (const auto& column:raw) {
        if (column.atom!=atom || column.l!=0 ||
            column.space!=ReferenceClass::ChemicalValence) {
            continue;
        }
        if (result==nullptr || column.n>result->n) result=&column;
    }
    return result;
}

PiOrbitalSelection select_occupied_projector_subspace(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& reference,
    const std::set<std::uint32_t>& active_atoms,
    const double degeneracy_tolerance,
    const std::set<std::size_t>& excluded) {
    PiOrbitalSelection result;
    const std::size_t target=reference.size();
    if (target==0u || wf.orbitals.empty()) return result;

    const bool has_beta=std::any_of(
        wf.orbitals.begin(),wf.orbitals.end(),
        [](const auto& mo){return mo.spin==Spin::Beta;});
    std::vector<Spin> spins{Spin::Alpha};
    if (has_beta) spins.push_back(Spin::Beta);
    result.weights.assign(wf.orbitals.size(),0.0);
    result.projection.assign(
        wf.orbitals.size(),std::vector<double>(target,0.0));
    std::vector<double> domain_weights(wf.orbitals.size(),0.0);
    for (const auto& candidate:wf.multicentre_candidates) {
        if (candidate.orbital_index>=domain_weights.size()) continue;
        for (std::size_t p=0;p<candidate.atoms.size() &&
             p<candidate.participation.size();++p) {
            if (active_atoms.count(candidate.atoms[p])!=0u) {
                domain_weights[candidate.orbital_index]+=candidate.participation[p];
            }
        }
    }
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& mo=wf.orbitals[i];
        if (mo.coefficients.size()!=wf.basis_count) continue;
        for (std::size_t p=0;p<target;++p) {
            const double amplitude=mo_s_dot(
                reference[p].coefficients,mo,
                wf.ao_overlap,wf.basis_count);
            result.projection[i][p]=amplitude;
            result.weights[i]+=amplitude*amplitude;
        }
    }

    auto groups=degeneracy_groups(wf,degeneracy_tolerance);
    std::map<Spin,std::vector<std::size_t>> selected_by_spin;
    for (const auto spin:spins) {
        std::vector<std::size_t> candidates;
        for (std::size_t g=0;g<groups.size();++g) {
            auto& group=groups[g];
            bool eligible=group.spin==spin;
            group.score=0.0;
            for (const auto index:group.orbitals) {
                if (index>=wf.orbitals.size() || excluded.count(index)!=0u ||
                    wf.orbitals[index].occupation<=0.25 ||
                    (!wf.orbitals[index].chemistry.valence_manifold &&
                     wf.orbitals[index].chemistry.valence_weight<0.05)) {
                    eligible=false;
                    break;
                }
                // The active-domain population selects the canonical source
                // span; the signed S-metric projector below verifies that the
                // selected span can actually resolve every equivalent bridge.
                group.score+=domain_weights[index]+0.10*result.weights[index];
            }
            if (eligible && group.orbitals.size()<=target) candidates.push_back(g);
        }

        const double negative=-std::numeric_limits<double>::infinity();
        std::vector<double> dp(target+1u,negative);
        std::vector<std::vector<std::size_t>> chosen(target+1u);
        dp[0]=0.0;
        for (const auto group_index:candidates) {
            const auto size=groups[group_index].orbitals.size();
            for (std::size_t used=target+1u;used-- > size;) {
                if (!std::isfinite(dp[used-size])) continue;
                const double trial=dp[used-size]+groups[group_index].score;
                if (trial>dp[used]+1.0e-12) {
                    dp[used]=trial;
                    chosen[used]=chosen[used-size];
                    chosen[used].push_back(group_index);
                }
            }
        }
        if (!std::isfinite(dp[target])) return {};
        auto& spin_selection=selected_by_spin[spin];
        for (const auto group_index:chosen[target]) {
            spin_selection.insert(
                spin_selection.end(),groups[group_index].orbitals.begin(),
                groups[group_index].orbitals.end());
        }
        std::sort(spin_selection.begin(),spin_selection.end());
        if (spin_selection.size()!=target) return {};
        result.orbitals.insert(
            result.orbitals.end(),spin_selection.begin(),spin_selection.end());
    }
    std::sort(result.orbitals.begin(),result.orbitals.end());

    double selected_weight=0.0;
    result.minimum_atom_coverage=std::numeric_limits<double>::infinity();
    result.minimum_pivot=std::numeric_limits<double>::infinity();
    for (const auto spin:spins) {
        const auto& selected=selected_by_spin[spin];
        double spin_weight=0.0;
        for (const auto index:selected) spin_weight+=result.weights[index];
        if (spin_weight<1.0e-5*static_cast<double>(target)) return {};
        selected_weight+=spin_weight;
        for (std::size_t p=0;p<target;++p) {
            double column_coverage=0.0;
            for (const auto index:selected) {
                const double amplitude=result.projection[index][p];
                column_coverage+=amplitude*amplitude;
            }
            result.minimum_atom_coverage=std::min(
                result.minimum_atom_coverage,column_coverage);
        }

        std::vector<std::vector<double>> residuals;
        residuals.reserve(target);
        for (const auto index:selected) residuals.push_back(result.projection[index]);
        for (std::size_t pivot=0;pivot<target;++pivot) {
            std::size_t best=pivot;
            double best_norm2=-1.0;
            for (std::size_t candidate=pivot;candidate<target;++candidate) {
                const double norm2=std::inner_product(
                    residuals[candidate].begin(),residuals[candidate].end(),
                    residuals[candidate].begin(),0.0);
                if (norm2>best_norm2) {
                    best_norm2=norm2;
                    best=candidate;
                }
            }
            if (best_norm2<1.0e-8) return {};
            std::swap(residuals[pivot],residuals[best]);
            result.minimum_pivot=std::min(result.minimum_pivot,best_norm2);
            const double inverse=1.0/std::sqrt(best_norm2);
            for (double& value:residuals[pivot]) value*=inverse;
            for (std::size_t candidate=pivot+1u;candidate<target;++candidate) {
                const double projection=std::inner_product(
                    residuals[pivot].begin(),residuals[pivot].end(),
                    residuals[candidate].begin(),0.0);
                for (std::size_t p=0;p<target;++p) {
                    residuals[candidate][p]-=
                        projection*residuals[pivot][p];
                }
            }
        }
    }
    if (result.minimum_atom_coverage<1.0e-6) return {};
    result.coverage=selected_weight/
        static_cast<double>(target*spins.size());
    return result;
}

const AtomicPReference* outer_valence_p_reference(
    const std::vector<AtomicPReference>& references,
    const std::uint32_t atom) {
    const auto found=std::find_if(
        references.begin(),references.end(),
        [atom](const auto& reference){return reference.atom==atom;});
    return found==references.end()?nullptr:&*found;
}

ReferenceColumn directed_three_centre_column(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw,
    const std::vector<AtomicPReference>& p_reference,
    const std::uint32_t atom,
    const std::array<double,3>& direction,
    const bool use_hybrid) {
    ReferenceColumn result;
    result.atom=atom;
    result.space=ReferenceClass::ChemicalValence;
    result.coefficients.assign(wf.basis_count,0.0);
    const auto* s=outer_valence_s_reference(raw,atom);
    const auto* p=outer_valence_p_reference(p_reference,atom);
    if (wf.atoms[atom].atomic_number<=2) {
        if (s==nullptr) return {};
        result=*s;
        result.label=reference_label(wf,atom,result.n,0)+" directed-s";
        return result;
    }
    if (p==nullptr) return {};
    result.n=p->principal_n;
    result.l=1;
    result.component=2;
    result.label=reference_label(wf,atom,result.n,1)+
        (use_hybrid?" directed-sp":" directed-p");
    const double p_scale=use_hybrid?std::sqrt(0.65):1.0;
    const double s_scale=use_hybrid?std::sqrt(0.35):0.0;
    if (use_hybrid && s==nullptr) return {};
    for (std::size_t axis=0;axis<3u;++axis) {
        const auto* source=p->component[axis];
        if (source==nullptr ||
            source->coefficients.size()!=result.coefficients.size()) {
            return {};
        }
        for (std::size_t mu=0;mu<result.coefficients.size();++mu) {
            result.coefficients[mu]+=
                p_scale*direction[axis]*source->coefficients[mu];
        }
    }
    if (use_hybrid) {
        for (std::size_t mu=0;mu<result.coefficients.size();++mu) {
            result.coefficients[mu]+=s_scale*s->coefficients[mu];
        }
    }
    return result;
}

using ThreeCentreAtoms=std::array<std::uint32_t,3>;

std::size_t three_centre_middle(const Wavefunction& wf,
                                const ThreeCentreAtoms& atoms) {
    std::size_t best=0u;
    double best_sum=std::numeric_limits<double>::infinity();
    for (std::size_t candidate=0;candidate<3u;++candidate) {
        double sum=0.0;
        for (std::size_t other=0;other<3u;++other) {
            if (other!=candidate) {
                sum+=distance_bohr(wf.atoms[atoms[candidate]],
                                   wf.atoms[atoms[other]]);
            }
        }
        if (sum<best_sum) {
            best_sum=sum;
            best=candidate;
        }
    }
    return best;
}

bool three_centre_geometry_supported(const Wavefunction& wf,
                                     const ThreeCentreAtoms& atoms,
                                     const std::size_t middle,
                                     bool& linear,
                                     bool& hydrogen_bridge) {
    std::array<std::size_t,2> terminal{};
    std::size_t out=0u;
    for (std::size_t i=0;i<3u;++i) {
        if (i!=middle) terminal[out++]=i;
    }
    const auto centre=atoms[middle];
    const auto left=atoms[terminal[0]];
    const auto right=atoms[terminal[1]];
    const auto a=unit_axis(wf.atoms[centre],wf.atoms[left]);
    const auto b=unit_axis(wf.atoms[centre],wf.atoms[right]);
    const double cosine=dot3(a,b);
    linear=cosine<=-0.88;

    const double d_left=distance_bohr(wf.atoms[centre],wf.atoms[left]);
    const double d_right=distance_bohr(wf.atoms[centre],wf.atoms[right]);
    const double d_terminal=distance_bohr(wf.atoms[left],wf.atoms[right]);
    hydrogen_bridge=wf.atoms[centre].atomic_number==1 &&
        d_terminal>1.15*std::max(d_left,d_right) &&
        std::max(d_left,d_right)<=1.40*std::max(1.0e-12,std::min(d_left,d_right));

    const double minimum=std::min({d_left,d_right,d_terminal});
    const double maximum=std::max({d_left,d_right,d_terminal});
    const bool all_hydrogen=
        wf.atoms[atoms[0]].atomic_number==1 &&
        wf.atoms[atoms[1]].atomic_number==1 &&
        wf.atoms[atoms[2]].atomic_number==1;
    const bool symmetric_h3=all_hydrogen && minimum>1.0e-12 &&
        maximum/minimum<=1.12;
    if (symmetric_h3) return true;
    if (hydrogen_bridge) return true;
    if (!linear) return false;

    // Treat the three atoms as a local unit rather than requiring the entire
    // molecule to be triatomic.  Positive, moderate centre--terminal support
    // admits substituted/spectator-bearing 3c4e units; ordinary multiple-bond
    // chains are rejected by their stronger centre--terminal pair orders.
    // A delocalised canonical density may give the remote terminal pair a
    // non-negligible positive Mayer index (notably I3-), so that value is
    // interpreted relative to both centre--terminal links and only after the
    // geometry has established well-separated collinear terminals.
    const double left_mayer=pair_mayer(wf,centre,left);
    const double right_mayer=pair_mayer(wf,centre,right);
    const double terminal_mayer=pair_mayer(wf,left,right);
    const bool separated_terminals=
        d_terminal>=1.65*std::max(d_left,d_right);
    const double weakest_central=std::min(left_mayer,right_mayer);
    return separated_terminals && left_mayer>=0.02 && right_mayer>=0.02 &&
        0.5*(left_mayer+right_mayer)<=0.95 && terminal_mayer>=0.0 &&
        terminal_mayer<0.55*weakest_central;
}

int three_centre_pair_support(const Wavefunction& wf,
                              const ThreeCentreAtoms& atoms) {
    int result=0;
    for (std::size_t i=0;i<3u;++i) {
        for (std::size_t j=i+1u;j<3u;++j) {
            if (pair_mayer(wf,atoms[i],atoms[j])>=0.02) ++result;
        }
    }
    return result;
}

bool has_three_centre_assignment(const Wavefunction& wf,
                                 const ThreeCentreAtoms& atoms) {
    return std::any_of(
        wf.multicentre_assignments.begin(),wf.multicentre_assignments.end(),
        [&](const auto& assignment){
            if (assignment.atoms.size()!=3u) return false;
            auto existing=assignment.atoms;
            std::sort(existing.begin(),existing.end());
            return std::equal(existing.begin(),existing.end(),atoms.begin());
        });
}

std::array<std::uint32_t,2> three_centre_terminals(
    const ThreeCentreAtoms& atoms,
    const std::size_t middle) {
    std::array<std::uint32_t,2> result{};
    std::size_t out=0u;
    for (std::size_t position=0;position<atoms.size();++position) {
        if (position!=middle) result[out++]=atoms[position];
    }
    if (result[1]<result[0]) std::swap(result[0],result[1]);
    return result;
}

ReferenceColumn hydrogen_bridge_bonding_reference(
    const Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw,
    const std::vector<AtomicPReference>& p_reference,
    const ThreeCentreAtoms& atoms,
    const std::size_t middle) {
    ReferenceColumn result;
    result.atom=atoms[middle];
    result.space=ReferenceClass::ChemicalValence;
    result.label="derived hydrogen-bridge bonding projector";
    result.coefficients.assign(wf.basis_count,0.0);

    std::vector<ReferenceColumn> components;
    components.reserve(3u);
    for (std::size_t position=0;position<atoms.size();++position) {
        const auto atom=atoms[position];
        std::array<double,3> direction{1.0,0.0,0.0};
        if (position==middle) {
            std::size_t terminal=position==0u?1u:0u;
            if (terminal==middle) terminal=2u;
            direction=unit_axis(wf.atoms[atom],wf.atoms[atoms[terminal]]);
        } else {
            direction=unit_axis(wf.atoms[atom],wf.atoms[atoms[middle]]);
        }
        auto component=directed_three_centre_column(
            wf,raw,p_reference,atom,direction,position!=middle);
        if (component.coefficients.size()!=wf.basis_count) return {};
        components.push_back(std::move(component));
    }

    // One projector represents one bridge-localised 3c bonding combination.
    // Several symmetry-equivalent bridges are orthogonalised and selected as
    // one joint active subspace below.  This is an analysis-only rotation: the
    // canonical MO coefficients, energies, symmetry labels and indices remain
    // untouched.
    const double scale=1.0/std::sqrt(static_cast<double>(components.size()));
    for (const auto& component:components) {
        for (std::size_t mu=0;mu<wf.basis_count;++mu) {
            result.coefficients[mu]+=scale*component.coefficients[mu];
        }
    }
    return result;
}

bool equivalent_hydrogen_bridge_group(
    const Wavefunction& wf,
    const std::vector<ThreeCentreAtoms>& bridges) {
    if (bridges.size()<2u) return false;
    std::array<double,4> reference{};
    bool have_reference=false;
    for (const auto& atoms:bridges) {
        const auto middle=three_centre_middle(wf,atoms);
        if (wf.atoms[atoms[middle]].atomic_number!=1) return false;
        const auto terminals=three_centre_terminals(atoms,middle);
        std::array<double,2> distances{
            distance_bohr(wf.atoms[atoms[middle]],wf.atoms[terminals[0]]),
            distance_bohr(wf.atoms[atoms[middle]],wf.atoms[terminals[1]])};
        std::array<double,2> mayer{
            pair_mayer(wf,atoms[middle],terminals[0]),
            pair_mayer(wf,atoms[middle],terminals[1])};
        if (mayer[0]<0.02 || mayer[1]<0.02) return false;
        std::sort(distances.begin(),distances.end());
        std::sort(mayer.begin(),mayer.end());
        const std::array<double,4> signature{
            distances[0],distances[1],mayer[0],mayer[1]};
        if (!have_reference) {
            reference=signature;
            have_reference=true;
            continue;
        }
        for (std::size_t i=0;i<signature.size();++i) {
            const double scale=std::max({
                std::abs(reference[i]),std::abs(signature[i]),1.0e-8});
            // A joint rotation is only valid for electronically and
            // geometrically equivalent channels.  Inequivalent bridges fall
            // back to independent analysis instead of being force-partitioned.
            if (std::abs(signature[i]-reference[i])/scale>0.12) return false;
        }
    }
    return true;
}

void derive_directed_three_centre_assignments(
    Wavefunction& wf,
    const std::vector<ReferenceColumn>& raw,
    const OrbitalChemistryOptions& options) {
    std::set<ThreeCentreAtoms> triples;
    for (const auto& candidate:wf.multicentre_candidates) {
        if (candidate.atoms.size()<3u || candidate.participation.size()<3u ||
            candidate.participation[0]+candidate.participation[1]+
                candidate.participation[2]<0.72) {
            continue;
        }
        ThreeCentreAtoms atoms{
            candidate.atoms[0],candidate.atoms[1],candidate.atoms[2]};
        std::sort(atoms.begin(),atoms.end());
        triples.insert(atoms);
    }

    // The density-only bond-analysis pass intentionally records permissive
    // three-atom *candidates*.  Its old strict occupation heuristic also left
    // provisional assignments behind, which could label ordinary tetrahedral
    // X-H/X-F fragments as 3c2e before geometry had been considered.  From
    // this point onward only geometry-qualified assignments (or an explicitly
    // identified shared bridge subspace) are authoritative.  Candidate
    // triples remain available below for full directed S-metric validation.
    wf.multicentre_assignments.erase(
        std::remove_if(
            wf.multicentre_assignments.begin(),
            wf.multicentre_assignments.end(),
            [](const MulticentreAssignment& assignment) {
                return assignment.provenance==DataProvenance::Derived &&
                    !assignment.geometry_qualified_framework &&
                    assignment.source_subspace_id.empty();
            }),
        wf.multicentre_assignments.end());

    // Canonical orbitals can delocalise several symmetry-equivalent bridges
    // over one joint subspace (the two B-H-B bridges in diborane are the
    // standard example).  Recover geometry-qualified hydrogen bridges from
    // atom-pair electronic evidence so a bridge is not lost merely because it
    // is fourth, rather than third, in every individual canonical MO.
    for (std::uint32_t hydrogen=0;hydrogen<wf.atoms.size();++hydrogen) {
        if (wf.atoms[hydrogen].atomic_number!=1) continue;
        std::vector<std::uint32_t> supported_terminals;
        for (std::uint32_t atom=0;atom<wf.atoms.size();++atom) {
            if (atom==hydrogen || wf.atoms[atom].atomic_number==1) continue;
            if (pair_mayer(wf,hydrogen,atom)>=0.02) {
                supported_terminals.push_back(atom);
            }
        }
        for (std::size_t i=0;i<supported_terminals.size();++i) {
            for (std::size_t j=i+1u;j<supported_terminals.size();++j) {
                ThreeCentreAtoms atoms{
                    hydrogen,supported_terminals[i],supported_terminals[j]};
                std::sort(atoms.begin(),atoms.end());
                const auto middle=three_centre_middle(wf,atoms);
                bool linear=false;
                bool hydrogen_bridge=false;
                if (atoms[middle]==hydrogen &&
                    three_centre_geometry_supported(
                        wf,atoms,middle,linear,hydrogen_bridge) &&
                    hydrogen_bridge) {
                    triples.insert(atoms);
                }
            }
        }
    }
    if (triples.empty()) return;

    const auto p_reference=atomic_valence_p_references(wf,raw);
    std::set<std::size_t> excluded;
    for (const auto& assignment:wf.multicentre_assignments) {
        for (const auto orbital:assignment.orbitals) excluded.insert(orbital);
    }

    std::map<std::array<std::uint32_t,2>,std::vector<ThreeCentreAtoms>>
        hydrogen_bridge_groups;
    for (const auto& atoms:triples) {
        const auto middle=three_centre_middle(wf,atoms);
        bool linear=false;
        bool hydrogen_bridge=false;
        if (three_centre_geometry_supported(
                wf,atoms,middle,linear,hydrogen_bridge) && hydrogen_bridge) {
            hydrogen_bridge_groups[three_centre_terminals(atoms,middle)]
                .push_back(atoms);
        }
    }

    std::set<ThreeCentreAtoms> jointly_handled;
    for (const auto& [terminals,bridges]:hydrogen_bridge_groups) {
        if (!equivalent_hydrogen_bridge_group(wf,bridges)) continue;
        if (std::any_of(bridges.begin(),bridges.end(),[&](const auto& atoms){
                return has_three_centre_assignment(wf,atoms) ||
                    three_centre_pair_support(wf,atoms)<2;
            })) {
            continue;
        }

        std::vector<ReferenceColumn> requested;
        requested.reserve(bridges.size());
        for (const auto& atoms:bridges) {
            const auto middle=three_centre_middle(wf,atoms);
            auto column=hydrogen_bridge_bonding_reference(
                wf,raw,p_reference,atoms,middle);
            if (column.coefficients.size()!=wf.basis_count) {
                requested.clear();
                break;
            }
            requested.push_back(std::move(column));
        }
        if (requested.size()!=bridges.size()) continue;
        auto reference=orthonormalise_reference(std::move(requested),wf);
        if (reference.size()!=bridges.size()) continue;
        std::set<std::uint32_t> active_atoms;
        for (const auto& atoms:bridges) {
            active_atoms.insert(atoms.begin(),atoms.end());
        }
        const auto selected=select_occupied_projector_subspace(
            wf,reference,active_atoms,
            options.degeneracy_tolerance_hartree,excluded);
        const bool has_beta=std::any_of(
            wf.orbitals.begin(),wf.orbitals.end(),
            [](const auto& mo){return mo.spin==Spin::Beta;});
        const std::size_t expected_orbitals=
            bridges.size()*(has_beta?2u:1u);
        if (selected.orbitals.size()!=expected_orbitals) continue;

        double joint_electrons=0.0;
        for (const auto orbital:selected.orbitals) {
            joint_electrons+=static_cast<double>(wf.orbitals[orbital].occupation);
        }
        const double expected_electrons=2.0*static_cast<double>(bridges.size());
        if (!near(joint_electrons,expected_electrons,0.20)) continue;

        const double confidence=std::clamp(
            0.50+0.30*selected.coverage+
            0.20*std::min(1.0,selected.minimum_atom_coverage),0.0,1.0);
        std::ostringstream shared_id;
        shared_id<<"equivalent-hydrogen-bridges:"
                 <<(terminals[0]+1u)<<':'<<(terminals[1]+1u)<<":source";
        for (const auto orbital:selected.orbitals) {
            shared_id<<':'<<(orbital+1u);
        }
        for (const auto& atoms:bridges) {
            MulticentreAssignment assignment;
            assignment.kind=MulticentreKind::ThreeCentreTwoElectron;
            assignment.atoms.assign(atoms.begin(),atoms.end());
            for (const auto orbital:selected.orbitals) {
                assignment.orbitals.push_back(
                    static_cast<std::uint32_t>(orbital));
            }
            assignment.electron_count=2.0;
            assignment.source_subspace_id=shared_id.str();
            assignment.source_subspace_electron_count=joint_electrons;
            assignment.source_subspace_fraction=
                1.0/static_cast<double>(bridges.size());
            assignment.confidence=confidence;
            assignment.geometry_qualified_framework=true;
            assignment.rationale=
                "geometry- and Mayer-qualified equivalent hydrogen bridges; "
                "full-rank joint bridge-localised S-metric projector; shared "
                "canonical source subspace with partitioned 2e channels; "
                "canonical orbitals unchanged; no molecule template";
            assignment.provenance=DataProvenance::Derived;
            wf.multicentre_assignments.push_back(std::move(assignment));
            jointly_handled.insert(atoms);
        }
        for (const auto orbital:selected.orbitals) excluded.insert(orbital);
    }

    for (const auto& atoms:triples) {
        if (jointly_handled.count(atoms)!=0u) continue;
        if (has_three_centre_assignment(wf,atoms) ||
            three_centre_pair_support(wf,atoms)<2) {
            continue;
        }
        const std::size_t middle=three_centre_middle(wf,atoms);
        bool linear=false;
        bool hydrogen_bridge=false;
        if (!three_centre_geometry_supported(
                wf,atoms,middle,linear,hydrogen_bridge)) {
            continue;
        }

        std::vector<ReferenceColumn> requested;
        requested.reserve(3u);
        for (std::size_t position=0;position<3u;++position) {
            const auto atom=atoms[position];
            std::array<double,3> direction{1.0,0.0,0.0};
            if (position==middle) {
                std::size_t terminal=position==0u?1u:0u;
                if (terminal==middle) terminal=2u;
                direction=unit_axis(wf.atoms[atom],wf.atoms[atoms[terminal]]);
            } else {
                direction=unit_axis(wf.atoms[atom],wf.atoms[atoms[middle]]);
            }
            const bool hybrid=hydrogen_bridge && !linear &&
                wf.atoms[atom].atomic_number>2;
            auto column=directed_three_centre_column(
                wf,raw,p_reference,atom,direction,hybrid);
            if (column.coefficients.size()!=wf.basis_count) {
                requested.clear();
                break;
            }
            requested.push_back(std::move(column));
        }
        if (requested.size()!=3u) continue;
        auto reference=orthonormalise_reference(std::move(requested),wf);
        if (reference.size()!=3u) continue;
        const auto selected=select_perpendicular_p_orbitals(
            wf,reference,options.degeneracy_tolerance_hartree,excluded);
        const bool has_beta=std::any_of(
            wf.orbitals.begin(),wf.orbitals.end(),
            [](const auto& mo){return mo.spin==Spin::Beta;});
        if (selected.orbitals.size()!=(has_beta?6u:3u)) continue;

        double electrons=0.0;
        for (const auto orbital:selected.orbitals) {
            electrons+=static_cast<double>(wf.orbitals[orbital].occupation);
        }
        MulticentreKind kind=MulticentreKind::Unclassified;
        if (near(electrons,2.0,0.20)) {
            kind=MulticentreKind::ThreeCentreTwoElectron;
        } else if (near(electrons,4.0,0.20)) {
            kind=MulticentreKind::ThreeCentreFourElectron;
        } else {
            continue;
        }

        MulticentreAssignment assignment;
        assignment.kind=kind;
        assignment.atoms.assign(atoms.begin(),atoms.end());
        for (const auto orbital:selected.orbitals) {
            assignment.orbitals.push_back(
                static_cast<std::uint32_t>(orbital));
            excluded.insert(orbital);
        }
        assignment.electron_count=electrons;
        assignment.confidence=std::clamp(
            0.50+0.30*selected.coverage+
            0.20*std::min(1.0,selected.minimum_atom_coverage),0.0,1.0);
        assignment.geometry_qualified_framework=true;
        assignment.rationale=
            "geometry-qualified three-centre framework; full-rank directed "
            "minimal-valence S-metric subspace; exact degeneracy-preserving "
            "occupation count; no molecule template";
        assignment.provenance=DataProvenance::Derived;
        wf.multicentre_assignments.push_back(std::move(assignment));
    }
}

void attach_multicentre_assignments(Wavefunction& wf) {
    for (const auto& assignment:wf.multicentre_assignments) {
        std::string label;
        if (assignment.kind==MulticentreKind::ThreeCentreTwoElectron) label="3c2e";
        else if (assignment.kind==MulticentreKind::ThreeCentreFourElectron) label="3c4e";
        else continue;
        std::vector<std::uint32_t> participating_atoms=assignment.atoms;
        double participating_electrons=assignment.electron_count;
        std::size_t channel_count=1u;
        if (!assignment.source_subspace_id.empty()) {
            participating_atoms.clear();
            participating_electrons=assignment.source_subspace_electron_count;
            channel_count=0u;
            for (const auto& channel:wf.multicentre_assignments) {
                if (channel.source_subspace_id!=assignment.source_subspace_id) continue;
                ++channel_count;
                participating_atoms.insert(
                    participating_atoms.end(),channel.atoms.begin(),channel.atoms.end());
            }
            std::sort(participating_atoms.begin(),participating_atoms.end());
            participating_atoms.erase(
                std::unique(participating_atoms.begin(),participating_atoms.end()),
                participating_atoms.end());
            std::ostringstream shared_label;
            shared_label<<channel_count<<"×"<<label<<" (shared "
                        <<participating_atoms.size()<<"c/"
                        <<std::max(0LL,std::llround(participating_electrons))
                        <<"e source)";
            label=shared_label.str();
        }
        for (const auto index:assignment.orbitals) {
            if (index>=wf.orbitals.size()) continue;
            auto& chemistry=wf.orbitals[index].chemistry;
            chemistry.multicentre_label=label;
            chemistry.multicentre_assignment_available=
                assignment.provenance!=DataProvenance::Unavailable;
            chemistry.multicentre_participating_atoms=participating_atoms.size();
            chemistry.multicentre_participating_electrons=participating_electrons;
            chemistry.multicentre_participating_atom_indices=participating_atoms;
            chemistry.multicentre_channel_count=channel_count;
            chemistry.multicentre_source_subspace_id=
                assignment.source_subspace_id;
            chemistry.multicentre_source_electron_count=
                assignment.source_subspace_electron_count;
            chemistry.multicentre_confidence=std::max(
                chemistry.multicentre_confidence,
                assignment.confidence);
            chemistry.family_symbol=family_symbol(chemistry.channel.dominant);
        }
    }
}

void generic_three_centre_fallback(Wavefunction& wf) {
    std::vector<std::size_t> selected;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].spin==Spin::Alpha &&
            wf.orbitals[i].chemistry.valence_manifold) selected.push_back(i);
    }
    if (selected.size()!=3u) return;

    std::set<std::uint32_t> atoms;
    double electrons=0.0;
    bool sigma=true;
    for (const auto i:selected) {
        electrons+=static_cast<double>(wf.orbitals[i].occupation);
        const auto& chemistry=wf.orbitals[i].chemistry;
        sigma=sigma &&
            chemistry.channel.sigma>=0.70 &&
            chemistry.channel.undetermined<=0.25;
        for (const auto& contribution:chemistry.ao_contributions) {
            if (contribution.weight>=0.05) atoms.insert(contribution.atom_index);
        }
    }
    if (atoms.size()!=3u || !sigma) return;
    int supported=0;
    std::vector<std::uint32_t> centres(atoms.begin(),atoms.end());
    for (std::size_t a=0;a<centres.size();++a) {
        for (std::size_t b=a+1u;b<centres.size();++b) {
            if (pair_mayer(wf,centres[a],centres[b])>=0.02) ++supported;
        }
    }
    if (supported<2) return;

    ThreeCentreAtoms triple{centres[0],centres[1],centres[2]};
    std::sort(triple.begin(),triple.end());
    const auto middle=three_centre_middle(wf,triple);
    bool linear=false;
    bool hydrogen_bridge=false;
    if (!three_centre_geometry_supported(
            wf,triple,middle,linear,hydrogen_bridge)) {
        return;
    }

    std::string label;
    if (near(electrons,2.0,0.20)) label="3c2e";
    else if (near(electrons,4.0,0.20)) label="3c4e";
    else return;
    for (const auto i:selected) {
        auto& chemistry=wf.orbitals[i].chemistry;
        if (!chemistry.multicentre_label.empty()) continue;
        chemistry.multicentre_label=label;
        chemistry.multicentre_participating_atoms=3u;
        chemistry.multicentre_participating_electrons=electrons;
        chemistry.multicentre_participating_atom_indices=centres;
        chemistry.multicentre_confidence=std::min(
            chemistry.confidence,0.45);
        chemistry.family_symbol="sigma";
    }
}

} // namespace

const char* chemistry_status_name(const ChemistryStatus value) noexcept {
    switch (value) {
        case ChemistryStatus::Determined: return "determined";
        case ChemistryStatus::Percentages: return "percentages";
        case ChemistryStatus::Undetermined: return "undetermined";
        case ChemistryStatus::NotApplicable: return "not-applicable";
        default: return "unavailable";
    }
}

const char* orbital_angular_family_name(
    const OrbitalAngularFamily value) noexcept {
    switch (value) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        case OrbitalAngularFamily::NotApplicable: return "not-applicable";
        default: return "mixed";
    }
}

const char* orbital_bonding_role_name(
    const OrbitalBondingRole value) noexcept {
    switch (value) {
        case OrbitalBondingRole::Bonding: return "bonding";
        case OrbitalBondingRole::Antibonding: return "antibonding";
        case OrbitalBondingRole::Nonbonding: return "nonbonding";
        case OrbitalBondingRole::NotApplicable: return "not-applicable";
        default: return "mixed";
    }
}

void derive_orbital_chemistry(
    Wavefunction& wf,
    const OrbitalChemistryOptions& options) {
    wf.delocalised_pi_assignments.clear();
    const std::size_t n=wf.basis_count;
    if (n==0u || wf.orbitals.empty() ||
        wf.ao_overlap.size()!=n*n) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="AO overlap metric unavailable";
        }
        return;
    }

    const auto raw=build_raw_reference(wf);
    const auto reference=orthonormalise_reference(raw,wf);
    if (reference.empty()) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="COV minimal atomic reference could not be built";
        }
        return;
    }

    std::size_t valence_rank=0;
    for (const auto& ref:reference) {
        if (ref.space==ReferenceClass::ChemicalValence) ++valence_rank;
    }
    const auto basis_atom=basis_atom_map(wf);
    if (basis_atom.size()!=n || valence_rank==0u) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="AO-to-atom map or chemical-valence rank unavailable";
        }
        return;
    }

    std::vector<std::vector<double>> amplitudes(
        wf.orbitals.size(),std::vector<double>(reference.size(),0.0));
    std::vector<std::vector<double>> atom_weight(wf.orbitals.size());

    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        auto& mo=wf.orbitals[i];
        auto& chemistry=mo.chemistry;
        chemistry=OrbitalChemistry{};
        chemistry.available=true;
        chemistry.method="COV FCHK S-metric minimal atomic-reference projection";
        atom_weight[i]=atomic_weights(wf,mo,basis_atom);

        std::map<std::tuple<std::uint32_t,int,int>,double> grouped;
        for (std::size_t r=0;r<reference.size();++r) {
            const double amplitude=mo_s_dot(
                reference[r].coefficients,mo,wf.ao_overlap,n);
            amplitudes[i][r]=amplitude;
            const double weight=amplitude*amplitude;
            switch (reference[r].space) {
                case ReferenceClass::DeepCore:
                    chemistry.deep_core_weight+=weight; break;
                case ReferenceClass::Semicore:
                    chemistry.semicore_weight+=weight; break;
                case ReferenceClass::ChemicalValence:
                    chemistry.valence_weight+=weight;
                    grouped[{reference[r].atom,reference[r].n,reference[r].l}]+=weight;
                    break;
            }
        }
        const double assigned=chemistry.deep_core_weight+
                              chemistry.semicore_weight+
                              chemistry.valence_weight;
        chemistry.unresolved_weight=std::clamp(1.0-assigned,0.0,1.0);
        for (const auto& [key,weight]:grouped) {
            OrbitalAOContribution contribution;
            contribution.atom_index=std::get<0>(key);
            contribution.principal_n=std::get<1>(key);
            contribution.angular_momentum=std::get<2>(key);
            contribution.label=reference_label(
                wf,contribution.atom_index,
                contribution.principal_n,
                contribution.angular_momentum);
            contribution.weight=weight;
            chemistry.ao_contributions.push_back(std::move(contribution));
        }
        std::sort(chemistry.ao_contributions.begin(),
                  chemistry.ao_contributions.end(),
                  [](const auto& a,const auto& b){return a.weight>b.weight;});
        chemistry.confidence=std::clamp(
            chemistry.valence_weight+
            0.5*(chemistry.deep_core_weight+chemistry.semicore_weight),
            0.0,1.0);
    }

    select_valence_manifold(
        wf,valence_rank,options.degeneracy_tolerance_hartree);

    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        auto& mo=wf.orbitals[i];
        auto& chemistry=mo.chemistry;
        for (std::uint32_t a=0;a<wf.atoms.size();++a) {
            for (std::uint32_t b=a+1u;b<wf.atoms.size();++b) {
                const double mayer=pair_mayer(wf,a,b);
                const bool geometric=
                    distance_bohr(wf.atoms[a],wf.atoms[b])<=
                    options.pair_distance_ceiling_bohr;
                const bool atoms_present=
                    a<atom_weight[i].size() && b<atom_weight[i].size() &&
                    atom_weight[i][a]>=options.pair_atom_weight_floor &&
                    atom_weight[i][b]>=options.pair_atom_weight_floor;
                if (std::abs(mayer)<options.pair_mayer_floor &&
                    !(geometric && atoms_present)) continue;

                double overlap_character=0.0;
                for (std::size_t mu=0;mu<n;++mu) {
                    if (basis_atom[mu]!=a) continue;
                    for (std::size_t nu=0;nu<n;++nu) {
                        if (basis_atom[nu]!=b) continue;
                        overlap_character+=
                            2.0*static_cast<double>(mo.coefficients[mu])*
                            wf.ao_overlap[mu*n+nu]*
                            static_cast<double>(mo.coefficients[nu]);
                    }
                }

                OrbitalPairInteraction interaction;
                interaction.atom_a=a;
                interaction.atom_b=b;
                interaction.atom_a_label=
                    wf.atoms[a].symbol+std::to_string(a+1u);
                interaction.atom_b_label=
                    wf.atoms[b].symbol+std::to_string(b+1u);
                interaction.total_mayer_index=mayer;
                interaction.overlap_character=overlap_character;
                interaction.occupied_overlap_contribution=
                    static_cast<double>(mo.occupation)*overlap_character;

                RawChannels channels;
                const auto axis=unit_axis(wf.atoms[a],wf.atoms[b]);
                add_atom_channels(a,axis,reference,amplitudes[i],channels);
                add_atom_channels(b,axis,reference,amplitudes[i],channels);
                interaction.channel=finish_channels(channels,options);
                interaction.bonding=pair_bonding(
                    overlap_character,atoms_present,options);
                chemistry.interactions.push_back(std::move(interaction));
            }
        }

        chemistry.channel=aggregate_channels(
            chemistry.interactions,options);
        chemistry.bonding=aggregate_bonding(
            chemistry.interactions);

        if (chemistry.interactions.empty()) {
            if (chemistry.deep_core_weight>=0.85) {
                chemistry.channel.status=ChemistryStatus::NotApplicable;
                chemistry.channel.dominant=OrbitalAngularFamily::NotApplicable;
                chemistry.channel.undetermined=0.0;
                chemistry.bonding.status=ChemistryStatus::NotApplicable;
                chemistry.bonding.dominant=OrbitalBondingRole::NotApplicable;
                chemistry.bonding.undetermined=0.0;
            } else {
                chemistry.note=
                    "No stable atom-pair interaction frame; chemistry remains UND";
            }
        }

        if (!chemistry.valence_manifold &&
            chemistry.valence_weight<0.10 &&
            chemistry.deep_core_weight<0.50) {
            chemistry.note=
                "Outside the selected minimal chemical-valence canonical manifold";
        }
    }

    // Directional valence-p projectors must first be block-orthogonal to all
    // lower-n p core/semicore functions on the same atom.  Without this step,
    // third-row valence references leak strongly into S 2p core MOs, so a
    // chemically complete SO2/thiourea p space can never pass the valence
    // selector even though the topology and canonical orbitals are correct.
    const auto directional_raw=
        core_orthogonalised_valence_p_reference(wf,raw);
    // Directed X--Y--X three-centre references must retain the full raw
    // valence-p span.  Same-atom core orthogonalisation is required for the
    // planar pi selector below, but applying it here removes the symmetry-
    // adapted terminal combinations of heavy 3c4e systems such as XeF2.
    derive_directed_three_centre_assignments(wf,raw,options);
    attach_multicentre_assignments(wf);
    generic_three_centre_fallback(wf);
    attach_planar_p_delocalised_families(wf,directional_raw,options);
}

} // namespace cov
