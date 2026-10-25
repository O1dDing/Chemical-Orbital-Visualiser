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

#include "cov/molden_parser.hpp"
#include "cov/density.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

struct RawPrimitive {
    double exponent = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
};

struct RawShell {
    std::size_t atom_index = 0;
    std::string type;
    std::vector<RawPrimitive> primitives;
};

struct BasisMarkerState {
    bool d_declared = false;
    bool f_declared = false;
    bool g_declared = false;
};

struct BasisConvention {
    bool pure_d = false;
    bool pure_f = false;
    bool pure_g = false;
};

enum class Section {
    None,
    Atoms,
    GTO,
    MO,
};

std::string trim(std::string s) {
    const auto is_ws = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), is_ws));
    s.erase(std::find_if_not(s.rbegin(), s.rend(), is_ws).base(), s.end());
    return s;
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool starts_with_ci(const std::string& s, const std::string& prefix) {
    if (s.size() < prefix.size()) {
        return false;
    }
    for (std::size_t i = 0; i < prefix.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(s[i])) !=
            std::tolower(static_cast<unsigned char>(prefix[i]))) {
            return false;
        }
    }
    return true;
}

double parse_fortran_double(std::string token) {
    for (char& c : token) {
        if (c == 'D' || c == 'd') {
            c = 'E';
        }
    }
    return std::stod(token);
}

std::uint8_t angular_momentum(const std::string& type) {
    if (type == "s") return 0;
    if (type == "p") return 1;
    if (type == "d") return 2;
    if (type == "f") return 3;
    if (type == "g") return 4;
    throw std::runtime_error("Unsupported Molden shell type: " + type);
}

bool is_section_header(const std::string& line) {
    const auto t = trim(line);
    return t.size() >= 3 && t.front() == '[' && t.back() == ']';
}

void apply_basis_marker(const std::string& line,
                        bool& pure_d,
                        bool& pure_f,
                        bool& pure_g,
                        BasisMarkerState& seen) {
    const std::string t = lower(trim(line));
    if (t == "[5d]") {
        pure_d = true;
        seen.d_declared = true;
    } else if (t == "[6d]") {
        pure_d = false;
        seen.d_declared = true;
    } else if (t == "[7f]") {
        pure_f = true;
        seen.f_declared = true;
    } else if (t == "[10f]") {
        pure_f = false;
        seen.f_declared = true;
    } else if (t == "[9g]") {
        pure_g = true;
        seen.g_declared = true;
    } else if (t == "[15g]") {
        pure_g = false;
        seen.g_declared = true;
    } else if (t == "[5d7f]") {
        pure_d = true;
        pure_f = true;
        seen.d_declared = true;
        seen.f_declared = true;
    } else if (t == "[5d10f]") {
        pure_d = true;
        pure_f = false;
        seen.d_declared = true;
        seen.f_declared = true;
    } else if (t == "[6d7f]") {
        pure_d = false;
        pure_f = true;
        seen.d_declared = true;
        seen.f_declared = true;
    } else if (t == "[6d10f]") {
        pure_d = false;
        pure_f = false;
        seen.d_declared = true;
        seen.f_declared = true;
    }
}

bool pure_for_l(const std::uint8_t l, const Wavefunction& wf) {
    if (l == 2) return wf.pure_d;
    if (l == 3) return wf.pure_f;
    if (l == 4) return wf.pure_g;
    return false;
}

std::uint32_t raw_basis_count(const std::vector<RawShell>& raw_shells,
                              const BasisConvention& convention) {
    std::uint32_t count = 0;
    for (const auto& shell : raw_shells) {
        if (shell.type == "s") {
            count += 1u;
        } else if (shell.type == "p") {
            count += 3u;
        } else if (shell.type == "sp") {
            count += 4u;
        } else if (shell.type == "d") {
            count += convention.pure_d ? 5u : 6u;
        } else if (shell.type == "f") {
            count += convention.pure_f ? 7u : 10u;
        } else if (shell.type == "g") {
            count += convention.pure_g ? 9u : 15u;
        }
    }
    return count;
}

std::uint32_t probe_first_mo_dimension(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        return 0u;
    }

    bool in_mo = false;
    std::uint32_t coefficient_count = 0u;
    std::uint32_t maximum_index = 0u;
    std::string line;
    while (std::getline(input, line)) {
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }
        if (!in_mo) {
            if (starts_with_ci(t, "[mo]")) {
                in_mo = true;
            }
            continue;
        }
        if (is_section_header(t)) {
            break;
        }
        if ((starts_with_ci(t, "sym=") || starts_with_ci(t, "ene=")) &&
            coefficient_count > 0u) {
            break;
        }

        std::istringstream css(t);
        int coefficient_index = 0;
        std::string coefficient_s;
        if (css >> coefficient_index >> coefficient_s && coefficient_index > 0) {
            ++coefficient_count;
            maximum_index = std::max(maximum_index,
                                     static_cast<std::uint32_t>(coefficient_index));
        }
    }

    // Only use a complete contiguous 1..N first MO as dimensional evidence.
    // The full parser remains responsible for duplicate and completeness checks.
    return coefficient_count == maximum_index ? maximum_index : 0u;
}

void infer_omitted_basis_markers(const std::filesystem::path& path,
                                 const std::vector<RawShell>& raw_shells,
                                 const BasisMarkerState& seen,
                                 Wavefunction& wf) {
    const std::uint32_t mo_dimension = probe_first_mo_dimension(path);
    if (mo_dimension == 0u) {
        return;
    }

    const BasisConvention declared{wf.pure_d, wf.pure_f, wf.pure_g};
    if (raw_basis_count(raw_shells, declared) == mo_dimension) {
        return;
    }

    const bool has_d = std::any_of(raw_shells.begin(), raw_shells.end(),
                                   [](const RawShell& shell) { return shell.type == "d"; });
    const bool has_f = std::any_of(raw_shells.begin(), raw_shells.end(),
                                   [](const RawShell& shell) { return shell.type == "f"; });
    const bool has_g = std::any_of(raw_shells.begin(), raw_shells.end(),
                                   [](const RawShell& shell) { return shell.type == "g"; });

    std::vector<BasisConvention> matches;
    for (int mask = 0; mask < 8; ++mask) {
        const BasisConvention candidate{
            (seen.d_declared || !has_d) ? wf.pure_d : ((mask & 1) != 0),
            (seen.f_declared || !has_f) ? wf.pure_f : ((mask & 2) != 0),
            (seen.g_declared || !has_g) ? wf.pure_g : ((mask & 4) != 0),
        };

        const auto duplicate = std::find_if(matches.begin(), matches.end(),
            [&](const BasisConvention& existing) {
                return existing.pure_d == candidate.pure_d &&
                       existing.pure_f == candidate.pure_f &&
                       existing.pure_g == candidate.pure_g;
            });
        if (duplicate != matches.end()) {
            continue;
        }
        if (raw_basis_count(raw_shells, candidate) == mo_dimension) {
            matches.push_back(candidate);
        }
    }

    // Some writers (observed with Multiwfn-generated Gaussian/def2 Molden files)
    // emit [5D] while omitting [7F], even though their MO vector length proves
    // that the f shell is pure spherical. Infer only when coefficient dimension
    // yields one unique convention for shell families that actually exist,
    // while respecting every explicit marker. Ambiguous/inconsistent files
    // remain strict and fail later as before.
    if (matches.size() == 1u) {
        wf.pure_d = matches.front().pure_d;
        wf.pure_f = matches.front().pure_f;
        wf.pure_g = matches.front().pure_g;
    }
}

std::string value_after_equals(const std::string& line) {
    const auto pos = line.find('=');
    if (pos == std::string::npos) {
        return {};
    }
    return trim(line.substr(pos + 1));
}

void finalise_mo(MolecularOrbital& mo,
                 bool& have_mo,
                 std::uint32_t& coefficient_count,
                 Wavefunction& wf,
                 const bool require_coefficients) {
    if (!have_mo) {
        return;
    }
    if (require_coefficients && mo.coefficients.size() != wf.basis_count) {
        throw std::runtime_error("Internal MO coefficient buffer size mismatch");
    }
    if (require_coefficients && coefficient_count != wf.basis_count) {
        throw std::runtime_error(
            "MO contains " + std::to_string(coefficient_count) +
            " coefficient entries, but the parsed basis contains " +
            std::to_string(wf.basis_count) +
            " functions. Check Cartesian/spherical shell convention or file completeness.");
    }
    if (mo.spin_provenance != DataProvenance::Unavailable) {
        mo.source_orbital_index = static_cast<std::size_t>(std::count_if(
            wf.orbitals.begin(), wf.orbitals.end(), [&](const auto& other) {
                return other.spin_provenance != DataProvenance::Unavailable && other.spin == mo.spin;
            }));
    }
    wf.orbitals.push_back(std::move(mo));
    mo = MolecularOrbital{};
    have_mo = false;
    coefficient_count = 0;
}

} // namespace

Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Unable to open Molden file: " + path.string());
    }

    Wavefunction wf;
    wf.source = WavefunctionSource::Molden;
    std::vector<RawShell> raw_shells;
    BasisMarkerState basis_markers;

    Section section = Section::None;
    bool atoms_in_angstrom = false;
    std::size_t current_atom = static_cast<std::size_t>(-1);

    std::string line;
    while (std::getline(input, line)) {
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }

        apply_basis_marker(t, wf.pure_d, wf.pure_f, wf.pure_g, basis_markers);

        if (starts_with_ci(t, "[atoms]")) {
            section = Section::Atoms;
            atoms_in_angstrom = lower(t).find("angs") != std::string::npos;
            continue;
        }
        if (starts_with_ci(t, "[gto]")) {
            section = Section::GTO;
            continue;
        }
        if (starts_with_ci(t, "[mo]")) {
            section = Section::MO;
            continue;
        }
        if (is_section_header(t)) {
            section = Section::None;
            continue;
        }

        if (section == Section::Atoms) {
            std::istringstream iss(t);
            Atom atom;
            int atom_index = 0;
            if (!(iss >> atom.symbol >> atom_index >> atom.atomic_number >> atom.x >> atom.y >> atom.z)) {
                throw std::runtime_error("Malformed [Atoms] line: " + t);
            }
            if (atoms_in_angstrom) {
                atom.x *= kAngstromToBohr;
                atom.y *= kAngstromToBohr;
                atom.z *= kAngstromToBohr;
            }
            wf.atoms.push_back(atom);
            if (wf.atoms.size() > options.max_atoms) {
                throw std::runtime_error(
                    "Molecule contains more than " + std::to_string(options.max_atoms) +
                    " atoms; this build is intentionally targeted at small/medium systems.");
            }
            continue;
        }

        if (section == Section::GTO) {
            std::istringstream iss(t);
            std::string first;
            iss >> first;
            const std::string first_l = lower(first);

            const bool shell_header =
                first_l == "s" || first_l == "p" || first_l == "sp" ||
                first_l == "d" || first_l == "f" || first_l == "g";

            if (!shell_header) {
                // Atom selector line: "<atom index> 0"
                try {
                    const int atom_index = std::stoi(first);
                    if (atom_index <= 0) {
                        throw std::runtime_error("Invalid atom index in [GTO]: " + t);
                    }
                    current_atom = static_cast<std::size_t>(atom_index - 1);
                } catch (const std::invalid_argument&) {
                    throw std::runtime_error("Malformed [GTO] atom selector: " + t);
                }
                continue;
            }

            if (current_atom == static_cast<std::size_t>(-1) ||
                current_atom >= wf.atoms.size()) {
                throw std::runtime_error("Shell encountered before a valid [GTO] atom selector");
            }

            int primitive_count = 0;
            double scale = 1.0;
            if (!(iss >> primitive_count)) {
                throw std::runtime_error("Malformed [GTO] shell header: " + t);
            }
            if (iss >> scale) {
                (void)scale;
            }
            if (primitive_count <= 0) {
                throw std::runtime_error("Non-positive primitive count in shell: " + t);
            }

            RawShell raw;
            raw.atom_index = current_atom;
            raw.type = first_l;
            raw.primitives.reserve(static_cast<std::size_t>(primitive_count));

            for (int i = 0; i < primitive_count; ++i) {
                std::string primitive_line;
                if (!std::getline(input, primitive_line)) {
                    throw std::runtime_error("Unexpected EOF inside [GTO] primitive block");
                }
                primitive_line = trim(primitive_line);
                if (primitive_line.empty()) {
                    --i;
                    continue;
                }
                std::istringstream pss(primitive_line);
                std::string exp_s, c1_s, c2_s;
                if (!(pss >> exp_s >> c1_s)) {
                    throw std::runtime_error("Malformed primitive line: " + primitive_line);
                }
                RawPrimitive primitive;
                primitive.exponent = parse_fortran_double(exp_s) * scale * scale;
                primitive.c1 = parse_fortran_double(c1_s);
                if (raw.type == "sp") {
                    if (!(pss >> c2_s)) {
                        throw std::runtime_error("SP shell primitive lacks p coefficient: " + primitive_line);
                    }
                    primitive.c2 = parse_fortran_double(c2_s);
                }
                raw.primitives.push_back(primitive);
            }
            raw_shells.push_back(std::move(raw));
        }
    }

    if (wf.atoms.empty()) {
        throw std::runtime_error("No atoms found in Molden file");
    }
    if (raw_shells.empty()) {
        throw std::runtime_error("No basis shells found in Molden file");
    }

    infer_omitted_basis_markers(path, raw_shells, basis_markers, wf);

    // Expand SP shells and assign basis offsets after all explicit markers and
    // any uniquely inferable omitted higher-shell convention are known.
    std::uint32_t basis_offset = 0;
    for (const RawShell& raw : raw_shells) {
        const auto append_shell = [&](const std::string& type, const bool use_second_coeff) {
            Shell shell;
            shell.atom_index = static_cast<std::uint32_t>(raw.atom_index);
            shell.angular_momentum = angular_momentum(type);
            shell.pure = pure_for_l(shell.angular_momentum, wf) ? 1u : 0u;
            shell.primitive_offset = static_cast<std::uint32_t>(wf.primitives.size());
            shell.primitive_count = static_cast<std::uint32_t>(raw.primitives.size());
            shell.basis_offset = basis_offset;

            for (const RawPrimitive& p : raw.primitives) {
                Primitive primitive;
                primitive.exponent = p.exponent;
                primitive.coefficient =
                    use_second_coeff ? p.c2 : p.c1;
                wf.primitives.push_back(primitive);
            }
            basis_offset += shell_basis_count(shell);
            wf.shells.push_back(shell);
        };

        if (raw.type == "sp") {
            append_shell("s", false);
            append_shell("p", true);
        } else {
            append_shell(raw.type, false);
        }
    }
    wf.basis_count = basis_offset;

    // Second streaming pass: MO metadata and coefficients.
    input.clear();
    input.seekg(0, std::ios::beg);
    section = Section::None;

    MolecularOrbital current_mo;
    bool have_mo = false;
    std::uint32_t coefficient_count = 0;
    std::vector<std::uint8_t> coefficient_seen(wf.basis_count, 0u);

    const auto begin_mo = [&]() {
        current_mo = MolecularOrbital{};
        current_mo.coefficients.assign(wf.basis_count, 0.0f);
        std::fill(coefficient_seen.begin(), coefficient_seen.end(), std::uint8_t{0});
        coefficient_count = 0;
        have_mo = true;
    };

    while (std::getline(input, line)) {
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }
        if (starts_with_ci(t, "[mo]")) {
            section = Section::MO;
            continue;
        }
        if (section != Section::MO) {
            continue;
        }
        if (is_section_header(t)) {
            break;
        }

        if (starts_with_ci(t, "sym=")) {
            // Sym= is the traditional Molden MO delimiter. Some writers put
            // Ene= first, so only close an existing MO after coefficients have
            // actually been read; otherwise this metadata belongs to the same MO.
            if (have_mo && coefficient_count > 0) {
                finalise_mo(current_mo, have_mo, coefficient_count, wf, true);
            }
            if (!have_mo) begin_mo();
            current_mo.symmetry = value_after_equals(t);
            if (!current_mo.symmetry.empty()) {
                current_mo.symmetry_provenance = DataProvenance::Producer;
            }
            continue;
        }

        if (starts_with_ci(t, "ene=")) {
            // Multiwfn and other valid Molden writers may omit Sym= entirely.
            // In that dialect each new Ene= after a complete coefficient block
            // is the next MO boundary. Without this boundary, coefficient index
            // 1 from the next orbital was falsely reported as a duplicate.
            if (have_mo && coefficient_count > 0) {
                finalise_mo(current_mo, have_mo, coefficient_count, wf, true);
            }
            if (!have_mo) begin_mo();
            current_mo.energy_hartree = parse_fortran_double(value_after_equals(t));
            continue;
        }

        if (!have_mo) {
            // Be permissive about metadata ordering while remaining strict about
            // coefficient completeness and duplicate indices inside one MO.
            if (starts_with_ci(t, "spin=") || starts_with_ci(t, "occup=")) {
                begin_mo();
            } else {
                continue;
            }
        }

        if (starts_with_ci(t, "spin=")) {
            current_mo.spin_source_text = value_after_equals(t);
            const auto spin = lower(trim(current_mo.spin_source_text));
            current_mo.spin_provenance = DataProvenance::Unavailable;
            if (spin == "alpha" || spin == "beta") {
                current_mo.spin = spin == "beta" ? Spin::Beta : Spin::Alpha;
                current_mo.spin_provenance = DataProvenance::Producer;
            }
        } else if (starts_with_ci(t, "occup=")) {
            current_mo.occupation =
                parse_fortran_double(value_after_equals(t));
            current_mo.occupation_provenance = DataProvenance::Producer;
        } else {
            std::istringstream css(t);
            int coefficient_index = 0;
            std::string coefficient_s;
            if (css >> coefficient_index >> coefficient_s) {
                if (coefficient_index <= 0 ||
                    static_cast<std::uint32_t>(coefficient_index) > wf.basis_count) {
                    throw std::runtime_error(
                        "MO coefficient index " + std::to_string(coefficient_index) +
                        " exceeds parsed basis count " + std::to_string(wf.basis_count) +
                        ". Check Cartesian/spherical shell convention.");
                }
                const std::size_t ci = static_cast<std::size_t>(coefficient_index - 1);
                if (coefficient_seen[ci]) {
                    throw std::runtime_error(
                        "Duplicate MO coefficient index " + std::to_string(coefficient_index));
                }
                coefficient_seen[ci] = 1u;
                ++coefficient_count;
                current_mo.coefficients[ci] = parse_fortran_double(coefficient_s);
            }
        }
    }
    finalise_mo(current_mo, have_mo, coefficient_count, wf, true);

    if (options.require_orbitals && wf.orbitals.empty()) {
        throw std::runtime_error("No molecular orbitals found in Molden file");
    }

    for (const auto& mo : wf.orbitals) {
        if (mo.coefficients.size() != wf.basis_count) {
            throw std::runtime_error("MO/basis consistency check failed");
        }
    }

    derive_molden_electron_counts(wf);

    return wf;
}

} // namespace cov
