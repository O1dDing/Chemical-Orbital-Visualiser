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

#include "cov/local_geometry.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/molecule_style.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

bool fchk_extension(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension == ".fch" || extension == ".fchk";
}

void append_input(const std::filesystem::path& input,
                  std::vector<std::filesystem::path>& files) {
    std::error_code ec;
    if (std::filesystem::is_regular_file(input, ec) && fchk_extension(input)) {
        files.push_back(std::filesystem::absolute(input));
        return;
    }
    ec.clear();
    if (!std::filesystem::is_directory(input, ec)) return;
    for (std::filesystem::recursive_directory_iterator it(
             input, std::filesystem::directory_options::skip_permission_denied, ec),
         end;
         it != end; it.increment(ec)) {
        if (ec) {
            ec.clear();
            continue;
        }
        if (it->is_regular_file(ec) && !ec && fchk_extension(it->path())) {
            files.push_back(std::filesystem::absolute(it->path()));
        }
        ec.clear();
    }
}

bool induced_subgraph_has_cycle(
    const std::vector<std::uint32_t>& atoms,
    const std::vector<cov::BondVisual>& bonds) {
    const std::set<std::uint32_t> atom_set(atoms.begin(),atoms.end());
    if (atom_set.size()<3u) return false;
    std::map<std::uint32_t,std::vector<std::uint32_t>> adjacency;
    for (const auto& bond:bonds) {
        if (atom_set.count(bond.atom_a)==0u ||
            atom_set.count(bond.atom_b)==0u) continue;
        adjacency[bond.atom_a].push_back(bond.atom_b);
        adjacency[bond.atom_b].push_back(bond.atom_a);
    }
    std::set<std::uint32_t> visited;
    const auto visit=[&](const std::uint32_t atom,
                         const std::uint32_t parent,
                         const auto& self)->bool {
        visited.insert(atom);
        for (const auto neighbour:adjacency[atom]) {
            if (neighbour==parent) continue;
            if (visited.count(neighbour)!=0u || self(neighbour,atom,self)) {
                return true;
            }
        }
        return false;
    };
    for (const auto atom:atom_set) {
        if (visited.count(atom)==0u &&
            visit(atom,std::numeric_limits<std::uint32_t>::max(),visit)) {
            return true;
        }
    }
    return false;
}

std::string validate_pi_assignments(
    const cov::Wavefunction& wf,
    const std::vector<cov::BondVisual>& bonds) {
    std::set<std::string> family_ids;
    std::map<std::uint32_t,std::string> orbital_owner;
    for (const auto& assignment:wf.delocalised_pi_assignments) {
        if (assignment.family_id.empty() ||
            !family_ids.insert(assignment.family_id).second) {
            return "missing or duplicate delocalised-pi family id";
        }
        const std::set<std::uint32_t> atom_set(
            assignment.atoms.begin(),assignment.atoms.end());
        if (atom_set.size()!=assignment.atoms.size() || atom_set.empty() ||
            *atom_set.rbegin()>=wf.atoms.size()) {
            return "invalid delocalised-pi family atoms";
        }
        std::set<std::uint32_t> orbital_set;
        double electrons=0.0;
        for (const auto orbital:assignment.orbitals) {
            if (orbital>=wf.orbitals.size() ||
                !orbital_set.insert(orbital).second) {
                return "invalid or duplicate delocalised-pi family orbital";
            }
            const auto [owner,inserted]=orbital_owner.emplace(
                orbital,assignment.family_id);
            if (!inserted && owner->second!=assignment.family_id) {
                return "canonical MO double-counted across pi families";
            }
            electrons+=static_cast<double>(wf.orbitals[orbital].occupation);
        }
        if (!std::isfinite(assignment.electron_count) ||
            std::abs(electrons-assignment.electron_count)>2.0e-3) {
            return "delocalised-pi occupation/electron-count mismatch";
        }
        if (assignment.orientation_channels.empty()) {
            return "delocalised-pi family has no orientation channel";
        }
        bool any_cyclic=false;
        for (std::size_t channel_index=0u;
             channel_index<assignment.orientation_channels.size();
             ++channel_index) {
            const auto& channel=assignment.orientation_channels[channel_index];
            const std::set<std::uint32_t> channel_atoms(
                channel.atoms.begin(),channel.atoms.end());
            if (channel_atoms.size()!=channel.atoms.size() ||
                channel_atoms.empty() ||
                !std::includes(atom_set.begin(),atom_set.end(),
                               channel_atoms.begin(),channel_atoms.end())) {
                return "invalid delocalised-pi orientation-channel atoms";
            }
            const double norm=std::sqrt(
                channel.direction[0]*channel.direction[0]+
                channel.direction[1]*channel.direction[1]+
                channel.direction[2]*channel.direction[2]);
            if (!std::isfinite(norm) || std::abs(norm-1.0)>1.0e-5 ||
                !std::isfinite(channel.coherence) ||
                channel.coherence<0.0 || channel.coherence>1.000001) {
                return "invalid delocalised-pi channel direction/coherence";
            }
            if (channel.cyclic &&
                !induced_subgraph_has_cycle(channel.atoms,bonds)) {
                return "cyclic pi claim has no cycle in structural topology";
            }
            any_cyclic=any_cyclic || channel.cyclic;
            for (std::size_t other=0u;other<channel_index;++other) {
                const auto& previous=assignment.orientation_channels[other];
                const bool share_atom=std::any_of(
                    channel.atoms.begin(),channel.atoms.end(),
                    [&](const auto atom){
                        return std::find(previous.atoms.begin(),
                                         previous.atoms.end(),atom)!=
                            previous.atoms.end();
                    });
                if (!share_atom) continue;
                const double direction_dot=std::abs(
                    channel.direction[0]*previous.direction[0]+
                    channel.direction[1]*previous.direction[1]+
                    channel.direction[2]*previous.direction[2]);
                if (direction_dot>0.35) {
                    return "shared-atom pi channels are not independent";
                }
            }
        }
        if (assignment.cyclic_topology!=any_cyclic) {
            return "pi family cyclic flag disagrees with channel topology";
        }
    }
    return {};
}

std::string validate_wavefunction(const cov::Wavefunction& wf) {
    if (wf.source != cov::WavefunctionSource::Fchk) return "source is not FCHK";
    if (wf.atoms.empty()) return "no atoms";
    if (wf.basis_count == 0u) return "zero basis count";
    if (wf.orbitals.empty()) return "no canonical orbitals";

    std::size_t shell_basis_total = 0u;
    for (const auto& shell : wf.shells) shell_basis_total += cov::shell_basis_count(shell);
    if (shell_basis_total != wf.basis_count) return "shell/basis dimension mismatch";

    double occupation_sum = 0.0;
    for (const auto& orbital : wf.orbitals) {
        if (!std::isfinite(orbital.energy_hartree) ||
            !std::isfinite(static_cast<double>(orbital.occupation))) {
            return "non-finite MO energy or occupation";
        }
        if (orbital.occupation < -1.0e-5f || orbital.occupation > 2.0001f) {
            return "out-of-range occupation";
        }
        if (orbital.coefficients.size() != wf.basis_count) {
            return "MO coefficient dimension mismatch";
        }
        occupation_sum += static_cast<double>(orbital.occupation);
    }
    const double electron_count = static_cast<double>(
        wf.alpha_electrons + wf.beta_electrons);
    if (std::abs(occupation_sum - electron_count) > 2.0e-3) {
        return "occupation/electron-count mismatch";
    }

    // A rectangular canonical-MO block cannot uniquely reconstruct S. Such a
    // producer file is still a valid raw-MO input; only a present matrix is
    // required to be complete and finite. S-dependent chemistry must remain
    // unavailable rather than turning the entire file into a parser failure.
    if (!wf.ao_overlap.empty()) {
        if (wf.ao_overlap.size() != wf.basis_count * wf.basis_count) {
            return "partial AO overlap matrix";
        }
        if (!std::all_of(wf.ao_overlap.begin(), wf.ao_overlap.end(),
                         [](double value) { return std::isfinite(value); })) {
            return "non-finite AO overlap";
        }
    }
    if (electron_count > 0.0 && wf.total_density_packed.empty()) {
        return "electron-bearing wavefunction has no total density";
    }

    const auto bonds = cov::analyse_bonds(wf);
    for (const auto& bond : bonds) {
        if (bond.atom_a >= wf.atoms.size() || bond.atom_b >= wf.atoms.size() ||
            bond.atom_a == bond.atom_b || !std::isfinite(bond.distance_bohr)) {
            return "invalid rendered bond graph";
        }
    }
    if (const auto failure=validate_pi_assignments(wf,bonds); !failure.empty()) {
        return failure;
    }
    for (const auto& geometry : cov::analyse_local_molecular_geometries(wf)) {
        if (!geometry.available() || geometry.centre_atom >= wf.atoms.size()) {
            return "invalid local geometry result";
        }
        for (const auto atom : geometry.neighbour_atoms) {
            if (atom >= wf.atoms.size() || atom == geometry.centre_atom) {
                return "invalid local geometry neighbour";
            }
        }
    }

    cov::MODiagramOptions diagram_options;
    diagram_options.selected_index = 0u;
    const auto diagram = cov::build_mo_diagram_data(wf, diagram_options);
    if (diagram.metadata.size() != wf.orbitals.size()) {
        return "MO diagram metadata lost canonical orbitals";
    }
    for (const auto index : diagram.selection.included_indices) {
        if (index >= wf.orbitals.size()) return "MO diagram selected invalid index";
    }
    return {};
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: cov_fchk_batch_regression <file-or-directory> [...]\n";
        return 2;
    }
    std::vector<std::filesystem::path> files;
    for (int i = 1; i < argc; ++i) append_input(argv[i], files);
    std::sort(files.begin(), files.end());
    files.erase(std::unique(files.begin(), files.end()), files.end());
    if (files.empty()) {
        std::cerr << "no FCHK/FCH files found\n";
        return 2;
    }

    std::size_t passed = 0u;
    std::size_t failed = 0u;
    for (const auto& path : files) {
        const auto start = std::chrono::steady_clock::now();
        try {
            cov::WavefunctionParseOptions options;
            options.max_atoms = 256u;
            options.require_orbitals = true;
            options.keep_density = true;
            options.reconstruct_density_if_missing = true;
            const auto wf = cov::parse_wavefunction(path, options);
            const std::string failure = validate_wavefunction(wf);
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            if (failure.empty()) {
                ++passed;
                std::cout << "PASS\t" << seconds << "\t" << path.string() << '\n';
            } else {
                ++failed;
                std::cout << "FAIL\t" << seconds << "\t" << path.string()
                          << "\t" << failure << '\n';
            }
        } catch (const std::exception& error) {
            ++failed;
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            std::cout << "FAIL\t" << seconds << "\t" << path.string()
                      << "\texception: " << error.what() << '\n';
        }
    }
    std::cout << "SUMMARY\tpassed=" << passed << "\tfailed=" << failed
              << "\ttotal=" << files.size() << '\n';
    return failed == 0u ? EXIT_SUCCESS : EXIT_FAILURE;
}
