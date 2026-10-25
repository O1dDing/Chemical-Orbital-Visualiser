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

#include "cov/orbital_tracking.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <set>
#include <string>
#include <vector>

namespace {

struct ManifestFrame {
    int frame = 0;
    std::filesystem::path fch_path;
};

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0u; index < line.size(); ++index) {
        const char character = line[index];
        if (character == '"') {
            if (quoted && index + 1u < line.size() && line[index + 1u] == '"') {
                field.push_back('"');
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (character == '\t' && !quoted) {
            fields.push_back(std::move(field));
            field.clear();
        } else if (character != '\r') {
            field.push_back(character);
        }
    }
    fields.push_back(std::move(field));
    return fields;
}

bool is_occupied(const cov::MolecularOrbital& orbital) {
    return orbital.occupation > 0.25f;
}

bool is_somo(const cov::MolecularOrbital& orbital) {
    return orbital.occupation > 0.25f && orbital.occupation < 1.75f;
}

struct SideCoverage {
    std::size_t total = 0u;
    std::size_t matched = 0u;
    std::size_t occupied = 0u;
    std::size_t occupied_matched = 0u;
    std::size_t somo = 0u;
    std::size_t somo_matched = 0u;
};

SideCoverage coverage(const cov::Wavefunction& wavefunction,
                      const std::vector<bool>& matched) {
    SideCoverage result;
    result.total = wavefunction.orbitals.size();
    for (std::size_t index = 0u; index < wavefunction.orbitals.size(); ++index) {
        if (matched[index]) ++result.matched;
        if (is_occupied(wavefunction.orbitals[index])) {
            ++result.occupied;
            if (matched[index]) ++result.occupied_matched;
        }
        if (is_somo(wavefunction.orbitals[index])) {
            ++result.somo;
            if (matched[index]) ++result.somo_matched;
        }
    }
    return result;
}

double fraction(const std::size_t numerator, const std::size_t denominator) {
    return denominator == 0u ? 1.0
                             : static_cast<double>(numerator) /
                                   static_cast<double>(denominator);
}

std::string unmatched_occupied_summary(
    const cov::Wavefunction& wavefunction,
    const std::vector<std::vector<std::size_t>>& groups) {
    std::ostringstream output;
    bool first = true;
    for (const auto& group : groups) {
        for (const auto index : group) {
            if (index >= wavefunction.orbitals.size() ||
                !is_occupied(wavefunction.orbitals[index])) continue;
            const auto& orbital = wavefunction.orbitals[index];
            if (!first) output << ',';
            first = false;
            output << index << '@' << orbital.energy_hartree;
            if (!orbital.symmetry.empty()) output << ':' << orbital.symmetry;
        }
    }
    return output.str();
}

std::string validate_transition(const cov::Wavefunction& from,
                                const cov::Wavefunction& to,
                                const cov::OrbitalTrackingResult& tracking,
                                SideCoverage& from_coverage,
                                SideCoverage& to_coverage) {
    if (!tracking.atom_mapping_compatible) return "incompatible atom mapping";
    if (tracking.composite_optimisation_truncated) {
        return "composite optimiser used conservative fallback";
    }
    std::vector<bool> matched_from(from.orbitals.size(), false);
    std::vector<bool> matched_to(to.orbitals.size(), false);
    std::vector<bool> accounted_from(from.orbitals.size(), false);
    std::vector<bool> accounted_to(to.orbitals.size(), false);
    for (const auto& match : tracking.matches) {
        if (match.from_members.size() != match.to_members.size() ||
            match.from_members.empty()) {
            return "dimension-changing subspace match";
        }
        if (!std::isfinite(match.similarity) || !std::isfinite(match.score)) {
            return "non-finite match metric";
        }
        for (const auto index : match.from_members) {
            if (index >= from.orbitals.size() || accounted_from[index]) {
                return "duplicate/invalid source member";
            }
            accounted_from[index] = true;
            matched_from[index] = true;
        }
        for (const auto index : match.to_members) {
            if (index >= to.orbitals.size() || accounted_to[index]) {
                return "duplicate/invalid target member";
            }
            accounted_to[index] = true;
            matched_to[index] = true;
        }
    }
    for (const auto& group : tracking.unmatched_from) {
        for (const auto index : group) {
            if (index >= from.orbitals.size() || accounted_from[index]) {
                return "duplicate/invalid unmatched source member";
            }
            accounted_from[index] = true;
        }
    }
    for (const auto& group : tracking.unmatched_to) {
        for (const auto index : group) {
            if (index >= to.orbitals.size() || accounted_to[index]) {
                return "duplicate/invalid unmatched target member";
            }
            accounted_to[index] = true;
        }
    }
    if (std::find(accounted_from.begin(), accounted_from.end(), false) !=
            accounted_from.end() ||
        std::find(accounted_to.begin(), accounted_to.end(), false) !=
            accounted_to.end()) {
        return "silent orbital loss outside matched/unmatched partitions";
    }

    from_coverage = coverage(from, matched_from);
    to_coverage = coverage(to, matched_to);
    const double total_coverage = std::min(
        fraction(from_coverage.matched, from_coverage.total),
        fraction(to_coverage.matched, to_coverage.total));
    const double occupied_coverage = std::min(
        fraction(from_coverage.occupied_matched, from_coverage.occupied),
        fraction(to_coverage.occupied_matched, to_coverage.occupied));
    if (total_coverage < 0.60) return "canonical-orbital coverage below 60%";
    // Occupied subspaces are more chemically important than the diffuse
    // virtual tail, but at an exact high-symmetry reaction frame canonical
    // orbitals can reorganise beyond what an atom/angular descriptor can
    // identify honestly. Require strong coverage and keep every remainder in
    // the explicit unmatched partition instead of inventing a false match.
    if (occupied_coverage < 0.80) return "occupied-orbital coverage below 80%";
    if (from_coverage.somo_matched != from_coverage.somo ||
        to_coverage.somo_matched != to_coverage.somo) {
        return "SOMO identity silently lost";
    }
    return {};
}

} // namespace

int main(int argc, char** argv) {
    std::cout << std::unitbuf;
    if (argc != 2) {
        std::cout << "usage: cov_scan_manifest_regression <manifest.tsv>\n";
        return 2;
    }

    const std::filesystem::path manifest = std::filesystem::absolute(argv[1]);
    std::ifstream input(manifest, std::ios::binary);
    if (!input) {
        std::cout << "FAIL\tmanifest\tcannot open " << manifest.string() << '\n';
        return 2;
    }

    std::string line;
    if (!std::getline(input, line)) {
        std::cout << "FAIL\tmanifest\tempty manifest\n";
        return 2;
    }
    const auto header = split_tsv(line);
    const auto find_column = [&](const std::string& name) {
        const auto iterator = std::find(header.begin(), header.end(), name);
        return iterator == header.end()
            ? header.size()
            : static_cast<std::size_t>(iterator - header.begin());
    };
    const std::size_t series_column = find_column("series_id");
    const std::size_t frame_column = find_column("frame");
    const std::size_t path_column = find_column("fch_path");
    if (series_column == header.size() || frame_column == header.size() ||
        path_column == header.size()) {
        std::cout << "FAIL\tmanifest\tmissing series_id/frame/fch_path column\n";
        return 2;
    }

    std::map<std::string, std::vector<ManifestFrame>> series;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tsv(line);
        if (fields.size() != header.size()) {
            std::cout << "FAIL\tmanifest\tmalformed TSV row\n";
            return 2;
        }
        ManifestFrame row;
        try {
            row.frame = std::stoi(fields[frame_column]);
        } catch (...) {
            std::cout << "FAIL\tmanifest\tinvalid frame index\n";
            return 2;
        }
        row.fch_path = fields[path_column];
        series[fields[series_column]].push_back(std::move(row));
    }

    std::size_t frames = 0u;
    std::size_t transitions = 0u;
    std::size_t passed = 0u;
    std::size_t failed = 0u;
    for (auto& [series_id, rows] : series) {
        std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
            return left.frame < right.frame;
        });
        frames += rows.size();
        if (rows.size() < 2u) {
            ++failed;
            std::cout << "FAIL\t" << series_id << "\tseries has fewer than two frames\n";
            continue;
        }
        cov::WavefunctionParseOptions options;
        options.max_atoms = 512u;
        options.require_orbitals = true;
        options.keep_density = true;
        options.reconstruct_density_if_missing = true;
        try {
            cov::Wavefunction previous = cov::parse_wavefunction(
                rows.front().fch_path, options);
            for (std::size_t index = 1u; index < rows.size(); ++index) {
                ++transitions;
                cov::Wavefunction current = cov::parse_wavefunction(
                    rows[index].fch_path, options);
                const auto tracking_start = std::chrono::steady_clock::now();
                const auto tracking = cov::track_orbital_subspaces(previous, current);
                const double tracking_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - tracking_start).count();
                SideCoverage from_coverage;
                SideCoverage to_coverage;
                const std::string failure = validate_transition(
                    previous, current, tracking, from_coverage, to_coverage);
                if (failure.empty()) {
                    ++passed;
                    std::cout << "PASS\t" << series_id << '\t'
                              << rows[index - 1u].frame << "->" << rows[index].frame
                              << "\ttotal=" << from_coverage.matched << '/'
                              << from_coverage.total << ":" << to_coverage.matched
                              << '/' << to_coverage.total << "\toccupied="
                              << from_coverage.occupied_matched << '/'
                              << from_coverage.occupied << ':'
                              << to_coverage.occupied_matched << '/'
                              << to_coverage.occupied << "\ttracking_s="
                              << tracking_seconds << '\n';
                } else {
                    ++failed;
                    std::cout << "FAIL\t" << series_id << '\t'
                              << rows[index - 1u].frame << "->" << rows[index].frame
                              << '\t' << failure << "\ttotal="
                              << from_coverage.matched << '/' << from_coverage.total
                              << ':' << to_coverage.matched << '/'
                              << to_coverage.total << "\toccupied="
                              << from_coverage.occupied_matched << '/'
                              << from_coverage.occupied << ':'
                              << to_coverage.occupied_matched << '/'
                              << to_coverage.occupied << "\tunmatched_occ="
                              << unmatched_occupied_summary(
                                     previous, tracking.unmatched_from)
                              << ':' << unmatched_occupied_summary(
                                     current, tracking.unmatched_to)
                              << '\n';
                }
                previous = std::move(current);
            }
        } catch (const std::exception& error) {
            ++failed;
            std::cout << "FAIL\t" << series_id << "\texception: "
                      << error.what() << '\n';
        }
    }
    std::cout << "SUMMARY\tseries=" << series.size() << "\tframes=" << frames
              << "\ttransitions=" << transitions << "\tpassed=" << passed
              << "\tfailed=" << failed << '\n';
    return failed == 0u ? EXIT_SUCCESS : EXIT_FAILURE;
}
