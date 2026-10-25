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

#include "cov/gaussian_log.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

std::string trim(std::string value) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool starts_with(const std::string& value, const char* prefix) {
    const std::string p(prefix);
    return value.size() >= p.size() && value.compare(0, p.size(), p) == 0;
}

std::vector<std::string> parenthesized_labels(const std::string& line) {
    std::vector<std::string> labels;
    std::size_t cursor = 0;
    while (true) {
        const auto open = line.find('(', cursor);
        if (open == std::string::npos) break;
        const auto close = line.find(')', open + 1u);
        if (close == std::string::npos) break;
        std::string label = trim(line.substr(open + 1u, close - open - 1u));
        if (!label.empty()) labels.push_back(std::move(label));
        cursor = close + 1u;
    }
    return labels;
}

std::string point_group_field(const std::string& line, const char* prefix,
                             std::string& field_text) {
    const auto pos = line.find(prefix);
    if (pos == std::string::npos) return {};
    const auto begin = pos + std::string(prefix).size();
    const auto end = line.find("NOp", begin);
    const std::string field = trim(line.substr(begin, end == std::string::npos
                                                        ? end : end - begin));
    field_text=field;
    // Validate Schoenflies syntax, not membership of the smaller display
    // catalogue. Gaussian's linear aliases D*H/C*V are valid producer values.
    // In particular, the next column name is never the missing field's value.
    static const std::regex syntax(
        R"((?:C(?:[1-9][0-9]*[HV]?|S|I|INFV|\*V)|D(?:[1-9][0-9]*[HD]?|INFH|\*H)|S[1-9][0-9]*|T[DH]?|OH?|IH?|KH?))",
        std::regex::icase);
    return std::regex_match(field, syntax) ? field : std::string{};
}

struct SymmetryBlock {
    std::vector<std::string> alpha;
    std::vector<std::string> beta;
    bool explicit_spin = false;
    std::size_t job_segment = 0;
    OrbitalSymmetrySourceRecord source;
};

struct TextRecord {
    std::string value;
    std::size_t job_segment = 0;
};

struct SpinSquaredRecord {
    std::size_t line = 0;
    std::size_t job_segment = 0;
    double before = 0.0;
    double after = 0.0;
};

struct StatusRecord {
    std::size_t line = 0;
    std::size_t job_segment = 0;
    bool positive = false;
    std::string detail;
};

std::optional<double> first_real_after(const std::string& line,
                                       const std::size_t offset) {
    if (offset >= line.size()) return std::nullopt;
    static const std::regex number(
        R"([+-]?(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[DdEe][+-]?[0-9]+)?)");
    std::smatch match;
    const std::string tail = line.substr(offset);
    if (!std::regex_search(tail, match, number)) return std::nullopt;
    std::string token = match.str();
    for (char& c : token) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    try {
        const double value = std::stod(token);
        return std::isfinite(value) ? std::optional<double>(value) : std::nullopt;
    } catch (...) {
        return std::nullopt;
    }
}

std::optional<SpinSquaredRecord> spin_squared_record(
    const std::string& line,
    const std::size_t line_index) {
    const std::string folded = lower(line);
    const bool spin_marker =
        folded.find("s**2") != std::string::npos ||
        folded.find("s^2") != std::string::npos ||
        folded.find("<s2>") != std::string::npos;
    if (!spin_marker) return std::nullopt;

    const auto before_phrase = folded.find("before annihilation");
    if (before_phrase == std::string::npos) return std::nullopt;
    const auto after_phrase = folded.find("after", before_phrase + 19u);
    if (after_phrase == std::string::npos) return std::nullopt;

    const auto before = first_real_after(line, before_phrase + 19u);
    const auto after = first_real_after(line, after_phrase + 5u);
    if (!before || !after || *before < 0.0 || *after < 0.0) {
        return std::nullopt;
    }
    return SpinSquaredRecord{line_index, 0u, *before, *after};
}

} // namespace

GaussianLogEnrichmentResult enrich_from_gaussian_log(
    Wavefunction& wavefunction,
    const std::filesystem::path& path) {
    GaussianLogEnrichmentResult result;
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        result.warning = "Unable to open Gaussian enrichment file: " + path.string();
        return result;
    }
    result.opened = true;
    wavefunction.point_group_source_records.clear();
    wavefunction.orbital_symmetry_source_records.clear();

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) lines.push_back(line);

    std::optional<TextRecord> point_group_detected;
    std::optional<TextRecord> point_group_used;
    SymmetryBlock last_block;
    std::optional<SpinSquaredRecord> last_spin_squared;
    std::optional<StatusRecord> last_scf_status;
    std::optional<StatusRecord> last_stability_status;
    std::size_t job_segment = 0u;

    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string t = trim(lines[i]);
        const std::string folded = lower(t);
        if (starts_with(t, "Entering Link 1")) ++job_segment;
        if (t.find("Full point group") != std::string::npos) {
            std::string field_text;
            const auto value = point_group_field(t, "Full point group",field_text);
            // An invalid later geometry record invalidates an earlier value
            // in the same job too; it must not silently resurrect that value.
            point_group_detected=TextRecord{value,job_segment};
            wavefunction.point_group_source_records.push_back(
                {"detected",lines[i],value,i+1u,job_segment,!value.empty(),field_text,
                 field_text.empty()?"missing-input":value.empty()?"invalid-input":"available"});
        }
        if (t.find("Largest Abelian subgroup") != std::string::npos) {
            std::string field_text;
            const auto value = point_group_field(t, "Largest Abelian subgroup",field_text);
            point_group_used=TextRecord{value,job_segment};
            wavefunction.point_group_source_records.push_back(
                {"used",lines[i],value,i+1u,job_segment,!value.empty(),field_text,
                 field_text.empty()?"missing-input":value.empty()?"invalid-input":"available"});
        }

        if (auto spin = spin_squared_record(t, i)) {
            spin->job_segment = job_segment;
            last_spin_squared = spin;
        }

        // Only explicit Gaussian outcomes are accepted.  A route containing
        // SCF=... or Stable=... is a request, not evidence of an outcome.
        if (starts_with(t, "SCF Done:")) {
            last_scf_status = StatusRecord{i, job_segment, true, t};
        } else if (folded.find("convergence failure -- run terminated") !=
                       std::string::npos ||
                   folded.find("scf has not converged") != std::string::npos ||
                   folded.find("scf failed to converge") != std::string::npos) {
            last_scf_status = StatusRecord{i, job_segment, false, t};
        }

        if (folded.find(
                "the wavefunction is stable under the perturbations considered") !=
            std::string::npos) {
            last_stability_status = StatusRecord{i, job_segment, true, t};
        } else if (folded.find("the wavefunction is unstable") !=
                       std::string::npos ||
                   (folded.find("the wavefunction has") != std::string::npos &&
                    folded.find("instability") != std::string::npos) ||
                   folded.find("wavefunction stability test failed") !=
                       std::string::npos) {
            last_stability_status = StatusRecord{i, job_segment, false, t};
        }

        // Deliberately exclude "Initial guess orbital symmetries:". In an
        // optimization or converged SCF output there may be several final
        // blocks; the last complete producer block wins.
        if (t != "Orbital symmetries:") continue;

        SymmetryBlock block;
        block.job_segment=job_segment;
        const auto source_path_utf8=path.u8string();
        block.source.source_path.assign(source_path_utf8.begin(),source_path_utf8.end());
        block.source.line_begin=i+1u;
        block.source.line_end=i+1u;
        block.source.job_segment=job_segment;
        if (point_group_detected && point_group_detected->job_segment==job_segment)
            block.source.detected_group_context=point_group_detected->value;
        if (point_group_used && point_group_used->job_segment==job_segment)
            block.source.abelian_group_context=point_group_used->value;
        Spin current_spin = Spin::Alpha;
        for (std::size_t j = i + 1u; j < lines.size() && j < i + 80u; ++j) {
            const std::string s = trim(lines[j]);
            if (starts_with(s,"Entering Link 1")) break;
            if (s == "Alpha Orbitals:" || s == "Alpha orbitals:") {
                current_spin = Spin::Alpha;
                block.explicit_spin = true;
                continue;
            }
            if (s == "Beta Orbitals:" || s == "Beta orbitals:") {
                current_spin = Spin::Beta;
                block.explicit_spin = true;
                continue;
            }
            if (starts_with(s, "The electronic state") ||
                starts_with(s, "Alpha  occ. eigenvalues") ||
                starts_with(s, "Beta  occ. eigenvalues") ||
                starts_with(s, "SCF Done:") ||
                s.find("Leave Link") != std::string::npos) {
                break;
            }

            const auto labels = parenthesized_labels(s);
            if (labels.empty()) continue;
            block.source.line_end=j+1u;
            auto& target = current_spin == Spin::Alpha ? block.alpha : block.beta;
            target.insert(target.end(), labels.begin(), labels.end());
        }

        if (!block.alpha.empty() || !block.beta.empty()) {
            last_block = std::move(block);
        }
    }

    const std::size_t final_job_segment=
        last_scf_status?last_scf_status->job_segment:job_segment;

    const bool detected_is_current=point_group_detected &&
        point_group_detected->job_segment==final_job_segment &&
        !point_group_detected->value.empty();
    const bool used_is_current=point_group_used &&
        point_group_used->job_segment==final_job_segment &&
        !point_group_used->value.empty();
    if (detected_is_current || used_is_current) {
        wavefunction.point_group_detected = detected_is_current
            ?point_group_detected->value:std::string{};
        wavefunction.point_group_used = used_is_current
            ?point_group_used->value:std::string{};
        wavefunction.point_group_provenance = DataProvenance::Producer;
        wavefunction.point_group_detected_provenance = detected_is_current
            ?DataProvenance::Producer:DataProvenance::Unavailable;
        wavefunction.point_group_used_provenance = used_is_current
            ?DataProvenance::Producer:DataProvenance::Unavailable;
        result.point_group_applied = true;
    }

    std::vector<std::size_t> alpha_indices;
    std::vector<std::size_t> beta_indices;
    for (std::size_t i = 0; i < wavefunction.orbitals.size(); ++i) {
        if (wavefunction.orbitals[i].spin == Spin::Beta) beta_indices.push_back(i);
        else alpha_indices.push_back(i);
    }

    const bool symmetry_segment_ok=last_block.job_segment==final_job_segment;
    bool symmetry_dimensions_ok = true;
    if (!last_block.beta.empty() || last_block.explicit_spin) {
        symmetry_dimensions_ok =
            last_block.alpha.size() == alpha_indices.size() &&
            last_block.beta.size() == beta_indices.size();
    } else if (!beta_indices.empty()) {
        // Never mirror one producer-reported spin block onto another.
        symmetry_dimensions_ok = false;
    } else {
        symmetry_dimensions_ok = last_block.alpha.size() == alpha_indices.size();
    }

    if (symmetry_segment_ok && symmetry_dimensions_ok &&
        !last_block.alpha.empty()) {
        auto source=last_block.source;
        source.orbital_indices=alpha_indices;
        source.orbital_indices.insert(source.orbital_indices.end(),beta_indices.begin(),beta_indices.end());
        source.labels=last_block.alpha;
        source.labels.insert(source.labels.end(),last_block.beta.begin(),last_block.beta.end());
        wavefunction.orbital_symmetry_source_records.push_back(std::move(source));
        for (std::size_t i = 0; i < alpha_indices.size(); ++i) {
            auto& mo = wavefunction.orbitals[alpha_indices[i]];
            mo.symmetry = last_block.alpha[i];
            mo.symmetry_provenance = DataProvenance::Producer;
            ++result.symmetry_labels_applied;
        }
        for (std::size_t i = 0; i < beta_indices.size(); ++i) {
            auto& mo = wavefunction.orbitals[beta_indices[i]];
            mo.symmetry = last_block.beta[i];
            mo.symmetry_provenance = DataProvenance::Producer;
            ++result.symmetry_labels_applied;
        }
    } else if (symmetry_segment_ok &&
               (!last_block.alpha.empty() || !last_block.beta.empty())) {
        result.warning =
            "Gaussian orbital-symmetry count does not match the loaded MO blocks; "
            "symmetry enrichment was left unapplied";
    }

    if (last_scf_status) {
        wavefunction.scf_convergence = last_scf_status->positive
            ? ScfConvergenceStatus::Converged
            : ScfConvergenceStatus::Failed;
        wavefunction.scf_convergence_provenance = DataProvenance::Producer;
        result.scf_diagnostic_applied = true;
    }

    const bool stability_is_current = last_stability_status &&
        (!last_scf_status ||
         (last_scf_status->positive &&
          last_stability_status->job_segment == last_scf_status->job_segment &&
          last_stability_status->line >= last_scf_status->line));
    if (stability_is_current) {
        wavefunction.stability = last_stability_status->positive
            ? WavefunctionStabilityStatus::Stable
            : WavefunctionStabilityStatus::Unstable;
        wavefunction.stability_provenance = DataProvenance::Producer;
        wavefunction.stability_detail = last_stability_status->detail;
        result.stability_diagnostic_applied = true;
    }

    // In an optimization a LOG can contain many S**2 reports.  Bind the final
    // one to the final successful SCF statement and reject stale values left
    // behind before a later SCF or in an earlier Link1 job. Standalone
    // diagnostic excerpts without an SCF status are still accepted because
    // the before/after pair is itself a complete producer statement.
    bool spin_is_current = last_spin_squared.has_value();
    if (last_spin_squared && last_scf_status) {
        const std::size_t separation = last_spin_squared->line >=
                                                last_scf_status->line
                                           ? last_spin_squared->line -
                                                 last_scf_status->line
                                           : std::numeric_limits<std::size_t>::max();
        spin_is_current = last_scf_status->positive &&
            last_spin_squared->job_segment == last_scf_status->job_segment &&
            separation <= 50u;
    }
    if (last_spin_squared && spin_is_current) {
        wavefunction.spin_squared_before_annihilation =
            last_spin_squared->before;
        wavefunction.spin_squared_after_annihilation =
            last_spin_squared->after;
        wavefunction.spin_squared_provenance = DataProvenance::Producer;
        result.spin_squared_applied = true;
    }

    if (result.applied() || !wavefunction.point_group_source_records.empty())
        wavefunction.enrichment_source = path.string();
    return result;
}

std::filesystem::path find_sibling_gaussian_log(
    const std::filesystem::path& wavefunction_path) {
    std::string filename = wavefunction_path.filename().string();
    const std::string folded = lower(filename);
    std::string root;
    const auto molden_pos = folded.find(".molden");
    if (molden_pos != std::string::npos) {
        root = filename.substr(0, molden_pos);
    } else {
        root = wavefunction_path.stem().string();
    }
    if (root.empty()) return {};

    const auto parent = wavefunction_path.parent_path();
    std::error_code ec;
    for (const char* extension : {".log", ".out"}) {
        const auto candidate = parent / (root + extension);
        if (std::filesystem::exists(candidate, ec) && !ec) return candidate;
        ec.clear();
    }

    // Gaussian projects commonly keep equivalent directory trees below
    // dedicated .chk and .log roots. Preserve the relative subdirectory and
    // probe that producer output only after the true sibling lookup above.
    // This is a filesystem convention lookup, not chemistry inference.
    std::filesystem::path relative;
    for (auto cursor = parent; !cursor.empty();) {
        if (lower(cursor.filename().string()) == ".chk") {
            const auto log_parent = cursor.parent_path() / ".log" / relative;
            for (const char* extension : {".log", ".out"}) {
                const auto candidate = log_parent / (root + extension);
                if (std::filesystem::exists(candidate, ec) && !ec) {
                    return candidate;
                }
                ec.clear();
            }
            break;
        }
        const auto next = cursor.parent_path();
        if (next == cursor) break;
        relative = cursor.filename() / relative;
        cursor = next;
    }
    return {};
}

} // namespace cov
