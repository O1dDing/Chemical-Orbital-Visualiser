#include "cov/coordination_geometry.hpp"
#include "cov/interaction_graph.hpp"
#include "cov/local_geometry.hpp"
#include "cov/model.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <locale>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using Row = std::map<std::string, std::string>;

struct AtomPair {
    std::uint32_t first = 0;
    std::uint32_t second = 0;
};

struct MulticentreExpectation {
    cov::MulticentreKind kind = cov::MulticentreKind::Unclassified;
    std::vector<std::uint32_t> atoms;
};

std::string trim(std::string value) {
    const auto not_space = [](const unsigned char character) {
        return !std::isspace(character);
    };
    const auto first = std::find_if(value.begin(), value.end(), not_space);
    const auto last = std::find_if(value.rbegin(), value.rend(), not_space).base();
    if (first >= last) return {};
    return std::string(first, last);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](const unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    return value;
}

std::string semantic_key(const std::string_view value) {
    std::string key;
    key.reserve(value.size());
    for (const unsigned char character : value) {
        if (std::isalnum(character)) {
            key.push_back(static_cast<char>(std::tolower(character)));
        }
    }
    return key;
}

std::string geometry_key(std::string value) {
    value = lower(trim(std::move(value)));
    for (const std::string_view modifier : {
             "constrained", "constraint", "idealized", "idealised", "ideal",
             "reference"}) {
        std::size_t position = 0u;
        while ((position = value.find(modifier, position)) != std::string::npos) {
            value.erase(position, modifier.size());
        }
    }
    return semantic_key(value);
}

bool unavailable(const std::string_view raw) {
    const std::string value = semantic_key(trim(std::string(raw)));
    return value.empty() || value == "na" || value == "null" ||
           value == "none" || value == "notapplicable";
}

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (character == '"') {
            if (quoted && index + 1u < line.size() && line[index + 1u] == '"') {
                field.push_back('"');
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (character == '\t' && !quoted) {
            fields.push_back(trim(std::move(field)));
            field.clear();
        } else if (character != '\r') {
            field.push_back(character);
        }
    }
    fields.push_back(trim(std::move(field)));
    return fields;
}

std::vector<std::string> split(const std::string_view value, const char separator) {
    std::vector<std::string> result;
    std::size_t begin = 0u;
    while (begin <= value.size()) {
        const std::size_t end = value.find(separator, begin);
        result.push_back(trim(std::string(value.substr(
            begin, end == std::string_view::npos ? value.size() - begin : end - begin))));
        if (end == std::string_view::npos) break;
        begin = end + 1u;
    }
    return result;
}

std::optional<long long> parse_integer(const std::string_view raw) {
    std::string value = trim(std::string(raw));
    if (value.empty()) return std::nullopt;
    bool positive_sign = value.front() == '+';
    if (positive_sign) value.erase(value.begin());
    if (value.empty()) return std::nullopt;
    long long parsed = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (result.ec != std::errc{} || result.ptr != value.data() + value.size()) {
        return std::nullopt;
    }
    return parsed;
}

std::optional<double> parse_real(const std::string_view raw) {
    std::string value=trim(std::string(raw));
    if (value.empty()) return std::nullopt;
    double parsed=0.0;
    // Apple's libc++ does not provide floating-point from_chars in every SDK.
    // The classic locale preserves the manifest's decimal-dot syntax.
    if (value.front()=='+') return std::nullopt;
    std::istringstream input(value);
    input.imbue(std::locale::classic());
    if (!(input >> std::noskipws >> parsed) || input.peek()!=std::char_traits<char>::eof() ||
        !std::isfinite(parsed)) {
        return std::nullopt;
    }
    return parsed;
}

std::optional<bool> parse_boolean(const std::string_view raw) {
    const std::string value = semantic_key(trim(std::string(raw)));
    if (value == "1" || value == "true" || value == "yes" || value == "y" ||
        value == "success" || value == "passed" || value == "complete" ||
        value == "completed" || value == "normal") {
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "n" ||
        value == "failed" || value == "failure" || value == "error" ||
        value == "abnormal") {
        return false;
    }
    return std::nullopt;
}

bool explicitly_not_ready(const std::string_view raw) {
    const std::string value = semantic_key(trim(std::string(raw)));
    static const std::set<std::string> states{
        "pending", "planned", "queued", "running", "notrun", "generated",
        "failed", "failure", "error", "aborted", "timeout", "missing",
        "skipped", "unrunnable",
    };
    return states.contains(value);
}

const std::string& cell(const Row& row, const std::string& name) {
    static const std::string empty;
    const auto found = row.find(name);
    return found == row.end() ? empty : found->second;
}

std::filesystem::path resolve_path(const std::filesystem::path& manifest,
                                   const std::string_view raw) {
    std::filesystem::path path(trim(std::string(raw)));
    if (path.empty() || path.is_absolute()) return path;
    return (manifest.parent_path() / path).lexically_normal();
}

std::optional<std::vector<AtomPair>> parse_pairs(const std::string& raw,
                                                 std::string& error) {
    std::vector<AtomPair> result;
    if (unavailable(raw)) return result;
    std::set<std::pair<std::uint32_t, std::uint32_t>> unique;
    for (std::string token : split(raw, ';')) {
        if (token.empty()) continue;
        const std::size_t annotation = token.find('@');
        if (annotation != std::string::npos) token.erase(annotation);
        const std::size_t dash = token.find('-');
        if (dash == std::string::npos || token.find('-', dash + 1u) != std::string::npos) {
            error = "invalid atom-pair token '" + token + "' (expected 1-2)";
            return std::nullopt;
        }
        const auto first = parse_integer(token.substr(0u, dash));
        const auto second = parse_integer(token.substr(dash + 1u));
        if (!first || !second || *first <= 0 || *second <= 0 || *first == *second) {
            error = "invalid one-based atom pair '" + token + "'";
            return std::nullopt;
        }
        const auto a = static_cast<std::uint32_t>(*first - 1);
        const auto b = static_cast<std::uint32_t>(*second - 1);
        const auto ordered = std::minmax(a, b);
        if (unique.emplace(ordered.first, ordered.second).second) {
            result.push_back({ordered.first, ordered.second});
        }
    }
    return result;
}

std::optional<std::vector<std::uint32_t>> parse_atom_set(const std::string& raw,
                                                         std::string& error) {
    std::vector<std::uint32_t> atoms;
    for (const auto& token : split(raw, '-')) {
        const auto index = parse_integer(token);
        if (!index || *index <= 0) {
            error = "invalid one-based multicentre atom '" + token + "'";
            return std::nullopt;
        }
        atoms.push_back(static_cast<std::uint32_t>(*index - 1));
    }
    std::sort(atoms.begin(), atoms.end());
    if (atoms.size() < 3u || std::adjacent_find(atoms.begin(), atoms.end()) != atoms.end()) {
        error = "multicentre atom set must contain at least three distinct atoms";
        return std::nullopt;
    }
    return atoms;
}

std::optional<std::vector<std::uint32_t>> parse_pi_atom_set(
    const std::string& raw,std::string& error) {
    std::vector<std::uint32_t> atoms;
    for (const auto& token:split(raw,'-')) {
        const auto index=parse_integer(token);
        if (!index || *index<=0) {
            error="invalid one-based pi atom '"+token+"'";
            return std::nullopt;
        }
        atoms.push_back(static_cast<std::uint32_t>(*index-1));
    }
    std::sort(atoms.begin(),atoms.end());
    if (atoms.size()<2u ||
        std::adjacent_find(atoms.begin(),atoms.end())!=atoms.end()) {
        error="pi atom set must contain at least two distinct atoms";
        return std::nullopt;
    }
    return atoms;
}

std::optional<std::vector<MulticentreExpectation>> parse_multicentre(
    const std::string& raw, bool& require_none, std::string& error) {
    require_none = semantic_key(raw) == "none";
    std::vector<MulticentreExpectation> result;
    if (require_none || unavailable(raw)) return result;
    for (const std::string& token : split(raw, ';')) {
        if (token.empty()) continue;
        const std::size_t colon = token.find(':');
        const std::string kind_text = semantic_key(token.substr(0u, colon));
        cov::MulticentreKind kind = cov::MulticentreKind::Unclassified;
        if (kind_text == "3c2e" || kind_text == "threecentretwoelectron") {
            kind = cov::MulticentreKind::ThreeCentreTwoElectron;
        } else if (kind_text == "3c4e" || kind_text == "threecentrefourelectron") {
            kind = cov::MulticentreKind::ThreeCentreFourElectron;
        } else {
            error = "unknown multicentre kind in '" + token + "'";
            return std::nullopt;
        }
        MulticentreExpectation expectation;
        expectation.kind = kind;
        if (colon != std::string::npos) {
            auto atoms = parse_atom_set(token.substr(colon + 1u), error);
            if (!atoms) return std::nullopt;
            expectation.atoms = std::move(*atoms);
        }
        result.push_back(std::move(expectation));
    }
    return result;
}

bool strong_pair(const cov::InteractionGraph& graph, const AtomPair pair) {
    return std::any_of(graph.edges.begin(), graph.edges.end(), [&](const auto& edge) {
        return cov::interaction_merges_fragments(edge.kind) &&
               ((edge.atom_a == pair.first && edge.atom_b == pair.second) ||
                (edge.atom_a == pair.second && edge.atom_b == pair.first));
    });
}

std::string one_based_pair(const AtomPair pair) {
    return std::to_string(pair.first + 1u) + "-" +
           std::to_string(pair.second + 1u);
}

bool geometry_matches(const cov::LocalGeometryEnvironment& environment,
                      const std::vector<std::string>& expected_keys) {
    const auto* descriptor = cov::coordination_geometry_descriptor(environment.geometry_id);
    if (descriptor == nullptr) return false;
    const std::vector<std::string> actual{
        semantic_key(descriptor->machine_id), semantic_key(descriptor->name),
        semantic_key(descriptor->point_group),
    };
    return std::any_of(expected_keys.begin(), expected_keys.end(),
                       [&](const auto& expected) {
                           return std::find(actual.begin(), actual.end(), expected) != actual.end();
                       });
}

std::vector<std::string> geometry_expectations(const std::string& raw) {
    std::vector<std::string> result;
    std::string token;
    for (const char character : raw) {
        if (character == '|' || character == '/') {
            const std::string key = geometry_key(token);
            if (!key.empty()) result.push_back(key);
            token.clear();
        } else {
            token.push_back(character);
        }
    }
    const std::string key = geometry_key(token);
    if (!key.empty()) result.push_back(key);
    return result;
}

void append_failure(std::vector<std::string>& failures, std::string failure) {
    failures.push_back(std::move(failure));
}

std::string join_failures(const std::vector<std::string>& failures) {
    std::ostringstream output;
    for (std::size_t index = 0u; index < failures.size(); ++index) {
        if (index != 0u) output << " | ";
        output << failures[index];
    }
    return output.str();
}

void validate_shared_multicentre_sources(
    const cov::Wavefunction& wavefunction,
    std::vector<std::string>& failures) {
    std::map<std::vector<std::uint32_t>,
             std::vector<const cov::MulticentreAssignment*>> repeated_orbitals;
    std::map<std::string,std::vector<const cov::MulticentreAssignment*>> groups;
    for (const auto& assignment:wavefunction.multicentre_assignments) {
        if (!assignment.orbitals.empty()) {
            repeated_orbitals[assignment.orbitals].push_back(&assignment);
        }
        if (!assignment.source_subspace_id.empty()) {
            groups[assignment.source_subspace_id].push_back(&assignment);
        }
    }
    for (const auto& [orbitals,channels]:repeated_orbitals) {
        if (channels.size()<2u) continue;
        const std::string id=channels.front()->source_subspace_id;
        if (id.empty() || std::any_of(
                channels.begin(),channels.end(),[&](const auto* channel) {
                    return channel->source_subspace_id!=id;
                })) {
            append_failure(failures,
                "repeated canonical multicentre source lacks one shared identity");
        }
    }
    for (const auto& [id,channels]:groups) {
        if (channels.size()<2u) {
            append_failure(failures,"shared multicentre source has one channel: "+id);
            continue;
        }
        const auto source_orbitals=channels.front()->orbitals;
        const double source_electrons=
            channels.front()->source_subspace_electron_count;
        std::vector<std::uint32_t> source_atoms;
        for (const auto* channel:channels) {
            source_atoms.insert(source_atoms.end(),channel->atoms.begin(),channel->atoms.end());
        }
        std::sort(source_atoms.begin(),source_atoms.end());
        source_atoms.erase(std::unique(source_atoms.begin(),source_atoms.end()),
                           source_atoms.end());
        double fraction_sum=0.0;
        double electron_sum=0.0;
        for (const auto* channel:channels) {
            if (channel->orbitals!=source_orbitals ||
                std::abs(channel->source_subspace_electron_count-
                         source_electrons)>1.0e-8 ||
                channel->source_subspace_fraction<=0.0 ||
                channel->source_subspace_fraction>1.0 ||
                std::abs(channel->electron_count-
                         source_electrons*channel->source_subspace_fraction)>0.05) {
                append_failure(failures,
                    "inconsistent shared multicentre source partition: "+id);
                break;
            }
            fraction_sum+=channel->source_subspace_fraction;
            electron_sum+=channel->electron_count;
        }
        if (std::abs(fraction_sum-1.0)>1.0e-8 ||
            std::abs(electron_sum-source_electrons)>0.05) {
            append_failure(failures,
                "shared multicentre source double-counts electrons: "+id);
        }
        for (const auto orbital_index:source_orbitals) {
            if (orbital_index>=wavefunction.orbitals.size()) {
                append_failure(failures,
                    "shared multicentre source has invalid orbital index: "+id);
                continue;
            }
            const auto& chemistry=wavefunction.orbitals[orbital_index].chemistry;
            if (chemistry.multicentre_channel_count!=channels.size() ||
                chemistry.multicentre_source_subspace_id!=id ||
                std::abs(chemistry.multicentre_source_electron_count-
                         source_electrons)>0.05 ||
                chemistry.multicentre_participating_atom_indices!=source_atoms ||
                std::abs(chemistry.multicentre_participating_electrons-
                         source_electrons)>0.05 ||
                chemistry.multicentre_label.find("shared")==std::string::npos) {
                append_failure(failures,
                    "shared multicentre source is not preserved in orbital chemistry: "+id);
                break;
            }
        }
    }
}

std::vector<std::string> validate_semantics(const Row& row,
                                            const std::filesystem::path& manifest,
                                            const cov::Wavefunction& wavefunction) {
    std::vector<std::string> failures;

    validate_shared_multicentre_sources(wavefunction,failures);

    if (!unavailable(cell(row,"shared_multicentre_channels"))) {
        const auto expected_channels=parse_integer(
            cell(row,"shared_multicentre_channels"));
        const auto expected_centres=parse_integer(
            cell(row,"shared_multicentre_centres"));
        const auto expected_electrons=parse_real(
            cell(row,"shared_multicentre_electrons"));
        if (!expected_channels || *expected_channels<2 ||
            !expected_centres || *expected_centres<3 ||
            !expected_electrons || *expected_electrons<=0.0) {
            append_failure(failures,"invalid shared multicentre expectation");
        } else {
            std::map<std::string,
                     std::vector<const cov::MulticentreAssignment*>> sources;
            for (const auto& assignment:wavefunction.multicentre_assignments) {
                if (!assignment.source_subspace_id.empty()) {
                    sources[assignment.source_subspace_id].push_back(&assignment);
                }
            }
            const bool matched=std::any_of(
                sources.begin(),sources.end(),[&](const auto& source) {
                    const auto& channels=source.second;
                    if (channels.size()!=static_cast<std::size_t>(*expected_channels)) {
                        return false;
                    }
                    const auto& orbitals=channels.front()->orbitals;
                    std::vector<std::uint32_t> atoms;
                    for (const auto* channel:channels) {
                        if (channel->orbitals!=orbitals) return false;
                        atoms.insert(atoms.end(),channel->atoms.begin(),channel->atoms.end());
                    }
                    std::sort(atoms.begin(),atoms.end());
                    atoms.erase(std::unique(atoms.begin(),atoms.end()),atoms.end());
                    return atoms.size()==static_cast<std::size_t>(*expected_centres) &&
                           std::abs(channels.front()->source_subspace_electron_count-
                                    *expected_electrons)<=0.05;
                });
            if (!matched) {
                append_failure(failures,
                    "required shared multicentre source was not derived");
            }
        }
    }

    if (!unavailable(cell(row, "charge"))) {
        const auto expected = parse_integer(cell(row, "charge"));
        if (!expected) {
            append_failure(failures, "invalid charge field");
        } else if (wavefunction.charge_provenance == cov::DataProvenance::Unavailable) {
            append_failure(failures, "FCHK has no producer charge");
        } else if (wavefunction.charge != *expected) {
            append_failure(failures, "charge expected " + std::to_string(*expected) +
                                         ", got " + std::to_string(wavefunction.charge));
        }
    }
    if (!unavailable(cell(row, "multiplicity"))) {
        const auto expected = parse_integer(cell(row, "multiplicity"));
        if (!expected || *expected <= 0) {
            append_failure(failures, "invalid multiplicity field");
        } else if (wavefunction.multiplicity_provenance == cov::DataProvenance::Unavailable) {
            append_failure(failures, "FCHK has no producer multiplicity");
        } else if (wavefunction.multiplicity != static_cast<std::uint32_t>(*expected)) {
            append_failure(failures, "multiplicity expected " + std::to_string(*expected) +
                                         ", got " + std::to_string(wavefunction.multiplicity));
        }
    }

    if (!unavailable(cell(row,"pi_channels"))) {
        const auto expected=parse_integer(cell(row,"pi_channels"));
        std::size_t actual=0u;
        for (const auto& assignment:wavefunction.delocalised_pi_assignments) {
            actual+=assignment.orientation_channels.size();
        }
        if (!expected || *expected<0) {
            append_failure(failures,"invalid pi_channels field");
        } else if (actual!=static_cast<std::size_t>(*expected)) {
            append_failure(failures,"pi channels expected "+
                std::to_string(*expected)+", got "+std::to_string(actual));
        }
    }
    if (!unavailable(cell(row,"pi_cyclic"))) {
        const auto expected=parse_boolean(cell(row,"pi_cyclic"));
        const bool actual=std::any_of(
            wavefunction.delocalised_pi_assignments.begin(),
            wavefunction.delocalised_pi_assignments.end(),
            [](const auto& assignment){return assignment.cyclic_topology;});
        if (!expected) append_failure(failures,"invalid pi_cyclic field");
        else if (actual!=*expected) {
            append_failure(failures,std::string("pi cyclic expected ")+
                (*expected?"true":"false")+", got "+
                (actual?"true":"false"));
        }
    }
    if (!unavailable(cell(row,"pi_atoms"))) {
        for (const auto& token:split(cell(row,"pi_atoms"),';')) {
            std::string error;
            const auto expected=parse_pi_atom_set(token,error);
            if (!expected) {
                append_failure(failures,error);
                continue;
            }
            const bool found=std::any_of(
                wavefunction.delocalised_pi_assignments.begin(),
                wavefunction.delocalised_pi_assignments.end(),
                [&](const auto& assignment){
                    auto atoms=assignment.atoms;
                    std::sort(atoms.begin(),atoms.end());
                    return atoms==*expected;
                });
            if (!found) append_failure(failures,"expected pi atom family not found: "+token);
        }
    }
    if (!unavailable(cell(row,"pi_electrons"))) {
        const auto expected=parse_real(cell(row,"pi_electrons"));
        double actual=0.0;
        for (const auto& assignment:wavefunction.delocalised_pi_assignments) {
            actual+=assignment.electron_count;
        }
        if (!expected) append_failure(failures,"invalid pi_electrons field");
        else if (std::abs(actual-*expected)>0.05) {
            append_failure(failures,"pi electrons expected "+
                std::to_string(*expected)+", got "+std::to_string(actual));
        }
    }

    const cov::InteractionGraph graph = cov::build_interaction_graph(wavefunction);
    std::string parse_error;
    const auto required = parse_pairs(cell(row, "required_bonds"), parse_error);
    if (!required) {
        append_failure(failures, parse_error);
    } else {
        for (const auto pair : *required) {
            if (pair.second >= wavefunction.atoms.size()) {
                append_failure(failures, "required bond " + one_based_pair(pair) +
                                             " exceeds atom count");
            } else if (!strong_pair(graph, pair)) {
                append_failure(failures, "missing required strong bond " +
                                             one_based_pair(pair));
            }
        }
    }
    parse_error.clear();
    const auto forbidden = parse_pairs(cell(row, "forbidden_bonds"), parse_error);
    if (!forbidden) {
        append_failure(failures, parse_error);
    } else {
        for (const auto pair : *forbidden) {
            if (pair.second >= wavefunction.atoms.size()) {
                append_failure(failures, "forbidden bond " + one_based_pair(pair) +
                                             " exceeds atom count");
            } else if (strong_pair(graph, pair)) {
                append_failure(failures, "forbidden strong bond present " +
                                             one_based_pair(pair));
            }
        }
    }

    if (!unavailable(cell(row, "fragment_count"))) {
        const auto expected = parse_integer(cell(row, "fragment_count"));
        if (!expected || *expected < 0) {
            append_failure(failures, "invalid fragment_count field");
        } else if (graph.fragment_analysis.fragments.size() !=
                   static_cast<std::size_t>(*expected)) {
            std::ostringstream detail;
            detail << "fragment_count expected " << *expected << ", got "
                   << graph.fragment_analysis.fragments.size() << "; Mayer";
            for (const auto& record : wavefunction.bond_orders) {
                detail << ' ' << (record.atom_a + 1u) << '-'
                       << (record.atom_b + 1u) << '=' << record.mayer_order;
            }
            detail << "; strong";
            for (const auto& edge : graph.edges) {
                if (!cov::interaction_merges_fragments(edge.kind)) continue;
                detail << ' ' << (edge.atom_a + 1u) << '-'
                       << (edge.atom_b + 1u) << ':'
                       << cov::interaction_kind_name(edge.kind);
            }
            append_failure(failures, detail.str());
        }
    }

    if (!unavailable(cell(row, "required_geometry"))) {
        const auto expected = geometry_expectations(cell(row, "required_geometry"));
        const auto environments = cov::analyse_local_molecular_geometries(wavefunction);
        std::size_t maximum_cn = 0u;
        for (const auto& environment : environments) {
            maximum_cn = std::max(maximum_cn, environment.coordination_number());
        }
        const bool local_match = std::any_of(
            environments.begin(), environments.end(), [&](const auto& environment) {
                return environment.coordination_number() == maximum_cn &&
                       geometry_matches(environment, expected);
            });
        const std::string used_point_group = semantic_key(wavefunction.point_group_used);
        const std::string detected_point_group = semantic_key(wavefunction.point_group_detected);
        const bool global_match = std::any_of(
            expected.begin(), expected.end(), [&](const auto& key) {
                return (!used_point_group.empty() && key == used_point_group) ||
                       (!detected_point_group.empty() && key == detected_point_group);
            });
        const bool cage_match = !graph.polyhedral_cages.empty() &&
            std::any_of(expected.begin(), expected.end(), [](const auto& key) {
                return key == "icosahedralcage" || key == "icosahedral";
            });
        if (!local_match && !global_match && !cage_match) {
            append_failure(failures, "required global/principal-local geometry '" +
                                         cell(row, "required_geometry") + "' not found");
        }
    }

    bool require_no_multicentre = false;
    parse_error.clear();
    const auto multicentre = parse_multicentre(cell(row, "multicentre_kind"),
                                               require_no_multicentre, parse_error);
    if (!multicentre) {
        append_failure(failures, parse_error);
    } else if (require_no_multicentre && !graph.multicentre_groups.empty()) {
        append_failure(failures, "unexpected strict multicentre assignment");
    } else {
        for (const auto& expected : *multicentre) {
            const bool matched = std::any_of(
                graph.multicentre_groups.begin(), graph.multicentre_groups.end(),
                [&](const auto& actual) {
                    if (actual.kind != expected.kind) return false;
                    if (expected.atoms.empty()) return true;
                    std::vector<std::uint32_t> atoms = actual.atoms;
                    std::sort(atoms.begin(), atoms.end());
                    return atoms == expected.atoms;
                });
            if (!matched) {
                std::string description = cov::multicentre_kind_name(expected.kind);
                if (!expected.atoms.empty()) {
                    description += ":";
                    for (std::size_t index = 0u; index < expected.atoms.size(); ++index) {
                        if (index != 0u) description += "-";
                        description += std::to_string(expected.atoms[index] + 1u);
                    }
                }
                std::ostringstream detail;
                detail << "missing strict multicentre assignment " << description
                       << "; candidates="
                       << wavefunction.multicentre_candidates.size();
                if (!expected.atoms.empty()) {
                    detail << "; pair-Mayer=";
                    for (std::size_t i = 0u; i < expected.atoms.size(); ++i) {
                        for (std::size_t j = i + 1u; j < expected.atoms.size(); ++j) {
                            double value = 0.0;
                            for (const auto& record : wavefunction.bond_orders) {
                                if ((record.atom_a == expected.atoms[i] &&
                                     record.atom_b == expected.atoms[j]) ||
                                    (record.atom_a == expected.atoms[j] &&
                                     record.atom_b == expected.atoms[i])) {
                                    value = record.mayer_order;
                                    break;
                                }
                            }
                            detail << (expected.atoms[i] + 1u) << '-'
                                   << (expected.atoms[j] + 1u) << ':' << value << ',';
                        }
                    }
                }
                detail << "; assignments="
                       << wavefunction.multicentre_assignments.size();
                for (const auto& assignment : wavefunction.multicentre_assignments) {
                    detail << " [" << cov::multicentre_kind_name(assignment.kind)
                           << " e=" << assignment.electron_count << " atoms=";
                    for (const auto atom : assignment.atoms) detail << (atom + 1u) << ',';
                    detail << " MOs=";
                    for (const auto orbital : assignment.orbitals) {
                        detail << (orbital + 1u) << ',';
                    }
                    detail << ']';
                }
                const std::size_t candidate_limit = std::min<std::size_t>(
                    6u, wavefunction.multicentre_candidates.size());
                for (std::size_t candidate_index = 0u;
                     candidate_index < candidate_limit; ++candidate_index) {
                    const auto& candidate =
                        wavefunction.multicentre_candidates[candidate_index];
                    detail << " [MO" << (candidate.orbital_index + 1u)
                           << " occ=" << candidate.occupation << " atoms=";
                    const std::size_t atom_limit = std::min<std::size_t>(
                        4u, candidate.atoms.size());
                    for (std::size_t atom_index = 0u;
                         atom_index < atom_limit; ++atom_index) {
                        if (atom_index != 0u) detail << ',';
                        detail << (candidate.atoms[atom_index] + 1u) << ':'
                               << candidate.participation[atom_index];
                    }
                    detail << ']';
                }
                append_failure(failures, detail.str());
            }
        }
    }

    (void)manifest;
    return failures;
}

} // namespace

int main(int argc, char** argv) {
    const bool dump_chemistry=argc==3 && std::string_view(argv[2])=="--dump-chemistry";
    if (argc != 2 && !dump_chemistry) {
        std::cout << "usage: cov_reference_manifest_regression <manifest.tsv> "
                     "[--dump-chemistry]\n";
        return 2;
    }

    const std::filesystem::path manifest = std::filesystem::absolute(argv[1]);
    std::ifstream input(manifest, std::ios::binary);
    if (!input) {
        std::cout << "FAIL\tmanifest\tcannot open " << manifest.string() << '\n';
        return 2;
    }

    std::string header_line;
    if (!std::getline(input, header_line)) {
        std::cout << "FAIL\tmanifest\tempty TSV\n";
        return 2;
    }
    if (header_line.size() >= 3u &&
        static_cast<unsigned char>(header_line[0]) == 0xefu &&
        static_cast<unsigned char>(header_line[1]) == 0xbbu &&
        static_cast<unsigned char>(header_line[2]) == 0xbfu) {
        header_line.erase(0u, 3u);
    }
    auto headers = split_tsv(header_line);
    std::set<std::string> unique_headers;
    for (auto& header : headers) {
        header = lower(trim(std::move(header)));
        if (header.empty() || !unique_headers.insert(header).second) {
            std::cout << "FAIL\tmanifest\tempty or duplicate header '" << header << "'\n";
            return 2;
        }
    }
    const std::vector<std::string> required_headers{
        "case_id", "path", "formula", "charge", "multiplicity", "family",
        "required_bonds", "forbidden_bonds", "fragment_count", "required_geometry",
        "multicentre_kind", "notes", "gjf_path", "log_path", "basis", "method",
        "runnable", "execution_state", "normal_termination", "formchk_state",
    };
    for (const auto& required : required_headers) {
        if (!unique_headers.contains(required)) {
            std::cout << "FAIL\tmanifest\tmissing required column " << required << '\n';
            return 2;
        }
    }

    std::size_t rows = 0u;
    std::size_t passed = 0u;
    std::size_t failed = 0u;
    std::size_t skipped = 0u;
    std::set<std::string> case_ids;
    std::string line;
    std::size_t line_number = 1u;
    while (std::getline(input, line)) {
        ++line_number;
        if (trim(line).empty() || trim(line).starts_with('#')) continue;
        ++rows;
        auto fields = split_tsv(line);
        if (fields.size() > headers.size()) {
            ++failed;
            std::cout << "FAIL\tline-" << line_number
                      << "\tfield count at most " << headers.size()
                      << ", got " << fields.size() << '\n';
            continue;
        }
        // Trailing columns are optional validation contracts.  Older rows may
        // omit them; an omitted trailing value is semantically unavailable.
        fields.resize(headers.size());
        Row row;
        for (std::size_t index = 0u; index < headers.size(); ++index) {
            row.emplace(headers[index], fields[index]);
        }
        const std::string case_id = trim(cell(row, "case_id"));
        if (case_id.empty() || !case_ids.insert(case_id).second) {
            ++failed;
            std::cout << "FAIL\tline-" << line_number
                      << "\tempty or duplicate case_id '" << case_id << "'\n";
            continue;
        }

        const auto runnable = parse_boolean(cell(row, "runnable"));
        if (!runnable) {
            ++failed;
            std::cout << "FAIL\t" << case_id << "\tinvalid runnable value '"
                      << cell(row, "runnable") << "'\n";
            continue;
        }
        if (!*runnable) {
            ++skipped;
            std::cout << "SKIP\t" << case_id << "\tnot runnable\n";
            continue;
        }
        const auto normal = parse_boolean(cell(row, "normal_termination"));
        if ((!unavailable(cell(row, "normal_termination")) && normal && !*normal) ||
            explicitly_not_ready(cell(row, "execution_state")) ||
            explicitly_not_ready(cell(row, "formchk_state"))) {
            ++skipped;
            std::cout << "SKIP\t" << case_id << "\tcalculation/FCHK not successful\n";
            continue;
        }

        const std::filesystem::path fchk_path = resolve_path(manifest, cell(row, "path"));
        std::error_code file_error;
        if (fchk_path.empty() || !std::filesystem::is_regular_file(fchk_path, file_error)) {
            ++failed;
            std::cout << "FAIL\t" << case_id << "\tsuccessful row has no FCHK file '"
                      << fchk_path.string() << "'\n";
            continue;
        }

        try {
            cov::WavefunctionParseOptions options;
            options.max_atoms = 2048u;
            options.require_orbitals = true;
            options.keep_density = true;
            options.reconstruct_density_if_missing = true;
            const std::filesystem::path log_path =
                resolve_path(manifest, cell(row, "log_path"));
            file_error.clear();
            if (!log_path.empty() && std::filesystem::is_regular_file(log_path, file_error)) {
                options.gaussian_log_path = log_path;
            }
            const cov::Wavefunction wavefunction = cov::parse_wavefunction(fchk_path, options);
            if (dump_chemistry) {
                std::cout<<"CHEM\t"<<case_id<<"\tPG="
                         <<wavefunction.point_group_detected<<"\tPI=";
                if (wavefunction.delocalised_pi_assignments.empty()) {
                    std::cout<<"none";
                }
                for (std::size_t assignment_index=0u;
                     assignment_index<wavefunction.delocalised_pi_assignments.size();
                     ++assignment_index) {
                    if (assignment_index!=0u) std::cout<<';';
                    const auto& assignment=
                        wavefunction.delocalised_pi_assignments[assignment_index];
                    std::cout<<'[';
                    for (std::size_t atom=0u;atom<assignment.atoms.size();++atom) {
                        if (atom!=0u) std::cout<<',';
                        std::cout<<(assignment.atoms[atom]+1u);
                    }
                    std::cout<<"]e="<<assignment.electron_count
                             <<",ch="<<assignment.orientation_channels.size()
                             <<",cyc="<<(assignment.cyclic_topology?1:0)
                             <<",MO=";
                    for (std::size_t orbital=0u;
                         orbital<assignment.orbitals.size();++orbital) {
                        if (orbital!=0u) std::cout<<',';
                        std::cout<<(assignment.orbitals[orbital]+1u);
                    }
                }
                std::cout<<"\tMC=";
                if (wavefunction.multicentre_assignments.empty()) {
                    std::cout<<"none";
                }
                for (std::size_t assignment_index=0u;
                     assignment_index<wavefunction.multicentre_assignments.size();
                     ++assignment_index) {
                    if (assignment_index!=0u) std::cout<<';';
                    const auto& assignment=
                        wavefunction.multicentre_assignments[assignment_index];
                    std::cout<<cov::multicentre_kind_name(assignment.kind)<<'[';
                    for (std::size_t atom=0u;atom<assignment.atoms.size();++atom) {
                        if (atom!=0u) std::cout<<',';
                        std::cout<<(assignment.atoms[atom]+1u);
                    }
                    std::cout<<"]e="<<assignment.electron_count<<",MO=";
                    for (std::size_t orbital=0u;
                         orbital<assignment.orbitals.size();++orbital) {
                        if (orbital!=0u) std::cout<<',';
                        std::cout<<(assignment.orbitals[orbital]+1u);
                    }
                }
                const auto graph=cov::build_interaction_graph(wavefunction);
                std::cout<<"\tEDGES=";
                for (std::size_t edge_index=0u;edge_index<graph.edges.size();++edge_index) {
                    if (edge_index!=0u) std::cout<<';';
                    const auto& edge=graph.edges[edge_index];
                    std::cout<<(edge.atom_a+1u)<<'-'<<(edge.atom_b+1u)<<':'
                             <<cov::interaction_kind_name(edge.kind);
                }
                std::cout<<'\n';
            }
            const auto failures = validate_semantics(row, manifest, wavefunction);
            if (failures.empty()) {
                ++passed;
                std::cout << "PASS\t" << case_id << "\t" << fchk_path.string() << '\n';
            } else {
                ++failed;
                std::cout << "FAIL\t" << case_id << "\t" << join_failures(failures) << '\n';
            }
        } catch (const std::exception& error) {
            ++failed;
            std::cout << "FAIL\t" << case_id << "\texception: " << error.what() << '\n';
        }
    }

    std::cout << "SUMMARY\trows=" << rows << "\tpassed=" << passed
              << "\tfailed=" << failed << "\tskipped=" << skipped << '\n';
    return failed == 0u ? EXIT_SUCCESS : EXIT_FAILURE;
}
