#include "cov/fchk_parser.hpp"

#include "cov/density.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

std::string field_context(const std::string& label,
                          const std::size_t line_number) {
    return " for field '" + label + "' near line " +
           std::to_string(line_number);
}

double parse_real(std::string token,
                  const std::string& label,
                  const std::size_t line_number) {
    for (char& c : token) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    std::size_t used = 0;
    try {
        const double value = std::stod(token, &used);
        if (used != token.size()) {
            throw std::runtime_error("Invalid FCHK real value" +
                                     field_context(label, line_number) +
                                     ": '" + token + "'");
        }
        return value;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Invalid FCHK real value" +
                                 field_context(label, line_number) +
                                 ": '" + token + "'");
    } catch (const std::out_of_range&) {
        throw std::runtime_error("Out-of-range FCHK real value" +
                                 field_context(label, line_number) +
                                 ": '" + token + "'");
    }
}

long long parse_integer(const std::string& token,
                        const std::string& label,
                        const std::size_t line_number) {
    std::size_t used = 0;
    try {
        const long long value = std::stoll(token, &used);
        if (used != token.size()) {
            throw std::runtime_error("Invalid FCHK integer value" +
                                     field_context(label, line_number) +
                                     ": '" + token + "'");
        }
        return value;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Invalid FCHK integer value" +
                                 field_context(label, line_number) +
                                 ": '" + token + "'");
    } catch (const std::out_of_range&) {
        throw std::runtime_error("Out-of-range FCHK integer value" +
                                 field_context(label, line_number) +
                                 ": '" + token + "'");
    }
}

struct Header {
    std::string label;
    char type = '\0';
    bool array = false;
    std::size_t count = 0;
    std::string scalar;
};

bool parse_header(const std::string& line, Header& header) {
    if (line.size() < 41u) return false;

    header = {};
    header.label = trim(line.substr(0, 40));
    if (header.label.empty()) return false;

    const std::string tail = line.substr(40);
    std::istringstream stream(tail);
    stream >> header.type;
    if (!stream || (header.type != 'I' && header.type != 'R' &&
                    header.type != 'C' && header.type != 'L')) {
        return false;
    }

    const auto npos = tail.find("N=");
    if (npos != std::string::npos) {
        header.array = true;
        std::istringstream count_stream(tail.substr(npos + 2u));
        unsigned long long count = 0;
        count_stream >> count;
        if (!count_stream || count > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("Invalid FCHK array length for field: " + header.label);
        }
        header.count = static_cast<std::size_t>(count);
    } else {
        std::string scalar;
        std::getline(stream, scalar);
        header.scalar = trim(scalar);
    }
    return true;
}

struct FchkRecords {
    std::unordered_map<std::string, long long> integers;
    std::unordered_map<std::string, double> reals;
    std::unordered_map<std::string, std::string> strings;
    std::unordered_map<std::string, std::vector<long long>> integer_arrays;
    std::unordered_map<std::string, std::vector<double>> real_arrays;
};

void reject_unsupported_wavefunction_record(const Header& header,
                                            const std::size_t line_number) {
    std::string label = header.label;
    std::transform(label.begin(), label.end(), label.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    // Unknown auxiliary fields remain legal. Only explicit wavefunction data
    // that cannot be represented by our real scalar orbitals is rejected;
    // e.g. an imaginary-frequency diagnostic is not a complex wavefunction.
    const bool orbital_data = label.find("mo coefficient") != std::string::npos ||
                              label.find("orbital coefficient") != std::string::npos ||
                              label.find("spinor coefficient") != std::string::npos ||
                              label.find("scf density") != std::string::npos ||
                              label.find("density matrix") != std::string::npos;
    const bool unsupported_representation = label.find("imaginary") != std::string::npos ||
                                            label.find("complex") != std::string::npos ||
                                            label.find("spinor") != std::string::npos;
    if (orbital_data && unsupported_representation) {
        throw std::runtime_error("Unsupported FCHK complex/spinor wavefunction data" +
                                 field_context(header.label, line_number));
    }
    // FCHK canonical occupations currently come from integer spin counts.
    // Do not silently discard an explicit contrary occupation specification.
    for (const auto* prefix : {"", "alpha ", "beta "}) {
        for (const auto* name : {"orbital occupations", "orbital occupancies",
                                 "orbital occupation numbers", "mo occupations",
                                 "mo occupancies", "mo occupation numbers"}) {
            if (label == std::string(prefix) + name) {
                throw std::runtime_error("Unsupported explicit FCHK orbital occupations" +
                                         field_context(header.label, line_number));
            }
        }
    }
}

FchkRecords read_records(std::ifstream& input) {
    FchkRecords records;
    std::string line;
    std::size_t line_number = 2u; // title and route-method lines were consumed
    while (std::getline(input, line)) {
        ++line_number;
        Header header;
        if (!parse_header(line, header)) continue;
        reject_unsupported_wavefunction_record(header, line_number);

        if (!header.array) {
            if (header.type == 'I') {
                records.integers[header.label] =
                    parse_integer(header.scalar, header.label, line_number);
            } else if (header.type == 'R') {
                records.reals[header.label] =
                    parse_real(header.scalar, header.label, line_number);
            } else {
                records.strings[header.label] = header.scalar;
            }
            continue;
        }

        // Gaussian character arrays are fixed-width A12 values, five values
        // per physical line.  Consume them even though COV does not currently
        // retain them: a continuation line can legitimately contain I/R/C/L
        // in column 41 and must never be mistaken for a new record header.
        if (header.type == 'C') {
            std::size_t values = 0;
            while (values < header.count && std::getline(input, line)) {
                ++line_number;
                values += (line.size() + 11u) / 12u;
            }
            if (values < header.count) {
                throw std::runtime_error(
                    "Unexpected end of FCHK character array '" + header.label +
                    "' near line " + std::to_string(line_number));
            }
            continue;
        }

        struct TokenAtLine {
            std::string value;
            std::size_t line = 0;
        };
        std::vector<TokenAtLine> tokens;
        tokens.reserve(header.count);
        while (tokens.size() < header.count && std::getline(input, line)) {
            ++line_number;
            std::istringstream values(line);
            std::string token;
            while (values >> token) {
                tokens.push_back({token, line_number});
                if (tokens.size() == header.count) break;
            }
        }
        if (tokens.size() != header.count) {
            throw std::runtime_error("Unexpected end of FCHK array: " + header.label);
        }

        if (header.type == 'L') {
            // Logical arrays have already been consumed token-for-token.
            continue;
        }
        if (header.type == 'I') {
            auto& values = records.integer_arrays[header.label];
            values.reserve(header.count);
            for (const auto& token : tokens) {
                values.push_back(parse_integer(
                    token.value, header.label, token.line));
            }
        } else {
            auto& values = records.real_arrays[header.label];
            values.reserve(header.count);
            for (const auto& token : tokens) {
                values.push_back(parse_real(
                    token.value, header.label, token.line));
            }
        }
    }
    return records;
}

const std::vector<long long>& require_int_array(const FchkRecords& records,
                                                const char* label) {
    const auto it = records.integer_arrays.find(label);
    if (it == records.integer_arrays.end()) {
        throw std::runtime_error(std::string("Required FCHK field is missing: ") + label);
    }
    return it->second;
}

const std::vector<double>& require_real_array(const FchkRecords& records,
                                              const char* label) {
    const auto it = records.real_arrays.find(label);
    if (it == records.real_arrays.end()) {
        throw std::runtime_error(std::string("Required FCHK field is missing: ") + label);
    }
    return it->second;
}

long long require_integer(const FchkRecords& records, const char* label) {
    const auto it = records.integers.find(label);
    if (it == records.integers.end()) {
        throw std::runtime_error(std::string("Required FCHK field is missing: ") + label);
    }
    return it->second;
}

const long long* find_integer(const FchkRecords& records, const char* label) {
    const auto it = records.integers.find(label);
    return it == records.integers.end() ? nullptr : &it->second;
}

const std::vector<double>* find_real_array(const FchkRecords& records,
                                           const char* label) {
    const auto it = records.real_arrays.find(label);
    return it == records.real_arrays.end() ? nullptr : &it->second;
}

std::string atomic_symbol(const int z) {
    static constexpr std::array<const char*, 119> symbols{{
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    }};
    if (z <= 0 || z >= static_cast<int>(symbols.size())) {
        return "X" + std::to_string(z);
    }
    return symbols[static_cast<std::size_t>(z)];
}

std::uint32_t checked_u32(const long long value, const char* label) {
    if (value < 0 || static_cast<unsigned long long>(value) >
                         std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error(std::string("FCHK integer is out of range: ") + label);
    }
    return static_cast<std::uint32_t>(value);
}

std::int32_t checked_i32(const long long value, const char* label) {
    if (value < std::numeric_limits<std::int32_t>::min() ||
        value > std::numeric_limits<std::int32_t>::max()) {
        throw std::runtime_error(std::string("FCHK integer is out of range: ") + label);
    }
    return static_cast<std::int32_t>(value);
}

void append_shell(Wavefunction& wf,
                  const std::uint32_t atom_index,
                  const std::uint8_t l,
                  const bool pure,
                  const std::vector<double>& exponents,
                  const std::vector<double>& coefficients,
                  const std::size_t begin,
                  const std::size_t count) {
    if (begin + count > exponents.size() || begin + count > coefficients.size()) {
        throw std::runtime_error("FCHK primitive/contraction array is shorter than shell metadata");
    }

    Shell shell;
    shell.atom_index = atom_index;
    shell.primitive_offset = static_cast<std::uint32_t>(wf.primitives.size());
    shell.primitive_count = static_cast<std::uint32_t>(count);
    shell.basis_offset = wf.basis_count;
    shell.angular_momentum = l;
    shell.pure = pure ? 1u : 0u;

    for (std::size_t p = 0; p < count; ++p) {
        Primitive primitive;
        primitive.exponent = exponents[begin + p];
        primitive.coefficient = coefficients[begin + p];
        wf.primitives.push_back(primitive);
    }

    wf.basis_count += shell_basis_count(shell);
    wf.shells.push_back(shell);
}

using BasisMapEntry = GaussianAoTransformEntry;

void append_basis_map(std::vector<BasisMapEntry>& basis_map,
                      const std::size_t source_offset,
                      const int shell_type) {
    const auto local = gaussian_ao_transform(shell_type, source_offset);
    basis_map.insert(basis_map.end(), local.begin(), local.end());
}

std::vector<double> transform_packed_density(const std::vector<double>& source,
                                             const std::vector<BasisMapEntry>& basis_map) {
    const std::size_t n = basis_map.size();
    const std::size_t expected = n * (n + 1u) / 2u;
    if (source.size() != expected) {
        throw std::runtime_error("FCHK packed density dimension does not match basis count");
    }

    auto packed = [](std::size_t i, std::size_t j) {
        if (j > i) std::swap(i, j);
        return i * (i + 1u) / 2u + j;
    };

    std::vector<double> transformed(expected, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const BasisMapEntry& mi = basis_map[i];
            const BasisMapEntry& mj = basis_map[j];
            transformed[packed(i, j)] =
                mi.coefficient_scale * mj.coefficient_scale *
                source[packed(mi.source_index, mj.source_index)];
        }
    }
    return transformed;
}

void set_global_purity_flags(Wavefunction& wf) {
    auto all_pure = [&](const std::uint8_t l) {
        bool seen = false;
        bool pure = true;
        for (const Shell& shell : wf.shells) {
            if (shell.angular_momentum != l) continue;
            seen = true;
            pure = pure && shell.pure != 0;
        }
        return seen && pure;
    };
    wf.pure_d = all_pure(2);
    wf.pure_f = all_pure(3);
    wf.pure_g = all_pure(4);
}

void append_orbital_set(Wavefunction& wf,
                        const std::vector<double>& energies,
                        const std::vector<double>& coefficients,
                        const std::vector<BasisMapEntry>& basis_map,
                        const Spin spin,
                        const std::uint32_t occupied_count,
                        const bool restricted) {
    const std::size_t basis = wf.basis_count;
    if (basis == 0) throw std::runtime_error("FCHK orbital data has zero basis functions");
    if (basis_map.size() != basis) {
        throw std::runtime_error("Internal FCHK basis-order map does not match basis count");
    }
    if (coefficients.size() != energies.size() * basis) {
        throw std::runtime_error(
            "FCHK MO coefficient dimension does not equal orbital_count * basis_count");
    }

    wf.orbitals.reserve(wf.orbitals.size() + energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        MolecularOrbital mo;
        mo.energy_hartree = energies[i];
        mo.spin = spin;
        mo.spin_provenance = DataProvenance::Producer;
        mo.spin_source_text = spin == Spin::Alpha ? "Alpha MO coefficients" : "Beta MO coefficients";
        mo.source_orbital_index = i;
        if (restricted) {
            const float alpha = i < wf.alpha_electrons ? 1.0f : 0.0f;
            const float beta = i < wf.beta_electrons ? 1.0f : 0.0f;
            mo.occupation = alpha + beta;
        } else {
            mo.occupation = i < occupied_count ? 1.0f : 0.0f;
        }
        mo.occupation_provenance = DataProvenance::Derived;
        mo.coefficients.resize(basis);
        const std::size_t offset = i * basis;
        mo.gaussian_source_coefficients.assign(coefficients.begin()+offset,
                                               coefficients.begin()+offset+basis);
        for (std::size_t internal = 0; internal < basis; ++internal) {
            const BasisMapEntry& map = basis_map[internal];
            mo.coefficients[internal] =
                map.coefficient_scale * coefficients[offset + map.source_index];
        }
        wf.orbitals.push_back(std::move(mo));
    }
}

void validate_density(const std::vector<double>& density,
                      const std::size_t basis,
                      const char* label) {
    if (density.empty()) return;
    const std::size_t expected = basis * (basis + 1u) / 2u;
    if (density.size() != expected) {
        throw std::runtime_error(std::string(label) +
                                 " dimension does not match packed AO basis dimension");
    }
}

} // namespace

Wavefunction parse_fchk(const std::filesystem::path& path,
                        const FchkParseOptions& options) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Unable to open FCHK file: " + path.string());
    }

    Wavefunction wf;
    wf.source = WavefunctionSource::Fchk;
    std::getline(input, wf.source_title);
    std::getline(input, wf.source_route);
    wf.source_title = trim(wf.source_title);
    wf.source_route = trim(wf.source_route);

    const FchkRecords records = read_records(input);

    const auto atom_count_ll = require_integer(records, "Number of atoms");
    if (atom_count_ll <= 0) throw std::runtime_error("FCHK Number of atoms must be positive");
    const std::size_t atom_count = static_cast<std::size_t>(atom_count_ll);
    if (atom_count > options.max_atoms) {
        throw std::runtime_error("FCHK atom count exceeds configured maximum of " +
                                 std::to_string(options.max_atoms));
    }

    const auto& atomic_numbers = require_int_array(records, "Atomic numbers");
    const auto& coordinates = require_real_array(records, "Current cartesian coordinates");
    if (atomic_numbers.size() != atom_count || coordinates.size() != atom_count * 3u) {
        throw std::runtime_error("FCHK atom-coordinate dimensions are inconsistent");
    }

    wf.atoms.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        Atom atom;
        atom.atomic_number = static_cast<int>(atomic_numbers[i]);
        atom.symbol = atomic_symbol(atom.atomic_number);
        atom.x = coordinates[3u * i + 0u];
        atom.y = coordinates[3u * i + 1u];
        atom.z = coordinates[3u * i + 2u];
        wf.atoms.push_back(std::move(atom));
    }

    wf.alpha_electrons = checked_u32(require_integer(records, "Number of alpha electrons"),
                                     "Number of alpha electrons");
    wf.beta_electrons = checked_u32(require_integer(records, "Number of beta electrons"),
                                    "Number of beta electrons");
    wf.electron_counts_provenance = DataProvenance::Producer;
    if (const auto* charge = find_integer(records, "Charge")) {
        wf.charge = checked_i32(*charge, "Charge");
        wf.charge_provenance = DataProvenance::Producer;
    }
    if (const auto* multiplicity = find_integer(records, "Multiplicity")) {
        if (*multiplicity <= 0) {
            throw std::runtime_error("FCHK Multiplicity must be positive");
        }
        wf.multiplicity = checked_u32(*multiplicity, "Multiplicity");
        wf.multiplicity_provenance = DataProvenance::Producer;
    }
    if (const auto* charges = find_real_array(records, "Mulliken Charges")) {
        if (charges->size() == atom_count &&
            std::all_of(charges->begin(), charges->end(),
                        [](const double value) { return std::isfinite(value); })) {
            wf.atomic_partial_charges = *charges;
            wf.atomic_partial_charge_scheme = "Mulliken";
            wf.atomic_partial_charge_provenance = DataProvenance::Producer;
        }
    }

    const auto& shell_types = require_int_array(records, "Shell types");
    const auto& primitive_counts = require_int_array(records, "Number of primitives per shell");
    const auto& shell_to_atom = require_int_array(records, "Shell to atom map");
    const auto& exponents = require_real_array(records, "Primitive exponents");
    const auto& contractions = require_real_array(records, "Contraction coefficients");
    const auto* sp_contractions = find_real_array(records, "P(S=P) Contraction coefficients");

    if (shell_types.size() != primitive_counts.size() ||
        shell_types.size() != shell_to_atom.size()) {
        throw std::runtime_error("FCHK shell metadata arrays have inconsistent lengths");
    }
    if (exponents.size() != contractions.size()) {
        throw std::runtime_error("FCHK primitive exponent/contraction arrays have inconsistent lengths");
    }

    std::vector<BasisMapEntry> basis_map;
    std::size_t primitive_cursor = 0;
    std::size_t source_basis_cursor = 0;
    for (std::size_t s = 0; s < shell_types.size(); ++s) {
        const long long count_ll = primitive_counts[s];
        if (count_ll <= 0) throw std::runtime_error("FCHK shell has no primitives");
        const std::size_t count = static_cast<std::size_t>(count_ll);
        if (primitive_cursor + count > exponents.size()) {
            throw std::runtime_error("FCHK shell primitive counts exceed primitive arrays");
        }

        const long long atom_map = shell_to_atom[s];
        if (atom_map <= 0 || static_cast<std::size_t>(atom_map) > wf.atoms.size()) {
            throw std::runtime_error("FCHK shell-to-atom map contains an invalid atom index");
        }
        const std::uint32_t atom_index = static_cast<std::uint32_t>(atom_map - 1);
        const int shell_type = static_cast<int>(shell_types[s]);

        if (shell_type == -1) {
            if (!sp_contractions || sp_contractions->size() != contractions.size()) {
                throw std::runtime_error(
                    "FCHK SP shell requires P(S=P) Contraction coefficients");
            }
            append_shell(wf, atom_index, 0u, false,
                         exponents, contractions, primitive_cursor, count);
            append_shell(wf, atom_index, 1u, false,
                         exponents, *sp_contractions, primitive_cursor, count);
            append_basis_map(basis_map, source_basis_cursor, shell_type);
            source_basis_cursor += 4u;
        } else {
            const int l = shell_type < 0 ? -shell_type : shell_type;
            if (l < 0 || l > 4) {
                throw std::runtime_error(
                    "FCHK shell angular momentum above g is not supported yet: type=" +
                    std::to_string(shell_type));
            }
            const bool pure = shell_type <= -2;
            append_shell(wf, atom_index, static_cast<std::uint8_t>(l), pure,
                         exponents, contractions, primitive_cursor, count);
            append_basis_map(basis_map, source_basis_cursor, shell_type);
            source_basis_cursor += shell_basis_count(wf.shells.back());
        }
        primitive_cursor += count;
    }
    if (primitive_cursor != exponents.size()) {
        throw std::runtime_error("FCHK primitive arrays contain unused entries");
    }

    const auto expected_basis = checked_u32(require_integer(records, "Number of basis functions"),
                                            "Number of basis functions");
    if (wf.basis_count != expected_basis || source_basis_cursor != expected_basis ||
        basis_map.size() != expected_basis) {
        throw std::runtime_error(
            "FCHK shell-derived basis count/map does not match Number of basis functions " +
            std::to_string(expected_basis));
    }
    set_global_purity_flags(wf);
    wf.gaussian_ao_transform = basis_map;

    const auto* alpha_energies = find_real_array(records, "Alpha Orbital Energies");
    const auto* alpha_coefficients = find_real_array(records, "Alpha MO coefficients");
    const auto* beta_energies = find_real_array(records, "Beta Orbital Energies");
    const auto* beta_coefficients = find_real_array(records, "Beta MO coefficients");

    if ((alpha_energies == nullptr) != (alpha_coefficients == nullptr)) {
        throw std::runtime_error("FCHK alpha orbital energies/coefficients are incomplete");
    }
    if ((beta_energies == nullptr) != (beta_coefficients == nullptr)) {
        throw std::runtime_error("FCHK beta orbital energies/coefficients are incomplete");
    }

    if (alpha_energies && alpha_coefficients) {
        const bool unrestricted = beta_energies && beta_coefficients;
        wf.orbital_occupation_model = unrestricted ? OrbitalOccupationModel::ExplicitSpin :
                                                    OrbitalOccupationModel::CanonicalShared;
        wf.orbital_occupation_model_provenance = DataProvenance::Producer;
        append_orbital_set(wf, *alpha_energies, *alpha_coefficients, basis_map,
                           Spin::Alpha, wf.alpha_electrons, !unrestricted);
        if (unrestricted) {
            append_orbital_set(wf, *beta_energies, *beta_coefficients, basis_map,
                               Spin::Beta, wf.beta_electrons, false);
        }
    }

    if (options.require_orbitals && wf.orbitals.empty()) {
        throw std::runtime_error("No molecular orbitals found in FCHK file");
    }

    if (options.keep_density) {
        if (const auto* total = find_real_array(records, "Total SCF Density")) {
            validate_density(*total, wf.basis_count, "Total SCF Density");
            wf.total_density_packed = transform_packed_density(*total, basis_map);
            wf.total_density_provenance = DataProvenance::Producer;
        }

        if (const auto* spin = find_real_array(records, "Spin SCF Density")) {
            validate_density(*spin, wf.basis_count, "Spin SCF Density");
            wf.spin_density_packed = transform_packed_density(*spin, basis_map);
            wf.spin_density_provenance = DataProvenance::Producer;
        }
        establish_density(wf, options.reconstruct_density_if_missing);
    }

    return wf;
}

} // namespace cov
