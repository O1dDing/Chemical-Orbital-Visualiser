#pragma once
#include "cov/model.hpp"
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cov {
enum class NboSpin { Total, Alpha, Beta };
struct NboSource {
    std::string path, block, producer_version, raw;
    std::size_t line_begin=0, line_end=0, analysis_segment=0;
};
struct NboPopulation {
    std::size_t atom=0; // producer one-based atom identity
    std::string symbol;
    NboSpin spin=NboSpin::Total;
    double charge=0, core=0, valence=0, rydberg=0, total=0;
    // NBO NPA includes the ECP replacement core in its printed population.
    // These are derived only when the archive reports effective charges.
    std::optional<double> effective_core_electrons, explicit_population;
    std::optional<double> spin_density;
    NboSource source;
};
struct NboNao {
    std::size_t id=0, atom=0;
    NboSpin spin=NboSpin::Total;
    std::string symbol, angular, type;
    double occupation=0;
    std::optional<double> energy_hartree; // NAO diagonal Fock energy, not canonical eigenvalue
    std::optional<double> spin_density; // total-spin NAO table prints Spin instead of Energy
    NboSource source;
};
struct NboLocalComponent {
    std::size_t atom=0;
    double percent=0, coefficient=0;
    std::string hybrid; // literal producer hybrid description, no renormalization
    NboSource source;
};
struct NboOrbital {
    std::size_t id=0, ordinal=0;
    NboSpin spin=NboSpin::Total;
    std::string kind, label;
    std::vector<std::size_t> atoms;
    double occupation=0;
    std::optional<double> diagonal_fock_hartree;
    std::optional<NboSource> energy_source;
    std::vector<NboLocalComponent> components;
    NboSource source;
};
struct NboE2 {
    std::size_t donor=0, acceptor=0;
    NboSpin spin=NboSpin::Total;
    double value=0, energy_gap_hartree=0, fock_hartree=0;
    std::string units="kcal/mol";
    std::optional<double> printing_threshold;
    NboSource source;
};
struct NboE2Section {
    NboSpin spin=NboSpin::Total;
    std::optional<double> printing_threshold;
    std::string units="kcal/mol", missing_reason;
    NboSource source;
};
struct NboWiberg {
    std::size_t atom_a=0, atom_b=0;
    NboSpin spin=NboSpin::Total;
    double value=0;
    NboSource source;
};
struct NboNlmo {
    std::size_t id=0;
    NboSpin spin=NboSpin::Total;
    double occupation=0, parent_percent=0;
    std::string parent_label;
    std::optional<std::size_t> parent_nbo; // matched by literal label, never by ordinal
    std::vector<NboLocalComponent> components;
    NboSource source;
};
struct NboMatrix {
    std::string kind; // AONBO: AO rows / NBO columns; NBOMO: NBO rows / MO columns
    NboSpin spin=NboSpin::Total;
    std::size_t rows=0, columns=0;
    std::vector<double> values; // row-major, always complete or rejected
    NboSource source;
};
struct NboArchive {
    std::vector<Atom> atoms;
    std::size_t basis_count=0;
    bool open_shell=false, density_is_bond_order=false;
    std::string fock_input_units="hartree"; // matrices are converted to hartree on read
    std::vector<int> centers, labels, ncomp, nprim, nptr;
    std::vector<double> exponents, cs, cp, cd, cf, cg;
    std::vector<NboMatrix> matrices;
    NboSource source;
};
struct NboCanonicalEvidence {
    NboSpin spin=NboSpin::Total;
    std::string coefficient_source;
    bool direct_fchk_coefficients=false, density_verified=false;
    std::string detail;
    std::vector<double> archive_to_fchk_phase; // per stored MO column; null columns have 0
};
struct NboNaoValidation {
    NboSpin spin=NboSpin::Total;
    bool available=false, direct_fchk_coefficients=false;
    std::string status="unavailable", detail;
    std::size_t effective_mo_columns=0;
    std::vector<std::string> column_status; // active or null_padding; never physical zero MOs
    std::optional<double> orthogonality_error, composition_error, projection_error;
    std::optional<double> normalization_error, nao_occupation_error, nao_nbo_composition_error;
    bool occupation_density_verified=false;
    std::optional<double> occupation_density_error;
    std::vector<double> canonical_occupations; // only populated after T f T^T density verification
};
struct NboNaoContribution {
    std::size_t nao_id=0, atom=0; // producer one-based identities
    std::string symbol, type, angular;
    std::optional<int> principal_n, angular_l; // literal producer shell label, never inferred
    double coefficient=0, weight=0;
    std::optional<double> electron_contribution;
    NboSource source;
};
struct NboNaoGroupContribution {
    std::size_t atom=0;
    std::string symbol, label;
    std::optional<int> principal_n, angular_l;
    double weight=0;
    std::optional<double> electron_contribution;
    std::vector<std::size_t> nao_ids;
};
struct NboMoDecomposition {
    std::size_t canonical_index=0, source_orbital_index=0; // zero-based; distinct identities
    NboSpin spin=NboSpin::Total;
    bool available=false;
    std::string status="unavailable", detail;
    std::optional<double> occupation, weight_sum, normalization_error;
    std::optional<double> projection_residual_norm; // AO-metric norm of canonical minus local projection
    NboSource matrix_source;
    std::vector<NboNaoContribution> rows;
    std::vector<NboNaoGroupContribution> atoms, shells;
};
struct NboAssociation {
    bool compatible=false;
    std::string status="not_checked", detail;
    // Numerical equivalence cannot prove a producer job/step identity. External
    // production manifests supply that evidence; filenames never do.
    std::string provenance_status="producer_step_unverified";
    double geometry_max_error_bohr=0, overlap_max_error=0, density_max_error=0, canonical_max_error=0;
    // NBO AO row -> literal Gaussian AO row, with C_G = scale * C_NBO.
    std::vector<std::size_t> gaussian_row;
    std::vector<double> coefficient_scale;
    std::vector<NboCanonicalEvidence> canonical_evidence;
};
struct NboDataset {
    std::string producer_version;
    std::vector<NboPopulation> populations;
    std::vector<NboNao> naos;
    std::vector<NboOrbital> orbitals;
    std::vector<NboE2> e2;
    std::vector<NboE2Section> e2_sections;
    std::vector<NboWiberg> wiberg;
    std::vector<NboNlmo> nlmos;
    std::vector<NboMatrix> matrices;
    std::vector<NboSource> cmo_summaries; // thresholded text is never a complete matrix
    std::vector<std::string> warnings;
    std::optional<NboArchive> archive;
    NboAssociation association;
    std::vector<NboNaoValidation> nao_validation;
    std::vector<NboMoDecomposition> mo_decompositions;
    NboSource source;
};
struct NboReadOptions {
    std::optional<std::size_t> analysis_segment; // zero-based; multiple analyses require explicit selection
    std::filesystem::path archive47, aonbo, nbomo, naomo, aonao, naonbo;
};
NboDataset read_nbo(const std::filesystem::path& output, const NboReadOptions& options={});
NboArchive read_nbo_archive(const std::filesystem::path& path);
// Explicit complete W-format matrix reader used by capability-isolated import.
std::vector<NboMatrix> read_nbo_matrix(const std::filesystem::path& path,const std::string& kind,std::size_t basis_count,bool open_shell);
std::vector<NboMatrix> read_nbo_matrix_rectangular(const std::filesystem::path& path,const std::string& kind,std::size_t rows,std::size_t columns,bool open_shell);
NboAssociation associate_nbo(NboDataset& dataset, const Wavefunction& canonical);
Wavefunction make_nbo_wavefunction(const NboDataset& dataset, const Wavefunction& canonical);
std::string serialize_nbo_json(const NboDataset& dataset);
const NboMoDecomposition* nbo_mo_decomposition(const NboDataset& dataset, std::size_t canonical_index) noexcept;
const char* nbo_spin_name(NboSpin spin) noexcept;
} // namespace cov
