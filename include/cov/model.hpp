#pragma once

#include "cov/gaussian_ao_transform.hpp"
#include "cov/topology_graph.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace cov {

constexpr double kAngstromToBohr = 1.8897261254578281;

enum class Spin : std::uint8_t {
    Alpha = 0,
    Beta = 1,
};

enum class WavefunctionSource : std::uint8_t {
    Unknown = 0,
    Molden = 1,
    Fchk = 2,
};

enum class DataProvenance : std::uint8_t {
    Unavailable = 0,
    Producer = 1,
    Derived = 2,
};

// One literal Gaussian point-group record. An absent field has no record;
// an empty/invalid printed field has a record with valid=false. Keep rejected
// text and Link1 identity without promoting it to producer symmetry metadata.
struct PointGroupSourceRecord {
    std::string field;
    std::string raw_record;
    std::string value;
    std::size_t line_number = 0;
    std::size_t job_segment = 0;
    bool valid = false;
    std::string field_text;
    std::string status;
};

struct OrbitalSymmetrySourceRecord {
    std::string source_path;
    std::size_t line_begin = 0;
    std::size_t line_end = 0;
    std::size_t job_segment = 0;
    std::string detected_group_context;
    std::string abelian_group_context;
    std::vector<std::size_t> orbital_indices;
    std::vector<std::string> labels;
};

// A derived irrep describes a specified MO subspace in a specified nuclear
// symmetry frame. It is independent of producer labels and local ligand-field
// interpretations. Axes are in the input Cartesian frame; they are absent
// unless axes_available is true. B1/B2 conventions can depend on this frame.
struct DerivedOrbitalSymmetryAssignment {
    std::string point_group;
    std::string label;
    std::vector<std::size_t> orbital_indices;
    double subspace_retention = std::numeric_limits<double>::quiet_NaN();
    double maximum_character_error = std::numeric_limits<double>::quiet_NaN();
    std::array<double,3> centre_bohr{};
    std::array<double,3> principal_axis{};
    std::array<double,3> secondary_axis{};
    bool axes_available = false;
    std::string axis_convention;
};

enum class NumericalStatus : std::uint8_t {
    NotComputed = 0,
    Available,
    MissingInput,
    InvalidInput,
    Failed,
};

// A format layout and an occupation interpretation are independent of MO
// display grouping. The shared integer interpretation is an explicit model
// assumption, not a producer statement about an unrestricted spin partition.
enum class OrbitalOccupationModel : std::uint8_t {
    Unspecified = 0,
    CanonicalShared,
    ExplicitSpin,
    SharedIntegerDeterminant,
    SharedFractionalUnresolved,
};

struct DensityMatrixDiagnostics {
    NumericalStatus status = NumericalStatus::NotComputed;
    DataProvenance input_provenance = DataProvenance::Unavailable;
    OrbitalOccupationModel occupation_model = OrbitalOccupationModel::Unspecified;
    DataProvenance model_provenance = DataProvenance::Unavailable;
    bool occupied_space_complete = false; // in the declared input model
    double occupation_sum = std::numeric_limits<double>::quiet_NaN();
    double expected_electron_count = std::numeric_limits<double>::quiet_NaN();
    double metric_trace = std::numeric_limits<double>::quiet_NaN();
    std::string detail;
};

struct DensityMatrixResult {
    std::vector<double> packed;
    DataProvenance provenance = DataProvenance::Unavailable;
    DensityMatrixDiagnostics diagnostics;
    bool available() const noexcept {
        return diagnostics.status == NumericalStatus::Available;
    }
};

struct DensityReconstruction {
    DensityMatrixResult total;
    DensityMatrixResult spin;
};

struct OrbitalMetricDiagnostics {
    Spin spin = Spin::Alpha;
    std::size_t orbital_count = 0;
    NumericalStatus status = NumericalStatus::NotComputed;
    double maximum_diagonal_error = std::numeric_limits<double>::quiet_NaN();
    double maximum_off_diagonal_error = std::numeric_limits<double>::quiet_NaN();
    double orthonormality_tolerance = 5.0e-6;
    std::string detail;
};

struct AoMetricDiagnostics {
    NumericalStatus status = NumericalStatus::NotComputed;
    std::size_t dimension = 0;
    bool finite = false;
    bool positive_semidefinite = false;
    double symmetry_error = std::numeric_limits<double>::quiet_NaN();
    double minimum_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    double maximum_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    double condition_number = std::numeric_limits<double>::quiet_NaN();
    double eigen_residual = std::numeric_limits<double>::quiet_NaN();
    double numerical_rank_tolerance = std::numeric_limits<double>::quiet_NaN();
    std::size_t numerical_rank = 0;
    std::vector<OrbitalMetricDiagnostics> orbital_blocks;
    std::string detail;
};

// Optional producer diagnostics are deliberately tri-state.  Absence means
// that no reliable statement was present in the input; it must not be
// presented as either success or failure.
enum class ScfConvergenceStatus : std::uint8_t {
    Unavailable = 0,
    Converged,
    Failed,
};

enum class WavefunctionStabilityStatus : std::uint8_t {
    Unavailable = 0,
    Stable,
    Unstable,
};

struct Atom {
    std::string symbol;
    int atomic_number = 0;
    double x = 0.0; // bohr
    double y = 0.0; // bohr
    double z = 0.0; // bohr

    // Producer-reported effective nuclear charge. This equals atomic_number
    // for all-electron calculations and is smaller when an ECP replaces core
    // electrons. Parsers set it explicitly; zero means "not reported".
    double nuclear_charge = 0.0;
};

struct Primitive {
    double exponent = 0.0;
    double coefficient = 0.0;
};

struct Shell {
    std::uint32_t atom_index = 0;
    std::uint32_t primitive_offset = 0;
    std::uint32_t primitive_count = 0;
    std::uint32_t basis_offset = 0;
    std::uint8_t angular_momentum = 0; // s=0, p=1, d=2, ...
    std::uint8_t pure = 0;             // 0 Cartesian, 1 real spherical
};

enum class ChemistryStatus : std::uint8_t {
    Determined = 0,
    Percentages,
    Undetermined,
    NotApplicable,
    Unavailable,
};

enum class OrbitalAngularFamily : std::uint8_t {
    Sigma = 0,
    Pi,
    Delta,
    Phi,
    Mixed,
    NotApplicable,
};

enum class OrbitalBondingRole : std::uint8_t {
    Bonding = 0,
    Antibonding,
    Nonbonding,
    Mixed,
    NotApplicable,
};

struct OrbitalChannelDistribution {
    double sigma = 0.0;
    double pi = 0.0;
    double delta = 0.0;
    double phi = 0.0;
    double undetermined = 1.0;
    OrbitalAngularFamily dominant = OrbitalAngularFamily::Mixed;
    ChemistryStatus status = ChemistryStatus::Unavailable;
};

struct OrbitalBondingDistribution {
    double bonding = 0.0;
    double antibonding = 0.0;
    double nonbonding = 0.0;
    double undetermined = 1.0;
    OrbitalBondingRole dominant = OrbitalBondingRole::Mixed;
    ChemistryStatus status = ChemistryStatus::Unavailable;
};

struct OrbitalAOContribution {
    std::uint32_t atom_index = 0;
    int principal_n = 0;
    int angular_momentum = 0;
    std::string label;
    double weight = 0.0;
};

struct OrbitalPairInteraction {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    std::string atom_a_label;
    std::string atom_b_label;

    // Total density-level Mayer index for the pair.
    double total_mayer_index = 0.0;

    // Per-unit-occupation Mulliken overlap character for this canonical MO.
    // This is representation/basis dependent and is not a universal bond order.
    double overlap_character = 0.0;
    double occupied_overlap_contribution = 0.0;

    OrbitalChannelDistribution channel;
    OrbitalBondingDistribution bonding;
};

struct OrbitalChemistry {
    bool available = false;
    bool valence_manifold = false;

    double deep_core_weight = 0.0;
    double semicore_weight = 0.0;
    double valence_weight = 0.0;
    double unresolved_weight = 1.0;

    OrbitalChannelDistribution channel;
    OrbitalBondingDistribution bonding;
    std::vector<OrbitalAOContribution> ao_contributions;
    std::vector<OrbitalPairInteraction> interactions;

    std::string multicentre_label;
    // True only when an actual non-unavailable MulticentreAssignment names
    // this MO. A generic label/population fallback remains a candidate.
    bool multicentre_assignment_available = false;
    std::string family_symbol;
    std::size_t multicentre_participating_atoms = 0;
    double multicentre_participating_electrons = 0.0;
    std::vector<std::uint32_t> multicentre_participating_atom_indices;
    std::size_t multicentre_channel_count = 0;
    std::string multicentre_source_subspace_id;
    double multicentre_source_electron_count = 0.0;
    double multicentre_confidence = 0.0;

    // A delocalised family is a property of a complete active subspace, not
    // of one canonical MO in isolation.  Repeating the stable identity and
    // complete membership on every member lets diagram/export code keep the
    // family intact even when only a filtered subset of levels is visible.
    std::size_t delocalised_participating_atoms = 0;
    double delocalised_participating_electrons = 0.0;
    std::vector<std::uint32_t> delocalised_participating_atom_indices;
    std::vector<std::uint32_t> delocalised_family_orbitals;
    std::string delocalised_family_id;
    std::size_t delocalised_orientation_channels = 0;
    bool delocalised_cyclic_topology = false;
    double delocalised_pi_confidence = 0.0;
    // Projection onto this member's oriented-p active family. This is not
    // the atom-pair-axis channel distribution above. NaN means unavailable.
    double delocalised_pi_weight = std::numeric_limits<double>::quiet_NaN();

    // FCHK-only canonical-MO evidence cannot always define a unique
    // donor/acceptor direction. In that case this remains "UND".
    std::string donor_acceptor = "UND";

    std::string method;
    double confidence = 0.0;
    std::string note;
};

struct MolecularOrbital {
    double energy_hartree = 0.0;
    double occupation = 0.0;
    Spin spin = Spin::Alpha;
    std::string symmetry;
    std::vector<double> coefficients;
    // Exact parsed Gaussian coefficients before the AO representation change.
    // Synthetic/Molden orbitals leave this Gaussian-specific field empty.
    std::vector<double> gaussian_source_coefficients;

    // Molden can carry occupation/symmetry explicitly; FCHK normally carries
    // electron counts rather than per-orbital occupations and usually does not
    // carry per-MO irreps. Keeping provenance prevents derived values from being
    // presented as producer data later in the UI/export layer.
    DataProvenance occupation_provenance = DataProvenance::Unavailable;
    DataProvenance spin_provenance = DataProvenance::Unavailable;
    std::string spin_source_text;
    // Index within the literal source spin block. Keeps canonical occupation
    // identity intact when an in-memory MO list is reordered.
    std::size_t source_orbital_index = std::numeric_limits<std::size_t>::max();
    DataProvenance symmetry_provenance = DataProvenance::Unavailable;

    // Derived after density/overlap/symmetry/bond analysis. Producer data is
    // never fabricated; unsupported conclusions remain percentage-valued or UND.
    OrbitalChemistry chemistry;
};

struct BondOrderRecord {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    double mayer_order = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

struct MulticentreCandidate {
    std::uint32_t orbital_index = 0;
    std::vector<std::uint32_t> atoms;
    std::vector<double> participation;
    double occupation = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

enum class MulticentreKind : std::uint8_t {
    Unclassified = 0,
    ThreeCentreTwoElectron,
    ThreeCentreFourElectron,
};

struct MulticentreAssignment {
    MulticentreKind kind = MulticentreKind::Unclassified;
    std::vector<std::uint32_t> atoms;
    std::vector<std::uint32_t> orbitals;
    double electron_count = 0.0;
    double confidence = 0.0;
    std::string rationale;
    DataProvenance provenance = DataProvenance::Unavailable;
    // True only when the three-centre active space passed an explicit
    // geometry-directed valence-subspace construction.  Generic population
    // candidates remain useful annotations, but must not take over the
    // compact MO diagram merely because they contain three orbitals.
    bool geometry_qualified_framework = false;
    // Equivalent local channels may be obtained by a derived rotation inside
    // one shared canonical source span.  In that case every channel carries
    // the same non-empty id, the total source electron count, and a partition
    // fraction.  This makes shared orbitals explicit and prevents downstream
    // consumers from summing their canonical occupations once per channel.
    std::string source_subspace_id;
    double source_subspace_electron_count = 0.0;
    double source_subspace_fraction = 1.0;
};

struct PiOrientationChannel {
    std::vector<std::uint32_t> atoms;
    std::array<double, 3> direction{0.0, 0.0, 1.0};
    double coherence = 0.0;
    bool cyclic = false;
};

enum class DelocalisedPiTopology : std::uint8_t {
    Unknown = 0,
    Path,
    Cycle,
    BranchedResonance,
    Spiro,
    HapticMetal,
    // Several disconnected but symmetry-equivalent pi subsystems whose
    // canonical MOs mix across the whole molecule. They are selected as one
    // direct-sum active space without claiming a fictitious covalent cycle.
    SymmetryDirectSum,
    // More than one electronic orientation channel without a proved,
    // channel-specific ring union. This does not assert molecular acyclicity.
    MultiChannel,
};

enum class PiTopologyGraphSource : std::uint8_t {
    Unavailable = 0,
    MayerDistanceModel,
    CovalentDistanceModel,
};

struct PiChannelRingWitness {
    std::array<std::size_t,2> channel_indices{};
    CycleUnionWitness ring_union;
};

// Structural graph and electronic-channel correspondence are distinct claims.
// Atom indices refer to the complete input atom list, including non-pi atoms
// that close a skeletal ring. These are model-selected edges, not input bonds.
struct PiTopologyGraphEvidence {
    PiTopologyGraphSource source = PiTopologyGraphSource::Unavailable;
    DataProvenance bond_order_provenance = DataProvenance::Unavailable;
    std::size_t atom_count = 0;
    double minimum_mayer_order = std::numeric_limits<double>::quiet_NaN();
    double maximum_covalent_radius_factor = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::pair<std::uint32_t,std::uint32_t>> edges;
    std::vector<std::uint32_t> examined_hubs;
    std::vector<std::uint32_t> unscoped_cycle_union_hubs;
    std::vector<PiChannelRingWitness> channel_ring_witnesses;
    bool channel_association_search_complete = false;
};

struct DelocalisedPiAssignment {
    std::string family_id;
    std::vector<std::uint32_t> atoms;
    std::vector<std::uint32_t> orbitals;
    double electron_count = 0.0;
    // A family can contain one planar channel, two degenerate perpendicular
    // channels (acetylene/O2/CO2), or several orthogonal channels meeting at a
    // shared sp centre (cumulenes). Individual canonical MOs are not assigned
    // to an arbitrary direction inside a degenerate subspace.
    std::vector<PiOrientationChannel> orientation_channels;
    // This topology is produced by the same pruned first-neighbour graph and
    // oriented-p analysis that generated the assignment.  Diagram policy must
    // consume it rather than rebuilding a second graph from raw all-pair
    // Mayer couplings (which contain legitimate through-bond terms).
    DelocalisedPiTopology topology = DelocalisedPiTopology::Unknown;
    bool cyclic_topology = false;
    PiTopologyGraphEvidence topology_graph;
    // For a unique branched centre, record the electronically occupied share
    // of its selected oriented-p column.  The normalisation is by that
    // column's selected canonical coverage, so 0..2 electrons has a stable
    // meaning across basis sizes.  It separates a strongly participating
    // resonance centre from a weak ligand-to-centre donor star without
    // element or molecule-name rules.
    double branch_centre_projected_occupation = 0.0;
    // Backward-compatible representative plane data. For a non-planar bundle,
    // plane_normal is only the first channel direction and plane_rms_bohr is 0;
    // consumers must inspect orientation_channels before making a plane claim.
    std::array<double, 3> plane_normal{0.0, 0.0, 1.0};
    double plane_rms_bohr = 0.0;
    double subspace_coverage = 0.0;
    double confidence = 0.0;
    std::string rationale;
    DataProvenance provenance = DataProvenance::Unavailable;
};

struct Wavefunction {
    std::vector<Atom> atoms;
    std::vector<Primitive> primitives;
    std::vector<Shell> shells;
    std::vector<MolecularOrbital> orbitals;
    std::uint32_t basis_count = 0;
    std::vector<GaussianAoTransformEntry> gaussian_ao_transform;
    bool pure_d = false;
    bool pure_f = false;
    bool pure_g = false;

    // Input/provenance metadata. These fields are intentionally orthogonal to
    // rendering so Molden and FCHK converge on the same Wavefunction model.
    WavefunctionSource source = WavefunctionSource::Unknown;
    std::uint32_t alpha_electrons = 0;
    std::uint32_t beta_electrons = 0;
    DataProvenance electron_counts_provenance = DataProvenance::Unavailable;
    OrbitalOccupationModel orbital_occupation_model = OrbitalOccupationModel::Unspecified;
    DataProvenance orbital_occupation_model_provenance = DataProvenance::Unavailable;
    // Gaussian FCHK carries these independently of alpha/beta counts.  Charge
    // zero is meaningful, so provenance rather than a magic value records
    // whether a producer actually supplied it.  Multiplicity zero remains the
    // harmless default for formats that do not report an electronic state.
    std::int32_t charge = 0;
    std::uint32_t multiplicity = 0;
    DataProvenance charge_provenance = DataProvenance::Unavailable;
    DataProvenance multiplicity_provenance = DataProvenance::Unavailable;
    // Optional producer atom-resolved population analysis. These values are
    // representation-dependent evidence for fragment/contact classification,
    // not formal oxidation states or integer atomic charges.
    std::vector<double> atomic_partial_charges;
    std::string atomic_partial_charge_scheme;
    DataProvenance atomic_partial_charge_provenance = DataProvenance::Unavailable;
    std::string source_title;
    std::string source_route;
    std::string enrichment_source;
    std::string point_group_detected;
    std::string point_group_used;
    DataProvenance point_group_provenance = DataProvenance::Unavailable;
    DataProvenance point_group_detected_provenance = DataProvenance::Unavailable;
    DataProvenance point_group_used_provenance = DataProvenance::Unavailable;
    std::vector<PointGroupSourceRecord> point_group_source_records;
    std::vector<OrbitalSymmetrySourceRecord> orbital_symmetry_source_records;
    std::vector<DerivedOrbitalSymmetryAssignment> derived_orbital_symmetry_assignments;

    // Optional Gaussian LOG/OUT diagnostics.  These fields are populated only
    // from explicit final producer statements; route keywords and inferred
    // spin values are never used as substitutes.
    ScfConvergenceStatus scf_convergence = ScfConvergenceStatus::Unavailable;
    DataProvenance scf_convergence_provenance = DataProvenance::Unavailable;
    WavefunctionStabilityStatus stability =
        WavefunctionStabilityStatus::Unavailable;
    DataProvenance stability_provenance = DataProvenance::Unavailable;
    std::string stability_detail;
    double spin_squared_before_annihilation = 0.0;
    double spin_squared_after_annihilation = 0.0;
    DataProvenance spin_squared_provenance = DataProvenance::Unavailable;

    // Packed lower-triangular AO density matrices when present or safely
    // reconstructable. Size is basis_count*(basis_count+1)/2.
    std::vector<double> total_density_packed;
    std::vector<double> spin_density_packed;
    DataProvenance total_density_provenance = DataProvenance::Unavailable;
    DataProvenance spin_density_provenance = DataProvenance::Unavailable;
    DensityMatrixDiagnostics total_density_diagnostics;
    DensityMatrixDiagnostics spin_density_diagnostics;
    // Preserve the actual producer matrices and the production MO reconstruction
    // separately. Serialize these once per load, never once per rendered frame.
    DensityReconstruction mo_density_reconstruction;

    // Analytic AO metric from the actual basis, independent of MO completeness.
    // The optional producer matrix remains separate for consistency checks.
    std::vector<double> ao_overlap;
    DataProvenance ao_overlap_provenance = DataProvenance::Unavailable;
    double ao_overlap_orthonormality_error = std::numeric_limits<double>::quiet_NaN();
    AoMetricDiagnostics ao_metric_diagnostics;
    std::vector<double> producer_ao_overlap;
    AoMetricDiagnostics producer_ao_metric_diagnostics;
    double producer_ao_overlap_basis_error = std::numeric_limits<double>::quiet_NaN();
    NumericalStatus producer_ao_overlap_basis_status = NumericalStatus::NotComputed;
    double producer_ao_overlap_basis_tolerance = 5.0e-6;

    // Provenance-aware analyses derived from density + overlap. Mayer bond
    // orders are pairwise evidence only. Multicentre assignments are made only
    // by the separate active-subspace classifier when its stricter evidence
    // requirements are met.
    std::vector<BondOrderRecord> bond_orders;
    DataProvenance bond_order_provenance = DataProvenance::Unavailable;
    std::vector<MulticentreCandidate> multicentre_candidates;
    std::vector<MulticentreAssignment> multicentre_assignments;
    std::vector<DelocalisedPiAssignment> delocalised_pi_assignments;
};

struct GridBox {
    float min_x = -5.0f;
    float min_y = -5.0f;
    float min_z = -5.0f;
    float max_x = 5.0f;
    float max_y = 5.0f;
    float max_z = 5.0f;
};

inline std::uint32_t cartesian_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>((l + 1u) * (l + 2u) / 2u);
}

inline std::uint32_t spherical_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>(2u * l + 1u);
}

inline std::uint32_t shell_basis_count(const Shell& shell) {
    return shell.pure ? spherical_basis_count(shell.angular_momentum)
                      : cartesian_basis_count(shell.angular_momentum);
}

inline const char* wavefunction_source_name(const WavefunctionSource source) noexcept {
    switch (source) {
        case WavefunctionSource::Molden: return "Molden";
        case WavefunctionSource::Fchk: return "FCHK";
        default: return "Unknown";
    }
}

inline const char* data_provenance_name(const DataProvenance provenance) noexcept {
    switch (provenance) {
        case DataProvenance::Producer: return "producer";
        case DataProvenance::Derived: return "derived";
        default: return "unavailable";
    }
}

inline const char* orbital_occupation_model_name(const OrbitalOccupationModel model) noexcept {
    switch (model) {
        case OrbitalOccupationModel::CanonicalShared: return "canonical-shared-source-identities";
        case OrbitalOccupationModel::ExplicitSpin: return "explicit-spin-orbital-occupations";
        case OrbitalOccupationModel::SharedIntegerDeterminant: return "shared-integer-determinant-assumption";
        case OrbitalOccupationModel::SharedFractionalUnresolved: return "shared-fractional-spin-unresolved";
        default: return "unspecified";
    }
}

inline const char* multicentre_kind_name(const MulticentreKind kind) noexcept {
    switch (kind) {
        case MulticentreKind::ThreeCentreTwoElectron: return "3c2e";
        case MulticentreKind::ThreeCentreFourElectron: return "3c4e";
        default: return "unclassified";
    }
}

} // namespace cov
