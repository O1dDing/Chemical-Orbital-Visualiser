#pragma once
#include "cov/nbo.hpp"
#include <map>

namespace cov {
enum class NboCapabilityState { Missing, Available, Rejected, Ambiguous, Unsupported };
struct NboCapability {
    std::string key;
    NboCapabilityState state=NboCapabilityState::Missing;
    std::string detail;
    std::vector<NboSource> sources;
    bool available() const noexcept { return state==NboCapabilityState::Available; }
};
enum class NboInputKind { Unknown, GaussianFchk, GaussianLog, NboReport, Archive, Matrix, Package };
struct NboIdentifiedInput {
    std::filesystem::path path;
    NboInputKind kind=NboInputKind::Unknown;
    std::string matrix_kind, detail;
    std::size_t analysis_count=0;
};
struct NboInputCandidate {
    std::string id, label;
    std::filesystem::path canonical, report, archive;
    std::optional<std::size_t> analysis_segment;
    std::map<std::string,std::filesystem::path> matrices;
    std::vector<std::string> diagnostics;
};
struct NboInputDiscovery {
    std::vector<NboIdentifiedInput> inputs;
    std::vector<NboInputCandidate> candidates;
    std::vector<std::string> diagnostics;
    bool selection_required=false;
};
// Discovery reads content only. It never runs producers or extracts executable
// packages. Multiple report/step candidates remain explicit user choices.
NboInputDiscovery discover_nbo_inputs(const std::vector<std::filesystem::path>& paths);

enum class NboOrbitalKind { Canonical, GaussianAO, NAO, PNAO, NHO, NBO, NLMO };
struct NboOrbitalRef {
    NboOrbitalKind kind=NboOrbitalKind::Canonical;
    NboSpin spin=NboSpin::Total;
    // Canonical: index in immutable Wavefunction::orbitals. GaussianAO:
    // literal Gaussian source AO index. Others: producer column/id minus one.
    std::size_t index=0;
};
bool operator==(const NboOrbitalRef& a,const NboOrbitalRef& b) noexcept;
const char* nbo_orbital_kind_name(NboOrbitalKind kind) noexcept;
const char* nbo_capability_state_name(NboCapabilityState state) noexcept;
// Centre-conditioned angular projections are independent evidence on each
// end of a local axis. Their fractions must never be summed as populations.
struct NboChannelEvidence {
    std::string status="insufficient_evidence", channel="unassigned", method, reason;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> axis_atoms;
    std::array<double,2> centre_fractions{};
    std::array<std::array<double,5>,2> conditional_angular_fractions{};
    double classification_threshold=0;
    double centre_threshold=0;
    NboSource source;
};
struct NboOrbitalDescriptor {
    NboOrbitalRef ref;
    std::string id, label, detail, energy_semantics;
    std::string display_label; // producer label remains unchanged in label
    std::vector<NboChannelEvidence> channels;
    std::vector<std::size_t> atoms; // ZERO-based canonical atom identities (raw dataset retains producer IDs)
    std::optional<double> occupation, energy_hartree;
    double metric_norm2=0;
    bool orthonormal_basis=false;
    std::vector<double> coefficients; // actual COV internal Gaussian AO convention
    NboSource source;
};
struct NboMoLink {
    NboOrbitalRef orbital;
    std::size_t canonical_index=0;
    double coefficient=0; // signed wavefunction coefficient, no display cutoff
    std::optional<double> weight; // squared coefficient only for orthonormal bases
    NboSource source;
};
enum class NboSelectionMode { Orbital, WeightedComponent, PartialSum, Combination, Overlay };
struct NboOrbitalTerm { NboOrbitalRef orbital; double coefficient=1; };
// Display derivation, separate from every producer orbital. A channel may map
// several original members onto one verified common spatial direction.
struct NboSpinSourceMember {
    std::string id,label;
    std::vector<NboOrbitalTerm> terms;
    std::optional<double> energy_hartree,occupation;
};
struct NboSpinChannelProjection {
    NboSpin spin=NboSpin::Total;
    std::vector<NboSpinSourceMember> members;
    std::vector<double> mapping; // common chi = sum mapping[i] * source member[i]
    std::vector<double> fock_matrix,density_matrix; // source-member basis, row-major
    std::optional<double> energy_hartree,occupation; // expectations on common chi
    std::string energy_status,occupation_status;
};
struct NboSpatialSpinInfo {
    std::string id,block_id,correspondence,detail,ro_method;
    std::size_t member_index=0,dimension=1;
    double tolerance=0,metric_error=0,spatial_residual=0;
    std::vector<NboSpinChannelProjection> channels;
    std::vector<NboSpinSourceMember> total_aliases; // provenance only, never counted again
    std::vector<double> total_alias_mapping; // same common chi in the Total alias basis
    double total_alias_residual=0,total_alias_metric_error=0;
    std::optional<double> energy_hartree,occupation;
    std::string energy_status,occupation_status;
};
std::string serialize_nbo_spatial_spin_json(const NboSpatialSpinInfo&);
struct NboOrbitalSelection {
    std::string dataset_id, label;
    std::string semantic_kind, group_id, source_id;
    NboSelectionMode mode=NboSelectionMode::Orbital;
    std::vector<NboOrbitalTerm> terms;
    std::optional<std::size_t> target_canonical_index;
    bool normalize=false; // raw components/partial sums retain their amplitudes
    std::optional<NboSpatialSpinInfo> spatial_spin;
};
std::string serialize_nbo_orbital_selection_json(const NboOrbitalSelection&);
struct NboStructureEvidence {
    std::string id, kind, label, detail, style;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> atoms; // ZERO-based canonical atoms; supports true multi-centre evidence
    std::vector<NboOrbitalRef> orbitals;
    std::optional<double> wiberg, value;
    std::optional<double> wiberg_alpha, wiberg_beta;
    std::string units, channel="unassigned"; // never infer sigma/pi/delta from ordinal
    std::vector<NboChannelEvidence> channels; // multiple associated local channels
    std::optional<unsigned> lewis_bond_count; // actual BD records; never rounded WBI
    NboSource source;
};
struct NboIntegration {
    std::string id, canonical_fingerprint;
    NboDataset dataset; // all successfully read raw evidence, including unavailable capabilities
    std::vector<NboCapability> capabilities;
    std::vector<NboOrbitalDescriptor> orbitals;
    std::vector<NboMoLink> links;
    std::vector<NboStructureEvidence> structure;
    std::vector<std::string> diagnostics;
};
// Layered association: bad advanced matrices disable their dependent capability
// without destroying canonical data or other independently verified capabilities.
NboIntegration integrate_nbo(const Wavefunction& canonical, const NboDataset& dataset);
NboIntegration read_nbo_integration(const Wavefunction& canonical, const NboInputCandidate& candidate);
std::string nbo_canonical_fingerprint(const Wavefunction& canonical);
const NboCapability* nbo_capability(const NboIntegration& data,const std::string& key) noexcept;
const NboOrbitalDescriptor* nbo_orbital(const NboIntegration& data,const NboOrbitalRef& ref) noexcept;
std::vector<NboMoLink> nbo_links_for_orbital(const NboIntegration& data,const NboOrbitalRef& ref);
std::vector<NboMoLink> nbo_links_for_mo(const NboIntegration& data,std::size_t canonical_index,NboOrbitalKind kind=NboOrbitalKind::NAO);
NboOrbitalSelection nbo_single_selection(const NboIntegration& data,const NboOrbitalRef& ref);
NboOrbitalSelection nbo_component_selection(const NboIntegration& data,const NboMoLink& link);
std::vector<NboOrbitalTerm> nbo_nlmo_components(const NboIntegration& data,const NboOrbitalRef& nlmo);
// Signed NAONHO column in the verified NAO basis; labels/atoms come from the
// referenced NAO descriptors. Sum coefficient squared by literal n/l to obtain
// hybrid composition without assigning shell identity from column order.
std::vector<NboOrbitalTerm> nbo_nho_components(const NboIntegration& data,const NboOrbitalRef& nho);
std::optional<NboOrbitalRef> nbo_nlmo_parent(const NboIntegration& data,const NboOrbitalRef& nlmo);
// Numerical fragment component of a selected MO; atoms are ZERO-based canonical
// identities. This is NOT an asserted SALC.
NboOrbitalSelection nbo_fragment_selection(const NboIntegration& data,std::size_t canonical_index,const std::vector<std::size_t>& atoms,NboOrbitalKind kind=NboOrbitalKind::NAO,bool normalize=false);
struct NboSelectionView {
    bool available=false;
    std::string status="unavailable", detail, label;
    Wavefunction wavefunction; // independent view; canonical is never mutated
    NboOrbitalSelection selection;
    std::vector<double> metric_norm2;
    std::vector<std::size_t> atoms;
    std::optional<double> reconstruction_error;
};
NboSelectionView make_nbo_selection_view(const NboIntegration& data,const Wavefunction& canonical,const NboOrbitalSelection& selection);
std::string serialize_nbo_integration_json(const NboIntegration& data);
std::string serialize_nbo_selection_json(const NboSelectionView& view);
} // namespace cov
