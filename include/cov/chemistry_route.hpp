#pragma once

#include "cov/interaction_graph.hpp"
#include "cov/nbo_integration.hpp"
#include "cov/pi_coupling.hpp"

#include <optional>
#include <string>
#include <vector>

namespace cov {

// A consumer asks for one capability at a time. A present report is evidence,
// not a licence to use values whose independent validation failed.
enum class RoutedStatus { Available, NotAnalysed, Insufficient, Rejected, Unsupported, NotApplicable,
                          NotReportedAboveThreshold };
enum class RoutedProvider { None, Legacy, Nbo };

const char* routed_status_name(RoutedStatus) noexcept;
const char* routed_provider_name(RoutedProvider) noexcept;

template<class T> struct RoutedResult {
    RoutedStatus status=RoutedStatus::NotAnalysed;
    RoutedProvider provider=RoutedProvider::None;
    std::string method, reason, fallback_reason;
    std::vector<NboSource> evidence;
    std::optional<T> value;
    [[nodiscard]] bool available() const noexcept {
        return status==RoutedStatus::Available && value.has_value();
    }
};

struct RoutedMoComposition {
    std::size_t canonical_index=0;
    std::vector<NboNaoContribution> rows;
    std::vector<NboNaoGroupContribution> atoms, shells;
    std::optional<double> retained_norm, residual_norm;
    std::optional<double> legacy_coverage, legacy_unresolved_fraction;
    bool complete=false;
    std::string semantics;
};
struct RoutedMoShellBalance {
    double core=0, valence=0, rydberg=0, other=0;
    bool complete=false;
};

// Producer-local E2/NLMO/multicentre facts mapped by actual transformation
// columns to a canonical MO. These are overlap relationships: a canonical MO
// is never re-labelled a localized donor or assigned a fraction of E(2).
struct RoutedLocalRelation {
    std::string id, kind, label, semantics;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> atoms;
    std::vector<NboOrbitalRef> orbitals;
    std::optional<double> printed_value;
    std::string units;
    NboSource source;
};
struct RoutedRelationProjection {
    std::string relation_id;
    std::vector<std::optional<double>> canonical_projection_weights;
};
// One source orbital's measured projections, shared by all local records that
// name it. Missing transforms remain null; a measured zero remains zero.
struct RoutedSourceOrbitalProjection {
    NboOrbitalRef orbital;
    std::vector<std::optional<double>> canonical_weights;
};

struct RoutedAnalysis {
    std::string canonical_fingerprint, integration_id;
    RoutedResult<std::vector<double>> total_atomic_charge;
    RoutedResult<std::vector<double>> atomic_spin;
    std::vector<RoutedResult<RoutedMoComposition>> mo_composition;
    std::vector<RoutedLocalRelation> local_relation_registry;
    std::vector<RoutedSourceOrbitalProjection> local_source_projections;
    std::vector<std::size_t> mo_relation_counts;
    bool local_relations_normalized=false;
    // Per-MO status metadata. In normalized routes the available value is an
    // empty placeholder; materialize the selected MO with routed_mo_relations.
    std::vector<RoutedResult<std::vector<RoutedRelationProjection>>> mo_relations;
    std::vector<RoutedResult<NboPiCoupling>> pi_couplings;
    RoutedResult<InteractionGraph> interaction_graph;
};

// Rebuild this value on canonical load and on every NBO reattachment. It owns
// its routed values and never modifies the canonical wavefunction or producer
// integration. A null integration means the independent legacy provider only.
RoutedAnalysis route_chemistry(const Wavefunction& canonical,
                               const NboIntegration* integration=nullptr);
inline std::optional<RoutedMoShellBalance> routed_mo_shell_balance(
    const RoutedAnalysis& routed,std::size_t index){
    if(index>=routed.mo_composition.size())return std::nullopt;
    const auto& result=routed.mo_composition[index];
    if(!result.available()||result.provider!=RoutedProvider::Nbo)return std::nullopt;
    RoutedMoShellBalance balance;balance.complete=result.value->complete;
    for(const auto& row:result.value->rows){
        if(row.type.rfind("Cor",0)==0)balance.core+=row.weight;
        else if(row.type.rfind("Val",0)==0)balance.valence+=row.weight;
        else if(row.type.rfind("Ryd",0)==0)balance.rydberg+=row.weight;
        else balance.other+=row.weight;
    }
    return balance;
}
RoutedResult<std::vector<RoutedRelationProjection>> routed_mo_relations(
    const RoutedAnalysis&,std::size_t canonical_index);
// v2 is normalized by default. The explicit expanded mode reproduces the v1
// relationship payload for consumers migrating from that schema.
std::string serialize_routed_analysis_json(const RoutedAnalysis&,bool include_projector_matrices=false,
    bool expand_relation_projections=false,
    const std::optional<std::vector<std::size_t>>& canonical_scope=std::nullopt);

// The scene is a typed identity, including zero-based canonical MO index zero.
// NBO inspection carries its exact signed selection rather than an NBO-set bool.
enum class ActiveOrbitalKind { Canonical, NboSet, Inspection };
struct ActiveOrbitalView {
    ActiveOrbitalKind kind=ActiveOrbitalKind::Canonical;
    // Inspection representation is independent of the orbital kind of its
    // first term: a SALC may be a signed NAO combination.
    std::string source_id, label, source_label, display_name_evidence;
    // Produced only by the typed naming serializer; keeps display scope and
    // source-subspace evidence without making analysis depend on the UI model.
    std::string display_name_metadata_json;
    std::string semantic_kind, group_id;
    NboSpin spin=NboSpin::Total;
    NboSpin source_spin=NboSpin::Total;
    std::string spin_semantics;
    std::optional<std::size_t> canonical_index, rendered_index;
    std::optional<NboOrbitalSelection> selection;
};
const char* active_orbital_kind_name(ActiveOrbitalKind) noexcept;
std::string serialize_active_orbital_view_json(const ActiveOrbitalView&);

} // namespace cov
