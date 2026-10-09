#pragma once
#include "cov/nbo_salc.hpp"

namespace cov {
struct NboRoCommonEnergyOrbital {
    std::size_t canonical_index=0, source_orbital_index=0; // both zero based
    Spin spin=Spin::Alpha;
    double source_energy_hartree=0, occupation=0;
    std::optional<double> alpha_energy_hartree, beta_energy_hartree, common_energy_hartree;
    bool available=false;
    std::string status="unavailable";
};
struct NboRoCommonEnergyModel {
    bool available=false;
    // Association proves numerical compatibility, not the original SCF job's
    // operator identity. No current manifest verifier certifies the latter.
    bool source_step_identity_verified=false;
    std::string status="unavailable", detail, dataset_id, canonical_fingerprint;
    std::string operator_semantics="associated_archive_spin_average_expectation";
    std::string provenance_status;
    NboRestrictedOpenShellEvidence restricted_open_shell;
    std::vector<NboSalcEnergyEvidence> spin_evidence;
    std::vector<NboRoCommonEnergyOrbital> orbitals;
    double canonical_linkage_error=0, metric_error=0, projection_residual=0;
    double hermiticity_error=0, density_error=0, maximum_offdiagonal_hartree=0;
};
// Evaluates <source MO|(F_alpha+F_beta)/2|source MO>. These are
// expectations, not re-canonicalized eigenvalues. Inputs are never modified.
NboRoCommonEnergyModel build_nbo_ro_common_energy(
    const Wavefunction&, const NboIntegration&, const NboSalcModel& raw);
std::string serialize_nbo_ro_common_energy_json(const NboRoCommonEnergyModel&);
} // namespace cov
