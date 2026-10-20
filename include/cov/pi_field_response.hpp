#pragma once
#include "cov/pi_coupling.hpp"
#include <span>

namespace cov {
// These roles describe a frozen d-manifold ligand-field response, not charge
// transfer or an SCF stabilization energy. One shared RO response is returned.
enum class PiFieldRole { Unavailable, Negligible, DonorDominant, AcceptorDominant, Mixed };
const char* pi_field_role_name(PiFieldRole) noexcept;
struct PiFieldResponseOptions {
    double numerical_tolerance=1e-8;
    // Explicit display resolution, not a universal chemical threshold.
    double response_resolution_hartree=1e-4;
    double minimum_step_overlap=0.98;
    // Majority correspondence to the TRACKED spectral subspace, never an
    // atomic metal-composition threshold or a criterion for the field sign.
    double minimum_group_spectral_fraction=0.5;
    std::size_t initial_steps=16, maximum_steps=512;
};
struct PiFieldResponseInput {
    // All columns and operators are in ONE orthonormal reference. The full
    // Fock matrix is retained, including the environment outside both scopes.
    NboMatrix fock, metal_d_basis, ligand_pi_basis, total_density;
    NboMatrix canonical_columns;
    std::vector<std::size_t> canonical_indices; // column -> immutable global MO
    std::vector<std::size_t> centre_atoms,ligand_atoms;
    std::string canonical_fingerprint,operator_semantics,scope_semantics;
    double operator_error_hartree=0;
    bool shared_spatial_verified=false,density_verified=false;
    // Occupied same-spin metal -> ligand acceptor evidence is a separate
    // qualification; it never votes on the joined field role.
    bool occupied_metal_backbond_evidence=false;
    NboMatrix occupied_backbond_donor_basis; // optional complete common-space span
};
struct PiFieldTracking {
    bool verified=false;
    std::string reason;
    std::size_t steps=0,rank=0;
    double minimum_step_overlap=1,minimum_d_projection=1;
    double endpoint_eigen_residual_hartree=0;
    std::vector<double> reference_energies_hartree,final_energies_hartree;
};
struct PiFieldGroupAssessment {
    bool available=false,applicable=false;
    PiFieldRole role=PiFieldRole::Unavailable;
    std::string reason,method="shared-spatial-full-environment-frozen-pi-field-v1";
    std::vector<std::size_t> members;
    std::size_t metal_support_rank=0;
    double d_fraction=0,responsive_d_fraction=0;
    double tracked_d_fraction=0,tracked_responsive_fraction=0;
    double mean_shift_hartree=0,trace_shift_hartree=0;
    double response_min_hartree=0,response_max_hartree=0;
    double gross_response_hartree=0,numerical_error_bound_hartree=0;
    bool occupied_metal_backbond_supported=false;
};
struct PiFieldResponse {
    bool available=false,full_environment_retained=false,shared_spatial_verified=false;
    PiFieldRole role=PiFieldRole::Unavailable;
    std::string id,reason,canonical_fingerprint,operator_semantics,scope_semantics;
    std::vector<std::size_t> centre_atoms,ligand_atoms,canonical_indices;
    std::size_t full_dimension=0,metal_rank=0,ligand_rank=0;
    double numerical_error_bound_hartree=0,response_resolution_hartree=0;
    double minimum_group_spectral_fraction=0.5;
    double mean_shift_hartree=0,trace_shift_hartree=0;
    double gross_response_hartree=0,nonadditivity_norm_hartree=0;
    double occupied_d_electrons=0;
    bool density_verified=false,occupied_metal_backbond_supported=false;
    PiFieldTracking tracking;
    // Small d-reference matrices; invariant under joint change of d basis.
    NboMatrix d_response,d_occupied_pi_response,d_vacant_pi_response;
    // Original columns projected into the SAME d reference. These do not
    // imply a metal-character display cutoff or replace full MO composition.
    NboMatrix canonical_d_coordinates;
    NboMatrix canonical_tracked_d_coordinates;
    NboMatrix d_backbond_projector;
    std::vector<PiFieldGroupAssessment> groups;
};
struct PiFieldResponseAnalysis {
    std::string status="insufficient_evidence",reason;
    std::vector<PiFieldResponse> responses;
};
PiFieldResponse assess_pi_field_response(const PiFieldResponseInput&,
    const PiFieldResponseOptions& options={});
PiFieldGroupAssessment assess_pi_field_group(const PiFieldResponse&,
    const std::vector<std::size_t>& canonical_members);
// Reuses only complete geometric/localized projectors from the channel
// producer. Existing direction/ordinary-display labels never select a role.
PiFieldResponseAnalysis analyse_pi_field_responses(const Wavefunction&,
    const NboIntegration&,std::span<const NboPiCoupling* const> channels,
    const PiFieldResponseOptions& options={});
std::string pi_field_response_json(const PiFieldResponse&,bool include_matrices=false);
std::string pi_field_group_assessment_json(const PiFieldGroupAssessment&);
std::string pi_field_response_analysis_json(const PiFieldResponseAnalysis&,bool include_matrices=false);
}
