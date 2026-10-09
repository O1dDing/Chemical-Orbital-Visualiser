#pragma once

#include "cov/model.hpp"
#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace cov {
struct AtomicOrbitalSymmetryOptions {
    double metric_tolerance = 2e-5;
    double squared_residual_tolerance = 4e-8;
    double degeneracy_tolerance_hartree = 2e-5;
    double relative_rank_tolerance = 1e-12;
};

struct AtomicAngularComponent {
    int angular_momentum = -1;
    std::string label;
    // Squared S norm / unchanged source squared S norm; never c_mu squared.
    double weight = 0;
    std::vector<double> ao_coefficients;
    std::vector<std::size_t> source_shell_indices;
};

struct AtomicOrbitalAssignment {
    std::size_t source_orbital_index = 0; // index in Wavefunction::orbitals
    bool angular_momentum_verified = false;
    int angular_momentum = -1;
    std::string label;
    int inversion_parity = 0; // +1 even, -1 odd; 0 unresolved
    // Deliberately never inferred from energies, electron counts, or ECP cores.
    int principal_n = 0;
    double source_squared_s_norm = std::numeric_limits<double>::quiet_NaN();
    double projection_residual_squared = std::numeric_limits<double>::quiet_NaN();
    bool decomposition_verified = false;
    double reconstruction_residual_squared = std::numeric_limits<double>::quiet_NaN();
    double weight_sum_error = std::numeric_limits<double>::quiet_NaN();
    std::vector<AtomicAngularComponent> components;
    std::string status;
    std::string detail;
    // A pure l column need not have its 2l+1 radial partners in the source set.
    bool containing_span_closed = false;
    bool partner_block_verified = false;
    std::vector<std::size_t> containing_members;
    std::size_t representation_multiplicity = 0;
    std::string partner_block_id;
    double generator_leakage_squared = std::numeric_limits<double>::quiet_NaN();
    double source_gram_error = std::numeric_limits<double>::quiet_NaN();
    std::string partner_status;
    // Order among verified complete copies, NOT the principal quantum number.
    // Zero when incomplete/mixed repeated copies or tied energies obstruct it.
    std::size_t radial_copy_ordinal = 0;
};

struct AtomicOrbitalSymmetryResult {
    bool applicable = false; // exactly one represented nuclear centre
    bool available = false; // analytic AO action and actual S pass their gates
    std::string symmetry_group; // SO(3), with inversion parity carried separately
    std::string status;
    std::string detail;
    std::size_t ao_metric_rank = 0;
    std::array<std::size_t,5> angular_reference_ranks{};
    double generator_metric_error = std::numeric_limits<double>::quiet_NaN();
    std::vector<AtomicOrbitalAssignment> orbitals;
};

// Read-only spatial analysis of actual scalar AO/MO coefficients. Cartesian
// degree L is resolved into L,L-2,..., not promoted wholesale to pure l=L.
// Three analytic infinitesimal rotation generators establish continuous
// closure. Spin blocks stay separate. This is neither a finite point-group
// assignment nor a many-electron atomic term or relativistic j assignment.
AtomicOrbitalSymmetryResult analyse_atomic_orbital_symmetry(
    const Wavefunction&, const AtomicOrbitalSymmetryOptions& = {});
} // namespace cov
