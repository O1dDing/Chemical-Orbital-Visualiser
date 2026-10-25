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

#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace cov {
enum class MetricSubspaceStatus { Available, MissingInput, InvalidInput, Failed };
struct MetricSubspaceDiagnostics {
    std::size_t supplied_columns=0;
    std::size_t numerical_rank=0;
    std::size_t unresolved_norm_columns=0;
    double rank_cutoff=0;
    double gram_inverse_residual=0;
    double orthonormality_residual=0;
    double retained_condition_number=0;
};
struct MetricSubspaceOverlap {
    MetricSubspaceStatus status=MetricSubspaceStatus::MissingInput;
    std::size_t basis_dimension=0;
    std::size_t metric_numerical_rank=0;
    bool basis_norm_preconditioned=false;
    // Spectral diagnostics below refer to the congruently transformed metric
    // with each individual AO normalized. Original S and coefficient data are
    // not changed. X and C undergo the same inverse basis transformation.
    double metric_rank_cutoff=0;
    double metric_symmetry_error=0;
    double metric_eigen_residual=0;
    double metric_minimum_eigenvalue=std::numeric_limits<double>::quiet_NaN();
    double subspace_overlap_trace=std::numeric_limits<double>::quiet_NaN();
    double mean_orbital_fraction=std::numeric_limits<double>::quiet_NaN();
    MetricSubspaceDiagnostics reference;
    MetricSubspaceDiagnostics orbitals;
    std::string detail;
};

class MetricSubspaceContext;
class MetricPreparedSubspace {
public:
    MetricPreparedSubspace()=default;
    MetricSubspaceDiagnostics diagnostics() const;
    MetricSubspaceStatus status() const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
    friend class MetricSubspaceContext;
};

// The context owns its normalized metric and source norm factors. Each space
// owns its prepared columns. Mutating a caller's input vector cannot affect a
// prepared result, and spaces from different contexts cannot be silently mixed.
class MetricSubspaceContext {
public:
    MetricSubspaceContext(std::size_t dimension,const std::vector<double>& metric);
    MetricSubspaceOverlap diagnostics() const;
    MetricPreparedSubspace prepare_space(std::size_t columns,const std::vector<double>& coefficients) const;
    MetricSubspaceOverlap compare(const MetricPreparedSubspace& reference,const MetricPreparedSubspace& orbitals) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
// Row-major S, X and C. The reported trace is tr(P_X P_C), and its mean is
// divided by the represented rank of C. It is invariant under well-conditioned
// changes of basis inside either space. The unweighted amount is not an
// electron population, nor an additive partition over overlapping centres.
// All numerical rank decisions and residuals are returned. A zero-rank target
// has trace zero and an undefined mean, which is preserved as NaN.
MetricSubspaceOverlap metric_subspace_overlap(
    std::size_t dimension,const std::vector<double>& metric,
    std::size_t reference_columns,const std::vector<double>& reference,
    std::size_t orbital_columns,const std::vector<double>& orbitals);
}
