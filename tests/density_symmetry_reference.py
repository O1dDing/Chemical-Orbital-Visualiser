# SPDX-License-Identifier: CPAL-1.0
#
# The contents of this file are subject to the Common Public Attribution License
# Version 1.0 (the “License”); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
# License Version 1.1 but Sections 14 and 15 have been added to cover use of software
# over a computer network and provide for limited attribution for the Original
# Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
#
# Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
# WARRANTY OF ANY KIND, either express or implied. See the License for the specific
# language governing rights and limitations under the License.
#
# The Original Code is Chemical Orbital Visualiser (COV).
#
# The Original Developer is O1dDing, who is also the Initial Developer.
#
# The Initial Developer of the Original Code is O1dDing. All portions of the code
# written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
#
# Contributor O1dDing.

"""Independent density-operator and MO-subspace symmetry diagnostics.

The Hilbert-Schmidt operator norm is not the real-space density L2 norm.
No scientific pass or irrep assignment is made by this module. Signed squared
differences, projection errors and arithmetic estimates remain explicit.
"""
from __future__ import annotations

from wavefunction_reference import np


def symmetric_roots(metric):
    metric = np.asarray(metric, dtype=float)
    if (metric.ndim != 2 or metric.shape[0] != metric.shape[1] or not len(metric)
            or not np.isfinite(metric).all()
            or not np.allclose(metric, metric.T, rtol=0, atol=1e-10)):
        raise ValueError('Finite nonempty symmetric metric required')
    values, vectors = np.linalg.eigh((metric+metric.T)/2)
    if values[0] <= 1e-10:
        raise ValueError('The supplied subspace has no stable positive metric')
    root = (vectors*np.sqrt(values)) @ vectors.T
    inverse = (vectors/np.sqrt(values)) @ vectors.T
    return root, inverse


def prepare_spin_densities(alpha, beta, overlap, occ_alpha, occ_beta, *, shared_coefficients):
    """Express spin densities in orthonormal alpha MO coordinates.

    For rectangular spaces the beta residual is measured directly in the AO
    metric. Its density-error upper bound is retained, never hidden by fitting.
    Original coefficient arrays are not modified.
    """
    alpha, beta, overlap = (np.asarray(value, dtype=float) for value in (alpha, beta, overlap))
    fa, fb = np.asarray(occ_alpha, dtype=float), np.asarray(occ_beta, dtype=float)
    if (alpha.ndim != 2 or beta.ndim != 2 or alpha.shape[0] != beta.shape[0]
            or overlap.shape != (alpha.shape[0], alpha.shape[0])
            or fa.shape != (alpha.shape[1],) or fb.shape != (beta.shape[1],)
            or any(not np.isfinite(value).all() for value in (alpha, beta, overlap, fa, fb))
            or np.any(fa < 0) or np.any(fb < 0)):
        raise ValueError('Finite compatible coefficients, AO metric and nonnegative occupations required')
    if not np.allclose(overlap, overlap.T, rtol=0, atol=1e-10):
        raise ValueError('The independent AO metric must be symmetric')
    cholesky = np.linalg.cholesky((overlap+overlap.T)/2)
    gram = alpha.T @ overlap @ alpha
    root, inverse = symmetric_roots(gram)
    qa = alpha @ inverse
    ha = (root*fa) @ root
    if shared_coefficients:
        if not np.array_equal(alpha, beta):
            raise ValueError('Shared spin coefficient provenance disagrees with the arrays')
        hb = (root*fb) @ root
        total = (root*(fa+fb)) @ root
        spin = (root*(fa-fb)) @ root
        change = np.eye(alpha.shape[1])
        relative_residual = 0.0
        bound = 0.0
        exact_zero_spin = bool(np.array_equal(fa, fb))
    else:
        beta_in_qa = qa.T @ overlap @ beta
        projected_beta = qa @ beta_in_qa
        residual = cholesky.T @ (beta-projected_beta)
        norms = np.linalg.norm(cholesky.T @ beta, axis=0)
        if np.any(norms <= 0):
            raise ValueError('A beta MO has zero AO-metric norm')
        relative_residual = float(np.max(np.linalg.norm(residual, axis=0)/norms))
        a = np.linalg.norm((cholesky.T @ projected_beta)*np.sqrt(fb))
        b = np.linalg.norm(residual*np.sqrt(fb))
        # ||B F E^T + E F B^T + E F E^T||HS <= 2 ||B sqrt(F)||F
        # ||E sqrt(F)||F + ||E sqrt(F)||F^2, in the independent AO metric.
        bound = float(2*a*b+b*b)
        hb = (beta_in_qa*fb) @ beta_in_qa.T
        total, spin = ha+hb, ha-hb
        change = inverse @ beta_in_qa
        exact_zero_spin = False
    return {'alpha_gram': gram, 'alpha_inverse_sqrt_gram': inverse,
            'beta_coefficients_in_alpha': change,
            'densities': {'alpha': ha, 'beta': hb, 'total': total, 'spin': spin},
            'beta_relative_column_residual': relative_residual,
            'beta_density_projection_error_bound_hs': bound,
            'exact_zero_spin_from_shared_occupations': exact_zero_spin,
            'density_scope': 'one-electron operators, not pointwise density norms'}


def density_operation_difference(source, target, operation, *, source_projection_bound=0., target_projection_bound=0.):
    """Measure ||U rho_source U^T - rho_target||HS, allowing subspace leakage.

    A = Q^T U Q may be a contraction. The projected direct difference is a lower
    bound for the represented source/target operators' full-space difference.
    The signed full-space expression avoids pretending negative cancellation
    is a measured physical zero. Arithmetic scales here are estimates, not
    rigorous floating-point interval bounds.
    """
    h, j, a = (np.asarray(value, dtype=float) for value in (source, target, operation))
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or j.shape != h.shape or a.shape != h.shape
            or not len(h) or any(not np.isfinite(value).all() for value in (h, j, a))
            or not np.allclose(h, h.T, rtol=0, atol=1e-10)
            or not np.allclose(j, j.T, rtol=0, atol=1e-10)):
        raise ValueError('Compatible finite symmetric densities and finite operation required')
    if min(source_projection_bound, target_projection_bound) < 0 or not np.isfinite(source_projection_bound+target_projection_bound):
        raise ValueError('Finite nonnegative projection error bounds required')
    projected = a @ h @ a.T
    source_norm2 = float(np.sum(h*h))
    target_norm2 = float(np.sum(j*j))
    cross = float(np.sum(j*projected))
    signed = source_norm2+target_norm2-2*cross
    direct = float(np.linalg.norm(projected-j))
    arithmetic = float(128*np.finfo(float).eps*len(h)*(
        source_norm2+target_norm2+2*abs(cross)+np.sum(projected*projected)))
    source_norm = float(np.sqrt(source_norm2))
    return {'source_hs_norm': source_norm, 'target_hs_norm': float(np.sqrt(target_norm2)),
            'full_space_hs_difference_squared_signed': signed,
            'arithmetic_squared_scale_estimate': arithmetic,
            'projected_direct_hs_difference': direct,
            'projected_direct_relative_difference': direct/source_norm if source_norm > 0 else None,
            'full_space_relative_norm_estimate': float(np.sqrt(max(signed, 0)))/source_norm if source_norm > 0 else None,
            'source_and_target_projection_error_bound_hs': source_projection_bound+target_projection_bound,
            'density_status': 'both_operators_exactly_zero' if not np.any(h) and not np.any(j) else 'measured',
            'strict_equality_claimed': False}


def energy_connected_blocks(energies, tolerance):
    energies = np.asarray(energies, dtype=float)
    if energies.ndim != 1 or not len(energies) or not np.isfinite(energies).all() or not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('Finite energies and nonnegative window required')
    order = np.argsort(energies, kind='stable')
    groups = [[int(order[0])]]
    for previous, current in zip(order[:-1], order[1:]):
        if energies[current]-energies[previous] <= tolerance:
            groups[-1].append(int(current))
        else:
            groups.append([int(current)])
    return [{'members_zero_based': indices,
             'energy_range_hartree': [float(energies[indices].min()), float(energies[indices].max())],
             'grouping_scope': 'connected diagnostic energy window; not an exact-degeneracy assertion'} for indices in groups]


def subspace_measurement(raw_operator, metric):
    _, inverse = symmetric_roots(metric)
    raw_operator = np.asarray(raw_operator, dtype=float)
    if raw_operator.shape != inverse.shape or not np.isfinite(raw_operator).all():
        raise ValueError('Finite compatible raw operator matrix required')
    projected = inverse @ raw_operator @ inverse
    values = np.linalg.svd(projected, compute_uv=False)
    dimension = len(projected)
    return {'dimension': dimension, 'character': float(np.trace(projected)),
            'min_principal_overlap': float(values.min()), 'max_principal_overlap': float(values.max()),
            'subspace_leakage_squared_signed': float(dimension-np.sum(values**2)),
            'arithmetic_squared_scale_estimate': float(64*np.finfo(float).eps*(dimension+np.sum(values**2))),
            'projected_unitarity_max_error': float(np.max(np.abs(projected.T @ projected-np.eye(dimension)))),
            'exact_irrep_claimed': False}
