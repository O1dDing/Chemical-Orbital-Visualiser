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

"""Independent finite-point electron-density evaluation from original FCHK MOs.

This module uses IOData/GBasis conventions without AO or MO fitting. Its point
statistics are neither continuous-space L2 integrals nor COV pass decisions.
"""
from __future__ import annotations

from wavefunction_reference import np
from gbasis.evals.eval import evaluate_basis


def density_points(coordinates, grid_shape=(13, 13, 13), padding=4.0, atom_radii=(0.1, 0.5, 1.5)):
    """Fixed Cartesian cloud plus explicit atom-near strata, all in bohr.

    Point identities, including coincidences, are retained. No quadrature
    weights are implied. Under a rigid covariance check, transform this very
    cloud along with the molecule instead of silently regenerating it.
    """
    xyz = np.asarray(coordinates, dtype=float)
    shape = np.asarray(grid_shape)
    radii = np.asarray(atom_radii, dtype=float)
    if (xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all()
            or shape.shape != (3,) or not np.issubdtype(shape.dtype, np.integer) or np.any(shape < 2)
            or not np.isfinite(padding) or padding <= 0 or radii.ndim != 1 or not len(radii)
            or not np.isfinite(radii).all() or np.any(radii <= 0)):
        raise ValueError('Finite coordinates, integer grid shape, and positive point-cloud lengths required')
    low, high = xyz.min(axis=0)-padding, xyz.max(axis=0)+padding
    fractional = np.indices(shape).reshape(3, -1).T/(shape-1)
    grid = low+fractional*(high-low)
    directions = np.concatenate((np.eye(3), -np.eye(3)), axis=0)
    blocks = [grid, xyz.copy()]
    strata = [{'name': 'cartesian_cloud', 'start': 0, 'stop': len(grid)},
              {'name': 'nuclear_centers', 'start': len(grid), 'stop': len(grid)+len(xyz)}]
    cursor = len(grid)+len(xyz)
    for radius in radii:
        block = (xyz[:, None, :]+radius*directions[None, :, :]).reshape(-1, 3)
        blocks.append(block)
        strata.append({'name': f'atom_offsets_{radius:g}_bohr', 'radius_bohr': float(radius),
                       'start': cursor, 'stop': cursor+len(block)})
        cursor += len(block)
    return np.concatenate(blocks), {'grid_shape': shape.tolist(), 'grid_low_bohr': low.tolist(),
        'grid_high_bohr': high.tolist(), 'atom_radii_bohr': radii.tolist(), 'strata': strata,
        'point_count': cursor, 'unit': 'bohr', 'quadrature_weights': None,
        'point_order': 'Cartesian C order, z fastest; then nuclei; then each radius, atom order, +x,+y,+z,-x,-y,-z',
        'scope': 'Explicit finite samples, with repeated coordinates retained as distinct sample identities'}


def spin_densities_at_points(basis, points, alpha, beta, occ_alpha, occ_beta, *, shared_coefficients=False, chunk_points=512):
    """Evaluate sum_i f_i |phi_i(r)|^2 separately for both stored spin states."""
    points, ca, cb, fa, fb = (np.asarray(value, dtype=float)
                            for value in (points, alpha, beta, occ_alpha, occ_beta))
    if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
            or ca.ndim != 2 or cb.ndim != 2 or ca.shape[0] != cb.shape[0] or not ca.shape[0]
            or fa.shape != (ca.shape[1],) or fb.shape != (cb.shape[1],)
            or any(not np.isfinite(x).all() for x in (points, ca, cb, fa, fb))
            or np.any(fa < 0) or np.any(fb < 0) or not isinstance(chunk_points, int) or chunk_points < 1):
        raise ValueError('Compatible finite coefficients, samples and nonnegative occupations required')
    if shared_coefficients and not np.array_equal(ca, cb):
        raise ValueError('Shared coefficient provenance disagrees with the arrays')
    occupied_a, occupied_b = fa > 0, fb > 0
    result = np.zeros((2, len(points)), dtype=float)
    for start in range(0, len(points), chunk_points):
        stop = min(start+chunk_points, len(points))
        ao = evaluate_basis(basis, points[start:stop])
        if ao.shape != (ca.shape[0], stop-start) or not np.isfinite(ao).all():
            raise ValueError('AO evaluation has an unexpected shape or a nonfinite value')
        if shared_coefficients:
            occupied = occupied_a | occupied_b
            phi2 = (ca[:, occupied].T @ ao)**2
            result[0, start:stop] = fa[occupied] @ phi2
            result[1, start:stop] = fb[occupied] @ phi2
        else:
            if occupied_a.any():
                result[0, start:stop] = fa[occupied_a] @ ((ca[:, occupied_a].T @ ao)**2)
            if occupied_b.any():
                result[1, start:stop] = fb[occupied_b] @ ((cb[:, occupied_b].T @ ao)**2)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ValueError('Electron density evaluation is nonfinite or negative')
    return result


def density_channels(alpha_beta):
    values = np.asarray(alpha_beta, dtype=float)
    if values.ndim != 2 or values.shape[0] != 2 or not values.shape[1] or not np.isfinite(values).all():
        raise ValueError('Two finite nonempty spin sample arrays required')
    return {'alpha': values[0], 'beta': values[1], 'total': values[0]+values[1], 'spin': values[0]-values[1]}


def sampled_difference(source, transformed):
    """Point-vector diagnostics with explicit zero/underflow scope and units."""
    a, b = np.asarray(source, dtype=float), np.asarray(transformed, dtype=float)
    if a.ndim != 1 or b.shape != a.shape or not len(a) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Equal finite nonempty sampled arrays required')
    delta = b-a
    norm, delta_norm = float(np.linalg.norm(a)), float(np.linalg.norm(delta))
    peak = float(np.max(np.abs(a)))
    return {'point_count': len(a), 'source_sample_norm': norm,
        'sample_difference_norm': delta_norm, 'relative_sample_norm_difference': delta_norm/norm if norm > 0 else None,
        'maximum_absolute_density_difference': float(np.max(np.abs(delta))),
        'maximum_absolute_density_difference_unit': 'electrons / bohr^3',
        'relative_peak_difference': float(np.max(np.abs(delta)))/peak if peak > 0 else None,
        'source_samples_exactly_zero': bool(not np.any(a)), 'transformed_samples_exactly_zero': bool(not np.any(b)),
        'maximum_difference_point_index': int(np.argmax(np.abs(delta))),
        'norm_scope': 'unweighted finite point vector; no continuous-space integral or all-space zero assertion',
        'formal_case_pass': False}
