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

"""Independent nuclear isometries and AO-metric operation measurements.

No COV labels, orbital energies, molecule names, or nearest geometry templates
enter the nuclear search. Reported groups always carry an explicit coordinate
tolerance and closure residual; they are not declarations of exact symmetry.
Row coordinates use x -> (x-origin) @ operation + origin throughout.
"""
from __future__ import annotations

from wavefunction_reference import (block_diag, inverse_sqrt_positive, np,
                                    overlap_integral, shell_rotation)
from scipy.optimize import linear_sum_assignment


def checked_operation(operation):
    matrix = np.asarray(operation, dtype=float)
    if (matrix.shape != (3, 3) or not np.isfinite(matrix).all()
            or not np.allclose(matrix.T @ matrix, np.eye(3), rtol=0, atol=1e-11)):
        raise ValueError("A finite orthogonal 3 by 3 operation is required")
    return matrix


def checked_geometry(coordinates, labels):
    xyz = np.asarray(coordinates, dtype=float)
    labels = tuple(labels)
    if (xyz.ndim != 2 or xyz.shape != (len(labels), 3) or not len(xyz)
            or not np.isfinite(xyz).all()):
        raise ValueError("Nonempty finite coordinates and one label per atom are required")
    # Hashable labels may include effective nuclear charge or basis identity.
    for label in labels:
        hash(label)
    return xyz, labels


def match_operation(coordinates, labels, operation, origin=None):
    """Element/identity-preserving bijection, never many-to-one nearest atoms."""
    xyz, labels = checked_geometry(coordinates, labels)
    matrix = checked_operation(operation)
    origin = xyz.mean(axis=0) if origin is None else np.asarray(origin, dtype=float)
    if origin.shape != (3,) or not np.isfinite(origin).all():
        raise ValueError("A finite three-dimensional origin is required")
    moved = (xyz-origin) @ matrix+origin
    permutation = np.full(len(xyz), -1, dtype=int)
    for label in dict.fromkeys(labels):
        indices = np.array([i for i, value in enumerate(labels) if value == label])
        distances = np.linalg.norm(moved[indices, None, :]-xyz[None, indices, :], axis=2)
        rows, cols = linear_sum_assignment(distances**2)
        permutation[indices[rows]] = indices[cols]
    if sorted(permutation.tolist()) != list(range(len(xyz))):
        raise ValueError("An atom bijection was not obtained")
    residual = np.linalg.norm(moved-xyz[permutation], axis=1)
    return {"permutation_zero_based": permutation.tolist(),
            "max_displacement_bohr": float(residual.max()),
            "rms_displacement_bohr": float(np.sqrt(np.mean(residual**2)))}


def _frame(first, second):
    a = first/np.linalg.norm(first)
    normal = np.cross(a, second/np.linalg.norm(second))
    normal /= np.linalg.norm(normal)
    # Subtracting nearly parallel vectors loses orthogonality near linearity.
    # Cross products retain an orthogonal frame without loosening the gate.
    b = np.cross(normal, a)
    b /= np.linalg.norm(b)
    return np.column_stack((a, b, np.cross(a, b)))


def group_closure(operations, matrix_tolerance):
    matrices = np.asarray(operations, dtype=float)
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3) or not len(matrices):
        raise ValueError("A nonempty operation set is required")
    if not np.isfinite(matrix_tolerance) or matrix_tolerance <= 0:
        raise ValueError("A positive finite matrix tolerance is required")
    for matrix in matrices:
        checked_operation(matrix)
    maximum = 0.0
    table = []
    for left in matrices:
        products = left @ matrices
        distances = np.linalg.norm(products[:, None, :, :]-matrices[None, :, :, :], axis=(2, 3))
        nearest = np.argmin(distances, axis=1)
        maximum = max(maximum, float(distances[np.arange(len(matrices)), nearest].max()))
        table.append(nearest.tolist())
    identity_error = float(np.linalg.norm(matrices-np.eye(3), axis=(1, 2)).min())
    inverse_error = max(float(np.linalg.norm(matrices-matrix.T, axis=(1, 2)).min())
                        for matrix in matrices)
    closed = max(maximum, identity_error, inverse_error) <= matrix_tolerance
    return {"closed_under_tolerance": bool(closed), "matrix_tolerance": matrix_tolerance,
            "max_product_residual": maximum, "identity_residual": identity_error,
            "max_inverse_residual": inverse_error,
            "multiplication_table": table if closed else None}


def _rotation_order(matrix, tolerance, maximum=512):
    product = np.eye(3)
    for order in range(1, maximum+1):
        product = product @ matrix
        if np.linalg.norm(product-np.eye(3)) <= tolerance:
            return order
    return None


def finite_group_name(operations, closure):
    """Classify a closed O(3) group from its operations, with no geometry fit."""
    if not closure["closed_under_tolerance"]:
        return None
    matrices = np.asarray(operations)
    tolerance = closure["matrix_tolerance"]*4
    proper = [m for m in matrices if np.linalg.det(m) > 0]
    improper = [m for m in matrices if np.linalg.det(m) < 0]
    orders = [_rotation_order(m, tolerance) for m in proper]
    if not orders or None in orders:
        return None
    size, n = len(proper), max(orders)
    inversion = any(np.linalg.norm(m+np.eye(3)) <= tolerance for m in improper)
    reflection = [m for m in improper if abs(np.trace(m)-1) <= tolerance
                  and np.linalg.norm(m @ m-np.eye(3)) <= tolerance]
    if size == 1:
        return "C1" if not improper else ("Ci" if inversion else "Cs" if reflection else None)
    if size == n:
        family = "C"
    elif size == 2*n:
        family = "D"
    elif (size, n) in ((12, 3), (24, 4), (60, 5)):
        family = {12: "T", 24: "O", 60: "I"}[size]
    else:
        return None
    if not improper:
        return family+str(n) if family in ("C", "D") else family
    if len(improper) != size:
        return None
    if family in ("T", "O", "I"):
        if family == "T":
            return "Th" if inversion else "Td"
        return family+"h" if inversion else None
    # A horizontal plane has its normal along a highest-order proper axis.
    horizontal = False
    for matrix, order in zip(proper, orders):
        if order != n:
            continue
        _, _, vt = np.linalg.svd(matrix-np.eye(3))
        axis = vt[-1]
        sigma_h = np.eye(3)-2*np.outer(axis, axis)
        horizontal |= any(np.linalg.norm(m-sigma_h) <= tolerance for m in reflection)
    if family == "D":
        return f"D{n}h" if horizontal else f"D{n}d"
    if horizontal:
        return f"C{n}h"
    return f"C{n}v" if reflection else f"S{2*n}"


def nuclear_isometries(coordinates, labels, *, geometry_tolerance, matrix_tolerance=1e-7):
    """Enumerate candidate finite isometries, or describe continuous freedom.

    Two noncollinear nuclear vectors determine an orthogonal map up to a
    normal sign. All identity-compatible target pairs are visited. Approximate
    operations are retained as observations even when they do not close.
    """
    xyz, labels = checked_geometry(coordinates, labels)
    if not np.isfinite(geometry_tolerance) or geometry_tolerance <= 0:
        raise ValueError("A positive finite geometry tolerance is required")
    if not np.isfinite(matrix_tolerance) or matrix_tolerance <= 0:
        raise ValueError("A positive finite matrix tolerance is required")
    origin = xyz.mean(axis=0)
    centered = xyz-origin
    _, singular, vt = np.linalg.svd(centered, full_matrices=True)
    linear_axis = vt[0]
    line_residual = np.linalg.norm(centered-np.outer(centered @ linear_axis, linear_axis), axis=1)
    radius = np.linalg.norm(centered, axis=1)
    result = {"coordinate_origin_bohr": origin.tolist(), "geometry_tolerance_bohr": geometry_tolerance,
              "matrix_tolerance": matrix_tolerance, "scope": "all_supplied_nuclei",
              "singular_values_bohr": singular.tolist(), "atom_count": len(xyz),
              "claim": "candidate symmetry under the stated tolerances; not exact physical symmetry"}
    if radius.max() <= geometry_tolerance:
        result.update({"kind": "continuous_spherical" if len(xyz) == 1 else "coincident_or_unresolved_nuclei",
                       "group_under_tolerance": "Kh" if len(xyz) == 1 else None,
                       "max_center_displacement_bohr": float(radius.max()), "operations": []})
        return result
    if line_residual.max() <= geometry_tolerance:
        inversion = match_operation(xyz, labels, -np.eye(3), origin)
        result.update({"kind": "continuous_linear", "axis": linear_axis.tolist(),
                       "max_line_displacement_bohr": float(line_residual.max()),
                       "inversion": inversion, "operations": [],
                       "group_under_tolerance": "Dinfh" if inversion["max_displacement_bohr"] <= geometry_tolerance else "Cinfv"})
        return result
    first = int(np.argmax(radius))
    perpendicular = centered-np.outer(centered @ (centered[first]/radius[first]), centered[first]/radius[first])
    second = int(np.argmax(np.linalg.norm(perpendicular, axis=1)))
    original_frame = _frame(centered[first], centered[second])
    anchor_distance = np.linalg.norm(centered[first]-centered[second])
    first_candidates = [i for i, label in enumerate(labels) if label == labels[first]
                        and abs(radius[i]-radius[first]) <= 2*geometry_tolerance]
    second_candidates = [i for i, label in enumerate(labels) if label == labels[second]
                         and abs(radius[i]-radius[second]) <= 2*geometry_tolerance]
    operations, observations = [np.eye(3)], [match_operation(xyz, labels, np.eye(3), origin)]
    # Deduplication is substantially tighter than the group closure threshold.
    deduplication_tolerance = min(matrix_tolerance*.01, 1e-10)
    for target_first in first_candidates:
        for target_second in second_candidates:
            if (target_first == target_second or
                    abs(np.linalg.norm(centered[target_first]-centered[target_second])-anchor_distance) > 4*geometry_tolerance):
                continue
            if np.linalg.norm(np.cross(centered[target_first], centered[target_second])) <= 1e-13*radius[first]*radius[second]:
                continue
            target_frame = _frame(centered[target_first], centered[target_second])
            for normal_sign in (1., -1.):
                operation = original_frame @ np.diag([1., 1., normal_sign]) @ target_frame.T
                if any(np.linalg.norm(operation-old) <= deduplication_tolerance for old in operations):
                    continue
                measurement = match_operation(xyz, labels, operation, origin)
                if measurement["max_displacement_bohr"] <= geometry_tolerance:
                    operations.append(operation)
                    observations.append(measurement)
    closure = group_closure(operations, matrix_tolerance)
    result.update({"kind": "finite_candidates", "group_under_tolerance": finite_group_name(operations, closure),
                   "operation_count": len(operations), "proper_operation_count": int(sum(np.linalg.det(m) > 0 for m in operations)),
                   "closure": closure, "anchor_atoms_zero_based": [first, second],
                   "operations": [{"matrix": m.tolist(), "determinant": float(np.linalg.det(m)), **obs}
                                  for m, obs in zip(operations, observations)]})
    return result


def pullback_basis(basis, operation, origin):
    """Construct A(x R) = A_moved(x) T, including improper operations."""
    matrix = checked_operation(operation)
    origin = np.asarray(origin, dtype=float)
    if origin.shape != (3,) or not np.isfinite(origin).all() or not len(basis):
        raise ValueError("Nonempty basis and finite origin required")
    moved, transforms = [], []
    for shell in basis:
        center = (shell.coord-origin) @ matrix.T+origin
        moved.append(type(shell)(shell.angmom, center, shell.coeffs.copy(), shell.exps.copy(),
                                 shell.coord_type, icenter=shell.icenter))
        transforms.append(shell_rotation(shell, matrix))
    return tuple(moved), block_diag(*transforms)


def ao_operation_integrals(basis, operation, origin):
    """Independent overlap S and pullback matrix K_ij = <AO_i|U_R AO_j>."""
    moved, transform = pullback_basis(basis, operation, origin)
    n = len(transform)
    combined = overlap_integral(tuple(basis)+moved)
    return combined[:n, :n], combined[:n, n:] @ transform


def subspace_operation(basis, coefficients, operation, origin):
    """Analytic AO-integral operator and leakage for a selected real subspace.

    Rectangular C and nonorthogonal Cartesian AOs are supported. No label is
    inferred from a trace alone; all singular values and closure errors remain.
    """
    overlap, pullback = ao_operation_integrals(basis, operation, origin)
    n = len(overlap)
    coefficients = np.asarray(coefficients, dtype=float)
    if (coefficients.ndim != 2 or coefficients.shape[0] != n or not coefficients.shape[1]
            or not np.isfinite(coefficients).all()):
        raise ValueError("Finite nonempty AO by subspace coefficient matrix required")
    gram = coefficients.T @ overlap @ coefficients
    normalizer = inverse_sqrt_positive(gram)
    raw = coefficients.T @ pullback @ coefficients
    projected = normalizer @ raw @ normalizer
    values = np.linalg.svd(projected, compute_uv=False)
    dimension = coefficients.shape[1]
    signed_loss = float(dimension-np.sum(values**2))
    return {"dimension": dimension, "projected_operator": projected.tolist(),
            "character": float(np.trace(projected)), "singular_values": values.tolist(),
            "subspace_leakage_squared_signed": signed_loss,
            "arithmetic_squared_scale": float(64*np.finfo(float).eps*(dimension+np.sum(values**2))),
            "projected_unitarity_max_error": float(np.max(np.abs(projected.T @ projected-np.eye(dimension)))),
            "input_gram_max_error": float(np.max(np.abs(gram-np.eye(dimension)))),
            "scope": "supplied AO subspace; no point-group or orbital-label acceptance implied"}
