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

"""Independent, rotation-aware wavefunction comparison in the AO metric.

This module reads Gaussian FCHK through IOData and uses GBasis integrals. It
does not modify an input wavefunction and does not use any COV implementation.
Coordinate alignment determines AO rotations; orbitals never determine them.
All quantities are diagnostic evidence, not an automatic chemical verdict.
"""
from __future__ import annotations

import math
from importlib.metadata import version
import os
from pathlib import Path
import sys
import time

# Set these before NumPy/BLAS initialization. A containing Job Object supplies
# the aggregate CPU/memory ceiling when this module is used by a batch worker.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
REFERENCE_DEPS = Path(r"F:\Dev\cov-validation-20260905\reference-deps")
sys.path.insert(0, str(REFERENCE_DEPS))

import numpy as np
from gbasis.integrals.overlap import overlap_integral
from gbasis.spherical import generate_transformation
from gbasis.wrappers import from_iodata
from iodata import load_one
from scipy.linalg import block_diag
from scipy.optimize import linear_sum_assignment


def rigid_alignment(old, new):
    """Return R for row coordinates (old-old_center) @ R + new_center.

    Atom correspondence must already have been established by the input
    identity gate. No permutation/reflection search or orbital fit is done.
    Linear/atomic geometries have undetermined rotations, recorded explicitly.
    """
    old, new = np.asarray(old, float), np.asarray(new, float)
    if old.shape != new.shape or old.ndim != 2 or old.shape[1] != 3 or not len(old):
        raise ValueError("Equal nonempty atom-coordinate arrays are required")
    center_old, center_new = old.mean(axis=0), new.mean(axis=0)
    a, b = old-center_old, new-center_new
    u, values, vt = np.linalg.svd(a.T @ b)
    signs = np.eye(3)
    signs[-1, -1] = 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0
    rotation = u @ signs @ vt
    displacement = a @ rotation-b
    cutoff = max(values[0]*1e-12, 1e-20)
    rank = int(np.count_nonzero(values > cutoff))
    return {"rotation": rotation, "center_old": center_old, "center_new": center_new,
            "coordinate_rank": rank, "rotation_underdetermined": rank < 2,
            "rms_displacement_bohr": float(np.sqrt(np.mean(np.sum(displacement**2, axis=1)))),
            "max_displacement_bohr": float(np.max(np.linalg.norm(displacement, axis=1)))}


def angular_double_factorial(powers):
    value = 1
    for power in powers:
        for factor in range(1, 2*int(power), 2):
            value *= factor
    return value


def cartesian_rotation(old_order, new_order, rotation):
    """Exact homogeneous-polynomial rotation of normalized Cartesian AOs.

    For AO row evaluations A, A_new(x @ R) = A_old(x) @ T. The ratio
    sqrt(prod((2*a-1)!!)/prod((2*b-1)!!)) retains primitive normalization.
    """
    old_order, new_order = np.asarray(old_order, int), np.asarray(new_order, int)
    rotation = np.asarray(rotation, float)
    if rotation.shape != (3, 3) or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12):
        raise ValueError("AO transform requires an orthogonal coordinate rotation")
    lookup = {tuple(powers): index for index, powers in enumerate(old_order)}
    transform = np.zeros((len(old_order), len(new_order)))
    for column, powers in enumerate(new_order):
        polynomial = {(0, 0, 0): 1.0}
        for new_axis, power in enumerate(powers):
            for _ in range(int(power)):
                expanded = {}
                for monomial, coefficient in polynomial.items():
                    for old_axis in range(3):
                        target = list(monomial)
                        target[old_axis] += 1
                        target = tuple(target)
                        expanded[target] = expanded.get(target, 0.0)+coefficient*rotation[old_axis, new_axis]
                polynomial = expanded
        for monomial, coefficient in polynomial.items():
            if monomial not in lookup:
                raise ValueError("Cartesian component orders do not span the same angular momentum")
            transform[lookup[monomial], column] = coefficient*math.sqrt(
                angular_double_factorial(monomial)/angular_double_factorial(powers))
    return transform


def shell_rotation(shell, rotation):
    order = shell.angmom_components_cart
    cart = cartesian_rotation(order, order, rotation)
    if shell.coord_type == "cartesian":
        angular = cart
    else:
        pure = generate_transformation(shell.angmom, order, shell.angmom_components_sph, "right")
        angular, _, rank, _ = np.linalg.lstsq(pure, cart @ pure, rcond=None)
        if rank != pure.shape[1] or np.max(np.abs(pure @ angular-cart @ pure)) > 1e-11:
            raise ValueError("Pure angular space is not closed under the coordinate rotation")
    # GBasis component ordering is segmented contraction first, component next.
    return block_diag(*([angular]*shell.num_seg_cont))


def aligned_basis(basis, alignment):
    rotation = alignment["rotation"]
    result, transforms = [], []
    for shell in basis:
        coordinate = (shell.coord-alignment["center_new"]) @ rotation.T+alignment["center_old"]
        result.append(type(shell)(shell.angmom, coordinate, shell.coeffs.copy(), shell.exps.copy(),
                                  shell.coord_type, icenter=shell.icenter))
        transforms.append(shell_rotation(shell, rotation))
    return tuple(result), block_diag(*transforms)


def metric_summary(coefficients, overlap):
    gram = coefficients.T @ overlap @ coefficients
    residual = gram-np.eye(gram.shape[0])
    return gram, {"n_ao": coefficients.shape[0], "n_mo": coefficients.shape[1],
                  "max_orthonormality_error": float(np.max(np.abs(residual))),
                  "rms_orthonormality_error": float(np.sqrt(np.mean(residual**2)))}


def inverse_sqrt_positive(matrix):
    eigenvalues, vectors = np.linalg.eigh((matrix+matrix.T)/2)
    if not len(eigenvalues) or eigenvalues[0] <= 1e-10:
        raise ValueError("MO subspace has no stable full-rank metric")
    return (vectors/np.sqrt(eigenvalues)) @ vectors.T


def principal_overlap(cross, gram_old, gram_new):
    """Principal overlaps with metric normalization; raw metric errors remain recorded."""
    if not len(gram_old) or not len(gram_new):
        return {"old_rank": len(gram_old), "new_rank": len(gram_new),
                "singular_values": [], "status": "empty_subspace"}
    normalized = inverse_sqrt_positive(gram_old) @ cross @ inverse_sqrt_positive(gram_new)
    values = np.linalg.svd(normalized, compute_uv=False)
    return {"old_rank": len(gram_old), "new_rank": len(gram_new),
            "singular_values": values.tolist(), "min_singular_value": float(values.min()),
            "max_singular_value": float(values.max()),
            "lost_squared_overlap": float(max(len(gram_old), len(gram_new))-np.sum(values**2)),
            "status": "measured"}


def matched_energy_blocks(cross, old_energy, new_energy, gram_old, gram_new, tolerance):
    """Collect connected near-energy blocks after overlap-based assignment.

    Assignment is for identifying candidate correspondences only. Acceptance
    of degeneracy is never inferred from a single assigned MO overlap.
    Actual members, both energy spans and full subspace overlaps are retained.
    """
    rows, columns = linear_sum_assignment(-np.abs(cross)**2)
    parent = list(range(len(rows)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # Adjacent energy gaps form connected diagnostic windows; their complete
    # span is recorded, so chained near-degeneracy never claims exact equality.
    for indices, energies in ((rows, old_energy), (columns, new_energy)):
        ordered = np.argsort(energies[indices])
        for a, b in zip(ordered[:-1], ordered[1:]):
            if abs(energies[indices[a]]-energies[indices[b]]) <= tolerance:
                parent[root(int(b))] = root(int(a))
    groups = {}
    for i, (a, b) in enumerate(zip(rows, columns)):
        group = groups.setdefault(root(i), [[], []])
        group[0].append(int(a)); group[1].append(int(b))
    result = []
    for old_members, new_members in groups.values():
        a, b = np.asarray(old_members), np.asarray(new_members)
        measurement = principal_overlap(cross[np.ix_(a, b)], gram_old[np.ix_(a, a)], gram_new[np.ix_(b, b)])
        result.append({"old_members_zero_based": a.tolist(), "new_members_zero_based": b.tolist(),
                       "old_energy_range_hartree": [float(old_energy[a].min()), float(old_energy[a].max())],
                       "new_energy_range_hartree": [float(new_energy[b].min()), float(new_energy[b].max())],
                       "mean_energy_delta_hartree": float(new_energy[b].mean()-old_energy[a].mean()),
                       "principal_overlap": measurement})
    return result


def density_matrices(molecule, transform=None):
    alpha = molecule.mo.coeffsa
    beta = molecule.mo.coeffsb
    if transform is not None:
        alpha, beta = transform @ alpha, transform @ beta
    pa = (alpha*molecule.mo.occsa) @ alpha.T
    pb = (beta*molecule.mo.occsb) @ beta.T
    return {"mo_total": pa+pb, "mo_spin": pa-pb}


def density_self_norm(density, overlap):
    product = density @ overlap
    return float(np.einsum("ij,ji->", product, product))


def density_difference(old, new, s_old, s_new, s_cross):
    first, second = density_self_norm(old, s_old), density_self_norm(new, s_new)
    mixed = float(np.einsum("ij,ji->", old @ s_cross, new @ s_cross.T))
    squared = first+second-2*mixed
    # Cancellation affects a difference of nearly equal density operators.
    # Preserve the signed squared result and its arithmetic floor as evidence.
    floor = 64*np.finfo(float).eps*(abs(first)+abs(second)+2*abs(mixed))
    return {"hilbert_schmidt_squared_signed": squared,
            "hilbert_schmidt_norm": float(math.sqrt(max(squared, 0))),
            "roundoff_squared_scale": floor,
            "relative_norm": float(math.sqrt(max(squared, 0)/first)) if first > 1e-20 else None,
            "old_electron_trace": float(np.einsum("ij,ji->", old, s_old)),
            "new_electron_trace": float(np.einsum("ij,ji->", new, s_new))}


def compare_wavefunctions(old_path, new_path, evidence_npz, energy_window=1e-5):
    """Compare a same-case reproduction; all orbitals, independent AO integrals.

    This function does not classify the physical state as correct, nor decide
    whether a newly converged SCF solution is preferable to the old solution.
    """
    start = time.monotonic()
    timings = {}
    old, new = load_one(str(old_path)), load_one(str(new_path))
    timings["read_fchk_seconds"] = time.monotonic()-start
    if not np.array_equal(old.atnums, new.atnums):
        raise ValueError("Atom identity/order differs; an independently verified mapping is required")
    alignment = rigid_alignment(old.atcoords, new.atcoords)
    old_basis, new_basis = tuple(from_iodata(old)), tuple(from_iodata(new))
    aligned, transform = aligned_basis(new_basis, alignment)
    timings["prepare_alignment_seconds"] = time.monotonic()-start-sum(timings.values())
    # Integrate both self metrics and the cross metric. Residual geometry or
    # basis rounding is thereby included, rather than assuming identical AOs.
    combined = overlap_integral(old_basis+aligned)
    timings["ao_integrals_seconds"] = time.monotonic()-start-sum(timings.values())
    no, nn = old.mo.coeffs.shape[0], new.mo.coeffs.shape[0]
    if combined.shape != (no+nn, no+nn):
        raise ValueError("AO/segmented-contraction ordering does not match the MO coefficient dimensions")
    so, sn, sx = combined[:no, :no], combined[no:, no:], combined[:no, no:]
    arrays = {"overlap_old": so, "overlap_new_aligned": sn, "overlap_cross": sx,
              "ao_rotation_new_to_old_axes": transform, "coordinate_rotation": alignment["rotation"]}
    result = {"scope": "Independent same-case reproduction evidence; chemical adjudication required",
              "reference_libraries": {"iodata": version("qc-iodata"), "gbasis": version("qc-gbasis"),
                                      "numpy": np.__version__, "scipy": version("scipy")},
              "energy_window_hartree": energy_window,
              "alignment": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in alignment.items()},
              "energy_old_hartree": old.energy, "energy_new_hartree": new.energy,
              "energy_delta_hartree": None if old.energy is None or new.energy is None else new.energy-old.energy,
              "spins": {}, "densities": {}}
    for suffix, spin in (("a", "alpha"), ("b", "beta")):
        co = getattr(old.mo, "coeffs"+suffix)
        cn = transform @ getattr(new.mo, "coeffs"+suffix)
        oo, on = getattr(old.mo, "occs"+suffix), getattr(new.mo, "occs"+suffix)
        go, mo = metric_summary(co, so)
        gn, mn = metric_summary(cn, sn)
        cross = co.T @ sx @ cn
        io, inn = np.flatnonzero(oo > 1e-8), np.flatnonzero(on > 1e-8)
        arrays[spin+"_cross_mo_overlap"] = cross
        arrays[spin+"_gram_old"], arrays[spin+"_gram_new"] = go, gn
        blocks = matched_energy_blocks(cross, getattr(old.mo, "energies"+suffix),
                                       getattr(new.mo, "energies"+suffix), go, gn, energy_window)
        result["spins"][spin] = {"old_metric": mo, "new_metric": mn,
            "old_electrons": float(oo.sum()), "new_electrons": float(on.sum()),
            "represented_space": principal_overlap(cross, go, gn),
            "occupied_space": principal_overlap(cross[np.ix_(io, inn)], go[np.ix_(io, io)], gn[np.ix_(inn, inn)]),
            "energy_blocks": blocks,
            "unmatched_old_count": max(0, len(oo)-len(on)), "unmatched_new_count": max(0, len(on)-len(oo))}
    timings["all_spin_subspaces_seconds"] = time.monotonic()-start-sum(timings.values())
    po, pn = density_matrices(old), density_matrices(new, transform)
    for name in po:
        result["densities"][name] = density_difference(po[name], pn[name], so, sn, sx)
    for molecule, prefix, metric, matrices, rotation in ((old, "old", so, po, None), (new, "new", sn, pn, transform)):
        for key, name in (("scf", "mo_total"), ("scf_spin", "mo_spin")):
            if key not in molecule.one_rdms:
                result["densities"][prefix+"_producer_"+key] = {"status": "not_supplied_by_reader"}
                continue
            producer = molecule.one_rdms[key]
            if rotation is not None:
                producer = rotation @ producer @ rotation.T
            result["densities"][prefix+"_producer_"+key] = {
                "status": "measured", "difference_from_mo_density":
                density_difference(producer, matrices[name], metric, metric, metric)}
    eigen_old, eigen_new = np.linalg.eigvalsh(so), np.linalg.eigvalsh(sn)
    result["ao_metric"] = {}
    for name, values in (("old", eigen_old), ("new_aligned", eigen_new)):
        result["ao_metric"][name] = {"min_eigenvalue": float(values[0]), "max_eigenvalue": float(values[-1]),
            "condition_number": float(values[-1]/values[0]) if values[0] > 0 else None,
            "nonpositive_eigenvalues": int(np.count_nonzero(values <= 0))}
    Path(evidence_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(evidence_npz, **arrays)
    result["matrix_evidence"] = str(evidence_npz)
    timings["density_conditioning_and_evidence_seconds"] = time.monotonic()-start-sum(timings.values())
    result["wall_seconds"] = time.monotonic()-start
    result["timings"] = timings
    return result
