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

"""Collect raw all-MO operator matrices without deciding orbital labels."""
from __future__ import annotations

from pathlib import Path
import time

from nuclear_symmetry_reference import (ao_operation_integrals, checked_operation,
                                        match_operation, np)
from wavefunction_reference import from_iodata, load_one, overlap_integral
from scipy.spatial.transform import Rotation


def operations_for_nuclear_evidence(nuclear):
    """Complete finite candidate list, or explicitly incomplete continuous probes."""
    if nuclear['kind'] == 'finite_candidates':
        return [{"id": f"operation-{index:03d}", "matrix": item['matrix'],
                 "source": "frozen_nuclear_candidate", "nuclear_displacement_bohr": item['max_displacement_bohr']}
                for index, item in enumerate(nuclear['operations'])]
    identity = {'id': 'identity', 'matrix': np.eye(3).tolist(), 'source': 'continuous_group_probe'}
    if nuclear['kind'] == 'continuous_linear':
        axis = np.asarray(nuclear['axis'])
        probe = np.eye(3)[int(np.argmin(np.abs(axis)))]
        perpendicular = np.cross(axis, probe)
        perpendicular /= np.linalg.norm(perpendicular)
        operators = [identity]
        for index, angle in enumerate((np.pi/4, np.pi/3, 2*np.pi/7)):
            operators.append({'id': f'axial-rotation-{index}',
                'matrix': Rotation.from_rotvec(axis*angle).as_matrix().tolist(),
                'source': 'continuous_group_probe', 'angle_radians': angle})
        operators.append({'id': 'vertical-reflection',
            'matrix': (np.eye(3)-2*np.outer(perpendicular, perpendicular)).tolist(),
            'source': 'continuous_group_probe'})
        if nuclear['group_under_tolerance'] == 'Dinfh':
            operators.append({'id': 'inversion', 'matrix': (-np.eye(3)).tolist(), 'source': 'continuous_group_probe'})
        return operators
    if nuclear['kind'] == 'continuous_spherical':
        operators = [identity, {'id': 'inversion', 'matrix': (-np.eye(3)).tolist(), 'source': 'continuous_group_probe'}]
        for index, axis in enumerate(np.eye(3)):
            operators.append({'id': f'rotation-axis-{index}',
                'matrix': Rotation.from_rotvec(axis*np.sqrt(2)).as_matrix().tolist(),
                'source': 'continuous_group_probe', 'angle_radians': float(np.sqrt(2))})
        return operators
    raise ValueError('Nuclear geometry does not provide a resolvable operation scope')


def collect_mo_operators(fchk, operators, origin, output, *, cached_overlap=None):
    """Store C^T S C and C^T K(R) C for every supplied MO and operation.

    No coefficient normalization, phase fitting, energy grouping or irrep
    selection happens during collection. The full raw matrices permit later
    independent block, whole-space and alpha/beta density analysis.
    """
    started = time.time()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    molecule = load_one(str(fchk))
    basis = from_iodata(molecule)
    coefficients = {'alpha': np.asarray(molecule.mo.coeffsa)}
    if molecule.mo.kind == 'unrestricted':
        coefficients['beta'] = np.asarray(molecule.mo.coeffsb)
    elif molecule.mo.kind != 'restricted':
        raise ValueError(f'Unsupported MO storage kind: {molecule.mo.kind}')
    overlap = overlap_integral(basis) if cached_overlap is None else np.asarray(cached_overlap)
    n = molecule.obasis.nbasis
    if overlap.shape != (n, n) or not np.isfinite(overlap).all():
        raise ValueError('Invalid independent AO overlap')
    if any(c.ndim != 2 or c.shape[0] != n or not np.isfinite(c).all() for c in coefficients.values()):
        raise ValueError('Invalid original coefficient matrix')
    labels = list(zip(molecule.atnums.tolist(), molecule.atcorenums.tolist()))
    identity_data = {'overlap': overlap, 'atomic_numbers': molecule.atnums, 'coordinates_bohr': molecule.atcoords,
                     'origin_bohr': np.asarray(origin), 'alpha_energies': molecule.mo.energiesa,
                     'beta_energies': molecule.mo.energiesb, 'alpha_occupations': molecule.mo.occsa,
                     'beta_occupations': molecule.mo.occsb}
    spins = {}
    for spin, c in coefficients.items():
        gram = c.T @ overlap @ c
        identity_data[f'{spin}_gram'] = gram
        spins[spin] = {'ao_count': c.shape[0], 'mo_count': c.shape[1],
                       'source_mos_one_based': list(range(1, c.shape[1]+1)),
                       'max_input_gram_error': float(np.max(np.abs(gram-np.eye(c.shape[1]))))}
    np.savez_compressed(output / 'identity.npz', **identity_data)
    records = []
    for index, operator in enumerate(operators):
        op_start = time.time()
        matrix = checked_operation(operator['matrix'])
        # The independent cached S is used for E only, with an arithmetic-level
        # identity test. Every nonidentity transformed basis is integrated.
        if np.linalg.norm(matrix-np.eye(3)) <= 1e-13:
            observed_overlap, pullback = overlap, overlap
        else:
            observed_overlap, pullback = ao_operation_integrals(basis, matrix, origin)
        delta = float(np.max(np.abs(observed_overlap-overlap)))
        matrices = {f'{spin}_raw_mo_operator': c.T @ pullback @ c for spin, c in coefficients.items()}
        if not all(np.isfinite(value).all() for value in matrices.values()):
            raise ValueError('Nonfinite MO operation matrix')
        filename = f'operation-{index:03d}.npz'
        np.savez_compressed(output / filename, **matrices)
        records.append({'id': operator['id'], 'matrix': matrix.tolist(), 'artifact': filename,
                        'source': operator['source'], 'overlap_recalculation_max_difference': delta,
                        'nuclear_matching': match_operation(molecule.atcoords, labels, matrix, origin),
                        'wall_seconds': time.time()-op_start})
    return {'status': 'all_mo_operators_collected', 'mo_kind': molecule.mo.kind,
            'spins': spins, 'restricted_beta_uses_alpha_coefficient_space': molecule.mo.kind == 'restricted',
            'operation_count': len(records), 'operations': records, 'identity_artifact': 'identity.npz',
            'started_epoch': started, 'finished_epoch': time.time(),
            'scope': 'raw all-MO integral measurements; no orbital labels or scientific pass assigned',
            'formal_case_pass': False}
