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

"""Strict Gaussian MO cube reader and independent full-grid comparisons.

Only same-grid, same-input evidence is comparable. Global orbital sign is
reported and may be reversed; no amplitude fitting, rotation, or interpolation
is performed. A passing grid is not a formal molecule acceptance verdict.
"""
from __future__ import annotations

import math
from pathlib import Path

from independent_grid_reference import np, sha256


def number(token):
    value = float(token.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise ValueError('Non-finite cube field')
    return value


def read_mo_cube(path):
    with Path(path).open(encoding='ascii') as stream:
        comments = [stream.readline().rstrip('\n'), stream.readline().rstrip('\n')]
        row = stream.readline().split()
        if len(row) not in (4, 5):
            raise ValueError('Malformed cube origin header')
        natom = int(row[0])
        if natom >= 0:
            raise ValueError('An orbital cube with explicit dataset identities is required')
        header_nval = int(row[4]) if len(row) == 5 else None
        if header_nval is not None and header_nval < 1:
            raise ValueError('Cube header NVAL must be positive')
        origin = np.array([number(value) for value in row[1:4]])
        shape, axes = [], []
        for _ in range(3):
            row = stream.readline().split()
            if len(row) != 4 or int(row[0]) < 2:
                raise ValueError('Positive Bohr axes with at least two points are required')
            shape.append(int(row[0]))
            axes.append([number(value) for value in row[1:]])
        axes = np.asarray(axes)
        if abs(np.linalg.det(axes)) < 1e-20:
            raise ValueError('Cube axes do not span three dimensions')
        atoms = []
        for _ in range(abs(natom)):
            row = stream.readline().split()
            if len(row) != 5:
                raise ValueError('Malformed cube atom record')
            atoms.append({'atomic_number': int(row[0]), 'nuclear_charge': number(row[1]),
                          'coordinates_bohr': [number(value) for value in row[2:]]})
        tokens = stream.read().split()
    if not tokens:
        raise ValueError('Missing orbital cube dataset header')
    datasets = int(tokens[0])
    if datasets < 1:
        raise ValueError('At least one orbital dataset is required')
    # Gaussian 16W's AMO/BMO=All writes the dataset count as NVAL. Other
    # orbital-cube writers leave NVAL absent or one and specify it only in
    # the following MO identity record. Both conventions are unambiguous.
    if header_nval not in (None, 1, datasets):
        raise ValueError('Cube header NVAL disagrees with its orbital dataset count')
    ids = [int(token) for token in tokens[1:1+datasets]]
    if len(ids) != datasets or len(set(ids)) != datasets:
        raise ValueError('Missing or duplicate cube orbital identity')
    count = math.prod(shape)
    body = tokens[1+datasets:]
    if len(body) != count*datasets:
        raise ValueError('Truncated or extra cube grid values')
    # Gaussian interleaves datasets at each point. Preserve C-order (z fastest).
    values = np.fromiter((number(token) for token in body), dtype=np.float64).reshape(count, datasets).T
    return {'path': str(Path(path).resolve()), 'sha256': sha256(path), 'comments': comments,
            'origin_bohr': origin, 'shape': shape, 'step_vectors_bohr': axes,
            'atoms': atoms, 'orbital_ids': ids, 'values': values}


def grid_metrics(actual, reference, thresholds):
    actual, reference = np.asarray(actual, float), np.asarray(reference, float)
    if actual.shape != reference.shape or not actual.size:
        raise ValueError('Nonempty, equal-sized grid arrays are required')
    if not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError('Non-finite orbital grid values')
    actual, reference = actual.ravel(), reference.ravel()
    scale = max(float(np.max(np.abs(actual))), float(np.max(np.abs(reference))))
    if scale == 0:
        return {'status': 'insufficient', 'pass': False, 'reason': 'Both sampled fields are zero',
                'points': actual.size}
    a, r = actual/scale, reference/scale
    aa, rr, ar = float(a@a), float(r@r), float(a@r)
    if aa == 0 or rr == 0:
        return {'status': 'insufficient', 'pass': False, 'reason': 'One sampled field has zero norm',
                'points': actual.size}
    phase = 1 if ar >= 0 else -1
    error = phase*a-r
    nrms = math.sqrt(float(error@error)/rr)
    cosine = min(1., abs(ar)/math.sqrt(aa*rr))
    relative_peak = float(np.max(np.abs(error)))/float(np.max(np.abs(r)))
    passed = (nrms <= thresholds['nrms_max'] and cosine >= thresholds['abs_cosine_min']
              and relative_peak <= thresholds['relative_peak_error_max'])
    return {'status': 'pass' if passed else 'fail', 'pass': passed, 'points': actual.size,
            'nrms': nrms, 'abs_cosine': cosine, 'relative_peak_error': relative_peak, 'phase': phase,
            'actual_peak': float(np.max(np.abs(actual))), 'reference_peak': float(np.max(np.abs(reference))),
            'max_absolute_error': float(np.max(np.abs(error)))*scale,
            'thresholds': dict(thresholds), 'amplitude_fitting': False}


def checked_reference_values(cube, record, directory):
    """Verify the common input, coordinates and reference artifact once per cube."""
    from iodata import load_one
    directory = Path(directory)
    input_path = Path(record['source_input'])
    if sha256(input_path) != record['input_sha256']:
        raise ValueError('Independent reference input identity changed')
    grid = record['grid']
    if cube['shape'] != grid['shape'] or not np.allclose(cube['origin_bohr'], grid['origin_bohr'], atol=1e-12, rtol=0):
        raise ValueError('Different cube/reference point grids')
    if not np.allclose(cube['step_vectors_bohr'], grid['step_vectors_bohr'], atol=1e-12, rtol=0):
        raise ValueError('Different cube/reference step vectors')
    mol = load_one(str(input_path))
    if len(cube['atoms']) != len(mol.atnums) or any(row['atomic_number'] != int(z) for row, z in zip(cube['atoms'], mol.atnums)):
        raise ValueError('Cube/reference atomic identities differ')
    coordinates = np.array([row['coordinates_bohr'] for row in cube['atoms']])
    if not np.allclose(coordinates, mol.atcoords, atol=6e-7, rtol=0):
        raise ValueError('Cube/reference atom coordinate frames differ beyond cube rounding precision')
    mo_path = directory / 'mo_values.npy'
    if sha256(mo_path) != record['artifacts'][mo_path.name]['sha256']:
        raise ValueError('Independent grid artifact identity changed')
    values = np.load(mo_path, mmap_mode='r', allow_pickle=False)
    if values.shape != (record['nmo'], math.prod(grid['shape'])):
        raise ValueError('Independent grid array dimensions disagree with its metadata')
    return values


def compare_cube_to_reference(cube, record, directory, column, expected_orbital, thresholds):
    """Compare one explicitly identified spin/MO column, after its batch barrier."""
    if cube['orbital_ids'] != [expected_orbital]:
        raise ValueError('Cube dataset does not contain exactly the requested orbital')
    if column < 0 or column >= record['nmo']:
        raise ValueError('MO column outside independent reference')
    values = checked_reference_values(cube, record, directory)
    result = grid_metrics(cube['values'][0], values[column], thresholds)
    result.update(cube_sha256=cube['sha256'], reference_input_sha256=record['input_sha256'],
                  orbital_identity=record['orbitals'][column], scientific_case_pass=False)
    return result


def compare_complete_spin_cube(cube, record, directory, spin, thresholds):
    """Compare every dataset in one explicit Gaussian AMO/BMO=All output.

    The caller must verify the frozen producer's input hash, command and cube
    hash before selecting spin. Neither dataset order nor comments alone are
    treated as sufficient provenance.
    """
    if spin not in ('alpha', 'beta'):
        raise ValueError('An explicit alpha or beta source is required')
    orbitals = [o for o in record['orbitals'] if o['spin'] == spin]
    if not orbitals or cube['orbital_ids'] != [o['source_mo_one_based'] for o in orbitals]:
        raise ValueError('Cube must contain every requested spin orbital exactly once, in source order')
    if cube['values'].shape != (len(orbitals), math.prod(cube['shape'])):
        raise ValueError('Cube dataset dimensions disagree with orbital identities')
    values = checked_reference_values(cube, record, directory)
    results = []
    for index, identity in enumerate(orbitals):
        column = identity['column_zero_based']
        result = grid_metrics(cube['values'][index], values[column], thresholds)
        result.update(cube_sha256=cube['sha256'], reference_input_sha256=record['input_sha256'],
                      orbital_identity=identity, scientific_case_pass=False)
        results.append(result)
    return results
