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

"""Collect full AO/MO grids with IOData/GBasis, independently of COV.

Outputs are reference evidence, never an acceptance verdict. The FCHK is read
without edits. A manifest records the complete point order and alpha/beta
column identity, so native texture layouts can be compared explicitly later.
Run this worker under the campaign's process-tree CPU and memory limiter.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import sys
import time

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
DEPS = Path(r'F:\Dev\cov-validation-20260905\reference-deps')
sys.path.insert(0, str(DEPS))
import numpy as np
from gbasis.evals.eval import evaluate_basis
from gbasis.wrappers import from_iodata
from iodata import load_one

from validation_process import atomic_json


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            result.update(block)
    return result.hexdigest()


def point_table(grid):
    original_shape = np.asarray(grid['shape'])
    shape = np.asarray(grid['shape'], dtype=np.int64)
    origin = np.asarray(grid['origin_bohr'], dtype=np.float64)
    if shape.shape != (3,) or np.any(shape < 2) or origin.shape != (3,) or not np.array_equal(shape, original_shape):
        raise ValueError('A three-dimensional grid with at least two points per axis is required')
    if 'step_vectors_bohr' in grid:
        axes = np.asarray(grid['step_vectors_bohr'], dtype=np.float64)
    else:
        axes = np.eye(3) * float(grid['step_bohr'])
    if axes.shape != (3, 3) or not np.isfinite(axes).all() or not np.isfinite(origin).all():
        raise ValueError('Non-finite or malformed grid coordinates')
    if abs(np.linalg.det(axes)) < 1e-20:
        raise ValueError('Grid step vectors must span three dimensions')
    indices = np.indices(shape).reshape(3, -1).T
    return origin + indices @ axes, shape.tolist(), axes.tolist()


def collect(input_path, specification, output, chunk_points=2048):
    started = time.monotonic()
    if output.exists():
        raise FileExistsError('Existing reference evidence is never overwritten: ' + str(output))
    if chunk_points < 1:
        raise ValueError('Positive chunk size required')
    specification_document = json.loads(specification.read_text(encoding='utf-8'))
    grid = specification_document.get('grid', specification_document)
    points, shape, axes = point_table(grid)
    input_digest = sha256(input_path)
    mol = load_one(str(input_path))
    if mol.mo is None or mol.mo.kind not in ('restricted', 'unrestricted'):
        raise ValueError('A restricted or unrestricted molecular orbital block is required')
    basis = tuple(from_iodata(mol))
    coeff = np.asarray(mol.mo.coeffs, dtype=np.float64)
    if coeff.ndim != 2 or not np.isfinite(coeff).all():
        raise ValueError('Missing or non-finite orbital coefficients')
    columns = [{'column_zero_based': i, 'spin': 'alpha' if mol.mo.kind == 'unrestricted' else 'shared_spatial',
                'source_mo_one_based': i+1} for i in range(mol.mo.norba)]
    if mol.mo.kind == 'unrestricted':
        columns.extend({'column_zero_based': mol.mo.norba+i, 'spin': 'beta',
                        'source_mo_one_based': i+1} for i in range(mol.mo.norbb))
    if len(columns) != coeff.shape[1]:
        raise ValueError('Orbital spin/column counts do not match the coefficient matrix')
    output.mkdir(parents=True)
    shutil.copy2(Path(__file__), output / 'collector-source.py')
    np.save(output / 'points_bohr.npy', points, allow_pickle=False)
    np.save(output / 'coefficients.npy', coeff, allow_pickle=False)
    np.save(output / 'energies_hartree.npy', mol.mo.energies, allow_pickle=False)
    np.save(output / 'occupations.npy', mol.mo.occs, allow_pickle=False)
    ao_path, mo_path = output / 'ao_values.partial.npy', output / 'mo_values.partial.npy'
    ao = np.lib.format.open_memmap(ao_path, mode='w+', dtype='<f8', shape=(coeff.shape[0], len(points)))
    mo = np.lib.format.open_memmap(mo_path, mode='w+', dtype='<f8', shape=(coeff.shape[1], len(points)))
    for start in range(0, len(points), chunk_points):
        stop = min(start+chunk_points, len(points))
        values = evaluate_basis(basis, points[start:stop])
        if values.shape != (coeff.shape[0], stop-start) or not np.isfinite(values).all():
            raise ValueError('Reference AO evaluation did not produce the expected finite grid')
        orbital_values = coeff.T @ values
        if not np.isfinite(orbital_values).all():
            raise ValueError('Non-finite MO grid values')
        ao[:, start:stop] = values
        mo[:, start:stop] = orbital_values
    ao.flush()
    mo.flush()
    del ao, mo
    ao_path.replace(output / 'ao_values.npy')
    mo_path.replace(output / 'mo_values.npy')
    # Re-read completed files; publication of collection.json is the barrier.
    for filename, expected in (('ao_values.npy', (coeff.shape[0], len(points))),
                               ('mo_values.npy', (coeff.shape[1], len(points)))):
        data = np.load(output / filename, mmap_mode='r', allow_pickle=False)
        if data.shape != expected or data.dtype != np.dtype('<f8'):
            raise ValueError('Saved grid shape or precision does not match its source')
        del data
    if sha256(input_path) != input_digest:
        raise ValueError('Input changed during reference collection')
    implementations = {str(path.relative_to(DEPS)): sha256(path)
                       for package in ('iodata', 'gbasis')
                       for path in sorted((DEPS / package).rglob('*.py'))}
    record = {'schema_version': 1, 'status': 'collected', 'scientific_pass': False,
              'source_input': str(input_path.resolve()), 'input_sha256': input_digest,
              'specification': str(specification.resolve()), 'specification_sha256': sha256(specification),
              'grid': {'shape': shape, 'origin_bohr': grid['origin_bohr'], 'step_vectors_bohr': axes,
                       'point_order': 'C order for (ix, iy, iz), iz fastest; explicit points_bohr.npy is authoritative',
                       'point_count': len(points), 'unit': 'bohr'},
              'orbitals': columns, 'mo_kind': mol.mo.kind, 'nao': coeff.shape[0], 'nmo': coeff.shape[1],
              'data_type': 'little-endian IEEE-754 float64',
              'grid_array_axes': ['basis_or_orbital_index', 'point_index'],
              'ao_convention_source': 'IOData FCHK conventions, passed unmodified to GBasis from_iodata',
              'units': {'ao_values': 'bohr^(-3/2)', 'mo_values': 'bohr^(-3/2)',
                        'coefficients': 'dimensionless', 'energies_hartree': 'hartree',
                        'occupations': 'electrons'},
              'libraries': {key: version(package) for key, package in
                            (('iodata', 'qc-iodata'), ('gbasis', 'qc-gbasis'), ('numpy', 'numpy'))},
              'implementation_sha256': sha256(Path(__file__)), 'library_implementations': implementations,
              'artifacts': {path.name: {'bytes': path.stat().st_size, 'sha256': sha256(path)}
                            for path in sorted(output.glob('*.npy'))},
              'wall_seconds': time.monotonic()-started,
              'scope': 'Independent full-grid reference collection; no COV scientific, UI, or full-space acceptance claim'}
    atomic_json(output / 'collection.json', record)
    print(json.dumps({'input': str(input_path), 'nao': record['nao'], 'nmo': record['nmo'],
                      'points': len(points), 'wall_seconds': record['wall_seconds'], 'status': 'collected'}), flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('input', 'specification', 'output'):
        parser.add_argument('--'+option, type=Path, required=True)
    parser.add_argument('--chunk-points', type=int, default=2048)
    args = parser.parse_args()
    collect(args.input, args.specification, args.output, args.chunk_points)


if __name__ == '__main__':
    main()
