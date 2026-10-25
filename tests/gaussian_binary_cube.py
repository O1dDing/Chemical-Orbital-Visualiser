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

"""Strict reader for the observed Gaussian 16W A.03 unformatted MO cube ABI.

This profile has little-endian 4-byte Fortran record markers, 8-byte integers
and IEEE float64 values. Unsupported ABIs are rejected, never guessed. The
original Gaussian file remains unchanged and retains its source hash.
"""
from __future__ import annotations

import math
from pathlib import Path
import struct

from independent_grid_reference import np, sha256


def read_binary_mo_cube(path, *, values_path=None):
    path = Path(path)
    size = path.stat().st_size
    with path.open('rb') as stream:
        def record(expected_bytes):
            prefix = stream.read(4)
            if len(prefix) != 4 or struct.unpack('<I', prefix)[0] != expected_bytes:
                raise ValueError('Unexpected Gaussian binary record length or unsupported ABI')
            if expected_bytes > size - stream.tell() - 4:
                raise ValueError('Truncated Gaussian binary record payload')
            payload = stream.read(expected_bytes)
            if stream.read(4) != prefix:
                raise ValueError('Mismatched Gaussian binary record footer')
            return payload

        raw_comments = [record(80), record(80)]
        comments = [text.decode('latin1').rstrip(' \x00') for text in raw_comments]
        atom_count, x, y, z, nval = struct.unpack('<q3dq', record(40))
        if atom_count >= 0 or nval < 1:
            raise ValueError('Explicit orbital datasets are required in a binary MO cube')
        origin = np.array([x, y, z])
        shape, axes = [], []
        for _ in range(3):
            count, x, y, z = struct.unpack('<q3d', record(32))
            if count < 2:
                raise ValueError('Positive Bohr axes with at least two points are required')
            shape.append(count)
            axes.append([x, y, z])
        axes = np.array(axes)
        if not np.isfinite(origin).all() or not np.isfinite(axes).all() or abs(np.linalg.det(axes)) < 1e-20:
            raise ValueError('Invalid binary cube coordinate frame')
        atom_data = record(abs(atom_count) * 40)
        atoms = []
        for start in range(0, len(atom_data), 40):
            atomic_number, charge, x, y, z = struct.unpack_from('<q4d', atom_data, start)
            if atomic_number < 1 or not all(math.isfinite(v) for v in (charge, x, y, z)):
                raise ValueError('Invalid binary cube atom identity or coordinates')
            atoms.append({'atomic_number': atomic_number, 'nuclear_charge': charge,
                          'coordinates_bohr': [x, y, z]})
        orbital_data = record((nval + 1) * 8)
        identities = np.frombuffer(orbital_data, dtype='<i8')
        if int(identities[0]) != nval:
            raise ValueError('Binary cube NVAL and explicit orbital count disagree')
        orbital_ids = [int(value) for value in identities[1:]]
        if min(orbital_ids) < 1 or len(set(orbital_ids)) != nval:
            raise ValueError('Invalid or duplicate binary cube orbital identity')
        row_bytes = shape[2] * nval * 8
        rows = shape[0] * shape[1]
        if size - stream.tell() != rows * (row_bytes + 8):
            raise ValueError('Truncated or extra binary cube rows')
        value_shape = (nval, math.prod(shape))
        if values_path is None:
            if math.prod(value_shape) * 8 > 1024**3:
                raise ValueError('Large binary grid requires an explicit memory-mapped values_path')
            values = np.empty(value_shape, dtype=np.float64)
        else:
            values_path = Path(values_path)
            if values_path.exists():
                raise FileExistsError('Never overwrite a prior grid artifact')
            values = np.lib.format.open_memmap(values_path, mode='w+', dtype=np.float64, shape=value_shape)
        for row in range(rows):
            numbers = np.frombuffer(record(row_bytes), dtype='<f8')
            if not np.isfinite(numbers).all():
                raise ValueError('Non-finite binary orbital grid value')
            begin = row * shape[2]
            values[:, begin:begin + shape[2]] = numbers.reshape(shape[2], nval).T
        if stream.read(1):
            raise ValueError('Trailing binary cube bytes')
    if isinstance(values, np.memmap):
        values.flush()
    values.flags.writeable = False
    return {'path': str(path.resolve()), 'sha256': sha256(path), 'comments': comments,
            'raw_comment_bytes_hex': [text.hex() for text in raw_comments],
            'origin_bohr': origin, 'shape': shape, 'step_vectors_bohr': axes,
            'atoms': atoms, 'orbital_ids': orbital_ids, 'values': values,
            'binary_profile': 'Gaussian16W-A03-LE-marker32-int64-real64'}


def write_single_orbital_cube(cube, dataset_index, destination):
    """Package one original dataset as text with round-trip float64 precision.

    This changes the storage format only. Callers retain the original binary
    file, dataset identity and hash, and independently verify numeric equality.
    """
    if not isinstance(dataset_index, int) or not 0 <= dataset_index < len(cube['orbital_ids']):
        raise ValueError('Dataset index outside the original Gaussian cube')
    destination = Path(destination)
    values = cube['values'][dataset_index]
    if values.size != math.prod(cube['shape']) or not np.isfinite(values).all():
        raise ValueError('Invalid dataset values')
    with destination.open('x', encoding='ascii', newline='\n') as stream:
        stream.write('Lossless packaging of original Gaussian binary MO data\n')
        stream.write('Parent SHA256 ' + cube['sha256'] + '\n')
        stream.write(str(-len(cube['atoms'])) + ' ' + ' '.join(f'{float(v):.17e}' for v in cube['origin_bohr']) + ' 1\n')
        for count, axis in zip(cube['shape'], cube['step_vectors_bohr']):
            stream.write(str(count) + ' ' + ' '.join(f'{float(v):.17e}' for v in axis) + '\n')
        for atom in cube['atoms']:
            stream.write(str(atom['atomic_number']) + ' ' + ' '.join(f'{float(v):.17e}' for v in
                         [atom['nuclear_charge'], *atom['coordinates_bohr']]) + '\n')
        stream.write('1 ' + str(cube['orbital_ids'][dataset_index]) + '\n')
        full = values.size // 6 * 6
        np.savetxt(stream, values[:full].reshape(-1, 6), fmt='%.17e')
        if full < values.size:
            stream.write(' '.join(f'{float(v):.17e}' for v in values[full:]) + '\n')
    return destination
