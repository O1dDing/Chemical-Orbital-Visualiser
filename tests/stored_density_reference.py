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

"""Audit raw Gaussian SCF density fields without IOData's ROHF suppression.

Density differences use the independently integrated AO metric. They are
one-electron operator norms, not pointwise-density or rendering errors.
No missing field, correlated density or diagnostic flag is a full COV verdict.
"""
from __future__ import annotations

from pathlib import Path
import re

from wavefunction_reference import np

SCALARS = {
    'Number of basis functions': 'I', 'Number of independent functions': 'I',
    'Number of electrons': 'I', 'Number of alpha electrons': 'I',
    'Number of beta electrons': 'I', 'Charge': 'I', 'Multiplicity': 'I',
}
ARRAYS = {'Alpha MO coefficients', 'Beta MO coefficients',
          'Total SCF Density', 'Spin SCF Density'}
HEADER = re.compile(r'^(.{40})\s+([IRCLH])\s+(.*?)\s*$')


def read_density_fields(path):
    """Read selected raw numeric fields, retaining other density declarations.

    This is a targeted reader, not a replacement for general FCHK validation.
    Numeric payloads of selected fields must have exactly their declared size;
    duplicate declarations and nonfinite/truncated payloads are rejected.
    """
    fields, declarations, other_densities = {}, {}, {}
    with Path(path).open(encoding='ascii', errors='strict') as stream:
        lines = iter(enumerate(stream, 1))
        try:
            _, title = next(lines)
            _, route = next(lines)
        except StopIteration as exc:
            raise ValueError('Missing two-line FCHK preamble') from exc
        for number, text in lines:
            match = HEADER.match(text.rstrip('\r\n'))
            if not match:
                continue
            name, kind, body = match.groups()
            name = name.strip()
            density = name.startswith(('Total ', 'Spin ')) and name.endswith(' Density')
            if name not in SCALARS and name not in ARRAYS and not density:
                continue
            if name in declarations:
                raise ValueError(f'Duplicate selected field: {name}')
            declaration = {'line': number, 'kind': kind}
            declarations[name] = declaration
            if name in SCALARS:
                if kind != SCALARS[name] or body.startswith('N='):
                    raise ValueError(f'Unexpected scalar field type: {name}')
                fields[name] = int(body)
                declaration['array_length'] = None
                continue
            if kind != 'R' or not body.startswith('N='):
                raise ValueError(f'Expected real array: {name}')
            count = int(body[2:].strip())
            if count < 0:
                raise ValueError(f'Negative array length: {name}')
            declaration['array_length'] = count
            if name not in ARRAYS:
                # Correlated density has a different meaning from an occupied
                # SCF MO sum. Only its declaration is inventoried here.
                other_densities[name] = dict(declaration)
                continue
            values = np.empty(count, dtype=float)
            cursor = 0
            while cursor < count:
                try:
                    _, row = next(lines)
                except StopIteration as exc:
                    raise ValueError(f'Truncated array: {name}') from exc
                tokens = row.split()
                if not tokens or len(tokens) > count-cursor:
                    raise ValueError(f'Invalid array payload length: {name}')
                try:
                    values[cursor:cursor+len(tokens)] = [
                        float(token.replace('D', 'E').replace('d', 'e')) for token in tokens]
                except ValueError as exc:
                    raise ValueError(f'Invalid numeric payload: {name}') from exc
                cursor += len(tokens)
            if not np.isfinite(values).all():
                raise ValueError(f'Nonfinite array: {name}')
            fields[name] = values
    return {'fields': fields, 'declarations': declarations,
            'other_density_declarations': other_densities,
            'title': title.strip(), 'method_header': route.rstrip()}


def unpack_lower_triangle(packed, dimension):
    """Gaussian stores rows of the lower triangle, with no off-diagonal factor."""
    packed = np.asarray(packed, dtype=float)
    if (not isinstance(dimension, (int, np.integer)) or dimension <= 0
            or packed.shape != (dimension*(dimension+1)//2,)
            or not np.isfinite(packed).all()):
        raise ValueError('A finite packed symmetric density of the exact AO dimension is required')
    rows, columns = np.tril_indices(dimension)
    result = np.empty((dimension, dimension), dtype=float)
    result[rows, columns] = packed
    result[columns, rows] = packed
    return result


def density_from_occupations(coefficients, occupations):
    coefficients, occupations = np.asarray(coefficients, float), np.asarray(occupations, float)
    if (coefficients.ndim != 2 or not coefficients.shape[0]
            or occupations.shape != (coefficients.shape[1],)
            or not np.isfinite(coefficients).all() or not np.isfinite(occupations).all()
            or np.any(occupations < 0)):
        raise ValueError('Finite compatible coefficients and nonnegative occupations required')
    selected = occupations != 0
    occupied = coefficients[:, selected]
    return (occupied*occupations[selected]) @ occupied.T


def density_record_metrics(stored, reconstructed, overlap, expected_trace):
    """Measure a direct matrix difference, without subtracting squared norms.

    If S=L L^T, ||L^T (Pstored-Pmo) L||F is the density-operator HS norm.
    This accepts rectangular MO-derived densities and signed spin matrices.
    AO entry/Frobenius differences are recorded separately and are not invariant
    to changes of AO coordinates. Positive definiteness comes from the actual S;
    no orbital-coefficient inversion, eigenvalue clipping or diagonal repair.
    """
    stored, reconstructed, overlap = (np.asarray(value, float) for value in (stored, reconstructed, overlap))
    if (stored.ndim != 2 or not len(stored) or stored.shape[0] != stored.shape[1]
            or reconstructed.shape != stored.shape or overlap.shape != stored.shape
            or not np.isfinite(expected_trace)
            or any(not np.isfinite(value).all() for value in (stored, reconstructed, overlap))
            or any(not np.allclose(value, value.T, rtol=0, atol=1e-10)
                   for value in (stored, reconstructed, overlap))):
        raise ValueError('Finite compatible symmetric densities and AO metric required')
    factor = np.linalg.cholesky(overlap)
    delta = stored-reconstructed
    physical_delta = factor.T @ delta @ factor
    physical_reference = factor.T @ reconstructed @ factor
    physical_stored = factor.T @ stored @ factor
    norm = float(np.linalg.norm(physical_reference))
    error = float(np.linalg.norm(physical_delta))
    ao_norm = float(np.linalg.norm(reconstructed))
    trace_stored = float(np.sum(stored*overlap.T))
    trace_reference = float(np.sum(reconstructed*overlap.T))
    index = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
    return {
        'expected_trace_electrons': float(expected_trace),
        'stored_trace_electrons': trace_stored,
        'mo_trace_electrons': trace_reference,
        'stored_trace_error_electrons': trace_stored-float(expected_trace),
        'mo_trace_error_electrons': trace_reference-float(expected_trace),
        'stored_minus_mo_trace_electrons': float(np.sum(delta*overlap.T)),
        'stored_hs_norm': float(np.linalg.norm(physical_stored)),
        'mo_hs_norm': norm,
        'hs_difference': error,
        'relative_hs_difference': error/norm if norm > 0 else None,
        'ao_frobenius_difference': float(np.linalg.norm(delta)),
        'relative_ao_frobenius_difference': float(np.linalg.norm(delta))/ao_norm if ao_norm > 0 else None,
        'ao_maximum_entry_difference': float(np.max(np.abs(delta))),
        'ao_maximum_difference_indices_zero_based': [int(x) for x in index],
        'stored_matrix_exactly_zero': bool(not np.any(stored)),
        'mo_matrix_exactly_zero': bool(not np.any(reconstructed)),
        'metric_scope': 'One-particle operator Hilbert-Schmidt norm in the independent AO metric',
        'relative_norm_reference': 'MO-reconstructed density operator',
        'formal_case_pass': False,
    }
