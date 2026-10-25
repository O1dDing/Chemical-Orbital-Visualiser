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

"""Negative controls for a reference reader: malformed files cannot pass."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from grid_reference_compare import compare_complete_spin_cube, grid_metrics, np, read_mo_cube
from independent_grid_reference import point_table

THRESHOLDS = {'nrms_max': 1e-4, 'abs_cosine_min': 1-1e-7, 'relative_peak_error_max': 1e-3}
HEADER = '''Test reference
Two orbital datasets at each grid point
-1 -1.0 -2.0 -3.0
2 0.5 0.0 0.0
2 0.0 0.5 0.0
2 0.0 0.0 0.5
1 1.0 0.0 0.0 0.0
'''


class GridReferenceTests(unittest.TestCase):
    def test_nonintegral_grid_counts_and_singular_coordinate_axes_are_rejected(self):
        for grid in ({'shape': [2, 2, 2.5], 'origin_bohr': [0, 0, 0], 'step_bohr': 1},
                     {'shape': [2, 2, 2], 'origin_bohr': [0, 0, 0],
                      'step_vectors_bohr': [[1, 0, 0], [0, 1, 0], [0, 0, 0]]}):
            with self.assertRaises(ValueError):
                point_table(grid)

    def read(self, body, header=HEADER):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'field.cube'
            path.write_text(header+body, encoding='ascii')
            return read_mo_cube(path)

    def test_multiple_dataset_point_interleaving_and_fortran_exponents(self):
        cube = self.read('2 7 12\n1D0 11 2 12 3 13 4 14 5 15 6 16 7 17 8 18\n')
        self.assertEqual(cube['orbital_ids'], [7, 12])
        np.testing.assert_array_equal(cube['values'], [np.arange(1, 9), np.arange(11, 19)])

    def test_gaussian_all_header_count_and_legacy_nval_are_consistent(self):
        body = '2 7 12\n' + '1 2 ' * 8
        for nval in (1, 2):
            cube = self.read(body, HEADER.replace('-1 -1.0 -2.0 -3.0', f'-1 -1.0 -2.0 -3.0 {nval}'))
            self.assertEqual(cube['values'].shape, (2, 8))
        for nval in (0, -1, 3):
            with self.assertRaises(ValueError):
                self.read(body, HEADER.replace('-1 -1.0 -2.0 -3.0', f'-1 -1.0 -2.0 -3.0 {nval}'))

    def test_complete_spin_uses_explicit_identities_and_all_columns(self):
        reference = np.array([np.arange(1., 9.), np.arange(2., 10.),
                              np.sin(np.arange(1., 9.)), np.cos(np.arange(1., 9.))])
        record = {'input_sha256': 'test', 'orbitals': [
            {'spin': spin, 'column_zero_based': offset + i, 'source_mo_one_based': i + 1}
            for spin, offset in (('alpha', 0), ('beta', 2)) for i in range(2)]}
        cube = {'orbital_ids': [1, 2], 'values': reference[:2], 'shape': [2, 2, 2], 'sha256': 'test'}
        with patch('grid_reference_compare.checked_reference_values', return_value=reference):
            alpha = compare_complete_spin_cube(cube, record, '.', 'alpha', THRESHOLDS)
            self.assertEqual(len(alpha), 2)
            self.assertTrue(all(r['pass'] for r in alpha))
            beta = compare_complete_spin_cube(cube, record, '.', 'beta', THRESHOLDS)
            self.assertTrue(all(not r['pass'] for r in beta))
            for spin, ids in (('unknown', [1, 2]), ('alpha', [1]), ('alpha', [1, 1]), ('alpha', [2, 1])):
                with self.assertRaises(ValueError):
                    compare_complete_spin_cube(dict(cube, orbital_ids=ids), record, '.', spin, THRESHOLDS)

    def test_truncation_extra_fields_invalid_tokens_and_duplicate_ids_are_rejected(self):
        for body in ('1 7\n1 2 3 4 5 6 7', '1 7\n1 2 3 4 5 6 7 8 9',
                     '1 7\n1 2 3 4 5 6 7 NaN', '1 7\n1 2 3 4 5 6 7 8 garbage',
                     '2 7 7\n'+'1 '*16):
            with self.subTest(body=body), self.assertRaises(ValueError):
                self.read(body)

    def test_density_and_unlabelled_units_cannot_be_mistaken_for_orbital_cube(self):
        for header in (HEADER.replace('-1 -1.0', '1 -1.0'), HEADER.replace('2 0.5', '-2 0.5')):
            with self.assertRaises(ValueError):
                self.read('1 7\n'+'1 '*8, header)

    def test_sign_is_allowed_but_scale_and_axis_permutation_fail(self):
        reference = np.arange(1., 9.).reshape(2, 2, 2)
        exact = grid_metrics(-reference, reference, THRESHOLDS)
        self.assertTrue(exact['pass'])
        self.assertEqual(exact['phase'], -1)
        self.assertFalse(grid_metrics(reference*1.01, reference, THRESHOLDS)['pass'])
        self.assertFalse(grid_metrics(reference.transpose(2, 1, 0), reference, THRESHOLDS)['pass'])

    def test_zero_and_nonfinite_fields_never_pass(self):
        self.assertEqual(grid_metrics(np.zeros(8), np.zeros(8), THRESHOLDS)['status'], 'insufficient')
        self.assertFalse(grid_metrics(np.ones(8), np.zeros(8), THRESHOLDS)['pass'])
        with self.assertRaises(ValueError):
            grid_metrics(np.array([float('inf')]), np.ones(1), THRESHOLDS)


if __name__ == '__main__':
    unittest.main()
