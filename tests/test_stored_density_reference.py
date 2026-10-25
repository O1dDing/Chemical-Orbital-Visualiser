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

import tempfile
from pathlib import Path
import unittest

from stored_density_reference import (density_from_occupations, density_record_metrics,
                                      np, read_density_fields, unpack_lower_triangle)


class StoredDensityReferenceTests(unittest.TestCase):
    def read_fixture(self, payload):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'raw.fchk'
            path.write_text('Synthetic density reference\nSP UHF STO-3G\n'+payload, encoding='ascii')
            return read_density_fields(path)

    @staticmethod
    def field(name, kind, body):
        return f'{name:40s}   {kind}   {body}\n'

    def test_raw_open_shell_fields_are_retained_and_correlated_density_is_separate(self):
        raw = self.read_fixture(
            self.field('Number of alpha electrons', 'I', '2')+
            self.field('Number of beta electrons', 'I', '1')+
            self.field('Total SCF Density', 'R', 'N= 3')+'2.0D+00 0.0D+00 1.0D+00\n'+
            self.field('Spin SCF Density', 'R', 'N= 3')+'0.0 0.0 1.0\n'+
            self.field('Total MP2 Density', 'R', 'N= 3')+'1.9 0.0 1.1\n')
        np.testing.assert_array_equal(raw['fields']['Total SCF Density'], [2, 0, 1])
        np.testing.assert_array_equal(raw['fields']['Spin SCF Density'], [0, 0, 1])
        self.assertNotIn('Total MP2 Density', raw['fields'])
        self.assertEqual(raw['other_density_declarations']['Total MP2 Density']['array_length'], 3)
        missing = self.read_fixture(self.field('Number of basis functions', 'I', '2'))
        self.assertNotIn('Total SCF Density', missing['fields'])

    def test_malformed_selected_fields_never_become_a_missing_or_zero_density(self):
        header = self.field('Total SCF Density', 'R', 'N= 3')
        for payload in (header+'1 2\n', header+'1 2 3 4\n', header+'1 NaN 3\n',
                        header+'1 2 3\n'+header+'1 2 3\n',
                        header+self.field('Number of basis functions', 'I', '2'),
                        self.field('Total SCF Density', 'I', 'N= 3')+'1 2 3\n'):
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                self.read_fixture(payload)

    def test_packed_order_off_diagonal_weight_and_rectangular_occupied_space(self):
        coefficients = np.array([[1., 0.], [0., .6], [0., .8]])
        alpha = density_from_occupations(coefficients, [1, 1])
        beta = density_from_occupations(coefficients, [1, 0])
        stored = unpack_lower_triangle([2, 0, .36, 0, .48, .64], 3)
        np.testing.assert_allclose(stored, alpha+beta, atol=1e-15)
        correct = density_record_metrics(stored, alpha+beta, np.eye(3), 3)
        self.assertLess(correct['hs_difference'], 1e-15)
        self.assertAlmostEqual(correct['stored_trace_electrons'], 3)
        damaged = stored.copy()
        damaged[1, 2] *= 2
        damaged[2, 1] *= 2
        self.assertGreater(density_record_metrics(damaged, alpha+beta, np.eye(3), 3)['relative_hs_difference'], .1)
        with self.assertRaises(ValueError):
            unpack_lower_triangle([1, 2, 3, 4], 3)

    def test_metric_result_is_invariant_to_consistent_nonorthogonal_ao_change(self):
        overlap = np.array([[1., .3], [.3, 1.2]])
        coefficients = np.linalg.inv(np.linalg.cholesky(overlap).T)
        reference = density_from_occupations(coefficients, [1, .5])
        stored = reference+np.array([[.002, -.001], [-.001, .003]])
        original = density_record_metrics(stored, reference, overlap, 1.5)
        transform = np.array([[1.4, .2], [-.1, .7]])
        inverse = np.linalg.inv(transform)
        changed = density_record_metrics(inverse@stored@inverse.T, inverse@reference@inverse.T,
                                         transform.T@overlap@transform, 1.5)
        for key in ('stored_trace_electrons', 'mo_trace_electrons', 'hs_difference', 'relative_hs_difference'):
            self.assertAlmostEqual(original[key], changed[key], places=13)
        wrong = density_record_metrics(stored[::-1, ::-1], reference, overlap, 1.5)
        self.assertGreater(wrong['hs_difference'], 10*original['hs_difference'])

    def test_signed_spin_zero_channel_and_small_direct_difference(self):
        reference = np.diag([1., -1.])
        stored = reference+np.diag([1e-9, -1e-9])
        result = density_record_metrics(stored, reference, np.eye(2), 0)
        self.assertAlmostEqual(result['relative_hs_difference'], 1e-9, delta=1e-15)
        self.assertAlmostEqual(result['stored_trace_electrons'], 0)
        zero = np.zeros((2, 2))
        empty = density_record_metrics(zero, zero, np.eye(2), 0)
        self.assertIsNone(empty['relative_hs_difference'])
        self.assertTrue(empty['mo_matrix_exactly_zero'])
        nonzero = density_record_metrics(np.diag([1e-7, 0]), zero, np.eye(2), 0)
        self.assertGreater(nonzero['hs_difference'], 0)
        self.assertIsNone(nonzero['relative_hs_difference'])

    def test_invalid_metric_and_occupation_inputs_are_rejected(self):
        for metric in (np.diag([1., 0]), np.diag([1., -1]), np.array([[1., .2], [.3, 1.]])):
            with self.subTest(metric=metric), self.assertRaises((ValueError, np.linalg.LinAlgError)):
                density_record_metrics(np.eye(2), np.eye(2), metric, 2)
        with self.assertRaises(ValueError):
            density_from_occupations(np.eye(2), [1, -1])
        with self.assertRaises(ValueError):
            density_from_occupations(np.eye(2), [1])


if __name__ == '__main__':
    unittest.main()
