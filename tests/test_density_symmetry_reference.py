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

"""Exact finite-Hilbert-space counterexamples for density symmetry analysis."""
import unittest

from density_symmetry_reference import (density_operation_difference, energy_connected_blocks,
                                       np, prepare_spin_densities, subspace_measurement)


class DensitySymmetryReferenceTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(186435)

    def test_operator_difference_matches_full_space_with_outside_subspace_leakage(self):
        rotation, _ = np.linalg.qr(self.rng.normal(size=(7, 7)))
        q, _ = np.linalg.qr(self.rng.normal(size=(7, 3)))
        change = np.array([[1.1, .3, -.2], [0., .8, .1], [.1, -.2, 1.2]])
        alpha = q @ change
        beta = alpha @ np.array([[.9, .2, 0.], [.1, 1.1, -.2], [.15, 0., .8]])
        fa, fb = np.array([1., .6, .1]), np.array([.7, .4, .2])
        original_alpha, original_beta = alpha.copy(), beta.copy()
        states = prepare_spin_densities(alpha, beta, np.eye(7), fa, fb, shared_coefficients=False)
        inverse = states['alpha_inverse_sqrt_gram']
        projected_operator = inverse @ (alpha.T @ rotation @ alpha) @ inverse
        pa, pb = (alpha*fa) @ alpha.T, (beta*fb) @ beta.T
        for name, density in (('alpha', pa), ('beta', pb), ('total', pa+pb), ('spin', pa-pb)):
            measured = density_operation_difference(states['densities'][name], states['densities'][name], projected_operator)
            direct = np.linalg.norm(rotation @ density @ rotation.T-density)
            self.assertAlmostEqual(measured['full_space_hs_difference_squared_signed'], direct**2, places=11)
            self.assertLessEqual(measured['projected_direct_hs_difference'], direct+1e-12)
            self.assertFalse(measured['strict_equality_claimed'])
        self.assertLess(states['beta_density_projection_error_bound_hs'], 1e-12)
        np.testing.assert_array_equal(alpha, original_alpha)
        np.testing.assert_array_equal(beta, original_beta)

    def test_total_density_can_be_symmetric_while_each_spin_is_exchanged(self):
        alpha = np.eye(2)
        beta = np.array([[0., 1.], [1., 0.]])
        states = prepare_spin_densities(alpha, beta, np.eye(2), [1., 0.], [1., 0.], shared_coefficients=False)
        exchange = beta
        densities = states['densities']
        total = density_operation_difference(densities['total'], densities['total'], exchange)
        spin = density_operation_difference(densities['spin'], densities['spin'], exchange)
        alpha_self = density_operation_difference(densities['alpha'], densities['alpha'], exchange)
        alpha_to_beta = density_operation_difference(densities['alpha'], densities['beta'], exchange)
        self.assertEqual(total['projected_direct_hs_difference'], 0.)
        self.assertAlmostEqual(alpha_self['full_space_hs_difference_squared_signed'], 2.)
        self.assertAlmostEqual(spin['full_space_hs_difference_squared_signed'], 8.)
        self.assertEqual(alpha_to_beta['projected_direct_hs_difference'], 0.)

    def test_near_equal_density_retains_direct_difference_and_cancellation_scale(self):
        angle = 1e-9
        operator = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        occupied = np.diag([1., 0.])
        measured = density_operation_difference(occupied, occupied, operator)
        self.assertAlmostEqual(measured['projected_direct_hs_difference']/(np.sqrt(2)*angle), 1., places=12)
        self.assertGreater(measured['arithmetic_squared_scale_estimate'], 2*angle**2)
        self.assertFalse(measured['strict_equality_claimed'])
        zeros = density_operation_difference(np.zeros((2, 2)), np.zeros((2, 2)), operator)
        self.assertEqual(zeros['density_status'], 'both_operators_exactly_zero')
        self.assertIsNone(zeros['full_space_relative_norm_estimate'])
        self.assertIsNone(zeros['projected_direct_relative_difference'])

    def test_nonorthogonal_ao_change_preserves_density_metrics(self):
        original = np.eye(4)[:, :3]
        beta = original @ np.array([[.8, .6, 0.], [-.6, .8, 0.], [0., 0., -1.]])
        overlap = np.eye(4)
        transform = np.array([[1.2, .3, 0., 0.], [0., .7, .1, 0.], [.2, 0., 1.1, .2], [0., .1, 0., .9]])
        inverse_ao = np.linalg.inv(transform)
        changed_s = transform.T @ overlap @ transform
        states = prepare_spin_densities(original, beta, overlap, [1., 1., 0.], [1., 0., 0.], shared_coefficients=False)
        changed = prepare_spin_densities(inverse_ao @ original, inverse_ao @ beta, changed_s,
                                          [1., 1., 0.], [1., 0., 0.], shared_coefficients=False)
        for name in states['densities']:
            np.testing.assert_allclose(states['densities'][name], changed['densities'][name], atol=2e-14)

    def test_truncated_beta_space_error_bound_and_shared_zero_are_explicit(self):
        alpha = np.eye(4)[:, :2]
        beta = np.eye(4)[:, [0, 2]]
        states = prepare_spin_densities(alpha, beta, np.eye(4), [1., 0.], [0., 1.], shared_coefficients=False)
        self.assertAlmostEqual(states['beta_relative_column_residual'], 1.)
        # The omitted occupied beta density is one normalized rank-one operator.
        self.assertAlmostEqual(states['beta_density_projection_error_bound_hs'], 1.)
        self.assertFalse(states['exact_zero_spin_from_shared_occupations'])
        shared = prepare_spin_densities(alpha, alpha, np.eye(4), [1., 0.], [1., 0.], shared_coefficients=True)
        self.assertTrue(shared['exact_zero_spin_from_shared_occupations'])
        np.testing.assert_array_equal(shared['densities']['spin'], np.zeros((2, 2)))
        with self.assertRaises(ValueError):
            prepare_spin_densities(alpha, beta, np.eye(4), [1., 0.], [1., 0.], shared_coefficients=True)

    def test_chained_energy_window_records_full_span_and_does_not_claim_irrep(self):
        blocks = energy_connected_blocks([-.3, -.5, -.5+7e-6, -.5+14e-6], 1e-5)
        self.assertEqual(blocks[0]['members_zero_based'], [1, 2, 3])
        self.assertGreater(blocks[0]['energy_range_hartree'][1]-blocks[0]['energy_range_hartree'][0], 1e-5)
        matrix = np.array([[0., -1.], [1., 0.]])
        whole = subspace_measurement(matrix, np.eye(2))
        partial = subspace_measurement(matrix[:1, :1], np.eye(1))
        self.assertAlmostEqual(whole['character'], partial['character'])
        self.assertAlmostEqual(whole['subspace_leakage_squared_signed'], 0.)
        self.assertAlmostEqual(partial['subspace_leakage_squared_signed'], 1.)
        self.assertFalse(whole['exact_irrep_claimed'])


if __name__ == '__main__':
    unittest.main()
