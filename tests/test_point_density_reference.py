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

import unittest

from point_density_reference import density_channels, density_points, sampled_difference, spin_densities_at_points, np
from gbasis.contractions import GeneralizedContractionShell


class PointDensityReferenceTests(unittest.TestCase):
    def shell(self, angular, center=None):
        return GeneralizedContractionShell(angular, np.zeros(3) if center is None else center,
            np.array([1.0]), np.array([0.8]), 'cartesian')

    def test_normalized_s_density_matches_closed_form_and_zero_beta(self):
        points = np.random.default_rng(2718).normal(size=(57, 3))
        c = np.ones((1, 1))
        actual = spin_densities_at_points([self.shell(0)], points, c, c, [1], [0], shared_coefficients=True, chunk_points=11)
        expected = (1.6/np.pi)**1.5*np.exp(-1.6*np.sum(points*points, axis=1))
        np.testing.assert_allclose(actual[0], expected, rtol=3e-15, atol=1e-16)
        self.assertFalse(np.any(actual[1]))
        zeros = sampled_difference(actual[1], actual[1])
        self.assertIsNone(zeros['relative_sample_norm_difference'])
        self.assertTrue(zeros['source_samples_exactly_zero'])

    def test_p_orbitals_exchange_spins_with_invariant_total_and_flipped_spin(self):
        shell = self.shell(1)
        powers = shell.angmom_components_cart
        ix = next(i for i, row in enumerate(powers) if np.array_equal(row, [1, 0, 0]))
        iy = next(i for i, row in enumerate(powers) if np.array_equal(row, [0, 1, 0]))
        ca, cb = np.eye(3)[:, [ix]], np.eye(3)[:, [iy]]
        points = np.array([[.2, .5, .1], [-.8, .3, .2], [.9, -.4, .7], [.1, .8, -.2]])
        r = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        original = spin_densities_at_points([shell], points, ca, cb, [1], [1])
        moved = spin_densities_at_points([shell], points @ r, ca, cb, [1], [1])
        normalizer = 3.2*(1.6/np.pi)**1.5*np.exp(-1.6*np.sum(points*points, axis=1))
        np.testing.assert_allclose(original[0], normalizer*points[:, 0]**2, atol=1e-16)
        np.testing.assert_allclose(original[1], normalizer*points[:, 1]**2, atol=1e-16)
        np.testing.assert_allclose(moved, original[::-1], atol=2e-16)
        a, b = density_channels(original), density_channels(moved)
        np.testing.assert_allclose(a['total'], b['total'], atol=2e-16)
        np.testing.assert_allclose(a['spin'], -b['spin'], atol=2e-16)
        self.assertGreater(sampled_difference(a['alpha'], b['alpha'])['relative_sample_norm_difference'], .5)
        self.assertAlmostEqual(sampled_difference(a['spin'], b['spin'])['relative_sample_norm_difference'], 2.)

    def test_sign_and_equal_occupation_gauge_preserve_point_density(self):
        rng = np.random.default_rng(1618)
        points = rng.normal(size=(31, 3))
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        c = np.eye(3)
        before = c.copy()
        a = spin_densities_at_points([self.shell(1)], points, c, c, [1, 1, 1], [.5, .5, .5], shared_coefficients=True)
        b = spin_densities_at_points([self.shell(1)], points, q, -q, [1, 1, 1], [.5, .5, .5])
        np.testing.assert_allclose(a, b, atol=3e-16, rtol=3e-15)
        np.testing.assert_array_equal(c, before)

    def test_rigid_translation_uses_the_same_point_identities(self):
        rng = np.random.default_rng(99)
        points, meta = density_points(np.array([[-.4, .2, 0.], [.5, .3, -.1]]), grid_shape=(3, 4, 5), atom_radii=(.2, .7))
        self.assertEqual(len(points), 60+2+24)
        self.assertEqual(meta['strata'][-1]['stop'], len(points))
        shift = rng.normal(size=3)
        a = spin_densities_at_points([self.shell(0)], points, [[1]], [[1]], [1], [1], shared_coefficients=True)
        b = spin_densities_at_points([self.shell(0, shift)], points+shift, [[1]], [[1]], [1], [1], shared_coefficients=True)
        np.testing.assert_allclose(a, b, atol=5e-16)

    def test_node_samples_are_not_promoted_to_all_space_zero(self):
        shell = self.shell(1)
        ix = next(i for i, row in enumerate(shell.angmom_components_cart) if np.array_equal(row, [1, 0, 0]))
        c = np.eye(3)[:, [ix]]
        nodes = np.array([[0., .3, .2], [0., -.8, 1.0]])
        rho = spin_densities_at_points([shell], nodes, c, c, [1], [0])
        result = sampled_difference(rho[0], rho[0])
        self.assertTrue(result['source_samples_exactly_zero'])
        self.assertIsNone(result['relative_sample_norm_difference'])
        self.assertIn('no continuous-space', result['norm_scope'])
        off_node = spin_densities_at_points([shell], nodes+np.array([.2, 0, 0]), c, c, [1], [0])
        self.assertTrue(np.all(off_node[0] > 0))

    def test_invalid_occupations_and_coefficient_provenance_are_rejected(self):
        with self.assertRaises(ValueError):
            spin_densities_at_points([self.shell(0)], [[0, 0, 0]], [[1]], [[1]], [-1], [0])
        with self.assertRaises(ValueError):
            spin_densities_at_points([self.shell(0)], [[0, 0, 0]], [[1]], [[-1]], [1], [1], shared_coefficients=True)
        with self.assertRaises(ValueError):
            density_points([[0, 0, 0]], grid_shape=(3.5, 3, 3))


if __name__ == '__main__':
    unittest.main()
