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

"""Mathematical checks for the independent reference, without formal-batch input."""
import unittest

from wavefunction_reference import (aligned_basis, cartesian_rotation, density_difference,
    matched_energy_blocks, np, principal_overlap, rigid_alignment, shell_rotation)
from gbasis.contractions import GeneralizedContractionShell
from gbasis.evals.eval import evaluate_basis
from gbasis.integrals.overlap import overlap_integral


class WavefunctionReferenceTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(314159)
        q, _ = np.linalg.qr(self.rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        self.rotation = q

    def shell(self, angular, kind, center=None):
        return GeneralizedContractionShell(angular, np.zeros(3) if center is None else center,
            np.array([0.6, 0.4]), np.array([0.8, 0.25]), kind)

    def test_all_cartesian_and_pure_shells_through_g_against_independent_evaluation(self):
        # Off-axis, different radii: no polynomial interpolation or MO fitting.
        points = self.rng.normal(size=(117, 3))
        for angular in range(5):
            for kind in ("cartesian", "spherical"):
                with self.subTest(angular=angular, kind=kind):
                    shell = self.shell(angular, kind)
                    transform = shell_rotation(shell, self.rotation)
                    a = evaluate_basis([shell], points).T
                    b = evaluate_basis([shell], points @ self.rotation).T
                    np.testing.assert_allclose(a @ transform, b, atol=2e-13, rtol=2e-12)
                    metric = overlap_integral([shell])
                    np.testing.assert_allclose(transform.T @ metric @ transform, metric, atol=2e-12)

    def test_cartesian_rotation_requires_primitive_normalization_ratios(self):
        shell = self.shell(2, "cartesian")
        transform = shell_rotation(shell, self.rotation)
        # Mixed Cartesian components carry a different angular normalizer.
        self.assertGreater(np.linalg.norm(transform.T @ transform-np.eye(6)), 0.1)
        with self.assertRaises(ValueError):
            cartesian_rotation(shell.angmom_components_cart, shell.angmom_components_cart,
                               self.rotation*1.001)

    def test_nuclear_alignment_recovers_rotation_and_translation(self):
        old = self.rng.normal(size=(9, 3))
        new = old @ self.rotation+np.array([2.1, -3.0, 0.7])
        result = rigid_alignment(old, new)
        np.testing.assert_allclose(result["rotation"], self.rotation, atol=1e-14)
        self.assertLess(result["max_displacement_bohr"], 1e-13)
        self.assertFalse(result["rotation_underdetermined"])

    def test_linear_geometry_records_rotation_ambiguity(self):
        old = np.array([[-1., 0., 0.], [1., 0., 0.]])
        result = rigid_alignment(old, old @ self.rotation)
        self.assertTrue(result["rotation_underdetermined"])
        self.assertEqual(result["coordinate_rank"], 1)
        self.assertAlmostEqual(np.linalg.det(result["rotation"]), 1.0)

    def test_transformed_basis_reproduces_values_in_original_coordinate_frame(self):
        centers = self.rng.normal(size=(5, 3))
        moved = centers @ self.rotation+np.array([2.1, -3., .7])
        alignment = rigid_alignment(centers, moved)
        old_basis = [self.shell(i, "cartesian" if i % 2 else "spherical", centers[i]) for i in range(5)]
        new_basis = [self.shell(i, "cartesian" if i % 2 else "spherical", moved[i]) for i in range(5)]
        aligned, transform = aligned_basis(new_basis, alignment)
        points = self.rng.normal(size=(53, 3))
        a = evaluate_basis(aligned, points).T
        b = evaluate_basis(new_basis, points @ self.rotation+np.array([2.1, -3., .7])).T
        np.testing.assert_allclose(a @ transform, b, atol=3e-13)
        np.testing.assert_allclose(a, evaluate_basis(old_basis, points).T, atol=3e-13)

    def test_signs_permutations_and_degenerate_rotation_preserve_subspaces(self):
        mix, _ = np.linalg.qr(self.rng.normal(size=(3, 3)))
        change = np.eye(6)
        change[1:4, 1:4] = mix
        change[:, 0] *= -1
        change = change[:, [4, 2, 0, 3, 5, 1]]
        old_e = np.array([-2., -.5, -.5, -.5, .3, .7])
        new_e = old_e[[4, 2, 0, 3, 5, 1]]
        full = principal_overlap(change, np.eye(6), np.eye(6))
        self.assertAlmostEqual(full["min_singular_value"], 1.0)
        groups = matched_energy_blocks(change, old_e, new_e, np.eye(6), np.eye(6), 1e-5)
        self.assertEqual(sorted(len(g["old_members_zero_based"]) for g in groups), [1, 1, 1, 3])
        for group in groups:
            self.assertAlmostEqual(group["principal_overlap"]["min_singular_value"], 1.0)
            self.assertAlmostEqual(group["mean_energy_delta_hartree"], 0.0)

    def test_truncated_spaces_and_a_changed_occupied_state_are_visible(self):
        unchanged = principal_overlap(np.eye(5)[:, :4], np.eye(5), np.eye(4))
        self.assertEqual(unchanged["old_rank"], 5)
        self.assertEqual(unchanged["new_rank"], 4)
        self.assertAlmostEqual(unchanged["lost_squared_overlap"], 1.0)
        changed = principal_overlap(np.diag([1., 0.8]), np.eye(2), np.eye(2))
        self.assertAlmostEqual(changed["min_singular_value"], 0.8)
        self.assertAlmostEqual(changed["lost_squared_overlap"], .36)

    def test_density_metric_is_invariant_to_ao_coordinate_change(self):
        matrix = self.rng.normal(size=(6, 6))
        metric = matrix @ matrix.T+np.eye(6)
        density = self.rng.normal(size=(6, 2))
        density = density @ density.T
        change = np.eye(6)+self.rng.normal(size=(6, 6))*.03
        inverse = np.linalg.inv(change)
        transformed_density = inverse @ density @ inverse.T
        new_metric = change.T @ metric @ change
        result = density_difference(density, transformed_density, metric, new_metric, metric @ change)
        self.assertLess(abs(result["hilbert_schmidt_squared_signed"]), result["roundoff_squared_scale"])
        self.assertAlmostEqual(result["old_electron_trace"], result["new_electron_trace"], places=10)


if __name__ == "__main__":
    unittest.main()
