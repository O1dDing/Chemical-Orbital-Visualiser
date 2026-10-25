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

"""Independent analytic shapes, group orbits and polynomial/grid cross-checks."""
import json
import unittest

from nuclear_symmetry_reference import (group_closure, match_operation, nuclear_isometries,
                                        np, pullback_basis, subspace_operation)
from gbasis.contractions import GeneralizedContractionShell
from gbasis.evals.eval import evaluate_basis
from scipy.spatial.transform import Rotation


class NuclearSymmetryReferenceTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(928301)
        self.rotation = Rotation.random(random_state=self.rng).as_matrix()

    def check_shape(self, xyz, labels, expected, count):
        result = nuclear_isometries(xyz, labels, geometry_tolerance=1e-9)
        result = json.loads(json.dumps(result, allow_nan=False))
        self.assertEqual(result['group_under_tolerance'], expected)
        self.assertEqual(result['operation_count'], count)
        self.assertTrue(result['closure']['closed_under_tolerance'])
        self.assertLess(result['closure']['max_product_residual'], 2e-12)
        permutation = self.rng.permutation(len(xyz))
        moved = np.asarray(xyz)[permutation] @ self.rotation+np.array([2.4, -3.7, .9])
        moved_labels = [labels[i] for i in permutation]
        transformed = nuclear_isometries(moved, moved_labels, geometry_tolerance=1e-9)
        self.assertEqual(transformed['group_under_tolerance'], expected)
        self.assertEqual(transformed['operation_count'], count)
        return result

    def test_regular_shapes_and_full_nuclear_scope(self):
        tetra = np.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]])
        self.check_shape(tetra, [1]*4, 'Td', 24)
        octa = np.concatenate((np.eye(3), -np.eye(3)))
        self.check_shape(octa, [1]*6, 'Oh', 48)
        phi = (1+np.sqrt(5))/2
        icosa = np.array([(0, a, b*phi) for a in (-1, 1) for b in (-1, 1)]
                        + [(a, b*phi, 0) for a in (-1, 1) for b in (-1, 1)]
                        + [(b*phi, 0, a) for a in (-1, 1) for b in (-1, 1)])
        self.check_shape(icosa, [1]*12, 'Ih', 120)
        self.check_shape(np.array([[1., 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]), [1]*4, 'D4h', 16)
        # One atom outside an ideal local octahedron breaks the whole-system group.
        whole = np.vstack((octa, [.27, .43, 2.13]))
        self.check_shape(whole, [1]*6+[9], 'C1', 1)

    def test_prism_and_antiprism_are_distinct_from_shape_names(self):
        for n in (3, 5, 7):
            angles = np.arange(n)*2*np.pi/n
            upper = np.column_stack((np.cos(angles), np.sin(angles), np.full(n, .61)))
            for offset, expected in ((0., f'D{n}h'), (np.pi/n, f'D{n}d')):
                lower = np.column_stack((np.cos(angles+offset), np.sin(angles+offset), np.full(n, -.61)))
                self.check_shape(np.vstack((upper, lower)), [1]*(2*n), expected, 4*n)

    def orbit(self, matrices):
        # Two generic orbits with distinct labels remove accidental symmetries.
        first = np.einsum('i,gij->gj', [1.17, .32, .83], matrices)
        second = np.einsum('i,gij->gj', [.48, 1.39, -.57], matrices)
        return np.vstack((first, second)), [1]*len(first)+[2]*len(second)

    def test_general_finite_families_from_independent_scipy_generators(self):
        for name in ('C2', 'C3', 'C5', 'D2', 'D3', 'D4', 'T', 'O', 'I'):
            with self.subTest(name=name):
                matrices = Rotation.create_group(name).as_matrix()
                xyz, labels = self.orbit(matrices)
                self.check_shape(xyz, labels, name, len(matrices))
        horizontal = np.diag([1., 1., -1.])
        vertical = np.diag([1., -1., 1.])
        for n in (2, 3, 4, 5):
            proper = Rotation.create_group(f'C{n}').as_matrix()
            for reflector, suffix in ((horizontal, 'h'), (vertical, 'v')):
                matrices = np.concatenate((proper, proper @ reflector))
                xyz, labels = self.orbit(matrices)
                self.check_shape(xyz, labels, f'C{n}{suffix}', 2*n)
        for n in (4, 6, 8):
            generator = Rotation.from_rotvec([0, 0, 2*np.pi/n]).as_matrix() @ horizontal
            matrices = np.array([np.linalg.matrix_power(generator, k) for k in range(n)])
            xyz, labels = self.orbit(matrices)
            self.check_shape(xyz, labels, f'S{n}', n)
        tetra = Rotation.create_group('T').as_matrix()
        xyz, labels = self.orbit(np.concatenate((tetra, -tetra)))
        self.check_shape(xyz, labels, 'Th', 24)

    def test_continuous_atomic_and_linear_freedom_and_element_identity(self):
        for xyz, labels, expected in (([[1., 2., 3.]], [8], 'Kh'),
                                      ([[-1., 0, 0], [1., 0, 0]], [1, 1], 'Dinfh'),
                                      ([[-1., 0, 0], [1., 0, 0]], [1, 9], 'Cinfv')):
            result = nuclear_isometries(xyz, labels, geometry_tolerance=1e-9)
            self.assertEqual(result['group_under_tolerance'], expected)
            self.assertEqual(result['operations'], [])
        result = match_operation([[-1., 0, 0], [1., 0, 0]], [1, 9], -np.eye(3))
        self.assertEqual(result['permutation_zero_based'], [0, 1])
        self.assertAlmostEqual(result['max_displacement_bohr'], 2.)
        unresolved = nuclear_isometries([[0., 0, 0], [1e-12, 0, 0]], [1, 1], geometry_tolerance=1e-9)
        self.assertEqual(unresolved['kind'], 'coincident_or_unresolved_nuclei')
        self.assertIsNone(unresolved['group_under_tolerance'])

    def test_bijection_closure_and_distorted_geometry_cannot_silently_pass(self):
        xyz = np.array([[0., 0, 0], [.05, 0, 0], [1.8, 0, 0]])
        result = match_operation(xyz, [1]*3, -np.eye(3), origin=[.1, 0, 0])
        self.assertEqual(sorted(result['permutation_zero_based']), [0, 1, 2])
        self.assertGreater(result['max_displacement_bohr'], 1.)
        c3 = Rotation.from_rotvec([0, 0, 2*np.pi/3]).as_matrix()
        self.assertFalse(group_closure([np.eye(3), c3], 1e-9)['closed_under_tolerance'])
        distorted = np.array([[1., .011, 0], [0, 1.03, 0], [-1., 0, .02], [0, -1., 0]])
        result = nuclear_isometries(distorted, [1]*4, geometry_tolerance=1e-7)
        self.assertEqual(result['group_under_tolerance'], 'C1')
        with self.assertRaises(ValueError):
            nuclear_isometries([[np.nan, 0, 0]], [1], geometry_tolerance=1e-7)
        with self.assertRaises(ValueError):
            match_operation(xyz, [1]*3, np.eye(3)*1.0001)

    def test_almost_collinear_but_resolvable_geometry_has_orthogonal_operations(self):
        xyz = np.array([[-1.1, 0., 0.], [.3, 1e-6, 0.], [1.7, 0., 0.]])
        # This remains a plane at the frozen threshold, not a linear molecule.
        result = self.check_shape(xyz, [1, 8, 9], 'Cs', 2)
        for record in result['operations']:
            matrix = np.asarray(record['matrix'])
            np.testing.assert_allclose(matrix.T @ matrix, np.eye(3), rtol=0, atol=1e-12)

    def shell(self, angular, kind, center):
        return GeneralizedContractionShell(angular, np.asarray(center, dtype=float),
            np.array([[.6, .2], [.4, .8]]), np.array([.8, .25]), kind)

    def test_proper_improper_pullbacks_through_g_match_independent_real_space_values(self):
        origin = np.array([.2, -.4, .7])
        points = self.rng.normal(size=(79, 3))
        for determinant in (1., -1.):
            matrix = self.rotation.copy()
            matrix[:, 0] *= determinant
            for angular in range(5):
                for kind in ('cartesian', 'spherical'):
                    with self.subTest(determinant=determinant, angular=angular, kind=kind):
                        basis = [self.shell(angular, kind, [.6, .3, -.1])]
                        moved, transform = pullback_basis(basis, matrix, origin)
                        actual = evaluate_basis(moved, points).T @ transform
                        expected = evaluate_basis(basis, (points-origin) @ matrix+origin).T
                        np.testing.assert_allclose(actual, expected, atol=4e-13, rtol=1e-11)

    def test_subspace_trace_alone_does_not_hide_leakage_and_is_gauge_invariant(self):
        basis = [GeneralizedContractionShell(1, np.zeros(3), np.array([1.]), np.array([.7]), 'cartesian')]
        c4 = Rotation.from_rotvec([0, 0, np.pi/2]).as_matrix()
        complete = subspace_operation(basis, np.eye(3), c4, [0, 0, 0])
        self.assertAlmostEqual(complete['character'], 1., places=12)
        self.assertLess(complete['projected_unitarity_max_error'], 1e-12)
        partial = subspace_operation(basis, np.eye(3)[:, :1], c4, [0, 0, 0])
        self.assertAlmostEqual(partial['character'], 0., places=12)
        self.assertAlmostEqual(partial['subspace_leakage_squared_signed'], 1., places=12)
        planar = subspace_operation(basis, np.eye(3)[:, :2], c4, [0, 0, 0])
        self.assertAlmostEqual(planar['character'], 0., places=12)
        self.assertLess(abs(planar['subspace_leakage_squared_signed']), 1e-12)
        mixed = np.eye(3)[:, :2] @ np.array([[.6, .8], [-.8, .6]])
        mixed[:, 0] *= -1
        changed = subspace_operation(basis, mixed, c4, [0, 0, 0])
        np.testing.assert_allclose(planar['singular_values'], changed['singular_values'], atol=1e-12)
        self.assertAlmostEqual(planar['character'], changed['character'], places=12)

    def test_cartesian_metric_and_rectangular_space_under_inversion(self):
        for angular, dimension in ((2, 6), (3, 10), (4, 15)):
            basis = [GeneralizedContractionShell(angular, np.zeros(3), np.array([1.]), np.array([.7]), 'cartesian')]
            result = subspace_operation(basis, np.eye(dimension)[:, :3], -np.eye(3), [0, 0, 0])
            self.assertAlmostEqual(result['character'], 3*(-1)**angular, places=11)
            self.assertLess(result['projected_unitarity_max_error'], 1e-11)


if __name__ == '__main__':
    unittest.main()
