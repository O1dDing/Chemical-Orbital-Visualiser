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

from decimal import Decimal, localcontext
import struct
import unittest

from coefficient_precision import measure_coefficient_rounding, np


class CoefficientPrecisionTests(unittest.TestCase):
    def test_nonorthogonal_cancellation_against_decimal_polynomial_integral(self):
        # This analytic two-function integral is computed independently with
        # Decimal arithmetic and struct's binary32 conversion, not Cholesky.
        c = np.array([[100000.003], [-99999.997]])
        s = np.array([[1.0, .99999999], [.99999999, 1.0]])
        untouched = c.copy()
        measured = measure_coefficient_rounding(c, s)
        with localcontext() as context:
            context.prec = 60
            a, b, overlap = map(Decimal.from_float, (float(c[0, 0]), float(c[1, 0]), float(s[0, 1])))
            qa = Decimal.from_float(struct.unpack('<f', struct.pack('<f', float(a)))[0])
            qb = Decimal.from_float(struct.unpack('<f', struct.pack('<f', float(b)))[0])
            da, db = qa-a, qb-b
            norm_squared = a*a + b*b + 2*a*b*overlap
            error_squared = da*da + db*db + 2*da*db*overlap
            expected = float((error_squared/norm_squared).sqrt())
        self.assertAlmostEqual(float(measured['relative_errors'][0])/expected, 1.0, places=8)
        self.assertTrue(np.array_equal(c, untouched))
        self.assertGreater(expected, 1e-5)

    def test_exact_binary32_values_and_rectangular_matrix(self):
        c = np.array([[1.0], [.5], [-2.0]])
        measured = measure_coefficient_rounding(c, np.diag([2., 3., 4.]))
        self.assertEqual(float(measured['error_norms'][0]), 0.0)
        self.assertAlmostEqual(float(measured['original_norms'][0]), np.sqrt(18.75))
        self.assertEqual(measured['relative_errors'].shape, (1,))

    def test_missing_or_invalid_metric_evidence_is_rejected(self):
        fixtures = [(np.zeros((2, 1)), np.eye(2)),
                    (np.ones((2, 1)), np.array([[1., 1.], [1., 1.]])),
                    (np.ones((2, 1)), np.array([[1., .1], [.2, 1.]])),
                    (np.ones((2, 1)), np.eye(3)),
                    (np.array([[np.nan]]), np.eye(1))]
        for c, s in fixtures:
            with self.subTest(c=c.tolist(), s=s.tolist()), self.assertRaises((ValueError, np.linalg.LinAlgError)):
                measure_coefficient_rounding(c, s)


if __name__ == '__main__':
    unittest.main()
