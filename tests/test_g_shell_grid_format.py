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

from g_shell_fixture_batch import grid_spec


class GaussianGridSyntaxTests(unittest.TestCase):
    def test_integer_valued_coordinates_are_serialized_as_real_tokens(self):
        text = grid_spec({'grid': {'origin_bohr': [-6, -6, -6], 'step_bohr': .375, 'shape': [33, 33, 33]}})
        rows = [row.split() for row in text.splitlines()]
        self.assertEqual([int(row[0]) for row in rows], [-1, -33, 33, 33])
        for row in rows:
            self.assertEqual(len(row), 4)
            self.assertTrue(all('.' in token for token in row[1:]))
        self.assertEqual([float(v) for v in rows[0][1:]], [-6., -6., -6.])
        self.assertEqual([[float(v) for v in row[1:]] for row in rows[1:]],
                         [[.375, 0., 0.], [0., .375, 0.], [0., 0., .375]])


if __name__ == '__main__':
    unittest.main()
