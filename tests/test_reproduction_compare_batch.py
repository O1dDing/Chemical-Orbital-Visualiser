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

import argparse
import json
from pathlib import Path
import tempfile
import unittest

from reproduction_compare_batch import require_barrier


class ComparisonBarrierTests(unittest.TestCase):
    def test_missing_partial_and_wrong_version_batches_cannot_start_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            args = argparse.Namespace(output=output)
            identity = {"case_set_identity": "frozen_original_cases"}
            self.assertFalse(require_barrier(args, identity))
            (output/"reproduction").mkdir()
            path = output/"reproduction/batch.json"
            for terminal, count in ((False, 273), (True, 272)):
                path.write_text(json.dumps(dict(identity, all_terminal=terminal, case_count=count)))
                self.assertFalse(require_barrier(args, identity))
            path.write_text(json.dumps({"case_set_identity": "other_cases", "all_terminal": True, "case_count": 273}))
            with self.assertRaisesRegex(ValueError, "case sets differ"):
                require_barrier(args, identity)
            path.write_text(json.dumps(dict(identity, all_terminal=True, case_count=273)))
            self.assertTrue(require_barrier(args, identity))


if __name__ == "__main__":
    unittest.main()
