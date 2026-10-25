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

from gaussian_reference_inputs import make_followup_input, make_initial_input, reference_purpose


class ReferenceInputTests(unittest.TestCase):
    def setUp(self):
        self.case = {"purpose": {"kind": "scan_frame", "reference_geometry": "fixed"},
                     "original_input": {"route": "#p PBE1PBE/def2SVP SP"},
                     "fchk_identity": {"Number of atoms": 2, "Atomic numbers": [1, 1],
                         "Current cartesian coordinates": [0, 0, 0, 0, 0, 3],
                         "Charge": 0, "Multiplicity": 1}}
        self.reference = {"case_id": "test", "purpose": reference_purpose(self.case),
                          "gaussian_basis_keyword": "Gen", "shell_flags": "5D 7F 9G"}

    def test_scan_stays_fixed_and_contains_force_evidence(self):
        result = make_initial_input(self.reference, self.case, "H 0\nS 1 1.0\n1 1\n****")
        self.assertIn(" Force ", result)
        self.assertNotIn("Opt=", result)
        self.assertIn("Integral=SuperFineGrid", result)
        self.assertNotIn("Guess=Read", result)
        repaired = make_followup_input(self.reference, "UPBE1PBE", "relax_after_stability")
        self.assertIn("UPBE1PBE/ChkBasis Force", repaired)
        self.assertNotIn("Opt=", repaired)

    def test_changed_unrestricted_state_is_retained_and_refrequenced(self):
        self.reference["purpose"] = {"reference_geometry": "optimize"}
        result = make_followup_input(self.reference, "UPBE1PBE", "relax_after_stability")
        self.assertIn("UPBE1PBE/ChkBasis Opt=", result)
        self.assertIn(" Freq ", result)
        self.assertIn("Guess=Read", result)

    def test_original_optimize_takes_precedence_over_name_classification(self):
        self.case["purpose"] = {"kind": "specified_geometry", "reference_geometry": "fixed"}
        self.case["original_input"]["route"] = "#p PBE1PBE/def2SVP Opt=(Tight,CalcFC) Freq"
        self.assertEqual(reference_purpose(self.case)["reference_geometry"], "optimize")

    def test_atom_has_no_internal_vibrations(self):
        self.case["fchk_identity"]["Number of atoms"] = 1
        purpose = reference_purpose(self.case)
        self.assertEqual(purpose["reference_geometry"], "single_atom")
        self.assertIn("no_internal_nuclear_coordinates", purpose["vibration_applicability"])


if __name__ == "__main__":
    unittest.main()
