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

from gaussian_source_events import extract_events, review_final_frequency_chain


class SourceEvidenceTests(unittest.TestCase):
    def review(self, text):
        return review_final_frequency_chain(extract_events(text.splitlines()), 1e-8)

    def test_changed_wavefunction_after_frequency_is_not_a_final_state_frequency(self):
        result = self.review(''' SCF Done: E(RPBE1PBE) = -10.125 A.U.
 Frequencies --  123.0 456.0 789.0
 Optimization completed.
 The wavefunction has an RHF -> UHF instability.
 SCF Done: E(UPBE1PBE) = -10.130 a.u.
 The wavefunction is stable under the perturbations considered.
 Normal termination of Gaussian 16 at Sun Sep 6
''')
        self.assertEqual(result['status'], 'post_frequency_instability_and_scf_state_change')
        self.assertAlmostEqual(result['energy_delta_hartree'], -0.005)
        self.assertEqual(result['last_frequency']['line'], 2)
        self.assertFalse(result['physical_pass'])

    def test_repeat_frequency_after_stability_removes_that_specific_gap(self):
        result = self.review(''' SCF Done: E(RHF) = -1.0 A.U.
 Frequencies -- 10.0
 The wavefunction has an internal instability.
 SCF Done: E(UHF) = -1.1 A.U.
 The wavefunction is stable under the perturbations considered.
 SCF Done: E(UHF) = -1.2 A.U.
 Optimization completed.
 Frequencies -- 20.0
 SCF Done: E(UHF) = -1.2 A.U.
 The wavefunction is already stable.
''')
        self.assertEqual(result['status'], 'no_explicit_post_frequency_state_change_at_log_precision')
        self.assertFalse(result['physical_pass'])

    def test_energy_change_without_instability_is_an_ambiguous_context_change(self):
        result = self.review(''' SCF Done: E(RHF) = -1.0 A.U.
 Frequencies -- 10.0
 SCF Done: E(RHF) = -1.1 A.U.
''')
        self.assertEqual(result['status'], 'post_frequency_scf_state_change_needs_context')
        self.assertEqual(result['later_instability_reports'], [])

    def test_missing_evidence_and_archive_text_are_not_zero_or_pass(self):
        self.assertEqual(self.review('')['status'], 'missing_scf_evidence')
        self.assertEqual(self.review('SCF Done: E(RHF) = -1.0')['status'], 'missing_frequency_evidence')
        self.assertEqual(self.review('Frequencies -- 12.0\nSCF Done: E(RHF) = -1.0')['status'],
                         'frequency_state_not_identified')
        self.assertEqual(extract_events([' archive\\SCF Done: E(RHF) = -1.0'])['events'], [])

    def test_fortran_numbers_spin_and_imaginary_frequencies_preserved(self):
        data = extract_events([' SCF Done: E(UHF) = -1.2D+01 A.U.',
                               ' Frequencies -- -1.5D+01 2.0E+02',
                               ' S**2 before annihilation 0.3154, after 0.0030'])
        self.assertEqual(data['events'][0]['energy_hartree'], -12.0)
        self.assertEqual(data['events'][1]['values_cm_inverse'], [-15.0, 200.0])
        self.assertEqual(data['events'][2]['before_annihilation'], 0.3154)

    def test_invalid_values_fail_instead_of_creating_false_evidence(self):
        for line in ['SCF Done: E(UHF) = NaN', 'Frequencies --',
                     'Frequencies -- 1.0 ******', 'SCF Done: E(UHF) = 1.0E999']:
            with self.subTest(line=line), self.assertRaises(ValueError):
                extract_events([line])


if __name__ == '__main__':
    unittest.main()
