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

import ctypes as ct
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest

from validation_process import physical_core_masks, run_tree
from gaussian_rebuild import input_text, source_identity_check


@unittest.skipUnless(os.name == "nt", "Windows Job Object contract")
class ProcessTreeTests(unittest.TestCase):
    def test_real_stdin_file_and_console_output_are_preserved(self):
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            input_path = directory / 'commands.txt'
            payload = b'-1 -6.000000 -6.000000 -6.000000\nq\n'
            input_path.write_bytes(payload)
            command = "import sys; print(sys.stdin.buffer.read().hex()); print('stderr-preserved',file=sys.stderr)"
            result = run_tree([sys.executable, '-c', command], directory, dict(os.environ),
                              physical_core_masks(1)[0], 1, 30, affinity=False, stdin_path=input_path)
            self.assertEqual(result['exit_code'], 0)
            self.assertEqual(Path(result['stdin_path']), input_path.resolve())
            output = (directory/'launcher.log').read_text()
            self.assertIn(payload.hex(), output)
            self.assertIn('stderr-preserved', output)
            self.assertEqual(input_path.read_bytes(), payload)

    def test_cpu_quota_applies_without_affinity_restriction(self):
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            child = (
                "import ctypes as c,json; from ctypes import wintypes as w; "
                "k=c.WinDLL('kernel32',use_last_error=True); "
                "k.QueryInformationJobObject.argtypes=[w.HANDLE,c.c_int,c.c_void_p,w.DWORD,c.c_void_p]; "
                "r=(w.DWORD*2)(); "
                "assert k.QueryInformationJobObject(None,15,c.byref(r),c.sizeof(r),None); "
                "open('rate.json','w').write(json.dumps(list(r)))"
            )
            result = run_tree([sys.executable, "-c", child], directory, dict(os.environ),
                              physical_core_masks(1)[0], 1, 15, affinity=False)
            flags, rate = json.loads((directory / "rate.json").read_text())
            self.assertEqual(flags, 5)
            self.assertEqual(rate, result["cpu_rate_hard_cap"])
            self.assertLessEqual(rate, 10000 / os.cpu_count())
            self.assertIsNone(result["cpu_mask"])

    def test_waits_for_descendants_and_inherits_affinity(self):
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            result = directory / "child.json"
            child = (
                "import ctypes as c,json,time; from ctypes import wintypes as w; "
                "time.sleep(1); k=c.WinDLL('kernel32',use_last_error=True); "
                "k.GetCurrentProcess.restype=w.HANDLE; "
                "k.GetProcessAffinityMask.argtypes=[w.HANDLE,c.c_void_p,c.c_void_p]; "
                "a=c.c_size_t(); b=c.c_size_t(); "
                "assert k.GetProcessAffinityMask(k.GetCurrentProcess(),c.byref(a),c.byref(b)); "
                "open('child.json','w').write(json.dumps({'mask':a.value}))"
            )
            parent = "import subprocess,sys; subprocess.Popen([sys.executable,'-c'," + repr(child) + "])"
            mask = physical_core_masks(1)[0]
            metadata = run_tree([sys.executable, "-c", parent], directory, dict(os.environ), mask, 1, 15)
            self.assertTrue(result.exists(), "Must await child even after launcher exits")
            self.assertEqual(json.loads(result.read_text())["mask"], mask)
            self.assertGreaterEqual(metadata["tree_process_count"], 2)
            self.assertGreaterEqual(metadata["wall_seconds"], 1)
            self.assertFalse(metadata["timed_out"])

    def test_timeout_terminates_descendant(self):
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            child = "import time; time.sleep(2); open('escaped.txt','w').write('bad')"
            parent = "import subprocess,sys; subprocess.Popen([sys.executable,'-c'," + repr(child) + "])"
            metadata = run_tree([sys.executable, "-c", parent], directory, dict(os.environ),
                                physical_core_masks(1)[0], 1, 0.25)
            time.sleep(2)
            self.assertTrue(metadata["timed_out"])
            self.assertFalse((directory / "escaped.txt").exists())


class ReproductionTests(unittest.TestCase):
    def setUp(self):
        self.case = {"case_id": "OLD-test", "original_input": {
            "basis": "GenECP", "extra": "H 0\nS 1 1.0\n 1.0 1.0\n****",
            "route": "#p UPBE1PBE/GenECP Opt=(Tight,CalcFC) Freq Guess=Mix NoSymm 5D 7F EmpiricalDispersion=GD3BJ"},
            "fchk_identity": {"method": "UPBE1PBE", "Charge": 0, "Multiplicity": 2,
                "Atomic numbers": [1], "Current cartesian coordinates": [0.1, 0.2, 0.3],
                "Number of alpha electrons": 1, "Number of beta electrons": 0,
                "Shell types": [0], "Shell to atom map": [1],
                "Number of primitives per shell": [1], "Number of basis functions": 1,
                "Nuclear charges": [1], "Primitive exponents": [1.0],
                "Contraction coefficients": [1.0]}}

    def test_preserves_basis_state_dispersion_and_fresh_guess(self):
        text = input_text(self.case)
        self.assertIn("UPBE1PBE/GenECP SP Units=Bohr", text)
        self.assertIn("Guess=Mix", text)
        self.assertIn("EmpiricalDispersion=GD3BJ", text)
        self.assertIn("%nprocshared=4", text)
        self.assertIn("%mem=24GB", text)
        self.assertNotIn("Opt=", text)
        self.assertNotIn("Freq", text)
        self.assertNotIn("Read", text)
        self.assertIn(self.case["original_input"]["extra"], text)
        self.assertIn("MaxCycle=2048", input_text(self.case, retry=True))

    def test_rejects_uninterpreted_extra_input(self):
        self.case["original_input"]["basis"] = "def2SVP"
        with self.assertRaises(ValueError):
            input_text(self.case)

    def test_identity_detects_ecp_and_shell_change(self):
        expected = self.case["fchk_identity"]
        self.assertTrue(source_identity_check(self.case, dict(expected))["input_identity_consistent"])
        changed = dict(expected, **{"Nuclear charges": [0]})
        self.assertIn("Nuclear charges", source_identity_check(self.case, changed)["mismatches"])
        changed = dict(expected, **{"Shell types": [-2]})
        self.assertIn("Shell types", source_identity_check(self.case, changed)["mismatches"])


if __name__ == "__main__":
    unittest.main()
