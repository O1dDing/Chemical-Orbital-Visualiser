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

import json
from pathlib import Path
import struct
import tempfile
import unittest

from native_full_grid_evidence import inspect_full_grid


class FullGridEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.metadata = {"nx": 2, "ny": 2, "nz": 2, "frame": 12, "generation": 3, "rendered_mo": 7,
            "full_grid": {"file": "probe.volume.f32", "scalar_type": "IEEE754-float32-little-endian",
                          "point_count": 8, "byte_count": 32}, "samples": [[0, 1.25], [7, -.5]]}
        (self.root/"probe.volume.f32").write_bytes(struct.pack("<8f", 1.25, 0, 0, 0, 0, 0, 0, -.5))
        self.state = {"frame": 12, "rendered_mo": 7, "rendered_generation": 3, "volume_generation": 3}
        self.write_state()
        header = bytearray(54)
        header[:2] = b"BM"
        struct.pack_into("<I", header, 2, 70)
        struct.pack_into("<ii", header, 18, 2, 2)
        (self.root/"probe.bmp").write_bytes(header+bytes(16))

    def write_state(self):
        (self.root/"probe.ui.json").write_text(json.dumps({"state": self.state}), encoding="utf-8")
        (self.root/"frames.jsonl").write_text(json.dumps(self.state)+"\n", encoding="utf-8")

    def inspect(self):
        path = self.root/"probe.volume.json"
        path.write_text(json.dumps(self.metadata), encoding="utf-8")
        return inspect_full_grid(path)

    def test_complete_capture_is_only_a_transport_pass(self):
        result = self.inspect()
        self.assertEqual(result["point_count"], 8)
        self.assertFalse(result["counts_as_scientific_pass"])

    def test_truncated_grid_is_rejected(self):
        (self.root/"probe.volume.f32").write_bytes(bytes(28))
        with self.assertRaisesRegex(ValueError, "missing points"):
            self.inspect()

    def test_wrong_frame_or_stale_generation_is_rejected(self):
        self.state["rendered_generation"] = 2
        self.write_state()
        with self.assertRaisesRegex(ValueError, "generation differs"):
            self.inspect()
        self.state["rendered_generation"] = 3
        self.state["frame"] = 11
        self.write_state()
        with self.assertRaisesRegex(ValueError, "different frame"):
            self.inspect()

    def test_sample_mismatch_or_missing_framebuffer_is_rejected(self):
        self.metadata["samples"][0][1] = 2.0
        with self.assertRaisesRegex(ValueError, "sample values"):
            self.inspect()
        self.metadata["samples"][0][1] = 1.25
        (self.root/"probe.bmp").write_bytes(b"BM")
        with self.assertRaisesRegex(ValueError, "framebuffer"):
            self.inspect()

    def test_duplicate_frame_log_is_rejected(self):
        (self.root/"frames.jsonl").write_text((json.dumps(self.state)+"\n")*2, encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unique"):
            self.inspect()


if __name__ == "__main__":
    unittest.main()
