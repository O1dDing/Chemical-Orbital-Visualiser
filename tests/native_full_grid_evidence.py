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

"""Check integrity and frame association of native volume_full evidence.

This is evidence transport validation. It is not a wavefunction accuracy or
UI correctness verdict; Gaussian/GBasis and image checks are still required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct


def file_hash(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024*1024), b""):
            value.update(data)
    return value.hexdigest()


def inspect_full_grid(metadata_path):
    metadata_path = Path(metadata_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    full = metadata["full_grid"]
    root = metadata_path.parent.resolve()
    binary = (root / full["file"]).resolve()
    if binary.parent != root:
        raise ValueError("Texture evidence must remain inside its capture directory")
    count = math.prod(metadata[axis] for axis in ("nx", "ny", "nz"))
    if count <= 0 or any(metadata[axis] < 2 for axis in ("nx", "ny", "nz")):
        raise ValueError("Invalid captured grid dimensions")
    if full["scalar_type"] != "IEEE754-float32-little-endian":
        raise ValueError("Unknown complete-grid scalar representation")
    if full["point_count"] != count or full["byte_count"] != count*4 or binary.stat().st_size != count*4:
        raise ValueError("Complete texture file is missing points or disagrees with its metadata")
    samples = metadata["samples"]
    seen = set()
    with binary.open("rb") as stream:
        for index, expected in samples:
            if index in seen or not 0 <= index < count:
                raise ValueError("Invalid or duplicated sample index")
            seen.add(index)
            stream.seek(index*4)
            actual, = struct.unpack("<f", stream.read(4))
            if not math.isfinite(actual) or actual != expected:
                raise ValueError("Sparse sample values do not match the complete texture")
    prefix = metadata_path.name.removesuffix(".volume.json")
    ui_path, bitmap_path = root/(prefix+".ui.json"), root/(prefix+".bmp")
    ui = json.loads(ui_path.read_text(encoding="utf-8"))
    state = ui["state"]
    if state["frame"] != metadata["frame"] or state["rendered_mo"] != metadata["rendered_mo"]:
        raise ValueError("Screen/UI capture belongs to a different frame or rendered orbital")
    if state["rendered_generation"] != metadata["generation"] or state["volume_generation"] != metadata["generation"]:
        raise ValueError("Renderer generation differs between screen/UI and texture evidence")
    with bitmap_path.open("rb") as stream:
        header = stream.read(54)
    if len(header) != 54 or header[:2] != b"BM":
        raise ValueError("Corresponding framebuffer is missing or is not a BMP")
    recorded_size, = struct.unpack_from("<I", header, 2)
    width, height = struct.unpack_from("<ii", header, 18)
    if recorded_size != bitmap_path.stat().st_size or width <= 0 or height <= 0:
        raise ValueError("Framebuffer file is truncated or has invalid dimensions")
    # Cross-check the independent per-frame log emitted by end_frame.
    matches = []
    with (root/"frames.jsonl").open(encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            if row["frame"] == metadata["frame"]:
                matches.append(row)
    if len(matches) != 1 or matches[0] != state:
        raise ValueError("A unique corresponding frame state is not present in the frame log")
    return {"status": "evidence_transport_pass", "counts_as_scientific_pass": False,
            "metadata": str(metadata_path), "metadata_sha256": file_hash(metadata_path),
            "texture": str(binary), "texture_sha256": file_hash(binary),
            "framebuffer": str(bitmap_path), "framebuffer_sha256": file_hash(bitmap_path),
            "ui_sha256": file_hash(ui_path), "frame": metadata["frame"],
            "rendered_mo": metadata["rendered_mo"], "generation": metadata["generation"],
            "grid_shape": [metadata[axis] for axis in ("nx", "ny", "nz")],
            "point_count": count, "matched_sparse_samples": len(seen),
            "scope": "Complete-file size, sparse-value agreement and same-frame association; independent numerical/image checks remain required"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = inspect_full_grid(arguments.metadata)
    encoded = json.dumps(result, indent=2)
    if arguments.output:
        arguments.output.write_text(encoded, encoding="utf-8")
    print(encoded)
