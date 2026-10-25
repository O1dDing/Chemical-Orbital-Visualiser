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

"""Independently parse production PNG/XML and compare against release exports.

Usage: python check_export_attribution.py CURRENT --version VERSION [--baseline BASE]
The fixture directories are produced by cov_export_attribution_smoke.
Uses Python's standard library only; no COV metadata or PNG code is reused.
"""
import argparse
import binascii
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET
import zlib


NAME = "Chemical Orbital Visualiser (COV)"
AUTHOR = "O1dDing"
URL = "https://github.com/O1dDing/Chemical-Orbital-Visualiser"
SVG = "{http://www.w3.org/2000/svg}"
RDF = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
XMP = "{http://ns.adobe.com/xap/1.0/}"
DOAP = "{http://usefulinc.com/ns/doap#}"
FOAF = "{http://xmlns.com/foaf/0.1/}"
KINDS = ("mo", "view", "focus", "aomo")


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def parse_png(path):
    data = path.read_bytes()
    require(data[:8] == b"\x89PNG\r\n\x1a\n", f"{path}: signature")
    offset, chunks, text = 8, [], {}
    while offset < len(data):
        require(offset + 12 <= len(data), f"{path}: truncated chunk")
        size = struct.unpack_from(">I", data, offset)[0]
        end = offset + size + 12
        require(end <= len(data), f"{path}: truncated payload")
        kind, payload = data[offset + 4:offset + 8], data[offset + 8:end - 4]
        crc = struct.unpack_from(">I", data, end - 4)[0]
        require(binascii.crc32(kind + payload) == crc, f"{path}: {kind!r} CRC")
        chunks.append((kind, payload, data[offset:end]))
        if kind == b"iTXt":
            keyword, rest = payload.split(b"\0", 1)
            require(1 <= len(keyword) <= 79, f"{path}: PNG text keyword length")
            require(rest[:2] == b"\0\0", f"{path}: iTXt compression fields")
            language, translated, value = rest[2:].split(b"\0", 2)
            language.decode("ascii")
            translated.decode("utf-8")
            key = keyword.decode("latin-1")
            require(key not in text, f"{path}: duplicate {key}")
            text[key] = value.decode("utf-8")
        offset = end
        if kind == b"IEND":
            break
    require(offset == len(data), f"{path}: trailing bytes")
    require(chunks[0][0] == b"IHDR" and chunks[-1][0] == b"IEND", f"{path}: chunk order")
    require(len(chunks[0][1]) == 13 and not chunks[-1][1], f"{path}: IHDR/IEND size")
    width, height, depth, colour, compression, filtering, interlace = struct.unpack(
        ">IIBBBBB", chunks[0][1])
    require((depth, colour, compression, filtering, interlace) == (8, 6, 0, 0, 0),
            f"{path}: expected non-interlaced RGBA8")
    idat = b"".join(payload for kind, payload, _ in chunks if kind == b"IDAT")
    decoder = zlib.decompressobj()
    raw = decoder.decompress(idat) + decoder.flush()
    require(decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail,
            f"{path}: invalid zlib stream")
    stride = width * 4
    require(len(raw) == height * (stride + 1), f"{path}: pixel stream size")
    pixels, previous = bytearray(), bytearray(stride)
    for y in range(height):
        start = y * (stride + 1)
        method, row = raw[start], bytearray(raw[start + 1:start + stride + 1])
        require(method <= 4, f"{path}: invalid PNG filter")
        for x in range(stride):
            a, b, c = (row[x - 4] if x >= 4 else 0), previous[x], (previous[x - 4] if x >= 4 else 0)
            if method == 1:
                predictor = a
            elif method == 2:
                predictor = b
            elif method == 3:
                predictor = (a + b) // 2
            elif method == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                predictor = a if pa <= pb and pa <= pc else b if pb <= pc else c
            else:
                predictor = 0
            row[x] = (row[x] + predictor) & 255
        pixels.extend(row)
        previous = row
    return chunks, text, (width, height, bytes(pixels)), idat


def check_svg(path, software, version):
    document = ET.parse(path).getroot()
    require(document.tag == SVG + "svg", f"{path}: SVG root")
    metadata = document.findall(SVG + "metadata")
    require(len(metadata) == 1, f"{path}: one metadata element")
    meta = metadata[0]
    require(meta.findtext(f".//{XMP}CreatorTool") == software, f"{path}: generating software")
    project = meta.find(f".//{DOAP}Project")
    require(project is not None and project.get(RDF + "about") == URL, f"{path}: software project")
    require(project.findtext(DOAP + "name") == NAME, f"{path}: software name")
    require(project.find(f"{DOAP}homepage").get(RDF + "resource") == URL, f"{path}: software URL")
    require(project.findtext(f"{DOAP}developer/{FOAF}Person/{FOAF}name") == AUTHOR,
            f"{path}: software developer")
    require(project.findtext(f"{DOAP}release/{DOAP}Version/{DOAP}revision") == version,
            f"{path}: software version")
    rendered = ET.tostring(meta, encoding="unicode")
    require("PRIVATE_INPUT_PATH_SENTINEL" not in rendered, f"{path}: input path in metadata")
    require("http://purl.org/dc/elements/1.1/" not in rendered, f"{path}: artwork attribution asserted")
    data = path.read_bytes()
    require(data.count(b"<metadata>") == 1, f"{path}: metadata serialization")
    return re.sub(rb"<metadata>.*?</metadata>", b"", data, flags=re.S)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    software = NAME + " " + args.version
    expected = {"Software": software, "Software Version": args.version,
                "Software Author": AUTHOR, "Software URL": URL}
    for kind in KINDS:
        stem = "fixture." + kind
        png = args.current / (stem + ".png")
        chunks, text, pixels, idat = parse_png(png)
        require(text == expected, f"{png}: software metadata mismatch: {text}")
        require(not any(key in text for key in ("Author", "Copyright")), f"{png}: artwork attribution")
        stripped = check_svg(args.current / (stem + ".svg"), software, args.version)
        if args.baseline:
            before, before_text, before_pixels, before_idat = parse_png(args.baseline / (stem + ".png"))
            require(not before_text, f"{stem}: baseline already contains attribution")
            require(pixels == before_pixels, f"{stem}: decoded pixels changed")
            require(idat == before_idat, f"{stem}: compressed pixel data changed")
            require([raw for tag, _, raw in chunks if tag != b"iTXt"] ==
                    [raw for _, _, raw in before], f"{stem}: existing PNG chunks changed")
            baseline_svg = args.baseline / (stem + ".svg")
            ET.parse(baseline_svg)
            require(stripped == baseline_svg.read_bytes(), f"{stem}: visible SVG bytes changed")
        print(f"PASS {kind}: valid PNG CRC/UTF-8/RGBA and SVG XML/software attribution" +
              ("; baseline pixels, IDAT, PNG chunks and visible SVG bytes identical" if args.baseline else ""))


if __name__ == "__main__":
    main()
