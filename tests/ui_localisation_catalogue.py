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

"""Check the four-language UI catalogues and their string-key callsites.

Run from any directory: python tests/ui_localisation_catalogue.py
This complements the compiled UI/glyph tests; it does not test rendering.
"""

import ast
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
LITERAL = r'"(?:\\.|[^"\\])*"'
ROW = re.compile(r'\{\s*(' + LITERAL + r')\s*,\s*(' + LITERAL
                 + r')\s*,\s*(' + LITERAL + r')\s*,\s*(' + LITERAL + r')\s*\}')
FORMAT = re.compile(r'%(?:\d+\$)?[-+ #0]*(?:\d+|\*)?(?:\.(?:\d+|\*))?'
                    r'(?:hh|ll|[hljztL])?[diuoxXfFeEgGaAcspn%]')


def source(name):
    return (ROOT / name).read_text(encoding="utf-8")


def rows(name):
    return [tuple(ast.literal_eval(value) for value in match.groups())
            for match in ROW.finditer(source(name))]


def main():
    paths = ["src/ui/ui.cpp", "src/ui/orbital_ui_text.cpp",
             "src/ui/nbo_aomo_text.cpp", "src/ui/nbo_aomo_hover.cpp"]
    checked = 0
    for path in paths:
        catalogue = rows(path)
        assert catalogue, f"No four-language entries found in {path}"
        for row in catalogue:
            assert all(row), f"Empty translation in {path}: {row[0]}"
            expected = FORMAT.findall(row[0])
            for value in row[1:]:
                assert FORMAT.findall(value) == expected, (
                    f"Inconsistent printf arguments in {path}: {row}")
            checked += 1

    # Most NBO/scene controls use inline localisers instead of a catalogue.
    inline = re.compile(r'(?:nbo_local|scene_text)\(\s*language\s*,\s*('
                        + LITERAL + r')\s*,\s*(' + LITERAL + r')\s*,\s*('
                        + LITERAL + r')\s*,\s*(' + LITERAL + r')\s*\)')
    inline_count = 0
    for path in [*ROOT.glob("src/ui/*.cpp"), ROOT / "src/main.cpp"]:
        for match in inline.finditer(path.read_text(encoding="utf-8")):
            row = tuple(ast.literal_eval(value) for value in match.groups())
            assert all(row), f"Empty inline translation in {path}: {row}"
            assert all(FORMAT.findall(value) == FORMAT.findall(row[0])
                       for value in row[1:]), f"Inline printf mismatch in {path}: {row}"
            inline_count += 1

    nbo_ui = source("src/ui/nbo_ui.cpp")
    label_block = nbo_ui.split("constexpr const char* labels[4][Count] = {", 1)[1].split("};", 1)[0]
    label_rows = [[ast.literal_eval(value) for value in re.findall(LITERAL, row)]
                  for row in re.findall(r'\{([^{}]*)\}', label_block)]
    assert len(label_rows) == 4 and all(len(row) == 29 for row in label_rows), (
        "NBO labels do not have all four languages and 29 entries")
    report_block = nbo_ui.split("ReportText report_text(Language language)", 1)[1].split(
        "std::string fmt", 1)[0]
    report_rows = [[ast.literal_eval(value) for value in re.findall(LITERAL, row)]
                   for row in re.findall(r'return\s*\{([^{}]*)\}', report_block)]
    assert len(report_rows) == 4 and all(len(row) == 27 for row in report_rows), (
        "NBO report text does not have all four languages and 27 entries")
    assert "std::string label=c.label" not in nbo_ui, "Backend candidate label reached UI"

    aomo_rows = rows("src/ui/nbo_aomo_text.cpp")
    keys = {row[0] for row in aomo_rows}
    assert len(keys) == len(aomo_rows), "Duplicate AOMO translation key"
    calls = re.compile(r'aomo_text\(\s*(?:language|state.language|view.language)\s*,\s*('
                       + LITERAL + ')')
    used = set()
    for path in [*ROOT.glob("src/ui/*.cpp"), ROOT / "src/main.cpp"]:
        for match in calls.finditer(path.read_text(encoding="utf-8")):
            key = ast.literal_eval(match.group(1))
            assert key in keys, f"English fallback at {path}: {key}"
            used.add(key)

    # These labels are selected indirectly rather than by a literal call.
    dynamic_keys = {
        "Valence d-pi channel", "Additional p-pi channel", "Additional radial channel",
        "Transverse p reference space", "Pi donation", "Pi back-donation",
        "Pi back-donation and donation", "Occupied-space mixing", "Pi mixing",
        "Pi donor field", "Pi acceptor field", "Mixed pi response",
        "Pi field below display resolution", "Current metal s / p / d",
        "Current metal f", "Other metal shells", "Ligand valence", "Other ligand space",
        "Core composition", "Uncovered composition", "Metal valence dominated",
        "Ligand valence dominated", "Additional space dominated", "Mixed composition",
        "The folder could not be read.", "The path could not be opened.",
    }
    assert dynamic_keys <= keys, f"Missing dynamic translations: {dynamic_keys - keys}"

    main_ui = source("src/main.cpp")
    assert 'ImGui::Text("Atom %zu' not in main_ui, "Untranslated atom tooltip"
    assert "status_detail = e.what()" not in main_ui, "Raw exception in the status line"
    assert "status_detail=e.what()" not in main_ui, "Raw exception in the status line"
    assert "status_detail=nbo_ui.error;" not in main_ui, "Raw NBO exception in status"
    for path in ["src/ui/nbo_aomo_hover.cpp", "src/ui/nbo_aomo_ui.cpp"]:
        text = source(path)
        assert not re.search(r'(?:lines\.push_back|ImGui::TextWrapped)\([^;\n]*'
                             r'(?:info\.detail|info\.energy_status|info\.occupation_status|'
                             r'node\.spatial_spin->energy_status)', text), (
            f"Untranslated spin metadata in normal UI: {path}")

    basic = {row[0]: row for row in rows("src/ui/ui.cpp")}
    assert basic["Molecule opacity"][1] == "分子不透明度"
    assert basic["Orbital opacity"][1] == "轨道不透明度"
    # Each new string uses the AOMO catalogue, all of which is included in the
    # font atlas. The compiled smoke test checks actual translated seed output.
    assert "for(const auto& row:rows){seed+=row[column(language)]" in source(
        "src/ui/nbo_aomo_text.cpp")
    assert "nbo_aomo_text_glyph_seed(language)" in source("src/ui/nbo_aomo_ui.cpp")
    print(f"UI localisation: {checked} four-language rows, {inline_count} inline localisers, "
          f"29 NBO labels, 27 report fields, {len(used)} literal keys, "
          f"{len(dynamic_keys)} indirect keys and regression checks passed")


if __name__ == "__main__":
    main()
