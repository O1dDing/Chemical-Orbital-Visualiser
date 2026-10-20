"""Prepare and inspect native, one-drop AOMO category acceptance runs.

This driver never launches COV. Every run directory is caller-owned and fresh;
native capture, rendering, and GPU scheduling remain with the orchestrator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CASES = ROOT / "work/aomo-cases/cases.json"
PACKAGES = ROOT / "work/aomo-validation/category-packages"
PLANS = ROOT / "work/aomo-validation/category-plans"
W_KEYS = ("aonbo", "nbomo", "sao", "naomo", "aonao", "naonbo", "aonho",
          "aonlmo", "aopnao", "naonho", "naonlmo", "nhonbo", "nbonlmo",
          "nlmomo", "aomo")
SOURCE_KEYS = ("fchk", "report", "archive47") + W_KEYS
DEST = {"fchk": "canonical.fchk", "report": "analysis.log", "archive47": "FILE.47"}
for _key, _number in zip(W_KEYS, (37, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                                   58, 59, 60, 61, 62), strict=True):
    DEST[_key] = f"FILE.{_number}"


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(value)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_cases() -> list[dict]:
    data = json.loads(CASES.read_text(encoding="utf-8"))
    if len(data) != 61 or len({x["case_id"] for x in data}) != 61:
        raise ValueError("Expected 61 distinct representative cases")
    return data


def select_cases(wanted: list[str], manifest: Path | None = None) -> list[dict]:
    cases = json.loads(manifest.read_text(encoding="utf-8")) if manifest else source_cases()
    if manifest and len({c["case_id"] for c in cases}) != len(cases):
        raise ValueError("Duplicate case IDs in explicit manifest")
    found = [c for c in cases if not wanted or c["case_id"] in wanted]
    if wanted and {c["case_id"] for c in found} != set(wanted):
        raise ValueError("Unknown case ID(s): " + ", ".join(sorted(set(wanted) -
                            {c["case_id"] for c in found})))
    return found


def next_dir(base: Path, prefix: str) -> Path:
    base.mkdir(parents=True, exist_ok=True)
    for number in range(1, 1000):
        path = base / f"{prefix}-{number:03d}"
        try:
            path.mkdir()
            return path
        except FileExistsError:
            pass
    raise RuntimeError(f"Exhausted attempt numbers in {base}")


def prepare(cases: list[dict]) -> list[dict]:
    records = []
    for case in cases:
        cid = case["case_id"]
        folder = next_dir(PACKAGES / cid, "package")
        try:
            files = []
            for key in SOURCE_KEYS:
                source = Path(case[key])
                if not source.is_file():
                    raise FileNotFoundError(f"{cid} {key}: {source}")
                target = folder / DEST[key]
                shutil.copy2(source, target)
                source_sha = digest(source)
                copied_sha = digest(target)
                if source_sha != copied_sha:
                    raise ValueError(f"Copy hash mismatch: {source}")
                if key == "fchk" and case.get("fchk_sha256") and source_sha != case["fchk_sha256"]:
                    raise ValueError(f"Producer FCHK hash mismatch: {cid}")
                files.append(dict(key=key, member=target.name, source=str(source),
                                  bytes=target.stat().st_size, sha256=copied_sha))
            # The importer resolves only contained relative members. No paths
            # or source provenance from outside this directory enter the drop.
            manifest = folder / "drop.covnbopkg"
            write_text(manifest, "COV_NBO_PACKAGE 1\n" +
                       "\n".join(x["member"] for x in files) + "\n")
            package = dict(schema="cov.aomo.product-package.v1", case_id=cid,
                           producer_manifest=case["selected_producer_manifest"],
                           fchk_relationship=case["fchk_relationship"],
                           charge=case["charge"], multiplicity=case["multiplicity"],
                           method=case["method"], basis=case["basis"],
                           covers=case["covers"], geometry_label=case.get("geometry_label"),
                           drag_manifest=str(manifest), files=files,
                           state="prepared_unexecuted")
            write_json(folder / "package.json", package)
            records.append(dict(case_id=cid, package=str(folder), bytes=sum(x["bytes"] for x in files)))
        except Exception as exc:
            write_json(folder / "prepare-error.json", {"case_id": cid, "error": repr(exc)})
            raise
    return records


def latest_package(cid: str) -> Path:
    dirs = sorted((PACKAGES / cid).glob("package-*"))
    good = [p for p in dirs if (p / "package.json").is_file()]
    if not good:
        raise FileNotFoundError(f"No complete package for {cid}")
    return good[-1]


def native_command(package: Path, plan: Path, output: Path) -> list[str]:
    return ["<cov_validation.exe>", str(package / "canonical.fchk"),
            "--validation-plan", str(plan), "--validation-output", str(output),
            "--validation-background"]


def quoted(value: str | Path) -> str:
    s = str(value).replace("\\", "/")
    if '"' in s or "\n" in s:
        raise ValueError("Unsafe native plan string")
    return f'"{s}"'


def capture_plan(case: dict) -> dict:
    cid = case["case_id"]
    package = latest_package(cid)
    folder = next_dir(PLANS / cid, "capture")
    plan = folder / "discovery.plan"
    output = folder / "native-output"
    lines = ["COV_VALIDATION 1", "window 2100 1250",
             "scene 0.32 0.65 0.35 2.2 0.03 96",
             "drop " + quoted(package / "drop.covnbopkg"),
             'wait "load"', 'seek "panel.diagram"']
    if case.get("format_variant") in ("pure", "cartesian"):
        lines += ['seek "aomo.show_rydberg"', 'click "aomo.show_rydberg"',
                  'seek "aomo.graph.canvas"']
    lines += ['capture "overall"', 'hover "aomo.graph.canvas"',
             'capture "graph"', 'hover "scene.viewport"',
             'capture "right-molecule"']
    write_text(plan, "\n".join(lines) + "\n")
    task = dict(schema="cov.aomo.product-run.v1", phase="capture", case_id=cid,
                package=str(package), plan=str(plan), output=str(output),
                native_command=native_command(package, plan, output),
                state="prepared_unexecuted", coverage=case["covers"],
                geometry_label=case.get("geometry_label"),
                multiplicity=case["multiplicity"])
    write_json(folder / "run.json", task)
    return task


def g_pan_capture_plan(case: dict) -> dict:
    """Expose actual high-angular-momentum graph rows by user-visible controls."""
    if case.get("format_variant") not in ("pure", "cartesian"):
        raise ValueError("g_pan_requires_a_g_shell_format_variant")
    cid = case["case_id"]
    package = latest_package(cid)
    folder = next_dir(PLANS / cid, "g-pan")
    plan = folder / "g-discovery.plan"
    output = folder / "native-output"
    lines = ["COV_VALIDATION 1", "window 2100 1250",
             "scene 0.32 0.65 0.35 2.2 0.03 96",
             "drop " + quoted(package / "drop.covnbopkg"),
             'wait "load"', 'seek "panel.diagram"',
             'capture "overall"',
             'seek "aomo.show_rydberg"', 'click "aomo.show_rydberg"',
             'seek "aomo.graph.canvas"']
    drag_count = 6 if case["format_variant"] == "cartesian" else 4
    lines += ['click "aomo.zoom.out"'] * 6
    lines += ['drag "aomo.graph.canvas" 0 -250'] * drag_count
    lines += ['hover "aomo.graph.canvas"', 'capture "graph"',
              'hover "scene.viewport"', 'capture "right-molecule"']
    write_text(plan, "\n".join(lines) + "\n")
    task = dict(schema="cov.aomo.product-run.v1", phase="capture", case_id=cid,
                package=str(package), plan=str(plan), output=str(output),
                native_command=native_command(package, plan, output),
                state="prepared_unexecuted", coverage=case["covers"],
                geometry_label=None, multiplicity=2,
                g_discovery=dict(show_rydberg=True, zoom_out_clicks=6,
                                 canvas_drag_dy=[-250] * drag_count))
    write_json(folder / "run.json", task)
    return task


def target_ids(ui: dict, prefix: str) -> list[str]:
    return [x["id"] for x in ui.get("targets", []) if x.get("visible") and
            x.get("id", "").startswith(prefix)]


def choose_interaction(ui: dict, integration: dict,
                       scene_override: str | None = None) -> dict:
    visible = {x["id"]: x for x in ui.get("targets", []) if x.get("visible")}
    targets = set(visible)
    if "aomo.graph.canvas" not in targets or "scene.viewport" not in targets:
        raise ValueError("Discovery lacks visible graph or molecule viewport")
    orbitals = integration.get("orbitals", [])
    canonical = {o["ref"]["index"]: o for o in orbitals
                 if o.get("ref", {}).get("kind") == "canonical"}
    mo_ids = [x for x in targets if re.fullmatch(r"aomo\.node\.canonical_mo:\d+", x)]
    available_mos = []
    for item in mo_ids:
        index = int(item.rsplit(":", 1)[1])
        orbital = canonical.get(index)
        if orbital:
            available_mos.append((item, orbital))
    if not available_mos:
        raise ValueError("No visible canonical MO node in discovery")
    valence_nodes = []
    for row in integration.get("dataset", {}).get("naos", []):
        if not str(row.get("type", "")).strip().lower().startswith("val"):
            continue
        node = f'aomo.node.NAO:{row["spin"]}:{row["id"] - 1}'
        if node in targets:
            valence_nodes.append((node, row))
    pairs = [(mid, mo, nid, nao) for mid, mo in available_mos
             for nid, nao in valence_nodes if nao["spin"] == mo["ref"]["spin"]]
    if not pairs:
        raise ValueError("No visible same-spin MO and valence NAO pair; stage another capture")
    # The actual central view may feature virtual orbitals (notably open-shell
    # complexes). Pick a visible same-spin pair near the rendered MO; this is
    # an interface action, never an occupied-orbital chemistry assertion.
    displayed_index = ui.get("state", {}).get("rendered_mo")
    if not isinstance(displayed_index, int):
        displayed_index = 0
    mo_id, mo, nao_id, nao = min(pairs, key=lambda row: (
        abs(row[1]["ref"]["index"] - displayed_index),
        -(row[3].get("occupation") or 0), row[1]["ref"]["index"]))
    beta_mos = [(mid, orbital) for mid, orbital in available_mos
                if orbital["ref"]["spin"] == "beta" and mid != mo_id]
    beta_id = min(beta_mos, key=lambda row: row[1]["ref"]["index"])[0] if beta_mos else None
    structure = integration.get("structure", [])
    bond = []
    for item in targets:
        match = re.fullmatch(r"scene\.bond\.(\d+)", item)
        if match and int(match.group(1)) < len(structure):
            evidence = structure[int(match.group(1))]
            if evidence.get("atoms"):
                bond.append((item, evidence))
    atoms = [x for x in targets if re.fullmatch(r"scene\.atom\.\d+", x)]
    def centre(item: dict) -> tuple[float, float]:
        r = item["rect"]
        return ((r[0] + r[2]) / 2, (r[1] + r[3]) / 2)
    margins = [(min(math.dist(centre(visible[item]), centre(visible[atom]))
                    for atom in atoms), item) for item, _ in bond] if atoms else []
    # Pick a bond only when its hit centre is separated from all atom hit
    # centres. Close atom/bond overlays can legitimately resolve to an atom.
    if scene_override:
        if scene_override not in targets or not re.fullmatch(r"scene\.(?:atom|bond)\.\d+", scene_override):
            raise ValueError("Override must be an actually visible scene atom/bond target")
        scene_id = scene_override
        scene_kind = "bond" if ".bond." in scene_id else "atom"
    elif margins and max(margins)[0] >= 24.0:
        scene_id = max(margins)[1]
        scene_kind = "bond"
    elif atoms:
        counts = Counter(a for e in structure for a in e.get("atoms", []))
        scene_id = max(atoms, key=lambda x: (
            min(math.dist(centre(visible[x]), centre(visible[y]))
                for y in atoms if y != x) if len(atoms) > 1 else 1e9,
            counts[int(x.rsplit(".", 1)[1])], x))
        scene_kind = "atom"
    elif bond:
        scene_id = sorted(bond, key=lambda row: row[0])[0][0]
        scene_kind = "bond"
    else:
        raise ValueError("No visible selectable atom or bond in molecule viewport")
    groups = sorted(x for x in targets if x.startswith("aomo.node.degenerate:"))
    if not groups and len(integration.get("orbitals", [])) > 400:
        groups = sorted(x for x in targets if x.startswith("aomo.node.atom:"))
    return dict(mo_node=mo_id, mo_ref=mo["ref"], nao_node=nao_id,
                nao_report_id=nao["id"], nao_spin=nao["spin"],
                nao_type=nao["type"], scene_node=scene_id,
                scene_kind=scene_kind, collapse_node=groups[0] if groups else None,
                beta_mo_node=beta_id)


def interaction_plan(case: dict, discovery: Path,
                     scene_override: str | None = None) -> dict:
    cid = case["case_id"]
    package = latest_package(cid)
    graph_file = discovery / "graph.ui.json"
    integration_files = sorted(discovery.glob("integration-*.json"),
                               key=lambda p: int(p.stem.split("-")[1]))
    if not graph_file.is_file() or not integration_files:
        raise FileNotFoundError("Discovery requires graph.ui.json and integration-*.json")
    ui = json.loads(graph_file.read_text(encoding="utf-8"))
    integration = json.loads(integration_files[-1].read_text(encoding="utf-8"))
    selection = choose_interaction(ui, integration, scene_override)
    if case["multiplicity"] > 1 and selection["mo_ref"]["spin"] != "beta" and \
            not selection["beta_mo_node"]:
        raise ValueError("Open-shell discovery has no visible beta MO; beta product coverage pending")
    folder = next_dir(PLANS / cid, "interaction")
    plan, output = folder / "interaction.plan", folder / "native-output"
    lines = ["COV_VALIDATION 1", "window 2100 1250",
             "scene 0.32 0.65 0.35 2.2 0.03 96",
             "drop " + quoted(package / "drop.covnbopkg"),
             'wait "load"', 'seek "panel.diagram"',
             'capture "overall"']
    if selection["collapse_node"]:
        lines += ["click " + quoted(selection["collapse_node"]),
                  'capture "group-collapsed"',
                  "click " + quoted(selection["collapse_node"]),
                  'capture "group-expanded"']
    # The NAO is selected while the discovered graph is still the active
    # central-MO view; changing MO can legitimately replace its NAO nodes.
    lines += ["click " + quoted(selection["nao_node"]),
              'capture "valence-nao"',
              'volume "valence-nao-field" "0"',
              "click " + quoted(selection["mo_node"]),
              'capture "central-mo"',
              'volume "central-mo-field" "0"',
              'hover "scene.viewport"',
              "click " + quoted(selection["scene_node"]),
              'capture "right-picked"']
    if selection["beta_mo_node"]:
        lines += ["click " + quoted(selection["beta_mo_node"]),
                  'capture "beta-mo"',
                  'volume "beta-mo-field" "0"']
    lines += ['export-name "category-export"',
              'click "aomo.export"','seek "aomo.export_options"','click "aomo.export_options"','seek "aomo.export_data"','click "aomo.export_data"','click "aomo.export_options"',
              'capture "exported"']
    write_text(plan, "\n".join(lines) + "\n")
    write_json(folder / "selected-ui.json", selection)
    expected = {"valence-nao-field": {"kind": "NAO", "index": selection["nao_report_id"] - 1},
                "central-mo-field": {"kind": "canonical", "index": selection["mo_ref"]["index"]}}
    if selection["beta_mo_node"]:
        expected["beta-mo-field"] = {"kind": "canonical",
                                     "index": int(selection["beta_mo_node"].rsplit(":", 1)[1])}
    write_json(folder / "expected-fields.json", expected)
    task = dict(schema="cov.aomo.product-run.v1", phase="interaction", case_id=cid,
                package=str(package), discovery=str(discovery),
                discovery_sha256={"graph": digest(graph_file),
                                  "integration": digest(integration_files[-1])},
                plan=str(plan), output=str(output),
                expected_fields=str(folder / "expected-fields.json"),
                native_command=native_command(package, plan, output),
                selected_ui=selection, multiplicity=case["multiplicity"],
                state="prepared_unexecuted")
    write_json(folder / "run.json", task)
    return task


def g_shell_interaction_plan(case: dict, discovery: Path) -> dict:
    """Choose a visible, report-labelled g NAO in each spin block."""
    cid = case["case_id"]
    if case.get("format_variant") not in ("pure", "cartesian"):
        raise ValueError("g_shell_interaction_requires_format_variant")
    package = latest_package(cid)
    graph_file = discovery / "graph.ui.json"
    integrations = sorted(discovery.glob("integration-*.json"),
                          key=lambda p: int(p.stem.split("-")[1]))
    if not graph_file.is_file() or not integrations:
        raise FileNotFoundError("g_shell_discovery_missing")
    ui = json.loads(graph_file.read_text(encoding="utf-8"))
    integration = json.loads(integrations[-1].read_text(encoding="utf-8"))
    visible = {t["id"] for t in ui.get("targets", []) if t.get("visible")}
    g_indices = sorted({int(o["ref"]["index"]) for o in integration.get("orbitals", [])
                        if o.get("ref", {}).get("kind") == "NAO" and
                        "g(" in o.get("label", "")})
    choices = [index for index in g_indices if all(
        f"aomo.node.NAO:{spin}:{index}" in visible for spin in ("alpha", "beta"))]
    if not choices:
        raise ValueError("no_same_index_visible_g_NAO_alpha_beta_pair")
    index = choices[0]
    labels = {o["ref"]["spin"]: o["label"] for o in integration["orbitals"]
              if o.get("ref", {}).get("kind") == "NAO" and
              o["ref"]["index"] == index and o["ref"]["spin"] in ("alpha", "beta")}
    if set(labels) != {"alpha", "beta"} or not all("g(" in v for v in labels.values()):
        raise ValueError("g_NAO_label_and_spin_not_proved")
    atoms = sorted(x for x in visible if re.fullmatch(r"scene\.atom\.\d+", x))
    if not atoms:
        raise ValueError("g_shell_molecule_atom_target_missing")
    picked_atom = atoms[0]
    folder = next_dir(PLANS / cid, "g-interaction")
    plan = folder / "g-interaction.plan"
    output = folder / "native-output"
    lines = ["COV_VALIDATION 1", "window 2100 1250",
             "scene 0.32 0.65 0.35 2.2 0.03 96",
             "drop " + quoted(package / "drop.covnbopkg"),
             'wait "load"', 'seek "panel.diagram"',
             'capture "overall"',
             'seek "aomo.show_rydberg"', 'click "aomo.show_rydberg"',
             'seek "aomo.graph.canvas"']
    drag_count = 6 if case["format_variant"] == "cartesian" else 4
    lines += ['click "aomo.zoom.out"'] * 6
    lines += ['drag "aomo.graph.canvas" 0 -250'] * drag_count
    lines += ['hover "aomo.graph.canvas"', 'capture "g-graph"',
              'click ' + quoted(f"aomo.node.NAO:alpha:{index}"),
              'volume "g-alpha-field" "0"',
              'click "scene.fit_component"', 'capture "g-alpha"',
              'click ' + quoted(f"aomo.node.NAO:beta:{index}"),
              'volume "g-beta-field" "0"',
              'click "scene.fit_component"', 'capture "g-beta"',
              'hover "scene.viewport"', 'click ' + quoted(picked_atom),
              'capture "right-picked"',
              'export-name "g-shell-export"', 'click "aomo.export"','seek "aomo.export_options"','click "aomo.export_options"','seek "aomo.export_data"','click "aomo.export_data"','click "aomo.export_options"',
              'capture "exported"']
    write_text(plan, "\n".join(lines) + "\n")
    expected = {"g-alpha-field": {"kind": "NAO", "index": index, "spin": "alpha"},
                "g-beta-field": {"kind": "NAO", "index": index, "spin": "beta"}}
    write_json(folder / "expected-fields.json", expected)
    task = dict(schema="cov.aomo.g-shell-product-run.v1", phase="g-interaction",
                case_id=cid, package=str(package), discovery=str(discovery),
                discovery_sha256={"graph": digest(graph_file),
                                  "integration": digest(integrations[-1])},
                plan=str(plan), output=str(output),
                expected_fields=str(folder / "expected-fields.json"),
                native_command=native_command(package, plan, output),
                g_nao_index=index, g_nao_labels=labels,
                scene_node=picked_atom, state="prepared_unexecuted")
    write_json(folder / "run.json", task)
    return task


def check_g_shell(run_file: Path) -> dict:
    task = json.loads(run_file.read_text(encoding="utf-8"))
    output = Path(task["output"])
    session = json.loads((output / "session.json").read_text(encoding="utf-8"))
    events = [json.loads(line) for line in (output / "events.jsonl").read_text(encoding="utf-8").splitlines()]
    selections = sorted((int(path.stem.split("-")[1]), json.loads(path.read_text(encoding="utf-8")))
                        for path in output.glob("selection-*.json"))
    expected = json.loads(Path(task["expected_fields"]).read_text(encoding="utf-8"))
    observed_fields = {}
    for name, wanted in expected.items():
        path = output / f"{name}.volume.json"
        if not path.is_file():
            observed_fields[name] = dict(passed=False, reason="volume_missing")
            continue
        volume = json.loads(path.read_text(encoding="utf-8"))
        prior = [selection for frame, selection in selections if frame < volume["frame"]]
        ref = prior[-1]["terms"][0]["orbital"] if prior and prior[-1].get("terms") else {}
        observed_fields[name] = dict(passed=all(ref.get(k) == wanted[k]
                                               for k in ("kind", "index", "spin")),
                                     observed=ref, expected=wanted)
    captures = {p.stem.removesuffix(".ui") for p in output.glob("*.ui.json")}
    required = {"overall", "g-graph", "g-alpha", "g-beta", "right-picked", "exported"}
    ui = json.loads((output / "g-graph.ui.json").read_text(encoding="utf-8"))
    g_nodes = {f"aomo.node.NAO:{spin}:{task['g_nao_index']}" for spin in ("alpha", "beta")}
    visible = {t["id"] for t in ui.get("targets", []) if t.get("visible")}
    actual_pick = [e.get("data", {}) for e in events if e.get("kind") == "scene.pick"]
    expected_pick = dict(kind=0, index=int(task["scene_node"].rsplit(".", 1)[1]))
    exports = {suffix: output / ("g-shell-export" + suffix) for suffix in
               (".aomo.json", ".aomo.csv", ".aomo.svg", ".aomo.png")}
    passed = (session.get("failed_commands") == 0 and required <= captures and
              g_nodes <= visible and all(x["passed"] for x in observed_fields.values()) and
              bool(actual_pick) and actual_pick[-1] == expected_pick and
              all(p.is_file() and p.stat().st_size > 0 for p in exports.values()) and
              all((output / f"{name}.bmp").is_file() for name in required))
    result = dict(schema="cov.aomo.g-shell-product-check.v1", case_id=task["case_id"],
                  status="product_observation_pass" if passed else "product_observation_fail",
                  output=str(output), executable_sha256=json.loads(
                      (run_file.parent / "process.json").read_text(encoding="utf-8"))[
                          "executable_sha256"],
                  required_captures=sorted(required), captures=sorted(captures),
                  visible_g_nodes=sorted(g_nodes & visible),
                  g_nao_labels=task["g_nao_labels"], fields=observed_fields,
                  expected_pick=expected_pick,
                  actual_pick=actual_pick[-1] if actual_pick else None,
                  export_files={k: dict(exists=p.is_file(), bytes=p.stat().st_size if p.is_file() else 0)
                                for k, p in exports.items()},
                  native_session=session)
    report = next_dir(PLANS / task["case_id"], "g-check") / "result.json"
    write_json(report, result)
    result["report"] = str(report)
    return result


def check_run(run_file: Path) -> dict:
    task = json.loads(run_file.read_text(encoding="utf-8"))
    output = Path(task["output"])
    if not (output / "session.json").is_file():
        result = dict(case_id=task["case_id"], phase=task["phase"],
                      status="not_executed_or_incomplete", output=str(output))
    else:
        session = json.loads((output / "session.json").read_text(encoding="utf-8"))
        events = [json.loads(line) for line in (output / "events.jsonl").read_text(encoding="utf-8").splitlines()]
        kinds = Counter(e.get("kind") for e in events)
        captures = {p.stem.removesuffix(".ui") for p in output.glob("*.ui.json")}
        integrations = sorted(output.glob("integration-*.json"))
        last = json.loads(integrations[-1].read_text(encoding="utf-8")) if integrations else None
        caps = {c["key"]: c["state"] for c in last.get("capabilities", [])} if last else {}
        selections = sorted(output.glob("selection-*.json"))
        selected = [json.loads(p.read_text(encoding="utf-8")) for p in selections]
        errors = {k: kinds[k] for k in ("input.package.error", "nbo.attach.error",
                   "aomo.selection.error") if kinds[k]}
        required = {"overall", "graph", "right-molecule"} if task["phase"] == "capture" else \
                   {"overall", "central-mo", "valence-nao", "right-picked", "exported"}
        export_files = {suffix: output / ("category-export" + suffix) for suffix in
                        (".aomo.json", ".aomo.csv", ".aomo.svg", ".aomo.png")}
        export_complete = all(p.is_file() and p.stat().st_size > 0
                              for p in export_files.values())
        picks = [e.get("data", {}) for e in events if e.get("kind") == "scene.pick"]
        expected_pick = None
        actual_pick = picks[-1] if picks else None
        pick_matches = False
        if task["phase"] == "interaction":
            scene_node = task["selected_ui"]["scene_node"]
            expected_pick = dict(kind=1 if ".bond." in scene_node else 0,
                                 index=int(scene_node.rsplit(".", 1)[1]))
            pick_matches = actual_pick == expected_pick
        observed = dict(captures=sorted(captures), required_captures=sorted(required),
                        integration_count=len(integrations),
                        source_association=caps.get("source_association"),
                        aomo=caps.get("aomo"), nao=caps.get("nao"),
                        selection_count=len(selections),
                        selection_kinds=[s.get("terms", [{}])[0].get("orbital", {}).get("kind")
                                         for s in selected if s.get("terms")],
                        scene_pick_count=kinds["scene.pick"],
                        expected_pick=expected_pick, actual_pick=actual_pick,
                        pick_matches=pick_matches,
                        nbo_export_count=kinds["nbo.export"],
                        export_actual_count=kinds["export.actual"], errors=errors)
        observed["aomo_export_files"] = {key: dict(exists=p.is_file(),
                                                bytes=p.stat().st_size if p.is_file() else 0,
                                                sha256=digest(p) if p.is_file() else None)
                                         for key, p in export_files.items()}
        basic = (session.get("failed_commands") == 0 and not errors and required <= captures
                 and caps.get("source_association") == "available"
                 and caps.get("aomo") == "available")
        if task["phase"] == "interaction":
            basic = basic and ("canonical" in observed["selection_kinds"] and
                               "NAO" in observed["selection_kinds"] and
                               pick_matches and export_complete)
            if task.get("selected_ui", {}).get("beta_mo_node"):
                beta_state = output / "beta-mo.ui.json"
                basic = basic and beta_state.is_file() and \
                    json.loads(beta_state.read_text(encoding="utf-8"))["state"].get("rendered_spin") == "beta"
            elif task.get("multiplicity", 1) > 1:
                central_state = output / "central-mo.ui.json"
                basic = basic and central_state.is_file() and \
                    json.loads(central_state.read_text(encoding="utf-8"))["state"].get("rendered_spin") == "beta"
        result = dict(case_id=task["case_id"], phase=task["phase"],
                      status="product_observation_pass" if basic else "product_observation_fail",
                      output=str(output), native_session=session, observed=observed,
                      checker_version=2,
                      scientific_status="requires_independent_case_reference")
    report = next_dir(PLANS / task["case_id"], "check") / "result.json"
    write_json(report, result)
    result["report"] = str(report)
    return result


def execute_run(run_file: Path, exe: Path, timeout: int) -> dict:
    """One bounded native process; caller serializes invocations across cases."""
    from validation_process import physical_core_masks, run_tree

    task = json.loads(run_file.read_text(encoding="utf-8"))
    output = Path(task["output"])
    if not exe.is_file():
        raise FileNotFoundError(exe)
    if output.exists() or (run_file.parent / "process.json").exists():
        raise FileExistsError(f"Native attempt already has output: {output}")
    executable_sha256 = digest(exe)
    expected = [str(exe.resolve())] + task["native_command"][1:]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="2", OMP_THREAD_LIMIT="2", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="2")
    result = run_tree(expected, run_file.parent, env,
                      sum(physical_core_masks(2)), 20, timeout, affinity=True,
                      on_started=lambda data: write_json(
                          run_file.parent / "process-start.json",
                          {**data, "executable_sha256": executable_sha256}))
    result["executable_sha256"] = executable_sha256
    result["executable_path"] = str(exe.resolve())
    write_json(run_file.parent / "process.json", result)
    return dict(case_id=task["case_id"], phase=task["phase"],
                exit_code=result.get("exit_code"), timed_out=result.get("timed_out"),
                wall_seconds=result.get("wall_seconds"),
                peak_tree_commit_bytes=result.get("peak_tree_commit_bytes"),
                process=str(run_file.parent / "process.json"), output=str(output))


def current_check(run_file: Path) -> dict:
    task = json.loads(run_file.read_text(encoding="utf-8"))
    files = sorted((PLANS / task["case_id"]).glob("check-*/result.json"))
    for path in reversed(files):
        row = json.loads(path.read_text(encoding="utf-8"))
        if row.get("phase") == task["phase"] and row.get("output") == task["output"] and \
                row.get("checker_version") == 2:
            return row
    return check_run(run_file)


def campaign(cases: list[dict], limit: int, exe: Path, timeout: int) -> dict:
    # Persistent append-only progress lets the root inspect a long serial queue
    # and preserves every failure if a later case needs a runner/core repair.
    logdir = PLANS / "campaign-001"
    logdir.mkdir(parents=True, exist_ok=True)
    logfile = logdir / "progress.jsonl"
    completed = 0
    skipped = 0
    def emit(row: dict) -> None:
        row = {"utc": datetime.now(timezone.utc).isoformat(), **row}
        with logfile.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, ensure_ascii=False) + "\n")
            stream.flush()
        print(json.dumps(row, ensure_ascii=False), flush=True)
    for case in cases:
        cid = case["case_id"]
        existing = sorted((PLANS / cid).glob("interaction-*/run.json"))
        done = False
        for run in reversed(existing):
            task = json.loads(run.read_text(encoding="utf-8"))
            gpu = Path(task["output"]) / "independent-gpu.json"
            if gpu.is_file() and json.loads(gpu.read_text(encoding="utf-8")).get("passed") and \
                    current_check(run).get("status") == "product_observation_pass":
                done = True
                break
        if done:
            skipped += 1
            continue
        if limit and completed >= limit:
            break
        try:
            capture_runs = sorted((PLANS / cid).glob("capture-*/run.json"))
            if not capture_runs:
                capture_plan(case)
                capture_runs = sorted((PLANS / cid).glob("capture-*/run.json"))
            capture = capture_runs[-1]
            capture_task = json.loads(capture.read_text(encoding="utf-8"))
            capture_output = Path(capture_task["output"])
            if not (capture_output / "session.json").is_file():
                emit(dict(case_id=cid, phase="capture", state="starting"))
                execute_run(capture, exe, timeout)
            capture_check = current_check(capture)
            if capture_check["status"] != "product_observation_pass":
                raise RuntimeError("capture product observation failed")
            usable = [r for r in existing if
                      (Path(json.loads(r.read_text(encoding="utf-8"))["output"]) / "session.json").is_file()
                      and json.loads(r.read_text(encoding="utf-8")).get("expected_fields")]
            if usable:
                interaction = usable[-1]
            else:
                interaction_plan(case, capture_output)
                interaction = sorted((PLANS / cid).glob("interaction-*/run.json"))[-1]
            task = json.loads(interaction.read_text(encoding="utf-8"))
            output = Path(task["output"])
            if not (output / "session.json").is_file():
                emit(dict(case_id=cid, phase="interaction", state="starting"))
                execute_run(interaction, exe, timeout)
            product = current_check(interaction)
            if product["status"] != "product_observation_pass":
                raise RuntimeError("interaction product observation failed")
            gpu_file = output / "independent-gpu.json"
            if gpu_file.is_file():
                gpu = json.loads(gpu_file.read_text(encoding="utf-8"))
            else:
                sys.path.insert(0, "E:/Dev/cov-validation-20260905/reference-deps")
                os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="2")
                from aomo_gpu_check import check as gpu_check
                gpu = gpu_check(Path(task["package"]), output,
                                json.loads(Path(task["expected_fields"]).read_text(encoding="utf-8")))
            if not gpu.get("passed"):
                raise RuntimeError("independent signed GPU field reference failed")
            completed += 1
            emit(dict(case_id=cid, phase="complete", state="product_and_gpu_pass",
                      interaction=str(interaction), output=str(output),
                      fields=len(gpu.get("fields", [])),
                      visual_status="pending_individual_image_review"))
        except Exception as exc:
            emit(dict(case_id=cid, phase="failed", state="stopped",
                      error=f"{type(exc).__name__}: {exc}"))
            raise
    return dict(completed_this_call=completed, skipped_prior_pass=skipped,
                progress_log=str(logfile))


def replay(cases: list[dict], limit: int, exe: Path, timeout: int,
           expected_exe_sha256: str) -> dict:
    """Run passed interaction plans on one frozen final binary, serially.

    Keep each result in a fresh attempt. Reuse package and UI discovery from
    the original run; only optional group screenshots are omitted to avoid
    another large image set. All interactions, fields, picks, and exports run.
    """
    actual_sha = digest(exe)
    if actual_sha.lower() != expected_exe_sha256.lower():
        raise RuntimeError(f"final_executable_sha256_mismatch:{actual_sha}")
    logdir = PLANS / "final-replay-001"
    logdir.mkdir(parents=True, exist_ok=True)
    logfile = logdir / "progress.jsonl"
    completed = 0
    skipped = 0
    def emit(row: dict) -> None:
        row = {"utc": datetime.now(timezone.utc).isoformat(),
               "executable_sha256": actual_sha, **row}
        with logfile.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, ensure_ascii=False) + "\n")
            stream.flush()
        print(json.dumps(row, ensure_ascii=False), flush=True)
    for case in cases:
        cid = case["case_id"]
        previous = sorted((PLANS / cid).glob("replay-*/run.json"))
        done = False
        for run in reversed(previous):
            task = json.loads(run.read_text(encoding="utf-8"))
            output = Path(task["output"])
            process_file = run.parent / "process.json"
            gpu_file = output / "independent-gpu.json"
            if not (process_file.is_file() and gpu_file.is_file()):
                continue
            process = json.loads(process_file.read_text(encoding="utf-8"))
            gpu = json.loads(gpu_file.read_text(encoding="utf-8"))
            if process.get("executable_sha256") == actual_sha and gpu.get("passed") and \
                    current_check(run).get("status") == "product_observation_pass":
                done = True
                break
        if done:
            skipped += 1
            continue
        if limit and completed >= limit:
            break
        try:
            sources = sorted((PLANS / cid).glob("interaction-*/run.json"))
            source_run = None
            for candidate in reversed(sources):
                source_task = json.loads(candidate.read_text(encoding="utf-8"))
                gpu_file = Path(source_task["output"]) / "independent-gpu.json"
                if gpu_file.is_file() and \
                        json.loads(gpu_file.read_text(encoding="utf-8")).get("passed") and \
                        current_check(candidate).get("status") == "product_observation_pass":
                    source_run = candidate
                    break
            if source_run is None:
                raise RuntimeError("no_passed_original_interaction_plan")
            source_task = json.loads(source_run.read_text(encoding="utf-8"))
            folder = next_dir(PLANS / cid, "replay")
            plan = folder / "interaction-replay.plan"
            original_lines = Path(source_task["plan"]).read_text(encoding="utf-8").splitlines()
            omitted = ['capture "group-collapsed"', 'capture "group-expanded"']
            kept_lines = [line for line in original_lines if line not in omitted]
            write_text(plan, "\n".join(kept_lines) + "\n")
            expected = folder / "expected-fields.json"
            write_text(expected, Path(source_task["expected_fields"]).read_text(encoding="utf-8"))
            output = folder / "native-output"
            task = dict(source_task, plan=str(plan), output=str(output),
                        expected_fields=str(expected),
                        native_command=native_command(Path(source_task["package"]), plan, output),
                        replay_source_run=str(source_run),
                        replay_source_plan_sha256=digest(Path(source_task["plan"])),
                        replay_omitted_optional_captures=omitted,
                        expected_executable_sha256=actual_sha,
                        state="prepared_unexecuted")
            write_json(folder / "run.json", task)
            emit(dict(case_id=cid, phase="replay", state="starting", run=str(folder / "run.json")))
            process = execute_run(folder / "run.json", exe, timeout)
            if process["exit_code"] != 0 or process["timed_out"]:
                raise RuntimeError("final_replay_native_process_failed")
            product = check_run(folder / "run.json")
            if product["status"] != "product_observation_pass":
                raise RuntimeError("final_replay_product_observation_failed")
            sys.path.insert(0, "E:/Dev/cov-validation-20260905/reference-deps")
            os.environ.update(OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="2")
            from aomo_gpu_check import check as gpu_check
            gpu = gpu_check(Path(task["package"]), output,
                            json.loads(expected.read_text(encoding="utf-8")))
            if not gpu.get("passed"):
                raise RuntimeError("final_replay_independent_signed_GPU_failed")
            completed += 1
            emit(dict(case_id=cid, phase="replay", state="product_and_gpu_pass",
                      run=str(folder / "run.json"), output=str(output),
                      fields=len(gpu.get("fields", []))))
        except Exception as exc:
            emit(dict(case_id=cid, phase="replay", state="stopped",
                      error=f"{type(exc).__name__}: {exc}"))
            raise
    return dict(completed_this_call=completed, skipped_prior_pass=skipped,
                executable_sha256=actual_sha, progress_log=str(logfile))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="phase", required=True)
    p = sub.add_parser("size")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--manifest", type=Path)
    p = sub.add_parser("prepare")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--manifest", type=Path)
    p = sub.add_parser("capture")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--manifest", type=Path)
    p = sub.add_parser("g-pan")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--manifest", type=Path, required=True)
    p = sub.add_parser("g-interaction")
    p.add_argument("case_id")
    p.add_argument("discovery", type=Path)
    p.add_argument("--manifest", type=Path, required=True)
    p = sub.add_parser("g-check")
    p.add_argument("run_file", type=Path)
    p = sub.add_parser("interaction")
    p.add_argument("case_id")
    p.add_argument("discovery", type=Path)
    p.add_argument("--scene-target")
    p.add_argument("--manifest", type=Path)
    p = sub.add_parser("check")
    p.add_argument("run_file", type=Path)
    p = sub.add_parser("execute")
    p.add_argument("run_file", type=Path)
    p.add_argument("--exe", type=Path, default=ROOT / "work/build-ON/cov_validation.exe")
    p.add_argument("--timeout", type=int, default=600)
    p = sub.add_parser("campaign")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--manifest", type=Path)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--exe", type=Path, default=ROOT / "work/build-ON/cov_validation.exe")
    p.add_argument("--timeout", type=int, default=600)
    p = sub.add_parser("replay")
    p.add_argument("--case", action="append", default=[])
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--exe", type=Path, default=ROOT / "work/build-ON/cov_validation.exe")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--expected-exe-sha256", required=True)
    sub.add_parser("summary")
    args = parser.parse_args()
    try:
        if args.phase == "size":
            cases = select_cases(args.case, args.manifest)
            bytes_total = sum(Path(c[key]).stat().st_size for c in cases for key in SOURCE_KEYS)
            result = dict(cases=len(cases), bytes=bytes_total, gib=bytes_total / 1024**3)
        elif args.phase == "prepare":
            result = prepare(select_cases(args.case, args.manifest))
        elif args.phase == "capture":
            result = [capture_plan(c) for c in select_cases(args.case, args.manifest)]
        elif args.phase == "g-pan":
            result = [g_pan_capture_plan(c) for c in select_cases(args.case, args.manifest)]
        elif args.phase == "g-interaction":
            result = g_shell_interaction_plan(select_cases([args.case_id], args.manifest)[0],
                                               args.discovery)
        elif args.phase == "g-check":
            result = check_g_shell(args.run_file)
        elif args.phase == "interaction":
            result = interaction_plan(select_cases([args.case_id], args.manifest)[0], args.discovery,
                                      args.scene_target)
        elif args.phase == "execute":
            result = execute_run(args.run_file, args.exe, args.timeout)
        elif args.phase == "campaign":
            result = campaign(select_cases(args.case, args.manifest), args.limit, args.exe, args.timeout)
        elif args.phase == "replay":
            result = replay(select_cases(args.case), args.limit, args.exe, args.timeout,
                            args.expected_exe_sha256)
        elif args.phase == "summary":
            cases = source_cases()
            entries = []
            for case in cases:
                cid = case["case_id"]
                packages = sorted((PACKAGES / cid).glob("package-*/package.json"))
                captures = sorted((PLANS / cid).glob("capture-*/run.json"))
                interactions = sorted((PLANS / cid).glob("interaction-*/run.json"))
                checks = sorted((PLANS / cid).glob("check-*/result.json"))
                entries.append(dict(case_id=cid, covers=case["covers"],
                                    geometry_label=case.get("geometry_label"),
                                    source_group=case["source_group"],
                                    mathematical_tags=case.get("mathematical_tags", []),
                                    fchk_relationship=case["fchk_relationship"],
                                    packages=[str(x) for x in packages],
                                    capture_runs=[str(x) for x in captures],
                                    interaction_runs=[str(x) for x in interactions],
                                    checks=[str(x) for x in checks],
                                    latest_check=json.loads(checks[-1].read_text(encoding="utf-8"))
                                    if checks else None))
            result = dict(schema="cov.aomo.product-category-summary.v1",
                          total_cases=len(entries), cases=entries,
                          status="prepared_or_observed_not_global_scientific_pass")
            out = next_dir(PLANS, "summary") / "summary.json"
            write_json(out, result)
            result = dict(path=str(out), cases=len(entries),
                          capture_ready=sum(bool(x["capture_runs"]) for x in entries),
                          interaction_ready=sum(bool(x["interaction_runs"]) for x in entries),
                          checks=sum(bool(x["checks"]) for x in entries))
        else:
            result = check_run(args.run_file)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except Exception as exc:
        print(f"aomo_product_cases: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
