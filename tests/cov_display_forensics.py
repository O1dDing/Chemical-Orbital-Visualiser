"""Collect final native UI draws without exports, fabricated text, or source changes.

Inventory is the default. --sweep or --all-nodes explicitly visits the source
nodes included in that inventory; folded headers are reported separately rather
than counted as inspected orbitals. Each case uses two isolated native sessions,
one inventory and one batch sweep. An unreachable target remains a failed attempt.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
import hashlib
import gzip
import json
import math
import ntpath
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import zipfile


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def compact_bundle(native: Path, destination: Path, minimum_bytes: int = 256) -> dict:
    """Losslessly intern repeated JSON containers, then gzip their shared DAG.

    Only actual view captures, session/actions, and the native plan enter the
    primary artifact. Source staging remains intact. SHA-256 addresses canonical
    decoded JSON; equal hashes must also have equal encoded values. Tagged
    containers avoid collisions with any source dictionary's keys. Containers
    referenced only once after parent deduplication are inlined again.
    """
    files = sorted(native.glob("*.view.json"))
    files += [p for p in (native / "session.json", native / "actions.jsonl", native / "plan.txt") if p.is_file()]
    if not files:
        raise ValueError("No actual native artifacts are available for a compact bundle")

    def read(path: Path) -> object:
        if path.suffix == ".jsonl":
            return read_jsonl(path)
        if path.suffix == ".txt":
            return path.read_text(encoding="utf-8")
        return json.loads(path.read_text(encoding="utf-8"))

    counts: dict[str, int] = {}

    def count(value: object) -> None:
        if not isinstance(value, (dict, list)):
            return
        data = canonical_bytes(value)
        if len(data) >= minimum_bytes:
            digest = hashlib.sha256(data).hexdigest()
            counts[digest] = counts.get(digest, 0) + 1
        children = value.values() if isinstance(value, dict) else value
        for child in children:
            count(child)

    for path in files:
        count(read(path))
    objects: list[dict] = []
    addresses: dict[str, int] = {}

    def encode(value: object) -> object:
        if not isinstance(value, (dict, list)):
            return value
        body = [0, {key: encode(child) for key, child in value.items()}] if isinstance(value, dict) else [1, [encode(child) for child in value]]
        data = canonical_bytes(value)
        digest = hashlib.sha256(data).hexdigest() if len(data) >= minimum_bytes else None
        if digest is not None and counts.get(digest, 0) > 1:
            if digest in addresses:
                index = addresses[digest]
                if canonical_bytes(objects[index]["value"]) != canonical_bytes(body):
                    raise ValueError("SHA-256 collision between unequal JSON objects")
            else:
                index = len(objects)
                addresses[digest] = index
                objects.append({"sha256": digest, "value": body})
            return [2, index]
        return body

    roots = {path.name: encode(read(path)) for path in files}
    references: dict[int, int] = {}
    visited: set[int] = set()

    def visit(value: object) -> None:
        if not isinstance(value, list):
            return
        tag, body = value
        if tag == 2:
            references[body] = references.get(body, 0) + 1
            if body not in visited:
                visited.add(body)
                visit(objects[body]["value"])
        else:
            children = body.values() if tag == 0 else body
            for child in children:
                visit(child)

    for value in roots.values():
        visit(value)
    kept = sorted(index for index, n in references.items() if n > 1)
    remap = {old: new for new, old in enumerate(kept)}

    def inline_singletons(value: object) -> object:
        if not isinstance(value, list):
            return value
        tag, body = value
        if tag == 2:
            return [2, remap[body]] if body in remap else inline_singletons(objects[body]["value"])
        return [0, {key: inline_singletons(child) for key, child in body.items()}] if tag == 0 else [1, [inline_singletons(child) for child in body]]

    bundle = {"schema": "cov.display-forensics.bundle.v1", "method": "canonical_json_sha256_shared_containers",
              "files": {name: inline_singletons(value) for name, value in roots.items()},
              "objects": [{"sha256": objects[index]["sha256"], "value": inline_singletons(objects[index]["value"])} for index in kept]}
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Reproducible gzip stream; no timestamp or private source filename header.
    with destination.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0) as stream:
            stream.write(canonical_bytes(bundle))
    return {"file": str(destination), "source_files": len(files), "view_files": sum(p.name.endswith(".view.json") for p in files),
            "source_bytes": sum(p.stat().st_size for p in files), "compressed_bytes": destination.stat().st_size,
            "shared_objects": len(kept), "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
            "lossless_decoded_json": True}


def load_bundle(path: Path) -> dict[str, object]:
    """Materialize independent decoded files and verify every shared SHA-256."""
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        bundle = json.load(stream)
    if bundle.get("schema") != "cov.display-forensics.bundle.v1":
        raise ValueError("Unsupported display bundle schema")
    objects = bundle.get("objects", [])
    decoded: dict[int, object] = {}
    active: set[int] = set()

    def decode(value: object) -> object:
        if not isinstance(value, list):
            return value
        if len(value) != 2 or value[0] not in (0, 1, 2):
            raise ValueError("Invalid tagged bundle value")
        tag, body = value
        if tag == 0:
            return {key: decode(child) for key, child in body.items()}
        if tag == 1:
            return [decode(child) for child in body]
        if not isinstance(body, int) or body < 0 or body >= len(objects) or body in active:
            raise ValueError("Invalid or cyclic shared bundle reference")
        if body not in decoded:
            active.add(body)
            decoded[body] = decode(objects[body]["value"])
            active.remove(body)
            actual = hashlib.sha256(canonical_bytes(decoded[body])).hexdigest()
            if actual != objects[body]["sha256"]:
                raise ValueError("Shared bundle object hash mismatch")
        return decoded[body]

    result = {name: copy.deepcopy(decode(value)) for name, value in bundle["files"].items()}
    for index in range(len(objects)):
        decode([2, index])
    return result


def quoted(value: str | Path) -> str:
    """Escape for C++ std::quoted, not a shell or JSON command string."""
    text = str(value).replace("\\", "/")
    if any(ord(c) < 32 for c in text):
        raise ValueError("Native plan strings cannot contain control characters")
    return '"' + text.replace('"', '\\"') + '"'


def artifact_id(value: str) -> str:
    readable = re.sub(r"[^A-Za-z0-9_-]", "-", value).strip("-")[:48] or "node"
    return readable + "-" + hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]


@dataclass
class Plan:
    commands: list[str] = field(default_factory=list)
    attempts: list[dict] = field(default_factory=list)

    def add(self, op: str, value: str | Path) -> None:
        self.commands.append(op + " " + quoted(value))

    def text(self) -> str:
        return "COV_VALIDATION 1\n" + "\n".join(self.commands) + "\n"


def inventory_plan(package: Path, capture: bool = False, language: str = "zh",
                   resolution: int = 48, preset: str = "current") -> Plan:
    plan = Plan(["window 2100 1250", f"scene 0.60 0.65 0.35 2.2 0.03 {resolution}"])
    plan.add("drop", package)
    plan.add("wait", "load")
    if language != "current":
        plan.add("click", "language")
        plan.add("key", "Home")
        for _ in range({"en": 0, "zh": 1, "ja": 2, "fr": 3}[language]):
            plan.add("key", "Down")
        plan.add("key", "Enter")
        plan.add("wait", "language-settle")
    if preset != "current":
        plan.add("seek", "aomo.preset")
        plan.add("click", "aomo.preset")
        plan.add("click", "aomo.preset." + preset)
        plan.add("wait", "preset-settle")
    plan.add("inspect", "inventory")
    if capture:
        plan.add("capture", "inventory-frame")
    return plan


def normalized_path(value: str | Path) -> str:
    text = str(value)
    if re.match(r"^[A-Za-z]:", text) or "\\" in text:
        return ntpath.normcase(ntpath.normpath(text))
    return str(Path(text).resolve())


def input_identity(view: dict, canonical: Path) -> dict:
    records = draw_records(view, "forensic.input")
    if not records:
        return {"matched": False, "failure": "Actual input identity record is missing"}
    actual = records[-1]
    matched = normalized_path(actual.get("canonical_path", "")) == normalized_path(canonical)
    result = dict(actual, matched=matched, expected_canonical_path=canonical.as_posix())
    if not matched:
        result["failure"] = "Loaded canonical path differs from the intended package"
    elif actual.get("error"):
        result.update(matched=False, failure="Native input error: " + str(actual["error"]))
    return result


def draw_records(view: dict, kind: str) -> list[dict]:
    return [item["data"] for item in view.get("draw_trace", [])
            if item.get("kind") == kind and isinstance(item.get("data"), dict)]


def node_inventory(view: dict) -> tuple[list[dict], dict]:
    records = draw_records(view, "forensic.aomo")
    kind = "forensic.aomo"
    if not records:
        records = [r for r in draw_records(view, "forensic.diagram") if "nodes" in r]
        kind = "forensic.diagram"
    if not records:
        raise ValueError("The actual inventory has no forensic node snapshot")
    snapshot = records[-1]
    seen: set[str] = set()
    nodes = []
    for node in snapshot.get("nodes", []):
        hit_id = node.get("hit_id")
        if not isinstance(hit_id, str) or not hit_id:
            raise ValueError("An included node has no stable semantic hit ID")
        if hit_id in seen:
            raise ValueError("Duplicate semantic node hit ID: " + hit_id)
        seen.add(hit_id)
        nodes.append(dict(node))
    return nodes, {"record": kind, "snapshot_id": snapshot.get("snapshot_id"),
                   "view": snapshot.get("view"), "node_count": len(nodes),
                   "group_header_count": sum(bool(n.get("group_header")) for n in nodes),
                   "scroll_clipped_count": sum(bool(n.get("scroll_clipped")) for n in nodes)}


def choose_nodes(nodes: list[dict], selected: list[str] | None = None,
                 limit: int | None = None) -> list[dict]:
    source_nodes = [n for n in nodes if not n.get("group_header")]
    if selected:
        by_hit = {n["hit_id"]: n for n in source_nodes}
        missing = sorted(set(selected) - set(by_hit))
        if missing:
            raise ValueError("Requested source hit IDs are absent from inventory: " + ", ".join(missing))
        source_nodes = [by_hit[hit] for hit in dict.fromkeys(selected)]
    return source_nodes[:limit] if limit is not None else source_nodes


def needs_choice(node: dict, inventory_nodes: list[dict]) -> bool:
    """Predict the existing overlap chooser from actual node rectangles only."""
    rect = node.get("logical_rect")
    if node.get("canonical_index") is None or not rect or len(rect) != 4:
        return False
    x, y = (rect[0] + rect[2]) / 2, (rect[1] + rect[3]) / 2
    count = 0
    for other in inventory_nodes:
        r = other.get("logical_rect")
        if other.get("canonical_index") is not None and r and len(r) == 4 and r[0] <= x <= r[2] and r[1] <= y <= r[3]:
            count += 1
    return count > 1


def sweep_plan(package: Path, nodes: list[dict], capture: bool = False,
               language: str = "zh", inventory_nodes: list[dict] | None = None,
               resolution: int = 48, expand_pi: bool = False, preset: str = "current") -> Plan:
    plan = inventory_plan(package, language=language, resolution=resolution, preset=preset)
    for node in nodes:
        hit = node["hit_id"]
        name = artifact_id(hit)
        start = len(plan.commands)
        # Native seek uses genuine pointer/wheel input, including clipped nodes.
        plan.add("seek", hit + ".reveal")
        plan.add("click", hit)
        choice = needs_choice(node, inventory_nodes or nodes)
        if choice:
            if not node.get("choice_hit_id"):
                raise ValueError("Overlapping node has no real chooser hit ID: " + hit)
            plan.add("click", node["choice_hit_id"])
        plan.add("wait", "selection-settle")
        plan.add("inspect", name + "-selected")
        plan.add("click", "diagram.details")
        if expand_pi:
            plan.add("expand", "details.pi.toggle")
        plan.add("inspect-details", name + "-details")
        if capture:
            plan.add("capture", name + "-details-frame")
        # inspect-details ends at the bottom; the lower close button is real UI.
        plan.add("click", "diagram.details.close.bottom")
        plan.attempts.append({"hit_id": hit, "node_id": node.get("id"), "artifact_id": name,
                              "overlap_choice": choice,
                              "pi_expansion_requested": expand_pi,
                              "command_start": start, "command_end": len(plan.commands),
                              "node": node})
    return plan


def read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def run_native(exe: Path, package: Path, plan: Plan, directory: Path,
               timeout: float) -> dict:
    directory.mkdir(parents=True, exist_ok=False)
    cwd = directory / "workdir"
    cwd.mkdir()
    plan_path = directory / "input.plan"
    plan_path.write_text(plan.text(), encoding="utf-8")
    native = directory / "native"
    command = [str(exe), str(package / "canonical.fchk"), "--validation-plan", str(plan_path),
               "--validation-output", str(native), "--validation-background", "--validation-forensic"]
    started = time.monotonic()
    result = {"command": command, "cwd": str(cwd), "native": str(native)}
    kwargs: dict = {}
    if os.name == "nt":
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        startup = subprocess.STARTUPINFO()
        startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        startup.wShowWindow = 0
        kwargs["startupinfo"] = startup
    with (directory / "native.log").open("wb") as log:
        try:
            process = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                                     timeout=timeout, check=False, **kwargs)
            result["returncode"] = process.returncode
        except subprocess.TimeoutExpired:
            result.update(returncode=None, error="Native session timed out", timed_out=True)
        except OSError as exc:
            result.update(returncode=None, error=str(exc))
    result["seconds"] = round(time.monotonic() - started, 3)
    try:
        result["session"] = json.loads((native / "session.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        result["session_error"] = str(exc)
    session = result.get("session", {})
    result["completed"] = (result.get("returncode") == 0 and
                           session.get("capture_status") == "completed" and
                           session.get("failed_commands") == 0)
    write_json(directory / "run.json", result)
    return result


def selection_match(node: dict, view: dict) -> bool:
    snapshots = draw_records(view, "forensic.aomo")
    if snapshots:
        current = snapshots[-1]
        if current.get("selected_side_node_id") == node.get("id"):
            return True
        canonical = node.get("canonical_index")
        if canonical is not None and node.get("lane") == 1:
            active = current.get("active_view") or {}
            return (current.get("focused_canonical_index") == canonical and
                    (active.get("canonical_index") == canonical or
                     view.get("state", {}).get("applied_mo") == canonical))
        return False
    canonical = node.get("canonical_index")
    current = draw_records(view, "forensic.diagram")
    return canonical is not None and bool(current) and current[-1].get("selected_canonical_index") == canonical


def rendered_runs(view: dict, details_only: bool = False) -> list[dict]:
    runs = view.get("rendered", {}).get("text_runs", [])
    if details_only:
        runs = [r for r in runs if "cov.orbital.details" in str(r.get("owner", ""))]
    return runs


def rendered_completeness(view: dict, details_only: bool = False) -> dict:
    rendered = view.get("rendered", {})
    failures = []
    if rendered.get("status") != "captured":
        failures.append("rendered glyph capture is " + str(rendered.get("status")))
    counters = rendered.get("counters", {})
    for key in ("ambiguous_glyphs", "unmatched_font_texture_quads", "callbacks_not_decoded", "invalid_index_groups"):
        if key not in counters:
            failures.append("rendered completeness counter is missing: " + key)
        elif counters[key] != 0:
            failures.append(f"{key}={counters[key]}; complete textual capture is unverified")
    clean_count = sum(clean_run(run) for run in rendered_runs(view, details_only))
    if details_only and clean_count == 0:
        failures.append("No clean, unclipped, uncovered run from the actual details window")
    return {"complete": not failures, "clean_run_count": clean_count, "failures": failures}


def duplicate_gap_pairs(view: dict) -> list[dict]:
    """Flag repeated displayed endpoint pairs; make no scientific inference."""
    by_pair: dict[tuple, list[dict]] = {}
    for gap in draw_records(view, "details.energy-gap"):
        key = (gap.get("gap_kind"), tuple(sorted(gap.get("lower_orbitals", []))),
               tuple(sorted(gap.get("upper_orbitals", []))),gap.get("channel_identity"))
        by_pair.setdefault(key, []).append(gap)
    return [{"observation": "duplicate_displayed_gap_endpoint_pair", "gap_kind": key[0],
             "lower_orbitals": list(key[1]), "upper_orbitals": list(key[2]),
             "channel_identity": key[3],
             "record_count": len(gaps), "labels": [g.get("label") for g in gaps]}
            for key, gaps in by_pair.items() if len(gaps) > 1]


COMPOSITION_BUCKETS = ("centre_current_s", "centre_current_p", "centre_current_d",
                       "centre_current_f", "centre_other", "ligand_valence",
                       "ligand_other", "other_atoms", "core", "unresolved")


def display_semantics(view: dict, require_pi_expanded: bool = False) -> dict:
    """Check the exact frozen graph and details records that produced this draw.

    These are mechanical invariants, not a substitute for chemical positive and
    negative controls. Missing complete source data is not reported as zero.
    """
    failures, checked = [], set()
    snapshots = draw_records(view, "forensic.aomo")
    snapshot = snapshots[-1] if snapshots else {}
    ledgers = [g.get("composition") for g in snapshot.get("group_audit", [])]
    ledgers += draw_records(view, "details.composition")
    for ledger in ledgers:
        if not ledger or not ledger.get("available") or not ledger.get("complete"):
            continue
        checked.add("complete-nao-normalization")
        # Older snapshots predate the explicit disconnected-atom bucket.
        values = [ledger.get(k, 0) if k == "other_atoms" else ledger.get(k)
                  for k in COMPOSITION_BUCKETS]
        if any(not isinstance(v, (int, float)) or not math.isfinite(v) or v < -1e-10 for v in values):
            failures.append("A complete NAO composition has an invalid exclusive bucket")
        else:
            # weight_sum is the measured source norm. Positive missing norm is
            # explicitly retained in unresolved, so exclusive buckets need not
            # equal weight_sum to tighter accuracy than that measured residual.
            error = ledger.get("normalization_error", 0)
            weight_sum = ledger.get("weight_sum", 0)
            if (not isinstance(error, (int, float)) or not math.isfinite(error) or
                    error < 0 or error > 1e-5 or
                    abs(weight_sum - 1) > error + 1e-10 or
                    abs(sum(values) - weight_sum) > error + 1e-10 or
                    abs(sum(values) - 1) > 1e-5 + 1e-10):
                failures.append("Complete NAO exclusive buckets do not close to the original norm")
    counts = snapshot.get("final_counts", {})
    groups = snapshot.get("group_audit", [])
    if groups and counts.get("counts_are_final"):
        checked.add("final-visible-counts")
        included = [g for g in groups if g.get("display_decision", {}).get("included")]
        members = {m for g in included for key in ("member_indices", "member_spin_counterparts")
                   for m in g.get(key, []) if isinstance(m, int)}
        if len(included) != counts.get("groups") or len(members) != counts.get("members"):
            failures.append("Final group/member counts differ from the frozen displayed membership")
        if counts.get("occupied_members", 0) + counts.get("empty_members", 0) != counts.get("members"):
            failures.append("Occupied/empty counts do not partition final members")
        if any(not g.get("display_decision", {}).get("reason_codes") for g in groups):
            failures.append("A retained or folded group has no recorded decision reason")
    detail_compositions = draw_records(view, "details.composition")
    if detail_compositions and not snapshot.get("selected_side_node_id"):
        focus = snapshot.get("focused_canonical_index")
        matching = [n for n in snapshot.get("nodes", []) if n.get("lane") == 1 and n.get("canonical_index") == focus]
        if matching and matching[0].get("composition"):
            checked.add("graph-details-same-composition")
            if matching[0]["composition"] != detail_compositions[-1]:
                failures.append("Graph and actual selected details use different composition ledgers")
    brief = snapshot.get("view", {}).get("preset") == 0
    for gap in draw_records(view, "details.energy-gap"):
        mode = gap.get("mode")
        if mode is None:
            continue  # Crystal-field and older records do not claim this check.
        checked.add("shared-mode-and-direction")
        if not mode.get("verified") or not mode.get("shared_mode_ids") or not mode.get("shared_fragment_contraction_norm_hartree", 0) > 0:
            failures.append("A displayed pi relation has no verified common fragment mode")
        if not mode.get("two_endpoint_relation"):
            failures.append("A multi-group mode was displayed as a two-endpoint relation")
        if mode.get("direction") in ("centre_to_ligand", "ligand_to_centre") and (
                not mode.get("direction_verified") or not mode.get("matched_edge_ids")):
            failures.append("Displayed direction lacks a linked ordered source relation")
        if brief and not gap.get("ordinary_display_eligible"):
            failures.append("A secondary or out-of-scope pi relation leaked into brief details")
    network_ids = set()
    for record in draw_records(view, "details.pi-network"):
        network = record.get("network", {})
        checked.add("one-network-per-common-mode")
        identity = (network.get("channel_id"), network.get("mode_id"))
        if identity in network_ids:
            failures.append("The same common-mode network was displayed more than once")
        network_ids.add(identity)
        if not network.get("verified") or not network.get("nodes") or not network.get("mode_id"):
            failures.append("A displayed network lacks a verified mode or actual member groups")
        if brief and not network.get("ordinary_display_eligible"):
            failures.append("An auxiliary network leaked into brief details")
        if network.get("direction") in ("centre_to_ligand", "ligand_to_centre", "bidirectional") and (
                not network.get("direction_verified") or not network.get("matched_edge_ids")):
            failures.append("A displayed network direction lacks linked ordered source evidence")
    target_ids = {t.get("id") for t in view.get("targets", [])}
    for record in draw_records(view, "details.pi.joined-summary"):
        assessment = record.get("assessment", {})
        checked.add("joined-occupied-pi-summary")
        if not assessment.get("available") or not assessment.get("applicable"):
            failures.append("A joined pi summary has no applicable group response")
        if "回馈" in record.get("label", "") and not assessment.get("occupied_metal_backbond_supported"):
            failures.append("A displayed back-donation claim has no occupied metal donor")
        if not assessment.get("members"):
            failures.append("A joined pi summary is detached from source group members")
    if require_pi_expanded and "details.pi.toggle.closed" in target_ids:
        failures.append("Applicable pi details remain closed; their text was not inspected")
    return {"checked": sorted(checked), "failures": sorted(set(failures)),
            "status": "failed" if failures else "passed" if checked else "not_applicable"}


BONDING_CAPTIONS = ("All-pair bonding character", "各原子对的综合成键作用", "全原子対の結合性",
                    "Caractère liant de l’ensemble des paires")


def observed_bonding_labels(view: dict) -> list[dict]:
    """Read the actually rendered caption/value row, using its final geometry."""
    runs = [r for r in rendered_runs(view, details_only=True) if clean_run(r)]
    result = []
    for run in runs:
        text = run.get("text", "")
        caption = next((c for c in BONDING_CAPTIONS if c in text), None)
        if caption is None:
            continue
        suffix = text.split(caption, 1)[1].lstrip(":： ")
        value_run = run
        if not suffix:
            bounds = run.get("bounds", [])
            if len(bounds) != 4:
                continue
            nearby = [r for r in runs if r is not run and r.get("owner") == run.get("owner")
                      and len(r.get("bounds", [])) == 4 and r["bounds"][0] >= bounds[2] - 1
                      and abs(r.get("baseline_origin_y", r["bounds"][1]) -
                              run.get("baseline_origin_y", bounds[1])) <= 2]
            if not nearby:
                continue
            value_run = min(nearby, key=lambda r: r["bounds"][0])
            suffix = value_run.get("text", "").strip()
        if suffix:
            result.append({"caption": caption, "value": suffix, "caption_bounds": run.get("bounds"),
                           "value_bounds": value_run.get("bounds"), "source": "actual_clean_rendered_runs"})
    return result


def evaluate_attempt(attempt: dict, native: Path, actions: list[dict],
                     canonical: Path | None = None) -> dict:
    result = {k: v for k, v in attempt.items() if k != "node"}
    by_command = {a.get("command"): a for a in actions}
    expected = range(attempt["command_start"], attempt["command_end"])
    result["actions"] = [by_command[i] for i in expected if i in by_command]
    failures = [f"command {i}: missing action" for i in expected if i not in by_command]
    failures += [f"command {a.get('command')}: {a.get('detail') or a.get('status')}"
                 for a in result["actions"] if a.get("status") != "executed"]
    name = attempt["artifact_id"]
    selected_path = native / (name + "-selected.view.json")
    pages = sorted(native.glob(name + "-details-p*.view.json"))
    result.update(selected_capture=selected_path.name if selected_path.is_file() else None,
                  details_pages=[p.name for p in pages], details_page_count=len(pages))
    selected = None
    textual_complete = True
    try:
        selected = json.loads(selected_path.read_text(encoding="utf-8"))
        if selected.get("state", {}).get("scene_matches_applied") is not True:
            failures.append("Selected capture scene does not match the applied identity")
        selected_text = rendered_completeness(selected)
        result["selected_text_completeness"] = selected_text
        textual_complete &= selected_text["complete"]
        failures.extend("Selected capture: " + item for item in selected_text["failures"])
        if canonical is not None:
            identity = input_identity(selected, canonical)
            result["input_identity"] = identity
            if not identity["matched"]:
                failures.append(identity["failure"])
        matched = selection_match(attempt["node"], selected)
        result["selection_matched"] = matched
        if not matched:
            failures.append("Selected actual UI identity does not match the requested source node")
    except (OSError, ValueError) as exc:
        failures.append("Selected inspection unavailable: " + str(exc))
        result["selection_matched"] = False
        textual_complete = False
    if not pages:
        failures.append("No actual details page was captured")
        textual_complete = False
    run_count = 0
    observations, bonding = [], []
    page_completeness, semantic_checks = [], []
    for page in pages:
        try:
            value = json.loads(page.read_text(encoding="utf-8"))
            if value.get("state", {}).get("scene_matches_applied") is not True:
                failures.append(page.name + ": scene does not match the applied identity")
            if canonical is not None:
                identity = input_identity(value, canonical)
                if not identity["matched"]:
                    failures.append(page.name + ": " + identity["failure"])
            detail = draw_records(value, "forensic.diagram.details")
            if not detail or not detail[-1].get("open") or not detail[-1].get("content_visible"):
                failures.append(page.name + ": details context is missing or closed")
            elif selected is not None:
                aomo = draw_records(selected, "forensic.aomo")
                diagram = draw_records(selected, "forensic.diagram")
                active = (aomo[-1] if aomo else diagram[-1] if diagram else {}).get("active_view")
                if active != detail[-1].get("active_view"):
                    failures.append(page.name + ": details active identity differs from selected capture")
            complete = rendered_completeness(value, details_only=True)
            semantics = display_semantics(value, attempt.get("pi_expansion_requested", False))
            semantic_checks.append(dict(semantics, page=page.name))
            failures.extend(page.name + ": " + item for item in semantics["failures"])
            page_completeness.append(dict(complete, page=page.name))
            textual_complete &= complete["complete"]
            failures.extend(page.name + ": " + item for item in complete["failures"])
            run_count += len(rendered_runs(value, details_only=True))
            observations.extend(dict(item, page=page.name) for item in duplicate_gap_pairs(value))
            bonding.extend(dict(item, page=page.name) for item in observed_bonding_labels(value))
        except (OSError, ValueError) as exc:
            failures.append(page.name + ": " + str(exc))
            textual_complete = False
    if pages and run_count == 0:
        failures.append("Details pages contain no decoded runs from the real details window")
    result.update(details_text_run_count=run_count, observed_duplicate_gap_pairs=observations,
                  textual_capture_complete=textual_complete, details_text_completeness=page_completeness,
                  display_semantics=semantic_checks,
                  observed_bonding_labels=bonding, failures=failures,
                  status="captured" if not failures else "failed")
    return result


def group_label_differences(nodes: list[dict], attempts: list[dict]) -> list[dict]:
    by_hit = {n["hit_id"]: n for n in nodes}
    groups: dict[str, dict[str, list[str]]] = {}
    for attempt in attempts:
        if attempt.get("status") != "captured":
            continue
        node = by_hit[attempt["hit_id"]]
        group = node.get("display_group_id") or node.get("subspace_id")
        if not group and node.get("row") is not None:
            group = "canonical-row:" + str(node["row"])
        if not group:
            continue
        values = sorted({item["value"] for item in attempt.get("observed_bonding_labels", [])})
        if values:
            groups.setdefault(group, {})[attempt["hit_id"]] = values
    return [{"observation": "differing_rendered_bonding_labels_within_displayed_group",
             "group_id": group, "members": members,
             "scientific_verdict": "not assessed; observed UI difference only"}
            for group, members in groups.items() if len(members) > 1 and
            len({tuple(values) for values in members.values()}) > 1]


def summarize_coverage(nodes: list[dict], attempts: list[dict], sweep: bool) -> dict:
    included = {n["hit_id"] for n in nodes if not n.get("group_header")}
    attempted = {a["hit_id"] for a in attempts}
    captured = {a["hit_id"] for a in attempts if a.get("status") == "captured"}
    return {"scope": "source_nodes_in_actual_initial_view; folded_headers_separate",
            "included_source_nodes": len(included), "included_group_headers": sum(bool(n.get("group_header")) for n in nodes),
            "attempted": len(attempted & included), "captured": len(captured & included),
            "failed": len((attempted - captured) & included), "not_attempted": sorted(included - attempted),
            "captured_fraction": len(captured & included) / len(included) if included else None,
            "complete": bool(sweep and included and captured >= included)}


def clean_run(run: dict) -> bool:
    return (run.get("clipping") == "unclipped" and
            run.get("occlusion") == "no_later_window_overlap" and not run.get("ambiguous"))


def write_rendered_report(native: Path, destination: Path) -> None:
    """Only decode recorded final draw runs; never replace them with metadata."""
    qualified_path = destination.with_name(destination.stem + "-qualified.txt")
    with destination.open("w", encoding="utf-8") as out, qualified_path.open("w", encoding="utf-8") as qualified:
        for path in sorted(native.glob("*.view.json")):
            try:
                view = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                out.write(f"\n{path.name}: unreadable: {exc}\n")
                continue
            details = "-details-p" in path.name
            runs = rendered_runs(view, details_only=details)
            rendered = view.get("rendered", {})
            out.write(f"\n=== {path.name} | {rendered.get('status')} | {len(runs)} actual runs ===\n")
            qualified.write(f"\n=== {path.name} | clipped, ambiguous, or possibly covered actual runs ===\n")
            for run in runs:
                target = out if clean_run(run) else qualified
                target.write(f"[{run.get('owner')} | {run.get('bounds')} | {run.get('clipping')} | "
                             f"{run.get('occlusion')}] {run.get('text', '')}\n")


def find_cases(library: Path, requested: list[str] | None) -> list[tuple[str, Path]]:
    library = library.resolve()
    base = library.parent if library.is_file() else library
    manifest = library if library.is_file() else library / "library.json"
    if manifest.is_file():
        rows = json.loads(manifest.read_text(encoding="utf-8")).get("cases", [])
        available = {r["case_id"]: (base / r.get("package", "packages/" + r["case_id"])).resolve() for r in rows}
    else:
        packages = library / "packages" if (library / "packages").is_dir() else library
        available = {p.name: p.resolve() for p in packages.iterdir() if p.is_dir() and (p / "canonical.fchk").is_file()}
    # Package folders can acquire readable suffixes without changing source IDs.
    # Resolve a requested stable ID only when its actual folder is unambiguous.
    packages = base / "packages" if (base / "packages").is_dir() else base
    for cid in requested or list(available):
        if cid in available and (available[cid] / "canonical.fchk").is_file():
            continue
        matches = [p.resolve() for p in packages.iterdir() if p.is_dir() and
                   p.name.startswith(cid + "-") and (p / "canonical.fchk").is_file()]
        if len(matches) == 1:
            available[cid] = matches[0]
    ids = list(dict.fromkeys(requested)) if requested else sorted(available)
    absent = [cid for cid in ids if cid not in available]
    if absent:
        raise ValueError("Case IDs absent from library: " + ", ".join(absent))
    result = [(cid, available[cid]) for cid in ids]
    for cid, package in result:
        if not (package / "canonical.fchk").is_file():
            raise FileNotFoundError(f"{cid}: canonical.fchk missing in {package}")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", nargs="+", help="Case IDs; comma-separated IDs also accepted")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--inventory", action="store_true", help="Inventory only (default)")
    modes.add_argument("--sweep", action="store_true", help="Inspect included source nodes in one batch per case")
    modes.add_argument("--all-nodes", action="store_true", help="Explicit sweep of all included source nodes")
    parser.add_argument("--selected", nargs="+", help="Only these stable hit IDs; implies a sweep")
    parser.add_argument("--limit", type=int, help="Maximum cases to collect")
    parser.add_argument("--node-limit", type=int, help="Maximum source nodes per case; coverage records omissions")
    parser.add_argument("--capture", action="store_true", help="Also preserve native BMP frame captures")
    parser.add_argument("--expand-pi", action="store_true", help="Use real input to open applicable pi details before each capture")
    parser.add_argument("--preset", choices=("current", "teaching", "research", "full"), default="current",
                        help="Select an available AO/MO view through its actual combo")
    parser.add_argument("--language", choices=("en", "zh", "ja", "fr", "current"), default="zh",
                        help="Select via the actual language combo (default zh)")
    parser.add_argument("--resolution", type=int, default=48, help="Actual scene mesh setting (default 48; no scientific surface acceptance)")
    parser.add_argument("--timeout", type=float, default=3600, help="Seconds allowed per native process")
    parser.add_argument("--zip", nargs="?", const="auto", help="ZIP all output after collection; original files retained")
    args = parser.parse_args(argv)
    if args.inventory and args.selected:
        parser.error("--selected implies a sweep and cannot accompany --inventory")
    if args.timeout <= 0 or any(v is not None and v <= 0 for v in (args.limit, args.node_limit)):
        parser.error("Timeout and limits must be positive")
    if args.resolution < 16 or args.resolution > 512:
        parser.error("Resolution must be in 16..512")
    exe = args.exe.resolve()
    if not exe.is_file():
        parser.error("Native validation executable does not exist")
    output = args.output.resolve()
    base = args.library.resolve()
    if base.is_file():
        base = base.parent
    if output == base or base in output.parents:
        parser.error("Output must be outside the immutable library")
    if output.exists():
        parser.error("Output already exists; use a new collection directory")
    requested = [cid for item in (args.cases or []) for cid in item.split(",") if cid] or None
    cases = find_cases(args.library, requested)
    if args.limit:
        cases = cases[:args.limit]
    output.mkdir(parents=True)
    sweep = bool(args.sweep or args.all_nodes or args.selected)
    summary = {"schema": "cov.display-forensics.collection.v1", "exe": str(exe),
               "library": str(args.library.resolve()), "mode": "sweep" if sweep else "inventory",
               "language": args.language, "scene_resolution": args.resolution, "requested_preset": args.preset,
               "cases": [], "scientific_verdict": "not assessed; actual UI evidence only"}
    any_failed = False
    for cid, package in cases:
        print(f"{cid}: collecting actual UI inventory", flush=True)
        folder = output / artifact_id(cid)
        folder.mkdir()
        case = {"case_id": cid, "package": str(package), "directory": folder.name, "attempts": []}
        nodes = []
        try:
            inv = run_native(exe, package, inventory_plan(package, args.capture, args.language, args.resolution,args.preset), folder / "inventory", args.timeout)
            case["inventory_run"] = inv
            view = json.loads((folder / "inventory/native/inventory.view.json").read_text(encoding="utf-8"))
            identity = input_identity(view, package / "canonical.fchk")
            case["input_identity"] = identity
            if not identity["matched"]:
                raise RuntimeError(identity["failure"])
            nodes, context = node_inventory(view)
            case["display_semantics"] = display_semantics(view)
            if case["display_semantics"]["failures"]:
                raise RuntimeError("; ".join(case["display_semantics"]["failures"]))
            case["inventory"] = context
            write_json(folder / "nodes.json", {"context": context, "nodes": nodes})
            write_rendered_report(folder / "inventory/native", folder / "inventory-rendered.txt")
            if not inv["completed"]:
                raise RuntimeError("Inventory native session did not complete successfully")
            if sweep:
                chosen = choose_nodes(nodes, args.selected, args.node_limit)
                plan = sweep_plan(package, chosen, args.capture, args.language, nodes, args.resolution,args.expand_pi,args.preset)
                print(f"{cid}: batch-inspecting {len(chosen)} source nodes (including clipped nodes)", flush=True)
                run = run_native(exe, package, plan, folder / "sweep", args.timeout)
                case["sweep_run"] = run
                actions = read_jsonl(folder / "sweep/native/actions.jsonl")
                case["attempts"] = [evaluate_attempt(a, folder / "sweep/native", actions, package / "canonical.fchk") for a in plan.attempts]
                write_rendered_report(folder / "sweep/native", folder / "sweep-rendered.txt")
                if not run["completed"]:
                    case["run_failure"] = "Sweep native session failed, timed out, or contains failed commands"
        except (OSError, ValueError, RuntimeError) as exc:
            case["error"] = str(exc)
        case["coverage"] = summarize_coverage(nodes, case["attempts"], sweep)
        case["observed_group_label_differences"] = group_label_differences(nodes, case["attempts"])
        native = folder / ("sweep/native" if (folder / "sweep/native").is_dir() else "inventory/native")
        if native.is_dir():
            try:
                case["display_bundle"] = compact_bundle(native, folder / (artifact_id(cid) + ".display.json.gz"))
                print(f"{cid}: primary lossless bundle {case['display_bundle']['compressed_bytes']} bytes", flush=True)
            except (OSError, ValueError) as exc:
                case["bundle_error"] = str(exc)
        case["status"] = "failed" if case.get("error") or case.get("run_failure") or case.get("bundle_error") or case["coverage"]["failed"] else "collected"
        any_failed |= case["status"] == "failed"
        write_json(folder / "case.json", case)
        summary["cases"].append(case)
        write_json(output / "collection.json", summary)
        print(f"{cid}: {case['status']}; {case['coverage']['captured']}/{case['coverage']['included_source_nodes']} source nodes captured", flush=True)
    summary["status"] = "failed" if any_failed else "collected"
    write_json(output / "collection.json", summary)
    report = ["Actual native display collection", "", f"Mode: {summary['mode']}",
              "Coverage applies to source nodes included in the initial actual view; folded headers are listed separately.",
              "Rendered text is decoded from final draw quads; raster legibility and scientific correctness are not assessed.", ""]
    report.insert(5, f"Actual scene resolution: {args.resolution}; no scientific isosurface acceptance is claimed.")
    for case in summary["cases"]:
        c = case["coverage"]
        report.append(f"{case['case_id']}: {case['status']}; included={c['included_source_nodes']}, attempted={c['attempted']}, captured={c['captured']}, failed={c['failed']}, not_attempted={len(c['not_attempted'])}, complete={c['complete']}")
        for attempt in case["attempts"]:
            if attempt["failures"]:
                report.append("  " + attempt["hit_id"] + ": " + "; ".join(attempt["failures"]))
        if case.get("error") or case.get("run_failure"):
            report.append("  " + str(case.get("error") or case.get("run_failure")))
        if case.get("display_bundle"):
            bundle = case["display_bundle"]
            report.append(f"  PRIMARY lossless display bundle: {Path(bundle['file']).name}; {bundle['compressed_bytes']} bytes from {bundle['source_bytes']} staging bytes; {bundle['shared_objects']} shared objects.")
        if case.get("bundle_error"):
            report.append("  Compact bundle failed: " + case["bundle_error"])
        for observation in case["observed_group_label_differences"]:
            report.append("  OBSERVED group label difference: " + json.dumps(observation, ensure_ascii=False))
        duplicate_count = sum(len(a.get("observed_duplicate_gap_pairs", [])) for a in case["attempts"])
        if duplicate_count:
            report.append(f"  OBSERVED repeated endpoint-pair records: {duplicate_count} page observations (see case.json).")
    (output / "REPORT.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    if args.zip:
        archive = output.with_suffix(".zip") if args.zip == "auto" else Path(args.zip).resolve()
        if archive.exists() or output in archive.parents:
            raise ValueError("Archive must be new and outside the collection directory")
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
            for path in sorted(output.rglob("*")):
                if path.is_file():
                    bundle.write(path, path.relative_to(output).as_posix())
        print("Compressed evidence: " + str(archive), flush=True)
    return 1 if any_failed else 0


if __name__ == "__main__":
    sys.exit(main())
