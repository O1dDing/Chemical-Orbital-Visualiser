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

"""Immutable input registry and resumable, bounded original-library recomputation.

This module prepares and collects Gaussian evidence. Scientific adjudication is
deliberately a separate post-batch stage. 'collected' is never a scientific pass.
Only the approved 273 registry entries are used; failed entries remain present.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import sys
import threading
import time

from validation_process import atomic_json, physical_core_masks, pin_current, run_tree

SCHEMA = 1
DEFAULT_COLLECTION = Path(r"F:\Codex\2026-09-05\branch-15\outputs\fchk-collection-20260905")
DEFAULT_INPUTS = Path(r"F:\CalChem\Gaussian\Files\.gjf")
DEFAULT_OUTPUT = Path(r"F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906")
DEFAULT_JOBS = Path(r"F:\Dev\cov-cycle-20260906\jobs")
DEFAULT_GAUSSIAN = Path(r"F:\CalChem\Gaussian\Gaussian 16 W")
BOHR_ANGSTROM = 0.529177210903
FIELDS = {"Number of atoms", "Charge", "Multiplicity", "Number of electrons",
          "Number of alpha electrons", "Number of beta electrons",
          "Number of basis functions", "Number of independent functions",
          "Atomic numbers", "Nuclear charges", "Current cartesian coordinates",
          "Shell types", "Number of primitives per shell", "Shell to atom map",
          "Primitive exponents", "Contraction coefficients", "P(S=P) Contraction coefficients",
          "Total Energy", "SCF Energy", "ECP-RNFroz"}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        while block := f.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def read_fchk(path):
    """Read selected fixed-width fields without retaining the large MO arrays."""
    values = {}
    with Path(path).open(encoding="ascii") as stream:
        title = next(stream).strip()
        header = next(stream).rstrip()
        for line in stream:
            if not line.strip():
                continue
            match = re.match(r"^(.{40})\s+([IRCLH])\s+(.*)$", line.rstrip())
            if not match:
                raise ValueError(f"Invalid FCHK header: {line[:100]}")
            name, kind, body = match.groups()
            name = name.strip()
            if body.startswith("N="):
                count = int(body[2:].strip())
                per_line = {"I": 6, "R": 5, "C": 5, "H": 9, "L": 72}[kind]
                rows = math.ceil(count / per_line)
                if name in FIELDS:
                    if kind not in "IR":
                        raise ValueError(f"Unexpected field type: {name}")
                    raw = " ".join(next(stream) for _ in range(rows)).replace("D", "E").split()
                    if len(raw) != count:
                        raise ValueError(f"Invalid FCHK array length: {name}")
                    parsed = [int(x) if kind == "I" else float(x) for x in raw]
                    if not all(math.isfinite(x) for x in parsed):
                        raise ValueError(f"Nonfinite FCHK data: {name}")
                    values[name] = parsed
                else:
                    for _ in range(rows):
                        next(stream)
                if name == "Beta MO coefficients":
                    values["has_beta_coefficients"] = True
            elif name in FIELDS:
                values[name] = int(body) if kind == "I" else float(body.replace("D", "E"))
    values.update(title=title, job_type=header[:10].strip(),
                  method=header[10:70].strip(), basis=header[70:].strip())
    required = ("Atomic numbers", "Current cartesian coordinates", "Charge", "Multiplicity",
                "Number of alpha electrons", "Number of beta electrons", "Shell types")
    if any(key not in values for key in required):
        raise ValueError(f"Required FCHK identity absent: {path}")
    return values


def read_gjf(path, atom_count):
    raw = Path(path).read_text(encoding="utf-8-sig", errors="strict")
    sections = re.split(r"(?im)^\s*--link1--\s*$", raw)
    lines = sections[0].splitlines()
    route_start = next(i for i, row in enumerate(lines) if row.lstrip().startswith("#"))
    end = route_start
    while end < len(lines) and lines[end].strip():
        end += 1
    route = " ".join(row.strip() for row in lines[route_start:end])
    index = end
    while index < len(lines) and not lines[index].strip():
        index += 1
    title = []
    while index < len(lines) and lines[index].strip():
        title.append(lines[index].strip())
        index += 1
    while index < len(lines) and not lines[index].strip():
        index += 1
    charge_mult = lines[index].split()
    if len(charge_mult) != 2:
        raise ValueError(f"Unexpected charge/multiplicity record in {path}")
    charge, mult = map(int, charge_mult)
    coords = lines[index + 1:index + 1 + atom_count]
    if len(coords) != atom_count or any(len(row.split()) != 4 for row in coords):
        raise ValueError(f"Non-Cartesian original input requires explicit support: {path}")
    end_coords = index + atom_count + 1
    if end_coords < len(lines) and lines[end_coords].strip():
        raise ValueError(f"Original atom count mismatch: {path}")
    extra = "\n".join(lines[end_coords:]).strip()
    method_basis = re.search(r"(?:^|\s)([A-Za-z][A-Za-z0-9+-]*)/([^\s]+)", route)
    if not method_basis:
        raise ValueError(f"No unambiguous method/basis token: {path}")
    return {"route": route, "title": " ".join(title), "charge": charge,
            "multiplicity": mult, "method": method_basis.group(1),
            "basis": method_basis.group(2), "extra": extra,
            "link_sections": len(sections)}


def keyword(route, name):
    match = re.search(r"(?i)(?<!\S)" + name + r"(?:\s*=\s*(\([^)]*\)|[^\s]+))?(?=\s|$)", route)
    return match.group().strip() if match else None


def reproduction_route(source, fields, retry=False):
    basis = source["basis"]
    # The FCHK method includes the actual R/U/RO state, including U singlets.
    parts = ["#p", fields["method"] + "/" + basis, "SP", "Units=Bohr"]
    route = source["route"]
    for key in ("5D", "6D", "7F", "10F", "9G", "15G", "EmpiricalDispersion",
                "Integral", "Int", "Symmetry", "NoSymm", "SCRF"):
        token = keyword(route, key)
        if token:
            parts.append(token)
    # Tighten iteration convergence without changing the Hamiltonian. The sole
    # retry extends the same XQC algorithm; no state, geometry or basis edits.
    parts += ["SCF=(VeryTight,XQC,MaxCycle=" + ("2048" if retry else "1024") + ")", "Pop=Full"]
    guess = keyword(route, "Guess")
    if guess and "mix" in guess.lower() and "read" not in guess.lower():
        parts.append("Guess=Mix")
    for forbidden in ("pseudo", "geom", "charge", "external", "oniom", "field", "iop"):
        token = keyword(route, forbidden)
        if token and forbidden != "iop":
            raise ValueError(f"Original option needs explicit reproduction support: {token}")
        # Known GFInput IOp(6/7=3) controls printing only, not the Hamiltonian.
        if token and re.sub(r"\s+", "", token).lower() not in ("iop(6/7=3)",):
            raise ValueError(f"Unreviewed internal option: {token}")
    return " ".join(parts)


def input_text(case, retry=False, cores=4, memory_gib=24):
    f, source = case["fchk_identity"], case["original_input"]
    route = reproduction_route(source, f, retry)
    coordinates = f["Current cartesian coordinates"]
    lines = ["%chk=job.chk", f"%mem={memory_gib}GB", f"%nprocshared={cores}", route, "",
             case["case_id"] + " independent reproduction at original final geometry", "",
             f'{f["Charge"]} {f["Multiplicity"]}']
    for index, z in enumerate(f["Atomic numbers"]):
        xyz = coordinates[3 * index:3 * index + 3]
        lines.append(str(z) + " " + " ".join(f"{x:.14e}" for x in xyz))
    lines.append("")
    general_basis = source["basis"].lower() in ("gen", "genecp")
    if general_basis:
        if "****" not in source["extra"]:
            raise ValueError("General basis missing from original input")
        lines.append(source["extra"])
        lines.append("")
    elif source["extra"]:
        raise ValueError("Uninterpreted additional input cannot be silently discarded")
    return "\n".join(lines) + "\n"


def purpose(relative_path, source, f):
    path = relative_path.lower()
    if "scan_master" in path:
        return {"kind": "scan_frame", "reference_geometry": "fixed",
                "basis": "Original registered scan frame; no minimum requirement"}
    if "geometry_matrix" in path:
        return {"kind": "coordination_geometry_candidate", "reference_geometry": "optimize",
                "basis": "Synthetic starting shell; target name is not a physical geometry verdict",
                "additional_reference": "original fixed geometry diagnostic"}
    fixed_markers = ("stretched", "twisted", "90deg", "square", "_planar", "eclipsed",
                     "separated", "negative")
    if any(token in path for token in fixed_markers):
        return {"kind": "specified_geometry", "reference_geometry": "fixed",
                "basis": "Explicit diagnostic geometry in registered input name",
                "needs_post_batch_purpose_review": True}
    if "nonplanar_start" in path:
        return {"kind": "equilibrium_candidate", "reference_geometry": "optimize",
                "basis": "Original nonplanar start is a convergence diagnostic"}
    return {"kind": "equilibrium_candidate", "reference_geometry": "optimize",
            "basis": "Candidate equilibrium at the specified charge and multiplicity",
            "competing_states_required": bool(f.get("has_beta_coefficients") or
                any(21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 for z in f["Atomic numbers"]))}


def freeze(args):
    path = args.output / "campaign.json"
    if path.exists():
        raise RuntimeError("Campaign already frozen; use status/resume, never replace its registry")
    input_index = {}
    for item in args.inputs.rglob("*"):
        if item.suffix.lower() in (".gjf", ".com"):
            input_index.setdefault(item.stem.casefold(), []).append(item)
    cases = []
    for directory in sorted((args.collection / "cases").glob("OLD-*")):
        identity = json.loads((directory / "source-identity.json").read_text(encoding="utf-8"))
        matches = input_index.get(Path(identity["input"]).stem.casefold(), [])
        if len(matches) != 1:
            raise ValueError(f'{directory.name}: expected one original GJF, found {len(matches)}')
        f = read_fchk(directory / "source.fch")
        gjf = read_gjf(matches[0], f["Number of atoms"])
        if (gjf["charge"], gjf["multiplicity"]) != (f["Charge"], f["Multiplicity"]):
            raise ValueError(f'{directory.name}: source charge/multiplicity mismatch')
        case = {"case_id": directory.name, "original_relative_path": identity["relative_path"],
                "original_fchk": str(directory / "source.fch"),
                "original_fchk_sha256": sha256(directory / "source.fch"),
                "original_gjf": str(matches[0]), "original_gjf_sha256": sha256(matches[0]),
                "original_logs": [str(directory / row["copy"]) for row in identity.get("logs", [])],
                "fchk_identity": f, "original_input": gjf,
                "purpose": purpose(identity["relative_path"], gjf, f)}
        # Validate every generated input before any job in this batch starts.
        case["reproduction_input_sha256"] = hashlib.sha256(input_text(case).encode("ascii")).hexdigest()
        case["case_identity"] = digest(case)
        cases.append(case)
    expected = [f"OLD-{i:03d}" for i in range(1, 274)]
    if [case["case_id"] for case in cases] != expected:
        raise ValueError("Original registry must contain OLD-001 through OLD-273 exactly once")
    for case in cases:
        directory = args.output / "inputs" / case["case_id"]
        directory.mkdir(parents=True, exist_ok=True)
        shutil.copy2(case["original_gjf"], directory / "original.gjf")
        (directory / "reproduction.gjf").write_text(input_text(case), encoding="ascii", newline="\n")
        atomic_json(directory / "identity.json", case)
    campaign = {"schema_version": SCHEMA, "created_epoch": time.time(),
                "name": "COV complete validation 2026-09-06", "original_case_count": 273,
                "external_case_count": 0, "pilot_case_id": "OLD-001", "cases": cases,
                "case_set_identity": digest([case["case_identity"] for case in cases]),
                "resource_policy": {"total_cores": 12, "total_memory_gib": 128,
                    "gaussian_cores": 4, "gaussian_memory_gib": 24,
                    "gaussian_tree_memory_limit_gib": 32, "gaussian_concurrent_jobs": 2,
                    "reserved_other_cores": 4, "reserved_other_memory_gib": 64},
                "reproduction_policy": {"max_attempts": 2, "timeout_hours_per_attempt": 24,
                    "first": "SCF VeryTight XQC MaxCycle=1024; original fresh guess policy",
                    "retry": "SCF convergence failure only; fresh guess; MaxCycle=2048",
                    "round_barrier": "All 273 terminal before scientific adjudication",
                    "collected_is_pass": False},
                "stages": {name: "pending" for name in ("input_reliability", "general_fixes",
                     "original_273_loop", "external_set", "complete_loop", "delivery")}}
    atomic_json(path, campaign)
    atomic_json(args.output / "progress.json", {"schema_version": 1,
                "case_set_identity": campaign["case_set_identity"], "updated_epoch": time.time(),
                "stage": "reproduction_prepared", "counts": {"pending": 273},
                "scientifically_passed": 0, "all_required_evidence_complete": False})
    print(json.dumps({"registry": str(path), "cases": len(cases),
                      "purpose_counts": dict(Counter(x["purpose"]["kind"] for x in cases))}))


def log_completion(path):
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    endings = list(re.finditer(r"Normal termination|Error termination", text))
    return {"normal": bool(endings and endings[-1].group() == "Normal termination"),
            "termination": endings[-1].group() if endings else "missing",
            "scf_convergence_failure": bool(re.search(
                r"Convergence failure|SCF has not converged|Convergence criterion not met", text, re.I)),
            "last_scf_energy": (re.findall(r"SCF Done:\s+E\([^)]+\)\s*=\s*([-+\d.DE]+)", text) or [None])[-1],
            "stability_messages": [row.strip() for row in text.splitlines()
                                   if "wavefunction is stable" in row.lower() or "instability" in row.lower()][-12:]}


def source_identity_check(case, result):
    original = case["fchk_identity"]
    keys = ("Atomic numbers", "Charge", "Multiplicity", "Number of alpha electrons",
            "Number of beta electrons", "Shell types", "Shell to atom map",
            "Number of primitives per shell", "Number of basis functions")
    mismatches = [key for key in keys if result.get(key) != original.get(key)]
    for key in ("Nuclear charges", "Primitive exponents", "Contraction coefficients",
                "P(S=P) Contraction coefficients"):
        a, b = original.get(key, []), result.get(key, [])
        if len(a) != len(b) or any(abs(x - y) > 2e-8 * max(1, abs(x)) for x, y in zip(a, b)):
            mismatches.append(key)
    a, b = original["Current cartesian coordinates"], result["Current cartesian coordinates"]
    distance_error = 0.0
    if len(a) == len(b):
        for i in range(0, len(a), 3):
            for j in range(0, i, 3):
                da = math.dist(a[i:i + 3], a[j:j + 3])
                db = math.dist(b[i:i + 3], b[j:j + 3])
                distance_error = max(distance_error, abs(da - db))
        if distance_error > 2e-6:
            mismatches.append("pair_distances_bohr")
    else:
        mismatches.append("Current cartesian coordinates")
    return {"input_identity_consistent": not mismatches, "mismatches": mismatches,
            "max_pair_distance_delta_bohr": distance_error,
            "note": "Identity gate only; energy/density/subspace comparison awaits batch barrier"}


def collect_case(case, args, mask, runner_identity):
    base = args.jobs / "R001" / case["case_id"] / "reproduction"
    base.mkdir(parents=True, exist_ok=True)
    state_path = base / "state.json"
    terminal = {"collected", "failed", "identity_mismatch", "timeout"}
    previous = json.loads(state_path.read_text()) if state_path.exists() else None
    if previous:
        if previous["case_identity"] != case["case_identity"] or previous["runner_identity"] != runner_identity:
            raise RuntimeError("Frozen case/runner identity changed; create a new round instead of reusing evidence")
        if previous["status"] in terminal:
            return previous
    state = previous or {"case_id": case["case_id"], "case_identity": case["case_identity"],
                         "runner_identity": runner_identity, "attempts": [], "status": "pending"}
    if state["status"] == "running":
        # The previous exclusive supervisor has exited; KILL_ON_JOB_CLOSE ended
        # its descendants. Preserve the interrupted attempt and use a fresh dir.
        state["attempts"][-1]["status"] = "interrupted"
    try:
        if sha256(Path(case["original_fchk"])) != case["original_fchk_sha256"]:
            raise ValueError("Frozen original FCHK changed")
        gjf_copy = args.output / "inputs" / case["case_id"] / "original.gjf"
        if sha256(gjf_copy) != case["original_gjf_sha256"]:
            raise ValueError("Frozen original Gaussian input changed")
        completed_attempts = [x for x in state["attempts"] if x["status"] != "interrupted"]
        retry = len(completed_attempts) > 0
        for _ in range(2 - len(completed_attempts)):
            number = len(state["attempts"]) + 1
            directory = base / f"attempt-{number:02d}"
            directory.mkdir(exist_ok=False)
            scratch = directory / "scratch"
            scratch.mkdir()
            gjf = directory / "job.gjf"
            gjf.write_text(input_text(case, retry=retry), encoding="ascii", newline="\n")
            attempt = {"number": number, "status": "running", "input_sha256": sha256(gjf),
                       "fresh_initial_guess": True, "directory": str(directory),
                       "predeclared_scf_retry": retry}
            state["attempts"].append(attempt)
            state["status"] = "running"
            atomic_json(state_path, state)
            env = dict(os.environ)
            env.update(GAUSS_EXEDIR=str(args.gaussian), GAUSS_SCRDIR=str(scratch),
                       OMP_NUM_THREADS="4", NCPUS="4", OMP_THREAD_LIMIT="4",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="4")

            def started(details):
                attempt.update(details)
                atomic_json(state_path, state)

            process = run_tree([args.gaussian / "g16.exe", gjf, directory / "job.log"],
                               directory, env, mask, 32, args.timeout_hours * 3600, started, affinity=False)
            attempt["process"] = process
            completion = log_completion(directory / "job.log")
            attempt["completion"] = completion
            if process["timed_out"]:
                attempt["status"] = state["status"] = "timeout"
                break
            if not completion["normal"] or process["exit_code"] != 0:
                attempt["status"] = "failed"
                atomic_json(state_path, state)
                if completion["scf_convergence_failure"] and not retry:
                    retry = True
                    continue
                state["status"] = "failed"
                break
            # G16W resolves a relative %chk against GAUSS_SCRDIR, whereas the
            # command's explicit log path is independent of that directory.
            chk = directory / "job.chk"
            if not chk.exists():
                chk = scratch / "job.chk"
            # The installed 16W utility takes two positional filenames. Unlike
            # the Unix utility, it interprets '-3' as a checkpoint filename.
            # Required output fields are checked after conversion below.
            conversion = run_tree([args.gaussian / "formchk.exe", chk, directory / "job.fch"],
                                  directory, env, mask, 32, 1800, affinity=False)
            attempt["formchk_process"] = conversion
            fchk_path = directory / "job.fch"
            if conversion["exit_code"] or conversion["timed_out"] or not fchk_path.exists():
                attempt["status"] = state["status"] = "failed"
                attempt["failure"] = "Formatted checkpoint conversion failed"
                break
            result = read_fchk(fchk_path)
            attempt["result_identity"] = result
            attempt["identity_gate"] = source_identity_check(case, result)
            attempt["status"] = state["status"] = (
                "collected" if attempt["identity_gate"]["input_identity_consistent"] else "identity_mismatch")
            attempt["fchk_sha256"] = sha256(fchk_path)
            state["result_fchk"] = str(fchk_path)
            state["result_log"] = str(directory / "job.log")
            break
        if state["status"] == "running":
            state["status"] = "failed"
            state["failure"] = "Predeclared SCF retries exhausted"
    except Exception as error:
        state["status"] = "failed"
        state["failure"] = f"{type(error).__name__}: {error}"
        if state["attempts"] and state["attempts"][-1]["status"] == "running":
            state["attempts"][-1]["status"] = "failed"
    state["finished_epoch"] = time.time()
    state["scientific_adjudication"] = "pending_whole_batch"
    atomic_json(state_path, state)
    atomic_json(args.output / "reproduction" / case["case_id"] / "result.json", state)
    return state


def run_reproduction(args):
    # Windows byte-range locking prevents two resumed coordinators from taking
    # the same resource allocation or executing the same pending case twice.
    import msvcrt
    args.jobs.mkdir(parents=True, exist_ok=True)
    lock = (args.jobs / "supervisor.lock").open("a+b")
    if lock.tell() == 0:
        lock.write(b"0")
        lock.flush()
    lock.seek(0)
    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    campaign = json.loads((args.output / "campaign.json").read_text(encoding="utf-8"))
    if len(campaign["cases"]) != 273:
        raise ValueError("Original batch registry count changed")
    cores = physical_core_masks(12)
    if len(cores) < 12:
        raise RuntimeError("The frozen resource layout needs 12 available physical cores")
    # Do not bind this parent: G16W inherits affinity and its legacy PGI startup
    # then faults. Each Gaussian process tree has a kernel CPU-rate hard cap.
    masks = [sum(cores[:4]), sum(cores[4:8])]
    runner_identity = digest({"runner": sha256(Path(__file__)),
                              "process": sha256(Path(__file__).with_name("validation_process.py")),
                              "case_set": campaign["case_set_identity"],
                              "gaussian_exe": sha256(args.gaussian / "g16.exe"),
                              "timeout_hours": args.timeout_hours})
    run_identity_path = args.jobs / "R001" / "runner-identity.json"
    if run_identity_path.exists():
        previous = json.loads(run_identity_path.read_text())
        if previous["identity"] != runner_identity:
            raise RuntimeError("Runner changed after batch freeze; finish/review the existing batch first")
    else:
        atomic_json(run_identity_path, {"identity": runner_identity,
                    "coordinator_pid": os.getpid(), "created_epoch": time.time(),
                    "gaussian_masks": masks, "other_work_mask": sum(cores[8:]),
                    "total_physical_cores": 12, "gaussian_tree_memory_total_gib": 64,
                    "reserved_other_memory_gib": 64,
                    "gaussian_cpu_limit": "4 runtime threads plus 4 logical CPU-equivalent job hard cap",
                    "other_work_limit": "4 cores includes the coordinator; reserve 1 for orchestration"})
        snapshots = run_identity_path.parent / "runner-source"
        snapshots.mkdir()
        for name in ("gaussian_rebuild.py", "validation_process.py"):
            shutil.copy2(Path(__file__).with_name(name), snapshots / name)
    pending = queue.Queue()
    for case in campaign["cases"]:
        pending.put(case)
    states = {}
    write_lock = threading.Lock()
    start = time.time()

    def write_progress():
        counts = Counter(value["status"] for value in states.values())
        counts["not_terminal"] = 273 - len(states)
        complete = len(states) == 273
        atomic_json(args.output / "progress.json", {
            "schema_version": 1, "case_set_identity": campaign["case_set_identity"],
            "runner_identity": runner_identity, "updated_epoch": time.time(),
            "stage": "reproduction_collected" if complete else "reproduction_running",
            "coordinator_pid": os.getpid(), "counts": dict(counts),
            "collection_barrier_reached": complete, "scientifically_passed": 0,
            "all_required_evidence_complete": False,
            "resource_allocation": {"gaussian_cores": 8, "gaussian_max_tree_memory_gib": 64,
                                    "other_work_cores": 4, "other_work_memory_gib": 64},
            "cases": {key: value["status"] for key, value in sorted(states.items())}})

    def worker(mask):
        while True:
            try:
                case = pending.get_nowait()
            except queue.Empty:
                return
            result = collect_case(case, args, mask, runner_identity)
            with write_lock:
                states[case["case_id"]] = result
                write_progress()
                print(json.dumps({"case": case["case_id"], "state": result["status"],
                                  "terminal": len(states), "total": 273}), flush=True)
            pending.task_done()

    write_progress()
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(worker, mask) for mask in masks]
        for future in futures:
            future.result()
    atomic_json(args.output / "reproduction" / "batch.json", {
        "case_count": len(states), "case_set_identity": campaign["case_set_identity"],
        "runner_identity": runner_identity, "all_terminal": len(states) == 273,
        "started_epoch": start, "finished_epoch": time.time(),
        "counts": dict(Counter(x["status"] for x in states.values())),
        "scientific_adjudication": "pending", "reference_recomputation": "pending"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("freeze", "run-reproduction", "status"))
    parser.add_argument("--collection", type=Path, default=DEFAULT_COLLECTION)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--jobs", type=Path, default=DEFAULT_JOBS)
    parser.add_argument("--gaussian", type=Path, default=DEFAULT_GAUSSIAN)
    parser.add_argument("--timeout-hours", type=float, default=24)
    args = parser.parse_args()
    if args.command == "freeze":
        freeze(args)
    elif args.command == "run-reproduction":
        run_reproduction(args)
    else:
        value = json.loads((args.output / "progress.json").read_text())
        value.pop("cases", None)
        print(json.dumps(value, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
