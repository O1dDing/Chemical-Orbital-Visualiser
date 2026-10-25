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

"""Collect every primary reference candidate before scientific adjudication.

Uses the same exclusive resource reservation as original-library reproduction.
Every calculation chain is resumable by immutable input and runner identity.
Normal termination is collection evidence, never the campaign's pass verdict.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import queue
import re
import shutil
import threading
import time

from gaussian_rebuild import DEFAULT_GAUSSIAN, DEFAULT_JOBS, DEFAULT_OUTPUT, digest, log_completion, read_fchk, sha256
from gaussian_reference_inputs import make_followup_input
from validation_process import atomic_json, physical_core_masks, run_tree


def stage_diagnostics(log):
    text = log.read_text(encoding="utf-8", errors="replace")
    frequencies = [float(value) for line in re.findall(r"Frequencies --\s*([^\r\n]+)", text)
                   for value in line.split()]
    stable = "wavefunction is stable under" in text.lower()
    unstable = bool(re.search(r"wavefunction has (?:an?|one or more).*instabilit", text, re.I))
    return {"optimization_completed": "Optimization completed." in text,
            "frequencies_cm1": frequencies,
            "reported_negative_frequencies_cm1": [f for f in frequencies if f < 0],
            "stability_reported": stable, "instability_encountered": unstable,
            "s2_records": re.findall(r"S\*\*2 before annihilation[^\r\n]*", text)[-3:],
            "interpretation": "Raw job diagnostics; physical verdict deferred until all candidates terminate"}


def identity_gate(reference, original, fields):
    problems = []
    for key in ("Atomic numbers", "Charge", "Multiplicity"):
        if fields[key] != original["fchk_identity"][key]:
            problems.append(key)
    for field, reference_key in (("Number of electrons", "explicit_electrons"),
                                 ("Number of alpha electrons", "alpha_electrons"),
                                 ("Number of beta electrons", "beta_electrons")):
        if fields[field] != reference[reference_key]:
            problems.append(field)
    ecp = reference["ecp_core_electrons_by_atomic_number"]
    expected = [z - ecp[str(z)] for z in fields["Atomic numbers"]]
    if fields["Nuclear charges"] != expected:
        problems.append("ECP nuclear charges")
    if any(abs(l) == 4 and l != reference["expected_g_shell_type"] for l in fields["Shell types"]):
        problems.append("g-shell representation")
    return {"consistent": not problems, "mismatches": problems}


def collect_stage(reference, original, args, case_base, name, mask, runner_identity,
                  checkpoint=None, method=None):
    base = case_base / name
    base.mkdir(exist_ok=True)
    state_path = base / "stage.json"
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if state["runner_identity"] != runner_identity or state["reference_identity"] != reference["identity"]:
            raise RuntimeError("Stage identity changed; never reuse a different build/input")
        if state["status"] in ("collected", "failed", "timeout", "identity_mismatch"):
            if state["status"] == "collected" and sha256(state["fchk"]) != state["fchk_sha256"]:
                raise ValueError("Previously collected checkpoint changed")
            return state
        if state["attempts"] and state["attempts"][-1]["status"] == "running":
            state["attempts"][-1]["status"] = "interrupted"
    else:
        state = {"name": name, "status": "pending", "runner_identity": runner_identity,
                 "reference_identity": reference["identity"], "attempts": []}
    finished_attempts = [item for item in state["attempts"] if item["status"] != "interrupted"]
    for attempt_index in range(len(finished_attempts), 2):
        directory = base / f'attempt-{len(state["attempts"])+1:02d}'
        directory.mkdir(exist_ok=False)
        (directory / "scratch").mkdir()
        retry = attempt_index > 0
        if name == "initial":
            source = args.references / "reference-inputs" / reference["case_id"] / "initial.gjf"
            if sha256(source) != reference["initial_input_sha256"]:
                raise ValueError("Frozen initial reference input changed")
            text = source.read_text(encoding="ascii")
            if retry:
                text = text.replace("MaxCycle=1024", "MaxCycle=2048")
        else:
            if not checkpoint or not method:
                raise ValueError("Reference-chain checkpoint/method absent")
            shutil.copy2(checkpoint, directory / "job.chk")
            stage = "stability" if name.startswith("stability") else "relax_after_stability"
            text = make_followup_input(reference, method, stage, retry=retry)
        (directory / "job.gjf").write_text(text, encoding="ascii", newline="\n")
        attempt = {"number": len(state["attempts"])+1, "status": "running", "directory": str(directory),
                   "input_sha256": sha256(directory / "job.gjf"), "predeclared_scf_retry": retry,
                   "starting_checkpoint": str(checkpoint) if checkpoint else None,
                   "starting_checkpoint_sha256": sha256(checkpoint) if checkpoint else None}
        state["attempts"].append(attempt)
        state["status"] = "running"
        atomic_json(state_path, state)
        env = dict(os.environ, GAUSS_EXEDIR=str(args.gaussian), GAUSS_SCRDIR=str(directory / "scratch"),
                   OMP_NUM_THREADS="4", NCPUS="4", OMP_THREAD_LIMIT="4", MKL_NUM_THREADS="4",
                   OPENBLAS_NUM_THREADS="1")

        def started(details):
            attempt.update(details)
            atomic_json(state_path, state)

        attempt["process"] = run_tree([args.gaussian / "g16.exe", directory / "job.gjf", directory / "job.log"],
                directory, env, mask, 32, args.timeout_hours*3600, started, affinity=False)
        completion = log_completion(directory / "job.log")
        attempt["completion"] = completion
        if attempt["process"]["timed_out"]:
            attempt["status"] = state["status"] = "timeout"
            break
        if attempt["process"]["exit_code"] or not completion["normal"]:
            attempt["status"] = "failed"
            if completion["scf_convergence_failure"] and not retry:
                atomic_json(state_path, state)
                continue
            state["status"] = "failed"
            break
        checkpoint_out = directory / "job.chk"
        if not checkpoint_out.exists():
            checkpoint_out = directory / "scratch" / "job.chk"
        attempt["formchk"] = run_tree([args.gaussian / "formchk.exe", checkpoint_out, directory / "job.fch"],
                                      directory, env, mask, 32, 1800, affinity=False)
        if attempt["formchk"]["exit_code"] or attempt["formchk"]["timed_out"]:
            attempt["status"] = state["status"] = "failed"
            attempt["failure"] = "FCHK conversion failed"
            break
        fields = read_fchk(directory / "job.fch")
        attempt["identity_gate"] = identity_gate(reference, original, fields)
        attempt["diagnostics"] = stage_diagnostics(directory / "job.log")
        attempt["status"] = state["status"] = "collected" if attempt["identity_gate"]["consistent"] else "identity_mismatch"
        state.update(checkpoint=str(checkpoint_out), fchk=str(directory / "job.fch"),
                     log=str(directory / "job.log"), fields=fields,
                     fchk_sha256=sha256(directory / "job.fch"), diagnostics=attempt["diagnostics"])
        break
    if state["status"] == "running":
        state["status"] = "failed"
    atomic_json(state_path, state)
    return state


def needs_relaxation(before, stable):
    """Predeclared job-chain rule, not a physical acceptance decision."""
    old, new = before["fields"], stable["fields"]
    return (stable["diagnostics"]["instability_encountered"] or
            old["method"] != new["method"] or
            abs(old["Total Energy"]-new["Total Energy"]) > 1e-8)


def collect_reference(reference, original, args, mask, runner_identity):
    base = args.jobs / args.round_name / reference["case_id"]
    base.mkdir(parents=True, exist_ok=True)
    result_path = base / "result.json"
    if result_path.exists():
        previous = json.loads(result_path.read_text())
        if previous["runner_identity"] != runner_identity:
            raise ValueError("Reference runner identity changed")
        if previous["status"] != "running":
            return previous
    result = {"case_id": reference["case_id"], "runner_identity": runner_identity,
              "reference_identity": reference["identity"], "status": "running", "stages": [],
              "physical_adjudication": "pending_whole_batch", "scientific_pass": False}
    atomic_json(result_path, result)
    try:
        current = collect_stage(reference, original, args, base, "initial", mask, runner_identity)
        result["stages"].append(str(base / "initial" / "stage.json"))
        if current["status"] != "collected":
            raise RuntimeError("Initial reference stage: " + current["status"])
        for cycle in range(3):
            name = f"stability-{cycle:02d}"
            stable = collect_stage(reference, original, args, base, name, mask, runner_identity,
                                   checkpoint=current["checkpoint"], method=current["fields"]["method"])
            result["stages"].append(str(base / name / "stage.json"))
            if stable["status"] != "collected":
                raise RuntimeError("Stability stage: " + stable["status"])
            if needs_relaxation(current, stable):
                if cycle == 2:
                    raise RuntimeError("Predeclared stability repair cycles exhausted")
                name = f"relax-after-stability-{cycle+1:02d}"
                current = collect_stage(reference, original, args, base, name, mask, runner_identity,
                                        checkpoint=stable["checkpoint"], method=stable["fields"]["method"])
                result["stages"].append(str(base / name / "stage.json"))
                if current["status"] != "collected":
                    raise RuntimeError("Reoptimization/refrequency stage: " + current["status"])
                continue
            if not stable["diagnostics"]["stability_reported"]:
                raise RuntimeError("Normal Gaussian termination without explicit stability evidence")
            result.update(status="candidate_collected", reference_fchk=current["fchk"],
                          reference_log=current["log"], reference_fchk_sha256=current["fchk_sha256"],
                          stability_fchk=stable["fchk"], stability_log=stable["log"],
                          stability_repair_cycles=cycle,
                          competing_states_pending=reference["competing_states_required"])
            break
    except Exception as error:
        result["status"] = "failed"
        result["failure"] = f"{type(error).__name__}: {error}"
    result["finished_epoch"] = time.time()
    atomic_json(result_path, result)
    atomic_json(args.references / "results" / reference["case_id"] / "result.json", result)
    return result


def run_batch(args):
    import msvcrt
    if not re.fullmatch(r"[A-Z][A-Z0-9-]{1,32}", args.round_name):
        raise ValueError("Invalid reference round name")
    manifest = json.loads((args.references / "reference-candidates.json").read_text(encoding="utf-8"))
    campaign = json.loads((args.output / "campaign.json").read_text(encoding="utf-8"))
    expected = [f"OLD-{i:03d}" for i in range(1, 274)]
    if [c["case_id"] for c in manifest["candidates"]] != expected:
        raise ValueError("Reference batch must preserve all 273 original cases")
    if manifest["original_case_set_identity"] != campaign["case_set_identity"]:
        raise ValueError("Reference/original case-set mismatch")
    args.jobs.mkdir(parents=True, exist_ok=True)
    lock = (args.jobs / "supervisor.lock").open("a+b")
    if lock.tell() == 0:
        lock.write(b"0")
        lock.flush()
    start = time.time()
    while True:
        lock.seek(0)
        try:
            msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError:
            if not args.wait_for_reproduction:
                raise RuntimeError("Another campaign owns the Gaussian resource reservation")
        else:
            barrier = args.output / "reproduction" / "batch.json"
            if barrier.exists():
                previous = json.loads(barrier.read_text())
                if (previous.get("all_terminal") and previous.get("case_count") == 273 and
                        previous.get("case_set_identity") == campaign["case_set_identity"]):
                    break
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            if not args.wait_for_reproduction:
                raise RuntimeError("Original reproduction has not reached its collection barrier")
        atomic_json(args.references / "progress.json", {"status": "waiting_for_original_reproduction",
                    "coordinator_pid": os.getpid(), "updated_epoch": time.time(),
                    "gaussian_jobs_running": 0, "scientific_passes": 0})
        time.sleep(30)
    cores = physical_core_masks(12)
    if len(cores) < 12:
        raise ValueError("Frozen resource policy expects 12 physical cores")
    masks = [sum(cores[:4]), sum(cores[4:8])]
    source_files = ("gaussian_reference_batch.py", "gaussian_reference_inputs.py",
                    "gaussian_rebuild.py", "validation_process.py")
    identity = digest({"source_files": {name: sha256(Path(__file__).with_name(name)) for name in source_files},
                       "reference_set": manifest["reference_set_identity"],
                       "gaussian_exe": sha256(args.gaussian / "g16.exe"),
                       "formchk_exe": sha256(args.gaussian / "formchk.exe"),
                       "timeout_hours": args.timeout_hours, "round_name": args.round_name})
    identity_path = args.jobs / args.round_name / "runner-identity.json"
    if identity_path.exists():
        if json.loads(identity_path.read_text())["identity"] != identity:
            raise ValueError("Reference runner changed; resume its frozen snapshot or create a new round")
    else:
        atomic_json(identity_path, {"identity": identity, "created_epoch": time.time(),
                    "reference_set_identity": manifest["reference_set_identity"],
                    "source_files": {name: sha256(Path(__file__).with_name(name)) for name in source_files},
                    "gaussian_concurrency": 2, "cores_per_job": 4,
                    "memory_limit_per_tree_gib": 32, "other_work_including_coordinator_cores": 4})
        snapshot = identity_path.parent / "runner-source"
        snapshot.mkdir()
        for name in source_files:
            shutil.copy2(Path(__file__).with_name(name), snapshot / name)
    pending = queue.Queue()
    originals = {case["case_id"]: case for case in campaign["cases"]}
    for reference in manifest["candidates"]:
        pending.put(reference)
    results = {}
    write_lock = threading.Lock()

    def progress():
        complete = len(results) == 273
        counts = dict(Counter(value["status"] for value in results.values()))
        counts["not_terminal"] = 273 - len(results)
        value = {"schema_version": 1, "runner_identity": identity,
                 "reference_set_identity": manifest["reference_set_identity"],
                 "updated_epoch": time.time(), "coordinator_pid": os.getpid(),
                 "stage": "primary_references_collected" if complete else "primary_references_running",
                 "counts": counts, "collection_barrier_reached": complete,
                 "scientifically_passed": 0, "all_required_evidence_complete": False,
                 "cases": {key: result["status"] for key, result in sorted(results.items())},
                 "remaining_requirements": ["post-batch input and physical adjudication",
                    "applicable competing electronic states and environments",
                    "COV general fixes and complete native/independent validation"]}
        atomic_json(args.references / "progress.json", value)
        atomic_json(args.output / "progress.json", value)

    def worker(mask):
        while True:
            try:
                reference = pending.get_nowait()
            except queue.Empty:
                return
            result = collect_reference(reference, originals[reference["case_id"]], args, mask, identity)
            with write_lock:
                results[reference["case_id"]] = result
                progress()
                print(json.dumps({"case": reference["case_id"], "state": result["status"],
                                  "terminal": len(results), "total": 273}), flush=True)
            pending.task_done()

    progress()
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(worker, mask) for mask in masks]
        for future in futures:
            future.result()
    atomic_json(args.references / "batch.json", {"all_terminal": len(results) == 273,
                "case_count": len(results), "runner_identity": identity,
                "reference_set_identity": manifest["reference_set_identity"],
                "counts": dict(Counter(result["status"] for result in results.values())),
                "started_including_queue_epoch": start, "finished_epoch": time.time(),
                "scientific_adjudication": "pending", "competing_states": "pending_when_applicable"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--references", type=Path, default=DEFAULT_OUTPUT / "references" / "REF-001")
    parser.add_argument("--jobs", type=Path, default=DEFAULT_JOBS)
    parser.add_argument("--gaussian", type=Path, default=DEFAULT_GAUSSIAN)
    parser.add_argument("--round-name", default="REF-001")
    parser.add_argument("--wait-for-reproduction", action="store_true")
    parser.add_argument("--timeout-hours", type=float, default=24)
    run_batch(parser.parse_args())
