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

"""Collect independent old/reproduction comparisons after the 273-case barrier.

The collection is immutable and exhaustive. Review flags request scientific
adjudication; neither a normal process exit nor absence of flags is a COV pass.
This batch never changes Gaussian inputs, COV algorithms, or acceptance limits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import sys
import time

from gaussian_rebuild import DEFAULT_OUTPUT, digest, sha256
from validation_process import atomic_json, run_tree

DEFAULT_ROOT = DEFAULT_OUTPUT / "reproduction-comparison" / "CMP-001"
DEFAULT_WORK = Path(r"F:\Dev\cov-cycle-20260906\reproduction-comparison\CMP-001")
DEPS = Path(r"F:\Dev\cov-validation-20260905\reference-deps")
SOURCES = ("reproduction_compare_batch.py", "wavefunction_reference.py", "gaussian_rebuild.py", "validation_process.py")
POLICY = {"all_273_terminal_before_start": True, "scientific_verdict": "manual_post_batch_adjudication",
          "cores_per_worker": 1, "workers": 2, "memory_gib_per_tree": 12,
          "timeout_hours_per_case": 4, "energy_window_hartree": 1e-5,
          "diagnostic_flags_not_acceptance_limits": {
              "geometry_displacement_bohr": 2e-6, "energy_delta_hartree": 1e-6,
              "mo_orthonormality_max_error": 5e-6, "occupied_min_principal_overlap": 1-1e-6,
              "energy_block_min_principal_overlap": 1-1e-5, "density_relative_hs_norm": 1e-4,
              "electron_trace_error": 5e-6}}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def freeze(args):
    if (args.root / "manifest.json").exists():
        raise ValueError("Comparison batch already frozen")
    campaign = read_json(args.output / "campaign.json")
    if len(campaign["cases"]) != 273:
        raise ValueError("Original registry must contain 273 fixed cases")
    source_hashes = {name: sha256(Path(__file__).with_name(name)) for name in SOURCES}
    # Include all independent reader/integral implementation source, alongside
    # package versions saved per result. Nothing is read from a COV calculation.
    library_hashes = {str(path.relative_to(DEPS)).replace("\\", "/"): sha256(path)
                      for package in ("iodata", "gbasis") for path in sorted((DEPS / package).rglob("*.py"))}
    manifest = {"schema_version": 1, "created_epoch": time.time(),
                "case_set_identity": campaign["case_set_identity"], "policy": POLICY,
                "source_hashes": source_hashes, "library_hashes": library_hashes,
                "cases": [{"case_id": case["case_id"], "case_identity": case["case_identity"]}
                          for case in campaign["cases"]]}
    manifest["identity"] = digest(manifest)
    snapshot = args.root / "runner-source"
    snapshot.mkdir(parents=True)
    for name in SOURCES:
        shutil.copy2(Path(__file__).with_name(name), snapshot / name)
    atomic_json(args.root / "manifest.json", manifest)
    print(json.dumps({"frozen_comparison": str(args.root), "identity": manifest["identity"]}))


def require_barrier(args, manifest):
    path = args.output / "reproduction" / "batch.json"
    if not path.exists():
        return False
    value = read_json(path)
    if value["case_set_identity"] != manifest["case_set_identity"]:
        raise ValueError("Reproduction and comparison case sets differ")
    return value["all_terminal"] and value["case_count"] == 273


def verify_runner(manifest):
    for name, expected in manifest["source_hashes"].items():
        if sha256(Path(__file__).with_name(name)) != expected:
            raise ValueError("Run from the frozen comparison source snapshot")
    for name, expected in manifest["library_hashes"].items():
        if sha256(DEPS / name) != expected:
            raise ValueError("Frozen independent reference library changed: " + name)


def review_flags(evidence, policy):
    limits = policy["diagnostic_flags_not_acceptance_limits"]
    flags = []
    if evidence["alignment"]["max_displacement_bohr"] > limits["geometry_displacement_bohr"]:
        flags.append("geometry_changed_beyond_rigid_alignment")
    energy = evidence["energy_delta_hartree"]
    if energy is None:
        flags.append("missing_total_energy")
    elif abs(energy) > limits["energy_delta_hartree"]:
        flags.append("fresh_scf_energy_difference")
    for spin, value in evidence["spins"].items():
        for version in ("old", "new"):
            if value[version+"_metric"]["max_orthonormality_error"] > limits["mo_orthonormality_max_error"]:
                flags.append(spin+"_"+version+"_mo_metric_residual")
        occupied = value["occupied_space"]
        if occupied["old_rank"] != occupied["new_rank"] or (
                occupied["status"] == "measured" and occupied["min_singular_value"] < limits["occupied_min_principal_overlap"]):
            flags.append(spin+"_occupied_space_difference")
        if value["unmatched_old_count"] or value["unmatched_new_count"]:
            flags.append(spin+"_represented_space_dimension_difference")
        if any(block["principal_overlap"]["min_singular_value"] < limits["energy_block_min_principal_overlap"]
               for block in value["energy_blocks"]):
            flags.append(spin+"_canonical_energy_block_difference")
    for name in ("mo_total", "mo_spin"):
        value = evidence["densities"][name]
        if value["relative_norm"] is not None and value["relative_norm"] > limits["density_relative_hs_norm"]:
            flags.append(name+"_density_operator_difference")
    return flags


def compare_one(args):
    manifest = read_json(args.root / "manifest.json")
    verify_runner(manifest)
    if not require_barrier(args, manifest):
        raise RuntimeError("Whole-original-batch terminal barrier has not been reached")
    case = read_json(args.output / "inputs" / args.case_id / "identity.json")
    original_entry = next(item for item in manifest["cases"] if item["case_id"] == args.case_id)
    if case["case_identity"] != original_entry["case_identity"]:
        raise ValueError("Frozen original identity changed")
    reproduction = read_json(args.output / "reproduction" / args.case_id / "result.json")
    result = {"case_id": args.case_id, "comparison_identity": manifest["identity"],
              "case_identity": case["case_identity"], "counts_as_cov_pass": False,
              "chemical_adjudication": "pending_whole_comparison_batch",
              "reproduction_collection_status": reproduction["status"], "review_flags": []}
    destination = args.root / "cases" / args.case_id
    destination.mkdir(parents=True, exist_ok=True)
    if reproduction["status"] != "collected":
        result.update(status="unavailable", reason="reproduction_not_collected_with_matching_identity")
    else:
        old_path, new_path = Path(case["original_fchk"]), Path(reproduction["result_fchk"])
        expected_new = reproduction["attempts"][-1]["fchk_sha256"]
        if sha256(old_path) != case["original_fchk_sha256"] or sha256(new_path) != expected_new:
            raise ValueError("Input evidence changed since collection")
        from wavefunction_reference import compare_wavefunctions
        result.update(original_fchk=str(old_path), original_fchk_sha256=case["original_fchk_sha256"],
                      reproduction_fchk=str(new_path), reproduction_fchk_sha256=expected_new)
        result["evidence"] = compare_wavefunctions(old_path, new_path, destination / "matrices.npz",
                                                   manifest["policy"]["energy_window_hartree"])
        result["review_flags"] = review_flags(result["evidence"], manifest["policy"])
        expected_electrons = case["fchk_identity"]["Number of electrons"]
        for version in ("old", "new"):
            trace = result["evidence"]["densities"]["mo_total"][version+"_electron_trace"]
            if abs(trace-expected_electrons) > manifest["policy"]["diagnostic_flags_not_acceptance_limits"]["electron_trace_error"]:
                result["review_flags"].append(version+"_electron_trace_residual")
        result["matrix_evidence_sha256"] = sha256(destination / "matrices.npz")
        result["status"] = "collected"
    atomic_json(destination / "comparison.json", result)


def collect_case(args, entry, identity):
    path = args.root / "cases" / entry["case_id"] / "result.json"
    if path.exists():
        previous = read_json(path)
        if previous["comparison_identity"] != identity:
            raise ValueError("Comparison result belongs to another frozen version")
        return previous
    work = args.work / entry["case_id"]
    work.mkdir(parents=True, exist_ok=True)
    result = dict(entry, comparison_identity=identity, status="failed", counts_as_cov_pass=False)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", PYTHONIOENCODING="utf-8")
    try:
        process = run_tree([sys.executable, Path(__file__), "worker", "--case-id", entry["case_id"],
                            "--output", args.output, "--root", args.root, "--work", args.work],
                           work, env, 1, POLICY["memory_gib_per_tree"],
                           POLICY["timeout_hours_per_case"]*3600, affinity=False)
        result["process"] = process
        comparison = path.with_name("comparison.json")
        if process["exit_code"] == 0 and comparison.exists():
            payload = read_json(comparison)
            if payload["comparison_identity"] != identity:
                raise ValueError("Worker returned evidence from another version")
            result.update(status=payload["status"], review_flags=payload["review_flags"],
                          evidence=str(comparison), evidence_sha256=sha256(comparison))
        else:
            result.update(status="timeout" if process["timed_out"] else "failed", failure="Inspect preserved launcher.log")
    except Exception as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    atomic_json(path, result)
    return result


def run_batch(args):
    import msvcrt
    manifest = read_json(args.root / "manifest.json")
    verify_runner(manifest)
    args.work.mkdir(parents=True, exist_ok=True)
    with (args.work / "supervisor.lock").open("a+b") as lock:
        if lock.tell() == 0:
            lock.write(b"0"); lock.flush()
        lock.seek(0)
        msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        while not require_barrier(args, manifest):
            if not args.wait_for_reproduction:
                raise RuntimeError("All 273 reproductions must terminate before numerical comparison")
            atomic_json(args.root / "progress.json", {"status": "waiting_for_all_original_reproductions",
                        "coordinator_pid": os.getpid(), "workers_running": 0, "scientific_passes": 0,
                        "updated_epoch": time.time()})
            time.sleep(30)
        while args.wait_for_fixture_release and not (args.wait_for_fixture_release / "compare-worker.json").exists():
            atomic_json(args.root / "progress.json", {"status": "waiting_for_auxiliary_allocation",
                        "coordinator_pid": os.getpid(), "workers_running": 0, "scientific_passes": 0,
                        "updated_epoch": time.time(), "original_reproduction_barrier_reached": True})
            time.sleep(15)
        records = []
        start = time.time()
        with ThreadPoolExecutor(max_workers=POLICY["workers"]) as pool:
            futures = [pool.submit(collect_case, args, entry, manifest["identity"]) for entry in manifest["cases"]]
            for future in as_completed(futures):
                records.append(future.result())
                atomic_json(args.root / "progress.json", {"status": "collecting_comparison_evidence",
                            "comparison_identity": manifest["identity"], "updated_epoch": time.time(),
                            "coordinator_pid": os.getpid(), "terminal_cases": len(records), "expected_cases": 273,
                            "counts": dict(Counter(row["status"] for row in records)), "scientific_passes": 0})
                print(json.dumps({"terminal": len(records), "case": records[-1]["case_id"], "status": records[-1]["status"]}), flush=True)
        records.sort(key=lambda row: row["case_id"])
        atomic_json(args.root / "batch.json", {"all_terminal": len(records) == 273, "case_count": len(records),
                    "comparison_identity": manifest["identity"], "case_set_identity": manifest["case_set_identity"],
                    "started_epoch": start, "finished_epoch": time.time(), "records": records,
                    "counts": dict(Counter(row["status"] for row in records)),
                    "chemical_adjudication": "pending", "counts_as_cov_pass": False})
        atomic_json(args.root / "progress.json", {"status": "comparison_collection_finished", "terminal_cases": len(records),
                    "scientific_passes": 0, "updated_epoch": time.time(), "chemical_adjudication": "pending"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("freeze", "run", "worker"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--case-id")
    parser.add_argument("--wait-for-reproduction", action="store_true")
    parser.add_argument("--wait-for-fixture-release", type=Path)
    args = parser.parse_args()
    if args.command == "freeze":
        freeze(args)
    elif args.command == "worker":
        compare_one(args)
    else:
        run_batch(args)


if __name__ == "__main__":
    main()
