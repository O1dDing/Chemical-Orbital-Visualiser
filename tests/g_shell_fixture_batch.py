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

"""Auxiliary pure/Cartesian g regression evidence using an existing species.

H2+ is OLD-039. These representation/basis variants do not count as new external
molecules or as passing formal cases. Gaussian produces every coefficient;
the source wavefunction is never edited to fit COV. All Gaussian/reference
collection finishes before the baseline CUDA component comparison starts.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from gaussian_rebuild import DEFAULT_GAUSSIAN, DEFAULT_OUTPUT, digest, log_completion, read_fchk, sha256
from validation_process import atomic_json, run_tree

DEFAULT_ROOT = DEFAULT_OUTPUT / "regression-fixtures" / "g-shell-v2"
DEFAULT_WORK = Path(r"F:\Dev\cov-cycle-20260906\g-shell-regression-v2")
DEFAULT_EXE = Path(r"F:\Dev\cov-native-validation-20260905-build\Release\cov_cuda_cube_reference.exe")


def gaussian_env(directory, cores=2):
    return dict(os.environ, GAUSS_EXEDIR=str(DEFAULT_GAUSSIAN), GAUSS_SCRDIR=str(directory),
                GAUSS_MEMDEF="4GB", GAUSS_PDEF=str(cores), OMP_NUM_THREADS=str(cores),
                NCPUS=str(cores), OMP_THREAD_LIMIT=str(cores), OPENBLAS_NUM_THREADS="1",
                MKL_NUM_THREADS=str(cores))


def prepare(args):
    if (args.root / "manifest.json").exists():
        raise RuntimeError("Fixture set is already frozen")
    source = json.loads((DEFAULT_OUTPUT / "inputs/OLD-039/identity.json").read_text(encoding="utf-8"))
    f = source["fchk_identity"]
    sys.path.insert(0, r"F:\Dev\cov-cycle-20260906\dependencies")
    import basis_set_exchange as bse
    basis = bse.get_basis("cc-pV5Z", elements=f["Atomic numbers"], version=1, fmt="gaussian94", header=False)
    u = [1/math.sqrt(14), 2/math.sqrt(14), 3/math.sqrt(14)]
    a, b, c = u
    rotation = [[a, -b, -c], [b, 1-b*b/(1+a), -b*c/(1+a)], [c, -b*c/(1+a), 1-c*c/(1+a)]]
    xyz = f["Current cartesian coordinates"]
    moved = [[sum(rotation[i][j]*xyz[3*n+j] for j in range(3)) for i in range(3)]
             for n in range(f["Number of atoms"])]
    fixtures = []
    for label, flags, shell_type in (("pure", "5D 7F", -4), ("cartesian", "6D 10F", 4)):
        directory = args.root / label
        directory.mkdir(parents=True, exist_ok=False)
        lines = ["%chk=job.chk", "%mem=4GB", "%nprocshared=2",
                 f"#p UPBE1PBE/Gen SP {flags} NoSymm SCF=(VeryTight,XQC,MaxCycle=1024) Integral=UltraFine",
                 "", "Existing H2+ species: off-axis g-shell convention regression", "",
                 f'{f["Charge"]} {f["Multiplicity"]}']
        for z, position in zip(f["Atomic numbers"], moved):
            lines.append(str(z) + " " + " ".join(f"{v*0.529177210903:.14e}" for v in position))
        lines.extend(["", basis.rstrip(), "", ""])
        (directory / "input.gjf").write_text("\n".join(lines), encoding="ascii", newline="\n")
        fixture = {"name": label, "source_case_id": "OLD-039", "source_case_identity": source["case_identity"],
                   "source_fchk": source["original_fchk"], "source_fchk_sha256": source["original_fchk_sha256"],
                   "species": "H2+ doublet", "basis": "cc-pV5Z", "basis_data_version": 1,
                   "basis_library": "BSE " + bse.version(), "coordinate_rotation": rotation,
                   "input_coordinates_unit": "angstrom", "gaussian_shell_flags": flags,
                   "expected_g_shell_type": shell_type, "input_sha256": sha256(directory / "input.gjf"),
                   "counts_as_new_external_molecule": False,
                   "purpose": "AO/MO convention and truncation regression; no equilibrium or global-state claim"}
        fixture["identity"] = digest(fixture)
        fixtures.append(fixture)
    atomic_json(args.root / "manifest.json", {"schema_version": 1, "created_epoch": time.time(),
                "fixtures": fixtures, "set_identity": digest(fixtures),
                "scope": "Auxiliary static numerical fixtures; not formal COV UI/renderer acceptance",
                "grid": {"origin_bohr": [-6, -6, -6], "step_bohr": 0.375, "shape": [33, 33, 33]},
                "thresholds": {"nrms_max": 1e-4, "abs_cosine_min": 1-1e-7, "relative_peak_error_max": 1e-3}})
    print("Prepared two g-shell representation fixtures from OLD-039", flush=True)


def compute_inputs(args, manifest):
    results = []
    for fixture in manifest["fixtures"]:
        public = args.root / fixture["name"]
        result_path = public / "calculation.json"
        if result_path.exists():
            result = json.loads(result_path.read_text())
            if result["fixture_identity"] != fixture["identity"]:
                raise ValueError("Fixture identity changed")
            results.append(result)
            continue
        directory = args.work / fixture["name"]
        directory.mkdir(parents=True, exist_ok=False)
        if sha256(public / "input.gjf") != fixture["input_sha256"]:
            raise ValueError("Frozen fixture input changed")
        shutil.copy2(public / "input.gjf", directory / "job.gjf")
        env = gaussian_env(directory)
        result = {"fixture_identity": fixture["identity"], "name": fixture["name"], "status": "failed"}
        result["gaussian"] = run_tree([DEFAULT_GAUSSIAN / "g16.exe", directory / "job.gjf", directory / "job.log"],
                                       directory, env, 5, 12, 3600, affinity=False)
        result["termination"] = log_completion(directory / "job.log")
        if result["gaussian"]["exit_code"] == 0 and result["termination"]["normal"]:
            result["formchk"] = run_tree([DEFAULT_GAUSSIAN / "formchk.exe", directory / "job.chk", directory / "job.fch"],
                                          directory, env, 5, 12, 300, affinity=False)
            if result["formchk"]["exit_code"] == 0:
                fields = read_fchk(directory / "job.fch")
                if fixture["expected_g_shell_type"] not in fields["Shell types"]:
                    raise ValueError("Requested g-shell representation was not produced")
                result["fchk_identity"] = fields
                result["fchk_sha256"] = sha256(directory / "job.fch")
                result["status"] = "collected"
                shutil.copy2(directory / "job.fch", public / "wavefunction.fch")
        shutil.copy2(directory / "job.log", public / "gaussian.log")
        atomic_json(result_path, result)
        results.append(result)
        print(json.dumps({"fixture": fixture["name"], "status": result["status"]}), flush=True)
    atomic_json(args.root / "input-collection.json", {"all_terminal": len(results) == 2,
                "results": [{"name": r["name"], "status": r["status"]} for r in results]})
    if any(result["status"] != "collected" for result in results):
        raise RuntimeError("Fixture calculations ended with errors; preserve and review before reference collection")


def grid_spec(manifest):
    grid = manifest["grid"]
    # Gaussian's free-format reader distinguishes integer and real tokens.
    # Coordinates and vector components must be real even when exactly zero.
    lines = ["-1 " + " ".join(f"{float(v):.17e}" for v in grid["origin_bohr"])]
    for axis in range(3):
        values = [grid["step_bohr"] if k == axis else 0 for k in range(3)]
        lines.append(str((-1 if axis == 0 else 1)*grid["shape"][axis]) + " " + " ".join(f"{float(v):.17e}" for v in values))
    return "\n".join(lines) + "\n"


def reference_worker(args):
    """Each producer owns a bounded tree; the caller reserves args.workers cores."""
    manifest = json.loads((args.root / "manifest.json").read_text())
    producer_hash = sha256(DEFAULT_GAUSSIAN / 'cubegen.exe')
    if manifest.get('cubegen_sha256', producer_hash) != producer_hash:
        raise ValueError('Frozen Gaussian cube producer identity changed')
    jobs = []
    for fixture in manifest["fixtures"]:
        directory = args.root / fixture["name"]
        calculation = json.loads((directory / "calculation.json").read_text())
        if calculation["status"] != "collected" or sha256(directory / "wavefunction.fch") != calculation["fchk_sha256"]:
            raise ValueError("Fixture collection identity not intact")
        f = calculation["fchk_identity"]
        count = f["Number of independent functions"]
        for spin, token, offset in (("alpha", "AMO", 0), ("beta", "BMO", count)):
            if spin == "beta" and not f.get("has_beta_coefficients"):
                continue
            for orbital in range(count):
                jobs.append({"fixture": fixture["name"], "spin": spin, "number": orbital+1,
                             "kind": f"{token}={orbital+1}", "internal_index": offset+orbital,
                             "input_sha256": calculation['fchk_sha256']})
    specification = grid_spec(manifest)

    def collect(job):
        root = args.root / job["fixture"]
        directory = root / "cubes" / f'{job["spin"]}-{job["number"]:04d}'
        directory.mkdir(parents=True, exist_ok=True)
        state_path = directory / "result.json"
        cube = directory / "orbital.cube"
        identity = digest({'job': job, 'grid_specification': specification, 'cubegen_sha256': producer_hash,
                           'runner_files': manifest.get('runner_files'),
                           'reference_round_identity': manifest.get('reference_round_identity')})
        if state_path.exists():
            result = json.loads(state_path.read_text())
            if result.get('job_identity') == identity:
                if result['status'] != 'collected' or (cube.is_file() and sha256(cube) == result['cube_sha256']):
                    return result  # Retain terminal failures; never silently retry them.
            rejected = dict(job, status='identity_mismatch', job_identity=identity,
                            reason='Existing record identity or cube content differs from frozen input')
            atomic_json(directory/'reuse-rejection.json', rejected)
            return rejected
        result = dict(job, status="failed", job_identity=identity)
        env = gaussian_env(directory, cores=1)
        env["GAUSS_MEMDEF"] = "2GB"
        start = time.monotonic()
        command = [str(DEFAULT_GAUSSIAN / "cubegen.exe"), "1", job["kind"], str(root / "wavefunction.fch"), str(cube), "-1", "h"]
        stdin_path = directory / "grid-specification.txt"
        stdin_path.write_text(specification, encoding="ascii", newline="\n")
        try:
            run = run_tree(command, directory, env, 1, 4, 180, affinity=False, stdin_path=stdin_path)
            log = directory / "launcher.log"
            result.update(process=run, exit_code=run['exit_code'],
                          console_log=str(log), console_sha256=sha256(log),
                          stdout=log.read_text(encoding="utf-8", errors="replace")[-4000:],
                          grid_specification_sha256=sha256(stdin_path))
            if run['exit_code'] == 0 and not run['timed_out'] and cube.exists():
                result.update(status="collected", cube_sha256=sha256(cube), cube=str(cube))
        except Exception as error:
            result["failure"] = f"{type(error).__name__}: {error}"
        result["wall_seconds"] = time.monotonic()-start
        atomic_json(state_path, result)
        return result

    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(collect, job) for job in jobs]
        for future in as_completed(futures):
            records.append(future.result())
    records.sort(key=lambda r: (r["fixture"], r["internal_index"]))
    atomic_json(args.root / "reference-grid-collection.json", {"all_terminal": len(records) == len(jobs),
                "expected_orbitals": len(jobs), "records": records,
                "source": "Installed Gaussian cubegen, every alpha/beta MO on the frozen grid"})
    print(json.dumps({"reference_orbitals": len(records), "failed": sum(r["status"] != "collected" for r in records)}))


def compare_baseline(args):
    collection = json.loads((args.root / "reference-grid-collection.json").read_text())
    if not collection["all_terminal"]:
        raise RuntimeError("Reference collection barrier not reached")
    if any(record["status"] != "collected" for record in collection["records"]):
        raise RuntimeError("Reference failures require review; do not generate a partial pass")
    manifest = json.loads((args.root/'manifest.json').read_text())
    if manifest.get('baseline_executable_sha256', sha256(args.exe)) != sha256(args.exe):
        raise RuntimeError('Frozen baseline executable changed')
    for record in collection['records']:
        if sha256(Path(record['cube'])) != record['cube_sha256']:
            raise RuntimeError('Collected cube identity changed')
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.pop("COV_REFERENCE_DIAGNOSTIC_ODD_M", None)
    result = {"executable": str(args.exe), "executable_sha256": sha256(args.exe),
              "scope": "Unmodified production CUDA evaluator component with full texture readback; no UI/frame coverage claim",
              "counts_as_formal_case_pass": False, "comparisons": []}
    for fixture in ("pure", "cartesian"):
        records = [r for r in collection["records"] if r["fixture"] == fixture]
        # Keep each Windows argv under its command-line length bound.
        for begin in range(0, len(records), 64):
            command = [str(args.exe), str(args.root / fixture / "wavefunction.fch")]
            for row in records[begin:begin+64]:
                command.extend([str(row["internal_index"]), row["cube"]])
            run = subprocess.run(command, capture_output=True, text=True, cwd=args.work, env=env, timeout=300)
            rows = [json.loads(line) for line in run.stdout.splitlines() if line.startswith('{"internal_index"')]
            result["comparisons"].append({"fixture": fixture, "offset": begin, "return_code": run.returncode,
                                          "expected_count": len(records[begin:begin+64]), "rows": rows,
                                          "stderr": run.stderr[-2000:]})
    atomic_json(args.root / "baseline-cuda-comparison.json", result)
    rows = [(group["fixture"], row) for group in result["comparisons"] for row in group["rows"]]
    print(json.dumps({"compared": len(rows), "failed": sum(not r["pass"] for _, r in rows),
                      "max_nrms": max(r["nrms"] for _, r in rows)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "collect", "reference-worker", "compare-worker"))
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--workers", type=int, choices=(1, 2, 3), default=1)
    args = parser.parse_args()
    if args.command != "prepare":
        frozen = json.loads((args.root / "manifest.json").read_text())
        for name, expected in frozen.get('runner_files', {}).items():
            if sha256(Path(__file__).parent/name) != expected:
                raise RuntimeError('Frozen reference runner identity changed: '+name)
        if 'reference_resource_policy' in frozen and args.workers != frozen['reference_resource_policy']['workers']:
            raise RuntimeError('Worker count differs from frozen resource allocation')
    if args.command == "prepare":
        prepare(args)
    elif args.command == "reference-worker":
        reference_worker(args)
    elif args.command == "compare-worker":
        compare_baseline(args)
    else:
        manifest = json.loads((args.root / "manifest.json").read_text())
        compute_inputs(args, manifest)
        # Producers enforce their own per-job tree limits. Wrapping them in a
        # second CPU-rate Job would multiply the nested quotas. The comparator
        # instead receives one bounded tree for its entire process group.
        for operation in ("reference-worker", "compare-worker"):
            logdir = args.work / operation
            logdir.mkdir(parents=True, exist_ok=True)
            command = [sys.executable, Path(__file__), operation, "--root", args.root,
                       "--work", args.work, "--exe", args.exe, "--workers", str(args.workers)]
            if operation == 'reference-worker':
                started = time.time()
                with (logdir/'launcher.log').open('ab') as log:
                    try:
                        run = subprocess.run(command, cwd=logdir, stdout=log, stderr=log,
                                             timeout=manifest.get('reference_resource_policy', {}).get('batch_timeout_seconds', 90600))
                        returncode, timed_out = run.returncode, False
                    except subprocess.TimeoutExpired:
                        # Closing the worker closes its non-inherited Job
                        # handles and terminates every unfinished producer.
                        returncode, timed_out = 124, True
                result = {'argv': list(map(str, command)), 'started_epoch': started,
                          'finished_epoch': time.time(), 'exit_code': returncode, 'timed_out': timed_out,
                          'producer_tree_limits': {'workers': args.workers, 'cores_each': 1,
                                                   'memory_gib_each': 4, 'timeout_seconds_each': 180}}
            else:
                result = run_tree(command, logdir, dict(os.environ), 1, 12, 3600, affinity=False)
            atomic_json(args.root / (operation + ".json"), result)
            if result["exit_code"]:
                raise RuntimeError(operation + " failed; inspect its preserved log")


if __name__ == "__main__":
    main()
