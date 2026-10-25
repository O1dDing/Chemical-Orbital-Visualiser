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

"""Windows Gaussian 16W A.03 / NBO7 one-submission template.

Default: read-only dry run. --run creates a NEW output directory and launches
bounded producers using the existing tests/validation_process.py supervisor.
Requires an existing Gaussian/NBO installation and the COV source tree.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import sys

MATRICES = {
    "AONBO": 37, "NBOMO": 49, "SAO": 50, "NAOMO": 51,
    "AONAO": 52, "NAONBO": 53, "AONHO": 54, "AONLMO": 55,
    "AOPNAO": 56, "NAONHO": 57, "NAONLMO": 58, "NHONBO": 59,
    "NBONLMO": 60, "NLMOMO": 61, "AOMO": 62,
}
KEYLIST = "PRINT=3 E2PERT=0.0 " + " ".join(
    f"{key}=W{number}" for key, number in MATRICES.items()
) + " ARCHIVE PLOT"
HEADINGS = {
    key: ("SAO" if key == "SAO" else
          f"{key[len(prefix):]}s in the {prefix} basis:")
    for key in MATRICES
    for prefix in [next((p for p in ("NLMO", "NAO", "NHO", "NBO", "AO")
                       if key.startswith(p)), "")]
}


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def make_input(recipe):
    required = ("title", "method", "basis", "charge", "multiplicity", "geometry_angstrom")
    for key in required:
        if key not in recipe:
            raise ValueError(f"Missing recipe field: {key}")
    if not isinstance(recipe["charge"], int) or not isinstance(recipe["multiplicity"], int):
        raise ValueError("Charge and multiplicity must be integers")
    if recipe["multiplicity"] < 1:
        raise ValueError("Multiplicity must be positive")
    for key in ("title", "method", "basis", "memory", "integral", "scf"):
        value = str(recipe.get(key, ""))
        if any(c in value for c in "\r\n%#$"):
            raise ValueError(f"Field must be one literal line: {key}")
    geometry = recipe["geometry_angstrom"].strip()
    if not geometry or re.search(r"(?im)^\s*(%|#|\$|--link1--)", geometry):
        raise ValueError("Supply Cartesian atom rows only")
    threads = int(recipe.get("nprocshared", 3))
    if not 1 <= threads <= 3:
        raise ValueError("This template permits 1 to 3 Gaussian threads")
    symmetry = recipe.get("symmetry", "")
    if symmetry not in ("", "NoSymm"):
        raise ValueError("Template symmetry field must be empty (Gaussian default) or NoSymm")
    route = (f"#p {recipe['method']}/{recipe['basis']} "
             f"SCF={recipe.get('scf', '(XQC,Tight)')} "
             f"Integral={recipe.get('integral', 'UltraFine')} "
             + (symmetry + " " if symmetry else "") + "Pop=NBO6Read")
    tail = recipe.get("basis_ecp_tail", "").strip()
    return (f"%nprocshared={threads}\n%mem={recipe.get('memory', '2GB')}\n"
            f"%chk=canonical.chk\n{route}\n\n{recipe['title']}\n\n"
            f"{recipe['charge']} {recipe['multiplicity']}\n{geometry}\n\n"
            + (tail + "\n\n" if tail else "")
            + f"$NBO {KEYLIST} $END\n\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gaussian-bin", type=Path, required=True, help="Gaussian 16W A.03 installation directory")
    parser.add_argument("--nbo-bin", type=Path, required=True, help="NBO7 i8 executable directory")
    parser.add_argument("--runner-dir", type=Path,
                        default=Path(__file__).resolve().parents[2] / "tests",
                        help="COV source-tree tests directory containing validation_process.py")
    parser.add_argument("--memory-limit-gib", type=int, default=20)
    parser.add_argument("--gaussian-timeout", type=int, default=900)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    recipe = json.loads(args.recipe.read_text(encoding="utf-8"))
    gaussian_input = make_input(recipe)
    target = args.output.resolve()
    if target.exists():
        raise FileExistsError(f"Fresh output directory required: {target}")
    executables = {"gaussian": args.gaussian_bin / "g16.exe",
                   "formchk": args.gaussian_bin / "formchk.exe",
                   "gennbo": args.nbo_bin / "gennbo.i8.exe",
                   "nbo": args.nbo_bin / "nbo7.i8.exe",
                   "interface": args.nbo_bin / "g16nbo.i8.exe"}
    bridge = args.gaussian_bin / "gaunbo6.bat"
    runner = args.runner_dir / "validation_process.py"
    for path in [*executables.values(), bridge, runner]:
        if not path.is_file():
            raise FileNotFoundError(path)
    bridge_text = bridge.read_text(encoding="ascii")
    if str(args.nbo_bin).replace("/", "\\").lower() not in bridge_text.lower():
        raise ValueError("Installed A.03 bridge does not name the selected NBO bin directory")
    if not 1 <= args.memory_limit_gib <= 20 or not 1 <= args.gaussian_timeout <= 1800:
        raise ValueError("Template process bounds: 1-20 GiB, Gaussian timeout 1-1800 seconds")
    preview = {"mode": "run" if args.run else "dry-run; no processes or output files",
               "target": str(target), "expected_gaussian_revision": "A.03",
               "recipe": recipe, "gaussian_input": gaussian_input,
               "executables": {key: str(value) for key, value in executables.items()},
               "installed_bridge": str(bridge), "supervisor": str(runner),
               "matrix_numbers": MATRICES, "pipeline": ["Gaussian + NBO", "formchk", "GenNBO", "package"]}
    if not args.run:
        print(json.dumps(preview, ensure_ascii=False, indent=2))
        return 0
    spec = importlib.util.spec_from_file_location("cov_one_job_validation_process", runner)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    masks = module.physical_core_masks(3)
    target.mkdir(parents=True, exist_ok=False)
    gaussian_dir, converted_dir, matrices_dir, package = [target / x for x in ("gaussian", "formchk", "matrices", "cov-package")]
    for folder in (gaussian_dir, converted_dir, matrices_dir):
        folder.mkdir()
    manifest = {**preview, "status": "started", "stages": {},
                "recipe_sha256": digest(args.recipe), "wrapper_sha256": digest(__file__),
                "bridge_sha256": digest(bridge), "supervisor_sha256": digest(runner),
                "executable_sha256": {key: digest(value) for key, value in executables.items()}}
    save(target / "job.json", manifest)
    env = dict(os.environ)
    env.update(GAUSS_EXEDIR=str(args.gaussian_bin),
               NBOEXE=str(executables["nbo"]), NBOMEM="512mb",
               OMP_NUM_THREADS="3", NCPUS="3", OMP_THREAD_LIMIT="3",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="3")
    env["PATH"] = str(args.gaussian_bin) + os.pathsep + env.get("PATH", "")

    def stage(name, argv, folder, timeout, affinity):
        local_env = dict(env, GAUSS_SCRDIR=str(folder))
        result = module.run_tree([str(x) for x in argv], folder, local_env,
                                 sum(masks), args.memory_limit_gib, timeout,
                                 affinity=affinity,
                                 on_started=lambda value: save(folder / "live.json", value))
        save(folder / "process.json", result)
        manifest["stages"][name] = result
        save(target / "job.json", manifest)
        if result["exit_code"] != 0 or result["timed_out"]:
            raise RuntimeError(f"Stage failed; retained outputs: {name}")

    try:
        inp = gaussian_dir / "job.gjf"
        inp.write_text(gaussian_input, encoding="ascii", newline="\n")
        log = gaussian_dir / "job.log"
        stage("gaussian", [executables["gaussian"], inp, log], gaussian_dir,
              args.gaussian_timeout, False)
        log_text = log.read_text(errors="replace")
        if "Normal termination of Gaussian" not in log_text or not re.search(r"(?:G16Rev\s*|Revision\s+)A\.03\b", log_text):
            raise RuntimeError("Expected Gaussian A.03 normal termination absent")
        chk = gaussian_dir / "canonical.chk"
        archive = gaussian_dir / "FILE.47"
        for path in (chk, archive):
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f"Missing Gaussian result: {path}")
        canonical_hash = digest(chk)
        archive_hash = digest(archive)
        fchk = converted_dir / "canonical.fchk"
        stage("formchk", [executables["formchk"], chk, fchk], converted_dir, 300, False)
        if not fchk.is_file() or not fchk.stat().st_size:
            raise RuntimeError("formchk output absent")
        original = archive.read_text(encoding="ascii")
        for section in ("GENNBO", "COORD", "BASIS", "CONTRACT", "OVERLAP", "DENSITY", "FOCK", "LCAOMO"):
            if not re.search(r"\$" + section + r"\b", original, flags=re.I):
                raise RuntimeError(f"Archive lacks ${section}; not synthesized")
        updated, count = re.subn(r"\$NBO\b.*?\$END", f"$NBO {KEYLIST} $END", original, count=1, flags=re.S | re.I)
        if count != 1 or len(re.findall(r"\$NBO\b", original, flags=re.I)) != 1:
            raise RuntimeError("Unique archive $NBO keylist required")
        matrix_input = matrices_dir / "aomo.47"
        matrix_input.write_text(updated, encoding="ascii", newline="\n")
        stage("gennbo", [executables["gennbo"], matrix_input.name], matrices_dir, 600, True)
        report = matrices_dir / "launcher.log"
        if not report.is_file() or not re.search(r"NBO\s+7\.0", report.read_text(errors="replace")):
            raise RuntimeError("GenNBO NBO7 report absent")
        for kind, number in MATRICES.items():
            path = matrices_dir / f"FILE.{number}"
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f"Requested matrix missing: {kind} {path}")
            if kind != "SAO" and HEADINGS[kind] not in path.read_text(errors="replace"):
                raise RuntimeError(f"Requested matrix heading mismatch: {kind} {path}")
        produced_archive = matrices_dir / "FILE.47"
        if not produced_archive.is_file() or not produced_archive.stat().st_size:
            raise RuntimeError("GenNBO ARCHIVE output absent")
        if digest(chk) != canonical_hash or digest(archive) != archive_hash:
            raise RuntimeError("Original Gaussian checkpoint/archive changed")
        package.mkdir()
        shutil.copy2(fchk, package / "canonical.fchk")
        shutil.copy2(report, package / "analysis.nbo")
        shutil.copy2(produced_archive, package / "FILE.47")
        for number in MATRICES.values():
            shutil.copy2(matrices_dir / f"FILE.{number}", package / f"FILE.{number}")
        members = ["canonical.fchk", "analysis.nbo", "FILE.47"] + [f"FILE.{n}" for n in MATRICES.values()]
        (package / "drop.covnbopkg").write_text("COV_NBO_PACKAGE 1\n" + "\n".join(members) + "\n", encoding="utf-8")
        manifest.update(status="producer_complete; COV import pending",
                        gaussian_checkpoint_sha256=canonical_hash,
                        preserved_gaussian_archive_sha256=archive_hash,
                        source_step="single Gaussian SP with NBO6Read",
                        reanalysis_step="GenNBO from that archived wavefunction; no new SCF",
                        gaussian_version_header=[x.strip() for x in log_text.splitlines() if "Revision" in x or "G16Rev" in x][:3],
                        nbo_version_header=[x.strip() for x in report.read_text(errors="replace").splitlines() if "NBO 7.0" in x or "NBO Version" in x or "Cite this program" in x][:5],
                        package_files={p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in package.iterdir()},
                        matrix_roles={key: f"FILE.{number}" for key, number in MATRICES.items()},
                        limits="Producer checks only; no COV import, spin completeness or matrix numeric acceptance implied")
        save(target / "job.json", manifest)
        # This metadata quotes producer banners. Keep it outside the input
        # folder so content-based discovery cannot treat it as another report.
        save(target / "cov-input-manifest.json", manifest)
        print(json.dumps({"status": manifest["status"], "package": str(package)}, ensure_ascii=False))
    except Exception as error:
        manifest.update(status="failed; outputs retained", error=str(error))
        save(target / "job.json", manifest)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
