"""Resource-bounded producer for the focused NBO matrix cases.

Every attempt is retained in a new directory. This module never edits the
initial ten-case producer data or the Gaussian/NBO installation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import math
from datetime import datetime, timezone

from validation_process import atomic_json, physical_core_masks, run_tree


WORK = Path(__file__).resolve().parents[2]
BASE = WORK / "nbo-cases"
CASES = WORK / "nbo-focus-cases"
COMPUTE = WORK / "nbo-focus-compute"
NBO_BIN = Path(r"D:\CalChem\NBO\bin")
GAUSSIAN_BIN = Path(r"D:\CalChem\Gaussian\Gaussian 16 W")
KEYLIST = "PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 NAOMO=W51 AONAO=W52 NAONBO=W53 ARCHIVE PLOT"
MATRIX_NUMBERS = {"aonbo": 37, "nbomo": 49, "naomo": 51, "aonao": 52, "naonbo": 53}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def job_env(directory: Path, *, gennbo: bool) -> dict:
    env = dict(os.environ)
    env.update(GAUSS_EXEDIR=str(GAUSSIAN_BIN), GAUSS_SCRDIR=str(directory),
               OMP_NUM_THREADS="3", NCPUS="3", OMP_THREAD_LIMIT="3",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="3",
               NBOEXE=str(NBO_BIN / "nbo7.i8.exe"), NBOMEM="512mb")
    env["PATH"] = str(GAUSSIAN_BIN) + os.pathsep + env.get("PATH", "")
    return env


def run_once(args: list[str | Path], directory: Path, label: str,
             timeout: int = 1800, *, gennbo: bool = False,
             env_overrides: dict[str, str] | None = None) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    record = directory / f"{label}-process.json"
    if record.exists():
        raise RuntimeError(f"attempt_already_exists:{record}")
    env = job_env(directory, gennbo=gennbo)
    env.update(env_overrides or {})
    result = run_tree([str(item) for item in args], directory,
                      env,
                      sum(physical_core_masks(3)), 20, timeout,
                      affinity=gennbo,
                      on_started=lambda data: atomic_json(directory / f"{label}-live.json", data))
    atomic_json(record, result)
    return result


def reanalyse_original(number: str, attempt: str = "matrix-01") -> Path:
    source = BASE / f"NBO-{number}" / ("analysis-03" if number == "10" else "analysis-02") / "FILE.47"
    if not source.is_file():
        raise FileNotFoundError(source)
    directory = CASES / f"NBO-{number}" / attempt
    if directory.exists():
        raise RuntimeError(f"attempt_directory_already_exists:{directory}")
    directory.mkdir(parents=True)
    source_text = source.read_text(encoding="ascii")
    modified, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", source_text, count=1)
    if count != 1:
        raise RuntimeError("source_archive_missing_empty_NBO_keylist")
    inp = directory / "focus.47"
    inp.write_text(modified, encoding="ascii", newline="\n")
    manifest = {"case": f"NBO-{number}", "kind": "GenNBO_reanalysis_of_existing_Gaussian_archive",
                "source": str(source), "source_sha256": sha256(source),
                "input": str(inp), "input_sha256": sha256(inp),
                "keylist": KEYLIST, "started_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(directory / "manifest.json", manifest)
    result = run_once([NBO_BIN / "gennbo.i8.exe", inp.name], directory, "gennbo", 300,
                      gennbo=True)
    manifest["result"] = {k: result.get(k) for k in ("exit_code", "timed_out", "peak_tree_commit_bytes", "core_count", "cpu_rate_hard_cap")}
    manifest["files"] = {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)}
                         for p in directory.glob("FILE.*") if p.is_file()}
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    atomic_json(directory / "manifest.json", manifest)
    print(json.dumps({"directory": str(directory), "result": manifest["result"],
                      "matrix_files": sorted(manifest["files"])}, indent=2))
    return directory


def full_formchk_original(number: str, attempt: str = "formchk3-02") -> Path:
    # The installed Windows A.03 formchk did not support the tested format-3
    # invocations. Four retained failures are in NBO-07/formchk3-01..04;
    # the first invocation also overwrote its source checkpoint, since restored.
    # Keep the named entry point for a clear failure, without starting a job.
    raise RuntimeError("formchk3_disabled_on_this_Gaussian_A03_installation; see NBO-07/formchk3-01..04")


def _stable_paths(case_id: str, analysis: Path, fchk: Path, report: Path,
                  *, resource_process: Path, charge: int, multiplicity: int,
                  geometry: str, model: str, provenance: str) -> dict:
    paths = {name: str((analysis / f"FILE.{number}").resolve())
             for name, number in MATRIX_NUMBERS.items()}
    paths.update(fchk=str(fchk.resolve()), report=str(report.resolve()),
                 archive47=str((analysis / "FILE.47").resolve()),
                 sao=str((analysis / "FILE.50").resolve()))
    result = json.loads(resource_process.read_text(encoding="utf-8"))
    return {"case_id": case_id, "charge": charge, "multiplicity": multiplicity,
            "geometry_angstrom": geometry, "model": model,
            "provenance": provenance, **paths,
            "resource": {"process_record": str(resource_process.resolve()),
                         "nprocshared": None if "gennbo" in resource_process.name else 3,
                         "gaussian_stage_nprocshared": 3,
                         "gaussian_mem_gib": 12,
                         "tree_memory_limit_gib": result.get("memory_limit_gib"),
                         "effective_core_count": result.get("core_count"),
                         "effective_core": result.get("core_count"),
                         "cpu_rate_hard_cap": result.get("cpu_rate_hard_cap"),
                         "peak_tree_commit_bytes": result.get("peak_tree_commit_bytes"),
                         "exit_code": result.get("exit_code"),
                         "timed_out": result.get("timed_out")}}


def publish_existing(number: str) -> dict:
    case_id = f"NBO-{number}"
    directory = CASES / case_id / "matrix-01"
    old = json.loads((BASE / case_id / "manifest.json").read_text(encoding="utf-8"))
    fchk_relationship = "same_Gaussian_NBO_checkpoint"
    if number == "07":
        fchk = BASE / case_id / "canonical" / "canonical.fchk"
        fchk_relationship = ("preceding_canonical_SCF; independently_tested_same_wavefunction_"
                             "columnwise_up_to_phase_with_archive_LCAOMO")
    else:
        fchk = BASE / case_id / "analysis-02" / "paired-canonical.fchk"
    record = _stable_paths(case_id, directory, fchk, directory / "launcher.log",
                           resource_process=directory / "gennbo-process.json",
                           charge=old["charge"], multiplicity=old["multiplicity"],
                           geometry=old["atoms_cartesian_angstrom"],
                           model=f"{old['method']}/{old['basis']}",
                           provenance="GenNBO_reanalysis_from_existing_Gaussian_NBO_archive")
    record["fchk_relationship"] = fchk_relationship
    record["source_archive47"] = str((BASE / case_id / "analysis-02" / "FILE.47").resolve())
    record["source_archive47_sha256"] = sha256(Path(record["source_archive47"]))
    record["gaussian_resource_process"] = str((BASE / case_id / "analysis-02" / "analysis-process.json").resolve())
    record["fchk_sha256"] = sha256(fchk)
    record["analysis_input_sha256"] = sha256(directory / "focus.47")
    record["keylist"] = KEYLIST
    record["source_step"] = "initial_scf" if number == "07" else "analysis-02"
    record["nbo_step"] = "analysis-02"
    record["reanalysis_step"] = "matrix-01"
    if number == "07":
        record["status"] = "numerically_column_compatible_with_per_column_phase"
        record["incomplete_same_step_fchk"] = str((BASE / case_id / "analysis-02" / "paired-canonical.fchk").resolve())
        record["incomplete_same_step_fchk_sha256"] = sha256(Path(record["incomplete_same_step_fchk"]))
        record["negative_fixture"] = {"kind": "incomplete_beta_FCHK",
                                      "path": record["incomplete_same_step_fchk"],
                                      "sha256": record["incomplete_same_step_fchk_sha256"]}
        record["column_comparison"] = str((directory / "reference.json").resolve())
    else:
        record["status"] = "same_analysis_step_paired_FCHK"
    atomic_json(CASES / case_id / "case.json", record)
    return record


def _carbonyl_geometry() -> str:
    return "C 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.128000"


def _nickel_carbonyl_geometry() -> str:
    # Fixed ideal tetrahedron, Ni-C 1.82 Å and radial C-O 1.15 Å.
    directions = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    rows = ["Ni 0.000000 0.000000 0.000000"]
    for sign in directions:
        rows.append("C " + " ".join(f"{1.82 * item / math.sqrt(3):.6f}" for item in sign))
        rows.append("O " + " ".join(f"{2.97 * item / math.sqrt(3):.6f}" for item in sign))
    return "\n".join(rows)


NEW_CASES = {
    "CO": {"case_id": "NBO-F11", "title": "carbon monoxide fixed geometry",
           "charge": 0, "multiplicity": 1, "geometry": _carbonyl_geometry(),
           "basis": "6-31G(d)", "ecp": None,
           "geometry_note": "Fixed C-O 1.128 A; interface model, not optimized ground state"},
    "NICO4": {"case_id": "NBO-F12", "title": "nickel tetracarbonyl fixed tetrahedral geometry",
              "charge": 0, "multiplicity": 1, "geometry": _nickel_carbonyl_geometry(),
              "basis": "GenECP", "ecp": "Ni LANL2DZ; C/O 6-31G(d)",
              "geometry_note": "Fixed ideal Td: Ni-C 1.82 A, radial C-O 1.15 A; interface model, not optimized ground state"}}


def _run_stage(args: list[str | Path], directory: Path, label: str,
               expected: list[Path], normal_log: Path | None = None) -> dict:
    record = directory / f"{label}-process.json"
    if record.exists():
        result = json.loads(record.read_text(encoding="utf-8"))
    else:
        result = run_once(args, directory, label)
    success = (result["exit_code"] == 0 and not result["timed_out"]
               and all(item.is_file() and item.stat().st_size > 0 for item in expected))
    if normal_log and normal_log.is_file():
        success &= "Normal termination of Gaussian" in normal_log.read_text(errors="replace")
    if not success:
        raise RuntimeError(f"stage_failed_or_incomplete:{record}")
    return result


def _reanalyse_new_archive(source: Path, directory: Path, case_id: str) -> dict:
    if not source.is_file():
        raise FileNotFoundError(source)
    directory.mkdir(parents=True, exist_ok=True)
    source_text = source.read_text(encoding="ascii")
    modified, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", source_text, count=1)
    if count != 1:
        raise RuntimeError(f"source_archive_missing_empty_NBO_keylist:{source}")
    inp = directory / "focus.47"
    if inp.exists() and inp.read_text(encoding="ascii") != modified:
        raise RuntimeError(f"reanalysis_input_mismatch:{inp}")
    inp.write_text(modified, encoding="ascii", newline="\n")
    matrices = [directory / f"FILE.{item}" for item in (37, 49, 50, 51, 52, 53, 47)]
    result = _run_stage([NBO_BIN / "gennbo.i8.exe", inp.name], directory, "gennbo",
                        [directory / "launcher.log", *matrices])
    atomic_json(directory / "manifest.json", {
        "case_id": case_id, "kind": "GenNBO_same_archive_matrix_reanalysis",
        "source_archive47": str(source.resolve()), "source_sha256": sha256(source),
        "input": str(inp.resolve()), "input_sha256": sha256(inp),
        "keylist": KEYLIST, "exit_code": result["exit_code"],
        "matrix_sha256": {p.name: sha256(p) for p in matrices}})
    return result


def produce_new(name: str) -> dict:
    spec = NEW_CASES[name]
    case_id = spec["case_id"]
    target = CASES / case_id
    canonical = target / "canonical"
    analysis = target / "analysis-01"
    canonical.mkdir(parents=True, exist_ok=True)
    analysis.mkdir(parents=True, exist_ok=True)
    basis_tail = ""
    if spec["ecp"]:
        basis_tail = "C O 0\n6-31G(d)\n****\nNi 0\nLANL2DZ\n****\n\nNi 0\nLANL2DZ\n\n"
    inp = canonical / "canonical.gjf"
    canonical_input = ("%nprocshared=3\n%mem=12GB\n%chk=canonical.chk\n"
                       f"#p PBE1PBE/{spec['basis']} SCF=(XQC,Tight) Integral=UltraFine NoSymm\n\n"
                       f"{spec['title']}\n\n{spec['charge']} {spec['multiplicity']}\n"
                       f"{spec['geometry']}\n\n{basis_tail}")
    if inp.exists() and inp.read_text(encoding="ascii") != canonical_input:
        raise RuntimeError(f"canonical_input_mismatch:{inp}")
    inp.write_text(canonical_input, encoding="ascii", newline="\n")
    canonical_chk = canonical / "canonical.chk"
    _run_stage([GAUSSIAN_BIN / "g16.exe", inp, canonical / "canonical.log"],
               canonical, "canonical", [canonical_chk, canonical / "canonical.log"],
               canonical / "canonical.log")
    _run_stage([GAUSSIAN_BIN / "formchk.exe", canonical_chk, canonical / "canonical.fchk"],
               canonical, "canonical-formchk", [canonical / "canonical.fchk"])
    analysis_chk = analysis / "analysis.chk"
    if not analysis_chk.exists():
        shutil.copy2(canonical_chk, analysis_chk)
    nbo_input = ("%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
                 "#p PBE1PBE/ChkBasis Geom=AllCheck Guess=Read SCF=(XQC,Tight) "
                 "Integral=UltraFine NoSymm Pop=NBO6Read\n\n"
                 f"$NBO {KEYLIST} $END\n\n")
    nbo_inp = analysis / "analysis.gjf"
    if nbo_inp.exists() and nbo_inp.read_text(encoding="ascii") != nbo_input:
        raise RuntimeError(f"analysis_input_mismatch:{nbo_inp}")
    nbo_inp.write_text(nbo_input, encoding="ascii", newline="\n")
    # Gaussian A.03/NBO7 interface reports NAONBO on lfn30 despite W53 in
    # the keylist. A subsequent GenNBO pass writes the full W37/49/50/51/52/53
    # set from this exact archive, avoiding a mixture of analysis versions.
    matrices = [analysis / f"FILE.{item}" for item in (30, 37, 49, 50, 51, 52, 47)]
    _run_stage([GAUSSIAN_BIN / "g16.exe", nbo_inp, analysis / "analysis.log"],
               analysis, "analysis", [analysis / "analysis.log", analysis_chk, *matrices],
               analysis / "analysis.log")
    paired = analysis / "paired-canonical.fchk"
    _run_stage([GAUSSIAN_BIN / "formchk.exe", analysis_chk, paired],
               analysis, "paired-formchk", [paired])
    matrix_dir = target / "matrix-01"
    _reanalyse_new_archive(analysis / "FILE.47", matrix_dir, case_id)
    record = _stable_paths(case_id, matrix_dir, paired, matrix_dir / "launcher.log",
                           resource_process=matrix_dir / "gennbo-process.json",
                           charge=spec["charge"], multiplicity=spec["multiplicity"],
                           geometry=spec["geometry"], model=f"PBE1PBE/{spec['basis']}",
                           provenance="GenNBO_reanalysis_from_same_Gaussian_NBO_step_archive")
    record.update(title=spec["title"], geometry_note=spec["geometry_note"],
                  ecp=spec["ecp"], fchk_relationship="same_Gaussian_NBO_checkpoint",
                  status="same_analysis_step_paired_FCHK",
                  source_step="analysis-01", nbo_step="analysis-01",
                  reanalysis_step="matrix-01",
                  fchk_sha256=sha256(paired), archive47_sha256=sha256(matrix_dir / "FILE.47"),
                  canonical_checkpoint_sha256=sha256(canonical_chk),
                  analysis_checkpoint_sha256=sha256(analysis_chk),
                  keylist=KEYLIST, analysis_input_sha256=sha256(matrix_dir / "focus.47"),
                  source_archive47=str((analysis / "FILE.47").resolve()),
                  source_archive47_sha256=sha256(analysis / "FILE.47"),
                  gaussian_report=str((analysis / "analysis.log").resolve()),
                  gaussian_resource_process=str((analysis / "analysis-process.json").resolve()))
    atomic_json(target / "case.json", record)
    return record


def publish_index() -> list[dict]:
    rows = []
    for path in sorted(CASES.glob("*/case.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        for key in ("fchk", "report", "archive47", "aonbo", "nbomo",
                    "naomo", "aonao", "naonbo"):
            target = Path(row[key])
            if not target.is_absolute() or not target.is_file() or target.stat().st_size == 0:
                raise RuntimeError(f"index_required_file_missing_or_empty:{row['case_id']}:{key}:{target}")
        reference = Path(row["aonbo"]).parent / "reference.json"
        if reference.is_file():
            reference_data = json.loads(reference.read_text(encoding="utf-8"))
            row["reference"] = str(reference.resolve())
            row["reference_status"] = reference_data["status"]
            row["electron_bookkeeping"] = reference_data["electron_bookkeeping"]
            row["ecp_core_derived"] = reference_data["electron_bookkeeping"]["ecp_core_derived"]
        atomic_json(path, row)
        rows.append(row)
    atomic_json(CASES / "cases.json", rows)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["reanalyse", "formchk3", "publish-existing", "produce", "index"])
    parser.add_argument("case", nargs="?")
    parser.add_argument("--attempt", default="matrix-01")
    args = parser.parse_args()
    if args.command == "formchk3":
        parser.error("formchk3 is disabled on this Gaussian A.03 installation; see NBO-07/formchk3-01..04")
    if args.command == "publish-existing":
        if args.case not in ("01", "07"):
            parser.error("publish-existing case must be 01 or 07")
        print(json.dumps(publish_existing(args.case), indent=2))
        return 0
    if args.command == "produce":
        if args.case not in NEW_CASES:
            parser.error("produce case must be CO or NICO4")
        print(json.dumps(produce_new(args.case), indent=2))
        return 0
    if args.command == "index":
        print(json.dumps({"case_count": len(publish_index()),
                          "path": str(CASES / "cases.json")}))
        return 0
    if args.case not in ("01", "07"):
        parser.error("case must be 01 or 07")
    directory = reanalyse_original(args.case, args.attempt)
    result = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))["result"]
    return 0 if result["exit_code"] == 0 and not result["timed_out"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
