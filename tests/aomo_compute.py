"""Produce complete NBO transformation matrices from a preserved FILE.47.

Each run is a new GenNBO attempt; no Gaussian or formchk process is launched.
The source archive, checkpoint, and earlier matrices are never modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
from datetime import datetime, timezone

from validation_process import atomic_json, physical_core_masks, run_tree


WORK = Path(__file__).resolve().parents[2]
NBO_BIN = Path(r"D:\CalChem\NBO\bin")
SOURCE = WORK / "nbo-cases"
DEST = WORK / "aomo-cases"
COMPUTE = WORK / "aomo-compute"
GAUSSIAN_BIN = Path(r"D:\CalChem\Gaussian\Gaussian 16 W")
NUMBER = {
    "AONBO": 37, "NBOMO": 49, "SAO": 50, "NAOMO": 51,
    "AONAO": 52, "NAONBO": 53, "AONHO": 54, "AONLMO": 55,
    "AOPNAO": 56, "NAONHO": 57, "NAONLMO": 58, "NHONBO": 59,
    "NBONLMO": 60, "NLMOMO": 61, "AOMO": 62,
}
KEYLIST = "PRINT=3 E2PERT=0.0 " + " ".join(
    f"{key}=W{number}" for key, number in NUMBER.items()
) + " ARCHIVE PLOT"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_archive(case: str, attempt: str) -> dict:
    if not re.fullmatch(r"\d{2}", case) or int(case) not in range(1, 11):
        raise ValueError("case_must_be_original_NBO_01_through_10")
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,50}", attempt):
        raise ValueError("attempt_must_be_short_lowercase_hyphenated_name")
    analysis = "analysis-03" if case == "10" else "analysis-02"
    source = SOURCE / f"NBO-{case}" / analysis / "FILE.47"
    if not source.is_file():
        raise FileNotFoundError(source)
    directory = DEST / f"NBO-{case}" / attempt
    if directory.exists():
        raise FileExistsError(f"new_attempt_required:{directory}")
    original = source.read_text(encoding="ascii")
    updated, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
    if count != 1:
        raise RuntimeError("source_archive_does_not_have_expected_empty_NBO_keylist")
    directory.mkdir(parents=True)
    input_path = directory / "aomo.47"
    input_path.write_text(updated, encoding="ascii", newline="\n")
    manifest = {
        "case_id": f"NBO-{case}", "attempt": attempt,
        "kind": "GenNBO_reanalysis_from_preserved_same_Gaussian_NBO_archive",
        "source_archive47": str(source.resolve()), "source_sha256": sha256(source),
        "input_archive47": str(input_path.resolve()), "input_sha256": sha256(input_path),
        "keylist": KEYLIST, "matrix_numbers": NUMBER,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(directory / "manifest.json", manifest)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="3", NCPUS="3", OMP_THREAD_LIMIT="3",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="3",
               NBOEXE=str(NBO_BIN / "nbo7.i8.exe"), NBOMEM="512mb",
               GAUSS_SCRDIR=str(directory))
    result = run_tree([str(NBO_BIN / "gennbo.i8.exe"), input_path.name],
                      directory, env, sum(physical_core_masks(3)), 20, 300,
                      affinity=True,
                      on_started=lambda data: atomic_json(directory / "gennbo-live.json", data))
    atomic_json(directory / "gennbo-process.json", result)
    expected = [directory / f"FILE.{number}" for number in NUMBER.values()]
    missing = [str(path) for path in expected if not path.is_file() or path.stat().st_size == 0]
    report = directory / "launcher.log"
    success = result["exit_code"] == 0 and not result["timed_out"] and not missing
    manifest.update(
        result={key: result.get(key) for key in ("exit_code", "timed_out", "core_count",
             "cpu_rate_hard_cap", "memory_limit_gib", "peak_tree_commit_bytes")},
        report=str(report.resolve()) if report.is_file() else None,
        missing_matrix_files=missing,
        files={path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
               for path in sorted(directory.glob("FILE.*")) if path.is_file()},
        status="producer_complete_unvalidated" if success else "producer_failed_preserved",
        finished_utc=datetime.now(timezone.utc).isoformat(),
    )
    atomic_json(directory / "manifest.json", manifest)
    return {"directory": str(directory.resolve()), "status": manifest["status"],
            "result": manifest["result"], "missing": missing,
            "matrix_files": sorted(manifest["files"])}


def reanalyse_seeded_archive(case_id: str, source_attempt: str, attempt: str) -> dict:
    if not re.fullmatch(r"OLD-\d{3}", case_id):
        raise ValueError("case_must_be_OLD_nnn")
    if not all(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,50}", word)
               for word in (source_attempt, attempt)):
        raise ValueError("invalid_attempt_name")
    prior = DEST / case_id / source_attempt
    source = prior / "analysis/FILE.47"
    if not source.is_file():
        raise FileNotFoundError(source)
    old_manifest = json.loads((prior / "manifest.json").read_text())
    if not (old_manifest["stages"].get("gaussian_nbo", {}).get("pass") or
            old_manifest["stages"].get("guessonly_nbo", {}).get("pass")):
        raise RuntimeError("prior_Gaussian_NBO_stage_not_complete")
    directory = DEST / case_id / attempt
    if directory.exists():
        raise FileExistsError(f"new_attempt_required:{directory}")
    original = source.read_text(encoding="ascii")
    updated, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
    if count != 1:
        raise RuntimeError("source_archive_has_unexpected_NBO_keylist")
    directory.mkdir(parents=True)
    inp = directory / "aomo.47"
    inp.write_text(updated, encoding="ascii", newline="\n")
    manifest = {"case_id": case_id, "kind": "GenNBO_reanalysis_of_seeded_Gaussian_archive",
                "source_archive47": str(source.resolve()), "source_sha256": sha256(source),
                "input_archive47": str(inp.resolve()), "input_sha256": sha256(inp),
                "frozen_fchk": old_manifest["frozen_fchk"],
                "frozen_fchk_sha256": old_manifest.get("frozen_fchk_sha256", old_manifest.get("frozen_sha256")),
                "canonical_fchk": old_manifest.get("canonical_fchk"),
                "fchk_relationship": ("canonical_checkpoint_snapshot_precedes_GuessOnly_analysis; "
                                      "requires_independent_same_wavefunction_identity_audit"
                                      if old_manifest["stages"].get("guessonly_nbo", {}).get("pass")
                                      else "requires_independent_column_or_subspace_compatibility_audit"),
                "keylist": KEYLIST, "matrix_numbers": NUMBER}
    atomic_json(directory / "manifest.json", manifest)
    result = _bounded_stage([str(NBO_BIN / "gennbo.i8.exe"), inp.name],
                            directory, "gennbo", gaussian=False, timeout=600)
    files = [directory / f"FILE.{number}" for number in NUMBER.values()]
    success = _stage_success(result, files)
    manifest.update(status="producer_complete_unvalidated" if success else "producer_failed_preserved",
                    result=str((directory / "gennbo-process.json").resolve()),
                    report=str((directory / "launcher.log").resolve()),
                    matrices={key: str((directory / f"FILE.{number}").resolve())
                              for key, number in NUMBER.items() if (directory / f"FILE.{number}").is_file()},
                    finished_utc=datetime.now(timezone.utc).isoformat())
    atomic_json(directory / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(directory.resolve()), "status": manifest["status"]}


def _bounded_stage(argv: list[str], directory: Path, label: str,
                   *, gaussian: bool, timeout: int) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    record = directory / f"{label}-process.json"
    if record.exists():
        raise RuntimeError(f"stage_already_attempted:{record}")
    env = dict(os.environ)
    env.update(GAUSS_EXEDIR=str(GAUSSIAN_BIN), GAUSS_SCRDIR=str(directory),
               OMP_NUM_THREADS="3", NCPUS="3", OMP_THREAD_LIMIT="3",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="3",
               NBOEXE=str(NBO_BIN / "nbo7.i8.exe"), NBOMEM="512mb")
    env["PATH"] = str(GAUSSIAN_BIN) + os.pathsep + env.get("PATH", "")
    result = run_tree(argv, directory, env, sum(physical_core_masks(3)), 20, timeout,
                      affinity=not gaussian,
                      on_started=lambda data: atomic_json(directory / f"{label}-live.json", data))
    atomic_json(record, result)
    return result


def _stage_success(result: dict, files: list[Path], *, normal_log: Path | None = None) -> bool:
    passed = result["exit_code"] == 0 and not result["timed_out"]
    passed = passed and all(path.is_file() and path.stat().st_size > 0 for path in files)
    if normal_log and normal_log.is_file():
        passed = passed and "Normal termination of Gaussian" in normal_log.read_text(errors="replace")
    return passed


def recover_fchk(case_id: str, attempt: str, *, native_symmetry: bool = False) -> dict:
    """Recreate a *new* Gaussian/NBO paired calculation from one frozen FCHK.

    The historical FCHK is source evidence, never a same-step pair for the new
    NBO matrices. Every stage and every unsuccessful attempt is retained.
    """
    ledger = json.loads((WORK / "aomo-coverage/coverage-ledger.json").read_text())
    row = next((item for item in ledger["cases"] if item["case_id"] == case_id), None)
    if row is None or not row["manifest_hash_match"]:
        raise ValueError(f"frozen_source_not_verified:{case_id}")
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,50}", attempt):
        raise ValueError("attempt_must_be_short_lowercase_hyphenated_name")
    source = Path(row["source_fchk"])
    source_sha = sha256(source)
    if source_sha != row["manifest_sha256"]:
        raise RuntimeError(f"frozen_fchk_hash_changed:{case_id}")
    target = DEST / case_id / attempt
    if target.exists():
        raise FileExistsError(f"new_attempt_required:{target}")
    target.mkdir(parents=True)
    imported = target / "imported"
    analysis = target / "analysis"
    transformed = target / "matrices"
    imported.mkdir()
    analysis.mkdir()
    transformed.mkdir()
    manifest = {"case_id": case_id, "attempt": attempt,
                "kind": "new_paired_Gaussian_NBO_analysis_seeded_by_frozen_FCHK",
                "frozen_fchk": str(source.resolve()), "frozen_fchk_sha256": source_sha,
                "frozen_manifest_sha256": row["manifest_sha256"],
                "method": row["method"], "basis": row["basis"],
                "charge": row["charge"], "multiplicity": row["multiplicity"],
                "fchk_relationship": "frozen_source_seed_only; not_same_analysis_step",
                "native_symmetry_retry": native_symmetry,
                "started_utc": datetime.now(timezone.utc).isoformat(), "stages": {}}
    atomic_json(target / "manifest.json", manifest)
    try:
        imported_chk = imported / "imported.chk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "unfchk.exe"), str(source),
                                 str(imported_chk)], imported, "unfchk", gaussian=True, timeout=120)
        manifest["stages"]["unfchk"] = {"result": str((imported / "unfchk-process.json").resolve()),
                                          "pass": _stage_success(result, [imported_chk])}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["unfchk"]["pass"]:
            raise RuntimeError("unfchk_failed")
        analysis_chk = analysis / "analysis.chk"
        shutil.copy2(imported_chk, analysis_chk)
        method = row["method"]
        if not method or not re.fullmatch(r"[RU]?[-A-Za-z0-9]+", method):
            raise RuntimeError(f"unsupported_method_token:{method}")
        inp = analysis / "analysis.gjf"
        gaussian_keylist = "PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT"
        symmetry_option = "" if native_symmetry else " NoSymm"
        inp.write_text("%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
                       f"#p {method}/ChkBasis Geom=AllCheck Guess=Read "
                       f"SCF=(XQC,Tight) Integral=UltraFine{symmetry_option} Pop=NBO6Read\n\n"
                       f"$NBO {gaussian_keylist} $END\n\n", encoding="ascii", newline="\n")
        log = analysis / "analysis.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(inp), str(log)],
                                analysis, "gaussian-nbo", gaussian=True, timeout=3600)
        archive = analysis / "FILE.47"
        manifest["stages"]["gaussian_nbo"] = {
            "result": str((analysis / "gaussian-nbo-process.json").resolve()),
            "pass": _stage_success(result, [analysis_chk, archive, log], normal_log=log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["gaussian_nbo"]["pass"]:
            raise RuntimeError("gaussian_nbo_failed")
        paired = analysis / "paired.fchk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "formchk.exe"),
                                 str(analysis_chk), str(paired)],
                                analysis, "formchk", gaussian=True, timeout=300)
        manifest["stages"]["formchk"] = {"result": str((analysis / "formchk-process.json").resolve()),
                                          "pass": _stage_success(result, [paired])}
        if manifest["stages"]["formchk"]["pass"]:
            fchk_text = paired.read_text(errors="replace")
            alpha_present = "Alpha MO coefficients" in fchk_text
            beta_present = "Beta MO coefficients" in fchk_text
            manifest["paired_fchk_spin_blocks"] = {"alpha": alpha_present, "beta": beta_present,
                                                   "beta_required": row["multiplicity"] > 1}
            if not alpha_present or (row["multiplicity"] > 1 and not beta_present):
                manifest["stages"]["formchk"]["pass"] = False
                manifest["stages"]["formchk"]["reason"] = "missing_required_MO_coefficient_spin_block"
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["formchk"]["pass"]:
            raise RuntimeError("formchk_failed")
        original = archive.read_text(encoding="ascii")
        modified, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
        if count != 1:
            raise RuntimeError("archive_has_unexpected_NBO_keylist")
        matrix_input = transformed / "aomo.47"
        matrix_input.write_text(modified, encoding="ascii", newline="\n")
        result = _bounded_stage([str(NBO_BIN / "gennbo.i8.exe"), matrix_input.name],
                                transformed, "gennbo", gaussian=False, timeout=600)
        files = [transformed / f"FILE.{number}" for number in NUMBER.values()]
        manifest["stages"]["gennbo"] = {"result": str((transformed / "gennbo-process.json").resolve()),
                                          "pass": _stage_success(result, files)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["gennbo"]["pass"]:
            raise RuntimeError("gennbo_failed")
        manifest.update(status="producer_complete_unvalidated", paired_fchk=str(paired.resolve()),
                        analysis_archive47=str(archive.resolve()), matrix_archive47=str((transformed / "FILE.47").resolve()),
                        report=str((transformed / "launcher.log").resolve()),
                        matrices={key: str((transformed / f"FILE.{number}").resolve())
                                  for key, number in NUMBER.items()},
                        matrix_sha256={key: sha256(transformed / f"FILE.{number}")
                                       for key, number in NUMBER.items()},
                        max_peak_tree_commit_bytes=max(json.loads(Path(stage["result"]).read_text())["peak_tree_commit_bytes"]
                            for stage in manifest["stages"].values()),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        manifest.update(status="producer_failed_preserved", failure=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    finally:
        manifest["frozen_fchk_sha256_after"] = sha256(source)
        atomic_json(target / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(target.resolve()),
            "status": manifest["status"], "stages": manifest["stages"],
            "failure": manifest.get("failure"),
            "frozen_source_unchanged": manifest["frozen_fchk_sha256_after"] == source_sha}


def two_stage_fchk(case_id: str, attempt: str) -> dict:
    """Pilot a preserved canonical SCF followed by a Guess=Only NBO step."""
    ledger = json.loads((WORK / "aomo-coverage/coverage-ledger.json").read_text())
    row = next((item for item in ledger["cases"] if item["case_id"] == case_id), None)
    if row is None or not row["manifest_hash_match"]:
        raise ValueError(f"frozen_source_not_verified:{case_id}")
    source = Path(row["source_fchk"])
    source_sha = sha256(source)
    if source_sha != row["manifest_sha256"]:
        raise RuntimeError(f"frozen_fchk_hash_changed:{case_id}")
    target = DEST / case_id / attempt
    if target.exists():
        raise FileExistsError(f"new_attempt_required:{target}")
    target.mkdir(parents=True)
    manifest = {"case_id": case_id, "attempt": attempt,
                "kind": "FCHK_seed_to_canonical_SCF_then_GuessOnly_NBO_pilot",
                "frozen_fchk": str(source.resolve()), "frozen_sha256": source_sha,
                "method": row["method"], "basis": row["basis"],
                "charge": row["charge"], "multiplicity": row["multiplicity"],
                "stages": {}, "started_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(target / "manifest.json", manifest)
    try:
        imported = target / "imported"
        canonical = target / "canonical"
        analysis = target / "analysis"
        imported.mkdir(); canonical.mkdir(); analysis.mkdir()
        imported_chk = imported / "imported.chk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "unfchk.exe"), str(source), str(imported_chk)],
                                imported, "unfchk", gaussian=True, timeout=120)
        manifest["stages"]["unfchk"] = {"result": str((imported / "unfchk-process.json").resolve()),
                                          "pass": _stage_success(result, [imported_chk])}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["unfchk"]["pass"]: raise RuntimeError("unfchk_failed")
        canonical_chk = canonical / "canonical.chk"
        shutil.copy2(imported_chk, canonical_chk)
        method = row["method"]
        canonical_inp = canonical / "canonical.gjf"
        canonical_inp.write_text("%nprocshared=3\n%mem=12GB\n%chk=canonical.chk\n"
            f"#p {method}/ChkBasis Geom=AllCheck Guess=Read SCF=(XQC,Tight) "
            "Integral=UltraFine NoSymm\n\n", encoding="ascii", newline="\n")
        canonical_log = canonical / "canonical.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(canonical_inp), str(canonical_log)],
                                canonical, "canonical", gaussian=True, timeout=3600)
        manifest["stages"]["canonical_scf"] = {
            "result": str((canonical / "canonical-process.json").resolve()),
            "pass": _stage_success(result, [canonical_chk, canonical_log], normal_log=canonical_log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["canonical_scf"]["pass"]: raise RuntimeError("canonical_scf_failed")
        canonical_fchk = canonical / "canonical.fchk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "formchk.exe"), str(canonical_chk), str(canonical_fchk)],
                                canonical, "formchk", gaussian=True, timeout=300)
        fchk_text = canonical_fchk.read_text(errors="replace") if canonical_fchk.is_file() else ""
        beta_present = "Beta MO coefficients" in fchk_text
        manifest["stages"]["canonical_formchk"] = {
            "result": str((canonical / "formchk-process.json").resolve()),
            "pass": _stage_success(result, [canonical_fchk]) and
                    (row["multiplicity"] == 1 or beta_present),
            "beta_coefficients_present": beta_present}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["canonical_formchk"]["pass"]: raise RuntimeError("canonical_formchk_incomplete")
        analysis_chk = analysis / "analysis.chk"
        shutil.copy2(canonical_chk, analysis_chk)
        keylist = "PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT"
        analysis_inp = analysis / "analysis.gjf"
        analysis_inp.write_text("%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
            f"#p {method}/ChkBasis Geom=AllCheck Guess=(Read,Only) "
            "Density=Current NoSymm Pop=NBO6Read\n\n"
            f"$NBO {keylist} $END\n\n", encoding="ascii", newline="\n")
        analysis_log = analysis / "analysis.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(analysis_inp), str(analysis_log)],
                                analysis, "guessonly-nbo", gaussian=True, timeout=1200)
        archive = analysis / "FILE.47"
        manifest["stages"]["guessonly_nbo"] = {
            "result": str((analysis / "guessonly-nbo-process.json").resolve()),
            "pass": _stage_success(result, [analysis_log, archive], normal_log=analysis_log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["guessonly_nbo"]["pass"]: raise RuntimeError("guessonly_nbo_failed")
        manifest.update(status="producer_pilot_complete_numerical_identity_unchecked",
                        canonical_fchk=str(canonical_fchk.resolve()),
                        analysis_archive47=str(archive.resolve()),
                        canonical_chk_sha256=sha256(canonical_chk),
                        analysis_chk_sha256=sha256(analysis_chk),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        manifest.update(status="producer_failed_preserved", failure=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    finally:
        manifest["frozen_fchk_sha256_after"] = sha256(source)
        atomic_json(target / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(target.resolve()),
            "status": manifest["status"], "stages": manifest["stages"],
            "failure": manifest.get("failure"),
            "frozen_source_unchanged": manifest["frozen_fchk_sha256_after"] == source_sha}


def _fchk_array(path: Path, label: str, numeric_type: str) -> list[float]:
    lines = path.read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if line.startswith(label):
            match = re.search(rf"\b{numeric_type}\s+N=\s*(\d+)", line)
            if not match:
                raise RuntimeError(f"invalid_FCHK_array:{label}")
            count = int(match.group(1))
            values = []
            for row in lines[i+1:]:
                values.extend(float(word.replace("D", "E")) for word in row.split())
                if len(values) >= count:
                    if len(values) != count:
                        raise RuntimeError(f"FCHK_array_overrun:{label}")
                    return values
            raise RuntimeError(f"truncated_FCHK_array:{label}")
    raise RuntimeError(f"missing_FCHK_array:{label}")


def _fchk_integer(path: Path, label: str) -> int:
    match = re.search(rf"^{re.escape(label)}\s+I\s+(-?\d+)\s*$",
                      path.read_text(errors="replace"), re.M)
    if not match:
        raise RuntimeError(f"missing_FCHK_integer:{label}")
    return int(match.group(1))


def fresh_contact(attempt: str) -> dict:
    """Fresh fixed-geometry OLD-173 calculation; old FCHK supplies geometry only."""
    case_id = "OLD-173"
    ledger = json.loads((WORK / "aomo-coverage/coverage-ledger.json").read_text())
    row = next(item for item in ledger["cases"] if item["case_id"] == case_id)
    source = Path(row["source_fchk"])
    source_sha = sha256(source)
    if source_sha != row["manifest_sha256"]:
        raise RuntimeError("frozen_FCHK_hash_changed")
    target = DEST / case_id / attempt
    if target.exists():
        raise FileExistsError(f"new_attempt_required:{target}")
    target.mkdir(parents=True)
    elements = {1:"H", 6:"C", 7:"N", 9:"F", 51:"Sb"}
    atnums = [int(z) for z in _fchk_array(source, "Atomic numbers", "I")]
    coords = _fchk_array(source, "Current cartesian coordinates", "R")
    if len(atnums) != row["atoms"] or len(coords) != 3 * row["atoms"]:
        raise RuntimeError("frozen_geometry_dimension_mismatch")
    geometry = "\n".join(f"{elements[z]} " + " ".join(
        f"{coords[3*i+j] * 0.529177210903:.10f}" for j in range(3))
        for i,z in enumerate(atnums))
    manifest = {"case_id": case_id, "attempt": attempt,
                "kind": "fresh_fixed_geometry_canonical_SCF_then_GuessOnly_NBO",
                "frozen_fchk": str(source.resolve()), "frozen_fchk_sha256": source_sha,
                "frozen_input_role": "atomic_sequence_charge_multiplicity_and_geometry_only",
                "geometry_angstrom": geometry,
                "method_basis": "PBE1PBE/def2SVP 5D 7F",
                "charge": row["charge"], "multiplicity": row["multiplicity"],
                "expected_nbasis": row["n_basis"], "expected_explicit_electrons": row["explicit_electrons"],
                "expected_effective_core_electrons": row["ecp_core_derived"],
                "stages": {}, "started_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(target / "manifest.json", manifest)
    try:
        canonical = target / "canonical"; analysis = target / "analysis"; matrices = target / "matrices"
        canonical.mkdir(); analysis.mkdir(); matrices.mkdir()
        canonical_chk = canonical / "canonical.chk"
        inp = canonical / "canonical.gjf"
        inp.write_text("%nprocshared=3\n%mem=12GB\n%chk=canonical.chk\n"
            "#p PBE1PBE/def2SVP 5D 7F NoSymm SCF=(XQC,Tight,MaxCycle=1024) Integral=UltraFine\n\n"
            "OLD-173 fresh fixed contact geometry\n\n"
            f"{row['charge']} {row['multiplicity']}\n{geometry}\n\n", encoding="ascii", newline="\n")
        log = canonical / "canonical.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(inp), str(log)],
                                canonical, "canonical", gaussian=True, timeout=3600)
        manifest["stages"]["canonical_scf"] = {"result": str((canonical / "canonical-process.json").resolve()),
            "pass": _stage_success(result, [canonical_chk, log], normal_log=log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["canonical_scf"]["pass"]: raise RuntimeError("canonical_scf_failed")
        canonical_fchk = canonical / "canonical.fchk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "formchk.exe"), str(canonical_chk), str(canonical_fchk)],
                                canonical, "formchk", gaussian=True, timeout=300)
        manifest["stages"]["formchk"] = {"result": str((canonical / "formchk-process.json").resolve()),
            "pass": _stage_success(result, [canonical_fchk])}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["formchk"]["pass"]: raise RuntimeError("formchk_failed")
        nbasis = _fchk_integer(canonical_fchk, "Number of basis functions")
        electrons = _fchk_integer(canonical_fchk, "Number of electrons")
        effective_core = sum(atnums) - row["charge"] - electrons
        manifest["actual_basis_and_ecp"] = {"nbasis": nbasis, "explicit_electrons": electrons,
                                             "effective_core_electrons": effective_core}
        if (nbasis, electrons, effective_core) != (row["n_basis"], row["explicit_electrons"],
                                                   row["ecp_core_derived"]):
            raise RuntimeError("fresh_model_basis_or_ECP_mismatch")
        analysis_chk = analysis / "analysis.chk"
        shutil.copy2(canonical_chk, analysis_chk)
        ai = analysis / "analysis.gjf"
        ai.write_text("%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
            "#p PBE1PBE/ChkBasis Geom=AllCheck Guess=(Read,Only) Density=Current NoSymm Pop=NBO6Read\n\n"
            "$NBO PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT $END\n\n",
            encoding="ascii", newline="\n")
        alog = analysis / "analysis.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(ai), str(alog)],
                                analysis, "guessonly-nbo", gaussian=True, timeout=1200)
        archive = analysis / "FILE.47"
        manifest["stages"]["guessonly_nbo"] = {
            "result": str((analysis / "guessonly-nbo-process.json").resolve()),
            "pass": _stage_success(result, [alog, archive], normal_log=alog)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["guessonly_nbo"]["pass"]: raise RuntimeError("guessonly_nbo_failed")
        original = archive.read_text(encoding="ascii")
        updated, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
        if count != 1: raise RuntimeError("archive_keylist_unexpected")
        minp = matrices / "aomo.47"
        minp.write_text(updated, encoding="ascii", newline="\n")
        result = _bounded_stage([str(NBO_BIN / "gennbo.i8.exe"), minp.name],
                                matrices, "gennbo", gaussian=False, timeout=600)
        files = [matrices / f"FILE.{number}" for number in NUMBER.values()]
        manifest["stages"]["gennbo"] = {"result": str((matrices / "gennbo-process.json").resolve()),
                                          "pass": _stage_success(result, files)}
        if not manifest["stages"]["gennbo"]["pass"]: raise RuntimeError("gennbo_failed")
        manifest.update(status="producer_complete_numerical_identity_unchecked",
            canonical_fchk=str(canonical_fchk.resolve()), analysis_archive47=str(archive.resolve()),
            matrix_archive47=str((matrices / "FILE.47").resolve()),
            report=str((matrices / "launcher.log").resolve()),
            matrices={key: str((matrices / f"FILE.{number}").resolve()) for key,number in NUMBER.items()},
            finished_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        manifest.update(status="producer_failed_preserved", failure=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    finally:
        manifest["frozen_fchk_sha256_after"] = sha256(source)
        atomic_json(target / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(target.resolve()),
            "status": manifest["status"], "failure": manifest.get("failure"),
            "actual_basis_and_ecp": manifest.get("actual_basis_and_ecp"),
            "stages": manifest["stages"],
            "frozen_source_unchanged": manifest["frozen_fchk_sha256_after"] == source_sha}


def explicit_geometry_contact_retry(attempt: str) -> dict:
    """Retry only the NBO step, avoiding OLD-173's failing Geom=AllCheck path."""
    case_id = "OLD-173"
    prior = DEST / case_id / "fresh-geometry-01"
    previous = json.loads((prior / "manifest.json").read_text())
    if not previous["stages"]["canonical_scf"]["pass"] or not previous["stages"]["formchk"]["pass"]:
        raise RuntimeError("fresh_canonical_not_complete")
    canonical_chk = prior / "canonical/canonical.chk"
    canonical_fchk = prior / "canonical/canonical.fchk"
    target = DEST / case_id / attempt
    if target.exists(): raise FileExistsError(f"new_attempt_required:{target}")
    analysis = target / "analysis"; matrices = target / "matrices"
    analysis.mkdir(parents=True); matrices.mkdir()
    analysis_chk = analysis / "analysis.chk"
    shutil.copy2(canonical_chk, analysis_chk)
    manifest = {"case_id": case_id, "attempt": attempt,
                "kind": "same_fresh_canonical_checkpoint_explicit_geometry_GuessOnly_NBO_retry",
                "canonical_fchk": str(canonical_fchk.resolve()),
                "canonical_fchk_sha256": sha256(canonical_fchk),
                "canonical_checkpoint_source": str(canonical_chk.resolve()),
                "canonical_checkpoint_source_sha256": sha256(canonical_chk),
                "analysis_checkpoint_copy_sha256": sha256(analysis_chk),
                "geometry_angstrom": previous["geometry_angstrom"],
                "method_basis": previous["method_basis"], "stages": {},
                "started_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(target / "manifest.json", manifest)
    try:
        inp = analysis / "analysis.gjf"
        inp.write_text("%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
            "#p PBE1PBE/ChkBasis Guess=(Read,Only) Density=Current NoSymm Pop=NBO6Read\n\n"
            "OLD-173 explicit contact geometry NBO analysis\n\n"
            f"{previous['charge']} {previous['multiplicity']}\n{previous['geometry_angstrom']}\n\n"
            "$NBO PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT $END\n\n",
            encoding="ascii", newline="\n")
        log = analysis / "analysis.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(inp), str(log)],
                                analysis, "explicit-nbo", gaussian=True, timeout=1200)
        archive = analysis / "FILE.47"
        manifest["stages"]["explicit_nbo"] = {"result": str((analysis / "explicit-nbo-process.json").resolve()),
            "pass": _stage_success(result, [log, archive], normal_log=log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["explicit_nbo"]["pass"]: raise RuntimeError("explicit_nbo_failed")
        original = archive.read_text(encoding="ascii")
        updated, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
        if count != 1: raise RuntimeError("archive_keylist_unexpected")
        minp = matrices / "aomo.47"
        minp.write_text(updated, encoding="ascii", newline="\n")
        result = _bounded_stage([str(NBO_BIN / "gennbo.i8.exe"), minp.name],
                                matrices, "gennbo", gaussian=False, timeout=600)
        files = [matrices / f"FILE.{number}" for number in NUMBER.values()]
        manifest["stages"]["gennbo"] = {"result": str((matrices / "gennbo-process.json").resolve()),
                                          "pass": _stage_success(result, files)}
        if not manifest["stages"]["gennbo"]["pass"]: raise RuntimeError("gennbo_failed")
        manifest.update(status="producer_complete_numerical_identity_unchecked",
            analysis_archive47=str(archive.resolve()),
            matrix_archive47=str((matrices / "FILE.47").resolve()),
            report=str((matrices / "launcher.log").resolve()),
            matrices={key: str((matrices / f"FILE.{number}").resolve()) for key,number in NUMBER.items()},
            finished_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        manifest.update(status="producer_failed_preserved", failure=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    finally:
        atomic_json(target / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(target.resolve()),
            "status": manifest["status"], "failure": manifest.get("failure"),
            "stages": manifest["stages"]}


def g_shell_variant(variant: str, attempt: str) -> dict:
    """Pair one preserved H2+ FCHK format variant with a new NBO analysis.

    The source is copied before unfchk. The canonical FCHK is saved before the
    Guess=Only NBO step, so its alpha/beta coefficient arrays remain intact.
    Pure and Cartesian runs are format variants of the same chemical species.
    """
    if variant not in ("pure", "cartesian"):
        raise ValueError("g_shell_variant_must_be_pure_or_cartesian")
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,50}", attempt):
        raise ValueError("attempt_must_be_short_lowercase_hyphenated_name")
    case_id = "G-PURE" if variant == "pure" else "G-CART"
    source = Path("E:/Codex/2026-09-05/branch-15/outputs/cov-complete-validation-20260906"
                  "/regression-fixtures/g-shell-v2") / variant / "wavefunction.fch"
    source_sha = sha256(source)
    basis = _fchk_integer(source, "Number of basis functions")
    expected_basis = 110 if variant == "pure" else 140
    if (basis != expected_basis or _fchk_integer(source, "Number of atoms") != 2 or
            _fchk_integer(source, "Charge") != 1 or
            _fchk_integer(source, "Multiplicity") != 2 or
            _fchk_integer(source, "Number of alpha electrons") != 1 or
            _fchk_integer(source, "Number of beta electrons") != 0):
        raise RuntimeError("g_shell_source_identity_mismatch")
    target = COMPUTE / "g-shell" / case_id / attempt
    if target.exists():
        raise FileExistsError(f"new_attempt_required:{target}")
    target.mkdir(parents=True)
    manifest = {
        "case_id": case_id, "variant": variant, "attempt": attempt,
        "kind": "same_H2plus_format_variant_FCHK_seeded_canonical_SCF_then_GuessOnly_NBO",
        "independent_chemical_species_increment": 0,
        "source_fchk": str(source.resolve()), "source_sha256": source_sha,
        "source_role": "seed_only_not_identical_canonical_wavefunction",
        "method_basis": "UPBE1PBE/cc-pV5Z", "charge": 1, "multiplicity": 2,
        "expected_nbasis": expected_basis, "expected_alpha_electrons": 1,
        "expected_beta_electrons": 0, "stages": {},
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(target / "manifest.json", manifest)
    try:
        imported = target / "imported"
        canonical = target / "canonical"
        analysis = target / "analysis"
        matrices = target / "matrices"
        for folder in (imported, canonical, analysis, matrices):
            folder.mkdir()
        copied_source = imported / "source.fchk"
        shutil.copy2(source, copied_source)
        if sha256(copied_source) != source_sha:
            raise RuntimeError("source_copy_hash_mismatch")
        imported_chk = imported / "imported.chk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "unfchk.exe"), str(copied_source),
                                 str(imported_chk)], imported, "unfchk", gaussian=True,
                                timeout=120)
        manifest["stages"]["unfchk"] = {
            "result": str((imported / "unfchk-process.json").resolve()),
            "pass": _stage_success(result, [imported_chk])}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["unfchk"]["pass"]:
            raise RuntimeError("unfchk_failed")

        canonical_chk = canonical / "canonical.chk"
        shutil.copy2(imported_chk, canonical_chk)
        canonical_inp = canonical / "canonical.gjf"
        canonical_inp.write_text(
            "%nprocshared=3\n%mem=12GB\n%chk=canonical.chk\n"
            "#p UPBE1PBE/ChkBasis Geom=AllCheck Guess=Read SCF=(XQC,Tight) "
            "Integral=UltraFine NoSymm\n\n", encoding="ascii", newline="\n")
        canonical_log = canonical / "canonical.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(canonical_inp),
                                 str(canonical_log)], canonical, "canonical",
                                gaussian=True, timeout=3600)
        manifest["stages"]["canonical_scf"] = {
            "result": str((canonical / "canonical-process.json").resolve()),
            "pass": _stage_success(result, [canonical_chk, canonical_log],
                                   normal_log=canonical_log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["canonical_scf"]["pass"]:
            raise RuntimeError("canonical_scf_failed")

        canonical_fchk = canonical / "canonical.fchk"
        result = _bounded_stage([str(GAUSSIAN_BIN / "formchk.exe"),
                                 str(canonical_chk), str(canonical_fchk)],
                                canonical, "formchk", gaussian=True, timeout=300)
        fchk_text = canonical_fchk.read_text(errors="replace") if canonical_fchk.is_file() else ""
        observed_basis = _fchk_integer(canonical_fchk, "Number of basis functions") if fchk_text else None
        observed_independent = (_fchk_integer(canonical_fchk, "Number of independent functions")
                                if fchk_text else None)
        observed_alpha = (_fchk_integer(canonical_fchk, "Number of alpha electrons")
                          if fchk_text else None)
        observed_beta = (_fchk_integer(canonical_fchk, "Number of beta electrons")
                         if fchk_text else None)
        observed_pure_d = (_fchk_integer(canonical_fchk, "Pure/Cartesian d shells")
                           if fchk_text else None)
        observed_pure_f = (_fchk_integer(canonical_fchk, "Pure/Cartesian f shells")
                           if fchk_text else None)
        expected_pure = 0 if variant == "pure" else 1
        manifest["stages"]["canonical_formchk"] = {
            "result": str((canonical / "formchk-process.json").resolve()),
            "pass": _stage_success(result, [canonical_fchk]) and
                    observed_basis == expected_basis and observed_independent == expected_basis and
                    observed_alpha == 1 and observed_beta == 0 and
                    observed_pure_d == expected_pure and observed_pure_f == expected_pure and
                    "Alpha MO coefficients" in fchk_text and
                    "Beta MO coefficients" in fchk_text,
            "observed_nbasis": observed_basis,
            "observed_nindependent": observed_independent,
            "observed_alpha_electrons": observed_alpha,
            "observed_beta_electrons": observed_beta,
            "observed_pure_cartesian_d": observed_pure_d,
            "observed_pure_cartesian_f": observed_pure_f,
            "alpha_coefficients_present": "Alpha MO coefficients" in fchk_text,
            "beta_coefficients_present": "Beta MO coefficients" in fchk_text}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["canonical_formchk"]["pass"]:
            raise RuntimeError("canonical_formchk_incomplete_or_basis_changed")

        analysis_chk = analysis / "analysis.chk"
        shutil.copy2(canonical_chk, analysis_chk)
        analysis_chk_copy_sha256 = sha256(analysis_chk)
        if analysis_chk_copy_sha256 != sha256(canonical_chk):
            raise RuntimeError("analysis_checkpoint_copy_hash_mismatch")
        analysis_inp = analysis / "analysis.gjf"
        analysis_inp.write_text(
            "%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
            "#p UPBE1PBE/ChkBasis Geom=AllCheck Guess=(Read,Only) "
            "Density=Current NoSymm Pop=NBO6Read\n\n"
            "$NBO PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT $END\n\n",
            encoding="ascii", newline="\n")
        analysis_log = analysis / "analysis.log"
        result = _bounded_stage([str(GAUSSIAN_BIN / "g16.exe"), str(analysis_inp),
                                 str(analysis_log)], analysis, "guessonly-nbo",
                                gaussian=True, timeout=1200)
        archive = analysis / "FILE.47"
        manifest["stages"]["guessonly_nbo"] = {
            "result": str((analysis / "guessonly-nbo-process.json").resolve()),
            "pass": _stage_success(result, [analysis_log, archive],
                                   normal_log=analysis_log)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["guessonly_nbo"]["pass"]:
            raise RuntimeError("guessonly_nbo_failed")

        original = archive.read_text(encoding="ascii")
        updated, count = re.subn(r"\$NBO\s+\$END", f"$NBO {KEYLIST} $END", original, count=1)
        if count != 1:
            raise RuntimeError("archive_keylist_unexpected")
        matrix_inp = matrices / "aomo.47"
        matrix_inp.write_text(updated, encoding="ascii", newline="\n")
        result = _bounded_stage([str(NBO_BIN / "gennbo.i8.exe"), matrix_inp.name],
                                matrices, "gennbo", gaussian=False, timeout=600)
        files = [matrices / f"FILE.{number}" for number in NUMBER.values()]
        manifest["stages"]["gennbo"] = {
            "result": str((matrices / "gennbo-process.json").resolve()),
            "pass": _stage_success(result, files)}
        atomic_json(target / "manifest.json", manifest)
        if not manifest["stages"]["gennbo"]["pass"]:
            raise RuntimeError("gennbo_failed_or_missing_matrix")

        manifest.update(
            status="producer_complete_numerical_identity_unchecked",
            canonical_fchk=str(canonical_fchk.resolve()),
            analysis_archive47=str(archive.resolve()),
            matrix_archive47=str((matrices / "FILE.47").resolve()),
            report=str((matrices / "launcher.log").resolve()),
            matrices={key: str((matrices / f"FILE.{number}").resolve())
                      for key, number in NUMBER.items()},
            canonical_chk_sha256=sha256(canonical_chk),
            analysis_chk_copy_sha256=analysis_chk_copy_sha256,
            analysis_chk_after_nbo_sha256=sha256(analysis_chk),
            max_peak_tree_commit_bytes=max(
                json.loads(Path(stage["result"]).read_text())["peak_tree_commit_bytes"]
                for stage in manifest["stages"].values()),
            finished_utc=datetime.now(timezone.utc).isoformat())
    except Exception as error:
        manifest.update(status="producer_failed_preserved", failure=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat())
    finally:
        manifest["source_sha256_after"] = sha256(source)
        atomic_json(target / "manifest.json", manifest)
    return {"case_id": case_id, "directory": str(target.resolve()),
            "status": manifest["status"], "failure": manifest.get("failure"),
            "source_unchanged": manifest["source_sha256_after"] == source_sha,
            "stages": manifest["stages"]}


def publish_index() -> dict:
    ledger = json.loads((WORK / "aomo-coverage/coverage-ledger.json").read_text())
    two_stage = json.loads((COMPUTE / "two-stage-identity.json").read_text())
    two_stage_pass = {item["case_id"] for item in two_stage["cases"] if item["pass"]}
    rows = []
    for source in ledger["cases"]:
        if not source["representative_for"]:
            continue
        case_id = source["case_id"]
        base = DEST / case_id
        if case_id == "OLD-173":
            selected = base / "explicit-geometry-01"
            production = json.loads((selected / "manifest.json").read_text())
            fchk = Path(production["canonical_fchk"])
            matrix_root = selected / "matrices"
            identity = COMPUTE / "old173-fresh-identity.json"
            relationship = "new_fixed_geometry_SCF_checkpoint_snapshot_before_GuessOnly_NBO"
            source_role = "geometry_atoms_charge_multiplicity_only"
            source_step = "fresh-geometry-01/canonical"
            nbo_step = "explicit-geometry-01/analysis"
        elif case_id in two_stage_pass:
            selected = base / "aomo-from-two-stage-01"
            production = json.loads((selected / "manifest.json").read_text())
            fchk = base / "two-stage-01/canonical/canonical.fchk"
            matrix_root = selected
            identity = COMPUTE / "two-stage-identity.json"
            relationship = "same_SCF_checkpoint_snapshot_before_GuessOnly_NBO; different_analysis_step"
            source_role = "frozen_FCHK_seed_only"
            source_step = "two-stage-01/canonical"
            nbo_step = "two-stage-01/analysis"
        else:
            selected = base / "fchk-seeded-01"
            production = json.loads((selected / "manifest.json").read_text())
            fchk = Path(production["paired_fchk"])
            matrix_root = selected / "matrices"
            identity = None
            relationship = "same_new_Gaussian_NBO_analysis_checkpoint"
            source_role = "frozen_FCHK_seed_only"
            source_step = "fchk-seeded-01/analysis"
            nbo_step = "fchk-seeded-01/analysis"
        if not production["status"].startswith("producer_complete_"):
            raise RuntimeError(f"selected_case_incomplete:{case_id}")
        if not fchk.is_file() or not (matrix_root / "FILE.47").is_file():
            raise RuntimeError(f"selected_pair_missing:{case_id}")
        matrix_paths = {key.lower(): str((matrix_root / f"FILE.{number}").resolve())
                        for key, number in NUMBER.items()}
        missing = [key for key, path in matrix_paths.items() if not Path(path).is_file()]
        if missing:
            raise RuntimeError(f"selected_matrices_missing:{case_id}:{missing}")
        attempts = sorted(path for path in base.glob("*/manifest.json") if path.is_file())
        failed = [str(path.resolve()) for path in attempts
                  if json.loads(path.read_text()).get("status") == "producer_failed_preserved"]
        record = {
            "case_id": case_id, "covers": source["representative_for"],
            "source_group": source["source_group"], "geometry_label": source["geometry_label"],
            "identity_key": source["identity_key"], "count_policy": source["count_policy"],
            "chemical_tags": source["chemical_tags"],
            "mathematical_tags": source["mathematical_tags"],
            "required_NBO_capabilities": source["required_NBO_capabilities"],
            "user_visible_acceptance": source["user_visible_acceptance"],
            "method": source["method"], "basis": source["basis"],
            "charge": source["charge"], "multiplicity": source["multiplicity"],
            "nbasis": source["n_basis"], "nindependent": source["n_independent"],
            "explicit_electrons": source["explicit_electrons"],
            "effective_core_electrons": source["ecp_core_derived"],
            "frozen_source_fchk": source["source_fchk"],
            "frozen_source_sha256": source["manifest_sha256"],
            "frozen_source_role": source_role,
            "fchk": str(fchk.resolve()), "fchk_sha256": sha256(fchk),
            "report": str((matrix_root / "launcher.log").resolve()),
            "archive47": str((matrix_root / "FILE.47").resolve()),
            "source_archive47": production.get("source_archive47", production.get("analysis_archive47")),
            **matrix_paths,
            "selected_producer_manifest": str((selected / "manifest.json").resolve()),
            "failed_attempt_manifests": failed,
            "source_step": source_step, "nbo_step": nbo_step,
            "matrix_step": matrix_root.name,
            "fchk_relationship": relationship,
            "independent_pair_identity_evidence": str(identity.resolve()) if identity else None,
            "producer_status": "complete_unvalidated",
            "scientific_status": "pending_independent_case_audit",
            "user_visible_status": "pending_product_acceptance",
        }
        rows.append(record)
    if len(rows) != 61:
        raise RuntimeError(f"expected_61_representatives_found_{len(rows)}")
    OUTPUT = DEST / "cases.json"
    atomic_json(OUTPUT, rows)
    summary = {"schema": "aomo_representative_production_summary_v1",
               "case_count": len(rows),
               "same_new_Gaussian_NBO_analysis_checkpoint": sum("same_new_Gaussian" in x["fchk_relationship"] for x in rows),
               "same_SCF_snapshot_different_analysis_step": sum("snapshot" in x["fchk_relationship"] for x in rows),
               "frozen_source_same_SCF_claims": 0,
               "scientific_status": "pending_independent_case_audit",
               "user_visible_status": "pending_product_acceptance",
               "cases_json": str(OUTPUT.resolve()),
               "geometry_label_count": len({x["geometry_label"] for x in rows if x["geometry_label"]}),
               "failed_attempt_count": sum(len(x["failed_attempt_manifests"]) for x in rows),
               "created_utc": datetime.now(timezone.utc).isoformat()}
    atomic_json(COMPUTE / "production-summary.json", summary)
    acceptance = {
        "schema": "aomo_representative_acceptance_run_index_v1",
        "scope": "61_deduplicated_source_cases_not_all_273",
        "frozen_corpus_case_count": 273,
        "frozen_corpus_FCHK_sha256_verified": 273,
        "selected_representative_count": len(rows),
        "production_status": "61_raw_packages_complete",
        "scientific_status": "pending_independent_per_case_checks",
        "product_status": "pending_real_viewer_interaction_and_visual_review",
        "required_first_gate": ["OLD-063", "OLD-068", "OLD-070"],
        "geometry_cases": {x["geometry_label"]: x["case_id"] for x in rows if x["geometry_label"]},
        "ligand_field_cases": [x["case_id"] for x in rows if x["source_group"] == "ligand_field"],
        "mathematical_boundary_cases": {tag: [x["case_id"] for x in rows if tag in x["mathematical_tags"]]
             for tag in sorted({tag for x in rows for tag in x["mathematical_tags"]})},
        "case_runs": [{
            "case_id": x["case_id"], "source_group": x["source_group"],
            "covers": x["covers"], "chemical_tags": x["chemical_tags"],
            "geometry_label": x["geometry_label"],
            "mathematical_tags": x["mathematical_tags"],
            "identity_key": x["identity_key"], "count_policy": x["count_policy"],
            "producer_manifest": x["selected_producer_manifest"],
            "fchk_relationship": x["fchk_relationship"],
            "independent_pair_identity_evidence": x["independent_pair_identity_evidence"],
            "numeric_acceptance": ["same_wavefunction_and_AO_order",
                "raw_AO_NAO_NHO_NBO_NLMO_MO_matrix_chain",
                "overlap_metric_effective_rank_and_null_columns",
                "canonical_coefficient_energy_occupation_spin_identity",
                "explicit_electron_and_ECP_core_bookkeeping",
                "orbital_field_and_signed_component_reconstruction"],
            "chemical_data_acceptance": x["required_NBO_capabilities"],
            "visible_acceptance": x["user_visible_acceptance"],
            "producer_status": x["producer_status"],
            "numeric_status": "pending",
            "chemical_status": "pending",
            "visible_status": "pending",
        } for x in rows],
    }
    atomic_json(WORK / "aomo-coverage/acceptance-run-index.json", acceptance)
    csv_path = WORK / "aomo-coverage/acceptance-run-index.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        columns = ["case_id", "source_group", "covers", "chemical_tags", "geometry_label",
                   "mathematical_tags", "count_policy", "fchk_relationship",
                   "producer_status", "numeric_status", "chemical_status", "visible_status"]
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for case in acceptance["case_runs"]:
            writer.writerow({key: ";".join(case[key]) if isinstance(case.get(key), list)
                             else case.get(key) for key in columns})
    gaps = {"schema": "aomo_acceptance_gaps_v1",
            "represented_source_case_count": len(rows),
            "unrepresented_selected_cases": [],
            "unrepresented_full_geometry_labels": [],
            "producer_gaps": [],
            "scientific_gaps": [
                "Independent raw-matrix, FCHK-identity, spin, ECP, effective-rank and field checks are pending on all 61 selected cases.",
                "The 14 checkpoint-snapshot cases have producer identity evidence but still need the common independent case audit.",
                "OLD-020 has NBAS 400 and 396 independent orbitals; its null columns require explicit effective-rank acceptance.",
                "Historical frozen FCHK files are source seeds or geometry sources, not declared paired with the new NBO analyses.",
            ],
            "product_gaps": [
                "The ordinary viewer one-drop, clickable whole AO-MO graph, right-side 3D, reverse picking and exports require visible end-to-end acceptance.",
                "NPA/Wiberg, directional NHO, E2, NLMO, multicentre and coordination overlays require per-class visual and semantic checks.",
                "Focus filtering, degenerate subspaces, spin switching and molecule replacement require stale-state and selection checks.",
            ],
            "corpus_format_boundaries": [
                "No Cartesian-g source in this 273-case corpus; use the existing g-shell fixture or record no positive corpus coverage.",
                "No pure RHF/ROHF source family in this 273-case corpus; present DFT R/U checks do not imply those method families passed.",
                "Non-Aufbau, fractional-occupation, spinor/complex and unsupported higher-angular-momentum boundaries need explicit rejection or separate fixtures, not invented positive chemistry cases.",
            ],
            "retained_failure_history": [path for x in rows for path in x["failed_attempt_manifests"]],
            "failed_attempt_manifest_count": summary["failed_attempt_count"],
            "status_rule": "producer_complete_does_not_imply_scientific_or_product_pass"}
    atomic_json(WORK / "aomo-coverage/gap-list.json", gaps)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", nargs="?", help="original NBO number 01..10, or OLD-nnn")
    parser.add_argument("--attempt", default="aomo-matrix-01")
    parser.add_argument("--recover-fchk", action="store_true")
    parser.add_argument("--reanalyse-seeded", metavar="SOURCE_ATTEMPT")
    parser.add_argument("--native-symmetry", action="store_true",
                        help="retry a seeded case without the NoSymm route option")
    parser.add_argument("--two-stage", action="store_true",
                        help="pilot full canonical FCHK followed by Guess=Only NBO analysis")
    parser.add_argument("--fresh-contact", action="store_true",
                        help="fresh fixed-geometry OLD-173 retry with verified basis and ECP")
    parser.add_argument("--explicit-contact-retry", action="store_true",
                        help="OLD-173 fresh checkpoint NBO retry with explicit geometry")
    parser.add_argument("--g-shell-variant", choices=("pure", "cartesian"),
                        help="seed one H2+ g-shell format variant from a preserved FCHK")
    parser.add_argument("--index", action="store_true", help="publish the 61-case production index")
    args = parser.parse_args()
    if args.index:
        result = publish_index()
        print(json.dumps(result, indent=2))
        return 0
    if args.case is None and args.g_shell_variant is None:
        parser.error("case is required unless --index is used")
    if args.g_shell_variant:
        if args.case is not None:
            parser.error("omit case when using --g-shell-variant")
        result = g_shell_variant(args.g_shell_variant, args.attempt)
    elif args.explicit_contact_retry:
        if args.case != "OLD-173":
            parser.error("--explicit-contact-retry is limited to OLD-173")
        result = explicit_geometry_contact_retry(args.attempt)
    elif args.fresh_contact:
        if args.case != "OLD-173":
            parser.error("--fresh-contact is limited to OLD-173")
        result = fresh_contact(args.attempt)
    elif args.two_stage:
        result = two_stage_fchk(args.case, args.attempt)
    elif args.reanalyse_seeded:
        result = reanalyse_seeded_archive(args.case, args.reanalyse_seeded, args.attempt)
    elif args.recover_fchk:
        result = recover_fchk(args.case, args.attempt, native_symmetry=args.native_symmetry)
    else:
        result = run_archive(args.case, args.attempt)
    print(json.dumps(result, indent=2))
    return 0 if result["status"].startswith("producer_complete_") or result["status"].startswith("producer_pilot_complete_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
