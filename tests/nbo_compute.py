"""Isolated, resource-bounded NBO integration calculations.

Only this module owns the new NBO computation directories. The existing
validation process supervisor remains unchanged and is imported read-only.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import re
import math

from validation_process import atomic_json, physical_core_masks, run_tree

ROOT = Path(__file__).resolve().parents[2]
NBO_BIN = Path(r"D:\CalChem\NBO\bin")
GAUSSIAN_BIN = Path(r"D:\CalChem\Gaussian\Gaussian 16 W")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gennbo_smoke() -> int:
    target = ROOT / "gennbo-smoke" / "ch3nh2-original"
    target.mkdir(parents=True, exist_ok=True)
    source = Path(r"D:\CalChem\NBO\ch3nh2.47")
    inp = target / "ch3nh2.47"
    if inp.exists() and sha256(inp) != sha256(source):
        raise RuntimeError("Smoke input exists with unexpected contents")
    if not inp.exists():
        shutil.copy2(source, inp)
    env = dict(os.environ)
    env.update(NBOEXE=str(NBO_BIN / "nbo7.i8.exe"), NBOMEM="512mb",
               GAUSS_SCRDIR=str(target))
    cores = physical_core_masks(3)
    args = [str(NBO_BIN / "gennbo.i8.exe"), inp.name]
    event_path = target / "runner-live.json"
    result = run_tree(args, target, env, sum(cores), 20, 300,
                      on_started=lambda data: atomic_json(event_path, data),
                      affinity=True)
    atomic_json(target / "runner-result.json", result)
    output = target / "launcher.log"
    shutil.copy2(output, target / "ch3nh2.nbo")
    atomic_json(target / "manifest.json", {
        "kind": "vendor_sample_smoke", "source": str(source),
        "source_sha256": sha256(source), "input_sha256": sha256(inp),
        "gennbo_sha256": sha256(NBO_BIN / "gennbo.i8.exe"),
        "nbo7_sha256": sha256(NBO_BIN / "nbo7.i8.exe"),
        "result": {k: result[k] for k in ("exit_code", "timed_out", "tree_process_count",
                                     "peak_tree_commit_bytes", "core_count")},
        "output_size": output.stat().st_size,
    })
    print(json.dumps({"target": str(target), "result": result["exit_code"],
                      "timed_out": result["timed_out"], "size": output.stat().st_size}))
    return 0 if result["exit_code"] == 0 and not result["timed_out"] else 1


def ensure_bridge() -> Path:
    """Install only the A.03-required batch bridge, preserving any existing file."""
    target = GAUSSIAN_BIN / "gaunbo6.bat"
    body = ("@echo off\r\n"
            "@setlocal\r\n"
            "set \"GAUNBO=g16nbo\"\r\n"
            "set \"INT=i8\"\r\n"
            "set \"NBOBIN=D:\\CalChem\\NBO\\bin\"\r\n"
            "set \"NBOMEM=512mb\"\r\n"
            "set \"INTERFACE=%NBOBIN%\\%GAUNBO%.%INT%.exe\"\r\n"
            "set \"NBOEXE=%NBOBIN%\\nbo7.%INT%.exe\"\r\n"
            "set \"inpfile=%~2\"\r\n"
            "set \"msgfile=%~4\"\r\n"
            "set \"fchkfile=%~5\"\r\n"
            "set \"matfile=%~6\"\r\n"
            "\"%INTERFACE%\"\r\n").encode("ascii")
    record_path = ROOT / "compute" / "bridge.json"
    previous_record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.exists() else {}
    if target.exists():
        if target.read_bytes() != body:
            raise RuntimeError(f"Existing bridge differs; refusing overwrite: {target}")
        created = bool(previous_record.get("created_by_nbo_compute")) if (
            previous_record.get("path") == str(target)
            and previous_record.get("sha256") == hashlib.sha256(body).hexdigest()
        ) else False
    else:
        with target.open("xb") as stream:
            stream.write(body)
        created = True
    record = {
        "path": str(target), "created_by_nbo_compute": created,
        "sha256": sha256(target), "bytes": len(body),
        "source_template": str(NBO_BIN / "gaunbo6.bat"),
        "source_template_sha256": sha256(NBO_BIN / "gaunbo6.bat"),
        "guide": r"D:\CalChem\NBO\install\INSTALL.G16W-A03",
        "rollback": "Remove only if SHA256 still matches this record and no run uses it.",
    }
    for key in ("first_created_utc", "creation_evidence", "metadata_repair"):
        if key in previous_record:
            record[key] = previous_record[key]
    if created and not previous_record:
        record["first_created_utc"] = datetime.fromtimestamp(
            target.stat().st_ctime, timezone.utc
        ).isoformat()
        record["creation_evidence"] = "Created with exclusive x mode after target did not exist."
    atomic_json(record_path, record)
    return target


METHYLAMINE = """C  0.745914  0.011106  0.000000
N -0.721743 -0.071848  0.000000
H  1.042059  1.060105  0.000000
H  1.129298 -0.483355  0.892539
H  1.129298 -0.483355 -0.892539
H -1.076988  0.386322 -0.827032
H -1.076988  0.386322  0.827032"""
METHYL = """C  0.000000  0.000000  0.000000
H  0.000000  1.073600  0.000000
H  0.929765 -0.536800  0.000000
H -0.929765 -0.536800  0.000000"""


def benzene_geometry() -> str:
    atoms = []
    for symbol, radius in (("C", 1.397), ("H", 2.482)):
        for i in range(6):
            angle = i * math.pi / 3
            atoms.append(f"{symbol} {radius * math.cos(angle):.6f} "
                         f"{radius * math.sin(angle):.6f} 0.000000")
    return "\n".join(atoms)


CASES = {
    "01": ("water", "O 0.000000 0.000000 0.110851\nH 0.000000 0.783837 -0.443405\nH 0.000000 -0.783837 -0.443405", 0, 1,
           "6-31G(d)", None, r"D:\CalChem\NBO\tests\gaussian\g16w-a03\h2o.gjf"),
    "02": ("methylamine", METHYLAMINE, 0, 1, "6-31G(d)", None,
           r"D:\CalChem\NBO\ch3nh2.47"),
    "03": ("formaldehyde", "C 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.210000\nH 0.940000 0.000000 -0.570000\nH -0.940000 0.000000 -0.570000", 0, 1,
           "6-31G(d)", None, "constructed planar C2v: C=O 1.21 A, C-H 1.10 A"),
    "04": ("acetylene", "H 0.000000 0.000000 -1.663000\nC 0.000000 0.000000 -0.601000\nC 0.000000 0.000000 0.601000\nH 0.000000 0.000000 1.663000", 0, 1,
           "6-31G(d)", None, "constructed linear: C-C 1.202 A, C-H 1.062 A"),
    "05": ("benzene", benzene_geometry(), 0, 1, "6-31G(d)", None,
           "constructed regular D6h ring: C radius 1.397 A, H radius 2.482 A"),
    "06": ("formate", "C 0.000000 0.000000 0.000000\nO 1.113000 0.568000 0.000000\nO -1.113000 0.568000 0.000000\nH 0.000000 -1.100000 0.000000", -1, 1,
           "6-31+G(d)", None, "constructed planar C2v: C-O 1.25 A, C-H 1.10 A"),
    "07": ("methyl radical", METHYL, 0, 2, "6-31G(d)", None,
           r"D:\CalChem\NBO\tests\gaussian\g16w-a03\ch3.gjf"),
    "08": ("ammonia borane", "B -0.800000 0.000000 0.000000\nN 0.800000 0.000000 0.000000\nH -1.203000 1.140000 0.000000\nH -1.203000 -0.570000 0.987269\nH -1.203000 -0.570000 -0.987269\nH 1.140000 0.940000 0.000000\nH 1.140000 -0.470000 0.814064\nH 1.140000 -0.470000 -0.814064", 0, 1,
           "6-31G(d)", None, "constructed staggered H3B-NH3: B-N 1.60 A; B-H ~1.21 A; N-H ~1.00 A"),
    "09": ("diborane", "B -0.885000 0.000000 0.000000\nB 0.885000 0.000000 0.000000\nH 0.000000 1.000000 0.000000\nH 0.000000 -1.000000 0.000000\nH -1.535000 0.000000 1.000000\nH -1.535000 0.000000 -1.000000\nH 1.535000 0.000000 1.000000\nH 1.535000 0.000000 -1.000000", 0, 1,
           "6-31G(d)", None, "constructed bridged D2h: B-B 1.77 A; B-Hbridge 1.34 A; B-Hterminal 1.19 A"),
    "10": ("copper dimer", "Cu 0.000000 0.000000 -1.110000\nCu 0.000000 0.000000 1.110000", 0, 1,
           "GenECP", "LANL2DZ", "constructed linear Cu-Cu 2.22 A, model-only singlet"),
}


def gauss_env(scratch: Path) -> dict:
    env = dict(os.environ)
    env.update(GAUSS_EXEDIR=str(GAUSSIAN_BIN), GAUSS_SCRDIR=str(scratch),
               OMP_NUM_THREADS="3", NCPUS="3", OMP_THREAD_LIMIT="3",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="3")
    env["PATH"] = str(GAUSSIAN_BIN) + os.pathsep + env.get("PATH", "")
    return env


def execute_stage(args: list[Path], directory: Path, name: str, timeout: int = 1800) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    cpu_mask = sum(physical_core_masks(3))
    result_path = directory / f"{name}-process.json"
    if result_path.exists():
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result["exit_code"] == 0 and not result["timed_out"]:
            return result
        raise RuntimeError(f"Prior failed stage retained: {result_path}")
    result = run_tree(args, directory, gauss_env(directory), cpu_mask, 20,
                      timeout, affinity=False,
                      on_started=lambda data: atomic_json(directory / f"{name}-live.json", data))
    atomic_json(result_path, result)
    return result


def gaussian_smoke(name: str) -> int:
    ensure_bridge()
    number = {"ch3nh2": "02", "ch3": "07"}.get(name, name)
    title, atoms, charge, mult, basis, ecp, geometry_source = CASES[number]
    case_id = "NBO-" + number
    target = ROOT / "nbo-cases" / case_id
    canonical = target / "canonical"
    analysis = target / ("analysis-03" if number == "10" else "analysis-02")
    for path in (canonical, analysis):
        path.mkdir(parents=True, exist_ok=True)
    method = "UPBE1PBE" if mult > 1 else "PBE1PBE"
    route = f"#p {method}/{basis} SCF=(XQC,Tight) Integral=UltraFine NoSymm"
    basis_tail = "Cu 0\nLANL2DZ\n****\n\nCu 0\nLANL2DZ\n\n" if ecp else ""
    canonical_input = (f"%nprocshared=3\n%mem=12GB\n%chk=canonical.chk\n"
                       f"{route}\n\n{title}; fixed geometry NBO integration\n\n"
                       f"{charge} {mult}\n{atoms}\n\n{basis_tail}")
    inp = canonical / "canonical.gjf"
    if inp.exists() and inp.read_text(encoding="ascii") != canonical_input:
        raise RuntimeError(f"Input mismatch: {inp}")
    inp.write_text(canonical_input, encoding="ascii", newline="\n")
    first = execute_stage([GAUSSIAN_BIN / "g16.exe", inp, canonical / "canonical.log"],
                          canonical, "canonical")
    log = canonical / "canonical.log"
    if first["exit_code"] or first["timed_out"] or "Normal termination of Gaussian" not in log.read_text(errors="replace"):
        raise RuntimeError(f"Canonical Gaussian failed: {log}")
    chk = canonical / "canonical.chk"
    if not chk.exists():
        chk = canonical / "canonical.chk"  # explicit failure below
    if not chk.exists():
        raise RuntimeError("Canonical checkpoint missing")
    conversion = execute_stage([GAUSSIAN_BIN / "formchk.exe", chk,
                                canonical / "canonical.fchk"], canonical, "formchk", 300)
    if conversion["exit_code"] or conversion["timed_out"] or not (canonical / "canonical.fchk").exists():
        raise RuntimeError("Canonical formchk failed")
    # A separate checkpoint protects the canonical SCF state from NBO link writes.
    nbo_chk = analysis / "analysis.chk"
    if not nbo_chk.exists():
        shutil.copy2(chk, nbo_chk)
    nbo_options = "PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT"
    nbo_basis = "ChkBasis" if ecp else basis
    nbo_input = (f"%nprocshared=3\n%mem=12GB\n%chk=analysis.chk\n"
                 f"#p {method}/{nbo_basis} Geom=AllCheck Guess=Read "
                 "SCF=(XQC,Tight) Integral=UltraFine NoSymm Pop=NBO6Read\n\n"
                 f"$NBO {nbo_options} $END\n\n")
    nbo_inp = analysis / "analysis.gjf"
    if nbo_inp.exists() and nbo_inp.read_text(encoding="ascii") != nbo_input:
        raise RuntimeError("NBO input mismatch")
    nbo_inp.write_text(nbo_input, encoding="ascii", newline="\n")
    second = execute_stage([GAUSSIAN_BIN / "g16.exe", nbo_inp, analysis / "analysis.log"],
                           analysis, "analysis")
    nbo_log = analysis / "analysis.log"
    normal = nbo_log.exists() and "Normal termination of Gaussian" in nbo_log.read_text(errors="replace")
    atomic_json(target / "manifest.json", {
        "case_id": case_id, "title": title, "charge": charge, "multiplicity": mult,
        "atoms_cartesian_angstrom": atoms, "geometry_source": geometry_source,
        "method": method, "basis": basis, "ecp": ecp,
        "nbo_options": nbo_options, "input_sha256": sha256(inp),
        "nbo_input_sha256": sha256(nbo_inp),
        "source_sha256": sha256(Path(__file__)),
        "g16_sha256": sha256(GAUSSIAN_BIN / "g16.exe"),
        "g16nbo_sha256": sha256(NBO_BIN / "g16nbo.i8.exe"),
        "nbo7_sha256": sha256(NBO_BIN / "nbo7.i8.exe"),
        "run_id": f"nbo-compute-{time.strftime('%Y%m%d')}",
        "steps": {"canonical_scf": first["exit_code"] == 0,
                  "canonical_formchk": conversion["exit_code"] == 0,
                  "nbo_link": second["exit_code"] == 0 and normal},
        "scope": "integration_check_only; no ground-state or scientific validation claim",
    })
    print(json.dumps({"target": str(target), "nbo_exit": second["exit_code"],
                      "normal": normal, "files": [p.name for p in analysis.iterdir()]}))
    return 0 if second["exit_code"] == 0 and normal else 1


def gennbo_reanalyse(number: str) -> int:
    case_id = "NBO-" + number
    analysis_name = "analysis-03" if number == "10" else "analysis-02"
    source = ROOT / "nbo-cases" / case_id / analysis_name / "FILE.47"
    if not source.exists():
        raise FileNotFoundError(source)
    directory = ROOT / "nbo-cases" / case_id / "gennbo-reanalysis"
    directory.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding="ascii")
    options = "PRINT=3 E2PERT=0.0 AONBO=W37 NBOMO=W49 SAO=W50 PLOT"
    modified, count = re.subn(r"\$NBO\s+\$END", f"$NBO {options} $END", text, count=1)
    if count != 1:
        raise RuntimeError("Archive does not have the expected empty $NBO keylist")
    inp = directory / "reanalysis.47"
    if inp.exists() and inp.read_text(encoding="ascii") != modified:
        raise RuntimeError("Reanalysis input differs")
    inp.write_text(modified, encoding="ascii", newline="\n")
    result_path = directory / "process.json"
    if result_path.exists():
        result = json.loads(result_path.read_text(encoding="utf-8"))
    else:
        env = gauss_env(directory)
        env.update(NBOEXE=str(NBO_BIN / "nbo7.i8.exe"), NBOMEM="512mb")
        result = run_tree([NBO_BIN / "gennbo.i8.exe", inp.name], directory, env,
                          sum(physical_core_masks(3)), 20, 300, affinity=True,
                          on_started=lambda data: atomic_json(directory / "live.json", data))
        atomic_json(result_path, result)
        shutil.copy2(directory / "launcher.log", directory / "reanalysis.nbo")
    output = directory / "reanalysis.nbo"
    success = (result["exit_code"] == 0 and not result["timed_out"]
               and output.exists() and "NBO 7.0" in output.read_text(errors="replace"))
    atomic_json(directory / "manifest.json", {
        "source_archive": str(source), "source_archive_sha256": sha256(source),
        "input_sha256": sha256(inp), "options": options,
        "result": "completed" if success else "failed",
        "output_sha256": sha256(output) if output.exists() else None,
    })
    print(json.dumps({"case": case_id, "success": success,
                      "output_size": output.stat().st_size if output.exists() else None}))
    return 0 if success else 1


def audit_case(number: str) -> dict:
    case_id = "NBO-" + number
    target = ROOT / "nbo-cases" / case_id
    analysis = target / ("analysis-03" if number == "10" else "analysis-02")
    canonical = target / "canonical"
    reanalysis = target / "gennbo-reanalysis"
    manifest_path = target / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("scientific_pass", None)
    expected_files = {
        "canonical_gjf": canonical / "canonical.gjf",
        "canonical_log": canonical / "canonical.log",
        "canonical_chk": canonical / "canonical.chk",
        "canonical_fchk": canonical / "canonical.fchk",
        "analysis_gjf": analysis / "analysis.gjf",
        "analysis_log": analysis / "analysis.log",
        "analysis_chk": analysis / "analysis.chk",
        "file47": analysis / "FILE.47",
        "aonbo37": analysis / "FILE.37",
        "nbomo49": analysis / "FILE.49",
        "sao50": analysis / "FILE.50",
        "basis31": analysis / "FILE.31",
        "gennbo_input": reanalysis / "reanalysis.47",
        "gennbo_output": reanalysis / "reanalysis.nbo",
    }
    files = {key: {"path": str(path), "size": path.stat().st_size,
                   "sha256": sha256(path)} if path.is_file() and path.stat().st_size else None
             for key, path in expected_files.items()}
    canonical_log = expected_files["canonical_log"].read_text(errors="replace")
    analysis_log = expected_files["analysis_log"].read_text(errors="replace")
    nbo_text = expected_files["gennbo_output"].read_text(errors="replace")
    fchk = expected_files["canonical_fchk"].read_text(errors="replace")
    basis_match = re.search(r"^Number of basis functions\s+I\s+(\d+)", fchk, re.M)
    nbasis = int(basis_match.group(1)) if basis_match else None
    checks = {
        "all_required_files_nonempty": all(files.values()),
        "canonical_scf_done": "SCF Done:" in canonical_log,
        "canonical_normal": "Normal termination of Gaussian" in canonical_log,
        "analysis_scf_done": "SCF Done:" in analysis_log,
        "analysis_normal": "Normal termination of Gaussian" in analysis_log,
        "gaussian_a03": "EM64W-G16RevA.03" in analysis_log,
        "nbo_7_0_10_gaussian": "NBO 7.0.10 (8-Feb-2021)" in analysis_log,
        "nbo_7_0_10_gennbo": "NBO 7.0.10 (8-Feb-2021)" in nbo_text,
        "npa": "Summary of Natural Population Analysis:" in analysis_log,
        "nbo_orbitals": "NATURAL BOND ORBITAL ANALYSIS" in analysis_log,
        "wiberg": "Wiberg bond index matrix" in analysis_log,
        "e2": "SECOND ORDER PERTURBATION THEORY ANALYSIS" in analysis_log,
        "basis_count_available": nbasis is not None,
        "archive_basis_count": bool(nbasis and
            f"NBAS={nbasis}" in expected_files["file47"].read_text(errors="replace")[:200]),
    }
    if number == "07":
        checks["spin_diagnostic"] = "S**2 before annihilation" in canonical_log
        checks["alpha_beta_nbo"] = ("ALPHA SPIN" in expected_files["aonbo37"].read_text(errors="replace")
                                    and "BETA  SPIN" in expected_files["aonbo37"].read_text(errors="replace"))
    if number == "10":
        checks["ecp_active"] = "Pseudopotential Parameters" in canonical_log and "ECPInt:" in analysis_log
    process_files = [canonical / "canonical-process.json", canonical / "formchk-process.json",
                     analysis / "analysis-process.json", reanalysis / "process.json"]
    processes = [json.loads(path.read_text(encoding="utf-8")) for path in process_files]
    checks["processes_exited_normally"] = all(p["exit_code"] == 0 and not p["timed_out"] for p in processes)
    checks["memory_limit_20_gib_each"] = all(p["memory_limit_gib"] == 20 and
                                           p["peak_tree_commit_bytes"] <= 20 * 1024**3 for p in processes)
    checks["cpu_cap_3_core_equivalent"] = all(p["core_count"] == 3 for p in processes)
    report = {"case_id": case_id, "nbasis": nbasis, "checks": checks,
              "producer_pipeline_complete": all(checks.values()),
              "files": files,
              "peak_tree_commit_bytes": [p["peak_tree_commit_bytes"] for p in processes],
              "process_tree_counts": [p["tree_process_count"] for p in processes],
              "scope": "producer_pipeline",
              "notes": "Artifact presence and producer logs only; matrix content and COV consumer acceptance are separate."}
    atomic_json(target / "audit.json", report)
    manifest["audit_path"] = str(target / "audit.json")
    manifest.pop("integration_complete", None)
    geometry_source = Path(manifest["geometry_source"])
    manifest["geometry_source_sha256"] = sha256(geometry_source) if geometry_source.is_file() else None
    manifest["geometry_source_kind"] = "vendor_file" if geometry_source.is_file() else "explicit_constructed_geometry"
    manifest["producer_versions"] = {"gaussian": "EM64W-G16RevA.03",
                                     "nbo": "NBO 7.0.10 (8-Feb-2021)"}
    spin_match = re.search(r"S\*\*2 before annihilation\s+([0-9.]+),\s+after\s+([0-9.]+)",
                           canonical_log)
    manifest["spin_diagnostic"] = ({"s2_before_annihilation": float(spin_match.group(1)),
                                    "s2_after_annihilation": float(spin_match.group(2))}
                                   if spin_match else {"status": "closed_shell_not_applicable"})
    manifest["producer_pipeline_complete"] = report["producer_pipeline_complete"]
    manifest["scope"] = "producer_pipeline_only; no COV consumer, ground-state, or scientific validation claim"
    atomic_json(manifest_path, manifest)
    return report


def paired_formchk(number: str) -> dict:
    target = ROOT / "nbo-cases" / ("NBO-" + number)
    analysis = target / ("analysis-03" if number == "10" else "analysis-02")
    chk = analysis / "analysis.chk"
    if not chk.is_file():
        raise FileNotFoundError(chk)
    fchk = analysis / "paired-canonical.fchk"
    result = execute_stage([GAUSSIAN_BIN / "formchk.exe", chk, fchk],
                           analysis, "paired-formchk", 300)
    completed = (result["exit_code"] == 0 and not result["timed_out"]
                 and fchk.is_file() and fchk.stat().st_size > 0)
    manifest_path = target / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["paired_fchk"] = {
        "path": str(fchk), "source_checkpoint": str(chk),
        "source_checkpoint_sha256": sha256(chk),
        "sha256": sha256(fchk) if completed else None,
        "status": "completed" if completed else "failed",
        "relationship": "Formatted from the same Gaussian/NBO analysis checkpoint; original canonical FCHK remains preserved",
    }
    atomic_json(manifest_path, manifest)
    return {"case_id": target.name, "completed": completed,
            "fchk": str(fchk), "bytes": fchk.stat().st_size if fchk.exists() else None}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command")
    args = parser.parse_args()
    if args.command == "gennbo-smoke":
        return gennbo_smoke()
    if args.command.startswith("gaussian-"):
        return gaussian_smoke(args.command.removeprefix("gaussian-"))
    if args.command.startswith("reanalyse-"):
        return gennbo_reanalyse(args.command.removeprefix("reanalyse-"))
    if args.command == "audit":
        reports = [audit_case(f"{i:02d}") for i in range(1, 11)]
        print(json.dumps({"complete": sum(r["producer_pipeline_complete"] for r in reports),
                          "cases": {r["case_id"]: [key for key, ok in r["checks"].items() if not ok]
                                    for r in reports}}))
        return 0 if all(r["producer_pipeline_complete"] for r in reports) else 1
    if args.command == "paired-formchk":
        reports = [paired_formchk(f"{i:02d}") for i in range(1, 11)]
        print(json.dumps(reports))
        return 0 if all(r["completed"] for r in reports) else 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
