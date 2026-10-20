"""Mutation-based NBO import contracts against the production C++ probe.

Never runs Gaussian, NBO, or the GUI. Inputs are copied into a new evidence
directory; originals are hashed before/after. A passing test requires the
expected diagnostic class AND the expected dataset/render availability.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def sources(root: Path, case: str) -> dict[str, Path]:
    case_dir = root / case
    successes = []
    for directory in sorted(case_dir.glob("analysis-*")):
        log = directory / "analysis.log"
        if log.is_file() and "Normal termination of Gaussian" in log.read_text(errors="replace"):
            if (directory / "paired-canonical.fchk").is_file():
                successes.append(directory)
    if len(successes) != 1:
        raise ValueError(f"Need one successful paired analysis for {case}, found {successes}")
    directory = successes[0]
    result = {"fchk": directory / "paired-canonical.fchk", "log": directory / "analysis.log",
              "archive": directory / "FILE.47", "aonbo": directory / "FILE.37", "nbomo": directory / "FILE.49"}
    for path in result.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    return result


def edit(path: Path, transform) -> None:
    path.write_text(transform(path.read_text(errors="replace")), encoding="utf-8")


def modify_coordinate(text: str) -> str:
    # Only the first atom's x coordinate changes; dimensions, basis and all
    # electronic matrices remain byte-for-byte the same.
    pattern = r"(?m)^(\s*\d+\s+\d+\s+)([+-]?\d+\.\d+)(\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s*)$"
    changed, count = re.subn(pattern, lambda m: m[1] + f"{float(m[2])+0.125:.6f}" + m[3], text, count=1)
    if count != 1:
        raise ValueError("coordinate mutation did not identify exactly one atom")
    return changed


def modify_density(text: str) -> str:
    pattern = r"(\$DENSITY\s+)([+-]?\d+\.\d+[EeDd][+-]\d+)"
    changed, count = re.subn(pattern, lambda m: m[1] + f"{float(m[2].replace('D','E'))+0.05:.12E}", text, count=1)
    if count != 1:
        raise ValueError("density mutation did not identify first density element")
    return changed


def remove_e2(text: str) -> str:
    begin = text.index("SECOND ORDER PERTURBATION")
    begin = text.rfind("\n", 0, begin) + 1
    end = text.index("NATURAL BOND ORBITALS (Summary)", begin)
    end = text.rfind("\n", begin, end) + 1
    return text[:begin] + text[end:]


def remove_nbo_summary(text: str) -> str:
    begin = text.index("NATURAL BOND ORBITALS (Summary)")
    begin = text.rfind("\n", 0, begin) + 1
    end = text.index("Total Lewis", begin)
    end = text.index("\n", end) + 1
    return text[:begin] + text[end:]


def swap_aonbo_columns(text: str, n: int) -> str:
    rows = text.splitlines(keepends=True)
    real_row = re.compile(r"^\s*(?:[+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?\s*)+$")
    indices, values = [], []
    for i, row in enumerate(rows):
        if real_row.fullmatch(row):
            indices.append(i)
            values.extend(row.split())
            if len(values) >= n*n:
                break
    if len(values) != n*n:
        raise ValueError("AONBO payload did not terminate at expected boundary")
    values[:2*n] = values[n:2*n] + values[:n]
    cursor = 0
    for i in indices:
        count = len(rows[i].split())
        rows[i] = " " + " ".join(values[cursor:cursor+count]) + "\n"
        cursor += count
    return "".join(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite evidence: {args.output}")
    exe = args.build / ("cov_nbo_probe.exe" if (args.build / "cov_nbo_probe.exe").is_file() else "cov_nbo_probe")
    if not exe.is_file():
        raise FileNotFoundError(exe)
    base = sources(args.cases, "NBO-02")
    open_base = sources(args.cases, "NBO-07")
    other = sources(args.cases, "NBO-01")
    originals = {str(p.resolve()): digest(p) for group in (base, open_base, other) for p in group.values()}
    args.output.mkdir(parents=True)
    results = []

    def run(name, expected_class, *, origin=base, mutate=None, omit=(), expected_exit=None, verify=None):
        directory = args.output / name
        directory.mkdir()
        local = {}
        for key, path in origin.items():
            local[key] = directory / path.name
            shutil.copyfile(path, local[key])
        if mutate:
            mutate(local)
        mutations = {key: {"path": str(path.resolve()), "sha256": digest(path),
                           "source": str(origin[key].resolve()), "source_sha256": digest(origin[key]),
                           "changed": digest(path) != digest(origin[key])} for key, path in local.items()}
        save(directory / "input-manifest.json", mutations)
        paths = [str(local[key].resolve()) if key not in omit else "-" for key in ("fchk", "log", "archive", "aonbo", "nbomo")]
        evidence_dir = directory / "probe"
        command = [str(exe.resolve()), *paths, str(evidence_dir.resolve())]
        process = subprocess.run(command, capture_output=True, text=True, errors="replace", timeout=120)
        (directory / "stdout.txt").write_text(process.stdout, encoding="utf-8")
        (directory / "stderr.txt").write_text(process.stderr, encoding="utf-8")
        dataset_path = evidence_dir / "dataset.json"
        dataset = json.loads(dataset_path.read_text(encoding="utf-8")) if dataset_path.is_file() else None
        rendered = (evidence_dir / "render-contract.json").is_file()
        diagnostic = process.stderr.strip()
        observed_class = dataset.get("association", {}).get("status") if dataset and process.returncode == 2 else diagnostic.split(":", 1)[0]
        assertions = {"expected_error_category": any(token in diagnostic for token in expected_class) if expected_class else process.returncode == 0,
                      "expected_exit": process.returncode in expected_exit if expected_exit is not None else True,
                      "no_render_on_rejection": not rendered if process.returncode != 0 else True}
        if verify:
            assertions.update(verify(dataset, rendered))
        result = {"name": name, "passed": all(assertions.values()), "expected_categories": list(expected_class),
                  "observed_category": observed_class, "exit_code": process.returncode, "diagnostic": diagnostic,
                  "dataset_available": dataset is not None, "render_available": rendered,
                  "assertions": assertions, "command": command, "manifest": str((directory / "input-manifest.json").resolve())}
        save(directory / "result.json", result)
        results.append(result)

    run("wrong_molecule_sidecars", ("incompatible_dimensions",),
        mutate=lambda p: [shutil.copyfile(other[k], p[k]) for k in ("log", "archive", "aonbo", "nbomo")], expected_exit=(2,))
    run("same_dimensions_changed_geometry", ("incompatible_geometry",),
        mutate=lambda p: edit(p["archive"], modify_coordinate), expected_exit=(2,))
    nbas = int(re.search(r"NBAS\s*=\s*(\d+)", base["archive"].read_text())[1])
    run("same_basis_wrong_orbital_columns", ("incompatible_nbomo", "incompatible_nbo_output"),
        mutate=lambda p: edit(p["aonbo"], lambda s: swap_aonbo_columns(s, nbas)), expected_exit=(2,))
    run("truncated_archive", ("Unterminated .47 block", ".47 matrix size mismatch"),
        mutate=lambda p: edit(p["archive"], lambda s: s[:len(s)//2]), expected_exit=(1,))
    for key in ("aonbo", "nbomo"):
        run("truncated_"+key, ("incomplete matrix payload",),
            mutate=lambda p, key=key: edit(p[key], lambda s: "\n".join(s.splitlines()[:8])+"\n"), expected_exit=(1,))
    run("wrong_spin_block_order", ("OPEN matrix requires exactly one alpha and beta block",), origin=open_base,
        mutate=lambda p: edit(p["aonbo"], lambda s: s.replace("ALPHA SPIN", "TEMP SPIN").replace("BETA  SPIN", "ALPHA SPIN").replace("TEMP SPIN", "BETA  SPIN")), expected_exit=(1,))
    run("same_basis_changed_density", ("incompatible_density",),
        mutate=lambda p: edit(p["archive"], modify_density), expected_exit=(2,))
    run("ambiguous_duplicate_analysis", ("Multiple NBO analyses",),
        mutate=lambda p: edit(p["log"], lambda s: s+"\n"+s), expected_exit=(1,))
    run("missing_e2_is_unknown", (), mutate=lambda p: edit(p["log"], remove_e2), expected_exit=(0,),
        verify=lambda d,r: {"dataset_readable": d is not None, "render_still_available": r,
                            "e2_absent": d is not None and not d["e2"],
                            "threshold_unknown": d is not None and all(x["printing_threshold"] is None for x in d["e2_sections"]),
                            "missing_reason_explicit": d is not None and all(x["missing_reason"] for x in d["e2_sections"])})
    run("missing_matrices_readable_not_renderable", ("Complete AONBO matrix",), omit=("aonbo", "nbomo"), expected_exit=(1,),
        verify=lambda d,r: {"dataset_readable": d is not None and bool(d["orbitals"]), "matrix_collection_empty": d is not None and not d["matrices"], "not_renderable": not r})
    run("nbomo_without_aonbo", ("missing_aonbo",), omit=("aonbo",), expected_exit=(2,))
    run("missing_archive_not_associated", ("missing_archive",), omit=("archive", "aonbo", "nbomo"), expected_exit=(2,),
        verify=lambda d,r: {"dataset_readable": d is not None and bool(d["orbitals"]), "association_rejected": d is not None and not d["association"]["compatible"]})
    cmo = "\n CMO: NBO Analysis of Canonical Molecular Orbitals\n Leading (> 5%) NBO Contributions to Molecular Orbitals\n MO 1 (occ): orbital energy = -1.0 a.u.\n 0.950*[ 1]: CR ( 1) C 1\n"
    run("thresholded_cmo_not_complete_matrix", ("Complete AONBO matrix",), omit=("aonbo", "nbomo"),
        mutate=lambda p: edit(p["log"], lambda s: s+cmo), expected_exit=(1,),
        verify=lambda d,r: {"cmo_summary_retained": d is not None and bool(d["cmo_summaries"]), "no_full_matrix_invented": d is not None and not d["matrices"], "threshold_warning": d is not None and any("thresholded" in w for w in d["warnings"]), "not_renderable": not r})
    run("printed_e2_threshold_preserved", (),
        mutate=lambda p: edit(p["log"], lambda s: re.sub(r"(Threshold for printing:\s*)0\.00(\s*kcal/mol)", r"\g<1>0.50\2", s)), expected_exit=(0,),
        verify=lambda d,r: {"threshold_is_half": d is not None and all(x["printing_threshold"] == 0.5 for x in d["e2_sections"]), "omissions_not_zero": d is not None and all("never" in x["missing_reason"].lower() for x in d["e2_sections"])})
    run("missing_nbo_energy_stays_null", (), mutate=lambda p: edit(p["log"], remove_nbo_summary), expected_exit=(0,),
        verify=lambda d,r: {"missing_energy_not_zero": d is not None and bool(d["orbitals"]) and all(x["diagonal_fock_hartree"] is None for x in d["orbitals"]),
                            "missing_energy_source_null": d is not None and all(x["energy_source"] is None for x in d["orbitals"]),
                            "render_independent_of_energy": r})
    preserved = all(Path(path).is_file() and digest(Path(path)) == sha for path, sha in originals.items())
    summary = {"schema": "cov.nbo.negative.v1", "passed": preserved and all(x["passed"] for x in results),
               "passed_count": sum(x["passed"] for x in results), "test_count": len(results),
               "originals_unchanged": preserved, "original_hashes": originals,
               "probe": {"path": str(exe.resolve()), "sha256": digest(exe)}, "tests": results,
               "scope": "Import/rejection and missing-data contracts only; no new chemistry, GUI, or scientific acceptance"}
    save(args.output / "summary.json", summary)
    print(json.dumps({k: summary[k] for k in ("passed", "passed_count", "test_count", "originals_unchanged")}))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Negative test harness failed: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(2)
