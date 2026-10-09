"""Independent raw-file algebra checks for the focused NBO matrix producer.

The input is the vendor archive, W-format matrices and the paired Gaussian
FCHK. No COV parser or exported COV representation is imported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np

from nbo_reference import (_archive_header, _archive_lcaomo, _archive_packed,
                           _fchk_electron_count, _fchk_packed, _read_sao,
                           _read_w_square, load_one)


WORK = Path(__file__).resolve().parents[2]
FOCUS = WORK / "nbo-focus-cases"
ORIGINAL = WORK / "nbo-cases"


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def measure(observed: float, tolerance: float, note: str = "") -> dict:
    return {"max_abs_residual": float(observed), "tolerance": tolerance,
            "status": "pass" if np.isfinite(observed) and observed <= tolerance else "fail",
            "note": note}


def _read_matrix(directory: Path, index: int, n: int, heading: str, unrestricted: bool):
    return _read_w_square(directory / f"FILE.{index}", n, heading, unrestricted)


def _fchk_coefficients(fchk: Path, n: int) -> dict[str, np.ndarray]:
    mol = load_one(str(fchk))
    if mol.mo is None:
        return {}
    coeff = np.asarray(mol.mo.coeffs, dtype=float)
    if coeff.shape[0] != n:
        raise ValueError(f"FCHK_AO_dimension_{coeff.shape[0]}_vs_{n}")
    result = {"alpha": coeff[:, :n]}
    if coeff.shape[1] >= 2 * n:
        result["beta"] = coeff[:, n:2*n]
    return result


def _fchk_integer_array(path: Path, label: str) -> list[int]:
    lines = path.read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if line.startswith(label):
            match = re.search(r"\bI\s+N=\s*(\d+)", line)
            if not match:
                raise ValueError(f"invalid_FCHK_integer_array:{label}")
            count = int(match.group(1))
            values: list[int] = []
            for row in lines[i+1:]:
                if not re.fullmatch(r"\s*[+\-\d\s]+", row):
                    break
                values.extend(int(word) for word in row.split())
                if len(values) >= count:
                    break
            if len(values) != count:
                raise ValueError(f"FCHK_{label}_expected_{count}_found_{len(values)}")
            return values
    raise ValueError(f"missing_FCHK_array:{label}")


def _fchk_scalar_integer(path: Path, label: str) -> int:
    match = re.search(rf"^{re.escape(label)}\s+I\s+(-?\d+)\s*$",
                      path.read_text(errors="replace"), re.M)
    if not match:
        raise ValueError(f"missing_FCHK_scalar:{label}")
    return int(match.group(1))


def _printed_population(report: Path) -> dict:
    text = report.read_text(errors="replace")
    # The open-shell total table appends a Spin column after its population.
    match = re.search(r"^\s*\* Total \*\s+[-+\d.]+\s+[-+\d.]+\s+[-+\d.]+\s+[-+\d.]+\s+([-+\d.]+)(?:\s+[-+\d.]+)?\s*$", text, re.M)
    ecp = re.search(r"\[(\d+) electrons found in the effective core potential\]", text)
    return {"printed_npa_total": float(match.group(1)) if match else None,
            "printed_ecp_core_count": int(ecp.group(1)) if ecp else 0}


def _printed_nao_labels(report: Path, nbasis: int) -> list[dict]:
    rows = []
    row_re = re.compile(r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(\d+)\s+(\S+)\s+"
                        r"(Cor|Val|Ryd)\(\s*(\d+)([spdfgh])\)\s+([+\-]?\d+(?:\.\d+)?)")
    for line in report.read_text(errors="replace").splitlines():
        match = row_re.match(line)
        if match:
            index, symbol, atom, angular, category, shell_n, shell_l, occupation = match.groups()
            rows.append({"index": int(index), "symbol": symbol, "atom": int(atom),
                         "angular": angular, "category": category,
                         "principal_n": int(shell_n), "shell_l": shell_l,
                         "printed_occupation": float(occupation), "raw": line})
            if len(rows) == nbasis:
                break
    return rows


def check(case: str) -> dict:
    case_id = f"NBO-{case}" if case in ("01", "07") else case
    meta = json.loads((FOCUS / case_id / "case.json").read_text(encoding="utf-8"))
    for key in ("fchk", "report", "archive47", "aonbo", "nbomo", "naomo", "aonao", "naonbo", "sao"):
        if meta.get(key) and not Path(meta[key]).is_absolute():
            meta[key] = str((FOCUS / case_id / meta[key]).resolve())
    directory = Path(meta["aonbo"]).parent
    fchk = Path(meta["fchk"])
    archive = Path(meta["archive47"])
    header_match = re.search(r"NBAS=(\d+)", archive.read_text(encoding="ascii")[:160])
    if not header_match:
        raise ValueError("missing_archive_NBAS")
    n = int(header_match.group(1))
    fchk_basis_count = _fchk_scalar_integer(fchk, "Number of basis functions")
    independent_count = _fchk_scalar_integer(fchk, "Number of independent functions")
    if fchk_basis_count != n:
        raise ValueError(f"archive_NBAS_{n}_vs_FCHK_basis_{fchk_basis_count}")
    if independent_count < n:
        raise ValueError(f"NBAS_gt_independent_MOs_{n}_vs_{independent_count}:rectangular_case_not_covered")
    header = _archive_header(archive, n)
    unrestricted = header["open_shell"]
    spins = ("alpha", "beta") if unrestricted else ("alpha",)
    # NBO writes one common AO->NAO basis for an unrestricted calculation;
    # spin-dependent transformations (W51/W53) carry explicit spin blocks.
    a_one = _read_matrix(directory, 52, n, "NAOs in the AO basis:", False)["alpha"]
    a = {spin: a_one for spin in spins}
    t = _read_matrix(directory, 51, n, "MOs in the NAO basis:", unrestricted)
    u = _read_matrix(directory, 53, n, "NBOs in the NAO basis:", unrestricted)
    b = _read_matrix(directory, 37, n, "NBOs in the AO basis:", unrestricted)
    v = _read_matrix(directory, 49, n, "MOs in the NBO basis:", unrestricted)
    s = _read_sao(directory / "FILE.50", n)
    c_arch = _archive_lcaomo(archive, n, unrestricted)
    c_fchk = _fchk_coefficients(fchk, n)
    original_fchk = ((ORIGINAL if case in ("01", "07") else FOCUS)
                     / case_id / "canonical" / "canonical.fchk")
    c_original = _fchk_coefficients(original_fchk, n) if original_fchk.is_file() else {}
    dens_arch = _archive_packed(archive, "DENSITY", n, len(spins))
    density = dict(zip(spins, dens_arch))
    report = {"case": case_id, "analysis": str(directory),
              "fchk": str(fchk), "archive": str(archive), "nbasis": n,
              "fchk_basis_count": fchk_basis_count,
              "fchk_independent_function_count": independent_count,
              "dimension_scope": "square_full_rank_case",
              "spins": list(spins), "dimensions": {}, "checks": {},
              "fchk_available_coefficients": sorted(c_fchk),
              "original_fchk": str(original_fchk),
              "original_fchk_column_comparison": {}}
    full_electrons = _fchk_scalar_integer(fchk, "Number of electrons")
    charge = _fchk_scalar_integer(fchk, "Charge")
    atomic_numbers = _fchk_integer_array(fchk, "Atomic numbers")
    z_sum = sum(atomic_numbers)
    ecp_core = z_sum - charge - full_electrons
    printed = _printed_population(Path(meta["report"]))
    report["electron_bookkeeping"] = {
        "fchk_electrons_explicit": full_electrons,
        "sum_atomic_numbers": z_sum, "charge": charge,
        "ecp_core_derived": ecp_core,
        "npa_printed_total": printed["printed_npa_total"],
        "nbo_printed_ecp_core": printed["printed_ecp_core_count"]}
    nao_labels = _printed_nao_labels(Path(meta["report"]), n)
    report["nao_labels"] = nao_labels
    report["nao_label_counts"] = {
        "rows": len(nao_labels), "zero_printed_occupation_rows": sum(
            row["printed_occupation"] == 0 for row in nao_labels),
        "categories": {category: sum(row["category"] == category for row in nao_labels)
                       for category in ("Cor", "Val", "Ryd")}}
    report["numerical_policy"] = {
        "W_format_decimal_places": 9,
        "matrix_absolute_tolerance": 1e-6,
        "density_absolute_tolerance": 1e-6,
        "selected_FCHK_phase_aligned_column_absolute_tolerance": 1e-5,
        "selected_FCHK_metric_tolerance": 1e-5,
        "NPA_printed_population_tolerance": 5e-4,
        "rationale": "W files print 9 decimal places; at NBAS 134, products sum up to 134 rounded terms. FCHK columns also round independently. Limits are above these rounding accumulations and below a chemically meaningful matrix discrepancy; no case-specific relaxation."}
    checks = report["checks"]
    checks["NAO_label_count"] = measure(abs(len(nao_labels) - n), 0)
    checks["NAO_label_indices"] = measure(
        0 if [row["index"] for row in nao_labels] == list(range(1, n+1)) else 1, 0)
    checks["ECP_core_vs_printed"] = measure(abs(ecp_core - printed["printed_ecp_core_count"]), 0)
    if printed["printed_npa_total"] is not None:
        checks["NPA_total_vs_all_electrons"] = measure(abs(printed["printed_npa_total"] - (z_sum-charge)), 5e-4)
        checks["NPA_minus_ECP_vs_explicit_electrons"] = measure(
            abs(printed["printed_npa_total"] - ecp_core - full_electrons), 5e-4)
    identity = np.eye(n)
    report["checks"]["overlap_positive_min_eigenvalue"] = {
        "value": float(np.linalg.eigvalsh(s)[0]),
        "status": "pass" if np.linalg.eigvalsh(s)[0] > 0 else "fail"}
    for spin in spins:
        report["dimensions"][spin] = {"AO_NAO": list(a[spin].shape),
                                      "NAO_MO": list(t[spin].shape),
                                      "NAO_NBO": list(u[spin].shape),
                                      "AO_NBO": list(b[spin].shape),
                                      "NBO_MO": list(v[spin].shape)}
        p = f"{spin}."
        checks = report["checks"]
        checks[p+"A_transpose_S_A"] = measure(max_abs(a[spin].T @ s @ a[spin] - identity), 1e-6)
        checks[p+"AO_NAO_times_NAO_MO_vs_archive_MO"] = measure(max_abs(a[spin] @ t[spin] - c_arch[spin]), 1e-6)
        checks[p+"AO_NAO_times_NAO_NBO_vs_AO_NBO"] = measure(max_abs(a[spin] @ u[spin] - b[spin]), 1e-6)
        checks[p+"NAO_NBO_times_NBO_MO_vs_NAO_MO"] = measure(max_abs(u[spin] @ v[spin] - t[spin]), 1e-6)
        checks[p+"AO_NBO_times_NBO_MO_vs_archive_MO"] = measure(max_abs(b[spin] @ v[spin] - c_arch[spin]), 1e-6)
        checks[p+"NAO_NBO_orthogonal"] = measure(max_abs(u[spin].T @ u[spin] - identity), 1e-6)
        checks[p+"NAO_MO_orthogonal"] = measure(max_abs(t[spin].T @ t[spin] - identity), 1e-6)
        checks[p+"AO_NBO_metric"] = measure(max_abs(b[spin].T @ s @ b[spin] - identity), 1e-6)
        checks[p+"NBO_MO_orthogonal"] = measure(max_abs(v[spin].T @ v[spin] - identity), 1e-6)
        checks[p+"NBO_MO_column_weights"] = measure(max_abs(np.sum(v[spin] ** 2, axis=0) - 1), 1e-6)
        checks[p+"NAO_NBO_column_weights"] = measure(max_abs(np.sum(u[spin] ** 2, axis=0) - 1), 1e-6)
        report.setdefault("null_columns", {})[spin] = {
            "AO_NAO": [int(i+1) for i, value in enumerate(np.linalg.norm(a[spin], axis=0)) if value < 1e-10],
            "NAO_MO": [int(i+1) for i, value in enumerate(np.linalg.norm(t[spin], axis=0)) if value < 1e-10],
            "NAO_NBO": [int(i+1) for i, value in enumerate(np.linalg.norm(u[spin], axis=0)) if value < 1e-10]}
        count = _fchk_electron_count(fchk, spin)
        report.setdefault("electron_counts", {})[spin] = count
        occ = np.r_[np.full(count, 1 if unrestricted else 2), np.zeros(n-count)]
        expected_density = (c_arch[spin] * occ) @ c_arch[spin].T
        checks[p+"archive_density_vs_occupied_MO"] = measure(max_abs(expected_density - density[spin]), 1e-6)
        explicit_count = count * (1 if unrestricted else 2)
        checks[p+"density_trace"] = measure(abs(float(np.trace(density[spin] @ s)) - explicit_count), 1e-6)
        if spin in c_fchk:
            selected = c_fchk[spin]
            metric = c_arch[spin].T @ s @ selected
            phases = np.where(np.diag(metric) < 0, -1.0, 1.0)
            aligned_error = max_abs(c_arch[spin] - selected * phases[None, :])
            offdiag_error = max_abs(metric - np.diag(np.diag(metric)))
            diagonal_error = float(np.max(np.abs(np.abs(np.diag(metric)) - 1.0)))
            checks[p+"archive_MO_vs_selected_FCHK_columnwise"] = {
                "status": "pass" if (aligned_error <= 1e-5 and offdiag_error <= 1e-5
                                     and diagonal_error <= 1e-5) else "fail",
                "max_abs_direct_coeff_residual": max_abs(c_arch[spin] - selected),
                "max_abs_phase_aligned_coeff_residual": aligned_error,
                "max_abs_metric_off_diagonal_overlap": offdiag_error,
                "max_abs_metric_diagonal_abs_deviation": diagonal_error,
                "min_abs_metric_diagonal_overlap": float(np.min(np.abs(np.diag(metric)))),
                "column_phase_alignment": [int(item) for item in phases],
                "tolerances": {"coeff": 1e-5, "metric_off_diagonal": 1e-5,
                               "metric_diagonal_abs_deviation": 1e-5},
                "fchk_relationship": meta["fchk_relationship"]}
        else:
            checks[p+"archive_MO_vs_selected_FCHK_columnwise"] = {"status": "unavailable",
                "reason": "selected_FCHK_missing_spin_coefficients"}
        if spin in c_original:
            reference = c_original[spin]
            direct = max_abs(c_arch[spin] - reference)
            overlap_orbitals = c_arch[spin].T @ s @ reference
            signs = np.where(np.diag(overlap_orbitals) < 0, -1.0, 1.0)
            signed = max_abs(c_arch[spin] - reference * signs[None, :])
            report["original_fchk_column_comparison"][spin] = {
                "max_abs_direct_coeff_residual": direct,
                "max_abs_phase_corrected_coeff_residual": signed,
                "column_phase_alignment": [int(item) for item in signs],
                "min_abs_metric_diagonal_overlap": float(np.min(np.abs(np.diag(overlap_orbitals)))),
                "max_abs_metric_off_diagonal_overlap": max_abs(overlap_orbitals - np.diag(np.diag(overlap_orbitals))),
                "columnwise_compatible_tolerance": 1e-5,
                "status": "compatible_after_column_phase" if signed <= 1e-5 else "incompatible"}
    scf = _fchk_packed(fchk, "Total SCF Density", n)
    if scf is not None:
        checks["FCHK_total_density"] = measure(max_abs(scf - sum(density.values())), 1e-6)
    if unrestricted:
        spin_density = _fchk_packed(fchk, "Spin SCF Density", n)
        if spin_density is not None:
            checks["FCHK_spin_density"] = measure(max_abs(spin_density - (density["alpha"] - density["beta"])), 1e-6)
    report["status"] = "pass" if all(item["status"] in ("pass", "unavailable")
                                       for item in checks.values()) else "fail"
    return report


def main() -> int:
    global FOCUS, ORIGINAL
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["01", "07", "NBO-F11", "NBO-F12"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--focus-root", type=Path, default=FOCUS)
    parser.add_argument("--original-root", type=Path, default=ORIGINAL)
    args = parser.parse_args()
    FOCUS, ORIGINAL = args.focus_root.resolve(), args.original_root.resolve()
    report = check(args.case)
    payload = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
