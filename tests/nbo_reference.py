"""Independent, read-only checks of raw Gaussian/NBO orbital representations.

Reads vendor W-format FILE.37/49/50 and FILE.47 directly. No COV parser or
exported render coefficients are used. The documented W-format matrix is the
first NBAS*NBAS values of each spin block; labels, occupations and plot data
that follow it are retained in the producer file but never parsed as matrix.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import os
import re
import sys

DEPS = Path(os.environ.get("COV_REFERENCE_DEPS",
                           r"E:\Dev\cov-validation-20260905\reference-deps"))
if not DEPS.is_dir():
    raise RuntimeError(f"Reference dependencies unavailable: {DEPS}")
sys.path.insert(0, str(DEPS))
import numpy as np
from gbasis.evals.eval import evaluate_basis
from gbasis.wrappers import from_iodata
from iodata import load_one


class ReferenceInputError(ValueError):
    """A required source is absent or cannot be interpreted without guessing."""


def _require_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise ReferenceInputError(f"missing_or_empty:{path}")
    return path


def _float_tokens(line: str) -> list[float] | None:
    words = line.split()
    if not words:
        return None
    if not all(re.fullmatch(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[EeDd][+-]?\d+)?", w) for w in words):
        return None
    return [float(w.replace("D", "E").replace("d", "e")) for w in words]


def _take_numeric_block(lines: list[str], start: int, count: int, label: str) -> np.ndarray:
    values: list[float] = []
    for line in lines[start:]:
        row = _float_tokens(line)
        if row is None:
            break
        values.extend(row)
        if len(values) >= count:
            break
    if len(values) != count:
        raise ReferenceInputError(f"{label}:expected_{count}_matrix_values_found_{len(values)}")
    result = np.asarray(values, dtype=np.float64)
    if not np.isfinite(result).all():
        raise ReferenceInputError(f"{label}:nonfinite_matrix_value")
    return result


def _read_w_square(path: Path, nbasis: int, heading: str, unrestricted: bool) -> dict[str, np.ndarray]:
    lines = _require_file(path).read_text(encoding="ascii", errors="strict").splitlines()
    if len(lines) < 4 or lines[1].strip() != heading or not lines[2].strip().startswith("---"):
        raise ReferenceInputError(f"{path}:unexpected_W_format_heading")
    if unrestricted:
        alpha = [i for i, line in enumerate(lines) if line.strip() == "ALPHA SPIN"]
        beta = [i for i, line in enumerate(lines) if line.strip() == "BETA  SPIN"]
        if len(alpha) != 1 or len(beta) != 1 or alpha[0] >= beta[0]:
            raise ReferenceInputError(f"{path}:expected_one_alpha_and_beta_block")
        starts = {"alpha": alpha[0] + 1, "beta": beta[0] + 1}
    else:
        if any("SPIN" in line for line in lines[:5]):
            raise ReferenceInputError(f"{path}:unexpected_spin_block")
        starts = {"alpha": 3}
    return {spin: _take_numeric_block(lines, start, nbasis * nbasis,
                                     f"{path.name}:{spin}").reshape((nbasis, nbasis), order="F")
            for spin, start in starts.items()}


def _read_sao(path: Path, nbasis: int) -> np.ndarray:
    lines = _require_file(path).read_text(encoding="ascii", errors="strict").splitlines()
    if len(lines) < 4 or lines[1].strip() != "AO overlap matrix:" or not lines[2].strip().startswith("---"):
        raise ReferenceInputError(f"{path}:unexpected_SAO_heading")
    packed = _take_numeric_block(lines, 3, nbasis * (nbasis + 1) // 2, path.name)
    overlap = np.zeros((nbasis, nbasis), dtype=np.float64)
    overlap[np.tril_indices(nbasis)] = packed
    overlap += np.tril(overlap, -1).T
    return overlap


def _archive_header(path: Path, nbasis: int) -> dict:
    header = _require_file(path).read_text(encoding="ascii", errors="strict")[:400]
    match = re.search(r"\$GENNBO\s+NATOMS=(\d+)\s+NBAS=(\d+)\s+([^$]+)\$END", header)
    if not match:
        raise ReferenceInputError(f"{path}:missing_GENNBO_header")
    if int(match.group(2)) != nbasis:
        raise ReferenceInputError(f"{path}:NBAS_mismatch_{match.group(2)}_vs_{nbasis}")
    is_open = "OPEN" in match.group(3).split()
    return {"natoms": int(match.group(1)), "nbasis": nbasis,
            "open_shell": is_open, "flags": match.group(3).split()}


def _archive_lcaomo(path: Path, nbasis: int, unrestricted: bool) -> dict[str, np.ndarray]:
    text = _require_file(path).read_text(encoding="ascii", errors="strict")
    blocks = re.findall(r"\$LCAOMO\s*(.*?)\$END", text, flags=re.S)
    if len(blocks) != 1:
        raise ReferenceInputError(f"{path}:expected_one_LCAOMO_datalist")
    words = blocks[0].split()
    expected = nbasis * nbasis * (2 if unrestricted else 1)
    if len(words) != expected:
        raise ReferenceInputError(f"{path}:LCAOMO_expected_{expected}_found_{len(words)}")
    try:
        values = np.asarray([float(word.replace("D", "E")) for word in words], dtype=np.float64)
    except ValueError as error:
        raise ReferenceInputError(f"{path}:LCAOMO_non_numeric:{error}") from error
    if not np.isfinite(values).all():
        raise ReferenceInputError(f"{path}:LCAOMO_nonfinite")
    result = {"alpha": values[:nbasis * nbasis].reshape((nbasis, nbasis), order="F")}
    if unrestricted:
        result["beta"] = values[nbasis * nbasis:].reshape((nbasis, nbasis), order="F")
    return result


def _fchk_electron_count(path: Path, spin: str) -> int:
    label = "Number of alpha electrons" if spin == "alpha" else "Number of beta electrons"
    match = re.search(rf"^{label}\s+I\s+(\d+)\s*$",
                      _require_file(path).read_text(errors="replace"), re.M)
    if not match:
        raise ReferenceInputError(f"{path}:missing_{label.replace(' ', '_')}")
    return int(match.group(1))


def _packed_to_symmetric(values: np.ndarray, nbasis: int) -> np.ndarray:
    matrix = np.zeros((nbasis, nbasis), dtype=np.float64)
    matrix[np.tril_indices(nbasis)] = values
    matrix += np.tril(matrix, -1).T
    return matrix


def _archive_packed(path: Path, name: str, nbasis: int, count: int = 1) -> list[np.ndarray]:
    text = _require_file(path).read_text(encoding="ascii", errors="strict")
    blocks = re.findall(rf"\${re.escape(name)}\s*(.*?)\$END", text, flags=re.S)
    if len(blocks) != 1:
        raise ReferenceInputError(f"{path}:expected_one_{name}_datalist")
    words = blocks[0].split()
    packed_count = nbasis * (nbasis + 1) // 2
    if len(words) != count * packed_count:
        raise ReferenceInputError(f"{path}:{name}_expected_{count * packed_count}_found_{len(words)}")
    try:
        values = np.asarray([float(word.replace("D", "E")) for word in words], dtype=np.float64)
    except ValueError as error:
        raise ReferenceInputError(f"{path}:{name}_non_numeric:{error}") from error
    if not np.isfinite(values).all():
        raise ReferenceInputError(f"{path}:{name}_nonfinite")
    return [_packed_to_symmetric(values[i * packed_count:(i+1) * packed_count], nbasis)
            for i in range(count)]


def _fchk_packed(path: Path, label: str, nbasis: int) -> np.ndarray | None:
    lines = _require_file(path).read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if line.startswith(label):
            match = re.search(r"\bR\s+N=\s*(\d+)\s*$", line)
            if not match:
                raise ReferenceInputError(f"{path}:{label}_invalid_header")
            count = nbasis * (nbasis + 1) // 2
            if int(match.group(1)) != count:
                raise ReferenceInputError(f"{path}:{label}_expected_{count}_found_{match.group(1)}")
            return _packed_to_symmetric(_take_numeric_block(lines, i + 1, count, label), nbasis)
    return None


def _check(observed: float, expected: float, tolerance: float, unit: str = "absolute") -> dict:
    return {"observed": float(observed), "expected": float(expected),
            "tolerance": tolerance, "unit": unit,
            "status": "pass" if np.isfinite(observed) and abs(observed - expected) <= tolerance else "fail"}


@dataclass
class CaseReference:
    case_dir: Path
    fchk_path: Path
    analysis_dir: Path
    fchk_source: str
    archive_header: dict
    archive_mo: dict[str, np.ndarray]
    archive_overlap: np.ndarray
    archive_density: dict[str, np.ndarray]
    mol: object
    basis: tuple
    overlap: np.ndarray
    aonbo: dict[str, np.ndarray]
    nbomo: dict[str, np.ndarray]
    mo_coeff: dict[str, np.ndarray]
    mo_occupations: dict[str, np.ndarray]
    nbo_coeff: dict[str, np.ndarray]
    mo_source: dict[str, str]

    @property
    def nbasis(self) -> int:
        return self.overlap.shape[0]

    @property
    def spins(self) -> tuple[str, ...]:
        return tuple(self.aonbo)

    def matrix_checks(self) -> dict:
        """Report matrix, AO metric, representation, and density discrepancies."""
        eye = np.eye(self.nbasis)
        checks: dict[str, dict] = {}
        checks["sao_symmetric"] = _check(np.max(np.abs(self.overlap - self.overlap.T)), 0, 2e-9)
        checks["sao_vs_archive_overlap"] = _check(
            np.max(np.abs(self.overlap - self.archive_overlap)), 0, 2e-8)
        checks["sao_positive_min_eigenvalue"] = {
            "observed": float(np.linalg.eigvalsh(self.overlap)[0]), "expected": ">0",
            "tolerance": 0.0, "status": "pass" if np.linalg.eigvalsh(self.overlap)[0] > 0 else "fail"}
        total_density = np.zeros_like(self.overlap)
        for spin in self.spins:
            mo = self.mo_coeff[spin]
            occ = self.mo_occupations[spin]
            nbo = self.nbo_coeff[spin]
            raw = self.aonbo[spin]
            t = self.nbomo[spin]
            checks[f"{spin}_nbomo_orthogonal"] = _check(np.max(np.abs(t.T @ t - eye)), 0, 1e-5)
            checks[f"{spin}_mo_ao_metric"] = _check(np.max(np.abs(mo.T @ self.overlap @ mo - eye)), 0, 1e-5)
            checks[f"{spin}_nbo_ao_metric"] = _check(np.max(np.abs(raw.T @ self.overlap @ raw - eye)), 0, 1e-5)
            checks[f"{spin}_aonbo_vs_selected_mo_nbomo"] = _check(np.max(np.abs(raw - nbo)), 0, 2e-4)
            checks[f"{spin}_aonbo_vs_selected_mo_nbomo"]["mo_source"] = self.mo_source[spin]
            if checks[f"{spin}_aonbo_vs_selected_mo_nbomo"]["status"] != "pass":
                checks[f"{spin}_aonbo_vs_selected_mo_nbomo"]["reason"] = (
                    "canonical_FCHK_MOs_do_not_reproduce_analysis_FILE37_orbitals; "
                    "possible_SCF_orbital_rotation; no orbitalwise grid acceptance")
            archive_nbo = self.archive_mo[spin] @ t.T
            checks[f"{spin}_aonbo_vs_archive_nbomo"] = _check(np.max(np.abs(raw - archive_nbo)), 0, 2e-4)
            checks[f"{spin}_fchk_vs_archive_mo"] = _check(
                np.max(np.abs(mo - self.archive_mo[spin])), 0, 2e-4)
            if self.mo_source[spin] != "paired_fchk":
                checks[f"{spin}_fchk_vs_archive_mo"]["status"] = "unknown"
                checks[f"{spin}_fchk_vs_archive_mo"]["reason"] = "paired_FCHK_has_no_beta_MO_coefficients"
                checks[f"{spin}_fchk_vs_archive_mo"]["observed"] = None
            density_mo = (mo * occ[np.newaxis, :]) @ mo.T
            checks[f"{spin}_archive_density"] = _check(
                np.max(np.abs(density_mo - self.archive_density[spin])), 0, 1e-5)
            density_nbo = nbo @ ((t * occ[np.newaxis, :]) @ t.T) @ nbo.T
            checks[f"{spin}_density_representation"] = _check(np.max(np.abs(density_mo - density_nbo)), 0, 1e-5)
            total_density += density_nbo
        scf = _fchk_packed(self.fchk_path, "Total SCF Density", self.nbasis)
        if scf is None or np.shape(scf) != total_density.shape:
            checks["fchk_scf_density"] = {"status": "unknown", "reason": "missing_or_shape_mismatched_scf_density"}
        else:
            checks["fchk_scf_density"] = _check(np.max(np.abs(scf - total_density)), 0, 1e-5)
        if len(self.spins) == 2:
            spin_density = _fchk_packed(self.fchk_path, "Spin SCF Density", self.nbasis)
            if spin_density is None:
                checks["fchk_spin_density"] = {"status": "unknown", "reason": "missing_Spin_SCF_Density"}
            else:
                alpha = (self.mo_coeff["alpha"] * self.mo_occupations["alpha"][np.newaxis, :]) @ self.mo_coeff["alpha"].T
                beta = (self.mo_coeff["beta"] * self.mo_occupations["beta"][np.newaxis, :]) @ self.mo_coeff["beta"].T
                checks["fchk_spin_density"] = _check(np.max(np.abs(spin_density - (alpha - beta))), 0, 1e-5)
        representation_checks = [item for key, item in checks.items()
                                 if not key.endswith("_fchk_vs_archive_mo")]
        return {"case_dir": str(self.case_dir), "fchk": str(self.fchk_path),
                "fchk_source": self.fchk_source, "nbasis": self.nbasis,
                "mo_source": self.mo_source,
                "spins": list(self.spins), "archive_header": self.archive_header,
                "checks": checks,
                "archive_representation_status": "pass" if all(
                    item["status"] == "pass" for item in representation_checks) else "not_pass",
                "status": "pass" if all(item["status"] == "pass" for item in checks.values()) else "not_pass"}

    def evaluate(self, points_bohr, orbital_set: str = "nbo", spin: str = "alpha",
                 source_one_based=None) -> np.ndarray:
        """Return [orbital, point] values on explicit bohr coordinates.

        orbital_set='nbo' uses FCHK canonical MOs and the full vendor NBOMO
        transform. 'aonbo' exposes the raw FILE.37 AO coefficients for a
        representation cross-check. Indices are 1-based source orbital numbers.
        """
        points = np.asarray(points_bohr, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ReferenceInputError("points_bohr:expected_finite_Nx3_array")
        if spin not in self.spins:
            raise ReferenceInputError(f"spin:{spin}:available_{self.spins}")
        arrays = {"nbo": self.nbo_coeff, "aonbo": self.aonbo,
                  "mo": self.mo_coeff, "archive_nbo": {
                      key: self.archive_mo[key] @ self.nbomo[key].T for key in self.spins}}
        if orbital_set not in arrays:
            raise ReferenceInputError(f"orbital_set:{orbital_set}:expected_nbo_aonbo_mo_or_archive_nbo")
        coeff = arrays[orbital_set][spin]
        if orbital_set == "nbo":
            disagreement = float(np.max(np.abs(self.nbo_coeff[spin] - self.aonbo[spin])))
            if disagreement > 2e-4:
                raise ReferenceInputError(
                    f"canonical_FCHK_x_NBOMO_disagrees_with_FILE37:{spin}:max_abs={disagreement:.6g}")
        if source_one_based is None:
            selected = np.arange(coeff.shape[1])
        else:
            items = [source_one_based] if isinstance(source_one_based, int) else list(source_one_based)
            if not items or any(not isinstance(i, int) or i < 1 or i > coeff.shape[1] for i in items):
                raise ReferenceInputError("source_one_based:out_of_range_or_empty")
            selected = np.asarray(items, dtype=int) - 1
        ao = evaluate_basis(self.basis, points)
        if ao.shape != (self.nbasis, len(points)) or not np.isfinite(ao).all():
            raise ReferenceInputError(f"AO_evaluation:unexpected_shape_or_nonfinite_{ao.shape}")
        return coeff[:, selected].T @ ao


def load_reference(case_dir: str | Path, fchk_source: str = "paired") -> CaseReference:
    """Load one producer case, refusing missing or ambiguous matrix blocks."""
    case_dir = Path(case_dir).resolve()
    analysis_dir = case_dir / ("analysis-03" if case_dir.name == "NBO-10" else "analysis-02")
    if fchk_source == "paired":
        fchk_path = _require_file(analysis_dir / "paired-canonical.fchk")
    elif fchk_source == "original":
        fchk_path = _require_file(case_dir / "canonical" / "canonical.fchk")
    else:
        raise ReferenceInputError(f"fchk_source:{fchk_source}:expected_paired_or_original")
    mol = load_one(str(fchk_path))
    if mol.mo is None or mol.mo.kind not in ("restricted", "unrestricted"):
        raise ReferenceInputError(f"{fchk_path}:unsupported_or_missing_MOs")
    nbasis = int(np.shape(mol.mo.coeffs)[0])
    archive = _archive_header(analysis_dir / "FILE.47", nbasis)
    unrestricted = archive["open_shell"]
    archive_mo = _archive_lcaomo(analysis_dir / "FILE.47", nbasis, unrestricted)
    archive_overlap = _archive_packed(analysis_dir / "FILE.47", "OVERLAP", nbasis)[0]
    density_blocks = _archive_packed(analysis_dir / "FILE.47", "DENSITY", nbasis,
                                     2 if unrestricted else 1)
    archive_density = {"alpha": density_blocks[0]}
    if unrestricted:
        archive_density["beta"] = density_blocks[1]
    aonbo = _read_w_square(analysis_dir / "FILE.37", nbasis, "NBOs in the AO basis:", unrestricted)
    nbomo = _read_w_square(analysis_dir / "FILE.49", nbasis, "MOs in the NBO basis:", unrestricted)
    overlap = _read_sao(analysis_dir / "FILE.50", nbasis)
    full_coeff = np.asarray(mol.mo.coeffs, dtype=np.float64)
    full_occ = np.asarray(mol.mo.occs, dtype=np.float64)
    if full_coeff.shape not in ((nbasis, nbasis * (2 if unrestricted else 1)),
                                (nbasis, nbasis)):
        raise ReferenceInputError(f"{fchk_path}:incomplete_MO_coefficients_{full_coeff.shape}")
    mo_coeff = {"alpha": full_coeff[:, :nbasis]}
    mo_occ = {"alpha": full_occ[:nbasis]}
    mo_source = {"alpha": "paired_fchk" if fchk_source == "paired" else "original_fchk"}
    if unrestricted:
        if full_coeff.shape[1] == 2 * nbasis:
            mo_coeff["beta"] = full_coeff[:, nbasis:]
            mo_occ["beta"] = full_occ[nbasis:]
            mo_source["beta"] = "paired_fchk" if fchk_source == "paired" else "original_fchk"
        else:
            mo_coeff["beta"] = archive_mo["beta"]
            nalpha = _fchk_electron_count(fchk_path, "alpha")
            nbeta = _fchk_electron_count(fchk_path, "beta")
            mo_occ["alpha"] = np.array([1.] * nalpha + [0.] * (nbasis - nalpha))
            mo_occ["beta"] = np.array([1.] * nbeta + [0.] * (nbasis - nbeta))
            mo_source["beta"] = "archive_LCAOMO; paired_FCHK_beta_missing"
    nbo_coeff = {spin: mo_coeff[spin] @ nbomo[spin].T for spin in aonbo}
    return CaseReference(case_dir, fchk_path, analysis_dir, fchk_source, archive,
                         archive_mo, archive_overlap, archive_density, mol,
                         tuple(from_iodata(mol)), overlap, aonbo, nbomo,
                         mo_coeff, mo_occ, nbo_coeff, mo_source)


def extract_raw_tables(case_dir: str | Path) -> dict:
    """Extract printed rows with line provenance for COV value comparisons.

    This is intentionally a narrow extraction of NBO's own tables. It retains
    raw lines and visible rounded values. An unparsed candidate is reported,
    never turned into a numeric zero. Wiberg print blocks are reassembled using
    their explicit column numbers, including benzene's wrapped columns.
    """
    case_dir = Path(case_dir).resolve()
    analysis = case_dir / ("analysis-03" if case_dir.name == "NBO-10" else "analysis-02")
    path = _require_file(analysis / "analysis.log")
    lines = path.read_text(errors="replace").splitlines()
    archive_path = _require_file(analysis / "FILE.47")
    basis_match = re.search(r"NBAS=(\d+)", archive_path.read_text(errors="replace")[:200])
    if not basis_match:
        raise ReferenceInputError(f"{archive_path}:missing_NBAS")
    archive = _archive_header(archive_path, int(basis_match.group(1)))
    natoms, nbasis = archive["natoms"], archive["nbasis"]
    open_shell = archive["open_shell"]
    starts = {
        "population": [i for i, line in enumerate(lines) if "NATURAL POPULATIONS:  Natural atomic orbital occupancies" in line],
        "nbo": [i for i, line in enumerate(lines) if "NATURAL BOND ORBITAL ANALYSIS" in line],
        "e2": [i for i, line in enumerate(lines) if "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX" in line],
        "summary": [i for i, line in enumerate(lines) if "NATURAL BOND ORBITALS (Summary):" in line],
    }
    population_scopes = ["total", "alpha", "beta"] if open_shell else ["total"]
    nbo_scopes = ["alpha", "beta"] if open_shell else ["alpha"]
    problems: list[str] = []
    if len(starts["population"]) != len(population_scopes):
        problems.append(f"population_sections:expected_{len(population_scopes)}_found_{len(starts['population'])}")
    if len(starts["nbo"]) != len(nbo_scopes):
        problems.append(f"nbo_sections:expected_{len(nbo_scopes)}_found_{len(starts['nbo'])}")
    if len(starts["e2"]) != len(nbo_scopes):
        problems.append(f"e2_sections:expected_{len(nbo_scopes)}_found_{len(starts['e2'])}")
    if len(starts["summary"]) != len(nbo_scopes):
        problems.append(f"summary_sections:expected_{len(nbo_scopes)}_found_{len(starts['summary'])}")
    summary_ends = [i for i, line in enumerate(lines) if "NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:" in line]
    all_heads = sorted(set(sum(starts.values(), []) + summary_ends)) + [len(lines)]
    def section(start: int) -> list[tuple[int, str]]:
        stop = next(i for i in all_heads if i > start)
        return [(i + 1, lines[i]) for i in range(start, stop)]

    out = {"case_dir": str(case_dir), "log": str(path), "natoms": natoms,
           "nbasis": nbasis, "population": {}, "orbitals": {}, "e2": {},
           "problems": problems}
    number = r"[-+]?(?:\d+\.\d+|\.\d+)(?:[EeDd][-+]?\d+)?"
    atom_row = re.compile(r"^\s*([A-Z][a-z]?)\s+(\d+)\s+(.+?)\s*$")
    nao_row = re.compile(rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+(\d+)\s+(\S+)\s+(Cor|Val|Ryd)\(.*?\)\s+({number})\s+({number})\s*$")
    wiberg_row = re.compile(r"^\s*(\d+)\.\s+([A-Z][a-z]?)\s+(.+?)\s*$")
    nbo_row = re.compile(rf"^\s*(\d+)\.\s+\(({number})\)\s+([A-Za-z0-9]+\*?)\s*(.*?)\s*$")
    component_row = re.compile(rf"^\s*\(\s*({number})%\)\s+({number})\*\s*([A-Z][a-z]?)\s+(\d+)\s*(.*?)\s*$")
    e2_row = re.compile(rf"^\s*(\d+)\.\s+(.+?)\s+(\d+)\.\s+(.+?)\s+({number})\s+({number})\s+({number})\s*$")
    summary_row = re.compile(rf"^\s*(\d+)\.\s+(.+?)\s+({number})\s+({number})(?:\s+.*)?$")

    for scope, begin in zip(population_scopes, starts["population"]):
        rows = section(begin)
        nao_header = next((line for _, line in rows[:8]
                           if "NAO Atom No" in line and "Occupancy" in line), "")
        nao_tail_kind = "spin_density" if re.search(r"Occupancy\s+Spin\s*$", nao_header) else "energy"
        summary_idx = next((j for j, (_, line) in enumerate(rows) if "Summary of Natural Population Analysis:" in line), None)
        wiberg_idx = next((j for j, (_, line) in enumerate(rows) if "Wiberg bond index matrix in the NAO basis:" in line), None)
        if summary_idx is None or wiberg_idx is None or summary_idx >= wiberg_idx:
            problems.append(f"{scope}:missing_NPA_or_Wiberg_heading")
            continue
        nao = []
        for line_no, line in rows[:summary_idx]:
            match = nao_row.match(line)
            if match:
                type_match = re.search(r"\b(?:Cor|Val|Ryd)\(\s*[^)]*\)", line)
                nao.append({"line": line_no, "index": int(match.group(1)),
                            "element": match.group(2), "atom": int(match.group(3)),
                            "angular": match.group(4), "kind": match.group(5),
                            "type": type_match.group(0) if type_match else None,
                            "occupation": float(match.group(6)),
                            "energy": float(match.group(7)) if nao_tail_kind == "energy" else None,
                            "spin_density": float(match.group(7)) if nao_tail_kind == "spin_density" else None,
                            "raw": line})
        atoms = []
        for line_no, line in rows[summary_idx:wiberg_idx]:
            match = atom_row.match(line)
            if match:
                words = match.group(3).split()
                if len(words) not in (5, 6) or not all(re.fullmatch(number, word) for word in words):
                    continue
                values = [float(word.replace("D", "E")) for word in words]
                atoms.append({"line": line_no, "element": match.group(1),
                              "atom": int(match.group(2)), "charge": values[0],
                              "core": values[1], "valence": values[2],
                              "rydberg": values[3], "population": values[4],
                              "spin_density": values[5] if len(values) == 6 else None,
                              "raw": line})
        matrix_rows: dict[int, dict[int, float]] = {}
        columns = []
        for line_no, line in rows[wiberg_idx + 1:]:
            stripped = line.strip()
            if stripped.startswith("Wiberg bond index, Totals by atom:"):
                break
            if stripped.startswith("Atom"):
                columns = [int(word) for word in stripped.split()[1:] if word.isdigit()]
                continue
            match = wiberg_row.match(line)
            if match and columns:
                raw_values = match.group(3).split()
                if len(raw_values) == len(columns) and all(_float_tokens(w) for w in raw_values):
                    index = int(match.group(1))
                    row = matrix_rows.setdefault(index, {})
                    for col, value in zip(columns, raw_values):
                        if col in row:
                            problems.append(f"{scope}:duplicate_Wiberg_cell_{index}_{col}")
                        row[col] = float(value)
        if len(nao) != nbasis:
            problems.append(f"{scope}:NAO_expected_{nbasis}_found_{len(nao)}")
        if len(atoms) != natoms:
            problems.append(f"{scope}:NPA_atoms_expected_{natoms}_found_{len(atoms)}")
        if len(matrix_rows) != natoms or any(len(matrix_rows.get(i, {})) != natoms for i in range(1, natoms + 1)):
            problems.append(f"{scope}:Wiberg_expected_{natoms}x{natoms}_found_"
                            f"{sum(len(row) for row in matrix_rows.values())}_cells")
        out["population"][scope] = {"nao_rows": nao, "npa_atoms": atoms,
                                     "wiberg": matrix_rows,
                                     "nao_tail_kind": nao_tail_kind,
                                     "counts": {"nao": len(nao), "npa_atoms": len(atoms),
                                                "wiberg_cells": sum(len(row) for row in matrix_rows.values())}}

    for scope, begin in zip(nbo_scopes, starts["nbo"]):
        rows = section(begin)
        orbitals = []
        current = None
        for line_no, line in rows:
            match = nbo_row.match(line)
            if match:
                description = match.group(4)
                ordinal_match = re.match(r"^\(\s*(\d+)\)\s*(.*)$", description)
                ordinal = int(ordinal_match.group(1)) if ordinal_match else None
                tail = ordinal_match.group(2) if ordinal_match else description
                atom_matches = list(re.finditer(r"([A-Z][a-z]?)\s+(\d+)", tail))
                atom_indices = [int(item.group(2)) for item in atom_matches]
                current = {"line": line_no, "index": int(match.group(1)),
                           "occupation": float(match.group(2)),
                           "kind": match.group(3), "ordinal": ordinal,
                           "description": description, "atoms": atom_indices,
                           "components": [], "raw": line}
                if len(atom_matches) == 1:
                    hybrid = (tail[:atom_matches[0].start()] + tail[atom_matches[0].end():]).strip()
                    current["components"].append({"line": line_no, "percentage": 100.0,
                                                  "coefficient": 1.0, "element": atom_matches[0].group(1),
                                                  "atom": atom_indices[0], "hybrid": hybrid,
                                                  "raw": line})
                orbitals.append(current)
                continue
            match = component_row.match(line)
            if match and current is not None:
                current["components"].append({
                    "line": line_no, "percentage": float(match.group(1)),
                    "coefficient": float(match.group(2)),
                    "element": match.group(3), "atom": int(match.group(4)),
                    "hybrid": match.group(5), "raw": line})
        if len(orbitals) != nbasis or sorted(item["index"] for item in orbitals) != list(range(1, nbasis + 1)):
            problems.append(f"{scope}:NBO_expected_indices_1_to_{nbasis}_found_{len(orbitals)}")
        out["orbitals"][scope] = orbitals

    for scope, begin in zip(nbo_scopes, starts["summary"]):
        orbitals_by_id = {item["index"]: item for item in out["orbitals"].get(scope, [])}
        seen = set()
        for line_no, line in section(begin):
            match = summary_row.match(line)
            if not match:
                continue
            index = int(match.group(1))
            if index not in orbitals_by_id or index in seen:
                continue
            seen.add(index)
            orbitals_by_id[index]["summary_occupation"] = float(match.group(3))
            orbitals_by_id[index]["diagonal_fock_hartree"] = float(match.group(4))
            orbitals_by_id[index]["summary_line"] = line_no
            orbitals_by_id[index]["summary_raw"] = line
        if len(seen) != nbasis:
            problems.append(f"{scope}:NBO_summary_expected_{nbasis}_found_{len(seen)}")

    for scope, begin in zip(nbo_scopes, starts["e2"]):
        rows = section(begin)
        interactions = []
        candidates = []
        threshold = next((float(match.group(1)) for _, line in rows
                          if (match := re.search(r"Threshold for printing:\s+([0-9.]+) kcal/mol", line))), None)
        for line_no, line in rows:
            match = e2_row.match(line)
            if match:
                interactions.append({"line": line_no, "donor": int(match.group(1)),
                                     "donor_label": match.group(2).strip(),
                                     "acceptor": int(match.group(3)),
                                     "acceptor_label": match.group(4).strip(),
                                     "e2_kcal_mol": float(match.group(5)),
                                     "delta_hartree": float(match.group(6)),
                                     "fock_hartree": float(match.group(7)), "raw": line})
            elif re.match(r"^\s*\d+\.\s+(?:CR|LP|BD|BD\*|RY|LV|3C)", line):
                candidates.append({"line": line_no, "raw": line})
        if candidates:
            problems.append(f"{scope}:unparsed_E2_candidates_{len(candidates)}")
        if not interactions:
            problems.append(f"{scope}:no_E2_rows")
        out["e2"][scope] = {"rows": interactions, "count": len(interactions),
                            "printed_threshold_kcal_mol": threshold,
                            "unparsed_candidates": candidates}
    out["status"] = "pass" if not problems else "not_pass"
    return out


def raw_table_checks(log_path: str | Path, dataset_json_object: dict) -> dict:
    """Compare every printed table row and rounded value with production JSON.

    The independent extraction uses only the original NBO log. The function
    reports every missing, unexpected, and unequal row; no missing entry is
    treated as zero. A CH3 total NPA row is distinct from alpha/beta rows.
    """
    log_path = Path(log_path).resolve()
    if log_path.name != "analysis.log":
        raise ReferenceInputError(f"expected_analysis_log:{log_path}")
    raw = extract_raw_tables(log_path.parent.parent)
    if log_path != Path(raw["log"]):
        raise ReferenceInputError(f"log_path_mismatch:{log_path}:{raw['log']}")
    details: list[dict] = [{"collection": "raw", "problem": item} for item in raw["problems"]]
    checks: dict[str, dict] = {}
    open_shell = len(raw["orbitals"]) == 2

    def spin_name(scope: str) -> str:
        return scope if open_shell else "total"

    def numeric_equal(left, right, tolerance: float = 1e-6) -> bool:
        return (isinstance(left, (int, float)) and isinstance(right, (int, float))
                and not isinstance(right, bool)
                and np.isfinite(right) and abs(float(left) - float(right)) <= tolerance)

    def row_source_equal(expected: dict, actual: dict) -> bool:
        source = actual.get("source") or {}
        return source.get("line_begin") == expected["line"] and source.get("raw") == expected["raw"]

    def compare(name: str, expected: list[dict], actual: object,
                key_fields: tuple[str, ...], text_fields: tuple[str, ...],
                numeric_fields: tuple[str, ...], *, source=True, tolerance=1e-6) -> None:
        before = len(details)
        if not isinstance(actual, list):
            details.append({"collection": name, "problem": "production_collection_missing_or_not_array"})
            actual = []
        def grouped(rows):
            mapping = {}
            for row in rows:
                try:
                    key = tuple(row[field] for field in key_fields)
                except (KeyError, TypeError):
                    details.append({"collection": name, "problem": "row_missing_identity", "row": row})
                    continue
                if key in mapping:
                    details.append({"collection": name, "problem": "duplicate_identity", "key": key})
                mapping[key] = row
            return mapping
        exp_map, got_map = grouped(expected), grouped(actual)
        for key in sorted(exp_map.keys() - got_map.keys()):
            details.append({"collection": name, "problem": "missing_production_row", "key": key,
                            "expected_line": exp_map[key].get("line")})
        for key in sorted(got_map.keys() - exp_map.keys()):
            details.append({"collection": name, "problem": "unexpected_production_row", "key": key,
                            "actual_line": (got_map[key].get("source") or {}).get("line_begin")})
        for key in sorted(exp_map.keys() & got_map.keys()):
            left, right = exp_map[key], got_map[key]
            for field in text_fields:
                if left.get(field) != right.get(field):
                    details.append({"collection": name, "problem": "text_mismatch",
                                    "key": key, "field": field, "expected": left.get(field), "actual": right.get(field)})
            for field in numeric_fields:
                if not numeric_equal(left.get(field), right.get(field), tolerance):
                    details.append({"collection": name, "problem": "numeric_mismatch",
                                    "key": key, "field": field, "expected": left.get(field),
                                    "actual": right.get(field), "tolerance": tolerance})
            if source and not row_source_equal(left, right):
                details.append({"collection": name, "problem": "source_line_mismatch",
                                "key": key, "expected_line": left.get("line"),
                                "actual_source": right.get("source")})
        checks[name] = {"status": "pass" if len(details) == before else "fail",
                        "expected_rows": len(expected), "production_rows": len(actual),
                        "difference_count": len(details) - before}

    population = []
    naos = []
    wiberg = []
    for scope, tables in raw["population"].items():
        for row in tables["npa_atoms"]:
            population.append({"spin": scope, "atom": row["atom"], "symbol": row["element"],
                               "charge": row["charge"], "core": row["core"],
                               "valence": row["valence"], "rydberg": row["rydberg"],
                               "total": row["population"],
                               "spin_density": row["spin_density"],
                               "line": row["line"], "raw": row["raw"]})
        for row in tables["nao_rows"]:
            naos.append({"spin": scope, "id": row["index"], "atom": row["atom"],
                         "angular": row["angular"], "type": row["type"],
                         "occupation": row["occupation"],
                         "diagonal_fock_hartree": row["energy"],
                         "spin_density": row["spin_density"],
                         "line": row["line"], "raw": row["raw"]})
        for atom_a, columns in tables["wiberg"].items():
            for atom_b, value in columns.items():
                wiberg.append({"spin": scope, "atom_a": int(atom_a), "atom_b": int(atom_b),
                               "value": value})
    compare("populations", population, dataset_json_object.get("populations"),
            ("spin", "atom"), ("symbol",), ("charge", "core", "valence", "rydberg", "total"))
    before_pop_spin = len(details)
    actual_pops = {(row.get("spin"), row.get("atom")): row
                   for row in dataset_json_object.get("populations", [])}
    for row in population:
        key = (row["spin"], row["atom"])
        actual = actual_pops.get(key)
        if actual is None:
            continue
        if row["spin_density"] is None:
            if actual.get("spin_density") is not None:
                details.append({"collection": "populations", "problem": "unprinted_spin_density_populated",
                                "key": key, "actual": actual.get("spin_density")})
        elif not numeric_equal(row["spin_density"], actual.get("spin_density")):
            details.append({"collection": "populations", "problem": "spin_density_mismatch",
                            "key": key, "expected": row["spin_density"],
                            "actual": actual.get("spin_density"), "tolerance": 1e-6})
    checks["populations"]["difference_count"] += len(details) - before_pop_spin
    if len(details) != before_pop_spin:
        checks["populations"]["status"] = "fail"
    compare("naos", naos, dataset_json_object.get("naos"),
            ("spin", "id"), ("atom", "angular", "type"),
            ("occupation",))
    # The total open-shell NAO table prints Spin in the last column, whereas
    # alpha and beta tables print Energy. Compare the meaning as well as value.
    before_nao_tail = len(details)
    actual_naos = {(row.get("spin"), row.get("id")): row
                   for row in dataset_json_object.get("naos", [])}
    for row in naos:
        key = (row["spin"], row["id"])
        actual = actual_naos.get(key)
        if actual is None:
            continue
        for field in ("diagonal_fock_hartree", "spin_density"):
            expected_value = row[field]
            actual_value = actual.get(field)
            if expected_value is None:
                if actual_value is not None:
                    details.append({"collection": "naos", "problem": "unprinted_field_populated",
                                    "key": key, "field": field, "actual": actual_value})
            elif not numeric_equal(expected_value, actual_value):
                details.append({"collection": "naos", "problem": "numeric_mismatch",
                                "key": key, "field": field, "expected": expected_value,
                                "actual": actual_value, "tolerance": 1e-6})
    checks["naos"]["difference_count"] += len(details) - before_nao_tail
    if len(details) != before_nao_tail:
        checks["naos"]["status"] = "fail"
    compare("wiberg", wiberg, dataset_json_object.get("wiberg"),
            ("spin", "atom_a", "atom_b"), (), ("value",), source=False, tolerance=1e-5)

    orbitals = []
    components = []
    for scope, rows in raw["orbitals"].items():
        spin = spin_name(scope)
        for row in rows:
            orbitals.append({"spin": spin, "id": row["index"], "ordinal": row["ordinal"],
                             "kind": row["kind"], "occupation": row["occupation"],
                             "diagonal_fock_hartree": row.get("diagonal_fock_hartree"),
                             "atoms": row["atoms"], "line": row["line"], "raw": row["raw"]})
            for ordinal, component in enumerate(row["components"]):
                components.append({"spin": spin, "orbital": row["index"],
                                   "component": ordinal, "atom": component["atom"],
                                   "percent": component["percentage"],
                                   "coefficient": component["coefficient"],
                                   "hybrid": component["hybrid"],
                                   "line": component["line"], "raw": component["raw"]})
    compare("orbitals", orbitals, dataset_json_object.get("orbitals"),
            ("spin", "id"), ("ordinal", "kind", "atoms"),
            ("occupation", "diagonal_fock_hartree"))
    production_components = []
    for orbital in dataset_json_object.get("orbitals", []):
        for ordinal, component in enumerate(orbital.get("components", [])):
            production_components.append({"spin": orbital.get("spin"), "orbital": orbital.get("id"),
                                          "component": ordinal, **component})
    compare("components", components, production_components,
            ("spin", "orbital", "component"), ("atom", "hybrid"),
            ("percent", "coefficient"))

    expected_e2 = []
    for scope, block in raw["e2"].items():
        spin = spin_name(scope)
        for ordinal, row in enumerate(block["rows"]):
            expected_e2.append({"spin": spin, "ordinal": ordinal,
                                "donor": row["donor"], "acceptor": row["acceptor"],
                                "value": row["e2_kcal_mol"],
                                "energy_gap_hartree": row["delta_hartree"],
                                "fock_hartree": row["fock_hartree"],
                                "printing_threshold": block["printed_threshold_kcal_mol"],
                                "line": row["line"], "raw": row["raw"]})
    production_e2 = []
    e2_ordinals = {}
    for row in dataset_json_object.get("e2", []):
        spin = row.get("spin")
        ordinal = e2_ordinals.get(spin, 0)
        e2_ordinals[spin] = ordinal + 1
        production_e2.append({"ordinal": ordinal, **row})
    compare("e2", expected_e2, production_e2,
            ("spin", "ordinal"), ("donor", "acceptor"),
            ("value", "energy_gap_hartree", "fock_hartree", "printing_threshold"))
    expected_sections = [{"spin": spin_name(scope),
                          "printing_threshold": block["printed_threshold_kcal_mol"]}
                         for scope, block in raw["e2"].items()]
    compare("e2_sections", expected_sections, dataset_json_object.get("e2_sections"),
            ("spin",), (), ("printing_threshold",), source=False)
    return {"status": "pass" if not details else "fail", "checks": checks,
            "details": details, "raw_extraction_status": raw["status"],
            "log": str(log_path)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args()
    try:
        report = load_reference(args.case_dir).matrix_checks()
    except Exception as error:
        report = {"case_dir": str(args.case_dir), "status": "not_pass",
                  "reason": f"{type(error).__name__}: {error}"}
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
