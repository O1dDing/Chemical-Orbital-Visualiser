"""Inventory the frozen 273 Gaussian cases for AO-MO/NBO coverage planning.

This is read-only with respect to the source corpus. It records file identity,
calculation metadata, category and evidence gaps; it never starts chemistry.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import re


WORK = Path(__file__).resolve().parents[2]
MANIFEST = WORK / "source/validation/view-snapshots-20260906/original-frozen-manifest.json"
OUTPUT = WORK / "aomo-coverage"
CASE_GROUPS = {
    "gaussian_validation_set": "baseline_chemistry",
    "molecular_ion_crosscheck_20260831": "ion_crosscheck",
    "molecular_ion_gold_opt_20260831": "ion_optimized_reference",
    "molecular_ion_master_20260831": "molecular_ion_master",
    "molecular_scan_master_20260831": "geometry_scan_frame",
    "(root)": "single_other",
    "tm_coordination_geometry_matrix_20260830": "coordination_geometry_sample",
    "tm_full_geometry_matrix_20260830": "full_coordination_geometry",
    "tm_ligand_field_matrix_20260829": "ligand_field",
}
SHELL_NAMES = {0:"s", 1:"p", -1:"SP", 2:"Cartesian_d", -2:"spherical_d",
               3:"Cartesian_f", -3:"spherical_f", 4:"Cartesian_g", -4:"spherical_g"}
REFERENCE_EXACT_SHA = set()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def migrated(path: str) -> Path:
    return Path(re.sub(r"^[Ff]:\\", "E:\\\\", path).replace("\\", "/"))


def parse_fchk(path: Path) -> dict:
    fields: dict[str, int] = {}
    arrays: dict[str, list[int]] = {}
    wanted_arrays = {"Atomic numbers", "Shell types"}
    wanted_scalars = {"Number of atoms", "Charge", "Multiplicity", "Number of electrons",
                      "Number of alpha electrons", "Number of beta electrons",
                      "Number of basis functions", "Number of independent functions"}
    active = None
    remaining = 0
    header = []
    coefficient_blocks = []
    orbital_energy_blocks = []
    occupations = []
    with path.open("r", errors="replace") as stream:
        for line_no, line in enumerate(stream):
            if line_no < 2:
                header.append(line.rstrip("\r\n"))
            if active:
                try:
                    values = [int(word) for word in line.split()]
                except ValueError:
                    raise ValueError(f"malformed_{active}_array:{path}")
                arrays[active].extend(values)
                remaining -= len(values)
                if remaining < 0:
                    raise ValueError(f"too_many_{active}_values:{path}")
                if remaining == 0:
                    active = None
                continue
            for key in wanted_scalars:
                if line.startswith(key) and re.search(r"\bI\s+-?\d+\s*$", line):
                    fields[key] = int(line.split()[-1])
                    break
            if "MO coefficients" in line and re.search(r"\bR\s+N=", line):
                coefficient_blocks.append(line.split("MO coefficients")[0].strip())
            if "Orbital Energies" in line and re.search(r"\bR\s+N=", line):
                orbital_energy_blocks.append(line.split("Orbital Energies")[0].strip())
            if "Orbital Occupation" in line or "MO Occupation" in line:
                occupations.append(line.strip())
            for key in wanted_arrays:
                if line.startswith(key):
                    match = re.search(r"\bI\s+N=\s*(\d+)", line)
                    if match:
                        active = key
                        remaining = int(match.group(1))
                        arrays[key] = []
                        if remaining == 0:
                            active = None
                    break
    if active:
        raise ValueError(f"truncated_{active}_array:{path}")
    required = wanted_scalars - {"Number of independent functions"}
    if not required.issubset(fields):
        raise ValueError(f"missing_fields_{sorted(required-fields.keys())}:{path}")
    method_basis = header[1].split() if len(header) > 1 else []
    shell_types = arrays.get("Shell types", [])
    atomic_numbers = arrays.get("Atomic numbers", [])
    if len(atomic_numbers) != fields["Number of atoms"]:
        raise ValueError(f"atomic_number_count_mismatch:{path}")
    nbas = fields["Number of basis functions"]
    nind = fields.get("Number of independent functions", nbas)
    ecp_core = sum(atomic_numbers) - fields["Charge"] - fields["Number of electrons"]
    return {
        "title": header[0].strip() if header else "",
        "job_kind": method_basis[0] if len(method_basis) > 0 else None,
        "method": method_basis[1] if len(method_basis) > 1 else None,
        "basis": method_basis[2] if len(method_basis) > 2 else None,
        "atom_count": fields["Number of atoms"], "atomic_numbers": atomic_numbers,
        "charge": fields["Charge"], "multiplicity": fields["Multiplicity"],
        "explicit_electrons": fields["Number of electrons"],
        "alpha_electrons": fields["Number of alpha electrons"],
        "beta_electrons": fields["Number of beta electrons"],
        "n_basis": nbas, "n_independent": nind,
        "ecp_core_derived": ecp_core,
        "shell_types": shell_types,
        "shell_labels": sorted({SHELL_NAMES.get(i, f"unsupported_{i}") for i in shell_types}),
        "coefficient_blocks": coefficient_blocks,
        "orbital_energy_blocks": orbital_energy_blocks,
        "explicit_occupation_fields": occupations,
    }


def identity_key(stem: str, group: str) -> tuple[str, str]:
    s = re.sub(r"^\d+_", "", stem.lower())
    s = re.sub(r"_(?:u?pbe1pbe|udft|pbe1pbe|dft|rhf|uhf|rohf)_(?:augccpvtz|def2tzvp|def2svp|def2svpd|6-31g.*)$", "", s)
    if group == "geometry_scan_frame":
        return re.sub(r"_\d{2}$", "", s), "scan_frame_not_independent_molecule"
    if "benzene_rotated" in s or "benzene_d6h" in s:
        return "benzene", "rotation_or_basis_variant_not_independent_molecule"
    if "butadiene" in s and ("twisted" in s or "trans" in s or "s_cis" in s or "s_trans" in s):
        return "butadiene", "conformer_not_independent_molecule"
    if "cp_anion" in s:
        return "cyclopentadienyl_anion", "basis_or_representation_variant_not_independent_molecule"
    return s, "identity_key_from_filename_candidate"


def chemical_tags(stem: str, group: str, meta: dict) -> list[str]:
    s = stem.lower()
    tags = set()
    if meta["charge"] > 0: tags.add("cation")
    if meta["charge"] < 0: tags.add("anion")
    if meta["multiplicity"] > 1: tags.add("open_shell")
    if any(x in s for x in ("h2o", "water", "methanol", "ammonia", "nh3", "hf", "lih", "h2_")):
        tags.update(("polar_or_sigma", "lone_pair_or_sigma"))
    if any(x in s for x in ("ethylene", "ethene", "formaldehyde", "so2", "ozone", "urea", "carbonate", "nitrate", "acrolein")):
        tags.add("double_bond_or_polar_pi")
    if any(x in s for x in ("acetylene", "propyne", "n2_", "co_reference", "_co_", "acetonitrile", "cn_minus")):
        tags.add("triple_bond_or_multiple_pi")
    if any(x in s for x in ("benzene", "butadiene", "naphthalene", "anthracene", "phenanthrene", "allyl", "cp_", "tropylium", "cyclobutadiene", "azulene", "hexatriene", "quinoline", "pyridine")):
        tags.add("conjugated_or_delocalized")
    if any(x in s for x in ("h3plus", "h3_plus", "diborane", "al2cl6", "ammonia_borane")):
        tags.add("three_center_two_electron_candidate")
    if any(x in s for x in ("xef2", "i3_minus", "cl3_minus", "hf2_minus", "pf5", "sf6", "xef4")):
        tags.add("hypervalent_or_three_center_four_electron_candidate")
    if any(x in s for x in ("dimer", "contact", "separated", "co2_ammonia", "zwitterion", "approach")):
        tags.add("noncovalent_or_ion_pair")
    if any(x in s for x in ("crco6", "nico4", "feco5", "feco", "ni_co4", "fecl", "mocn", "crcn")):
        tags.add("metal_ligand_donation_backbonding_candidate")
    if "re2cl8" in s:
        tags.add("metal_metal_delta_candidate")
    if group in ("coordination_geometry_sample", "full_coordination_geometry", "ligand_field"):
        tags.add("coordination_or_ligand_field")
    if group == "geometry_scan_frame": tags.add("geometry_scan_frame")
    if not tags: tags.add("general_gaussian_molecule")
    return sorted(tags)


def math_tags(meta: dict, stem: str) -> list[str]:
    labels = set(meta["shell_labels"])
    tags = set()
    for label in labels:
        if label == "SP": tags.add("merged_SP_shell")
        if label.startswith("Cartesian_"): tags.add(label)
        if label.startswith("spherical_"): tags.add(label)
        if label.startswith("unsupported_"): tags.add("unsupported_shell_type")
    if meta["ecp_core_derived"] > 0: tags.add("ECP")
    if meta["ecp_core_derived"] == 0: tags.add("all_electron_by_count")
    if meta["n_independent"] < meta["n_basis"]: tags.add("NBAS_gt_NMO_or_linear_dependence")
    if meta["multiplicity"] > 1: tags.add("open_shell_spin_blocks")
    if "Beta" not in meta["coefficient_blocks"] and meta["multiplicity"] > 1:
        tags.add("beta_coefficients_absent")
    b = (meta.get("basis") or "").lower()
    if "aug" in b or "+" in b or "svpd" in b:
        tags.add("diffuse_basis_candidate")
    if meta["explicit_occupation_fields"]:
        tags.add("explicit_occupation_fields_present")
    return sorted(tags)


def required_capabilities(tags: list[str]) -> list[str]:
    t = set(tags)
    needs = {"NAOMO", "AONAO", "NAONBO", "signed_AO_NAO_MO_links", "NPA_Wiberg"}
    if "open_shell" in t: needs.add("alpha_beta_SOMO_spin_population")
    if "double_bond_or_polar_pi" in t or "triple_bond_or_multiple_pi" in t:
        needs.add("sigma_pi_bond_order_overlay")
    if "conjugated_or_delocalized" in t: needs.add("fragment_pi_degeneracy_subspace")
    if "coordination_or_ligand_field" in t or "metal_ligand_donation_backbonding_candidate" in t:
        needs.update(("metal_ligand_grouping", "NHO_hybrid_direction", "E2_direction_if_printed"))
    if "metal_metal_delta_candidate" in t: needs.add("delta_channel_and_noninteger_bond_semantics")
    if "three_center_two_electron_candidate" in t or "hypervalent_or_three_center_four_electron_candidate" in t:
        needs.add("multicentre_or_hypervalent_representation")
    if "noncovalent_or_ion_pair" in t: needs.add("weak_interaction_without_forced_bond")
    if "geometry_scan_frame" in t: needs.add("input_switch_stale_overlay_clear")
    return sorted(needs)


def visible_acceptance(tags: list[str]) -> list[str]:
    t = set(tags)
    checks = ["one_drop_matches_source_or_explains_missing_NBO", "whole_AO_MO_graph_clicks_to_right_3D"]
    if "open_shell" in t: checks.append("alpha_beta_switch_preserves_SOMO_and_spin_identity")
    if "double_bond_or_polar_pi" in t: checks.append("visible_double_bond_and_sigma_pi_drilldown")
    if "triple_bond_or_multiple_pi" in t: checks.append("visible_triple_bond_and_two_pi_channels")
    if "conjugated_or_delocalized" in t: checks.append("fragment_and_degenerate_subspace_expand_without_false_unique_orbital")
    if "coordination_or_ligand_field" in t: checks.append("metal_ligand_geometry_grouping_and_pickable_overlay")
    if "three_center_two_electron_candidate" in t or "hypervalent_or_three_center_four_electron_candidate" in t:
        checks.append("multicentre_view_without_forced_integer_two_centre_bonds")
    if "noncovalent_or_ion_pair" in t: checks.append("contact_selection_without_forced_covalent_bond")
    if "geometry_scan_frame" in t: checks.append("step_change_clears_old_graph_and_orbital_field")
    return checks


def geometry_label(stem: str, group: str) -> str | None:
    if group != "full_coordination_geometry": return None
    match = re.search(r"_CN\d+_([^_]+)_", stem)
    return match.group(1) if match else None


REPRESENTATIVE_CLASSES = {
    "simple_sigma_H2":"OLD-040", "water_lone_pair":"OLD-063",
    "tetrahedral_degeneracy_CH4":"OLD-065", "benzene_pi":"OLD-005",
    "Cp_anion_pi":"OLD-019", "CO2_two_pi":"OLD-061", "O2_triplet":"OLD-049",
    "H3plus_3c2e":"OLD-001", "XeF2_3c4e":"OLD-002",
    "TiF6_oct_ligand_field":"OLD-003", "ZnCl4_tetra_ECP":"OLD-004",
    "square_planar_d8":"OLD-248", "formaldehyde_double":"OLD-068",
    "acetylene_triple":"OLD-070", "butadiene_conjugation":"OLD-008",
    "naphthalene_fused":"OLD-015", "butadiene_twisted":"OLD-016",
    "allyl_radical":"OLD-014", "diborane_3c2e":"OLD-010",
    "CrCO6_backbond":"OLD-012", "NiCO4_backbond":"OLD-272",
    "Re2Cl8_delta":"OLD-013", "water_dimer_noncovalent":"OLD-131",
    "ion_pair_contact":"OLD-173", "diffuse_allyl_anion":"OLD-032",
    "scan_H2_dissociation":"OLD-177", "scan_butadiene_torsion":"OLD-192",
    "scan_HF2_asymmetric":"OLD-221", "scan_ion_approach":"OLD-228",
    "ion_pair_contact_alternate_same_molecule":"OLD-228",
    "math_Cartesian_d_f_and_NBAS_gt_NMO":"OLD-020",
    "math_merged_SP_shell":"OLD-233",
    "math_spherical_d_f":"OLD-013",
    "math_ECP_explicit_electron_count":"OLD-002",
    "math_diffuse_anion":"OLD-032",
    "math_alpha_beta_blocks":"OLD-049",
}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    frozen = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if frozen.get("case_count") != 273 or len(frozen["cases"]) != 273:
        raise RuntimeError("unexpected_frozen_case_count")
    rows = []
    errors = []
    for source_case in frozen["cases"]:
        case_id = source_case["case_id"]
        relative = source_case["relative_path"].replace("\\", "/")
        group_raw = relative.split("/")[0] if "/" in relative else "(root)"
        group = CASE_GROUPS.get(group_raw, "unclassified_source_group")
        stem = Path(relative).stem
        fchk = migrated(source_case["input"])
        available = fchk.is_file()
        size_match = available and fchk.stat().st_size == source_case["source_size"]
        actual_sha = sha256(fchk) if available else None
        hash_match = actual_sha == source_case["sha256"]
        logs = [migrated(path) for path in source_case.get("log_candidates", [])]
        log = next((path for path in logs if path.is_file()), None)
        meta = parse_fchk(fchk) if available else {}
        identity, count_policy = identity_key(stem, group)
        chem = chemical_tags(stem, group, meta) if meta else []
        math = math_tags(meta, stem) if meta else []
        geom = geometry_label(stem, group)
        representative_for = [name for name, selected in REPRESENTATIVE_CLASSES.items() if selected == case_id]
        if group == "full_coordination_geometry": representative_for.append(f"geometry_{geom}")
        if group == "ligand_field": representative_for.append(f"ligand_field_{stem}")
        row = {
            "case_id": case_id, "relative_path": relative, "source_group": group,
            "source_fchk": str(fchk), "fchk_available": available,
            "fchk_size_bytes": fchk.stat().st_size if available else None,
            "manifest_size_match": size_match, "manifest_sha256": source_case["sha256"],
            "actual_sha256": actual_sha, "manifest_hash_match": hash_match,
            "source_log": str(log) if log else None, "log_available": log is not None,
            "source_checkpoint": None, "checkpoint_available": False,
            "source_archive47": None, "same_source_NBO_available": False,
            "title": meta.get("title"), "job_kind": meta.get("job_kind"),
            "method": meta.get("method"), "basis": meta.get("basis"),
            "atoms": meta.get("atom_count"), "atomic_numbers": meta.get("atomic_numbers"),
            "charge": meta.get("charge"), "multiplicity": meta.get("multiplicity"),
            "explicit_electrons": meta.get("explicit_electrons"),
            "ecp_core_derived": meta.get("ecp_core_derived"),
            "n_basis": meta.get("n_basis"), "n_independent": meta.get("n_independent"),
            "shell_types": meta.get("shell_types"), "shell_labels": meta.get("shell_labels"),
            "coefficient_blocks": meta.get("coefficient_blocks"),
            "explicit_occupation_fields": meta.get("explicit_occupation_fields"),
            "geometry_label": geom, "identity_key": identity, "count_policy": count_policy,
            "chemical_tags": chem, "mathematical_tags": math,
            "representative_for": representative_for,
            "required_NBO_capabilities": required_capabilities(chem),
            "user_visible_acceptance": visible_acceptance(chem),
            "status": "Gaussian_source_verified_NBO_missing" if hash_match else "source_missing_or_hash_mismatch",
            "scientific_acceptance": "not_assessed_for_AOMO_NBO",
        }
        if not hash_match: errors.append(case_id)
        rows.append(row)
    by_group = Counter(row["source_group"] for row in rows)
    geometry = {row["geometry_label"]: row["case_id"] for row in rows if row["geometry_label"]}
    if len(geometry) != 26:
        raise RuntimeError(f"expected_26_geometry_labels_found_{len(geometry)}")
    if errors:
        raise RuntimeError(f"source_hash_failures:{errors}")
    tag_counts = Counter(tag for row in rows for tag in row["chemical_tags"])
    math_counts = Counter(tag for row in rows for tag in row["mathematical_tags"])
    summary = {
        "schema": "aomo_coverage_summary_v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_manifest": str(MANIFEST), "frozen_manifest_sha256": sha256(MANIFEST),
        "case_count": len(rows), "group_counts": dict(sorted(by_group.items())),
        "fchk_available": sum(row["fchk_available"] for row in rows),
        "fchk_sha256_verified": sum(row["manifest_hash_match"] for row in rows),
        "logs_available": sum(row["log_available"] for row in rows),
        "missing_logs": [row["case_id"] for row in rows if not row["log_available"]],
        "same_source_NBO_available_in_frozen_inputs": 0,
        "scan_frame_count": by_group["geometry_scan_frame"],
        "scan_families": sorted({row["identity_key"] for row in rows if row["source_group"] == "geometry_scan_frame"}),
        "full_coordination_geometry_labels": dict(sorted(geometry.items())),
        "chemical_tag_counts": dict(sorted(tag_counts.items())),
        "mathematical_tag_counts": dict(sorted(math_counts.items())),
        "coverage_status": "inventory_complete; AOMO_NBO_class_acceptance_not_yet_run",
    }
    representatives = []
    for row in rows:
        if row["representative_for"]:
            representatives.append({"case_id": row["case_id"],
                                    "covers": row["representative_for"],
                                    "identity_key": row["identity_key"],
                                    "source_fchk": row["source_fchk"],
                                    "source_log": row["source_log"],
                                    "required_NBO_capabilities": row["required_NBO_capabilities"],
                                    "status": "candidate_Gaussian_source_only; NBO_not_produced"})
    (OUTPUT / "coverage-ledger.json").write_text(json.dumps({"schema":"aomo_coverage_ledger_v1",
        "cases":rows},indent=2)+"\n",encoding="utf-8")
    with (OUTPUT / "coverage-ledger.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=["case_id","relative_path","source_group","identity_key",
            "count_policy","geometry_label","method","basis","charge","multiplicity","n_basis",
            "n_independent","ecp_core_derived","fchk_available","manifest_hash_match","log_available",
            "chemical_tags","mathematical_tags","representative_for","required_NBO_capabilities",
            "user_visible_acceptance","status"])
        writer.writeheader()
        for row in rows:
            writer.writerow({k:";".join(row[k]) if isinstance(row.get(k),list) else row.get(k)
                             for k in writer.fieldnames})
    (OUTPUT / "coverage-summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    (OUTPUT / "representative-plan.json").write_text(json.dumps({"schema":"aomo_representative_plan_v1",
        "cases":representatives},indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"cases":len(rows),"groups":summary["group_counts"],
        "logs_available":summary["logs_available"],"geometry_labels":len(geometry),
        "representative_case_count":len(representatives),
        "scan_families":len(summary["scan_families"])}))


if __name__ == "__main__":
    main()
