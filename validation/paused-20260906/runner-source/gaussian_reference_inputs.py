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

"""Prepare traceable primary high-quality reference candidates for all old cases.

The generated inputs are candidates, not physical ground-state certificates.
Stability, competing states, environmental models and post-batch adjudication
remain separate requirements. Basis/ECP data are obtained verbatim from BSE.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time

from gaussian_rebuild import DEFAULT_OUTPUT, BOHR_ANGSTROM, digest, sha256
from validation_process import atomic_json

DEFAULT_BSE = Path(r"F:\Dev\cov-cycle-20260906\dependencies")


def get_basis_library(path):
    sys.path.insert(0, str(path))
    import basis_set_exchange as bse
    if bse.version() != "0.12":
        raise ValueError("The approved basis-data snapshot uses BSE 0.12")
    return bse


def reference_purpose(case):
    result = dict(case["purpose"])
    route = case["original_input"]["route"].lower()
    # A planar/nonplanar START is not a frozen final geometry. Original
    # optimization instructions take precedence over a filename heuristic.
    if "opt=" in route or "opt " in route:
        if result["kind"] != "scan_frame":
            result = {"kind": "equilibrium_candidate", "reference_geometry": "optimize",
                      "basis": "Original input requested optimization; final FCHK geometry is the new starting point"}
    if case["fchk_identity"]["Number of atoms"] == 1:
        result.update(kind="isolated_atom_or_ion", reference_geometry="single_atom",
                      vibration_applicability="not_applicable_no_internal_nuclear_coordinates")
    result["electronic_state"] = "specified_input_state; competing states remain to be examined when required"
    return result


def make_initial_input(reference, case, basis_text, cores=4, memory_gib=24, retry=False):
    fields = case["fchk_identity"]
    kind = reference["purpose"]["reference_geometry"]
    task = ("Opt=(VeryTight,CalcFC,MaxCycles=512) Freq" if kind == "optimize" else
            "SP" if kind == "single_atom" else "Force")
    method = "UPBE1PBE" if fields.get("has_beta_coefficients") else "RPBE1PBE"
    route = (f"#p {method}/{reference['gaussian_basis_keyword']} {task} Units=Bohr "
             f"{reference['shell_flags']} EmpiricalDispersion=GD3BJ "
             f"SCF=(VeryTight,XQC,MaxCycle={2048 if retry else 1024}) "
             "Integral=SuperFineGrid NoSymm")
    if fields.get("has_beta_coefficients") and fields["Multiplicity"] == 1:
        route += " Guess=Mix"
    lines = ["%chk=job.chk", f"%mem={memory_gib}GB", f"%nprocshared={cores}", route, "",
             reference["case_id"] + " PBE0-D3BJ primary reference candidate", "",
             f'{fields["Charge"]} {fields["Multiplicity"]}']
    xyz = fields["Current cartesian coordinates"]
    for i, z in enumerate(fields["Atomic numbers"]):
        lines.append(str(z) + " " + " ".join(f"{v:.14e}" for v in xyz[3*i:3*i+3]))
    lines.extend(["", basis_text.rstrip(), "", ""])
    return "\n".join(lines)


def make_followup_input(reference, method, stage, cores=4, memory_gib=24, retry=False):
    """Read only a checkpoint created within this new reference chain.

    The actual R/U method after Stable=Opt is retained. An R-to-U instability
    must not be projected back into R by a hard-coded follow-up route.
    """
    if method not in ("RPBE1PBE", "UPBE1PBE", "ROPBE1PBE"):
        raise ValueError("Unexpected method in newly generated reference checkpoint")
    if stage == "stability":
        task = "Stable=Opt"
    elif stage == "relax_after_stability":
        kind = reference["purpose"]["reference_geometry"]
        task = ("Opt=(VeryTight,CalcFC,MaxCycles=512) Freq" if kind == "optimize" else
                "SP" if kind == "single_atom" else "Force")
    else:
        raise ValueError("Unknown reference-chain stage")
    return (f"%chk=job.chk\n%mem={memory_gib}GB\n%nprocshared={cores}\n"
            f"#p {method}/ChkBasis {task} Geom=AllCheck Guess=Read "
            f"{reference['shell_flags']} EmpiricalDispersion=GD3BJ "
            f"SCF=(VeryTight,XQC,MaxCycle={2048 if retry else 1024}) "
            "Integral=SuperFineGrid NoSymm\n\n")


def prepare(args):
    destination = args.reference_output
    path = destination / "reference-candidates.json"
    if path.exists():
        raise RuntimeError("Reference candidates already frozen; preserve the existing version")
    campaign = json.loads((args.output / "campaign.json").read_text(encoding="utf-8"))
    annotations = json.loads(args.annotations.read_text(encoding="utf-8")) if args.annotations else {}
    bse = get_basis_library(args.bse)
    basis_root = destination / "basis-data"
    basis_root.mkdir(parents=True, exist_ok=True)
    all_elements = sorted({z for case in campaign["cases"] for z in case["fchk_identity"]["Atomic numbers"]})
    for name in ("def2-TZVPP", "def2-TZVPPD"):
        data = bse.get_basis(name, elements=all_elements, version=1)
        atomic_json(basis_root / (name + "-v1.json"), data)
        (basis_root / (name + "-v1.bib")).write_text(
            bse.get_references(name, elements=all_elements, version=1, fmt="bib"), encoding="utf-8")
    references = []
    for case in campaign["cases"]:
        fields = case["fchk_identity"]
        annotation = annotations.get(case["case_id"], {})
        diffuse = fields["Charge"] < 0 or annotation.get("contains_formal_anionic_fragment", False)
        name = "def2-TZVPPD" if diffuse else "def2-TZVPP"
        basis = bse.get_basis(name, elements=sorted(set(fields["Atomic numbers"])), version=1)
        core_by_z = {int(z): data.get("ecp_electrons", 0) for z, data in basis["elements"].items()}
        ecp = any(core_by_z.values())
        ne = sum(z - core_by_z[z] for z in fields["Atomic numbers"]) - fields["Charge"]
        spin_difference = fields["Multiplicity"] - 1
        if ne < spin_difference or (ne - spin_difference) % 2:
            raise ValueError(f"Electron/spin parity invalid in new reference: {case['case_id']}")
        # Preserve the existing Cartesian representation test as a distinct
        # input representation; all other reference inputs use pure d/f/g.
        cartesian = any(l >= 2 for l in fields["Shell types"])
        reference = {"case_id": case["case_id"], "original_case_identity": case["case_identity"],
                     "purpose": reference_purpose(case), "chemical_annotation": annotation,
                     "method": "PBE0-D3(BJ)", "basis_name": name, "basis_version": 1,
                     "basis_library": "Basis Set Exchange 0.12",
                     "basis_reason": "Net anion or documented anionic fragment" if diffuse else "Default quality reference",
                     "ecp_core_electrons_by_atomic_number": core_by_z,
                     "explicit_electrons": ne, "alpha_electrons": (ne + spin_difference) // 2,
                     "beta_electrons": (ne - spin_difference) // 2,
                     "gaussian_basis_keyword": "GenECP" if ecp else "Gen",
                     # 16W A.03 rejects standalone 9G/15G route keywords. Its
                     # verified 5D/7F versus 6D/10F convention also selects g.
                     "shell_flags": "6D 10F" if cartesian else "5D 7F",
                     "expected_g_shell_type": 4 if cartesian else -4,
                     "stability_required": True,
                     "competing_states_required": bool(fields.get("has_beta_coefficients") or
                         any(21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 for z in fields["Atomic numbers"])),
                     "is_verified_reference": False}
        directory = destination / "reference-inputs" / case["case_id"]
        directory.mkdir(parents=True, exist_ok=False)
        text = bse.get_basis(name, elements=sorted(set(fields["Atomic numbers"])),
                             version=1, fmt="gaussian94", header=False)
        (directory / "basis.gbs").write_text(text, encoding="ascii", newline="\n")
        (directory / "initial.gjf").write_text(make_initial_input(reference, case, text), encoding="ascii", newline="\n")
        reference["basis_sha256"] = sha256(directory / "basis.gbs")
        reference["initial_input_sha256"] = sha256(directory / "initial.gjf")
        reference["identity"] = digest(reference)
        atomic_json(directory / "reference.json", reference)
        references.append(reference)
    result = {"schema_version": 1, "created_epoch": time.time(),
              "original_case_set_identity": campaign["case_set_identity"],
              "reference_set_identity": digest([row["identity"] for row in references]),
              "candidates": references,
              "policy": {
                  "first_guess": "New atomic/Harris guess; U singlets may use fresh Guess=Mix; no original checkpoint",
                  "stages": ["fresh Opt+Freq or fixed-geometry Force or atomic SP", "Stable=Opt",
                             "if instability found: reoptimize and recompute frequencies (or Force/SP at fixed geometry), then recheck stability"],
                  "maximum_stability_repair_cycles": 2,
                  "maximum_attempts_per_stage": 2,
                  "scf_retry": "Only declared SCF convergence failures; extend XQC from 1024 to 2048 cycles",
                  "maximum_hours_per_stage_attempt": 24,
                  "geometry_iteration_limit": 512,
                  "terminal_result": "candidate_collected or explicit failure; never automatic scientific pass",
                  "physical_review": "All primary candidates collected before review; competing states, environments and other methods remain required when applicable"}}
    atomic_json(path, result)
    print(json.dumps({"cases": len(references), "basis_counts": dict(Counter(x["basis_name"] for x in references)),
                      "geometry_tasks": dict(Counter(x["purpose"]["reference_geometry"] for x in references)),
                      "with_ecp": sum(x["gaussian_basis_keyword"] == "GenECP" for x in references),
                      "path": str(path)}, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference-output", type=Path, default=DEFAULT_OUTPUT / "references" / "REF-001")
    parser.add_argument("--bse", type=Path, default=DEFAULT_BSE)
    parser.add_argument("--annotations", type=Path)
    prepare(parser.parse_args())
