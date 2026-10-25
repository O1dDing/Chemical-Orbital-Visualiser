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

"""Pure, testable scheduling/input/progress policy. Never starts Gaussian.

CPU counts are physical-core-derived budgets, not a claim of exclusive pinning.
Gaussian 16W A.03 cannot inherit restricted affinity on the validated machine.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re

VERSION = '3.0-adaptive-physical'
METHODS = ('RPBE1PBE', 'UPBE1PBE', 'ROPBE1PBE')
TASK = 'Opt=(VeryTight,CalcFC,MaxCycles=512)'
ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn '
            'Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd '
            'In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
            'Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu '
            'Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def historical_reference_digest(reference: dict) -> str:
    """REF-001 hashed integer atomic-number keys before JSON stringified them."""
    row = dict(reference)
    row.pop('identity', None)
    row['ecp_core_electrons_by_atomic_number'] = {
        int(k): v for k, v in row['ecp_core_electrons_by_atomic_number'].items()}
    return digest(row)


def core_budget(physical: int, requested: int | None = None) -> int:
    if isinstance(physical, bool) or not isinstance(physical, int) or physical < 1:
        raise ValueError('Physical core count must be a positive integer, not logical threads')
    # 16 and 8 are reviewed host profiles. Other hosts retain at least one
    # physical core (two for >= 9); never infer capacity from SMT threads.
    reserve = 2 if physical >= 9 else (1 if physical >= 2 else 0)
    ceiling = min(14, max(1, physical - reserve))
    if requested is not None and (isinstance(requested, bool) or not isinstance(requested, int) or not 1 <= requested <= ceiling):
        raise ValueError(f'CPU budget must be between 1 and {ceiling} on this host')
    return ceiling if requested is None else requested


def preferred_cores(nbasis: int, open_shell: bool = False, physical: int = 16) -> int:
    if isinstance(nbasis, bool) or not isinstance(nbasis, int) or nbasis < 1:
        raise ValueError('Basis-function count unavailable; do not guess from old small-basis FCHK')
    base = 4 if nbasis < 300 else 5 if nbasis < 450 else 7
    budget = core_budget(physical)
    # All REF cases already use hybrid exchange and strict numerical settings.
    # Open-shell/convergence flags remain scientific provenance, not arbitrary
    # overrides to the common dimension tiers or a measured speedup claim.
    return min(budget, max(1, math.ceil(base * budget / 14)))


def grants(requests: list[int], available: int) -> list[int]:
    """FIFO whole-request admission. An unfittable head blocks later requests.

    Never downgrade a heavy request. The coordinator rotates a finished phase
    behind waiting cases, preventing short stages from repeatedly cutting in.
    """
    if available < 0 or any(n < 1 or n > 14 for n in requests):
        raise ValueError('Invalid resource request')
    result = [0] * len(requests)
    for i, n in enumerate(requests):
        if n > available:
            break
        result[i] = n
        available -= n
    return result


def basis_count(basis: str, numbers: list[int], pure: bool) -> int:
    """Count contracted Gaussian94 shells, not primitives; ignore the ECP block.

    This strict parser supports the element-labelled BSE files frozen by REF-001.
    Other syntax is rejected instead of producing a misleading resource estimate.
    """
    counts: dict[str, int] = {}
    for block in basis.split('****'):
        lines = [x.strip() for x in block.splitlines() if x.strip() and not x.lstrip().startswith('!')]
        if not lines:
            continue
        header = re.fullmatch(r'([A-Za-z]{1,2})\s+0', lines[0])
        if not header:
            raise ValueError('Unsupported Gaussian94 basis block')
        symbol = header[1].capitalize()
        if symbol not in ELEMENTS:
            raise ValueError('Unknown basis element')
        # ECP sections follow the last **** and have e.g. Re-ECP 4 60.
        if len(lines) > 1 and not re.fullmatch(r'(?:[SPDFGHI]|SP)\s+\d+\s+\S+', lines[1], re.I):
            continue
        if symbol in counts:
            raise ValueError('Duplicate orbital basis block')
        count, i = 0, 1
        while i < len(lines):
            match = re.fullmatch(r'(SP|[SPDFGHI])\s+(\d+)\s+\S+', lines[i], re.I)
            if not match:
                raise ValueError('Unsupported basis shell header: ' + lines[i])
            label, primitives = match[1].upper(), int(match[2])
            if primitives < 1 or i + primitives >= len(lines):
                raise ValueError('Truncated orbital basis')
            for row in lines[i+1:i+1+primitives]:
                values = row.replace('D', 'E').replace('d', 'e').split()
                if len(values) != (3 if label == 'SP' else 2):
                    raise ValueError('Unsupported generalized contraction')
                if not all(math.isfinite(float(v)) for v in values):
                    raise ValueError('Nonfinite basis entry')
            l = 'SPDFGHI'.find(label)
            count += 4 if label == 'SP' else (2*l+1 if pure else (l+1)*(l+2)//2)
            i += primitives + 1
        counts[symbol] = count
    if not numbers or any(z < 1 or z >= len(ELEMENTS) or ELEMENTS[z] not in counts for z in numbers):
        raise ValueError('Orbital basis missing an atom element')
    return sum(counts[ELEMENTS[z]] for z in numbers)


def resource_header(cores: int, memory: int = 24) -> str:
    if not 1 <= cores <= 14 or memory < 1:
        raise ValueError('Invalid Link 0 resources')
    return f'%chk=job.chk\n%mem={memory}GB\n%nprocshared={cores}\n'


def initial_part(text: str, part: str, cores: int, memory: int = 24, retry: bool = False) -> str:
    """Split only the known frozen route; do not loosen scientific thresholds."""
    if '--link1--' in text.lower() or text.lower().count('#p ') != 1:
        raise ValueError('Only the frozen single-job REF-001 input is supported')
    expected = '%chk=job.chk\n%mem=24GB\n%nprocshared=4\n'
    text = text.replace('\r\n', '\n')
    if not text.startswith(expected):
        raise ValueError('Unexpected frozen Link 0 input')
    result = resource_header(cores, memory) + text[len(expected):]
    if part == 'opt':
        if TASK + ' Freq' not in result:
            raise ValueError('Optimization route does not match frozen protocol')
        result = result.replace(TASK + ' Freq', TASK, 1)
    elif part not in ('sp', 'force'):
        raise ValueError('Invalid initial part')
    if retry:
        result = result.replace('MaxCycle=1024', 'MaxCycle=2048')
    return result


def followup(reference: dict, method: str, part: str, cores: int,
             memory: int = 24, retry: bool = False, restart_opt: bool = False) -> str:
    if method not in METHODS or part not in ('opt', 'freq', 'stability', 'sp', 'force'):
        raise ValueError('Unsupported method or stage; never project a U solution back into R')
    header = resource_header(cores, memory)
    if restart_opt:
        if part != 'opt':
            raise ValueError('Opt=Restart cannot restart a frequency or stability calculation')
        # Model chemistry and optimization state come from the validated checkpoint.
        # Do not invent Freq=Restart for an analytic frequency calculation.
        return header + '#p Opt=(Restart,VeryTight,MaxCycles=512)\n\n'
    task = {'opt': TASK, 'freq': 'Freq', 'stability': 'Stable=Opt', 'sp': 'SP', 'force': 'Force'}[part]
    return (header + f'#p {method}/ChkBasis {task} Geom=AllCheck Guess=Read '
            f"{reference['shell_flags']} EmpiricalDispersion=GD3BJ "
            f'SCF=(VeryTight,XQC,MaxCycle={2048 if retry else 1024}) '
            'Integral=SuperFineGrid NoSymm\n\n')


def route_signature(text: str) -> str:
    """Separate science/restart intent from memory/thread/path changes."""
    return digest('\n'.join(x for x in text.replace('\r\n', '\n').splitlines()
                            if not x.lstrip().startswith('%')))


FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][-+]?\d+)?'


def log_progress(text: str) -> dict:
    """Report observations, never a fictitious SCF/geometry percent or ETA."""
    result: dict = {'optimization_completed': 'Optimization completed.' in text,
                    'normal_termination': 'Normal termination of Gaussian' in text,
                    'error_termination': 'Error termination' in text}
    steps = re.findall(r'Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', text)
    if steps:
        result['optimization_step'] = int(steps[-1][0])
        result['maximum_steps'] = int(steps[-1][1])
    basis = re.findall(r'NBasis=\s*(\d+)', text)
    if basis:
        result['nbasis'] = int(basis[-1])
    scf = re.findall(r'Cycle\s+(\d+)\s+Pass\s+\d+', text)
    result['last_scf_cycle'] = int(scf[-1]) if scf else None
    tables = []
    rows = re.findall(r'(Maximum\s+Force|RMS\s+Force|Maximum\s+Displacement|RMS\s+Displacement)'
                      r'\s+(' + FLOAT + r')\s+(' + FLOAT + r')\s+(YES|NO)', text)
    for label, value, threshold, passed in rows:
        label = ' '.join(label.split())
        if label == 'Maximum Force':
            tables.append([])
        if tables:
            v, t = (float(x.replace('D', 'E').replace('d', 'e')) for x in (value, threshold))
            tables[-1].append({'item': label, 'value': v, 'printed_threshold': t,
                               'ratio_to_printed_threshold': v/t if t else None, 'passed': passed == 'YES'})
    complete = [t for t in tables if len(t) == 4]
    result['last_convergence_tables'] = complete[-3:]
    result['convergence_note'] = 'Ratios use rounded printed thresholds; four tests are not a time percentage.'
    result['phase_observed'] = ('finished' if result['normal_termination'] else
                                'post_optimization' if result['optimization_completed'] else
                                'optimization' if steps or rows else 'SCF_or_initialization')
    result['eta_seconds'] = None
    return result


@dataclass(frozen=True)
class Stage:
    name: str
    part: str
    cycle: int
    source: dict | None = None


def next_stage(reference: dict, records: dict[str, dict]) -> Stage | str:
    """Same two-repair-cycle scientific chain as REF-001, with Opt/Freq split."""
    geometry = reference['purpose']['reference_geometry']
    part = 'opt' if geometry == 'optimize' else 'sp' if geometry == 'single_atom' else 'force'
    source = None
    for cycle in range(3):
        name = f'{part}-{cycle:02d}'
        if name not in records:
            return Stage(name, part, cycle, source)
        current = records[name]
        if part == 'opt':
            name = f'freq-{cycle:02d}'
            if name not in records:
                return Stage(name, 'freq', cycle, current)
            current = records[name]
        name = f'stability-{cycle:02d}'
        if name not in records:
            return Stage(name, 'stability', cycle, current)
        stable = records[name]
        before, after = current['fields'], stable['fields']
        relax = (stable['diagnostics']['instability_encountered'] or before['method'] != after['method'] or
                 abs(before['Total Energy'] - after['Total Energy']) > 1e-8)
        if not relax:
            return 'candidate_collected' if stable['diagnostics']['stability_reported'] else 'missing_stability_evidence'
        source = stable
    return 'stability_repair_cycles_exhausted'
