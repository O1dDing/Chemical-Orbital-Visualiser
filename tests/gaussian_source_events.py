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

"""Extract ordered source-log evidence, without claiming physical validity.

No molecule names, orbital indices, or fitted thresholds enter this parser.
The resulting event timeline is suitable for review only after a frozen batch
has finished. Equality of printed SCF energies does not prove equal densities.
"""
from __future__ import annotations

import math
import re
from collections.abc import Iterable


_REAL = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?'
_SCF = re.compile(r'^\s*SCF Done:\s+E\(([^)]+)\)\s*=\s*(\S+)')
_FREQUENCY = re.compile(r'^\s*Frequencies\s+--\s*(.*?)\s*$')
_SPIN = re.compile(r'^\s*S\*\*2 before annihilation\s+(' + _REAL +
                   r'),?\s+after\s+(' + _REAL + r')\s*$')


def _real(token: str) -> float:
    if re.fullmatch(_REAL, token) is None:
        raise ValueError(f'Invalid real token in Gaussian evidence: {token!r}')
    value = float(token.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise ValueError('Non-finite Gaussian evidence')
    return value


def extract_events(lines: Iterable[str]) -> dict:
    events: list[dict] = []
    total_lines = 0
    for total_lines, line in enumerate(lines, 1):
        text = line.strip()
        scf = _SCF.match(line)
        frequency = _FREQUENCY.match(line)
        spin = _SPIN.match(line)
        event = None
        if scf:
            event = {'kind': 'scf', 'method_label': scf.group(1),
                     'energy_hartree': _real(scf.group(2))}
        elif frequency:
            tokens = frequency.group(1).split()
            if not tokens:
                raise ValueError('Empty frequency row')
            event = {'kind': 'frequencies', 'values_cm_inverse': [_real(t) for t in tokens]}
        elif spin:
            event = {'kind': 'spin_expectation', 'before_annihilation': _real(spin.group(1)),
                     'after_annihilation': _real(spin.group(2))}
        elif text == 'Optimization completed.':
            # Gaussian can also print this during frequency bookkeeping.
            # It is recorded literally, never promoted to proof of a new optimization.
            event = {'kind': 'optimization_complete_message'}
        elif text.startswith('The wavefunction has ') and 'instability' in text.lower():
            event = {'kind': 'instability_report'}
        elif text == 'The wavefunction is stable under the perturbations considered.':
            event = {'kind': 'stable_report'}
        elif text == 'The wavefunction is already stable.':
            event = {'kind': 'already_stable_report'}
        elif text.startswith('Normal termination of Gaussian '):
            event = {'kind': 'normal_termination'}
        elif text.startswith('Error termination '):
            event = {'kind': 'error_termination'}
        if event is not None:
            events.append({'line': total_lines, **event, 'text': text})
    return {'line_count': total_lines, 'events': events,
            'scope': 'Printed event provenance only; not a physical or COV acceptance verdict'}


def review_final_frequency_chain(extracted: dict, energy_change_hartree: float) -> dict:
    """Identify explicit changes after the last frequency, preserving uncertainty.

    This narrow check is applied only to matching equilibrium source logs by
    the batch reviewer. Missing printed events are evidence gaps, not zeros.
    """
    if not math.isfinite(energy_change_hartree) or energy_change_hartree <= 0:
        raise ValueError('A positive finite frozen energy-change threshold is required')
    events = extracted['events']
    scfs = [e for e in events if e['kind'] == 'scf']
    frequencies = [e for e in events if e['kind'] == 'frequencies']
    if not scfs:
        return {'status': 'missing_scf_evidence', 'physical_pass': False}
    if not frequencies:
        return {'status': 'missing_frequency_evidence', 'final_scf': scfs[-1],
                'physical_pass': False}
    last_frequency = frequencies[-1]
    before = [e for e in scfs if e['line'] < last_frequency['line']]
    if not before:
        return {'status': 'frequency_state_not_identified', 'last_frequency': last_frequency,
                'final_scf': scfs[-1], 'physical_pass': False}
    frequency_scf, final_scf = before[-1], scfs[-1]
    delta = final_scf['energy_hartree'] - frequency_scf['energy_hartree']
    label_changed = final_scf['method_label'] != frequency_scf['method_label']
    energy_changed = abs(delta) > energy_change_hartree
    after = [e for e in events if e['line'] > last_frequency['line']]
    later_instabilities = [e for e in after if e['kind'] == 'instability_report']
    status = 'no_explicit_post_frequency_state_change_at_log_precision'
    if energy_changed or label_changed:
        status = ('post_frequency_instability_and_scf_state_change' if later_instabilities
                  else 'post_frequency_scf_state_change_needs_context')
    return {'status': status, 'frequency_scf': frequency_scf, 'last_frequency': last_frequency,
            'final_scf': final_scf, 'energy_delta_hartree': delta,
            'method_label_changed': label_changed, 'energy_change_exceeds_frozen_threshold': energy_changed,
            'later_instability_reports': later_instabilities,
            'later_stable_reports': [e for e in after if e['kind'] in ('stable_report', 'already_stable_report')],
            'physical_pass': False,
            'limitations': ['Printed energy equality does not prove wavefunction equality',
                            'Optimization messages alone do not establish a final-state stationary point',
                            'One log may omit evidence saved in a different calculation']}
