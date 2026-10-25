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

"""Explicit, isolated local Gaussian acceptance. NEVER runs in hosted CI.

Run only via resume.py native-acceptance under the same host/work leases.
No formal case is calculated/relabelled, no OLD-018 checkpoint is consumed.
Calibration uses disposable distorted copies of two small frozen inputs.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path
import re
import shutil
import time


def distances(xyz):
    atoms = [xyz[i:i+3] for i in range(0, len(xyz), 3)]
    return [math.dist(a, b) for i, a in enumerate(atoms) for b in atoms[i+1:]]


def compare(a, b):
    energy = abs(a['fields']['Total Energy'] - b['fields']['Total Energy'])
    xyz = max((abs(x-y) for x, y in zip(distances(a['fields']['Current cartesian coordinates']),
                                       distances(b['fields']['Current cartesian coordinates']))), default=0.)
    freq1, freq2 = sorted(a['diagnostics']['frequencies_cm1']), sorted(b['diagnostics']['frequencies_cm1'])
    freq = max((abs(x-y) for x, y in zip(freq1, freq2)), default=0.)
    same = a['fields']['method'] == b['fields']['method'] and len(freq1) == len(freq2) and bool(freq1)
    return {'energy_delta_ha': energy, 'pair_distance_max_delta_bohr': xyz,
            'frequency_max_delta_cm1': freq, 'same_method': same,
            'passed': same and energy <= 1e-7 and xyz <= 1e-4 and freq <= .2}


class ObservingBackend:
    SUPPORTS_RAM_PAUSE = True

    def __init__(self, real, operation):
        self.real, self.operation = real, operation
        self.events = []
        self.triggered = False
        self.held_at = None
        self.hold_requested = None
        self.hold_released = None
        self.capture_started = None
        self.cpu_at_hold = None
        self.cpu_delta_while_held = None
        self.attempts = []

    def host_info(self):
        return self.real.host_info()

    def run_tree(self, argv, cwd, cores, memory, timeout, **kw):
        if Path(argv[0]).name.lower() != 'g16.exe':
            return self.real.run_tree(argv, cwd, cores, memory, timeout, **kw)
        self.held_at = self.hold_requested = self.hold_released = None
        self.cpu_at_hold = self.cpu_delta_while_held = None
        from resume import text
        first_event = len(self.events)
        t0 = time.monotonic()
        original_tick = kw.get('on_tick', lambda _: None)
        original_cancel = kw.get('cancel', lambda: False)
        def tick(value):
            self.events.append(dict(value))
            if value.get('execution_state') == 'ram_paused':
                if self.held_at is None:
                    self.held_at = time.monotonic()
                    self.cpu_at_hold = value['tree_cpu_seconds']
                self.cpu_delta_while_held = value['tree_cpu_seconds'] - self.cpu_at_hold
            original_tick(value)
        def hold():
            # Delay very briefly so this tests actual running arithmetic, not
            # merely a process left suspended at CreateProcess time.
            now = time.monotonic()
            if self.operation == 'ram_pause' and now - t0 >= .2:
                if self.hold_requested is None:
                    self.hold_requested = now
                if self.held_at is None or now - self.held_at < 2:
                    self.triggered = True
                    return True
                self.hold_released = now
            return False
        def cancel():
            if original_cancel():
                return True
            if self.operation == 'rwf_interrupt':
                log = text(Path(cwd) / 'job.log')
                # Require derivative/response work, not an empty initialized RWF.
                matches = re.findall(r'Enter[^\r\n]*[\\/]l(\d+)\.exe', log, re.I)
                rwf = Path(cwd) / 'job.rwf'
                if (any(n in ('701', '1002', '1110') for n in matches) and
                        rwf.exists() and rwf.stat().st_size > 4096 and
                        'Normal termination of Gaussian' not in log):
                    self.triggered = True
                    return True
            return False
        kw.update(on_tick=tick, hold=hold, cancel=cancel)
        result = self.real.run_tree(argv, cwd, cores, memory, timeout, **kw)
        if self.operation == 'ram_pause':
            self.attempts.append({'directory': str(cwd), 'hold_acknowledged': self.held_at is not None,
                'ack_latency_seconds': self.held_at-self.hold_requested if self.held_at is not None else None,
                'cpu_delta_during_hold': self.cpu_delta_while_held,
                'pause_diagnostics': [x['pause_diagnostic'] for x in self.events[first_event:]
                                      if x.get('execution_state') == 'pausing'],
                'pause_failures': result.get('pause_failures', []),
                'active_processes_at_return': result.get('active_processes_at_return')})
        return result


def run_acceptance(parent, host, capability='all'):
    from resume import Engine, Paused, atomic, read, sha, text
    from policy import next_stage, core_budget, digest, preferred_cores
    from fast_checkpoint import persistent_input
    root = parent.root / 'native-acceptance' / str(time.time_ns())
    root.mkdir(parents=True, exist_ok=False)
    report = {'schema': 1, 'runtime_identity': parent.identity, 'binaries': parent.data.binaries,
              'formal_cases_modified': False, 'scientific_pass': False, 'trials': {},
              'started_epoch': time.time(), 'capabilities': {}}
    cores1, cores2 = 1, min(2, core_budget(host['physical_cores']))
    if host['available_gib'] < 40:
        raise RuntimeError('Insufficient memory for isolated 24/32 GiB Gaussian acceptance job')
    evidence = []
    # Water / NO exercise restricted and unrestricted PBE0, without touching any
    # formal input or result. Exact BSE basis and numerical settings are retained.
    for case_id in ('OLD-063', 'OLD-057'):
        cores1 = preferred_cores(parent.data.sizes[case_id], physical=host['physical_cores'])
        cores2 = cores1
        trial = report['trials'][case_id] = {}
        data = copy.copy(parent.data)
        data.candidates = copy.deepcopy(parent.data.candidates)
        data.references = root / case_id / 'reference-copy'
        src = parent.data.references / 'reference-inputs' / case_id
        dst = data.references / 'reference-inputs' / case_id
        shutil.copytree(src, dst)
        original = text(dst / 'initial.gjf')
        lines = original.splitlines()
        n = parent.data.originals[case_id]['fchk_identity']['Number of atoms']
        coordinate_start = next(i+1 for i, row in enumerate(lines) if re.fullmatch(r'\s*-?\d+\s+\d+\s*', row))
        for i in range(coordinate_start, coordinate_start+n):
            fields = lines[i].split()
            lines[i] = fields[0] + ' ' + ' '.join(f'{float(v)*1.12:.14e}' for v in fields[1:])
        (dst / 'initial.gjf').write_text('\n'.join(lines)+'\n', encoding='ascii')
        trial['probe_input_sha256'] = sha(dst / 'initial.gjf')
        trial['purpose'] = 'isolated native acceptance, 12% distorted initial geometry; not a formal reference attempt'
        def create_engine(label, backend, **extra):
            config = dict(parent.config, work_root=str(root / case_id / label),
                          timeout_hours=2, disk_reserve_gib=1, opt_step_checkpoints=False,
                          freq_recovery='replay', **extra)
            engine = Engine(config, data, backend)
            engine._native_probe = True
            engine.bind(); engine.set_mode('run')
            return engine
        def phase(engine):
            engine._records_cache.pop(case_id, None)
            return next_stage(data.candidates[case_id], engine.records(case_id))
        try:
            baseline = create_engine('baseline', parent.backend)
            for _ in range(2):
                baseline.calculate(case_id, phase(baseline), cores1)
            baseline_freq = baseline.records(case_id)['freq-00']
            trial['method'] = baseline_freq['fields']['method']
            # RAM hold/resume of the same Gaussian attempt.
            observer = ObservingBackend(parent.backend, 'ram_pause')
            paused = create_engine('ram-pause', observer)
            paused.calculate(case_id, phase(paused), cores1)
            paused.calculate(case_id, phase(paused), cores1)
            ram_cmp = compare(baseline_freq, paused.records(case_id)['freq-00'])
            ram_cmp['triggered'] = observer.triggered
            ram_cmp['ack_latency_seconds'] = (observer.held_at-observer.hold_requested
                                              if observer.held_at is not None else None)
            ram_cmp['cpu_delta_during_hold'] = observer.cpu_delta_while_held
            ram_cmp['attempts'] = observer.attempts
            ram_cmp['passed'] = bool(ram_cmp['passed'] and observer.triggered and observer.held_at is not None
                                    and len(observer.attempts) == 2 and all(
                                        x['hold_acknowledged'] and not x['pause_failures'] and
                                        x['cpu_delta_during_hold'] is not None and x['cpu_delta_during_hold'] < .1 and
                                        x['active_processes_at_return'] == 0 for x in observer.attempts))
            trial['ram_pause'] = ram_cmp
            if capability == 'ram_pause':
                for path in (root / case_id).rglob('job.log'):
                    evidence.append({'path': str(path.resolve()), 'sha256': sha(path)})
                atomic(root / 'report.json', report)
                continue
            try:
                # Repeated producer-controlled l103 segments; refuse ambiguous KJob
                # producer output rather than certify a hard-kill as a safe boundary.
                segmented = create_engine('opt-segments', parent.backend)
                segmented.config['opt_step_checkpoints'] = True
                checkpointed = 0
                for segment in range(64):
                    result = segmented.calculate(case_id, phase(segmented), cores1 if segment == 0 else cores2)
                    if result == 'checkpointed':
                        checkpointed += 1
                    else:
                        break
                if phase(segmented).part != 'freq':
                    raise RuntimeError('Opt segmentation did not finish within the native trial budget')
                segmented.calculate(case_id, phase(segmented), cores2)
                opt_cmp = compare(baseline_freq, segmented.records(case_id)['freq-00'])
                opt_cmp.update(checkpointed_segments=checkpointed, cores_before=cores1, cores_after=cores2)
                opt_cmp['passed'] = bool(opt_cmp['passed'] and checkpointed >= 1)
                trial['opt_l103_segments'] = opt_cmp
            except Exception as error:
                trial['opt_l103_segments'] = {'passed': False, 'error': str(error)}
            # Copy the completed Opt receipt as an immutable donor. Interrupt only
            # a standalone analytic frequency, then invoke #p Restart from RWF.
            observer = ObservingBackend(parent.backend, 'rwf_interrupt')
            frequency = create_engine('rwf-restart', observer)
            atomic(frequency.jobs / case_id / 'opt-00/stage.json', baseline.records(case_id)['opt-00'])
            try:
                frequency.calculate(case_id, phase(frequency), cores1)
            except Paused:
                pass
            else:
                raise RuntimeError('Frequency completed before interrupt trigger; RWF recovery NOT tested')
            observer.operation = 'normal'
            frequency.config['freq_recovery'] = 'rwf'
            frequency.set_mode('run')
            frequency.calculate(case_id, phase(frequency), cores2)
            rwf_cmp = compare(baseline_freq, frequency.records(case_id)['freq-00'])
            rwf_cmp['interrupt_triggered'] = observer.triggered
            rwf_cmp['passed'] = bool(rwf_cmp['passed'] and observer.triggered)
            trial['analytic_rwf_restart'] = rwf_cmp
        except Exception as error:
            trial['error'] = f'{type(error).__name__}: {error}'
        for path in (root / case_id).rglob('job.log'):
            evidence.append({'path': str(path.resolve()), 'sha256': sha(path)})
        atomic(root / 'report.json', report)
    for capability in ('ram_pause', 'opt_l103_segments', 'analytic_rwf_restart'):
        methods = sorted({v['method'] for v in report['trials'].values()
                          if v.get(capability, {}).get('passed')})
        report['capabilities'][capability] = {'passed': bool(methods), 'methods': methods}
    report['finished_epoch'] = time.time()
    atomic(root / 'report.json', report)
    evidence.append({'path': str((root / 'report.json').resolve()), 'sha256': sha(root / 'report.json')})
    receipt = {'runtime_identity': parent.identity, 'binaries': parent.data.binaries,
               'capabilities': report['capabilities'], 'evidence': evidence}
    atomic(parent.root / 'native-capabilities.json', receipt)
    print('Native acceptance report:', root / 'report.json')
    print('Capabilities:', report['capabilities'])
    required = report['capabilities'] if capability == 'all' else ('ram_pause',)
    if any(not report['capabilities'][c]['passed'] for c in required):
        raise RuntimeError('Some native capabilities are NOT accepted; see report, do not override gate')
