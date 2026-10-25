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

"""REF-001 runtime v2: immutable attempts, checkpoint recovery, bounded scheduling.

Read README.md before migration. Never hot-patch a running legacy coordinator.
Only explicit run/import commands start local executables. No GitHub writes.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import hashlib
import importlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid

from fast_checkpoint import (cold_snapshot, restore_snapshot, persistent_input,
                             rwf_restart_input, cooperative_opt_stop)
from state_store import ControlStore, atomic_json as atomic, read_json as read, lease, occupied

from policy import (VERSION, METHODS, Stage, basis_count, core_budget, digest, followup,
                    grants, historical_reference_digest, initial_part, log_progress, next_stage, preferred_cores, route_signature)

HERE = Path(__file__).resolve().parent
RUNTIME_DIRECTORY = 'runtime-v2-unified'
FROZEN = HERE.parent / 'paused-20260906' / 'runner-source'


def runtime_directory(config):
    value = config.get('runtime_directory', RUNTIME_DIRECTORY)
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError('runtime_directory must be one directory name')
    return value


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def text(path: Path) -> str:
    return path.read_text(encoding='ascii', errors='replace') if path.is_file() else ''


class Paused(Exception):
    pass


class NeedsReview(Exception):
    pass


class FrozenData:
    def __init__(self, config: dict):
        expected = read(FROZEN / 'identity.json')['files']
        for name, checksum in expected.items():
            if sha(FROZEN / name) != checksum:
                raise ValueError('Frozen source identity changed: ' + name)
        sys.path.insert(0, str(FROZEN))
        self.rebuild = importlib.import_module('gaussian_rebuild')
        self.runner = importlib.import_module('gaussian_reference_batch')
        self.source = Path(config['data_root']).resolve()
        self.work = Path(config['work_root']).resolve()
        if self.source == self.work or self.source.is_relative_to(self.work) or self.work.is_relative_to(self.source):
            raise ValueError('Archive data and mutable work must be separate directories')
        self.gaussian = Path(config['gaussian']).resolve()
        self.references = self.work / 'references/REF-001'
        self.manifest = read(self.references / 'reference-candidates.json')
        campaign = read(self.source / 'campaign/campaign.json')
        barrier = read(self.source / 'campaign/reproduction/batch.json')
        if not (barrier['all_terminal'] and barrier['case_count'] == 273 and
                barrier['case_set_identity'] == campaign['case_set_identity'] == self.manifest['original_case_set_identity']):
            raise ValueError('Invalid original reproduction barrier or case set')
        self.candidates = {r['case_id']: r for r in self.manifest['candidates']}
        if list(self.candidates) != [f'OLD-{i:03d}' for i in range(1, 274)]:
            raise ValueError('Reference identity must contain all 273 ordered cases')
        self.originals = {r['case_id']: r for r in campaign['cases']}
        self.binaries = {name: sha(self.gaussian / name) for name in ('g16.exe', 'formchk.exe')}
        identity = self.rebuild.digest({'source_files': expected, 'reference_set': self.manifest['reference_set_identity'],
            'gaussian_exe': self.binaries['g16.exe'], 'formchk_exe': self.binaries['formchk.exe'],
            'timeout_hours': 24, 'round_name': 'REF-001'})
        self.parent_identity = read(self.work / 'jobs/REF-001/runner-identity.json')['identity']
        if identity != self.parent_identity:
            raise ValueError('Gaussian binaries or frozen dataset differ from REF-001; no silent cross-version checkpoint reuse')
        row_ids = []
        self.sizes = {}
        for case_id, ref in self.candidates.items():
            row = dict(ref)
            recorded = row.pop('identity')
            if historical_reference_digest(row) != recorded:
                raise ValueError('Reference row identity changed: ' + case_id)
            row_ids.append(recorded)
            folder = self.references / 'reference-inputs' / case_id
            if sha(folder / 'initial.gjf') != ref['initial_input_sha256'] or sha(folder / 'basis.gbs') != ref['basis_sha256']:
                raise ValueError('Frozen input/basis changed: ' + case_id)
            if ref['original_case_identity'] != self.originals[case_id]['case_identity']:
                raise ValueError('Original case identity changed: ' + case_id)
            self.sizes[case_id] = basis_count(text(folder / 'basis.gbs'),
                self.originals[case_id]['fchk_identity']['Atomic numbers'], ref['shell_flags'] == '5D 7F')
        if self.rebuild.digest(row_ids) != self.manifest['reference_set_identity']:
            raise ValueError('Reference manifest digest mismatch')
        # Fine-grained restart receipts must not survive changing a link binary
        # while keeping the same g16 launcher. Hash the installed link executables
        # and DLLs without including license files or uploading binary contents.
        for binary in sorted(self.gaussian.iterdir()):
            if binary.is_file() and (re.fullmatch(r'l\d+\.exe', binary.name, re.I) or binary.suffix.lower() == '.dll'):
                self.binaries[binary.name.lower()] = sha(binary)

    def validate(self, case_id: str, fchk: Path) -> dict:
        fields = self.rebuild.read_fchk(fchk)
        gate = self.runner.identity_gate(self.candidates[case_id], self.originals[case_id], fields)
        if not gate['consistent'] or fields['method'] not in METHODS:
            raise NeedsReview('Checkpoint identity/state mismatch: ' + repr(gate))
        if fields['Number of basis functions'] != self.sizes[case_id]:
            raise NeedsReview('Reference basis dimension mismatch; never reuse original small-basis checkpoint')
        return fields


class Engine:
    def __init__(self, config: dict, data, backend):
        self.config, self.data, self.backend = config, data, backend
        self.work = Path(config['work_root']).resolve()
        self.root = self.work / runtime_directory(config)
        self.jobs = self.root / 'jobs'
        self.control = self.root / 'control.json'
        self.active: dict[str, dict] = {}
        self.mutex = threading.Lock()
        self._records_cache: dict[str, dict] = {}
        self.identity = digest({'version': '3.0-adaptive-physical', 'parent': data.parent_identity,
                                'files': {n: sha(HERE / n) for n in ('resume.py', 'policy.py', 'windows_job.py', 'windows_pause.py', 'fast_checkpoint.py', 'native_acceptance.py', 'state_store.py', 'process_control.py', 'migration.py', 'resource_scheduler.py')},
                                'binaries': data.binaries})

    def bind(self):
        target = self.root / 'binding.json'
        current = {'runtime_identity': self.identity, 'parent_runner_identity': self.data.parent_identity,
                   'reference_set_identity': self.data.manifest['reference_set_identity'],
                   'binaries': self.data.binaries, 'version': '3.0-adaptive-physical'}
        if target.exists() and read(target) != current:
            raise ValueError('Runtime implementation/inputs changed. Review migration; do not overwrite v2 evidence')
        atomic(target, current)

    def mode(self) -> str:
        return ControlStore(self.root).read()['mode']

    def set_mode(self, mode: str):
        set_control(self.root, mode)

    def interrupt_requested(self, created_epoch: float, control_sequence=None) -> bool:
        state = read(self.control, {})
        if control_sequence is not None:
            return int(state.get('interrupt_sequence', 0)) > control_sequence
        # A subsequent Resume must not erase a save-stop already requested for
        # this active attempt. New attempts have a later creation timestamp.
        return float(state.get('interrupt_epoch', 0)) >= created_epoch

    def producer(self, case_id, phase):
        return {'runtime_identity': self.identity, 'reference_identity': self.data.candidates[case_id]['identity'],
                'binaries': self.data.binaries, 'case_id': case_id, 'stage': phase.name}

    def require_capability(self, capability, method):
        if getattr(self, '_native_probe', False):
            return  # only the isolated native acceptance harness sets this attribute
        path = self.root / 'native-capabilities.json'
        if not path.exists():
            raise NeedsReview(capability + ': run native-acceptance first; no unverified Gaussian restart')
        receipt = read(path)
        if receipt.get('runtime_identity') != self.identity or receipt.get('binaries') != self.data.binaries:
            raise NeedsReview('Native capability belongs to another runtime/Gaussian installation')
        proof = receipt.get('capabilities', {}).get(capability, {})
        if not proof.get('passed') or method not in proof.get('methods', []):
            raise NeedsReview('Native capability not accepted for ' + capability + '/' + method)
        for record in receipt.get('evidence', []):
            if sha(Path(record['path'])) != record['sha256']:
                raise NeedsReview('Native acceptance evidence changed')
        if not receipt.get('evidence'):
            raise NeedsReview('Native acceptance evidence absent')

    def records(self, case_id: str) -> dict:
        if case_id in self._records_cache:
            return self._records_cache[case_id]
        rows = {}
        for path in (self.jobs / case_id).glob('*/stage.json'):
            row = read(path)
            if row.get('runtime_identity') != self.identity:
                raise ValueError('Stage runtime identity mismatch')
            if row.get('status') == 'collected':
                for kind in ('checkpoint', 'fchk', 'log'):
                    if sha(Path(row[kind])) != row[kind + '_sha256']:
                        raise ValueError('Completed stage artifact changed: ' + str(path))
                rows[path.parent.name] = row
        self._records_cache[case_id] = rows
        return rows

    def probe(self, case_id: str, checkpoint: Path, destination: Path, cores: int = 1) -> dict:
        """Always called on a cold attempt: after tree exit or under legacy lock.

        Parse-valid is not a promise that Gaussian can restart its internal Opt
        history. A rejected Opt=Restart remains an explicit needs_review result.
        """
        destination.mkdir(parents=True, exist_ok=True)
        (destination / 'scratch').mkdir(exist_ok=True)
        copy = destination / 'job.chk'
        if checkpoint.resolve() != copy.resolve():
            before = sha(checkpoint)
            shutil.copy2(checkpoint, copy)
            if sha(copy) != before or sha(checkpoint) != before:
                raise NeedsReview('Checkpoint changed during cold copy; another writer may still exist')
        if not copy.is_file() or not copy.stat().st_size:
            raise NeedsReview('Checkpoint absent/empty; preserving attempt without guessing a restart')
        fchk = destination / 'probe.fch'
        process = self.backend.run_tree([self.data.gaussian / 'formchk.exe', copy, fchk],
                                       destination, cores, 32, 1800)
        if process['exit_code'] or process.get('interrupted'):
            raise NeedsReview('formchk rejected the interrupted checkpoint; original attempt retained')
        fields = self.data.validate(case_id, fchk)
        return {'checkpoint': str(copy), 'checkpoint_sha256': sha(copy), 'fchk': str(fchk),
                'fchk_sha256': sha(fchk), 'fields': fields, 'probe_process': process,
                'validity': 'parse_and_identity_valid; restart_history_not_certified'}

    def receipt(self, case_id: str, phase: Stage, directory: Path, provenance: dict, cores: int, *, cooperative=False) -> dict:
        log = directory / 'job.log'
        contents = text(log)
        if (not cooperative and 'Normal termination of Gaussian' not in contents) or 'Error termination' in contents:
            raise NeedsReview('Gaussian did not terminate normally')
        diagnostics = self.data.runner.stage_diagnostics(log)
        if phase.part == 'opt' and not diagnostics['optimization_completed']:
            raise NeedsReview('Normal termination without optimization completion')
        if phase.part == 'freq' and not diagnostics['frequencies_cm1']:
            raise NeedsReview('Normal termination without frequency evidence')
        chk = directory / 'job.chk'
        if not chk.exists():
            chk = directory / 'scratch/job.chk'
        probe = self.probe(case_id, chk, directory / 'validated', cores)
        result = dict(probe, name=phase.name, status='collected', runtime_identity=self.identity,
                      reference_identity=self.data.candidates[case_id]['identity'],
                      log=str(log), log_sha256=sha(log), diagnostics=diagnostics,
                      provenance=provenance, collected_epoch=time.time(), scientific_pass=False)
        atomic(directory.parent / 'stage.json', result)
        self._records_cache.pop(case_id, None)
        return result

    def interrupted_source(self, case_id: str, phase: Stage) -> dict | None:
        base = self.jobs / case_id / phase.name
        attempts = sorted(base.glob('attempt-*/attempt.json'))
        for path in reversed(attempts):
            attempt = read(path)
            if (attempt['runtime_identity'] != self.identity or
                    attempt.get('reference_identity') != self.data.candidates[case_id]['identity'] or
                    attempt.get('binaries') != self.data.binaries):
                raise NeedsReview('Interrupted attempt identity mismatch')
            if sha(path.parent / 'job.gjf') != attempt['input_sha256']:
                raise NeedsReview('Interrupted attempt input changed; checkpoint science is not trusted')
            if attempt['status'] in ('running', 'interrupted', 'timeout', 'checkpointed'):
                # A running record with no OS-held coordinator lease is historical.
                checkpoint = path.parent / 'job.chk'
                if not checkpoint.exists():
                    checkpoint = path.parent / 'scratch/job.chk'
                if checkpoint.exists():
                    if attempt.get('checkpoint_snapshot'):
                        snap = attempt['checkpoint_snapshot']
                        manifest_path = Path(snap['directory']) / 'manifest.json'
                        if sha(manifest_path) != snap['manifest_sha256']:
                            raise NeedsReview('Saved snapshot manifest changed')
                        saved = read(manifest_path)['files'].get('job.chk', {})
                        if sha(checkpoint) != saved.get('sha256'):
                            raise NeedsReview('Raw interrupted checkpoint changed after cold preservation')
                    return {'checkpoint': str(checkpoint), 'log': str(path.parent / 'job.log'),
                            'origin': str(path), 'status': attempt['status'],
                            'snapshot': attempt.get('checkpoint_snapshot')}
                # Frequency can restart its entire analytic Hessian from the
                # completed opt checkpoint. Other missing checkpoints need review.
                if phase.part != 'freq':
                    raise NeedsReview('Interrupted checkpoint missing; no silent restart from initial geometry')
                # Keep the interruption visible. Auto/RWF mode must not silently
                # become a complete frequency replay merely because CHK vanished.
                return {'checkpoint': str(checkpoint), 'log': str(path.parent / 'job.log'),
                        'origin': str(path), 'status': attempt['status'],
                        'snapshot': attempt.get('checkpoint_snapshot')}
        imported = base / 'resume-source.json'
        return read(imported) if imported.exists() else None

    def calculate(self, case_id: str, phase: Stage, cores: int) -> str:
        base = self.jobs / case_id / phase.name
        base.mkdir(parents=True, exist_ok=True)
        reference = self.data.candidates[case_id]
        for retry in (False, True):
            if self.mode() != 'run':
                raise Paused('between_stages')
            attempts = sorted(base.glob('attempt-*'))
            directory = base / f'attempt-{len(attempts)+1:04d}'
            source = self.interrupted_source(case_id, phase)
            from migration import validate_source
            source_producer = validate_source(source, self, case_id, phase) if source else self.producer(case_id, phase)
            directory.mkdir(exist_ok=False)
            (directory / 'scratch').mkdir()
            strategy = 'fresh_reference_input'
            method = None
            donor = phase.source
            if source and 'Normal termination of Gaussian' in text(Path(source['log'])):
                shutil.copy2(source['log'], directory / 'job.log')
                shutil.copy2(source['checkpoint'], directory / 'job.chk')
                self.receipt(case_id, phase, directory,
                    {'strategy': 'finished_gaussian_recovered_before_collection', 'source': source}, cores)
                return 'collected'
            rwf_restored = False
            if source and phase.part == 'freq' and self.config.get('freq_recovery', 'replay') in ('rwf', 'auto'):
                method = phase.source['fields']['method'] if phase.source else None
                self.require_capability('analytic_rwf_restart', method)
                snapshot = source.get('snapshot')
                if not snapshot:
                    raise NeedsReview('Interrupted Freq has no cold RWF snapshot; select reviewed replay explicitly')
                restore_snapshot(snapshot, directory, source_producer, require_rwf=True)
                probe = self.probe(case_id, directory / 'job.chk', directory / 'recovery', cores)
                if probe['fields']['method'] != method:
                    raise NeedsReview('RWF restart checkpoint changed R/U method')
                strategy = 'analytic_frequency_rwf_restart'
                rwf_restored = True
            if source and phase.part != 'freq':
                if source.get('snapshot'):
                    isolated = directory / 'restored-source'
                    isolated.mkdir()
                    restore_snapshot(source['snapshot'], isolated, source_producer)
                    source = dict(source, checkpoint=str(isolated / 'job.chk'))
                probe = self.probe(case_id, Path(source['checkpoint']), directory / 'recovery', cores)
                donor = probe
                strategy = 'opt_restart' if phase.part == 'opt' else 'checkpoint_guess_replay'
                # An old combined Opt+Freq may have finished Opt before interruption.
                if phase.part == 'opt' and 'Optimization completed.' in text(Path(source['log'])):
                    diagnostics = self.data.runner.stage_diagnostics(Path(source['log']))
                    shutil.copy2(source['log'], directory / 'job.log')
                    row = dict(probe, name=phase.name, status='collected', runtime_identity=self.identity,
                               reference_identity=reference['identity'], log=str(directory / 'job.log'),
                               log_sha256=sha(directory / 'job.log'), diagnostics=diagnostics,
                               provenance={'strategy': 'completed_opt_in_interrupted_combined_stage', 'source': source},
                               scientific_pass=False)
                    atomic(base / 'stage.json', row)
                    self._records_cache.pop(case_id, None)
                    return 'collected'
            elif source and phase.part == 'freq' and not rwf_restored:
                strategy = 'analytic_frequency_replay_from_completed_opt'
            if rwf_restored:
                input_text = rwf_restart_input(cores)
            elif donor:
                method = donor['fields']['method']
                checkpoint = Path(donor['checkpoint'])
                if sha(checkpoint) != donor['checkpoint_sha256']:
                    raise NeedsReview('Starting checkpoint identity changed')
                shutil.copy2(checkpoint, directory / 'job.chk')
                input_text = followup(reference, method, phase.part, cores, retry=retry,
                                      restart_opt=strategy == 'opt_restart')
                if strategy == 'fresh_reference_input':
                    strategy = 'completed_stage_checkpoint'
            else:
                if phase.part not in ('opt', 'force', 'sp') or phase.cycle != 0:
                    raise NeedsReview('Missing preceding stage')
                input_text = initial_part(text(self.data.references / 'reference-inputs' / case_id / 'initial.gjf'),
                                          phase.part, cores, retry=retry)
            preference = self.config.get('opt_step_checkpoints', False)
            segmented = phase.part == 'opt' and preference is not False
            if segmented:
                effective_method = method or re.search(r'#p\s+(\w+)/', input_text)[1]
                try:
                    self.require_capability('opt_l103_segments', effective_method)
                except NeedsReview:
                    if preference == 'auto':
                        segmented = False  # RAM pause still works; no untested KJob injection
                    else:
                        raise
                if segmented and len(attempts) >= 1024:
                    raise NeedsReview('Optimization segment budget exhausted; no unbounded restart loop')
            if not rwf_restored:
                input_text = persistent_input(input_text, opt_segments=segmented)
            (directory / 'job.gjf').write_text(input_text, encoding='ascii', newline='\n')
            attempt = {'status': 'running', 'runtime_identity': self.identity, 'reference_identity': reference['identity'],
                       'phase': phase.name, 'strategy': strategy, 'cores': cores, 'input_memory_gib': 24,
                       'scheduling_policy': VERSION,
                       'frozen_initial_input_sha256': reference.get('initial_input_sha256'),
                       'tree_memory_gib': 32, 'input_sha256': sha(directory / 'job.gjf'),
                       'route_signature': route_signature(input_text), 'scf_retry': retry,
                       'source': source, 'starting_checkpoint_sha256': sha(directory / 'job.chk') if (donor or rwf_restored) else None,
                        'binaries': self.data.binaries, 'created_epoch': time.time(),
                        'control_sequence': int(ControlStore(self.root).read().get('sequence', 0))}
            atomic(directory / 'attempt.json', attempt)
            def started(value):
                attempt['process_started'] = value
                atomic(directory / 'attempt.json', attempt)
            def tick(value):
                with self.mutex:
                    self.active[case_id] = dict(value, stage=phase.name, cores=cores, log=str(directory / 'job.log'))
            if self.mode() != 'run':
                attempt['status'] = 'paused_before_launch'
                atomic(directory / 'attempt.json', attempt)
                raise Paused('before_launch')
            process = self.backend.run_tree([self.data.gaussian / 'g16.exe', directory / 'job.gjf', directory / 'job.log'],
                directory, cores, 32, float(self.config.get('timeout_hours', 72))*3600,
                cancel=lambda: self.interrupt_requested(attempt['created_epoch'], attempt['control_sequence']), on_started=started, on_tick=tick,
                **({'hold': lambda: self.mode() == 'hold',
                    'on_hold_error': lambda message: print('RAM PAUSE FAILED (rolled back): ' + message, flush=True),
                    'on_timeout': lambda: self.set_mode('hold')}
                   if getattr(self.backend, 'SUPPORTS_RAM_PAUSE', False) else {}))
            attempt['process'] = process
            # All file capture happens AFTER run_tree reports its entire Job exited.
            if process.get('interrupted'):
                with self.mutex:
                    self.active[case_id] = {'stage': phase.name, 'cores': cores,
                        'execution_state': 'saving_cold_checkpoint', 'active_processes': 0,
                        'log': str(directory / 'job.log'), 'pause_is_disk_checkpoint': False}
                try:
                    attempt['checkpoint_snapshot'] = cold_snapshot(directory, self.producer(case_id, phase),
                        writers_exited=process.get('active_processes_at_return', 0) == 0,
                        reserve_bytes=int(self.config.get('disk_reserve_gib', 20))*2**30, include_scratch=phase.part == 'freq')
                except Exception as error:
                    attempt['checkpoint_snapshot_error'] = str(error)
                atomic(directory / 'attempt.json', attempt)
            if process.get('interrupted'):
                attempt['status'] = 'timeout' if process['stop_reason'] == 'timeout' else 'interrupted'
                try:
                    chk = directory / 'job.chk'
                    if not chk.exists():
                        chk = directory / 'scratch/job.chk'
                    attempt['recovery_probe'] = self.probe(case_id, chk, directory / 'saved', cores)
                except Exception as error:
                    attempt['recovery_probe_error'] = str(error)
                atomic(directory / 'attempt.json', attempt)
                if process['stop_reason'] == 'timeout':
                    self.set_mode('pause')
                raise Paused(process['stop_reason'])
            log = text(directory / 'job.log')
            if segmented and cooperative_opt_stop(log, input_text):
                # Expected producer-controlled exit, not a failed scientific calculation.
                # An Opt checkpoint is a resume point, not a completed candidate.
                if 'Optimization completed.' in log:
                    row = self.receipt(case_id, phase, directory,
                        {'strategy': 'cooperative_opt_complete', 'attempt': str(directory / 'attempt.json')},
                        cores, cooperative=True)
                    attempt.update(status='collected', checkpoint_sha256=row['checkpoint_sha256'])
                    atomic(directory / 'attempt.json', attempt)
                    return 'collected'
                probe = self.probe(case_id, directory / 'job.chk', directory / 'saved', cores)
                if donor and probe['checkpoint_sha256'] == donor['checkpoint_sha256']:
                    raise NeedsReview('Cooperative segment made no checkpoint progress')
                attempt['checkpoint_snapshot'] = cold_snapshot(directory, self.producer(case_id, phase),
                    writers_exited=process.get('active_processes_at_return', 0) == 0,
                    reserve_bytes=int(self.config.get('disk_reserve_gib', 20))*2**30, include_scratch=phase.part == 'freq')
                attempt.update(status='checkpointed', recovery_probe=probe,
                               note='Opt not converged; intentional L103 boundary, not a scientific pass')
                atomic(directory / 'attempt.json', attempt)
                return 'checkpointed'
            if process['exit_code'] or 'Normal termination of Gaussian' not in log or 'Error termination' in log:
                attempt['status'] = 'failed'
                atomic(directory / 'attempt.json', attempt)
                scf_failure = bool(re.search(r'Convergence failure|SCF has not converged', log, re.I))
                if not retry and scf_failure and strategy != 'opt_restart':
                    continue
                raise NeedsReview('Gaussian failure; no unbounded retry: ' + str(directory / 'job.log'))
            row = self.receipt(case_id, phase, directory, {'strategy': strategy, 'attempt': str(directory / 'attempt.json')}, cores)
            attempt.update(status='collected', checkpoint_sha256=row['checkpoint_sha256'], fchk_sha256=row['fchk_sha256'])
            atomic(directory / 'attempt.json', attempt)
            return 'collected'
        raise NeedsReview('SCF retry budget exhausted')

    def import_legacy(self) -> dict:
        """Under both locks only. Copy records/checkpoints, never rewrite v1."""
        imported, pending = 0, 0
        for case_id, reference in self.data.candidates.items():
            old = self.work / 'jobs/REF-001' / case_id
            geometry = reference['purpose']['reference_geometry']
            part = 'opt' if geometry == 'optimize' else 'sp' if geometry == 'single_atom' else 'force'
            for cycle in range(3):
                mappings = [(('initial' if cycle == 0 else f'relax-after-stability-{cycle:02d}'),
                             [f'{part}-{cycle:02d}'] + ([f'freq-{cycle:02d}'] if part == 'opt' else [])),
                            (f'stability-{cycle:02d}', [f'stability-{cycle:02d}'])]
                for old_name, names in mappings:
                    old_path = old / old_name / 'stage.json'
                    if not old_path.exists():
                        continue
                    stage = read(old_path)
                    if stage.get('runner_identity') != self.data.parent_identity or stage.get('reference_identity') != reference['identity']:
                        raise ValueError('Legacy stage identity mismatch: ' + str(old_path))
                    if stage['status'] in ('failed', 'timeout', 'identity_mismatch'):
                        atomic(self.jobs / case_id / 'review.json', {'status': 'needs_review',
                               'reason': 'legacy_' + stage['status'], 'legacy_stage': str(old_path)})
                        continue
                    if stage['status'] == 'collected':
                        if sha(Path(stage['fchk'])) != stage['fchk_sha256']:
                            raise ValueError('Legacy FCHK checksum mismatch')
                        for name in names:
                            target = self.jobs / case_id / name
                            if (target / 'stage.json').exists():
                                continue
                            folder = target / ('imported-' + uuid.uuid4().hex)
                            folder.mkdir(parents=True, exist_ok=False)
                            probe = self.probe(case_id, Path(stage['checkpoint']), folder)
                            # Match the binary checkpoint to its recorded formatted representation.
                            recorded_fields = self.data.validate(case_id, Path(stage['fchk']))
                            for key in ('method', 'Total Energy', 'Current cartesian coordinates', 'Number of basis functions'):
                                if probe['fields'][key] != recorded_fields[key]:
                                    raise NeedsReview('Legacy chk/FCHK mismatch: ' + key)
                            shutil.copy2(stage['log'], folder / 'job.log')
                            diagnostics = self.data.runner.stage_diagnostics(folder / 'job.log')
                            if not log_progress(text(folder / 'job.log'))['normal_termination']:
                                raise NeedsReview('Collected legacy stage lacks normal termination')
                            if name.startswith('opt-') and not diagnostics['optimization_completed']:
                                raise NeedsReview('Legacy optimization completion evidence absent')
                            if name.startswith('freq-') and not diagnostics['frequencies_cm1']:
                                raise NeedsReview('Legacy frequency evidence absent')
                            row = dict(probe, status='collected', name=name, runtime_identity=self.identity,
                                reference_identity=reference['identity'], diagnostics=diagnostics,
                                log=str(folder / 'job.log'), log_sha256=sha(folder / 'job.log'),
                                provenance={'legacy_stage': str(old_path), 'legacy_stage_sha256': sha(old_path),
                                            'strategy': 'verified_legacy_import'}, scientific_pass=False)
                            atomic(target / 'stage.json', row)
                            imported += 1
                    elif stage.get('attempts'):
                        # Only the first missing geometry part receives the combined-stage checkpoint.
                        target = self.jobs / case_id / names[0]
                        if (target / 'stage.json').exists() or (target / 'resume-source.json').exists():
                            continue
                        attempt = stage['attempts'][-1]
                        original = Path(attempt['directory'])
                        raw_chk = original / 'job.chk'
                        if not raw_chk.exists():
                            raw_chk = original / 'scratch/job.chk'
                        if not raw_chk.exists():
                            atomic(self.jobs / case_id / 'review.json', {'status': 'needs_review',
                                   'reason': 'legacy checkpoint missing', 'legacy_stage': str(old_path)})
                            continue
                        if sha(original / 'job.gjf') != attempt['input_sha256']:
                            raise NeedsReview('Legacy interrupted input checksum mismatch')
                        folder = target / ('legacy-interrupted-' + uuid.uuid4().hex)
                        folder.mkdir(parents=True, exist_ok=False)
                        probe = self.probe(case_id, raw_chk, folder)
                        shutil.copy2(original / 'job.log', folder / 'job.log')
                        atomic(target / 'resume-source.json', dict(probe, log=str(folder / 'job.log'),
                               legacy_stage=str(old_path), legacy_stage_sha256=sha(old_path), status='interrupted'))
                        pending += 1
        summary = {'imported_parts': imported, 'interrupted_sources': pending,
                   'legacy_rewritten': False, 'created_epoch': time.time(), 'scientific_passes': 0}
        atomic(self.root / 'import-summary.json', summary)
        return summary

    def calculate_reserved(self, case_id, phase, cores, reservation):
        if reservation is None:  # Test doubles do not start native executables.
            return self.calculate(case_id, phase, cores)
        try:
            with reservation.scope():
                return self.calculate(case_id, phase, cores)
        finally:
            reservation.release()

    def run(self, selected: list[str], host: dict):
        budget = core_budget(host['physical_cores'], self.config.get('core_budget'))
        memory_budget = min(128, int(host['total_gib']) - 16)
        if memory_budget < 32:
            raise ValueError('At least 48 GiB visible system memory is required for the retained 24/32 GiB job envelope')
        slots = min(int(self.config.get('workers', 2)), 2, budget, memory_budget // 32)
        self.set_mode('run')
        state = {'runtime_identity': self.identity, 'selected_cases': selected,
                 'physical_core_budget': budget, 'host': host, 'max_parallel_jobs': slots,
                 'cases': {}, 'active': {}, 'started_epoch': time.time(), 'scientific_passes': 0,
                 'coordinator_running': True}
        running = {}
        order = list(selected)
        excluded = set()
        previous_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, lambda *_: self.set_mode('interrupt'))
        pool = ThreadPoolExecutor(max_workers=slots)
        try:
            while True:
                for case_id, item in list(running.items()):
                    future, cores, phase = item
                    if not future.done():
                        continue
                    try:
                        outcome = future.result()
                        self._records_cache.pop(case_id, None)
                        state['cases'][case_id] = 'checkpointed' if outcome == 'checkpointed' else 'stage_collected'
                        print(f'{case_id} {phase.name}: {outcome}', flush=True)
                    except Paused as error:
                        state['cases'][case_id] = 'paused_' + str(error)
                        print(f'{case_id} {phase.name}: paused ({error})', flush=True)
                    except Exception as error:
                        excluded.add(case_id)
                        state['cases'][case_id] = 'needs_review'
                        print(f'{case_id} {phase.name}: needs_review: {error}', flush=True)
                        atomic(self.jobs / case_id / 'review.json', {'status': 'needs_review',
                               'stage': phase.name, 'reason': str(error), 'epoch': time.time()})
                    finally:
                        with self.mutex:
                            self.active.pop(case_id, None)
                        del running[case_id]
                        order.remove(case_id)
                        order.append(case_id)
                ready = []
                for case_id in order:
                    if case_id in running or case_id in excluded:
                        continue
                    if (self.jobs / case_id / 'review.json').exists():
                        state['cases'][case_id] = 'needs_review'
                        excluded.add(case_id)
                        continue
                    phase = next_stage(self.data.candidates[case_id], self.records(case_id))
                    if isinstance(phase, str):
                        state['cases'][case_id] = phase
                        excluded.add(case_id)
                        if phase != 'candidate_collected':
                            atomic(self.jobs / case_id / 'review.json', {'status': 'needs_review', 'reason': phase})
                        else:
                            atomic(self.jobs / case_id / 'result.json', {'status': phase, 'runtime_identity': self.identity,
                                  'case_id': case_id, 'scientific_pass': False,
                                  'competing_states_pending': self.data.candidates[case_id]['competing_states_required'],
                                  'stages': sorted(self.records(case_id)), 'finished_epoch': time.time()})
                        continue
                    ready.append((case_id, phase))
                mode = self.mode()
                if mode == 'run' and ready:
                    free = budget - sum(item[1] for item in running.values())
                    candidates = ready[:max(0, slots-len(running))]
                    asks = [preferred_cores(self.data.sizes[c],
                            self.data.originals[c]['fchk_identity'].get('has_beta_coefficients', False) or
                            self.data.candidates[c]['alpha_electrons'] != self.data.candidates[c]['beta_electrons'],
                            host['physical_cores'])
                            for c, _ in candidates]
                    if any(n > budget for n in asks):
                        raise ValueError('Configured budget cannot fit a whole workload tier; increase it before running')
                    allocation = grants(asks, free)
                    available_memory = self.backend.host_info()['available_gib']
                    for (case_id, phase), cores in zip(candidates, allocation):
                        if cores == 0:
                            continue
                        if available_memory < 40:
                            state['admission_warning'] = 'Waiting for 32 GiB job envelope plus 8 GiB free-memory reserve'
                            continue
                        if shutil.disk_usage(self.work).free < int(self.config.get('disk_reserve_gib', 20))*2**30:
                            self.set_mode('pause')
                            state['admission_warning'] = 'Low disk space; no new job dispatched'
                            break
                        available_memory -= 32
                        reservation = None
                        if hasattr(self.backend, 'new_reservation'):
                            from resource_scheduler import AdmissionDenied
                            try:
                                reservation = self.backend.new_reservation(cores,
                                    label=case_id + '/' + phase.name, budget=budget,
                                    metadata={'runtime_identity': self.identity,
                                        'reference_identity': self.data.candidates[case_id]['identity'],
                                        'nbasis': self.data.sizes[case_id], 'phase': phase.name})
                            except AdmissionDenied as error:
                                state['admission_warning'] = str(error)
                                break  # FIFO: never bypass a blocked heavy head.
                        with self.mutex:
                            self.active[case_id] = {'stage': phase.name, 'cores': cores, 'status': 'launching'}
                        print(f'{case_id} {phase.name}: starting with {cores} cores', flush=True)
                        try:
                            future = pool.submit(self.calculate_reserved, case_id, phase, cores, reservation)
                        except BaseException:
                            if reservation:
                                reservation.release()
                            raise
                        running[case_id] = (future, cores, phase)
                        state['cases'][case_id] = 'running'
                with self.mutex:
                    state['active'] = dict(self.active)
                state.update(mode=self.mode(), heartbeat_epoch=time.time(), coordinator_pid=os.getpid(),
                             counts=dict(Counter(state['cases'].values())),
                             no_active_writer_trees=not bool(running),
                             can_close_coordinator_without_killing_jobs=not bool(running),
                             ram_pause_warning='RAM pause is not a disk save; do not power off or close coordinator')
                atomic(self.root / 'status.json', state)
                if (mode == 'shutdown' and not running) or (not ready and not running):
                    break
                time.sleep(getattr(self, 'poll_seconds', .25))
        except BaseException:
            # Signal cancellation BEFORE executor.shutdown waits for its futures.
            self.set_mode('interrupt')
            raise
        finally:
            pool.shutdown(wait=True)
            signal.signal(signal.SIGINT, previous_handler)
            state.update(finished_epoch=time.time(), active={}, coordinator_running=False)
            atomic(self.root / 'status.json', state)

def set_control(root: Path, mode: str):
    return ControlStore(root).set(mode)


def report(work: Path, directory=RUNTIME_DIRECTORY):
    root = work / directory
    status = read(root / 'status.json', {})
    status['work_lease_held'] = occupied(work / 'jobs/supervisor.lock')
    status['control'] = read(root / 'control.json', None)
    # A coordinator's last selection is not the full 273-case census.
    counts = Counter()
    for i in range(1, 274):
        folder = root / 'jobs' / f'OLD-{i:03d}'
        if (folder / 'result.json').is_file():
            counts['candidate_collected'] += 1
        elif (folder / 'review.json').is_file():
            counts['needs_review'] += 1
        elif any(folder.glob('*/stage.json')) or any(folder.glob('*/attempt-*')) or any(folder.glob('*/resume-source.json')):
            counts['partial'] += 1
        else:
            counts['not_started'] += 1
    status['last_selection_counts'] = status.get('counts', {})
    status['counts'] = dict(counts)
    status['scheduling_policy'] = VERSION
    if status.get('heartbeat_epoch'):
        status['heartbeat_age_seconds'] = time.time() - status['heartbeat_epoch']
    status['live_note'] = 'History is not live process state. Check OS lease, heartbeat age and owned Job accounting together.'
    print(json.dumps(status, indent=2, ensure_ascii=False))
    for case_id, active in status.get('active', {}).items():
        log = Path(active.get('log', ''))
        if log.is_file():
            print(case_id, json.dumps(log_progress(text(log)), indent=2, ensure_ascii=False))


def menu(config_path: Path, config: dict):
    work = Path(config['work_root'])
    while True:
        print('\nCOV REF-001 统一运行时 2.2\n'
              '1 检查数据及物理核预算（不算）\n'
              '2 导入 2.1 结果和检查点（旧协调器必须已退出）\n'
              '3 开始队列（自动分配核心）\n'
              '4 快速内存暂停（保留现场，不是落盘保存；不可关机/关闭计算窗口）\n'
              '5 查看状态、核心分配、最近三组几何收敛表\n'
              '6 撤销暂停，恢复全部可用计算槽位\n'
              '8 停止并冷保存 CHK/RWF（会丢失未保存工作；等保存回执后才能关机）\n'
              '9 当前阶段收尾后退出协调器\n'
              '10 优化检查点/子阶段边界暂停（比内存暂停慢，但可落盘）\n'
              '0 仅退出菜单')
        choice = input('> ').strip()
        try:
            if choice == '0':
                return
            if choice in ('1', '2', '3'):
                action = {'1': 'check', '2': 'import-runtime', '3': 'run'}[choice]
                argv = [sys.executable, '-X', 'utf8', str(HERE / 'resume.py'), '--config', str(config_path), action]
                if choice == '3':
                    limit = int(input('最多选取多少个未完成案例 [273]：').strip() or '273')
                    argv += ['--limit', str(limit)]
                    subprocess.Popen(argv, creationflags=subprocess.CREATE_NEW_CONSOLE)
                else:
                    subprocess.run(argv, check=True)
            elif choice in ('4', '6', '8', '9', '10'):
                mode = {'4': 'hold', '6': 'run', '8': 'interrupt', '9': 'shutdown', '10': 'pause'}[choice]
                if choice == '8' and input('确认中断并尝试保存？输入 SAVE：').strip() != 'SAVE':
                    continue
                set_control(work / runtime_directory(config), mode)
                print('v2 control written:', mode, '(no effect on an old v1 coordinator)')
                if choice == '6' and not occupied(work / 'jobs/supervisor.lock'):
                    print('No coordinator owns this work directory. Use 3 to start it.')
            elif choice == '5':
                report(work, runtime_directory(config))
        except (Exception, KeyboardInterrupt) as error:
            print('ERROR:', error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=HERE / 'config.json')
    parser.add_argument('command', choices=('menu', 'check', 'import-legacy', 'import-runtime', 'run', 'pause', 'hold', 'resume', 'interrupt', 'shutdown', 'status', 'retry-reviewed', 'native-acceptance'))
    parser.add_argument('--source', type=Path)
    parser.add_argument('--capability', choices=('all', 'ram_pause'), default='all')
    parser.add_argument('--case', action='append')
    parser.add_argument('--limit', type=int, default=273)
    parser.add_argument('--reason')
    parser.add_argument('--workers', type=int, choices=(1, 2))
    args = parser.parse_args()
    config = read(args.config)
    if args.workers is not None:
        config['workers'] = args.workers
    if config.get('workers', 2) not in (1, 2):
        raise ValueError('workers must be 1 or 2')
    if config.get('freq_recovery', 'replay') not in ('replay', 'rwf', 'auto'):
        raise ValueError('freq_recovery must be replay, rwf or auto')
    if config.get('opt_step_checkpoints', False) != 'auto' and not isinstance(config.get('opt_step_checkpoints', False), bool):
        raise ValueError('opt_step_checkpoints must be boolean or auto')
    import math
    hours = float(config.get('timeout_hours', 72))
    if not math.isfinite(hours) or not 0 < hours <= 168:
        raise ValueError('timeout_hours must be finite and between 0 and 168')
    if not 1 <= int(config.get('disk_reserve_gib', 20)) <= 10000:
        raise ValueError('Invalid disk reserve')
    work = Path(config['work_root']).resolve()
    if args.command == 'menu':
        return menu(args.config.resolve(), config)
    if args.command in ('pause', 'hold', 'resume', 'interrupt', 'shutdown'):
        set_control(work / runtime_directory(config), {'resume': 'run'}.get(args.command, args.command))
        print('Command written; only a v2 coordinator reacts. Legacy v1 is not hot-patched.')
        return
    if args.command == 'status':
        return report(work, runtime_directory(config))
    data = FrozenData(config)
    import windows_job as backend
    host = backend.host_info()
    budget = core_budget(host['physical_cores'], config.get('core_budget'))
    if args.command == 'check':
        print(json.dumps({'host': host, 'core_budget': budget, 'verified_inputs': len(data.candidates),
                         'basis_counts': data.sizes, 'calculation_started': False}, indent=2))
        return
    if not 1 <= args.limit <= 273:
        raise ValueError('limit must be between 1 and 273')
    selected = list(dict.fromkeys(args.case or list(data.candidates)))
    if any(c not in data.candidates for c in selected):
        raise ValueError('Unknown case ID')
    engine = Engine(config, data, backend)
    global_lock = Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'COV/gaussian-runtime-v2.lock'
    with ExitStack() as stack:
        stack.enter_context(lease(global_lock))
        stack.enter_context(lease(work / 'jobs/supervisor.lock'))
        foreign = backend.existing_gaussian(data.gaussian)
        if foreign:
            raise RuntimeError('Foreign/legacy Gaussian still exists; refusing admission: ' + repr(foreign))
        backend.full_startup_affinity()
        engine.bind()
        if args.command == 'native-acceptance':
            from native_acceptance import run_acceptance
            return run_acceptance(engine, host, capability=args.capability)
        if args.command == 'retry-reviewed':
            if not args.case or not args.reason or not args.reason.strip():
                raise ValueError('retry-reviewed requires explicit --case and a nonempty --reason')
            for case_id in selected:
                review = engine.jobs / case_id / 'review.json'
                if not review.exists():
                    raise ValueError('No review hold for ' + case_id)
                archived = review.with_name('review-' + str(time.time_ns()) + '.json')
                record = read(review)
                atomic(archived, dict(record, operator_retry_reason=args.reason, approved_epoch=time.time()))
                review.unlink()
            print('Review holds archived. No calculation started; use run explicitly.')
            return
        if args.command == 'import-legacy':
            print(json.dumps(engine.import_legacy(), indent=2))
            return
        if args.command == 'import-runtime':
            from migration import import_runtime
            result = import_runtime(engine, args.source or work / 'runtime-v2-fastpause')
            print(json.dumps({k: v for k, v in result.items() if k != 'source_evidence'}, indent=2))
            return
        if not (engine.root / 'import-summary.json').exists():
            raise RuntimeError('Import idle runtime evidence with import-runtime (or import-legacy for v1) first')
        selected = [c for c in selected if next_stage(data.candidates[c], engine.records(c)) != 'candidate_collected'][:args.limit]
        if not selected:
            print('No uncollected cases in selection; no calculation started.')
            return
        # Prevent unattended Windows sleep, not display-off. Restore on exit.
        k, _ = backend.api()
        k.SetThreadExecutionState.argtypes = [ctypes_dword()]
        if not k.SetThreadExecutionState(0x80000001):
            raise RuntimeError('Failed to inhibit system sleep')
        try:
            engine.run(selected, host)
        finally:
            k.SetThreadExecutionState(0x80000000)


def ctypes_dword():
    from ctypes import wintypes
    return wintypes.DWORD


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {type(error).__name__}: {error}', file=sys.stderr)
        sys.exit(1)
