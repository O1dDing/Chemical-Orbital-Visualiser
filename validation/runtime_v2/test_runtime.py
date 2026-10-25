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

"""Policy/scheduler tests use a fake Gaussian. Windows Job tests use Python only.

Neither test suite is evidence of Gaussian scientific/native restart acceptance.
"""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import threading
import time
import unittest
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import policy as p
import resume as r


REFERENCE = {'case_id': 'OLD-018', 'identity': 'science', 'purpose': {'reference_geometry': 'optimize'},
             'shell_flags': '5D 7F', 'alpha_electrons': 42, 'beta_electrons': 39,
             'competing_states_required': True}
INPUT = ('%chk=job.chk\n%mem=24GB\n%nprocshared=4\n'
         '#p UPBE1PBE/Gen Opt=(VeryTight,CalcFC,MaxCycles=512) Freq Units=Bohr 5D 7F '
         'EmpiricalDispersion=GD3BJ SCF=(VeryTight,XQC,MaxCycle=1024) Integral=SuperFineGrid NoSymm\n\n'
         'test\n\n3 4\n24 0 0 0\n\nCr 0\nS 1 1.0\n1.0 1.0\n****\n\n')
FIELDS = {'method': 'UPBE1PBE', 'Total Energy': -1.0, 'Current cartesian coordinates': [0., 0., 0.],
          'Number of basis functions': 502}


def diag(log):
    value = r.text(log)
    return {'optimization_completed': 'Optimization completed.' in value,
            'frequencies_cm1': [100.0] if 'Frequencies --' in value else [],
            'stability_reported': 'wavefunction is stable under' in value,
            'instability_encountered': 'wavefunction has an instability' in value}


def record(unstable=False, energy=-1., method='UPBE1PBE'):
    return {'fields': dict(FIELDS, method=method, **{'Total Energy': energy}),
            'diagnostics': {'stability_reported': not unstable, 'instability_encountered': unstable}}


class PolicyTests(unittest.TestCase):
    def test_physical_budgets(self):
        for physical, expected in ((1, 1), (2, 1), (3, 2), (4, 3), (6, 5), (8, 7), (12, 10), (16, 14), (32, 14)):
            with self.subTest(physical=physical):
                self.assertEqual(p.core_budget(physical), expected)
        for invalid in (0, -1, True, 8.5):
            with self.assertRaises(ValueError):
                p.core_budget(invalid)
        for invalid in (0, 8, 14, True, 3.5):
            with self.assertRaises(ValueError):
                p.core_budget(8, invalid)
        self.assertEqual(p.core_budget(16, 1), 1)

    def test_all_allocations_bounded(self):
        for budget in range(1, 15):
            for a in range(1, 15):
                for b in range(1, 15):
                    out = p.grants([a, b], budget)
                    self.assertLessEqual(sum(out), budget)
                    self.assertTrue(0 <= out[0] <= a and 0 <= out[1] <= b)
        self.assertEqual(p.grants([10, 10], 14), [10, 0])
        self.assertEqual(p.grants([10, 2], 14), [10, 2])
        self.assertEqual(p.grants([10, 10], 6), [0, 0])

    def test_workload_classes(self):
        for n, shell, expected in ((30, False, 4), (100, False, 4), (250, False, 4),
                                   (502, True, 7), (526, False, 7), (2000, True, 7)):
            self.assertEqual(p.preferred_cores(n, shell), expected)
        with self.assertRaises(ValueError):
            p.preferred_cores(0, False)

    def test_basis_counts_shells_not_primitives(self):
        b = 'H 0\nS 3 1.0\n3.0 0.1\n2.0 0.2\n1.0 0.3\nP 1 1.0\n1.0 1.0\nD 1 1.0\n1.0 1.0\n****\n'
        self.assertEqual(p.basis_count(b, [1, 1], True), 18)
        self.assertEqual(p.basis_count(b, [1], False), 10)
        b += 'Re 0\nRe-ECP 4 60\nignored ECP potentials\n'
        self.assertEqual(p.basis_count(b, [1], True), 9)
        b = 'C 0\nSP 1 1.0\n1.0 0.5 0.4\n****\n'
        self.assertEqual(p.basis_count(b, [6, 6], True), 8)

    def test_bad_basis_rejected(self):
        for b in ('', 'H 0\nS 2 1.0\n1.0 1.0\n****', 'H 0\nS 1 1.0\nnan 1.0\n****',
                  'H 0\nS 1 1.0\n1 1\n****\nH 0\nS 1 1.0\n1 1\n****'):
            with self.assertRaises(ValueError):
                p.basis_count(b, [1], True)

    def test_split_and_resource_independence(self):
        a, b = p.initial_part(INPUT, 'opt', 1), p.initial_part(INPUT, 'opt', 14)
        self.assertNotIn(' Freq', a)
        for keyword in ('VeryTight', 'CalcFC', 'MaxCycles=512', 'SuperFineGrid', 'NoSymm', 'GD3BJ', '5D 7F'):
            self.assertIn(keyword, a)
        self.assertEqual(p.route_signature(a), p.route_signature(b))
        self.assertIn('%nprocshared=14', b)
        self.assertIn('MaxCycle=2048', p.initial_part(INPUT, 'opt', 4, retry=True))
        with self.assertRaises(ValueError):
            p.initial_part(INPUT.replace('VeryTight', 'Loose'), 'opt', 4)

    def test_restart_does_not_fake_freq_restart(self):
        opt = p.followup(REFERENCE, 'UPBE1PBE', 'opt', 6, restart_opt=True)
        self.assertIn('Opt=(Restart,VeryTight', opt)
        self.assertNotIn('Freq', opt)
        freq = p.followup(REFERENCE, 'UPBE1PBE', 'freq', 6)
        self.assertIn('UPBE1PBE/ChkBasis Freq Geom=AllCheck Guess=Read', freq)
        self.assertNotIn('Restart', freq)
        with self.assertRaises(ValueError):
            p.followup(REFERENCE, 'UPBE1PBE', 'freq', 6, restart_opt=True)
        with self.assertRaises(ValueError):
            p.followup(REFERENCE, 'B3LYP', 'opt', 6)

    def test_chain_has_separate_frequency_boundary(self):
        rows = {}
        self.assertEqual(p.next_stage(REFERENCE, rows).name, 'opt-00')
        rows['opt-00'] = record()
        self.assertEqual(p.next_stage(REFERENCE, rows).name, 'freq-00')
        rows['freq-00'] = record()
        self.assertEqual(p.next_stage(REFERENCE, rows).name, 'stability-00')
        rows['stability-00'] = record()
        self.assertEqual(p.next_stage(REFERENCE, rows), 'candidate_collected')

    def test_repair_two_cycles_not_three(self):
        rows = {}
        for i in range(3):
            rows[f'opt-{i:02d}'] = record()
            rows[f'freq-{i:02d}'] = record()
            rows[f'stability-{i:02d}'] = record(unstable=True)
            if i < 2:
                self.assertEqual(p.next_stage(REFERENCE, rows).name, f'opt-{i+1:02d}')
        self.assertEqual(p.next_stage(REFERENCE, rows), 'stability_repair_cycles_exhausted')

    def test_r_to_u_and_energy_change_require_relaxation(self):
        rows = {'opt-00': record(method='RPBE1PBE'), 'freq-00': record(method='RPBE1PBE'),
                'stability-00': record()}
        self.assertEqual(p.next_stage(REFERENCE, rows).name, 'opt-01')
        rows['freq-00'] = record()
        rows['stability-00'] = record(energy=-1.1)
        self.assertEqual(p.next_stage(REFERENCE, rows).name, 'opt-01')

    def test_atoms_and_scan_frames_never_optimized(self):
        for geometry, part in (('single_atom', 'sp'), ('fixed', 'force')):
            reference = copy.deepcopy(REFERENCE)
            reference['purpose']['reference_geometry'] = geometry
            first = p.next_stage(reference, {})
            self.assertEqual(first.part, part)
            self.assertEqual(p.next_stage(reference, {first.name: record()}).part, 'stability')

    def test_progress_no_percent_or_fabricated_eta(self):
        log = ('Step number 37 out of a maximum of 512\nNBasis= 502\nCycle 11 Pass 1\n'
               'Maximum Force 0.000255 0.000002 NO\nRMS     Force 0.000071 0.000001 NO\n'
               'Maximum Displacement 0.212566 0.000006 NO\nRMS     Displacement 0.042942 0.000004 NO\n')
        result = p.log_progress(log)
        self.assertEqual(result['nbasis'], 502)
        self.assertEqual(result['optimization_step'], 37)
        self.assertEqual(len(result['last_convergence_tables'][0]), 4)
        self.assertEqual(result['last_convergence_tables'][0][0]['ratio_to_printed_threshold'], 127.50000000000001)
        self.assertIsNone(result['eta_seconds'])
        self.assertNotIn('percent', result)
        self.assertEqual(p.log_progress(log + 'Optimization completed.\n')['phase_observed'], 'post_optimization')

    def test_atomic_writes_and_lock_collision(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r.atomic(root / 'state.json', {'a': 1})
            r.atomic(root / 'state.json', {'a': 2})
            self.assertEqual(r.read(root / 'state.json'), {'a': 2})
            self.assertFalse(list(root.glob('*.tmp')))
            with r.lease(root / 'supervisor.lock'):
                self.assertTrue(r.occupied(root / 'supervisor.lock'))
            self.assertFalse(r.occupied(root / 'supervisor.lock'))


class FakeData:
    def __init__(self, root, count=2):
        self.parent_identity = 'legacy'
        self.binaries = {'g16.exe': 'fake', 'formchk.exe': 'fake'}
        self.manifest = {'reference_set_identity': 'reference'}
        self.gaussian = root / 'Gaussian with spaces'
        self.references = root / 'work/references/REF-001'
        self.candidates = {}
        self.originals = {}
        self.sizes = {}
        self.runner = SimpleNamespace(stage_diagnostics=diag)
        for i in range(1, count+1):
            case = f'OLD-{i:03d}'
            self.candidates[case] = dict(REFERENCE, case_id=case)
            self.originals[case] = {'fchk_identity': {'has_beta_coefficients': True}}
            self.sizes[case] = 502
            folder = self.references / 'reference-inputs' / case
            folder.mkdir(parents=True)
            (folder / 'initial.gjf').write_text(INPUT)

    def validate(self, case_id, path):
        fields = r.read(path)
        if fields.get('Number of basis functions') != 502:
            raise r.NeedsReview('Wrong basis')
        return fields


class FakeBackend:
    def __init__(self):
        self.calls = []
        self.interrupt_next = False
        self.fail_next = False
        self.hold_case = None
        self.release = threading.Event()
        self.started = threading.Event()
        self.cancel_case = None
        self.completion_hook = lambda cwd: None

    def host_info(self):
        return {'available_gib': 128, 'total_gib': 192, 'physical_cores': 16}

    def run_tree(self, argv, cwd, cores, memory, timeout, cancel=lambda: False,
                 on_started=lambda x: None, on_tick=lambda x: None):
        argv = [Path(a) for a in argv]
        if argv[0].name == 'formchk.exe':
            r.read(argv[1])  # reject corrupt mock binary
            shutil_copy(argv[1], argv[2])
            return {'exit_code': 0, 'interrupted': False}
        self.calls.append({'cwd': cwd, 'cores': cores, 'input': r.text(argv[1])})
        on_started({'pid': 12345, 'cores': cores})
        on_tick({'active_processes': 1, 'elapsed_seconds': 0.01})
        self.started.set()
        if self.hold_case and self.hold_case in str(cwd):
            deadline = time.monotonic() + 10
            while not self.release.is_set() and not cancel() and time.monotonic() < deadline:
                time.sleep(.01)
        r.atomic(cwd / 'job.chk', dict(FIELDS))
        if self.interrupt_next or cancel():
            self.interrupt_next = False
            (cwd / 'job.log').write_text('Berny optimization\nStep number 12 out of a maximum of 512\n')
            return {'exit_code': 123, 'interrupted': True, 'stop_reason': 'operator_interrupt'}
        if self.fail_next:
            self.fail_next = False
            (cwd / 'job.log').write_text('Error termination: synthetic test failure\n')
            return {'exit_code': 1, 'interrupted': False}
        contents = r.text(argv[1])
        log = ('Optimization completed.\n' if 'Opt=' in contents else '')
        if ' Freq ' in contents:
            log += 'Frequencies -- 100.0\n'
        if 'Stable=Opt' in contents:
            log += 'The wavefunction is stable under the perturbations considered.\n'
        log += 'Normal termination of Gaussian 16\n'
        (cwd / 'job.log').write_text(log)
        self.completion_hook(cwd)
        return {'exit_code': 0, 'interrupted': False, 'wall_seconds': .01, 'tree_cpu_seconds': .01}


def shutil_copy(a, b):
    import shutil
    shutil.copy2(a, b)


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.data = FakeData(self.root)
        self.backend = FakeBackend()
        self.config = {'work_root': str(self.root / 'work'), 'timeout_hours': 1, 'disk_reserve_gib': 1}
        self.engine = r.Engine(self.config, self.data, self.backend)
        self.engine.poll_seconds = .01
        self.engine.bind()
        self.engine.set_mode('run')

    def tearDown(self):
        self.tmp.cleanup()

    def phase(self, case='OLD-001'):
        self.engine._records_cache.pop(case, None)
        return p.next_stage(self.data.candidates[case], self.engine.records(case))

    def test_success_and_stage_reuse(self):
        for _ in range(3):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertEqual(self.phase(), 'candidate_collected')
        self.assertEqual(len(self.backend.calls), 3)
        self.assertTrue(all(not x.get('scientific_pass') for x in self.engine.records('OLD-001').values()))

    def test_pause_does_not_start_gaussian(self):
        self.engine.set_mode('pause')
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertFalse(self.backend.calls)

    def test_interrupted_opt_changes_resources_without_losing_geometry(self):
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 4)
        self.engine.calculate('OLD-001', self.phase(), 14)
        self.assertIn('Opt=(Restart,VeryTight', self.backend.calls[-1]['input'])
        self.assertIn('%nprocshared=14', self.backend.calls[-1]['input'])
        attempts = sorted((self.engine.jobs / 'OLD-001/opt-00').glob('attempt-*/attempt.json'))
        self.assertEqual(len(attempts), 2)
        self.assertEqual(r.read(attempts[0])['status'], 'interrupted')
        self.assertEqual(r.read(attempts[1])['status'], 'collected')

    def test_analytic_frequency_replay_preserves_completed_opt(self):
        self.engine.calculate('OLD-001', self.phase(), 7)
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.engine.calculate('OLD-001', self.phase(), 14)
        self.assertIn(' Freq ', self.backend.calls[-1]['input'])
        self.assertNotIn('Opt=', self.backend.calls[-1]['input'])
        self.assertNotIn('Freq=Restart', self.backend.calls[-1]['input'])
        self.assertEqual(len([a for a in self.backend.calls if 'Opt=' in a['input']]), 1)

    def test_corrupt_interrupted_checkpoint_not_silently_restarted(self):
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 7)
        chk = self.backend.calls[0]['cwd'] / 'job.chk'
        chk.write_text('corrupted')
        with self.assertRaises(Exception):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertEqual(len(self.backend.calls), 1)
        self.assertEqual(chk.read_text(), 'corrupted')

    def test_interrupted_input_tampering_rejected(self):
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 7)
        (self.backend.calls[0]['cwd'] / 'job.gjf').write_text('modified method')
        with self.assertRaises(r.NeedsReview):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertEqual(len(self.backend.calls), 1)

    def test_missing_interrupted_checkpoint_needs_review(self):
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.engine.calculate('OLD-001', self.phase(), 7)
        (self.backend.calls[0]['cwd'] / 'job.chk').unlink()
        with self.assertRaises(r.NeedsReview):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertEqual(len(self.backend.calls), 1)

    def test_completed_artifact_tampering_rejected(self):
        self.engine.calculate('OLD-001', self.phase(), 7)
        rows = self.engine.records('OLD-001')
        Path(rows['opt-00']['checkpoint']).write_text('tampered')
        self.engine._records_cache.clear()
        with self.assertRaises(ValueError):
            self.engine.records('OLD-001')

    def test_runtime_identity_change_rejected(self):
        value = r.read(self.engine.root / 'binding.json')
        value['runtime_identity'] = 'bad'
        r.atomic(self.engine.root / 'binding.json', value)
        with self.assertRaises(ValueError):
            self.engine.bind()

    def test_failure_not_infinite_retry(self):
        self.backend.fail_next = True
        with self.assertRaises(r.NeedsReview):
            self.engine.calculate('OLD-001', self.phase(), 7)
        self.assertEqual(len(self.backend.calls), 1)

    def test_schedule_small_host_respects_seven_core_budget(self):
        with patch('resume.shutil.disk_usage', return_value=SimpleNamespace(free=100*2**30)):
            self.engine.run(['OLD-001', 'OLD-002'], {'physical_cores': 8, 'total_gib': 96})
        self.assertEqual(r.read(self.engine.root / 'status.json')['physical_core_budget'], 7)
        self.assertTrue(all(x['cores'] == 4 for x in self.backend.calls))
        self.assertEqual(r.read(self.engine.jobs / 'OLD-001/result.json')['status'], 'candidate_collected')

    def test_waiting_case_precedes_completed_cases_next_phase(self):
        third = 'OLD-003'
        self.data.candidates[third] = dict(REFERENCE, case_id=third)
        self.data.originals[third] = self.data.originals['OLD-001']
        self.data.sizes[third] = 502
        folder = self.data.references / 'reference-inputs' / third
        folder.mkdir(); (folder / 'initial.gjf').write_text(INPUT)
        self.engine.config['workers'] = 1
        with patch('resume.shutil.disk_usage', return_value=SimpleNamespace(free=100*2**30)):
            self.engine.run(['OLD-001','OLD-002',third], {'physical_cores':16,'total_gib':192})
        starts = [(next(x for x in c['cwd'].parts if x.startswith('OLD-')), c['cwd'].parent.name) for c in self.backend.calls]
        self.assertEqual(starts[:3], [('OLD-001','opt-00'),('OLD-002','opt-00'),('OLD-003','opt-00')])
        self.assertTrue(all(c['cores']==7 for c in self.backend.calls))

    def test_cancel_pause_refills_idle_worker_while_other_job_still_runs(self):
        self.backend.hold_case = 'OLD-001'
        paused = threading.Event()
        observed, errors = [], []
        def hook(cwd):
            if 'OLD-002' in str(cwd) and 'opt-00' in str(cwd):
                self.engine.set_mode('pause')
                paused.set()
        self.backend.completion_hook = hook
        def controller():
            try:
                self.assertTrue(paused.wait(3))
                time.sleep(.15)
                self.assertEqual(len(self.backend.calls), 2)
                self.engine.set_mode('run')
                deadline = time.monotonic() + 3
                while time.monotonic() < deadline:
                    if any('OLD-002' in str(c['cwd']) and 'freq-00' in str(c['cwd']) for c in self.backend.calls):
                        self.assertFalse(self.backend.release.is_set())
                        observed.append(True)
                        break
                    time.sleep(.01)
            except BaseException as error:
                errors.append(error)
            finally:
                self.backend.release.set()
        thread = threading.Thread(target=controller)
        thread.start()
        with patch('resume.shutil.disk_usage', return_value=SimpleNamespace(free=100*2**30)):
            self.engine.run(['OLD-001', 'OLD-002'], {'physical_cores': 16, 'total_gib': 192})
        thread.join(5)
        self.assertFalse(errors)
        self.assertTrue(observed)
        self.assertEqual(len(self.backend.calls), 6)

    def test_legacy_completed_import_is_nondestructive_and_idempotent(self):
        folder = self.root / 'work/jobs/REF-001/OLD-001/initial/attempt-01'
        folder.mkdir(parents=True)
        r.atomic(folder / 'job.chk', dict(FIELDS))
        r.atomic(folder / 'job.fch', dict(FIELDS))
        (folder / 'job.log').write_text('Optimization completed.\nFrequencies -- 100.0\nNormal termination of Gaussian 16\n')
        original = {'status': 'collected', 'runner_identity': 'legacy', 'reference_identity': 'science',
                    'checkpoint': str(folder / 'job.chk'), 'fchk': str(folder / 'job.fch'),
                    'fchk_sha256': r.sha(folder / 'job.fch'), 'log': str(folder / 'job.log')}
        r.atomic(folder.parent / 'stage.json', original)
        before = {f: r.sha(f) for f in folder.parent.rglob('*') if f.is_file()}
        summary = self.engine.import_legacy()
        self.assertEqual(summary['imported_parts'], 2)
        self.assertEqual(self.engine.import_legacy()['imported_parts'], 0)
        self.assertEqual(before, {f: r.sha(f) for f in before})
        self.assertEqual(self.phase().name, 'stability-00')


class RepositoryFixtureTests(unittest.TestCase):
    def test_frozen_reference_digest_numeric_key_order(self):
        ref = dict(REFERENCE, ecp_core_electrons_by_atomic_number={1: 0, 7: 0, 24: 0})
        ref.pop('identity')
        old = p.digest(ref)
        decoded = json.loads(json.dumps(ref))
        self.assertNotEqual(p.digest(decoded), old)
        self.assertEqual(p.historical_reference_digest(decoded), old)

    def test_all_273_real_inputs_when_repository_is_available(self):
        root = r.FROZEN.parent
        manifest = root / 'reference-candidates.json'
        if not manifest.exists():
            self.skipTest('Full repository fixture not in this local patch workspace; CI must run this test')
        refs = r.read(manifest)['candidates']
        campaign = r.read(root / 'campaign.json')
        cases = {c['case_id']: c for c in campaign['cases']}
        self.assertEqual(len(refs), 273)
        sizes = {}
        for ref in refs:
            case = ref['case_id']
            folder = root / 'reference-inputs' / case
            with self.subTest(case=case):
                self.assertEqual(p.historical_reference_digest(ref), ref['identity'])
                self.assertEqual(r.sha(folder / 'initial.gjf'), ref['initial_input_sha256'])
                self.assertEqual(r.sha(folder / 'basis.gbs'), ref['basis_sha256'])
                n = p.basis_count(r.text(folder / 'basis.gbs'),
                    cases[case]['fchk_identity']['Atomic numbers'], ref['shell_flags'] == '5D 7F')
                self.assertGreater(n, 0)
                sizes[case] = n
        self.assertEqual(sizes['OLD-018'], 502)
        print('Real frozen reference basis sizes:', json.dumps(sizes, sort_keys=True))


@unittest.skipUnless(os.name == 'nt', 'Native Windows Job Object tests; not Gaussian tests')
class NativeWindowsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(r.FROZEN))
        import windows_job
        cls.backend = windows_job
        cls.backend.full_startup_affinity()

    def test_topology_not_smt_count(self):
        host = self.backend.host_info()
        self.assertLessEqual(host['physical_cores'], host['logical_processors'])
        self.assertLessEqual(p.core_budget(host['physical_cores']), 14)

    def test_normal_process_with_unicode_and_space_paths(self):
        with tempfile.TemporaryDirectory(prefix='COV 空格 ') as name:
            folder = Path(name)
            (folder / 'scratch').mkdir()
            result = self.backend.run_tree([sys.executable, '-c', 'import os; print(os.getenv("OMP_NUM_THREADS"))'],
                                           folder, 1, 1, 30)
            self.assertEqual(result['exit_code'], 0)
            self.assertEqual(result['active_processes_at_return'], 0)
            self.assertIn('1', (folder / 'launcher.log').read_text())

    def test_interrupt_waits_for_entire_owned_tree(self):
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name)
            (folder / 'scratch').mkdir()
            code = ('import subprocess,sys,time; '
                    'subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"]); time.sleep(60)')
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', code], folder, 1, 1, 30,
                                           cancel=lambda: time.monotonic()-start > 2)
            self.assertTrue(result['interrupted'])
            self.assertEqual(result['active_processes_at_return'], 0)
            self.assertGreaterEqual(result['tree_process_count'], 2)
            self.assertLess(result['wall_seconds'], 15)


if __name__ == '__main__':
    unittest.main(verbosity=2)
