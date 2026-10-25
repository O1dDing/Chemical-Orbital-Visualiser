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

"""Regression tests for quick RAM pause, cold save and gated native restarts.

NativeWindowsPauseTests runs Python processes, NOT licensed Gaussian.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fast_checkpoint as cp
import resume as r
from policy import Stage, next_stage, route_signature
from test_runtime import FakeBackend, FakeData, INPUT, FIELDS


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name)
        for name, data in (('job.chk', b'chk'), ('job.rwf', b'rwf'), ('job.gjf', b'input'), ('job.log', b'log')):
            (self.path / name).write_bytes(data)
        self.producer = {'runtime': 'r', 'case': 'case', 'binaries': 'b'}

    def tearDown(self):
        self.tmp.cleanup()

    def test_ram_pause_cannot_be_used_as_disk_snapshot(self):
        with self.assertRaises(ValueError):
            cp.cold_snapshot(self.path, self.producer, writers_exited=False)
        self.assertFalse(list(self.path.glob('snapshot-*')))

    def test_complete_snapshot_identity_and_independent_restore(self):
        snap = cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        dst = self.path / 'new'; dst.mkdir()
        cp.restore_snapshot(snap, dst, self.producer, require_rwf=True)
        self.assertEqual((dst / 'job.chk').read_bytes(), b'chk')
        (dst / 'job.rwf').write_bytes(b'new iteration')
        self.assertEqual((Path(snap['directory']) / 'job.rwf').read_bytes(), b'rwf')
        self.assertFalse((dst / 'job.log').exists())
        self.assertFalse((dst / 'job.gjf').exists())
        self.assertFalse(r.read(Path(snap['directory']) / 'manifest.json')['restart_validated'])

    def test_tampered_snapshot_member_rejected(self):
        snap = cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        (Path(snap['directory']) / 'job.rwf').write_bytes(b'bad')
        dst = self.path / 'new'; dst.mkdir()
        with self.assertRaises(ValueError):
            cp.restore_snapshot(snap, dst, self.producer, require_rwf=True)
        self.assertFalse((dst / 'job.chk').exists())

    def test_other_science_or_binary_rejected(self):
        snap = cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        dst = self.path / 'new'; dst.mkdir()
        with self.assertRaises(ValueError):
            cp.restore_snapshot(snap, dst, {'runtime': 'new'}, require_rwf=True)

    def test_opt_copy_small_and_raw_scratch_retained(self):
        snap = cp.cold_snapshot(self.path, self.producer, writers_exited=True, include_scratch=False)
        self.assertTrue((self.path / 'job.rwf').is_file())
        self.assertNotIn('job.rwf', r.read(Path(snap['directory']) / 'manifest.json')['files'])
        dst = self.path / 'new'; dst.mkdir()
        with self.assertRaises(ValueError):
            cp.restore_snapshot(snap, dst, self.producer, require_rwf=True)

    def test_missing_space_never_deletes_original(self):
        before = {p.name: p.read_bytes() for p in self.path.iterdir()}
        with patch('fast_checkpoint.shutil.disk_usage', return_value=SimpleNamespace(free=1)):
            with self.assertRaises(OSError):
                cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.path.iterdir()})

    def test_hot_file_mutation_during_copy_rejected(self):
        original = cp.shutil.copyfileobj
        def mutate(src, dst, *args):
            original(src, dst, *args)
            with (self.path / Path(src.name).name).open('ab') as file:
                file.write(b'live writer')
        with patch('fast_checkpoint.shutil.copyfileobj', side_effect=mutate):
            with self.assertRaises(ValueError):
                cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        self.assertTrue(all(p.name.endswith('.partial') for p in self.path.glob('snapshot-*')))

    def test_path_injection_manifest_rejected(self):
        snap = cp.cold_snapshot(self.path, self.producer, writers_exited=True)
        manifest = Path(snap['directory']) / 'manifest.json'
        data = r.read(manifest); data['files']['../outside'] = {'bytes': 0, 'sha256': ''}
        r.atomic(manifest, data); snap['manifest_sha256'] = r.sha(manifest)
        dst = self.path / 'new'; dst.mkdir()
        with self.assertRaises(ValueError):
            cp.restore_snapshot(snap, dst, self.producer)

    def test_named_rwf_does_not_change_scientific_route(self):
        augmented = cp.persistent_input(INPUT)
        self.assertIn('%RWF=job.rwf', augmented)
        self.assertIn('%Save', augmented)
        self.assertNotIn('%NoSave', augmented)
        self.assertEqual(route_signature(INPUT), route_signature(augmented))
        with self.assertRaises(ValueError):
            cp.persistent_input('%NoSave\n' + INPUT)

    def test_analytic_restart_is_general_restart_not_freq_restart(self):
        result = cp.rwf_restart_input(14)
        self.assertIn('#p Restart\n\n', result)
        self.assertNotIn('Freq=Restart', result)
        self.assertNotIn('Opt=', result)

    def test_controlled_opt_stop_not_a_generic_failure(self):
        log = 'Step number 2\nMaximum Force\nLeave Link 103 at time\nLeave Link 103 at time\n'
        input_text = cp.persistent_input(INPUT, opt_segments=True)
        self.assertTrue(cp.cooperative_opt_stop(log, input_text))
        for extra in ('Error termination', 'Convergence failure', 'Erroneous write', 'Leave Link 502 at time'):
            self.assertFalse(cp.cooperative_opt_stop(log + extra, input_text))
        self.assertFalse(cp.cooperative_opt_stop(log, INPUT))

    def test_held_duration_not_charged_to_active_timeout(self):
        clock = cp.ActiveClock(0)
        clock.set_held(True, 10)
        self.assertEqual(clock.active(10000), 10)
        clock.set_held(False, 10000)
        self.assertEqual(clock.active(10010), 20)
        self.assertEqual(clock.paused(10010), 9990)
        clock.set_held(False, 20000)
        self.assertEqual(clock.paused(20000), 9990)


class FastEngineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.data = FakeData(self.root)
        self.backend = FakeBackend()
        self.config = {'work_root': str(self.root / 'work'), 'timeout_hours': 1, 'disk_reserve_gib': 1}
        self.e = r.Engine(self.config, self.data, self.backend)
        self.e.bind(); self.e.set_mode('run')

    def tearDown(self):
        self.tmp.cleanup()

    def phase(self):
        self.e._records_cache.clear()
        return next_stage(self.data.candidates['OLD-001'], self.e.records('OLD-001'))

    def test_hold_blocks_new_dispatch_without_killing_anything(self):
        self.e.set_mode('hold')
        with self.assertRaises(r.Paused):
            self.e.calculate('OLD-001', self.phase(), 4)
        self.assertFalse(self.backend.calls)

    def test_resume_does_not_erase_unacknowledged_stop(self):
        before = time.time() - 1
        self.e.set_mode('interrupt')
        self.e.set_mode('run')
        self.assertTrue(self.e.interrupt_requested(before))
        self.assertFalse(self.e.interrupt_requested(time.time()+1))

    def test_no_opt_segments_without_native_receipt(self):
        self.e.config['opt_step_checkpoints'] = True
        with self.assertRaises(r.NeedsReview):
            self.e.calculate('OLD-001', self.phase(), 4)
        self.assertFalse(self.backend.calls)

    def test_native_receipt_must_have_real_evidence(self):
        r.atomic(self.e.root / 'native-capabilities.json', {'runtime_identity': self.e.identity,
                 'binaries': self.data.binaries, 'capabilities': {'opt_l103_segments':
                 {'passed': True, 'methods': ['UPBE1PBE']}}, 'evidence': []})
        with self.assertRaises(r.NeedsReview):
            self.e.require_capability('opt_l103_segments', 'UPBE1PBE')

    def test_interrupted_freq_rwf_gated_no_silent_full_replay(self):
        self.e.calculate('OLD-001', self.phase(), 4)
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.e.calculate('OLD-001', self.phase(), 4)
        count = len(self.backend.calls)
        self.e.config['freq_recovery'] = 'rwf'
        with self.assertRaises(r.NeedsReview):
            self.e.calculate('OLD-001', self.phase(), 4)
        self.assertEqual(len(self.backend.calls), count)

    def test_cold_save_receipt_not_completed_candidate(self):
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.e.calculate('OLD-001', self.phase(), 4)
        path = sorted((self.e.jobs / 'OLD-001/opt-00').glob('attempt-*/attempt.json'))[-1]
        row = r.read(path)
        self.assertEqual(row['status'], 'interrupted')
        self.assertIn('checkpoint_snapshot', row)
        self.assertFalse((path.parent.parent / 'stage.json').exists())
        self.assertIsNotNone(self.e.interrupted_source('OLD-001', self.phase())['snapshot'])

    def accept_mock_capability(self, capability):
        proof = self.root / 'mock-proof.txt'; proof.write_text('Mock evidence, not native acceptance')
        r.atomic(self.e.root / 'native-capabilities.json', {'runtime_identity': self.e.identity,
            'binaries': self.data.binaries, 'capabilities': {capability:
            {'passed': True, 'methods': ['UPBE1PBE']}},
            'evidence': [{'path': str(proof), 'sha256': r.sha(proof)}]})

    def test_gated_rwf_restart_uses_saved_scratch_not_new_freq(self):
        self.e.calculate('OLD-001', self.phase(), 4)
        original_run = self.backend.run_tree
        def run(argv, cwd, *args, **kw):
            if Path(argv[0]).name == 'g16.exe':
                inp = r.text(Path(argv[1]))
                if '#p Restart' in inp:
                    self.assertEqual((cwd / 'job.rwf').read_text(), 'saved derivative state')
                else:
                    (cwd / 'job.rwf').write_text('saved derivative state')
                result = original_run(argv, cwd, *args, **kw)
                if '#p Restart' in inp:
                    (cwd / 'job.log').write_text('Frequencies -- 100.0\nNormal termination of Gaussian 16\n')
                return result
            return original_run(argv, cwd, *args, **kw)
        self.backend.run_tree = run
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.e.calculate('OLD-001', self.phase(), 4)
        self.accept_mock_capability('analytic_rwf_restart')
        self.e.config['freq_recovery'] = 'rwf'
        self.e.calculate('OLD-001', self.phase(), 6)
        self.assertIn('#p Restart', self.backend.calls[-1]['input'])
        self.assertNotIn(' Freq ', self.backend.calls[-1]['input'])
        self.assertEqual(self.phase().part, 'stability')
        self.assertEqual(len(self.backend.calls), 3)

    def test_checkpointed_opt_segment_is_not_candidate_and_resumes(self):
        self.accept_mock_capability('opt_l103_segments')
        self.e.config['opt_step_checkpoints'] = True
        original_run = self.backend.run_tree
        segmented_calls = []
        def run(argv, cwd, *args, **kw):
            result = original_run(argv, cwd, *args, **kw)
            if Path(argv[0]).name == 'g16.exe':
                segmented_calls.append(cwd)
                if len(segmented_calls) == 1:
                    (cwd / 'job.log').write_text('Step number 2\nMaximum Force\n'
                        'Leave Link 103 at time\nLeave Link 103 at time\n')
                    result['exit_code'] = 1
            return result
        self.backend.run_tree = run
        self.assertEqual(self.e.calculate('OLD-001', self.phase(), 4), 'checkpointed')
        self.assertEqual(self.phase().part, 'opt')
        self.assertFalse((self.e.jobs / 'OLD-001/result.json').exists())
        self.e.calculate('OLD-001', self.phase(), 6)
        self.assertIn('Opt=(Restart,VeryTight', self.backend.calls[-1]['input'])
        self.assertEqual(self.phase().part, 'freq')

    def test_missing_interrupted_chk_does_not_silently_replay_freq_in_auto(self):
        self.e.calculate('OLD-001', self.phase(), 4)
        self.backend.interrupt_next = True
        with self.assertRaises(r.Paused):
            self.e.calculate('OLD-001', self.phase(), 4)
        attempt = self.backend.calls[-1]['cwd']
        (attempt / 'job.chk').unlink()
        self.accept_mock_capability('analytic_rwf_restart')
        self.e.config['freq_recovery'] = 'auto'
        with self.assertRaises(ValueError):
            self.e.calculate('OLD-001', self.phase(), 4)
        self.assertEqual(len(self.backend.calls), 2)

    def test_old_v2_work_directory_not_overwritten(self):
        old = self.root / 'work/runtime-v2/binding.json'; old.parent.mkdir(parents=True)
        old.write_text('old immutable record')
        self.e.bind()
        self.assertEqual(old.read_text(), 'old immutable record')
        self.assertNotEqual(old.parent, self.e.root)


@unittest.skipUnless(os.name == 'nt', 'Real Windows process pause/resume tests (not Gaussian)')
class NativeWindowsPauseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(r.FROZEN))
        import windows_job
        cls.backend = windows_job
        cls.backend.full_startup_affinity()

    def test_live_parent_child_pause_resume_and_foreign_process_untouched(self):
        with tempfile.TemporaryDirectory(prefix='COV pause ') as tmp:
            path = Path(tmp); (path / 'scratch').mkdir()
            worker = path / 'worker.py'
            worker.write_text('import sys,time,subprocess\nfrom pathlib import Path\n'
                'p=Path(sys.argv[1])\n'
                'if len(sys.argv)>2: subprocess.Popen([sys.executable,__file__,str(p.with_name("child.txt"))])\n'
                'for i in range(35):\n p.write_text(str(i)); time.sleep(.07)\n')
            foreign = subprocess.Popen([sys.executable, str(worker), str(path / 'foreign.txt')])
            want = threading.Event(); held = threading.Event(); release = threading.Event()
            samples = []; errors = []
            def tick(value):
                samples.append(value)
                if value.get('execution_state') == 'ram_paused':
                    held.set()
            def controller():
                try:
                    deadline = time.monotonic()+8
                    while not (path / 'child.txt').exists() and time.monotonic() < deadline:
                        time.sleep(.02)
                    want.set()
                    self.assertTrue(held.wait(6))
                    time.sleep(.2)
                    a = (path/'owned.txt').read_text(); b = (path/'child.txt').read_text()
                    f = (path/'foreign.txt').read_text()
                    time.sleep(.7)
                    self.assertEqual(a, (path/'owned.txt').read_text())
                    self.assertEqual(b, (path/'child.txt').read_text())
                    self.assertNotEqual(f, (path/'foreign.txt').read_text())
                except BaseException as error:
                    errors.append(error)
                finally:
                    want.clear(); release.set()
            thread = threading.Thread(target=controller); thread.start()
            started = time.monotonic()
            try:
                result = self.backend.run_tree([sys.executable, worker, path/'owned.txt', 'spawn'],
                    path, 1, 1, 12, hold=want.is_set, on_tick=tick,
                    cancel=lambda: time.monotonic()-started > 15)
            finally:
                want.clear(); thread.join(8)
                foreign.wait(timeout=10)
            self.assertFalse(errors, repr(errors) + ' process=' + json.dumps(result))
            self.assertEqual(result['exit_code'], 0)
            self.assertEqual(result['active_processes_at_return'], 0)
            self.assertGreater(result['ram_paused_seconds'], .5)
            self.assertGreaterEqual(result['tree_process_count'], 2)
            self.assertEqual((path/'owned.txt').read_text(), '34')
            print('Windows RAM-pause metrics:', json.dumps(result))

    def test_pause_longer_than_timeout_does_not_kill_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path/'scratch').mkdir()
            want = threading.Event(); held = threading.Event()
            def controller():
                time.sleep(.15); want.set()
                held.wait(5); time.sleep(3); want.clear()
            thread = threading.Thread(target=controller); thread.start()
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', 'import time\nfor i in range(10): time.sleep(.08)'],
                path, 1, 1, 2.5, hold=want.is_set,
                on_tick=lambda v: held.set() if v.get('execution_state') == 'ram_paused' else None,
                cancel=lambda: time.monotonic()-start > 10)
            thread.join(6)
            self.assertTrue(held.is_set())
            self.assertFalse(result['interrupted'])
            self.assertGreater(result['wall_seconds'], 3)
            self.assertLess(result['active_wall_seconds'], 2.5)

    def test_interrupt_from_held_state_leaves_no_descendants(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path/'scratch').mkdir()
            held = threading.Event()
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', 'import time; time.sleep(30)'],
                path, 1, 1, 15, hold=lambda: time.monotonic()-start > .2,
                on_tick=lambda v: held.set() if v.get('execution_state') == 'ram_paused' else None,
                cancel=lambda: held.is_set() or time.monotonic()-start > 8)
            self.assertTrue(held.is_set())
            self.assertTrue(result['interrupted'])
            self.assertEqual(result['stop_reason'], 'operator_interrupt')
            self.assertEqual(result['active_processes_at_return'], 0)

    def test_active_timeout_can_park_instead_of_terminate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path/'scratch').mkdir()
            want = threading.Event(); held = threading.Event()
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', 'import time;time.sleep(20)'],
                path, 1, 1, .25, hold=want.is_set, on_timeout=want.set,
                on_tick=lambda v: held.set() if v.get('execution_state') == 'ram_paused' else None,
                cancel=lambda: held.is_set() or time.monotonic()-start > 8)
            self.assertTrue(want.is_set())
            self.assertTrue(held.is_set())
            self.assertEqual(result['stop_reason'], 'operator_interrupt')

    def check_early_deadline(self, timeout):
        with self.subTest(timeout=timeout), tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path / 'scratch').mkdir()
            want = threading.Event(); held = threading.Event()
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', 'import time; time.sleep(20)'],
                path, 1, 1, timeout, hold=want.is_set, on_timeout=want.set,
                on_tick=lambda v: held.set() if v['execution_state'] == 'ram_paused' else None,
                cancel=lambda: held.is_set() or time.monotonic()-start > 8)
            self.assertTrue(held.is_set(), result)
            self.assertEqual(result['pause_failures'], [])
            self.assertEqual(result['active_processes_at_return'], 0)
            self.assertLess(result['wall_seconds'], 8)
            print('Early pause trial:', json.dumps(result))

    def test_explicit_resume_grants_another_bounded_active_interval(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path / 'scratch').mkdir()
            want = threading.Event(); acknowledgements = []
            def tick(value):
                if value['execution_state'] == 'ram_paused':
                    acknowledgements.append(value['active_timeout_count'])
                    if len(acknowledgements) == 1:
                        want.clear()
            start = time.monotonic()
            result = self.backend.run_tree([sys.executable, '-c', 'import time; time.sleep(30)'],
                path, 1, 1, .25, hold=want.is_set, on_timeout=want.set, on_tick=tick,
                cancel=lambda: len(acknowledgements) >= 2 or time.monotonic()-start > 8)
            self.assertEqual(acknowledgements[:2], [1, 2])
            self.assertEqual(result['active_processes_at_return'], 0)

    def test_started_callback_failure_cleans_suspended_owned_process(self):
        import ctypes as ct
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path / 'scratch').mkdir()
            kernel, _ = self.backend.api()
            handles = []
            def started(value):
                handle = kernel.OpenProcess(0x100000, False, value['pid'])
                self.assertTrue(handle)
                handles.append(handle)
                raise RuntimeError('Injected callback failure')
            try:
                with self.assertRaisesRegex(RuntimeError, 'Injected callback failure'):
                    self.backend.run_tree([sys.executable, '-c', 'import time; time.sleep(30)'],
                        path, 1, 1, 10, on_started=started)
                self.assertEqual(kernel.WaitForSingleObject(handles[0], 1000), 0)
            finally:
                for handle in handles:
                    kernel.CloseHandle(handle)


for _index, _timeout in enumerate((.005, .025, .05, .10, .25) * 4):
    def _trial(self, timeout=_timeout):
        self.check_early_deadline(timeout)
    setattr(NativeWindowsPauseTests, f'test_early_deadline_{_index:02d}', _trial)

if __name__ == '__main__':
    unittest.main(verbosity=2)
