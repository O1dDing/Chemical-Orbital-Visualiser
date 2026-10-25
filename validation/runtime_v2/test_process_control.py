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

"""Deterministic ownership, cancellation and persistence regression tests."""
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest import mock

from process_control import PauseTransaction, PauseDeadline
import state_store as state


class Clock:
    now = 0.
    def __call__(self): return self.now
    def sleep(self, seconds): self.now += seconds


class Threads:
    def __init__(self, clock):
        self.clock = clock
        self.processes = {1}
        self.live = {(1, 11, 100), (1, 12, 101)}
        self.counts = {key: 0 for key in self.live}
        self.handles = {}
        self.serial = 0
        self.suspend_calls = []
        self.resume_calls = []
        self.fail_suspend = None
        self.fail_resume = set()
        self.pin_seconds = 0.

    def pids(self, guard): guard('pids'); return set(self.processes)
    def threads(self, pids, guard):
        for key in sorted(self.live):
            guard('enumeration')
            yield key[:2]
    def pin(self, pid, tid, guard):
        self.clock.sleep(self.pin_seconds)
        guard('pin')
        key = next((k for k in self.live if k[:2] == (pid, tid)), None)
        if key is None: return None
        self.serial += 1
        self.handles[self.serial] = key
        return key, self.serial
    def exited(self, handle): return self.handles[handle] not in self.live
    def suspend(self, handle):
        key = self.handles[handle]
        if key == self.fail_suspend: raise OSError('injected suspension failure')
        self.counts[key] = self.counts.get(key, 0) + 1
        self.suspend_calls.append(key)
        return True
    def resume(self, handle):
        key = self.handles[handle]
        if key in self.fail_resume: return False
        self.counts[key] -= 1
        self.resume_calls.append(key)
        return True
    def close(self, handle): del self.handles[handle]


class PauseTests(unittest.TestCase):
    def setUp(self):
        self.clock = Clock()
        self.api = Threads(self.clock)
        self.pause = PauseTransaction(self.api, clock=self.clock, sleep=self.clock.sleep)

    def test_each_owned_increment_is_undone_once(self):
        self.api.counts[(1, 11, 100)] = 2
        self.assertEqual(self.pause.hold(), 2)
        self.assertEqual(self.pause.hold(refresh=True), 2)
        self.assertEqual(len(self.api.suspend_calls), 2)
        self.pause.resume()
        self.assertEqual(self.api.counts, {(1, 11, 100): 2, (1, 12, 101): 0})
        self.assertFalse(self.api.handles)
        self.assertEqual(self.pause.state, 'running')

    def test_missing_process_threads_prevent_false_pause_confirmation(self):
        self.api.processes.add(2)
        with self.assertRaises(PauseDeadline): self.pause.hold(timeout=.2)
        self.assertFalse(self.pause.held)
        self.assertTrue(all(n == 0 for n in self.api.counts.values()))

    def test_deadline_is_checked_inside_owned_thread_iteration(self):
        self.api.live = {(1, i, i) for i in range(100, 200)}
        self.api.pin_seconds = .03
        with self.assertRaises(PauseDeadline): self.pause.hold(timeout=.1)
        self.assertLess(self.clock.now, .14)
        self.assertLess(len(self.api.suspend_calls), 4)
        self.assertTrue(all(n == 0 for n in self.api.counts.values()))
        self.assertFalse(self.api.handles)

    def test_cancellation_interrupts_enumeration(self):
        def slow_scan(pids, guard):
            for _ in range(1000):
                self.clock.sleep(.01)
                guard('enumeration')
            return iter([])
        self.api.threads = slow_scan
        self.assertEqual(self.pause.hold(cancel=lambda: self.clock.now >= .08), 0)
        self.assertLess(self.clock.now, .12)
        self.assertEqual(self.pause.state, 'running')

    def test_cancel_after_partial_pause_rolls_back(self):
        self.api.pin_seconds = .03
        self.assertEqual(self.pause.hold(cancel=lambda: bool(self.api.suspend_calls)), 0)
        self.assertTrue(self.api.suspend_calls)
        self.assertTrue(all(n == 0 for n in self.api.counts.values()))
        self.assertFalse(self.api.handles)

    def test_failure_and_callback_failure_roll_back(self):
        self.api.fail_suspend = (1, 12, 101)
        with self.assertRaises(OSError): self.pause.hold()
        self.assertEqual(self.api.counts[(1, 11, 100)], 0)
        self.api.fail_suspend = None
        with self.assertRaisesRegex(RuntimeError, 'callback'):
            self.pause.hold(on_progress=lambda _: (_ for _ in ()).throw(RuntimeError('callback')))
        self.assertTrue(all(n == 0 for n in self.api.counts.values()))

    def test_exit_and_reused_thread_id_do_not_reuse_old_handle(self):
        self.pause.hold()
        self.api.live.remove((1, 11, 100))
        self.api.live.add((1, 11, 200))
        self.api.counts[(1, 11, 200)] = 3
        self.pause.hold(refresh=True)
        self.assertNotIn((1, 11, 100), self.pause.handles)
        self.pause.resume()
        self.assertEqual(self.api.counts[(1, 11, 200)], 3)
        self.assertNotIn((1, 11, 100), self.api.resume_calls)

    def test_failed_rollback_retains_handle_for_owner_cleanup(self):
        self.api.fail_suspend = (1, 12, 101)
        self.api.fail_resume.add((1, 11, 100))
        with self.assertRaisesRegex(RuntimeError, 'undo owned suspend'):
            self.pause.hold()
        self.assertTrue(self.pause.held)
        self.assertEqual(self.pause.state, 'resume_failed')
        self.assertEqual(len(self.pause.handles), 1)

    def test_empty_job_is_exited_not_paused(self):
        self.api.processes.clear()
        self.api.live.clear()
        self.assertEqual(self.pause.hold(), 0)
        self.assertFalse(self.pause.held)
        self.assertEqual(self.pause.state, 'exited')

    def test_deadline_after_successful_suspend_keeps_rollback_pair(self):
        original = self.api.suspend
        def slow(handle):
            result = original(handle)
            self.clock.sleep(.2)
            return result
        self.api.suspend = slow
        with self.assertRaises(PauseDeadline): self.pause.hold(timeout=.1)
        self.assertEqual(len(self.api.suspend_calls), 1)
        self.assertEqual(self.api.suspend_calls, self.api.resume_calls)


class StoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / 'state.json'

    def test_reader_and_replacement_share_one_transaction_boundary(self):
        state.atomic_json(self.path, {'version': 1})
        opened, release = threading.Event(), threading.Event()
        rows, errors = [], []
        real_read = Path.read_text
        def slow_read(path, *args, **kwargs):
            value = real_read(path, *args, **kwargs)
            opened.set()
            if not release.wait(3): raise TimeoutError('test reader release missing')
            return value
        def reader():
            try: rows.append(state.read_json(self.path))
            except BaseException as error: errors.append(error)
        written = threading.Event()
        def writer():
            try:
                state.atomic_json(self.path, {'version': 2})
                written.set()
            except BaseException as error: errors.append(error)
        with mock.patch.object(Path, 'read_text', slow_read):
            thread = threading.Thread(target=reader)
            thread.start()
            write_thread = None
            try:
                self.assertTrue(opened.wait(2))
                write_thread = threading.Thread(target=writer)
                write_thread.start()
                if state.os.name == 'nt':
                    self.assertFalse(written.wait(.05))
            finally:
                release.set()
                thread.join(3)
                if write_thread is not None: write_thread.join(3)
        self.assertFalse(errors)
        self.assertTrue(written.is_set())
        self.assertEqual(rows, [{'version': 1}])
        self.assertEqual(state.read_json(self.path), {'version': 2})

    def test_failed_replace_preserves_old_document(self):
        state.atomic_json(self.path, {'version': 1})
        with mock.patch.object(state.os, 'replace', side_effect=PermissionError('persistent denial')):
            with self.assertRaises(PermissionError): state.atomic_json(self.path, {'version': 2})
        self.assertEqual(state.read_json(self.path), {'version': 1})
        self.assertFalse(list(self.root.glob('*.tmp')))

    def test_invalid_or_missing_document_remains_an_error(self):
        with self.assertRaises(FileNotFoundError): state.read_json(self.path)
        self.path.write_text('{invalid')
        with self.assertRaises(json.JSONDecodeError): state.read_json(self.path)
        state.atomic_json(self.root / 'control.json', {'mode': 'invalid'})
        with self.assertRaises(ValueError): state.ControlStore(self.root).read()

    def test_concurrent_commands_are_serialized_and_interrupt_is_sticky(self):
        control = state.ControlStore(self.root)
        interrupted = control.set('interrupt')
        errors = []
        def commands():
            try:
                for _ in range(10): control.set('run')
            except BaseException as error: errors.append(error)
        threads = [threading.Thread(target=commands) for _ in range(4)]
        for thread in threads: thread.start()
        for thread in threads: thread.join(5)
        self.assertFalse(errors)
        self.assertFalse(any(t.is_alive() for t in threads))
        current = control.read()
        self.assertEqual(current['sequence'], 41)
        self.assertEqual(current['interrupt_epoch'], interrupted['interrupt_epoch'])
        self.assertEqual(current['interrupt_sequence'], 1)
        self.assertEqual(len(list((self.root / 'commands').glob('*.json'))), 41)

    def test_commands_and_readers_in_separate_processes(self):
        import subprocess
        import sys
        script = ('import sys; from pathlib import Path; from state_store import ControlStore; '
                  's=ControlStore(Path(sys.argv[1])); '
                  '[(s.set("pause"),s.read()) for _ in range(15)]')
        state.ControlStore(self.root).set('interrupt')
        children = [subprocess.Popen([sys.executable, '-B', '-c', script, str(self.root)],
            cwd=Path(__file__).resolve().parent) for _ in range(3)]
        try:
            self.assertEqual([p.wait(timeout=8) for p in children], [0, 0, 0])
        finally:
            for process in children:
                if process.poll() is None:
                    process.terminate(); process.wait(timeout=3)
        current = state.ControlStore(self.root).read()
        self.assertEqual(current['sequence'], 46)
        self.assertEqual(current['interrupt_sequence'], 1)

    def test_independent_readers_racing_document_replacement(self):
        import subprocess
        import sys
        state.atomic_json(self.path, {'sequence': 0})
        script = ('import sys; from state_store import read_json; '
                  '[read_json(sys.argv[1]) for _ in range(1500)]')
        children = [subprocess.Popen([sys.executable, '-B', '-c', script, str(self.path)],
            cwd=Path(__file__).resolve().parent) for _ in range(2)]
        try:
            for i in range(600):
                state.atomic_json(self.path, {'sequence': i+1})
            self.assertEqual([p.wait(timeout=8) for p in children], [0, 0])
        finally:
            for process in children:
                if process.poll() is None:
                    process.terminate(); process.wait(timeout=3)
        self.assertEqual(state.read_json(self.path), {'sequence': 600})


if __name__ == '__main__':
    unittest.main(verbosity=2)
