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

"""Non-computing checks for desktop startup and explicit dispatch boundaries."""
from pathlib import Path
import json
import os
import tempfile
import time
import types
import unittest
from unittest.mock import Mock

from helper_backend import HelperBackend
from helper_deployment import fingerprint, make_plan, migrate_frozen, protected_roots, validate_target


class HelperTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.home = Path(self.temp.name).resolve()
        self.bundle = self.home / 'runtime/releases/example'
        (self.bundle / 'validation/runtime_v2').mkdir(parents=True)
        (self.bundle / 'validation/runtime_v2/resume.py').write_text('# fixture\n')
        self.config = self.home / 'runtime/config.json'
        self.config.write_text(json.dumps({'work_root': str(self.home / 'work'), 'runtime_directory': 'runtime-example'}))
        (self.home / 'runtime/active.json').write_text(json.dumps({'bundle': str(self.bundle), 'config': str(self.config)}))
        self.status = {'control': {'mode': 'pause'}, 'work_lease_held': False, 'counts': {'candidate_collected': 103}}
        self.runner = Mock(return_value=types.SimpleNamespace(returncode=0, stdout=json.dumps(self.status), stderr=''))
        self.launcher = Mock(return_value=types.SimpleNamespace(pid=123))
        self.backend = HelperBackend(self.home, runner=self.runner, launcher=self.launcher)

    def test_construction_never_starts_or_controls_work(self):
        self.assertFalse(self.backend.armed)
        self.runner.assert_not_called()
        self.launcher.assert_not_called()
        self.assertFalse((self.home / 'work').exists())

    def test_startup_status_uses_only_status_and_pinned_config(self):
        self.assertEqual(self.backend.status()['control']['mode'], 'pause')
        argv = self.runner.call_args.args[0]
        self.assertEqual(argv[-1], 'status')
        self.assertEqual(argv[-2], str(self.config))
        self.launcher.assert_not_called()

    def test_status_retains_progress_after_json_document(self):
        self.runner.return_value.stdout += '\nOLD-001 {"progress": 2}'
        self.assertIn('OLD-001', self.backend.status()['progress_text'])

    def test_all_mutating_actions_require_in_memory_enable(self):
        for command in ('hold', 'pause', 'resume', 'interrupt', 'shutdown'):
            with self.assertRaises(RuntimeError):
                self.backend.control(command)
        with self.assertRaises(RuntimeError):
            self.backend.start(1, 1)
        self.runner.assert_not_called()
        self.launcher.assert_not_called()

    def test_new_window_does_not_inherit_manual_enable(self):
        self.backend.arm(True)
        other = HelperBackend(self.home, runner=self.runner, launcher=self.launcher)
        self.assertFalse(other.armed)

    def test_invalid_selection_cannot_dispatch(self):
        self.backend.arm(True)
        for args in ((0, 1, ''), (274, 1, ''), (1, 3, ''), (1, 1, '--help'), (1, 1, 'OLD-999')):
            with self.assertRaises(ValueError):
                self.backend.start(*args)
        self.runner.assert_not_called()
        self.launcher.assert_not_called()

    def test_control_with_no_coordinator_writes_nothing(self):
        self.backend.arm(True)
        with self.assertRaises(RuntimeError):
            self.backend.control('resume')
        self.assertEqual([c.args[0][-1] for c in self.runner.call_args_list], ['status'])

    def test_explicit_start_passes_selected_cases_without_shell(self):
        self.backend.arm(True)
        self.backend.start(2, 1, 'old-001,OLD-003 OLD-001')
        call = self.launcher.call_args
        self.assertEqual(call.args[0][-8:], ['--limit', '2', '--workers', '1', '--case', 'OLD-001', '--case', 'OLD-003'])
        self.assertNotIn('shell', call.kwargs)
        self.assertFalse((self.home / 'work').exists())

    def test_existing_coordinator_prevents_second_start(self):
        self.backend.arm(True)
        self.status['work_lease_held'] = True
        self.runner.return_value.stdout = json.dumps(self.status)
        with self.assertRaises(RuntimeError):
            self.backend.start(1, 1)
        self.launcher.assert_not_called()

    def test_check_requires_explicit_non_computing_receipt(self):
        with self.assertRaises(RuntimeError):
            self.backend.check()
        self.runner.return_value.stdout = json.dumps({'calculation_started': False, 'verified_inputs': 273})
        self.assertEqual(self.backend.check()['verified_inputs'], 273)

    def test_bundle_cannot_escape_deployment(self):
        pointer = {'bundle': str(self.home.parent), 'config': str(self.config)}
        (self.home / 'runtime/active.json').write_text(json.dumps(pointer))
        with self.assertRaises(ValueError):
            HelperBackend(self.home)

    @unittest.skipUnless(os.name == 'nt', 'Hidden Tk window check requires Windows')
    def test_hidden_window_startup_and_close_are_read_only(self):
        import tkinter as tk
        from helper_ui import HelperWindow
        root = tk.Tk()
        root.withdraw()
        window = HelperWindow(root, self.backend)
        deadline = time.monotonic() + 3
        while window.busy and time.monotonic() < deadline:
            root.update()
            time.sleep(.01)
        self.assertFalse(window.busy)
        self.assertIn('空闲', window.banner.get())
        self.assertFalse(window.enabled.get())
        self.assertTrue(all('disabled' in b.state() for b in window.action_buttons))
        window.close()
        self.assertEqual([c.args[0][-1] for c in self.runner.call_args_list], ['status'])
        self.launcher.assert_not_called()


class DeploymentBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.home, self.fast = self.root / 'Resume', self.root / 'FastPause-2.1'
        self.home.mkdir()
        self.fast.mkdir()
        self.source = self.fast / 'validation/paused-20260906'
        self.source.mkdir(parents=True)
        (self.source / 'input.gjf').write_bytes(b'original frozen input\n')
        (self.source / 'job.chk').write_bytes(bytes(range(256)))
        self.destination = self.home / 'recovery/frozen-reference-20260910'
        self.plan = {'home': str(self.home), 'fast': str(self.fast), 'migration': {
            'source': str(self.source), 'destination': str(self.destination), 'files': fingerprint(self.source)}}

    def test_only_named_retired_root_can_be_removed(self):
        roots = protected_roots(self.home, self.fast)
        self.assertEqual(validate_target(self.fast, [self.home, self.fast], roots), self.fast)
        for path in (self.home, self.root):
            with self.assertRaises(ValueError):
                validate_target(path, [self.home, self.fast], roots)

    def test_calculation_directories_cannot_be_cleanup_targets(self):
        for path in protected_roots(self.home, self.fast):
            path.mkdir(exist_ok=True)
            (path / 'result.json').write_text('{}')
            with self.assertRaises(ValueError):
                validate_target(path / 'result.json', [self.home, self.fast], protected_roots(self.home, self.fast))

    def test_migration_preserves_source_and_destination_bytes(self):
        result = migrate_frozen(self.plan)
        self.assertTrue(result['verified'])
        self.assertEqual(fingerprint(self.source), self.plan['migration']['files'])
        self.assertEqual(fingerprint(self.destination), self.plan['migration']['files'])
        self.assertTrue(migrate_frozen(self.plan)['verified'])

    def test_migration_never_overwrites_different_existing_material(self):
        self.destination.mkdir(parents=True)
        (self.destination / 'job.chk').write_bytes(b'unique different checkpoint')
        with self.assertRaises(ValueError):
            migrate_frozen(self.plan)
        self.assertEqual((self.destination / 'job.chk').read_bytes(), b'unique different checkpoint')

    def test_migration_refuses_source_mutation_or_changed_destination(self):
        (self.source / 'job.chk').write_bytes(b'changed')
        with self.assertRaises(ValueError):
            migrate_frozen(self.plan)
        self.plan['migration']['destination'] = str(self.root / 'elsewhere')
        with self.assertRaises(ValueError):
            migrate_frozen(self.plan)

    def test_unclassified_old_result_blocks_cleanup_plan(self):
        old = self.home / 'runtime/releases/old'
        new = self.home / 'runtime/releases/new'
        for bundle in (old, new):
            bundle.mkdir(parents=True)
            (bundle / 'manifest.json').write_text(json.dumps({'files': {}}))
        (self.home / 'runtime/active.json').write_text(json.dumps({'bundle': str(old)}))
        unexpected = self.fast / 'validation/runtime_v2/result.json'
        unexpected.parent.mkdir()
        unexpected.write_text('{"unique_result": true}')
        with self.assertRaisesRegex(ValueError, 'Unclassified'):
            make_plan(self.home, self.fast, {'bundle': str(new)})
        self.assertTrue(unexpected.exists())


if __name__ == '__main__':
    unittest.main()
