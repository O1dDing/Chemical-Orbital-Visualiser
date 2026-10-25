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

"""Real Windows children exercise the opt-in finite build-driver policy."""
from __future__ import annotations
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from validation_process import atomic_json, run_tree


@unittest.skipUnless(os.name == 'nt', 'Windows Job Objects are required')
class FiniteDriverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        specified = os.environ.get('COV_FINITE_DRIVER_TEST_ROOT')
        cls.root = Path(specified).resolve() if specified else Path(
            tempfile.mkdtemp(prefix='cov-finite-driver-')).resolve()
        cls.root.mkdir(parents=True, exist_ok=True)
        cls.program = cls.root / 'driver and child.py'
        cls.program.write_text('''import json, os, pathlib, subprocess, sys, time
directory = pathlib.Path(sys.argv[2])
if sys.argv[1] == 'child':
    (directory/'child-ready.json').write_text(json.dumps({'pid':os.getpid()}))
    time.sleep(1.5 if sys.argv[3] == 'work' else 300)
    (directory/'child-completed.json').write_text(json.dumps({'pid':os.getpid(), 'epoch':time.time()}))
else:
    child = subprocess.Popen([sys.executable, '-X', 'utf8', __file__, 'child', str(directory), sys.argv[3]],
                             stdin=subprocess.DEVNULL)
    children = [(child,directory)]
    if len(sys.argv)>5:
        worker_directory = directory/'undeclared worker'
        worker_directory.mkdir()
        worker_env = dict(os.environ,PYTHONHOME=sys.base_prefix,
                          PATH=sys.base_prefix+os.pathsep+os.environ.get('PATH',''))
        worker = subprocess.Popen([sys.argv[5],'-X','utf8',__file__,'child',str(worker_directory),'work'],
                                  stdin=subprocess.DEVNULL,env=worker_env)
        children.append((worker,worker_directory))
    while not all((path/'child-ready.json').exists() for _,path in children):
        if any(process.poll() is not None for process,_ in children):
            raise SystemExit(91)
        time.sleep(.01)
    raise SystemExit(int(sys.argv[4]))
''', encoding='utf-8')
        # Same basename at a different full path; one test runs it as work that
        # must finish normally while the declared helper is cleaned separately.
        copy_directory = cls.root / 'different image directory'
        copy_directory.mkdir(exist_ok=True)
        cls.other_image_path = copy_directory / Path(sys.executable).name
        shutil.copyfile(sys.executable, cls.other_image_path)
        for pattern in ('python3*.dll', 'vcruntime*.dll'):
            for library in Path(sys.base_prefix).glob(pattern):
                shutil.copyfile(library, copy_directory/library.name)

    def run_driver(self, name, mode='work', root_code=0, helpers=(), worker_image=None):
        directory = self.root / name
        directory.mkdir(exist_ok=False)
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', PYTHONIOENCODING='utf-8')
        def observe(event):
            with (directory/'phase-events.jsonl').open('a', encoding='utf-8') as stream:
                stream.write(json.dumps(event, allow_nan=False)+'\n')
        command = [sys.executable, '-X', 'utf8', str(self.program),
                   'driver', str(directory), mode, str(root_code)]
        if worker_image is not None:
            command.append(str(worker_image))
        result = run_tree(command,
                          directory, env, 1, 4, 120, affinity=False,
                          finite_driver_helpers=helpers, on_phase=observe)
        atomic_json(directory/'process.json', result)
        self.assertFalse(result['timed_out'], f'Real control timed out: {directory}')
        self.assertEqual(result['exit_code'], root_code)
        self.assertGreaterEqual(result['tree_process_count'], 2)
        events = result['lifecycle_events']
        names = [event['phase'] for event in events]
        self.assertLess(names.index('root_exit_observed'), names.index('tree_empty'))
        self.assertLess(names.index('tree_empty'), names.index('cleanup_complete'))
        self.assertGreater(result['root_exit_observation']['active_job_processes'], 0,
                           'The child must actually outlive its root for this control to count')
        return result, directory

    def assert_targeted_cleanup(self, result, directory):
        self.assertEqual(result['completion_policy'], 'finite_driver_declared_helpers')
        self.assertFalse((directory/'child-completed.json').exists())
        child = json.loads((directory/'child-ready.json').read_text())
        records = result['terminated_declared_helpers']
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['pid'], child['pid'])
        self.assertTrue(records[0]['verified_current_job_member'])
        self.assertEqual(os.path.normcase(records[0]['image_path']), os.path.normcase(sys.executable))
        self.assertGreater(records[0]['creation_filetime'], 0)
        self.assertEqual(records[0]['termination_exit_code'], 125)
        names = [event['phase'] for event in result['lifecycle_events']]
        self.assertLess(names.index('root_exit_observed'), names.index('finite_driver_cleanup_enter'))
        self.assertLess(names.index('finite_driver_helper_termination_requested'), names.index('tree_empty'))

    def test_default_waits_for_scientific_child_after_root_exit(self):
        result, directory = self.run_driver('default waits')
        self.assertEqual(result['completion_policy'], 'entire_tree')
        self.assertEqual(result['terminated_declared_helpers'], [])
        self.assertTrue((directory/'child-completed.json').exists())

    def test_declared_helper_cleanup_preserves_success_and_outside_process(self):
        # The sentinel uses the exact same image but is outside run_tree's Job.
        sentinel = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(300)'],
                                    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL,
                                    creationflags=subprocess.CREATE_NO_WINDOW)
        try:
            result, directory = self.run_driver('success and outside sentinel', mode='hold',
                                                helpers=(Path(sys.executable),))
            self.assert_targeted_cleanup(result, directory)
            self.assertIsNone(sentinel.poll(), 'A same-image process outside this Job must remain alive')
            atomic_json(directory/'outside-sentinel.json', {'pid':sentinel.pid,
                        'alive_after_inner_job_empty':True, 'cleanup':'owned Popen handle in finally'})
        finally:
            if sentinel.poll() is None:
                sentinel.terminate()
            sentinel.wait(timeout=10)

    def test_same_basename_at_another_path_does_not_authorize_cleanup(self):
        result, directory = self.run_driver('undeclared same basename', helpers=(self.other_image_path,))
        self.assertEqual(result['terminated_declared_helpers'], [])
        self.assertTrue((directory/'child-completed.json').exists())
        blocked = [event for event in result['lifecycle_events']
                   if event['phase'] == 'finite_driver_cleanup_deferred']
        self.assertTrue(blocked)
        self.assertTrue(any(member['status'] == 'undeclared_image'
                            for event in blocked for member in event['members']))

    def test_failed_driver_remains_failed_after_declared_helper_cleanup(self):
        result, directory = self.run_driver('failed driver', mode='hold', root_code=7,
                                            helpers=(Path(sys.executable),))
        self.assert_targeted_cleanup(result, directory)

    def test_undeclared_worker_finishes_after_targeted_helper_cleanup(self):
        result, directory = self.run_driver('helper plus continuing worker', mode='hold',
                                            helpers=(Path(sys.executable),), worker_image=self.other_image_path)
        self.assert_targeted_cleanup(result, directory)
        worker = json.loads((directory/'undeclared worker/child-completed.json').read_text())
        cleanup = next(event for event in result['lifecycle_events']
                       if event['phase']=='finite_driver_helper_termination_requested')
        empty = next(event for event in result['lifecycle_events'] if event['phase']=='tree_empty')
        self.assertLess(cleanup['epoch'], worker['epoch'])
        self.assertLessEqual(worker['epoch'], empty['epoch'])
        blocked = [member for event in result['lifecycle_events']
                   if event['phase']=='finite_driver_cleanup_deferred' for member in event['members']]
        self.assertTrue(any(member['pid']==worker['pid'] and member['status']=='undeclared_image'
                            for member in blocked))

    def test_exit_259_is_an_exited_root_not_a_running_process(self):
        result, directory = self.run_driver('real exit 259', mode='hold', root_code=259,
                                            helpers=(Path(sys.executable),))
        self.assert_targeted_cleanup(result, directory)


if __name__ == '__main__':
    unittest.main()
