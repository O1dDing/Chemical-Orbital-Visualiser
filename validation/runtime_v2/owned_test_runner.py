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

"""Run Python regressions inside a separately owned, bounded Windows Job.

No Gaussian is launched by this runner. A wedged test cannot outlive its guard.
The deliberate watchdog self-test only terminates its own worker/descendant.
"""
from pathlib import Path
import argparse
import ctypes as ct
from ctypes import wintypes as wt
import json
import os
import subprocess
import sys
import time
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.append(str(HERE.parent / 'paused-20260906/runner-source'))


def worker(folder, names, hang=False):
    from state_store import atomic_json
    limit = time.monotonic() + 10
    while not (folder / 'admitted').exists():
        if time.monotonic() >= limit:
            raise TimeoutError('No guard admission')
        time.sleep(.01)
    if hang:
        subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])
        time.sleep(60)
        return 1
    class Result(unittest.TextTestResult):
        def startTest(self, test):
            atomic_json(folder / 'current-test.json', {'sequence': self.testsRun + 1, 'test': test.id()})
            super().startTest(test)
    suite = unittest.defaultTestLoader.loadTestsFromNames(names)
    result = unittest.TextTestRunner(verbosity=2, resultclass=Result).run(suite)
    return 0 if result.wasSuccessful() else 1


def guarded(folder, names, deadline=20, hang=False):
    from state_store import atomic_json, read_json
    import windows_job
    folder.mkdir(parents=True, exist_ok=False)
    k, v = windows_job.api()
    job = v.require(k.CreateJobObjectW(None, None))
    process = None
    started = time.monotonic()
    try:
        limits = v.EXTENDEDLIMIT()
        limits.BasicLimitInformation.LimitFlags = 0x2000
        v.require(k.SetInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits)))
        with (folder / 'tests.log').open('wb') as log:
            argv = [sys.executable, '-B', '-X', 'utf8', str(Path(__file__).resolve()),
                    '--worker', '--output', str(folder)] + (['--self-test'] if hang else []) + names
            process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                       creationflags=0x8)
            v.require(k.AssignProcessToJobObject(job, int(process._handle)))
            (folder / 'admitted').write_text('Owned by test guard\n')
            sequence, test_started, timed_out = None, time.monotonic(), False
            while process.poll() is None:
                current = folder / 'current-test.json'
                value = read_json(current, None)
                if value is not None:
                    if value['sequence'] != sequence:
                        sequence, test_started = value['sequence'], time.monotonic()
                if time.monotonic()-test_started > deadline or time.monotonic()-started > 600:
                    timed_out = True
                    v.require(k.TerminateJobObject(job, 124))
                    break
                time.sleep(.05)
            process.wait(timeout=5)
            windows_job._drain_owned_tree(k, v, job)
            accounting = v.ACCOUNTING()
            v.require(k.QueryInformationJobObject(job, 1, ct.byref(accounting), ct.sizeof(accounting), None))
            last_test = read_json(folder / 'current-test.json', None)
            if last_test is not None:
                sequence = last_test['sequence']
            result = {'exit_code': process.returncode, 'watchdog_timeout': timed_out,
                'wall_seconds': time.monotonic()-started, 'tests_started': sequence,
                'remaining_owned_processes': accounting.ActiveProcesses,
                'self_test': hang, 'test_log': str(folder / 'tests.log')}
            result['passed'] = accounting.ActiveProcesses == 0 and (
                timed_out and accounting.TotalProcesses >= 2 if hang else not timed_out and process.returncode == 0)
            atomic_json(folder / 'summary.json', result)
            return result
    finally:
        k.CloseHandle(job)
        if process is not None and process.poll() is None:
            process.wait(timeout=5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--deadline', type=float, default=20)
    parser.add_argument('tests', nargs='*', default=[])
    args = parser.parse_args()
    names = args.tests or ['test_runtime', 'test_fast_pause', 'test_process_control', 'test_migration', 'test_entry', 'test_helper', 'test_resource_scheduler']
    if args.worker:
        return worker(args.output, names, args.self_test)
    if os.name != 'nt':
        return 0 if unittest.TextTestRunner(verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromNames(names)).wasSuccessful() else 1
    result = guarded(args.output.resolve(), names, args.deadline, args.self_test)
    print(json.dumps(result, indent=2))
    if not result['passed']:
        print((args.output / 'tests.log').read_text(encoding='utf-8', errors='replace'))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
