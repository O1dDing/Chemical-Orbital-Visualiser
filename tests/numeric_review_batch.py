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

"""Bounded, resumable full review of an already frozen COV collection.

The numerical reviewer and its thresholds belong to the collection snapshot.
This orchestrator only schedules it and records terminal results. No Gaussian
work is launched; failures never remove a case or abort the remaining cases.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import os
import queue
import shutil
import sys
import threading
import time
import traceback

PYTHON = Path(r'F:\Dev\Python312\python.exe')
DEPS = Path(r'F:\Dev\cov-validation-20260905\reference-deps')


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def identities(directory):
    return {str(p.relative_to(directory)): {'bytes': p.stat().st_size, 'sha256': sha(p)}
            for p in sorted(directory.rglob('*')) if p.is_file()
            and '__pycache__' not in p.parts and p.suffix != '.pyc'}


def verify_files(directory, records):
    for name, record in records.items():
        path = directory / name
        if path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
            raise RuntimeError('Artifact identity changed: ' + str(path))


def barrier_check(root, manifest):
    barrier = read(root / 'collection-complete.json')
    expected = {row['case_id'] for row in manifest['cases']}
    if (not barrier['all_terminal'] or barrier['terminal_cases'] != len(expected)
            or manifest['case_count'] != len(expected)
            or {row['case_id'] for row in barrier['cases']} != expected
            or barrier['round_identity'] != manifest['round_identity']):
        raise RuntimeError('Full matching collection barrier is required before review')


def controller(root, review_root):
    sys.path.insert(0, str(root / 'runner-source'))
    from validation_process import atomic_json, physical_core_masks, run_tree
    import msvcrt
    manifest = read(root / 'manifest.json')
    barrier_check(root, manifest)
    locks = []
    try:
        # Keep the collection itself immutable, and reject duplicate reviewers.
        for path in (root / 'coordinator.lock', review_root / 'coordinator.lock'):
            stream = path.open('a+b')
            stream.seek(0)
            if stream.read(1) == b'':
                stream.write(b'0'); stream.flush()
            stream.seek(0)
            try:
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            except Exception:
                stream.close()
                raise
            locks.append(stream)
        for name, wanted in manifest['runner_sources'].items():
            if sha(root / 'runner-source' / name) != wanted:
                raise RuntimeError('Frozen reviewer dependency changed: ' + name)
        freeze = review_root / 'manifest.json'
        if freeze.exists():
            identity = read(freeze)
            if (identity['round_identity'] != manifest['round_identity']
                    or identity['collection_manifest_sha256'] != sha(root / 'manifest.json')
                    or identity['collection_barrier_sha256'] != sha(root / 'collection-complete.json')
                    or identity['orchestrator_sha256'] != sha(Path(__file__))):
                raise RuntimeError('Review resume identity mismatch')
            verify_files(DEPS, identity['reference_dependencies'])
            verify_files(PYTHON.parent, identity['python_runtime'])
        else:
            runtime_names = ('python.exe', 'python312.dll', 'python3.dll')
            runtime = {name: {'bytes': (PYTHON.parent / name).stat().st_size,
                              'sha256': sha(PYTHON.parent / name)}
                       for name in runtime_names if (PYTHON.parent / name).is_file()}
            identity = {'schema': 1, 'round_identity': manifest['round_identity'],
                        'collection_manifest_sha256': sha(root / 'manifest.json'),
                        'collection_barrier_sha256': sha(root / 'collection-complete.json'),
                        'orchestrator_sha256': sha(Path(__file__)),
                        'reviewer_sha256': manifest['runner_sources']['review_numeric_fix_case.py'],
                        'reference_dependencies_root': str(DEPS),
                        'reference_dependencies': identities(DEPS),
                        'python_runtime': runtime, 'python_version': sys.version,
                        'case_ids': [case['case_id'] for case in manifest['cases']],
                        'resources': {'workers': 4, 'cores_per_worker': 3,
                                      'gib_per_worker': 24, 'total_cores': 12, 'total_gib': 128},
                        'case_timeout_seconds': 3600, 'created_epoch': time.time(),
                        'formal_case_passes': 0}
            atomic_json(freeze, identity)
        workers = identity['resources']['workers']
        threads = identity['resources']['cores_per_worker']
        masks = physical_core_masks(workers * threads)
        if len(masks) != workers * threads:
            raise RuntimeError('Cannot allocate the frozen CPU resource slots')
        slots = queue.Queue()
        for index in range(workers):
            slots.put(sum(masks[index * threads:(index + 1) * threads]))
        mutex = threading.RLock()
        results = {}
        running = {}
        start = time.perf_counter()
        epoch = time.time()

        def update_progress():
            with mutex:
                counts = {status: sum(row['status'] == status for row in results.values())
                          for status in ('numeric_subset_pass', 'numeric_subset_fail', 'review_error')}
                atomic_json(review_root / 'progress.json', {
                    'status': 'reviewing', 'round_identity': identity['round_identity'],
                    'started_epoch': epoch, 'updated_epoch': time.time(),
                    'expected_cases': len(identity['case_ids']), 'terminal_cases': len(results),
                    'counts': counts, 'running': dict(running), 'formal_case_passes': 0,
                    'unified_adjudication': 'pending all numerical reviews terminal'})

        def run_case(case_id):
            mask = slots.get()
            attempt_start = time.perf_counter()
            case_root = review_root / 'cases' / case_id
            case_root.mkdir(parents=True, exist_ok=True)
            attempt = case_root / f'attempt-{len(list(case_root.glob("attempt-*"))) + 1:03d}'
            attempt.mkdir()
            process_dir = attempt / 'process'
            process_dir.mkdir()
            output = attempt / 'result'
            record = {'case_id': case_id, 'round_identity': identity['round_identity'],
                      'review_manifest_sha256': sha(freeze), 'attempt': str(attempt.relative_to(review_root))}
            try:
                env = dict(os.environ, COV_CPU_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
                           OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
                           PYTHONDONTWRITEBYTECODE='1')
                def on_started(data):
                    atomic_json(process_dir / 'started.json', data)
                    with mutex:
                        running[case_id] = {'pid': data['pid'], 'started_epoch': data['started_epoch'],
                                            'directory': str(attempt)}
                    update_progress()
                process = run_tree([PYTHON, '-u', root / 'runner-source' / 'review_numeric_fix_case.py',
                                    '--round', root, '--case', case_id, '--output', output],
                                   process_dir, env, mask, identity['resources']['gib_per_worker'],
                                   identity['case_timeout_seconds'], on_started=on_started)
                atomic_json(process_dir / 'process.json', process)
                record['process_exit_code'] = process['exit_code']
                record['process_wall_seconds'] = process['wall_seconds']
                report_path = output / 'review.json'
                if process['timed_out']:
                    raise RuntimeError('Independent-review process tree timed out; partial evidence retained')
                if not report_path.exists():
                    raise RuntimeError('Independent reviewer did not produce a terminal report')
                report = read(report_path)
                if (report['round_identity'] != identity['round_identity'] or report['case_id'] != case_id
                        or report['reviewer_sha256'] != identity['reviewer_sha256']):
                    raise RuntimeError('Numerical result identity differs')
                wanted_exit = {'numeric_subset_pass': 0, 'numeric_subset_fail': 2}[report['status']]
                if process['exit_code'] != wanted_exit:
                    raise RuntimeError('Reviewer process and numerical result status disagree')
                record.update(status=report['status'], review=str(report_path.relative_to(review_root)),
                              check_counts={status: sum(row['status'] == status for row in report['checks'])
                                            for status in ('pass', 'fail', 'insufficient')})
            except Exception as error:
                record.update(status='review_error', error=f'{type(error).__name__}: {error}')
                (attempt / 'exception.txt').write_text(traceback.format_exc(), encoding='utf-8')
            finally:
                slots.put(mask)
            record['artifacts'] = identities(attempt)
            record['wall_seconds_with_packaging'] = time.perf_counter() - attempt_start
            record['formal_case_pass'] = False
            atomic_json(case_root / 'terminal.json', record)
            with mutex:
                results[case_id] = record
                running.pop(case_id, None)
            update_progress()
            print(json.dumps({key: record[key] for key in ('case_id', 'status', 'wall_seconds_with_packaging')}), flush=True)

        pending = []
        for case_id in identity['case_ids']:
            terminal = review_root / 'cases' / case_id / 'terminal.json'
            if terminal.exists():
                result = read(terminal)
                if (result['round_identity'] != identity['round_identity']
                        or result['review_manifest_sha256'] != sha(freeze)):
                    raise RuntimeError('Prior terminal review identity differs')
                verify_files(review_root / result['attempt'], result['artifacts'])
                results[case_id] = result
            else:
                pending.append(case_id)
        update_progress()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(run_case, case_id) for case_id in pending]
            for future in as_completed(futures):
                future.result()
        # The shared reference installation is identified before and after the
        # complete batch. A change invalidates the batch, never an automatic pass.
        if identities(DEPS) != identity['reference_dependencies']:
            raise RuntimeError('Independent reference environment changed during review')
        verify_files(PYTHON.parent, identity['python_runtime'])
        for name, wanted in manifest['runner_sources'].items():
            if sha(root / 'runner-source' / name) != wanted:
                raise RuntimeError('Frozen reviewer dependency changed during review')
        if len(results) != len(identity['case_ids']):
            raise RuntimeError('Missing terminal numerical reviews')
        final = {'status': 'all_numeric_reviews_terminal', 'all_terminal': True,
                 'round_identity': identity['round_identity'], 'review_manifest_sha256': sha(freeze),
                 'expected_cases': len(identity['case_ids']), 'terminal_cases': len(results),
                 'counts': {status: sum(row['status'] == status for row in results.values())
                            for status in ('numeric_subset_pass', 'numeric_subset_fail', 'review_error')},
                 'wall_seconds_this_session': time.perf_counter() - start, 'finished_epoch': time.time(),
                 'formal_case_passes': 0, 'unified_adjudication': 'ready; not yet performed',
                 'cases': [results[key] for key in sorted(results)]}
        atomic_json(review_root / 'review-complete.json', final)
        atomic_json(review_root / 'progress.json', {k: v for k, v in final.items() if k != 'cases'})
        print(json.dumps({k: v for k, v in final.items() if k != 'cases'}), flush=True)
        return 0 if all(row['status'] == 'numeric_subset_pass' for row in results.values()) else 2
    finally:
        for stream in reversed(locks):
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            stream.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--round', type=Path, required=True)
    parser.add_argument('--review-id', default='review-v1')
    parser.add_argument('--controller', action='store_true')
    args = parser.parse_args()
    if not args.review_id or Path(args.review_id).name != args.review_id or args.review_id in ('.', '..'):
        raise ValueError('Review ID must be one local directory name')
    root = args.round.resolve()
    review_root = root / args.review_id
    if args.controller:
        return controller(root, review_root)
    barrier_check(root, read(root / 'manifest.json'))
    review_root.mkdir(parents=True, exist_ok=True)
    frozen = review_root / 'orchestrator.py'
    if frozen.exists():
        if sha(frozen) != sha(Path(__file__)):
            raise RuntimeError('Existing review uses a different frozen orchestrator')
    else:
        shutil.copy2(Path(__file__), frozen)
    sys.path.insert(0, str(root / 'runner-source'))
    from validation_process import atomic_json, physical_core_masks, run_tree
    directory = review_root / ('supervisor-' + str(time.time_ns()))
    directory.mkdir()
    env = dict(os.environ, COV_CPU_THREADS='3', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    result = run_tree([PYTHON, '-u', frozen, '--round', root, '--review-id', args.review_id, '--controller'],
                      directory, env, sum(physical_core_masks(12)), 128, 7 * 86400)
    atomic_json(directory / 'process.json', result)
    print(json.dumps({'exit': result['exit_code'], 'seconds': result['wall_seconds'], 'evidence': str(directory)}), flush=True)
    return result['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())
