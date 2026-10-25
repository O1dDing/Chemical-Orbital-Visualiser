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

"""Verify/restore archived evidence and explicitly resume a bounded REF-001 subset.

No calculation starts when this file is opened, imported, verified or restored.
The run command uses the frozen scientific functions and checks their identities.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import queue
import shutil
import sys
import threading
import time
import zipfile

HERE = Path(__file__).resolve().parent
OLD_CAMPAIGN = r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906'
OLD_JOBS = r'F:\Dev\cov-cycle-20260906\jobs'
OLD_COLLECTION = r'F:\Codex\2026-09-05\branch-15\outputs\fchk-collection-20260905'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def contained(root, relative):
    if '\\' in relative or ':' in relative:
        raise ValueError('Archive names must use portable relative paths')
    root = Path(root).resolve()
    target = (root / relative).resolve()
    if target == root or not target.is_relative_to(root):
        raise ValueError('Archive entry escapes the requested directory')
    return target


def verify_assets(args):
    manifest = read(args.manifest)
    checked = 0
    for archive in manifest['archives']:
        path = args.assets / archive['name']
        if path.stat().st_size != archive['bytes'] or sha(path) != archive['sha256']:
            raise ValueError('Archive checksum mismatch: ' + archive['name'])
        with zipfile.ZipFile(path) as bundle:
            expected = {record['path']: record for record in archive['files']}
            if len(bundle.namelist()) != len(expected) or set(bundle.namelist()) != set(expected):
                raise ValueError('Archive membership mismatch')
            for relative, record in expected.items():
                target = contained(args.data_root, relative)
                digest = hashlib.sha256()
                size = 0
                output = None
                temporary = target.with_name(target.name + '.restoring')
                if args.extract:
                    if target.exists():
                        if sha(target) != record['sha256']:
                            raise ValueError('Refusing to overwrite different restored evidence: ' + relative)
                    else:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        output = temporary.open('wb')
                try:
                    with bundle.open(relative) as source:
                        for block in iter(lambda: source.read(1024 * 1024), b''):
                            digest.update(block)
                            size += len(block)
                            if output:
                                output.write(block)
                finally:
                    if output:
                        output.close()
                if size != record['bytes'] or digest.hexdigest() != record['sha256']:
                    raise ValueError('Member checksum mismatch: ' + relative)
                if output:
                    os.replace(temporary, target)
                checked += 1
        print(json.dumps({'verified_archive': archive['name'], 'verified_files': checked}), flush=True)
    print(json.dumps({'all_verified': True, 'files': checked, 'calculation_started': False}))


def relocated(value, mappings):
    if isinstance(value, dict):
        return {key: relocated(item, mappings) for key, item in value.items()}
    if isinstance(value, list):
        return [relocated(item, mappings) for item in value]
    if isinstance(value, str):
        for old, new in mappings:
            if value == old or value.startswith(old + '\\') or value.startswith(old + '/'):
                remainder = value[len(old):].lstrip('\\/').replace('\\', '/')
                return str(Path(new) / remainder)
    return value


def prepare_runtime(data, work):
    """Rebase only mutable runtime copies; never rewrite archived evidence."""
    data, work = data.resolve(), work.resolve()
    if data == work or work.is_relative_to(data) or data.is_relative_to(work):
        raise ValueError('Use separate archive-data and mutable-work directories')
    reference_manifest = data / 'campaign/references/REF-001/reference-candidates.json'
    identity = read(reference_manifest)['reference_set_identity']
    binding = work / 'restore-binding.json'
    if binding.exists():
        old = read(binding)
        if old['reference_set_identity'] != identity or old['data_root'] != str(data):
            raise ValueError('Working directory belongs to different archived data')
        return
    if work.exists() and any(work.iterdir()):
        raise ValueError('Unbound working directory must be empty')
    work.mkdir(parents=True, exist_ok=True)
    shutil.copytree(data / 'jobs/REF-001', work / 'jobs/REF-001')
    shutil.copytree(data / 'campaign/references/REF-001', work / 'references/REF-001')
    mappings = [(OLD_JOBS, work / 'jobs'), (OLD_CAMPAIGN + r'\references\REF-001', work / 'references/REF-001'),
                (OLD_CAMPAIGN, data / 'campaign'), (OLD_COLLECTION, data / 'original-collection')]
    for base in (work / 'jobs', work / 'references'):
        for path in base.rglob('*.json'):
            atomic(path, relocated(read(path), mappings))
    atomic(binding, {'reference_set_identity': identity, 'data_root': str(data),
                     'created_epoch': time.time(), 'archived_evidence_rewritten': False,
                     'path_relocation_only': True})


class PauseRequested(BaseException):
    pass


def run_incremental(args):
    if os.name != 'nt':
        raise RuntimeError('This frozen Gaussian 16W supervisor requires Windows')
    data, work = args.data_root.resolve(), args.work_root.resolve()
    # Preparation is explicit and separate so verification never launches Gaussian.
    prepare_runtime(data, work)
    source = HERE / 'runner-source'
    expected = read(source / 'identity.json')
    for name, digest in expected['files'].items():
        if sha(source / name) != digest:
            raise ValueError('Frozen runner source changed: ' + name)
    sys.path.insert(0, str(source))
    import gaussian_reference_batch as runner
    import gaussian_rebuild as rebuild
    from validation_process import physical_core_masks
    manifest = read(work / 'references/REF-001/reference-candidates.json')
    campaign = read(data / 'campaign/campaign.json')
    barrier = read(data / 'campaign/reproduction/batch.json')
    if not (barrier['all_terminal'] and barrier['case_count'] == 273 and
            barrier['case_set_identity'] == campaign['case_set_identity'] == manifest['original_case_set_identity']):
        raise ValueError('Original 273 reproduction barrier or case identity is invalid')
    computed_identity = rebuild.digest({'source_files': expected['files'],
        'reference_set': manifest['reference_set_identity'],
        'gaussian_exe': sha(args.gaussian / 'g16.exe'), 'formchk_exe': sha(args.gaussian / 'formchk.exe'),
        'timeout_hours': 24, 'round_name': 'REF-001'})
    runner_identity = read(work / 'jobs/REF-001/runner-identity.json')['identity']
    if computed_identity != runner_identity:
        raise ValueError('Gaussian executable/source/input identity differs. Preserve REF-001 and create a separately reviewed round.')
    for reference in manifest['candidates']:
        directory = work / 'references/REF-001/reference-inputs' / reference['case_id']
        if (sha(directory / 'initial.gjf') != reference['initial_input_sha256'] or
                sha(directory / 'basis.gbs') != reference['basis_sha256']):
            raise ValueError('Prepared reference input changed: ' + reference['case_id'])
    completed_stages = 0
    for path in (work / 'jobs/REF-001').glob('OLD-*/*/stage.json'):
        stage = read(path)
        if stage['status'] == 'collected':
            if sha(stage['fchk']) != stage['fchk_sha256'] or not Path(stage['checkpoint']).is_file():
                raise ValueError('Restored completed stage is missing or changed: ' + str(path))
            completed_stages += 1
    if args.command == 'check':
        print(json.dumps({'identities_verified': True, 'prepared_reference_inputs_verified': 273,
                          'completed_stages_verified': completed_stages, 'calculation_started': False}))
        return
    pause_path = work / 'PAUSE.json'
    if pause_path.exists():
        raise RuntimeError('Pause remains requested. Inspect progress, then explicitly remove PAUSE.json to resume.')
    candidates = {item['case_id']: item for item in manifest['candidates']}
    originals = {item['case_id']: item for item in campaign['cases']}
    chosen = []
    for case_id in args.case or list(candidates):
        if case_id not in candidates:
            raise ValueError('Unknown case ID: ' + case_id)
        result_path = work / 'jobs/REF-001' / case_id / 'result.json'
        previous = read(result_path) if result_path.exists() else {}
        if previous.get('status') == 'candidate_collected':
            continue
        if previous.get('status') not in (None, 'running'):
            if args.case:
                raise ValueError('A recorded failure needs post-batch adjudication and a new reviewed attempt/round: ' + case_id)
            continue
        chosen.append(case_id)
        if len(chosen) == args.limit:
            break
    if not chosen:
        print(json.dumps({'selected_cases': [], 'calculation_started': False}))
        return
    # Keep the verified Gaussian Fortran startup convention even if the caller
    # itself inherited a narrower affinity. Every calculation still receives its
    # own 4-core CPU-rate Job cap before its first instruction runs.
    import ctypes as ct
    from ctypes import wintypes as wt
    kernel = ct.WinDLL('kernel32', use_last_error=True)
    kernel.GetCurrentProcess.restype = wt.HANDLE
    kernel.GetProcessAffinityMask.argtypes = [wt.HANDLE, ct.POINTER(ct.c_size_t), ct.POINTER(ct.c_size_t)]
    kernel.SetProcessAffinityMask.argtypes = [wt.HANDLE, ct.c_size_t]
    current_mask, system_mask = ct.c_size_t(), ct.c_size_t()
    current_process = kernel.GetCurrentProcess()
    if not kernel.GetProcessAffinityMask(current_process, ct.byref(current_mask), ct.byref(system_mask)):
        raise ct.WinError(ct.get_last_error())
    if not kernel.SetProcessAffinityMask(current_process, system_mask):
        raise ct.WinError(ct.get_last_error())
    import msvcrt
    lock_path = work / 'jobs/supervisor.lock'
    with lock_path.open('a+b') as lock:
        if lock.tell() == 0:
            lock.write(b'0')
            lock.flush()
        lock.seek(0)
        msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        # A separate orchestration identity does not rewrite the frozen scientific identity.
        invocation = work / 'incremental-runs' / str(time.time_ns())
        state = {'runner_identity': runner_identity, 'orchestration_sha256': sha(__file__),
                 'selected_cases': chosen, 'workers': args.workers, 'cpu_cores_per_job': 4,
                 'tree_memory_gib_per_job': 32, 'input_memory_gib': 24,
                 'started_epoch': time.time(), 'cases': {}, 'scientific_passes': 0}
        atomic(invocation / 'progress.json', state)
        args.references = work / 'references/REF-001'
        args.jobs = work / 'jobs'
        args.output = data / 'campaign'
        args.round_name = 'REF-001'
        args.timeout_hours = 24
        original_stage = runner.collect_stage
        original_run = runner.run_tree

        def guarded_stage(*a, **kw):
            if pause_path.exists():
                raise PauseRequested()
            return original_stage(*a, **kw)

        def guarded_run(argv, *a, **kw):
            if Path(argv[0]).name.lower() == 'g16.exe' and pause_path.exists():
                raise PauseRequested()
            return original_run(argv, *a, **kw)

        runner.collect_stage, runner.run_tree = guarded_stage, guarded_run
        pending = queue.Queue()
        for case_id in chosen:
            pending.put(case_id)
        write_lock = threading.Lock()
        cores = physical_core_masks(12)
        if len(cores) < 12:
            raise ValueError('This frozen resource policy requires at least 12 physical cores')

        def worker(mask):
            while not pause_path.exists():
                try:
                    case_id = pending.get_nowait()
                except queue.Empty:
                    return
                try:
                    result = runner.collect_reference(candidates[case_id], originals[case_id], args, mask, runner_identity)
                except PauseRequested:
                    result = {'case_id': case_id, 'status': 'paused_between_stages', 'scientific_pass': False}
                with write_lock:
                    state['cases'][case_id] = result['status']
                    state['updated_epoch'] = time.time()
                    atomic(invocation / 'progress.json', state)
                    print(json.dumps({'case': case_id, 'status': result['status']}), flush=True)
                pending.task_done()

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(worker, sum(cores[4*i:4*i+4])) for i in range(args.workers)]
            try:
                for future in futures:
                    future.result()
            except KeyboardInterrupt:
                atomic(pause_path, {'requested_epoch': time.time(), 'reason': 'operator interrupt; finish current stage'})
                for future in futures:
                    future.result()
        state.update(finished_epoch=time.time(), pause_requested=pause_path.exists(),
                     selection_all_terminal=len(state['cases']) == len(chosen) and
                     all(x in ('candidate_collected', 'failed') for x in state['cases'].values()),
                     full_273_reference_batch_complete=False)
        atomic(invocation / 'progress.json', state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    verify = sub.add_parser('verify')
    verify.add_argument('--manifest', type=Path, required=True)
    verify.add_argument('--assets', type=Path, required=True)
    verify.add_argument('--data-root', type=Path, required=True)
    verify.add_argument('--extract', action='store_true')
    prepare = sub.add_parser('prepare')
    prepare.add_argument('--data-root', type=Path, required=True)
    prepare.add_argument('--work-root', type=Path, required=True)
    run = sub.add_parser('run')
    run.add_argument('--data-root', type=Path, required=True)
    run.add_argument('--work-root', type=Path, required=True)
    run.add_argument('--gaussian', type=Path, required=True)
    run.add_argument('--case', action='append')
    run.add_argument('--limit', type=int, default=1)
    run.add_argument('--workers', type=int, choices=(1, 2), default=1)
    check = sub.add_parser('check')
    check.add_argument('--data-root', type=Path, required=True)
    check.add_argument('--work-root', type=Path, required=True)
    check.add_argument('--gaussian', type=Path, required=True)
    check.set_defaults(limit=1, workers=1, case=None)
    pause = sub.add_parser('pause')
    pause.add_argument('--work-root', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'verify':
        verify_assets(args)
    elif args.command == 'prepare':
        prepare_runtime(args.data_root, args.work_root)
        print(json.dumps({'prepared': True, 'calculation_started': False}))
    elif args.command == 'pause':
        atomic(args.work_root / 'PAUSE.json', {'requested_epoch': time.time(), 'reason': 'operator requested stage-boundary pause'})
        print(json.dumps({'pause_requested': True, 'current_stage_finishes_before_exit': True}))
    else:
        if args.limit < 1 or args.limit > 273:
            parser.error('--limit must be between 1 and 273')
        run_incremental(args)


if __name__ == '__main__':
    main()
