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

"""Prepare an immutable bundle, then activate compatible local entry points.

All replaced files are copied and hashed before replacement. No data, result,
checkpoint or source archive is removed. Rollback refuses to overwrite later edits.
"""
from contextlib import ExitStack
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import sys
import time

from state_store import atomic_json, read_json, lease, file_transaction
from fast_checkpoint import file_hash

HERE = Path(__file__).resolve().parent


def bundle_files():
    files = {}
    for directory in (HERE, HERE.parent / 'paused-20260906', HERE.parent / 'runtime-integration-20260913'):
        if not directory.is_dir():
            continue
        for path in directory.rglob('*'):
            if path.is_file() and path.name != 'config.json' and not any(
                    part in ('__pycache__', '.git') for part in path.parts) and path.suffix != '.pyc':
                files[path.relative_to(HERE.parent).as_posix()] = file_hash(path)
    return files


def prepare(home, config):
    files = bundle_files()
    identity = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    bundle = home / 'runtime/releases' / identity
    for relative, checksum in files.items():
        target = bundle / 'validation' / relative
        if target.exists():
            if file_hash(target) != checksum:
                raise ValueError('Existing immutable bundle changed: ' + str(target))
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(HERE.parent / relative, target)
            if file_hash(target) != checksum:
                raise ValueError('Bundle copy failed verification')
    atomic_json(bundle / 'manifest.json', {'bundle_identity': identity, 'files': files})
    central = home / 'runtime/config.json'
    desired = read_json(config)
    if central.exists() and read_json(central) != desired:
        raise ValueError('Existing shared configuration differs; review it explicitly')
    if not central.exists():
        atomic_json(central, desired)
    return {'bundle_identity': identity, 'bundle': str(bundle), 'config': str(central)}


def stub(home, module, legacy=False):
    pointer = str(home / 'runtime/active.json')
    return ("# Compatible entry point for the COV unified runtime.\n"
        "import json, runpy, sys\nfrom pathlib import Path\n"
        f"_pointer = json.loads(Path({pointer!r}).read_text(encoding='utf-8'))\n"
        "_runtime = Path(_pointer['bundle']) / 'validation/runtime_v2'\n"
        "sys.path.insert(0, str(_runtime))\n"
        "if __name__ == '__main__':\n"
        + ("    sys.argv = [str(_runtime / 'entry.py'), _pointer['config'], '--legacy'] + sys.argv[1:]\n"
           "    runpy.run_path(str(_runtime / 'entry.py'), run_name='__main__')\n" if legacy else
           f"    sys.argv = [str(_runtime / {module!r})] + "
           "(['--config', _pointer['config']] if " + repr(module) + " == 'resume.py' and '--config' not in sys.argv else []) + sys.argv[1:]\n"
           f"    runpy.run_path(str(_runtime / {module!r}), run_name='__main__')\n")
        + "else:\n"
        + f"    globals().update(runpy.run_path(str(_runtime / {module!r}), run_name=__name__))\n")


def activate(home, fast, prepared):
    config = read_json(prepared['config'])
    work = Path(config['work_root']).resolve()
    import resume
    data = resume.FrozenData(config)
    import windows_job
    with ExitStack() as stack:
        stack.enter_context(lease(Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'COV/gaussian-runtime-v2.lock'))
        stack.enter_context(lease(work / 'jobs/supervisor.lock'))
        stack.enter_context(lease(home / 'helper-operation.lock'))
        if windows_job.existing_gaussian(data.gaussian):
            raise RuntimeError('Gaussian still active; no entry-point replacement')
        engine = resume.Engine(config, data, windows_job)
        binding = read_json(engine.root / 'binding.json')
        if binding['runtime_identity'] != engine.identity or not (engine.root / 'import-summary.json').exists():
            raise ValueError('Prepare, validate and import the new runtime before activating entries')
        for method in ('RPBE1PBE', 'UPBE1PBE'):
            engine.require_capability('ram_pause', method)
        original = fast / 'validation/runtime_v2'
        changes = {home / 'cov_resume_menu.py': stub(home, 'resume.py', legacy=True).encode('utf-8')}
        for path in original.glob('*.py'):
            if (HERE / path.name).exists():
                changes[path] = stub(home, path.name).encode('utf-8')
        changes[original / 'config.json'] = (json.dumps(config, indent=2) + '\n').encode()
        note = ('# COV 统一运行时 2.2\n\n此入口使用 Resume/runtime/active.json 指向的同一套运行时和配置。\n'
                '操作说明见当前 bundle/validation/runtime_v2/README.md。\n'
                '原入口及代码已逐文件保存在 Resume/runtime/backups；历史结果和恢复材料保持原位。\n')
        for name in ('README.md', 'FAST_PAUSE.md', 'WORK_PROMPT.md'):
            changes[original / name] = note.encode('utf-8')
        changes[home / 'README-先读我.md'] = note.encode('utf-8')
        backup = home / 'runtime/backups' / str(time.time_ns())
        backup.mkdir(parents=True)
        rows = []
        for number, (path, content) in enumerate(changes.items()):
            saved = backup / f'{number:03d}-{path.name}'
            before = None
            if path.exists():
                before = file_hash(path); shutil.copy2(path, saved)
                if file_hash(saved) != before:
                    raise ValueError('Entry backup mismatch')
            rows.append({'path': str(path), 'backup': str(saved), 'before_sha256': before,
                         'after_sha256': hashlib.sha256(content).hexdigest()})
        active = home / 'runtime/active.json'
        previous = read_json(active) if active.exists() else None
        receipt = {'prepared': prepared, 'previous_active': previous, 'replacements': rows, 'data_removed': False}
        atomic_json(backup / 'replacement-receipt.json', receipt)
        # Publish the complete, verified pointer before any stub can use it.
        atomic_json(active, prepared)
        for path, content in changes.items():
            temporary = path.with_name(path.name + '.' + str(time.time_ns()) + '.tmp')
            with temporary.open('xb') as stream:
                stream.write(content); stream.flush(); os.fsync(stream.fileno())
            try:
                with file_transaction(path):
                    os.replace(temporary, path)
            finally:
                temporary.unlink(missing_ok=True)
        for row in rows:
            if file_hash(row['path']) != row['after_sha256']:
                raise ValueError('Entry replacement mismatch')
        receipt['verified'] = True
        atomic_json(backup / 'replacement-receipt.json', receipt)
        return {'active': prepared, 'rollback_receipt': str(backup / 'replacement-receipt.json'), 'replaced_files': len(rows)}


def rollback(home, receipt_path):
    receipt = read_json(receipt_path)
    config = read_json(receipt['prepared']['config'])
    work = Path(config['work_root']).resolve()
    import resume
    data = resume.FrozenData(config)
    import windows_job
    with ExitStack() as stack:
        for lock in (Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'COV/gaussian-runtime-v2.lock',
                     work / 'jobs/supervisor.lock', home / 'helper-operation.lock'):
            stack.enter_context(lease(lock))
        if windows_job.existing_gaussian(data.gaussian):
            raise RuntimeError('Gaussian still active; no rollback')
        for row in receipt['replacements']:
            if file_hash(row['path']) not in (row['after_sha256'], row['before_sha256']):
                raise ValueError('Entry changed since activation: ' + row['path'])
            if row['before_sha256'] is None or file_hash(row['backup']) != row['before_sha256']:
                raise ValueError('No verified original to restore')
        for row in receipt['replacements']:
            shutil.copy2(row['backup'], row['path'])
            if file_hash(row['path']) != row['before_sha256']:
                raise ValueError('Rollback copy mismatch')
        if receipt['previous_active']:
            atomic_json(home / 'runtime/active.json', receipt['previous_active'])
        return {'restored_files': len(receipt['replacements']), 'results_retained': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'activate', 'rollback'))
    parser.add_argument('--home', type=Path, required=True)
    parser.add_argument('--fast', type=Path)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    if args.command == 'rollback':
        if not args.receipt:
            parser.error('rollback requires --receipt')
        print(json.dumps(rollback(args.home.resolve(), args.receipt.resolve()), indent=2))
        return
    prepared = prepare(args.home.resolve(), args.config.resolve())
    if args.command == 'activate':
        if not args.fast:
            parser.error('activate requires --fast')
        prepared = activate(args.home.resolve(), args.fast.resolve(), prepared)
    print(json.dumps(prepared, indent=2))


if __name__ == '__main__':
    main()
