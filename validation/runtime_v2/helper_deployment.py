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

"""Deploy the desktop Helper and retire replaced program files to Recycle Bin.

Calculation directories and frozen recovery dependencies are never cleanup targets.
An explicit plan pins every retired file and validates resolved paths before apply.
"""
from contextlib import ExitStack
from pathlib import Path
import argparse
import ctypes
from ctypes import wintypes
import hashlib
import json
import os
import shutil
import stat
import sys
import time

from fast_checkpoint import file_hash
from install_runtime import prepare
from state_store import atomic_json, lease, read_json


def fingerprint(path):
    """Reject reparse points before traversing or retiring a directory."""
    result = {}
    todo = [path]
    while todo:
        current = todo.pop()
        info = current.lstat()
        if getattr(info, 'st_file_attributes', 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT or current.is_symlink():
            raise ValueError('Reparse point cannot be retired: ' + str(current))
        if current.is_dir():
            todo.extend(current.iterdir())
        else:
            result[current.relative_to(path).as_posix()] = {'bytes': info.st_size, 'sha256': file_hash(current)}
    return result


def validate_target(path, homes, protected):
    # Validate both the literal ancestry and the final resolved absolute target.
    path = Path(path).absolute()
    if not any(path != root and path.is_relative_to(root) for root in homes) and path != homes[1]:
        raise ValueError('Cleanup target is outside the named deployment directories')
    for parent in (path, *path.parents):
        if parent.exists() and (parent.is_symlink() or getattr(parent.lstat(), 'st_file_attributes', 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT):
            raise ValueError('Cleanup ancestry contains a reparse point')
        if parent in homes:
            break
    resolved = path.resolve(strict=True)
    if resolved != path or (not any(resolved != root and resolved.is_relative_to(root) for root in homes) and resolved != homes[1]):
        raise ValueError('Resolved cleanup target escaped its deployment')
    if any(resolved.is_relative_to(root) or root.is_relative_to(resolved) for root in protected):
        raise ValueError('Cleanup target intersects protected calculation material')
    return resolved


def protected_roots(home, fast):
    return [home / name for name in ('assets', 'data', 'work', 'OLD-018-SALVAGE', 'tools', 'recovery')]


def make_plan(home, fast, prepared):
    homes, protected = [home, fast], protected_roots(home, fast)
    previous = read_json(home / 'runtime/active.json')
    old_bundle, new_bundle = Path(previous['bundle']).resolve(), Path(prepared['bundle']).resolve()
    if old_bundle == new_bundle:
        raise ValueError('This deployment is already active; inspect its existing receipt')
    new_manifest = read_json(new_bundle / 'manifest.json')['files']
    retired_bundles = []
    for bundle in (home / 'runtime/releases').iterdir():
        if bundle == new_bundle:
            continue
        validate_target(bundle, homes, protected)
        old_manifest = read_json(bundle / 'manifest.json')['files']
        expected_paths = {'validation/' + relative for relative in old_manifest} | {'manifest.json'}
        if set(fingerprint(bundle)) != expected_paths:
            raise ValueError('Unclassified material in an inactive source bundle: ' + str(bundle))
        # Frozen input bytes must remain available in the replacement bundle.
        for relative, digest in old_manifest.items():
            if file_hash(bundle / 'validation' / relative) != digest:
                raise ValueError('Inactive bundle identity changed: ' + relative)
            if relative.startswith('paused-20260906/') and new_manifest.get(relative) != digest:
                raise ValueError('Replacement changed a frozen calculation dependency: ' + relative)
        retired_bundles.append(bundle)
    paths = [(home / name, 'Replaced by the shared desktop Helper and configuration') for name in (
        'COV-REF001-Resume-Helper.zip', 'cov_resume_menu.py', 'START.cmd', 'README-先读我.md',
        'local-settings.json', 'TEST-STATUS.txt', 'helper-setup-verified.json')]
    paths += [(home / 'runtime/backups', 'Superseded program backups; originals remain recoverable in Recycle Bin'),
              (fast, 'Retired deployment; frozen calculation material is migrated and verified before retirement')]
    paths += [(bundle, 'Superseded or inactive source bundle; frozen dependencies verified in the replacement') for bundle in retired_bundles]
    rows = []
    for path, reason in paths:
        if path.exists():
            exact = validate_target(path, homes, protected)
            rows.append({'path': str(exact), 'reason': reason, 'files': fingerprint(exact)})
    frozen = fast / 'validation/paused-20260906'
    frozen_files = fingerprint(frozen)
    # Unknown additions are not silently classified as obsolete program files.
    allowed_modules = {'resume', 'policy', 'windows_job', 'windows_pause', 'native_acceptance',
                       'fast_checkpoint', 'test_runtime', 'test_fast_pause'}
    for relative in fingerprint(fast):
        path = Path(relative)
        if relative.startswith('validation/paused-20260906/') or relative.startswith('.github/'):
            continue
        if relative == 'validation/README.md':
            continue
        if path.parent.as_posix() == 'validation/runtime_v2' and (
                path.stem in allowed_modules and path.suffix == '.py' or
                path.name in ('START.cmd', 'config.json', 'config.example.json', 'README.md', 'FAST_PAUSE.md', 'WORK_PROMPT.md', '.gitignore', '.gitattributes')):
            continue
        if path.parent.as_posix() == 'validation/runtime_v2/__pycache__' and path.name.split('.')[0] in allowed_modules and path.suffix == '.pyc':
            continue
        raise ValueError('Unclassified material in retired deployment: ' + relative)
    return {'schema': 1, 'home': str(home), 'fast': str(fast), 'protected_roots': [str(p) for p in protected],
            'previous_active': previous, 'prepared': prepared, 'cleanup': rows,
            'migration': {'source': str(frozen), 'destination': str(home / 'recovery/frozen-reference-20260910'),
                          'files': frozen_files}}


def recycle_metadata(drive):
    found = set()
    try:
        folders = list((Path(drive + os.sep) / '$Recycle.Bin').iterdir())
    except OSError:
        return found
    for folder in folders:
        try:
            found.update(folder.glob('$I*'))
        except OSError:
            pass
    return found


def recycle(path):
    if os.name != 'nt':
        raise RuntimeError('Deployment cleanup requires the Windows Recycle Bin')
    before = recycle_metadata(path.drive)
    class Operation(ctypes.Structure):
        _fields_ = [('hwnd', wintypes.HWND), ('wFunc', wintypes.UINT), ('pFrom', wintypes.LPCWSTR),
                    ('pTo', wintypes.LPCWSTR), ('fFlags', wintypes.WORD), ('fAnyOperationsAborted', wintypes.BOOL),
                    ('hNameMappings', ctypes.c_void_p), ('lpszProgressTitle', wintypes.LPCWSTR)]
    source = ctypes.create_unicode_buffer(str(path) + '\0\0')
    operation = Operation()
    operation.wFunc = 3  # FO_DELETE
    operation.pFrom = ctypes.cast(source, wintypes.LPCWSTR)
    operation.fFlags = 0x40 | 0x10 | 0x4 | 0x400  # ALLOWUNDO, NOCONFIRMATION, SILENT, NOERRORUI
    shell = ctypes.WinDLL('shell32', use_last_error=True)
    shell.SHFileOperationW.argtypes = [ctypes.POINTER(Operation)]
    shell.SHFileOperationW.restype = ctypes.c_int
    result = shell.SHFileOperationW(ctypes.byref(operation))
    if result or operation.fAnyOperationsAborted or path.exists():
        raise OSError(result, 'Recycle Bin operation did not complete: ' + str(path))
    started = time.monotonic()
    # Explorer may publish its $I/$R directory receipt just after the shell call
    # completes. Wait only for that exact original path, with a bounded deadline.
    while time.monotonic() - started < 5:
        for metadata in recycle_metadata(path.drive) - before:
            try:
                data = metadata.read_bytes()
                version = int.from_bytes(data[:8], 'little')
                offset = 28 if version == 2 else 24
                original = data[offset:].decode('utf-16-le').rstrip('\0')
                if Path(original) == path:
                    content = metadata.with_name('$R' + metadata.name[2:])
                    if content.exists():
                        return {'metadata': str(metadata), 'content': str(content), 'original_path': original,
                                'receipt_wait_seconds': time.monotonic() - started}
            except (OSError, UnicodeError):
                continue
        time.sleep(.05)
    raise RuntimeError('File was retired but its Recycle Bin receipt could not be resolved: ' + str(path))


def entry_files(home, python):
    pythonw = python.with_name('pythonw.exe')
    if not pythonw.is_file():
        raise FileNotFoundError(pythonw)
    entry = home / 'COV-Helper.pyw'
    stub = ("# Shared desktop entry; opening this file only reads status.\n"
            "from pathlib import Path\nimport json, runpy, sys\n"
            "home = Path(__file__).resolve().parent\n"
            "pointer = json.loads((home / 'runtime/active.json').read_text(encoding='utf-8'))\n"
            "runtime = Path(pointer['bundle']) / 'validation/runtime_v2'\n"
            "sys.path.insert(0, str(runtime))\n"
            "sys.argv = [str(runtime / 'helper_ui.py'), '--home', str(home)] + sys.argv[1:]\n"
            "try:\n    runpy.run_path(str(runtime / 'helper_ui.py'), run_name='__main__')\n"
            "except Exception as error:\n"
            "    import tkinter as tk\n    from tkinter import messagebox\n"
            "    root = tk.Tk(); root.withdraw()\n"
            "    messagebox.showerror('COV Helper', str(error), parent=root)\n"
            "    root.destroy()\n    raise\n")
    launch = '@echo off\r\nstart "" "' + str(pythonw) + '" -B -X utf8 "' + str(entry) + '"\r\nexit /b\r\n'
    guide = ('# COV Helper 2.2.1\n\n双击 START.cmd 打开可见窗口。打开后只读取状态，计算操作每次默认关闭。\n'
             '需要手动运行时，先启用窗口中的计算控制，再选择案例上限、槽位和案例 ID；打开面板本身不会恢复队列。\n'
             '刷新状态和检查部署均不启动计算。关闭面板也不会停止已经运行的协调器。\n\n'
             'Resume 是唯一程序与配置目录；runtime/active.json 指定当前代码包，runtime/config.json 指定原位计算材料。\n'
             'data、work、assets、OLD-018-SALVAGE 与 tools 中的冻结恢复依赖保持原位。\n'
             '原独立目录中的冻结计算材料已完整迁入 recovery/frozen-reference-20260910；Resume 是唯一部署目录。\n\n'
             '清理记录位于 runtime/deployments。旧文件在 Windows 回收站，回执列出原路径、回收站位置和摘要。\n'
             '恢复旧部署前关闭 Helper 并确认没有协调器或 Gaussian；不要将回收站中的旧入口覆盖到仍在使用的新部署上。\n'
             '候选收集完成不等于科学验收通过。RAM 暂停不是磁盘保存，关机前须取得冷保存回执。\n')
    return {entry: stub.encode('utf-8'), home / 'START.cmd': launch.encode('utf-8'),
            home / 'README-先读我.md': guide.encode('utf-8')}


def migrate_frozen(plan):
    home, fast = Path(plan['home']).resolve(), Path(plan['fast']).resolve()
    migration = plan['migration']
    source, destination = Path(migration['source']), Path(migration['destination'])
    if source != fast / 'validation/paused-20260906' or destination != home / 'recovery/frozen-reference-20260910':
        raise ValueError('Frozen migration does not match the named deployment boundaries')
    if fingerprint(source) != migration['files']:
        raise ValueError('Frozen source changed after planning')
    if destination.resolve() != destination or not destination.is_relative_to(home):
        raise ValueError('Frozen migration destination escaped the Resume deployment')
    for parent in (destination, *destination.parents):
        if parent.exists() and (parent.is_symlink() or getattr(parent.lstat(), 'st_file_attributes', 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT):
            raise ValueError('Frozen migration destination contains a reparse point')
        if parent == home:
            break
    for relative, identity in migration['files'].items():
        target = destination / relative
        if not target.resolve().is_relative_to(destination) or Path(relative).is_absolute():
            raise ValueError('Invalid frozen migration member')
        if target.exists():
            if file_hash(target) != identity['sha256']:
                raise ValueError('Frozen migration would overwrite different data: ' + relative)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / relative, target)
    if fingerprint(destination) != migration['files']:
        raise ValueError('Frozen destination does not exactly match the original')
    atomic_json(home / 'recovery/frozen-reference-20260910.manifest.json', migration)
    return dict(migration, verified=True, file_count=len(migration['files']))


def apply_plan(plan, receipt_path):
    home, fast = Path(plan['home']).resolve(), Path(plan['fast']).resolve()
    prepared = plan['prepared']
    protected = protected_roots(home, fast)
    if plan['protected_roots'] != [str(p) for p in protected]:
        raise ValueError('Protected calculation boundaries changed')
    if read_json(home / 'runtime/active.json') != plan['previous_active']:
        raise ValueError('Active deployment changed after planning')
    config = read_json(prepared['config'])
    import resume
    import windows_job
    with ExitStack() as stack:
        for lock in (Path(os.environ['LOCALAPPDATA']) / 'COV/gaussian-runtime-v2.lock',
                     Path(config['work_root']) / 'jobs/supervisor.lock', home / 'helper-operation.lock'):
            stack.enter_context(lease(lock))
        data = resume.FrozenData(config)
        if windows_job.existing_gaussian(data.gaussian):
            raise RuntimeError('Gaussian is running; no deployment changes were made')
        engine = resume.Engine(config, data, windows_job)
        binding = read_json(engine.root / 'binding.json')
        if binding['runtime_identity'] != engine.identity:
            raise ValueError('Helper deployment must preserve the accepted runtime identity')
        if read_json(engine.root / 'control.json')['mode'] != 'pause':
            raise ValueError('Deployment requires the existing paused control state')
        for method in ('RPBE1PBE', 'UPBE1PBE'):
            engine.require_capability('ram_pause', method)
        manifest = read_json(Path(prepared['bundle']) / 'manifest.json')
        for name, digest in manifest['files'].items():
            if file_hash(Path(prepared['bundle']) / 'validation' / name) != digest:
                raise ValueError('Prepared deployment changed')
        for row in plan['cleanup']:
            path = validate_target(Path(row['path']), [home, fast], protected)
            if fingerprint(path) != row['files']:
                raise ValueError('Cleanup target changed after planning: ' + str(path))
        files = entry_files(home, Path(sys.executable))
        receipt = dict(plan, helper_version='2.2.1', runtime_identity=engine.identity,
                       completed=False, retired=[], new_files=[], calculation_started=False)
        atomic_json(receipt_path, receipt)
        receipt['migration'] = migrate_frozen(plan)
        atomic_json(receipt_path, receipt)
        atomic_json(home / 'runtime/active.json', prepared)
        for row in plan['cleanup']:
            path = validate_target(Path(row['path']), [home, fast], protected)
            if fingerprint(path) != row['files']:
                raise ValueError('Cleanup target changed immediately before retirement: ' + str(path))
            location = recycle(path)
            if fingerprint(Path(location['content'])) != row['files']:
                raise ValueError('Retired file verification failed: ' + str(path))
            receipt['retired'].append(dict(row, recycle=location))
            atomic_json(receipt_path, receipt)
        for path, content in files.items():
            if path.exists():
                raise ValueError('New entry unexpectedly exists: ' + str(path))
            path.write_bytes(content)
            receipt['new_files'].append({'path': str(path), 'sha256': hashlib.sha256(content).hexdigest()})
        launcher = {'helper_version': '2.2.1', 'python': str(Path(sys.executable)),
                    'pythonw': str(Path(sys.executable).with_name('pythonw.exe')),
                    'entry': str(home / 'COV-Helper.pyw'), 'startup_mode': 'read_only', 'auto_start': False,
                    'runtime_identity': engine.identity}
        atomic_json(home / 'runtime/launcher.json', launcher)
        receipt.update(completed=True, completed_epoch=time.time(), launcher=launcher)
        atomic_json(receipt_path, receipt)
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'apply'))
    parser.add_argument('--home', type=Path, required=True)
    parser.add_argument('--fast', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    home, fast = args.home.resolve(), args.fast.resolve()
    if args.command == 'plan':
        prepared = prepare(home, home / 'runtime/config.json')
        plan = make_plan(home, fast, prepared)
        atomic_json(args.plan, plan)
        print(json.dumps({'prepared': prepared, 'protected_roots': plan['protected_roots'],
                          'cleanup': [{'path': r['path'], 'files': len(r['files']), 'reason': r['reason']} for r in plan['cleanup']]}, indent=2))
    else:
        if args.receipt is None:
            parser.error('apply requires --receipt')
        plan = read_json(args.plan)
        if Path(plan['home']) != home or Path(plan['fast']) != fast:
            raise ValueError('Plan belongs to a different deployment')
        result = apply_plan(plan, args.receipt)
        print(json.dumps({'completed': result['completed'], 'helper_version': result['helper_version'],
                          'retired_targets': len(result['retired']), 'runtime_identity': result['runtime_identity'],
                          'calculation_started': result['calculation_started'], 'receipt': str(args.receipt)}, indent=2))


if __name__ == '__main__':
    main()
