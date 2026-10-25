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

"""Durable runtime state and coordinator leases, shared by every entry point."""
from __future__ import annotations

from contextlib import contextmanager
import json
import hashlib
import os
from pathlib import Path
import time
import uuid

MODES = frozenset(('run', 'pause', 'hold', 'interrupt', 'shutdown'))


@contextmanager
def file_transaction(path):
    """Serialize Windows readers and replacement without writing beside frozen data."""
    if os.name != 'nt':
        yield
        return
    import ctypes as ct
    from ctypes import wintypes as wt
    kernel = ct.WinDLL('kernel32', use_last_error=True)
    kernel.CreateMutexW.argtypes = [ct.c_void_p, wt.BOOL, wt.LPCWSTR]
    kernel.CreateMutexW.restype = wt.HANDLE
    kernel.WaitForSingleObject.argtypes = [wt.HANDLE, wt.DWORD]
    kernel.WaitForSingleObject.restype = wt.DWORD
    kernel.ReleaseMutex.argtypes = [wt.HANDLE]
    kernel.ReleaseMutex.restype = wt.BOOL
    kernel.CloseHandle.argtypes = [wt.HANDLE]
    # Python 3.12 realpath opens the final file with share mode zero on Windows.
    # Resolving that file BEFORE the mutex races with another process's replace.
    # Resolve only its parent; lock naming must never open the mutable document.
    path = Path(path)
    canonical = path.parent.resolve() / path.name
    identity = hashlib.sha256(str(canonical).casefold().encode('utf-8')).hexdigest()
    handle = kernel.CreateMutexW(None, False, 'Local\\COV.RuntimeState.' + identity)
    if not handle:
        raise ct.WinError(ct.get_last_error())
    try:
        waited = kernel.WaitForSingleObject(handle, 5000)
        if waited == 258:
            raise TimeoutError('Runtime state transaction is still busy: ' + str(path))
        if waited not in (0, 0x80):
            raise ct.WinError(ct.get_last_error())
        try:
            yield
        finally:
            if not kernel.ReleaseMutex(handle):
                raise ct.WinError(ct.get_last_error())
    finally:
        kernel.CloseHandle(handle)


_MISSING = object()


def read_json(path, default=_MISSING):
    with file_transaction(path):
        try:
            source = Path(path).read_text(encoding='utf-8-sig')
        except FileNotFoundError:
            if default is not _MISSING:
                return default
            raise
    return json.loads(source)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        with file_transaction(path):
            os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class LeaseOccupied(RuntimeError):
    pass


@contextmanager
def lease(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+b') as stream:
        if stream.tell() == 0:
            stream.write(b'0')
            stream.flush()
        stream.seek(0)
        try:
            if os.name == 'nt':
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise LeaseOccupied(f'Coordinator still owns {path}; keep the lock and use the active coordinator') from error
        yield


def occupied(path):
    try:
        with lease(path):
            return False
    except LeaseOccupied:
        return True


class ControlStore:
    """control.json is the committed command; commands/ records request history."""
    def __init__(self, root):
        self.root = Path(root)
        self.path = self.root / 'control.json'

    def read(self):
        try:
            value = read_json(self.path)
        except FileNotFoundError:
            return {'mode': 'pause', 'sequence': 0, 'interrupt_epoch': 0}
        if not isinstance(value, dict) or value.get('mode') not in MODES:
            raise ValueError('Invalid control state; refusing to continue')
        return value

    def set(self, mode):
        if mode not in MODES:
            raise ValueError('Unknown control command')
        deadline = time.monotonic() + 5
        while True:
            try:
                with lease(self.root / 'control.lock'):
                    previous = self.read()
                    now = time.time()
                    sequence = int(previous.get('sequence', 0)) + 1
                    record = {'schema': 2, 'mode': mode, 'sequence': sequence,
                        'requested_epoch': now, 'request_id': uuid.uuid4().hex,
                        'interrupt_epoch': now if mode == 'interrupt' else previous.get('interrupt_epoch', 0),
                        'interrupt_sequence': sequence if mode == 'interrupt' else previous.get('interrupt_sequence', 0)}
                    atomic_json(self.root / 'commands' / (str(time.time_ns()) + '.json'), record)
                    atomic_json(self.path, record)
                    return record
            except LeaseOccupied:
                if time.monotonic() >= deadline:
                    raise
                time.sleep(.02)
