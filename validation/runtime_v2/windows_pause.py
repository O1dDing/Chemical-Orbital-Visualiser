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

"""Win32 thread handles owned by one Job; transactions live in process_control."""
from __future__ import annotations
import ctypes as ct
from ctypes import wintypes as wt
import os
from pathlib import Path
from process_control import PauseTransaction


class THREADENTRY32(ct.Structure):
    _fields_ = [(name, wt.DWORD) for name in ('dwSize', 'cntUsage', 'th32ThreadID', 'th32OwnerProcessID')] + [
        ('tpBasePri', wt.LONG), ('tpDeltaPri', wt.LONG), ('dwFlags', wt.DWORD)]


class OwnedThreads:
    def __init__(self, kernel, job):
        self.k, self.job = kernel, job
        self.infrastructure = {}
        self.membership_events = []
        self.unresolved = {}
        signatures = {
            'OpenProcess': ([wt.DWORD, wt.BOOL, wt.DWORD], wt.HANDLE),
            'OpenThread': ([wt.DWORD, wt.BOOL, wt.DWORD], wt.HANDLE),
            'IsProcessInJob': ([wt.HANDLE, wt.HANDLE, ct.POINTER(wt.BOOL)], wt.BOOL),
            'GetProcessIdOfThread': ([wt.HANDLE], wt.DWORD),
            'GetThreadTimes': ([wt.HANDLE] + [ct.POINTER(wt.FILETIME)] * 4, wt.BOOL),
            'QueryFullProcessImageNameW': ([wt.HANDLE, wt.DWORD, wt.LPWSTR, ct.POINTER(wt.DWORD)], wt.BOOL),
            'SuspendThread': ([wt.HANDLE], wt.DWORD), 'ResumeThread': ([wt.HANDLE], wt.DWORD),
            'WaitForSingleObject': ([wt.HANDLE, wt.DWORD], wt.DWORD),
            'CreateToolhelp32Snapshot': ([wt.DWORD, wt.DWORD], wt.HANDLE),
            'Thread32First': ([wt.HANDLE, ct.POINTER(THREADENTRY32)], wt.BOOL),
            'Thread32Next': ([wt.HANDLE, ct.POINTER(THREADENTRY32)], wt.BOOL),
            'CloseHandle': ([wt.HANDLE], wt.BOOL),
        }
        for name, (args, result) in signatures.items():
            fun = getattr(kernel, name)
            fun.argtypes, fun.restype = args, result

    def pids(self, guard):
        capacity = 16
        while capacity <= 65536:
            guard('query Job process list')
            class IDS(ct.Structure):
                _fields_ = [('assigned', wt.DWORD), ('listed', wt.DWORD), ('pids', ct.c_size_t * capacity)]
            ids = IDS()
            ok = self.k.QueryInformationJobObject(self.job, 3, ct.byref(ids), ct.sizeof(ids), None)
            error = ct.get_last_error()
            guard('read Job process list')
            if ok and ids.listed >= ids.assigned:
                result = set()
                self.infrastructure = {}
                self.unresolved = {}
                for pid in ids.pids[:ids.listed]:
                    guard(f'OpenProcess for Job member pid={pid}')
                    process = self.k.OpenProcess(0x1000 | 0x100000, False, pid)
                    if not process:
                        error = ct.get_last_error()
                        current = IDS()
                        if not self.k.QueryInformationJobObject(self.job, 3, ct.byref(current), ct.sizeof(current), None):
                            raise ct.WinError(ct.get_last_error())
                        if pid not in current.pids[:current.listed]:
                            self.membership_events.append({'pid': int(pid), 'open_error': error,
                                'confirmed_absent_from_fresh_job_list': True})
                            continue
                        # A still-listed process may be starting/exiting. Keep it
                        # in required coverage, but do not suspend unidentified
                        # threads. Stable coverage cannot acknowledge this scan.
                        self.unresolved[int(pid)] = {'operation': 'OpenProcess', 'error': error}
                        result.add(int(pid))
                        continue
                    try:
                        image = ct.create_unicode_buffer(32768)
                        size = wt.DWORD(len(image))
                        guard(f'QueryFullProcessImageNameW for Job member pid={pid}')
                        if not self.k.QueryFullProcessImageNameW(process, 0, image, ct.byref(size)):
                            error = ct.get_last_error()
                            if self.exited(process):
                                continue
                            self.unresolved[int(pid)] = {'operation': 'QueryFullProcessImageNameW', 'error': error}
                            result.add(int(pid))
                            continue
                        # Windows may create this infrastructure for descendants
                        # even when the initial compute process has no console.
                        # Keep the OS console server responsive; it is not a
                        # computational thread or a disk-checkpoint writer.
                        host = Path(os.environ['SystemRoot']) / 'System32/conhost.exe'
                        if Path(image.value) == host:
                            self.infrastructure[int(pid)] = image.value
                        else:
                            result.add(int(pid))
                    finally:
                        self.close(process)
                return result
            if not ok and error != 234:
                raise ct.WinError(error)
            capacity = max(capacity * 2, int(ids.assigned) + 8)
        raise RuntimeError('Job process list too large to pause safely')

    def threads(self, pids, guard):
        guard('CreateToolhelp32Snapshot')
        snapshot = self.k.CreateToolhelp32Snapshot(0x4, 0)
        if snapshot == ct.c_void_p(-1).value:
            raise ct.WinError(ct.get_last_error())
        found = []
        try:
            guard('Thread32First')
            entry = THREADENTRY32()
            entry.dwSize = ct.sizeof(entry)
            ok = self.k.Thread32First(snapshot, ct.byref(entry))
            error = ct.get_last_error()
            while ok:
                guard('Thread32Next')
                if entry.th32OwnerProcessID in pids and entry.th32OwnerProcessID not in self.unresolved:
                    found.append((int(entry.th32OwnerProcessID), int(entry.th32ThreadID)))
                entry.dwSize = ct.sizeof(entry)
                ok = self.k.Thread32Next(snapshot, ct.byref(entry))
                error = ct.get_last_error()
            if error != 18:
                raise ct.WinError(error)
            return found
        finally:
            self.close(snapshot)

    def pin(self, pid, tid, guard):
        guard('OpenThread')
        handle = self.k.OpenThread(0x2 | 0x40 | 0x100000, False, tid)
        if not handle:
            if ct.get_last_error() == 87:
                return None
            raise ct.WinError(ct.get_last_error())
        retained = False
        try:
            guard('GetProcessIdOfThread')
            if self.k.GetProcessIdOfThread(handle) != pid:
                return None
            process = self.k.OpenProcess(0x1000, False, pid)
            if not process:
                if self.exited(handle):
                    return None
                raise ct.WinError(ct.get_last_error())
            try:
                guard('IsProcessInJob')
                member = wt.BOOL()
                if not self.k.IsProcessInJob(process, self.job, ct.byref(member)):
                    raise ct.WinError(ct.get_last_error())
                if not member.value:
                    return None
                image = ct.create_unicode_buffer(32768)
                size = wt.DWORD(len(image))
                guard('QueryFullProcessImageNameW')
                if not self.k.QueryFullProcessImageNameW(process, 0, image, ct.byref(size)):
                    if self.exited(handle):
                        return None
                    raise ct.WinError(ct.get_last_error())
                # Unknown console infrastructure fails closed. The exact Windows
                # System32 host is classified out before thread enumeration.
                if Path(image.value).name.casefold() in ('conhost.exe', 'openconsole.exe'):
                    raise RuntimeError('Unclassified console host in owned computation scope')
            finally:
                self.close(process)
            guard('GetThreadTimes')
            times = [wt.FILETIME() for _ in range(4)]
            if not self.k.GetThreadTimes(handle, *(ct.byref(t) for t in times)):
                if self.exited(handle):
                    return None
                raise ct.WinError(ct.get_last_error())
            created = (times[0].dwHighDateTime << 32) | times[0].dwLowDateTime
            retained = True
            return (pid, tid, created), handle
        finally:
            if not retained:
                self.close(handle)

    def exited(self, handle):
        result = self.k.WaitForSingleObject(handle, 0)
        if result == 0xFFFFFFFF:
            raise ct.WinError(ct.get_last_error())
        return result == 0

    def suspend(self, handle):
        if self.k.SuspendThread(handle) != 0xFFFFFFFF:
            return True
        if self.exited(handle):
            return False
        raise ct.WinError(ct.get_last_error())

    def resume(self, handle):
        return self.k.ResumeThread(handle) != 0xFFFFFFFF

    def close(self, handle):
        self.k.CloseHandle(handle)


class JobPause(PauseTransaction):
    def __init__(self, kernel, job):
        super().__init__(OwnedThreads(kernel, job))
