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

"""Bounded Windows process trees for the COV/Gaussian validation campaign.

Processes start suspended, join a Job Object, and only then execute. Affinity
and committed-memory limits apply to every descendant (including Gaussian
links). Completion waits for the entire tree, not only the launcher. No GUI
window or shell command construction is required.
"""
from __future__ import annotations

import ctypes as ct
from ctypes import wintypes as wt
import json
import os
from pathlib import Path
import subprocess
import time

GIB = 1024 ** 3


def atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def kernel():
    if os.name != "nt":
        raise RuntimeError("The production process supervisor requires Windows")
    return ct.WinDLL("kernel32", use_last_error=True)


def require(value):
    if not value:
        raise ct.WinError(ct.get_last_error())
    return value


class STARTUPINFO(ct.Structure):
    _fields_ = [("cb", wt.DWORD), ("lpReserved", wt.LPWSTR),
                ("lpDesktop", wt.LPWSTR), ("lpTitle", wt.LPWSTR),
                ("dwX", wt.DWORD), ("dwY", wt.DWORD),
                ("dwXSize", wt.DWORD), ("dwYSize", wt.DWORD),
                ("dwXCountChars", wt.DWORD), ("dwYCountChars", wt.DWORD),
                ("dwFillAttribute", wt.DWORD), ("dwFlags", wt.DWORD),
                ("wShowWindow", wt.WORD), ("cbReserved2", wt.WORD),
                ("lpReserved2", ct.POINTER(ct.c_byte)),
                ("hStdInput", wt.HANDLE), ("hStdOutput", wt.HANDLE),
                ("hStdError", wt.HANDLE)]


class PROCESSINFO(ct.Structure):
    _fields_ = [("hProcess", wt.HANDLE), ("hThread", wt.HANDLE),
                ("dwProcessId", wt.DWORD), ("dwThreadId", wt.DWORD)]


class BASICLIMIT(ct.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ct.c_int64),
                ("PerJobUserTimeLimit", ct.c_int64), ("LimitFlags", wt.DWORD),
                ("MinimumWorkingSetSize", ct.c_size_t),
                ("MaximumWorkingSetSize", ct.c_size_t),
                ("ActiveProcessLimit", wt.DWORD), ("Affinity", ct.c_size_t),
                ("PriorityClass", wt.DWORD), ("SchedulingClass", wt.DWORD)]


class IOCounters(ct.Structure):
    _fields_ = [(name, ct.c_uint64) for name in
                ("ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                 "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class EXTENDEDLIMIT(ct.Structure):
    _fields_ = [("BasicLimitInformation", BASICLIMIT), ("IoInfo", IOCounters),
                ("ProcessMemoryLimit", ct.c_size_t), ("JobMemoryLimit", ct.c_size_t),
                ("PeakProcessMemoryUsed", ct.c_size_t), ("PeakJobMemoryUsed", ct.c_size_t)]


class ACCOUNTING(ct.Structure):
    _fields_ = [(name, ct.c_int64) for name in
                ("TotalUserTime", "TotalKernelTime", "ThisPeriodTotalUserTime",
                 "ThisPeriodTotalKernelTime")] + [(name, wt.DWORD) for name in
                ("TotalPageFaultCount", "TotalProcesses", "ActiveProcesses",
                 "TotalTerminatedProcesses")]


class COREINFO(ct.Structure):
    _fields_ = [("ProcessorMask", ct.c_size_t), ("Relationship", wt.DWORD),
                ("Union", ct.c_uint64 * 2)]


def physical_core_masks(limit=12):
    """One logical processor per physical core; this avoids counting SMT as cores."""
    k = kernel()
    k.GetLogicalProcessorInformation.argtypes = [ct.c_void_p, ct.POINTER(wt.DWORD)]
    size = wt.DWORD()
    k.GetLogicalProcessorInformation(None, ct.byref(size))
    if not size.value:
        raise ct.WinError(ct.get_last_error())
    buf = ct.create_string_buffer(size.value)
    require(k.GetLogicalProcessorInformation(buf, ct.byref(size)))
    cores = []
    for item in (COREINFO * (size.value // ct.sizeof(COREINFO))).from_buffer(buf):
        if item.Relationship == 0:
            # Selecting one sibling gives a hard, conservative thread ceiling.
            cores.append(item.ProcessorMask & -item.ProcessorMask)
    if len(cores) < limit:
        limit = len(cores)
    if not limit:
        raise RuntimeError("No processor cores detected")
    return sorted(cores)[:limit]


def pin_current(mask):
    k = kernel()
    k.GetCurrentProcess.restype = wt.HANDLE
    k.SetProcessAffinityMask.argtypes = [wt.HANDLE, ct.c_size_t]
    require(k.SetProcessAffinityMask(k.GetCurrentProcess(), mask))


def run_tree(args, cwd: Path, env: dict, cpu_mask: int, memory_gib: int,
             timeout_seconds: float, on_started=None, affinity=True):
    """Run an argv list under kernel-enforced limits and record real tree time."""
    if cpu_mask <= 0 or memory_gib <= 0 or timeout_seconds <= 0:
        raise ValueError("Resource and timeout limits must be positive")
    k = kernel()
    k.CreateJobObjectW.argtypes = [ct.c_void_p, wt.LPCWSTR]
    k.CreateJobObjectW.restype = wt.HANDLE
    k.SetInformationJobObject.argtypes = [wt.HANDLE, ct.c_int, ct.c_void_p, wt.DWORD]
    k.QueryInformationJobObject.argtypes = [wt.HANDLE, ct.c_int, ct.c_void_p,
                                          wt.DWORD, ct.c_void_p]
    k.AssignProcessToJobObject.argtypes = [wt.HANDLE, wt.HANDLE]
    k.TerminateJobObject.argtypes = [wt.HANDLE, wt.UINT]
    k.GetExitCodeProcess.argtypes = [wt.HANDLE, ct.POINTER(wt.DWORD)]
    k.CloseHandle.argtypes = [wt.HANDLE]
    k.ResumeThread.argtypes = [wt.HANDLE]
    k.ResumeThread.restype = wt.DWORD
    k.CreateProcessW.argtypes = [wt.LPCWSTR, wt.LPWSTR, ct.c_void_p, ct.c_void_p,
                                wt.BOOL, wt.DWORD, ct.c_void_p, wt.LPCWSTR,
                                ct.POINTER(STARTUPINFO), ct.POINTER(PROCESSINFO)]
    job = require(k.CreateJobObjectW(None, None))
    process = PROCESSINFO()
    import msvcrt
    # Some Windows Fortran runtimes dereference the standard handles during
    # startup even when the application receives explicit input/output paths.
    # Give the hidden process real handles instead of the null GUI defaults.
    stdin_file = open(os.devnull, "rb")
    console_file = (Path(cwd) / "launcher.log").open("ab")
    started = time.time()
    begin = time.monotonic()
    timed_out = False
    try:
        limits = EXTENDEDLIMIT()
        # AFFINITY | JOB_MEMORY | KILL_ON_JOB_CLOSE; children cannot break away.
        limits.BasicLimitInformation.LimitFlags = 0x200 | 0x2000
        if affinity:
            limits.BasicLimitInformation.LimitFlags |= 0x10
            limits.BasicLimitInformation.Affinity = cpu_mask
        limits.JobMemoryLimit = memory_gib * GIB
        require(k.SetInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits)))
        cpu_rate = None
        if not affinity:
            # Gaussian 16W A.03's PGI runtime faults before reading its input
            # under a restricted inherited affinity. A kernel hard CPU-rate
            # limit plus the input/runtime thread count bounds its work without
            # triggering that startup defect. All links inherit the job limit.
            cpu_rate = max(1, min(10000, int(10000 * cpu_mask.bit_count() / os.cpu_count())))
            rate = (wt.DWORD * 2)(0x1 | 0x4, cpu_rate)
            require(k.SetInformationJobObject(job, 15, ct.byref(rate), ct.sizeof(rate)))
        info = STARTUPINFO()
        info.cb = ct.sizeof(info)
        info.dwFlags = 1 | 0x100  # STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES
        info.wShowWindow = 0
        info.hStdInput = msvcrt.get_osfhandle(stdin_file.fileno())
        info.hStdOutput = info.hStdError = msvcrt.get_osfhandle(console_file.fileno())
        os.set_handle_inheritable(info.hStdInput, True)
        os.set_handle_inheritable(info.hStdOutput, True)
        command = ct.create_unicode_buffer(subprocess.list2cmdline([str(x) for x in args]))
        environment = ct.create_unicode_buffer("\0".join(
            f"{key}={value}" for key, value in sorted(env.items(), key=lambda x: x[0].upper())
        ) + "\0\0")
        require(k.CreateProcessW(str(args[0]), command, None, None, True,
                                 0x4 | 0x400 | 0x10, environment, str(cwd),
                                 ct.byref(info), ct.byref(process)))
        try:
            require(k.AssignProcessToJobObject(job, process.hProcess))
        except Exception:
            k.TerminateProcess.argtypes = [wt.HANDLE, wt.UINT]
            k.TerminateProcess(process.hProcess, 125)
            raise
        if on_started:
            on_started({"pid": process.dwProcessId, "started_epoch": started,
                        "cpu_mask": cpu_mask if affinity else None,
                        "cpu_rate_hard_cap": cpu_rate, "memory_limit_gib": memory_gib})
        if k.ResumeThread(process.hThread) == 0xFFFFFFFF:
            raise ct.WinError(ct.get_last_error())
        accounting = ACCOUNTING()
        while True:
            require(k.QueryInformationJobObject(job, 1, ct.byref(accounting),
                                                ct.sizeof(accounting), None))
            if accounting.ActiveProcesses == 0:
                break
            if time.monotonic() - begin > timeout_seconds:
                require(k.TerminateJobObject(job, 124))
                timed_out = True
                # Wait for every link to be gone before releasing its resources.
            time.sleep(0.25 if timed_out else 0.5)
        exit_code = wt.DWORD()
        require(k.GetExitCodeProcess(process.hProcess, ct.byref(exit_code)))
        require(k.QueryInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits), None))
        return {"argv": [str(x) for x in args], "cwd": str(cwd),
                "pid": process.dwProcessId, "started_epoch": started,
                "finished_epoch": time.time(), "wall_seconds": time.monotonic() - begin,
                "exit_code": exit_code.value, "timed_out": timed_out,
                "cpu_mask": cpu_mask if affinity else None,
                "cpu_rate_hard_cap": cpu_rate, "core_count": cpu_mask.bit_count(),
                "memory_limit_gib": memory_gib,
                "peak_tree_commit_bytes": limits.PeakJobMemoryUsed,
                "tree_process_count": accounting.TotalProcesses,
                "tree_cpu_seconds": (accounting.TotalUserTime + accounting.TotalKernelTime) / 1e7}
    finally:
        # A supervisor crash or exception kills unfinished descendants, so they
        # cannot retain a resource allocation invisibly on resume.
        remaining = ACCOUNTING()
        if k.QueryInformationJobObject(job, 1, ct.byref(remaining), ct.sizeof(remaining), None):
            if remaining.ActiveProcesses:
                k.TerminateJobObject(job, 125)
                while (k.QueryInformationJobObject(job, 1, ct.byref(remaining), ct.sizeof(remaining), None)
                       and remaining.ActiveProcesses):
                    time.sleep(0.05)
        k.CloseHandle(job)
        if process.hThread:
            k.CloseHandle(process.hThread)
        if process.hProcess:
            k.CloseHandle(process.hProcess)
        stdin_file.close()
        console_file.close()
