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


def _job_process_ids(k, job):
    """Read actual Job membership; process ancestry is not a substitute."""
    capacity = 32
    while True:
        class PROCESSIDS(ct.Structure):
            _fields_ = [("assigned", wt.DWORD), ("listed", wt.DWORD),
                        ("pids", ct.c_size_t * capacity)]
        members = PROCESSIDS()
        if k.QueryInformationJobObject(job, 3, ct.byref(members), ct.sizeof(members), None):
            if members.listed == members.assigned:
                return [int(pid) for pid in members.pids[:members.listed]]
        elif ct.get_last_error() != 234:  # ERROR_MORE_DATA: membership grew.
            raise ct.WinError(ct.get_last_error())
        capacity = max(capacity * 2, members.assigned)


def _cleanup_declared_helpers(k, job, allowed, phase, previous_deferred):
    """Keep held handles throughout identity checks and any targeted cleanup.

    The caller has observed the root handle signaled. Unknown or unreadable
    members are never terminated here and retain the Job limits. In particular,
    an attached console host may exit naturally only after its helper exits.
    Even after targeted cleanup, the caller waits for the actual Job to be empty.
    """
    held = []
    observations = []
    terminated = []
    try:
        try:
            pids = _job_process_ids(k, job)
        except OSError as error:
            pids = []
            observations.append({"status": "membership_query_failed", "error": str(error)})
        for pid in pids:
            entry = {"pid": pid}
            observations.append(entry)
            handle = k.OpenProcess(0x1000 | 0x100000 | 0x1, False, pid)
            if not handle:
                entry.update(status="process_open_failed", winerror=ct.get_last_error())
                continue
            held.append((handle, entry))
            try:
                wait = k.WaitForSingleObject(handle, 0)
                if wait == 0:
                    entry["status"] = "already_exited"
                    continue
                if wait != 258:
                    raise ct.WinError(ct.get_last_error())
                member = wt.BOOL()
                require(k.IsProcessInJob(handle, job, ct.byref(member)))
                entry["verified_current_job_member"] = bool(member.value)
                if not member.value:
                    entry["status"] = "not_current_job_member"
                    continue
                length = wt.DWORD(32768)
                image = ct.create_unicode_buffer(length.value)
                require(k.QueryFullProcessImageNameW(handle, 0, image, ct.byref(length)))
                created, exited, kernel_time, user_time = (wt.FILETIME() for _ in range(4))
                require(k.GetProcessTimes(handle, ct.byref(created), ct.byref(exited),
                                         ct.byref(kernel_time), ct.byref(user_time)))
                entry.update(image_path=image.value,
                             creation_filetime=(created.dwHighDateTime << 32) | created.dwLowDateTime)
                image_key = os.path.normcase(os.path.normpath(image.value))
                entry["status"] = "declared_helper" if image_key in allowed else "undeclared_image"
            except OSError as error:
                entry.update(status="identity_query_failed", error=str(error))
        deferred = [row for row in observations if row["status"] not in ("declared_helper", "already_exited")]
        signature = None
        if deferred:
            signature = json.dumps(observations, sort_keys=True)
            if signature != previous_deferred:
                phase("finite_driver_cleanup_deferred", members=observations)
        helpers = [(handle, row) for handle, row in held if row["status"] == "declared_helper"]
        if helpers:
            # Evidence records the decision before acting. A failed observer is
            # retained by phase(); only the already verified handles are used.
            phase("finite_driver_cleanup_enter", helpers=[dict(row) for _, row in helpers])
        for handle, row in helpers:
            if k.WaitForSingleObject(handle, 0) == 0:
                continue
            member = wt.BOOL()
            if not k.IsProcessInJob(handle, job, ct.byref(member)) or not member.value:
                phase("finite_driver_helper_cleanup_failed", pid=row["pid"], reason="membership_recheck_failed")
                continue
            if k.TerminateProcess(handle, 125):
                record = {**row, "termination_exit_code": 125}
                terminated.append(record)
                phase("finite_driver_helper_termination_requested", **record)
            else:
                failure = ct.get_last_error()
                if k.WaitForSingleObject(handle, 0) != 0:
                    phase("finite_driver_helper_cleanup_failed", pid=row["pid"], winerror=failure)
        return terminated, signature
    finally:
        for handle, _ in held:
            k.CloseHandle(handle)


def run_tree(args, cwd: Path, env: dict, cpu_mask: int, memory_gib: int,
             timeout_seconds: float, on_started=None, affinity=True, *, stdin_path: Path | None = None,
             on_phase=None, finite_driver_helpers=()):
    """Run an argv list under kernel-enforced limits and record real tree time.

    Lifecycle events use an integer performance-counter origin at entry. The
    legacy started_epoch/wall_seconds fields retain their existing meaning.
    Optional on_phase receives an event snapshot before execution proceeds;
    its elapsed time and any ordinary exception are recorded separately. A
    failed observer does not remove process limits or interrupt cleanup.

    finite_driver_helpers is an explicit opt-in for finite build drivers, using
    existing absolute executable paths declared safe to stop after the finite
    root exits. Only signaled-root, image-verified current Job members may be
    terminated. Other members retain their limits and must finish normally or
    reach the job deadline. Root failures remain failures. The default waits for every descendant and
    must be used for scientific launchers whose children perform the real work.
    """
    if cpu_mask <= 0 or memory_gib <= 0 or timeout_seconds <= 0:
        raise ValueError("Resource and timeout limits must be positive")
    allowed_helpers = set()
    for helper in finite_driver_helpers:
        helper = Path(helper)
        if not helper.is_absolute() or not helper.is_file():
            raise ValueError("Finite-driver helpers require existing absolute executable paths")
        allowed_helpers.add(os.path.normcase(os.path.normpath(str(helper.resolve()))))
    lifecycle_clock = time.get_clock_info("perf_counter")
    if not lifecycle_clock.monotonic:
        raise RuntimeError("Lifecycle evidence requires a monotonic performance counter")
    lifecycle_begin_ns = time.perf_counter_ns()
    lifecycle = []

    def phase(name, **detail):
        elapsed_ns = time.perf_counter_ns() - lifecycle_begin_ns
        event = {"sequence": len(lifecycle), "phase": name,
                 "elapsed_nanoseconds": elapsed_ns, "elapsed_seconds": elapsed_ns / 1e9,
                 "epoch": time.time(), "observer_called": on_phase is not None, **detail}
        lifecycle.append(event)
        observer_begin_ns = time.perf_counter_ns()
        if on_phase is not None:
            try:
                on_phase(dict(event))
            except Exception as error:
                event["observer_error"] = f"{type(error).__name__}: {error}"
            observer_ns = time.perf_counter_ns() - observer_begin_ns
        else:
            # This is an absent callback, not a below-resolution measurement.
            observer_ns = 0
        event["observer_wall_nanoseconds"] = observer_ns
        event["observer_wall_seconds"] = observer_ns / 1e9
        after_ns = time.perf_counter_ns() - lifecycle_begin_ns
        event["after_observer_elapsed_nanoseconds"] = after_ns
        event["after_observer_elapsed_seconds"] = after_ns / 1e9

    phase("supervisor_enter")
    k = kernel()
    k.CreateJobObjectW.argtypes = [ct.c_void_p, wt.LPCWSTR]
    k.CreateJobObjectW.restype = wt.HANDLE
    k.SetInformationJobObject.argtypes = [wt.HANDLE, ct.c_int, ct.c_void_p, wt.DWORD]
    k.QueryInformationJobObject.argtypes = [wt.HANDLE, ct.c_int, ct.c_void_p,
                                          wt.DWORD, ct.c_void_p]
    k.AssignProcessToJobObject.argtypes = [wt.HANDLE, wt.HANDLE]
    k.TerminateJobObject.argtypes = [wt.HANDLE, wt.UINT]
    k.GetExitCodeProcess.argtypes = [wt.HANDLE, ct.POINTER(wt.DWORD)]
    k.WaitForSingleObject.argtypes = [wt.HANDLE, wt.DWORD]
    k.WaitForSingleObject.restype = wt.DWORD
    k.OpenProcess.argtypes = [wt.DWORD, wt.BOOL, wt.DWORD]
    k.OpenProcess.restype = wt.HANDLE
    k.IsProcessInJob.argtypes = [wt.HANDLE, wt.HANDLE, ct.POINTER(wt.BOOL)]
    k.QueryFullProcessImageNameW.argtypes = [wt.HANDLE, wt.DWORD, wt.LPWSTR, ct.POINTER(wt.DWORD)]
    k.GetProcessTimes.argtypes = [wt.HANDLE] + [ct.POINTER(wt.FILETIME)] * 4
    k.TerminateProcess.argtypes = [wt.HANDLE, wt.UINT]
    k.CloseHandle.argtypes = [wt.HANDLE]
    k.ResumeThread.argtypes = [wt.HANDLE]
    k.ResumeThread.restype = wt.DWORD
    k.CreateProcessW.argtypes = [wt.LPCWSTR, wt.LPWSTR, ct.c_void_p, ct.c_void_p,
                                wt.BOOL, wt.DWORD, ct.c_void_p, wt.LPCWSTR,
                                ct.POINTER(STARTUPINFO), ct.POINTER(PROCESSINFO)]
    phase("job_creation_enter")
    job = require(k.CreateJobObjectW(None, None))
    phase("job_created")
    process = PROCESSINFO()
    import msvcrt
    # Some Windows Fortran runtimes dereference the standard handles during
    # startup even when the application receives explicit input/output paths.
    # Give the hidden process real handles instead of the null GUI defaults.
    phase("standard_stream_open_enter")
    stdin_file = open(os.devnull if stdin_path is None else stdin_path, "rb")
    console_file = (Path(cwd) / "launcher.log").open("ab")
    phase("standard_streams_opened")
    started = time.time()
    begin = time.monotonic()
    timed_out = False
    root_exit_observed = None
    cleaned_helpers = []
    deferred_helpers = None
    try:
        phase("resource_limits_enter")
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
        phase("resource_limits_installed")
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
        phase("process_creation_enter")
        require(k.CreateProcessW(str(args[0]), command, None, None, True,
                                 0x4 | 0x400 | 0x10, environment, str(cwd),
                                 ct.byref(info), ct.byref(process)))
        phase("process_created_suspended", pid=process.dwProcessId)
        try:
            phase("job_assignment_enter")
            require(k.AssignProcessToJobObject(job, process.hProcess))
            phase("job_assigned")
        except Exception:
            k.TerminateProcess.argtypes = [wt.HANDLE, wt.UINT]
            k.TerminateProcess(process.hProcess, 125)
            raise
        if on_started:
            phase("started_callback_enter")
            on_started({"pid": process.dwProcessId, "started_epoch": started,
                        "cpu_mask": cpu_mask if affinity else None,
                        "cpu_rate_hard_cap": cpu_rate, "memory_limit_gib": memory_gib})
            phase("started_callback_returned")
        else:
            phase("started_callback_absent")
        phase("resume_thread_enter")
        if k.ResumeThread(process.hThread) == 0xFFFFFFFF:
            raise ct.WinError(ct.get_last_error())
        phase("thread_resumed")
        accounting = ACCOUNTING()
        phase("tree_monitor_enter")
        while True:
            require(k.QueryInformationJobObject(job, 1, ct.byref(accounting),
                                                ct.sizeof(accounting), None))
            if root_exit_observed is None:
                root_wait = k.WaitForSingleObject(process.hProcess, 0)
                if root_wait == 0:
                    root_code = wt.DWORD()
                    require(k.GetExitCodeProcess(process.hProcess, ct.byref(root_code)))
                    root_exit_observed = {"exit_code": root_code.value,
                                          "active_job_processes": accounting.ActiveProcesses}
                    phase("root_exit_observed", **root_exit_observed)
                elif root_wait != 258:
                    raise ct.WinError(ct.get_last_error())
            if accounting.ActiveProcesses == 0:
                phase("tree_empty")
                break
            if time.monotonic() - begin > timeout_seconds:
                if not timed_out:
                    phase("execution_deadline_observed", active_processes=accounting.ActiveProcesses)
                require(k.TerminateJobObject(job, 124))
                timed_out = True
                # Wait for every link to be gone before releasing its resources.
            elif allowed_helpers and root_exit_observed is not None:
                cleaned, deferred_helpers = _cleanup_declared_helpers(
                    k, job, allowed_helpers, phase, deferred_helpers)
                cleaned_helpers.extend(cleaned)
            time.sleep(0.25 if timed_out else 0.5)
        exit_code = wt.DWORD()
        require(k.GetExitCodeProcess(process.hProcess, ct.byref(exit_code)))
        require(k.QueryInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits), None))
        phase("final_accounting_collected", exit_code=exit_code.value)
        return {"argv": [str(x) for x in args], "cwd": str(cwd),
                "stdin_path": str(Path(stdin_path).resolve()) if stdin_path is not None else None,
                "console_log": str((Path(cwd) / "launcher.log").resolve()),
                "pid": process.dwProcessId, "started_epoch": started,
                "finished_epoch": time.time(), "wall_seconds": time.monotonic() - begin,
                "exit_code": exit_code.value, "timed_out": timed_out,
                "completion_policy": "finite_driver_declared_helpers" if allowed_helpers else "entire_tree",
                "finite_driver_helpers": sorted(allowed_helpers),
                "root_exit_observation": root_exit_observed,
                "terminated_declared_helpers": cleaned_helpers,
                "cpu_mask": cpu_mask if affinity else None,
                "cpu_rate_hard_cap": cpu_rate, "core_count": cpu_mask.bit_count(),
                "memory_limit_gib": memory_gib,
                "peak_tree_commit_bytes": limits.PeakJobMemoryUsed,
                "tree_process_count": accounting.TotalProcesses,
                "tree_cpu_seconds": (accounting.TotalUserTime + accounting.TotalKernelTime) / 1e7,
                "lifecycle_schema_version": 2,
                "lifecycle_clock": "perf_counter_ns from run_tree entry; integer ticks and observer time retained",
                "lifecycle_clock_origin_nanoseconds": lifecycle_begin_ns,
                "lifecycle_clock_info": {"name": "perf_counter", "implementation": lifecycle_clock.implementation,
                    "resolution_seconds": lifecycle_clock.resolution, "monotonic": lifecycle_clock.monotonic,
                    "adjustable": lifecycle_clock.adjustable,
                    "resolution_is_not_an_accuracy_guarantee": True},
                "legacy_wall_seconds_clock_info": vars(time.get_clock_info("monotonic")),
                "epoch_clock_info": vars(time.get_clock_info("time")),
                "lifecycle_events": lifecycle}
    finally:
        phase("cleanup_enter")
        # A supervisor crash or exception kills unfinished descendants, so they
        # cannot retain a resource allocation invisibly on resume.
        remaining = ACCOUNTING()
        if k.QueryInformationJobObject(job, 1, ct.byref(remaining), ct.sizeof(remaining), None):
            if remaining.ActiveProcesses:
                phase("unfinished_tree_cleanup", active_processes=remaining.ActiveProcesses)
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
        phase("cleanup_complete")
