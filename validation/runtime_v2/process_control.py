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

"""Pause transaction independent of Win32 and persistent command storage."""
from __future__ import annotations
import time


class PauseCancelled(Exception):
    pass


class PauseDeadline(TimeoutError):
    pass


class PauseTransaction:
    def __init__(self, threads, *, clock=time.monotonic, sleep=time.sleep):
        self.threads_api = threads
        self.clock, self.sleep = clock, sleep
        self.handles = {}
        self.held = False
        self.state = 'running'
        self.last_diagnostic = {}

    def hold(self, timeout=5.0, refresh=False, cancel=lambda: False, on_progress=lambda _: None):
        if self.held and not refresh:
            return len(self.handles)
        if timeout <= 0:
            raise ValueError('Pause deadline must be positive')
        started = self.clock()
        deadline = started + timeout
        next_cancel = started
        stable, previous, scans = 0, None, 0
        operation = 'start'

        def guard(where):
            nonlocal operation, next_cancel
            operation = where
            now = self.clock()
            if now >= deadline:
                raise PauseDeadline('Owned thread scan exceeded the pause deadline at ' + where)
            if now >= next_cancel:
                next_cancel = now + .02
                if cancel():
                    raise PauseCancelled()

        self.state = 'pausing'
        try:
            while True:
                guard('job membership')
                pids = self.threads_api.pids(guard)
                if not pids:
                    self.close_after_exit()
                    return 0
                live = set()
                covered = set()
                for pid, tid in self.threads_api.threads(pids, guard):
                    guard('pin owned thread')
                    token = self.threads_api.pin(pid, tid, guard)
                    if token is None:
                        continue
                    key, handle = token
                    retained = False
                    try:
                        if self.threads_api.exited(handle):
                            continue
                        if key not in self.handles:
                            guard('suspend owned thread')
                            if not self.threads_api.suspend(handle):
                                continue
                            # Record success before any callback or deadline check.
                            self.handles[key] = handle
                            retained = True
                        guard('confirm suspension')
                        live.add(key)
                        covered.add(pid)
                    finally:
                        if not retained:
                            self.threads_api.close(handle)
                for key, handle in list(self.handles.items()):
                    guard('retire exited thread')
                    if self.threads_api.exited(handle):
                        self.threads_api.close(handle)
                        del self.handles[key]
                        live.discard(key)
                covered = {key[0] for key in live}
                observed = (frozenset(pids), frozenset(live))
                stable = stable + 1 if observed == previous and covered == pids else 0
                previous = observed
                scans += 1
                self.last_diagnostic = {'operation': operation, 'scans': scans,
                    'elapsed_seconds': self.clock()-started, 'owned_processes': sorted(pids),
                    'covered_processes': sorted(covered), 'suspended_thread_handles': len(self.handles)}
                self.last_diagnostic['unsuspended_infrastructure'] = dict(getattr(self.threads_api, 'infrastructure', {}))
                self.last_diagnostic['membership_events'] = list(getattr(self.threads_api, 'membership_events', []))
                self.last_diagnostic['unresolved_processes'] = dict(getattr(self.threads_api, 'unresolved', {}))
                on_progress(self.last_diagnostic)
                guard('confirm stable membership')
                if stable >= 2 and live:
                    self.held = True
                    self.state = 'ram_paused'
                    return len(self.handles)
                self.sleep(min(.04, max(0., deadline-self.clock())))
        except BaseException as error:
            self.last_diagnostic = dict(self.last_diagnostic, operation=operation,
                elapsed_seconds=self.clock()-started, error_type=type(error).__name__, error=str(error))
            self.resume()
            if isinstance(error, PauseCancelled):
                self.state = 'running'
                return 0
            self.state = 'pause_failed'
            raise

    def resume(self):
        self.state = 'resuming'
        errors = []
        for key, handle in list(self.handles.items()):
            try:
                if not self.threads_api.exited(handle) and not self.threads_api.resume(handle):
                    if not self.threads_api.exited(handle):
                        errors.append(key)
                        continue
                self.threads_api.close(handle)
                del self.handles[key]
            except Exception as error:
                errors.append((key, str(error)))
        self.held = bool(self.handles)
        if errors:
            self.state = 'resume_failed'
            raise RuntimeError('Could not undo owned suspend counts: ' + repr(errors))
        self.state = 'running'

    def close_after_exit(self):
        for handle in self.handles.values():
            self.threads_api.close(handle)
        self.handles.clear()
        self.held = False
        self.state = 'exited'
