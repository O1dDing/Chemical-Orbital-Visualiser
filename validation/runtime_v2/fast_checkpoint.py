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

"""Checkpoint preservation and opt segment policy, independent of process launch.

An immutable copy + hash proves saved bytes, not Gaussian restart validity.
All snapshots here MUST be taken after the whole owned writer tree exits.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import uuid
from state_store import atomic_json as write_atomic, read_json


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def persistent_input(text, opt_segments=False):
    # Named local files keep a recoverable per-attempt RWF (not an anonymous
    # scratch pathname). Save prevents normal-end cleanup of restart evidence.
    if re.search(r'(?im)^\s*%(?:rwf|int|d2e|kjob|nosave|errorsave|save)\b', text):
        raise ValueError('Unexpected pre-existing scratch/stop directives')
    prefix = '%RWF=job.rwf\n%Int=job.int\n%D2E=job.d2e\n%Save\n'
    if opt_segments:
        prefix += '%KJob L103 2\n'
    return prefix + text


def rwf_restart_input(cores, memory=24):
    if not 1 <= cores <= 14 or memory < 1:
        raise ValueError('Invalid resources for RWF restart')
    return (f'%RWF=job.rwf\n%Int=job.int\n%D2E=job.d2e\n%Save\n%Chk=job.chk\n'
            f'%Mem={memory}GB\n%NProcShared={cores}\n#p Restart\n\n')


def cold_snapshot(directory, producer, *, writers_exited, reserve_bytes=0, include_scratch=True):
    """Copy only explicit local restart files; never snapshot a suspended writer."""
    directory = Path(directory).resolve()
    if not writers_exited:
        raise ValueError('Cold checkpoint requires zero active writers; RAM pause is not a save')
    names = ('job.chk', 'job.gjf', 'job.log') + (('job.rwf', 'job.int', 'job.d2e') if include_scratch else ())
    files = [directory / n for n in names
             if (directory / n).is_file()]
    if not (directory / 'job.chk').is_file():
        raise ValueError('No checkpoint to preserve')
    if any(p.is_symlink() or not p.resolve().is_relative_to(directory) for p in files):
        raise ValueError('Restart files must be regular files inside their owned attempt')
    required = sum(p.stat().st_size for p in files)
    if shutil.disk_usage(directory).free < required + reserve_bytes:
        raise OSError('Insufficient free space for cold copy; original restart files retained untouched')
    stage = directory / ('snapshot-' + uuid.uuid4().hex + '.partial')
    stage.mkdir()
    rows = {}
    for source in files:
        before = (source.stat().st_size, source.stat().st_mtime_ns, file_hash(source))
        target = stage / source.name
        with source.open('rb') as src, target.open('xb') as dst:
            shutil.copyfileobj(src, dst, 4*1024*1024)
            dst.flush(); os.fsync(dst.fileno())
        after = (source.stat().st_size, source.stat().st_mtime_ns, file_hash(source))
        if before != after or file_hash(target) != before[2]:
            raise ValueError('Source changed or copy checksum mismatch; partial snapshot not accepted')
        rows[source.name] = {'bytes': before[0], 'sha256': before[2]}
    receipt = {'schema': 1, 'files': rows, 'producer': producer, 'created_epoch': time.time(),
               'writers_exited': True, 'restart_validated': False,
               'note': 'Immutable bytes only; formchk/identity and native restart acceptance are separate'}
    write_atomic(stage / 'manifest.json', receipt)
    final = stage.with_suffix('')
    os.replace(stage, final)
    return {'directory': str(final), 'manifest_sha256': file_hash(final / 'manifest.json')}


def restore_snapshot(snapshot, destination, producer, *, require_rwf=False):
    folder = Path(snapshot['directory']).resolve()
    if folder.name.endswith('.partial') or file_hash(folder / 'manifest.json') != snapshot['manifest_sha256']:
        raise ValueError('Incomplete or modified checkpoint snapshot')
    receipt = read_json(folder / 'manifest.json')
    if receipt['producer'] != producer or not receipt.get('writers_exited'):
        raise ValueError('Snapshot producer/runtime/scientific identity mismatch')
    if require_rwf and not receipt['files'].get('job.rwf', {}).get('bytes'):
        raise ValueError('Named nonempty RWF absent; cannot pretend analytic Freq can restart from CHK alone')
    if not receipt['files'].get('job.chk', {}).get('bytes'):
        raise ValueError('Snapshot has no nonempty CHK')
    destination = Path(destination)
    allowed = {'job.chk', 'job.rwf', 'job.int', 'job.d2e', 'job.gjf', 'job.log'}
    for name, row in receipt['files'].items():
        if name not in allowed or (folder / name).is_symlink():
            raise ValueError('Unexpected snapshot member')
        src = folder / name
        if src.stat().st_size != row['bytes'] or file_hash(src) != row['sha256']:
            raise ValueError('Checkpoint snapshot member changed: ' + name)
    bytes_needed = sum(row['bytes'] for n, row in receipt['files'].items() if n not in ('job.log', 'job.gjf'))
    if shutil.disk_usage(destination).free < bytes_needed + 2**30:
        raise OSError('Insufficient space for independent RWF restart copy')
    for name in receipt['files']:
        if name in ('job.log', 'job.gjf'):
            continue  # each new attempt owns its new input and output
        target = destination / name
        if target.exists():
            raise ValueError('Refusing to overwrite an attempt restart file')
        with (folder / name).open('rb') as src, target.open('xb') as dst:
            shutil.copyfileobj(src, dst, 4*1024*1024)
            dst.flush(); os.fsync(dst.fileno())
        if file_hash(target) != receipt['files'][name]['sha256']:
            raise ValueError('Restart copy checksum mismatch')
    return receipt


def cooperative_opt_stop(log, input_text):
    """Strict detector, additionally gated by same-installation native acceptance.

    KJob is predeclared, NOT a hot-edited input or a timing-based kill after a
    log line. The two l103 returns and optimizer evidence distinguish the stop
    from an unrelated error. Unrecognized producer output fails closed.
    """
    if not re.search(r'(?im)^%KJob\s+L103\s+2\s*$', input_text):
        return False
    if re.search(r'Convergence failure|SCF has not converged|Erroneous write|No space left|segmentation|Error termination', log, re.I):
        return False
    returned = re.findall(r'Leave Link\s+(\d+)\s+at', log)
    return (returned.count('103') == 2 and returned[-1:] == ['103'] and
            'Maximum Force' in log and 'Step number' in log)


class ActiveClock:
    def __init__(self, now):
        self.started = now
        self.held_since = None
        self.held_seconds = 0.0

    def set_held(self, held, now):
        if held and self.held_since is None:
            self.held_since = now
        elif not held and self.held_since is not None:
            self.held_seconds += now - self.held_since
            self.held_since = None

    def paused(self, now):
        return self.held_seconds + (now - self.held_since if self.held_since is not None else 0.0)

    def active(self, now):
        return max(0.0, now - self.started - self.paused(now))
