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

"""Content-addressed independent mathematical references; never cache verdicts.

Each generation is immutable. A truncated/corrupt generation is retained and
recomputed; publication uses a unique temporary pointer and atomic replacement.
Concurrent misses may compute twice, but cannot expose an incomplete NPY file.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import time
import uuid

import numpy as np

SCHEMA = 1
KINDS = frozenset(('overlap', 'sampled_mo', 'full_mo'))


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode('utf-8')


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            value.update(block)
    return value.hexdigest()


def array_identity(value):
    value = np.ascontiguousarray(value)
    if value.dtype.hasobject:
        raise ValueError('Object arrays are not mathematical references')
    return {'shape': list(value.shape), 'dtype': value.dtype.str,
            'sha256': hashlib.sha256(memoryview(value).cast('B')).hexdigest()}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with temp.open('xb') as stream:
        stream.write(canonical(value) + b'\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


class ReferenceCache:
    def __init__(self, root, source_sha256, engine):
        if not re.fullmatch('[0-9a-f]{64}', source_sha256):
            raise ValueError('An exact source FCHK hash is required')
        if not engine or not engine.get('reference_code') or not engine.get('environment'):
            raise ValueError('Independent code and environment identities are required')
        self.root = Path(root).resolve()
        self.common = {'schema': SCHEMA, 'source_fchk_sha256': source_sha256,
                       'engine': engine}
        self.events = []

    def key(self, kind, specification):
        if kind not in KINDS:
            raise ValueError('Only independent overlap/MO values may be cached')
        contract = dict(self.common, kind=kind, specification=specification)
        return hashlib.sha256(canonical(contract)).hexdigest(), contract

    def _read(self, folder, key, contract):
        pointer = json.loads((folder / 'current.json').read_text(encoding='utf-8'))
        generation = pointer['generation']
        if not re.fullmatch('generation-[0-9a-f]{32}', generation):
            raise ValueError('Invalid cache generation')
        directory = folder / generation
        record_path = directory / 'record.json'
        if sha(record_path) != pointer['record_sha256']:
            raise ValueError('Cache metadata checksum differs')
        record = json.loads(record_path.read_text(encoding='utf-8'))
        if record['key'] != key or record['contract'] != contract:
            raise ValueError('Cache identity differs')
        path = directory / 'reference.npy'
        if path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
            raise ValueError('Cache data checksum differs')
        array = np.load(path, mmap_mode='r', allow_pickle=False)
        if array.dtype.hasobject or list(array.shape) != record['shape'] or array.dtype.str != record['dtype']:
            raise ValueError('Cache array layout differs')
        del array
        return path, record

    def materialize(self, kind, specification, writer):
        """writer(path) creates an independent NPY, optionally via a memmap."""
        start = time.perf_counter()
        key, contract = self.key(kind, specification)
        folder = self.root / 'objects' / key[:2] / key
        reason = None
        if (folder / 'current.json').exists():
            try:
                path, record = self._read(folder, key, contract)
                self.events.append({'kind': kind, 'key': key, 'status': 'hit',
                                    'reference_sha256': record['sha256'],
                                    'seconds': time.perf_counter() - start})
                return path
            except (OSError, ValueError, KeyError, EOFError) as error:
                reason = f'{type(error).__name__}: {error}'
        directory = folder / ('generation-' + uuid.uuid4().hex)
        directory.mkdir(parents=True)
        path = directory / 'reference.npy'
        try:
            writer(path)
            # An invalid array or interrupted writer can never become current.
            array = np.load(path, mmap_mode='r', allow_pickle=False)
            if array.dtype.hasobject:
                raise ValueError('Object arrays are forbidden')
            record = {'key': key, 'contract': contract, 'shape': list(array.shape),
                      'dtype': array.dtype.str, 'bytes': path.stat().st_size,
                      'sha256': sha(path)}
            del array
            with path.open('r+b') as stream:
                os.fsync(stream.fileno())
            atomic_json(directory / 'record.json', record)
            atomic_json(folder / 'current.json', {'generation': directory.name,
                        'record_sha256': sha(directory / 'record.json')})
        except Exception as error:
            atomic_json(directory / 'failed.json', {'error': f'{type(error).__name__}: {error}',
                        'published': False, 'prior_corruption': reason})
            raise
        self.events.append({'kind': kind, 'key': key,
                            'status': 'recomputed_corrupt' if reason else 'miss',
                            'prior_corruption': reason,
                            'reference_sha256': record['sha256'],
                            'seconds': time.perf_counter() - start})
        return path

    def array(self, kind, specification, compute):
        def writer(path):
            np.save(path, np.asarray(compute()), allow_pickle=False)
        return np.load(self.materialize(kind, specification, writer),
                       mmap_mode='r', allow_pickle=False)


def from_review_environment(source_sha256, reviewer_path):
    """The bounded batch supervisor hashes the dependencies before/after batch.

    Missing provenance disables caching. The environment document is supplied
    by the supervisor, and its exact hash is passed independently to each child.
    Production EXEs, case names, pass results and thresholds are not cache keys.
    """
    root = os.environ.get('COV_REFERENCE_CACHE_ROOT')
    environment = os.environ.get('COV_REFERENCE_ENV_MANIFEST')
    if not root and not environment:
        return None
    if not root or not environment or sha(environment) != os.environ.get('COV_REFERENCE_ENV_SHA256'):
        raise ValueError('Missing or changed independent reference environment')
    manifest = json.loads(Path(environment).read_text(encoding='utf-8'))
    cache_config = manifest['reference_cache']
    if Path(root).resolve() != Path(cache_config['root']).resolve():
        raise ValueError('Cache root differs from frozen review configuration')
    engine = {'reference_code': {
                'reviewer': sha(reviewer_path), 'cache_module': sha(Path(__file__))},
              'environment': {key: manifest[key] for key in
                ('reference_dependencies', 'python_runtime', 'python_version')},
              'coordinate_convention': 'existing CUDA float interpolation; points_at source pinned by reviewer hash',
              'ao_convention': 'IOData/GBasis source basis; AO transform retained in existing checker',
              'arithmetic': 'numpy float64; OMP/OpenBLAS/MKL/NumExpr each 1 thread'}
    return ReferenceCache(root, source_sha256, engine)
