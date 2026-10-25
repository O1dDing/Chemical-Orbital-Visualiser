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

"""Verified, non-destructive import of idle Runtime 2.1 evidence.

Imported receipts retain their original producer and immutable artifact paths.
Only the consumer runtime identity changes. Native capability receipts never move.
The caller must hold the global and work coordinator leases throughout import.
"""
from __future__ import annotations
import copy
from pathlib import Path
import re

from state_store import atomic_json, read_json
from fast_checkpoint import file_hash


def evidence(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': file_hash(path)}


def verify_evidence(rows):
    for row in rows:
        if file_hash(Path(row['path'])) != row['sha256']:
            raise ValueError('Imported source changed: ' + row['path'])


def validate_source(source, engine, case_id, phase):
    imported = source.get('runtime_import')
    if not imported:
        return engine.producer(case_id, phase)
    verify_evidence(imported['evidence'])
    producer = imported['producer']
    expected = dict(engine.producer(case_id, phase), runtime_identity=imported['source_runtime_identity'])
    if producer != expected or imported['consumer_runtime_identity'] != engine.identity:
        raise ValueError('Imported recovery producer/reference/binaries mismatch')
    return producer


def import_runtime(engine, source):
    from policy import next_stage
    source = Path(source).resolve()
    if source == engine.root.resolve() or source.parent != engine.work:
        raise ValueError('Import source must be a separate runtime under the locked work root')
    binding_path = source / 'binding.json'
    binding = read_json(binding_path)
    if (binding.get('version') not in ('2.1-fastpause', '2.2-unified', '3.0-adaptive-physical') or
            binding.get('parent_runner_identity') != engine.data.parent_identity or
            binding.get('reference_set_identity') != engine.data.manifest['reference_set_identity'] or
            binding.get('binaries') != engine.data.binaries):
        raise ValueError('Source runtime science or Gaussian installation does not match')
    origin = binding['runtime_identity']
    binding_proof = evidence(binding_path)
    summary_path = engine.root / 'import-summary.json'
    previous = read_json(summary_path) if summary_path.exists() else None
    if previous:
        if previous.get('source_binding') != binding_proof:
            raise ValueError('Runtime already imported from another source')
        verify_evidence(previous['source_evidence'])
        return previous
    plan, proofs = {}, [binding_proof]
    counts = {'collected_parts': 0, 'recovery_sources': 0, 'review_holds': 0, 'candidate_collected': 0}
    for case_id, reference in engine.data.candidates.items():
        old_case, new_case = source / 'jobs' / case_id, engine.jobs / case_id
        records = {}
        if (new_case / 'result.json').exists() or any(new_case.glob('*/attempt-*')):
            raise ValueError('Import requires a destination without new calculation attempts')
        for old_phase in sorted(old_case.glob('*')):
            if not old_phase.is_dir() or not re.fullmatch(r'(opt|freq|stability|force|sp)-\d+', old_phase.name):
                continue
            stage_path = old_phase / 'stage.json'
            if stage_path.exists():
                row = read_json(stage_path)
                if row.get('runtime_identity') != origin or row.get('reference_identity') != reference['identity']:
                    raise ValueError('Source stage identity mismatch: ' + str(stage_path))
                if row.get('status') == 'collected':
                    stage_proofs = [evidence(stage_path)]
                    for kind in ('checkpoint', 'fchk', 'log'):
                        proof = evidence(row[kind])
                        if proof['sha256'] != row[kind + '_sha256']:
                            raise ValueError('Source stage artifact changed: ' + row[kind])
                        stage_proofs.append(proof)
                    row = copy.deepcopy(row)
                    row.update(runtime_identity=engine.identity, scientific_pass=False,
                               original_producer_runtime_identity=row.get('original_producer_runtime_identity', origin),
                               runtime_import={'source_binding': binding_proof, 'source_stage': stage_proofs[0],
                                               'previous_import': row.get('runtime_import')})
                    plan[new_case / old_phase.name / 'stage.json'] = row
                    proofs.extend(stage_proofs)
                    records[old_phase.name] = row
                    counts['collected_parts'] += 1
                    continue
            recovery = None
            local_proofs = [binding_proof]
            for attempt_path in reversed(sorted(old_phase.glob('attempt-*/attempt.json'))):
                attempt = read_json(attempt_path)
                if attempt.get('status') not in ('running', 'interrupted', 'timeout', 'checkpointed'):
                    continue
                if (attempt.get('runtime_identity') != origin or attempt.get('binaries') != engine.data.binaries or
                        attempt.get('reference_identity') != reference['identity']):
                    raise ValueError('Source attempt identity mismatch')
                input_proof = evidence(attempt_path.parent / 'job.gjf')
                if input_proof['sha256'] != attempt['input_sha256']:
                    raise ValueError('Source attempt input changed')
                checkpoint = attempt_path.parent / 'job.chk'
                if not checkpoint.exists():
                    checkpoint = attempt_path.parent / 'scratch/job.chk'
                recovery = {'checkpoint': str(checkpoint), 'log': str(attempt_path.parent / 'job.log'),
                            'origin': str(attempt_path), 'status': attempt['status'],
                            'snapshot': attempt.get('checkpoint_snapshot')}
                local_proofs.extend([evidence(attempt_path), input_proof])
                break
            legacy_source = old_phase / 'resume-source.json'
            if recovery is None and legacy_source.exists():
                recovery = read_json(legacy_source)
                local_proofs.append(evidence(legacy_source))
            if recovery is not None:
                producer = {'runtime_identity': origin, 'reference_identity': reference['identity'],
                            'binaries': engine.data.binaries, 'case_id': case_id, 'stage': old_phase.name}
                inherited = recovery.get('runtime_import')
                if inherited:
                    verify_evidence(inherited['evidence'])
                    if (inherited['consumer_runtime_identity'] != origin or
                        inherited['producer'] != dict(producer, runtime_identity=inherited['source_runtime_identity'])):
                        raise ValueError('Transitive recovery provenance mismatch')
                    producer = inherited['producer']
                    local_proofs.extend(inherited['evidence'])
                for kind in ('checkpoint', 'log'):
                    if Path(recovery[kind]).is_file():
                        local_proofs.append(evidence(recovery[kind]))
                snapshot = recovery.get('snapshot')
                if snapshot:
                    manifest_path = Path(snapshot['directory']) / 'manifest.json'
                    proof = evidence(manifest_path)
                    receipt = read_json(manifest_path)
                    if proof['sha256'] != snapshot['manifest_sha256'] or receipt['producer'] != producer:
                        raise ValueError('Source cold snapshot producer changed')
                    local_proofs.append(proof)
                    for name, member in receipt['files'].items():
                        if Path(name).name != name:
                            raise ValueError('Invalid snapshot member')
                        item = evidence(manifest_path.parent / name)
                        if item['sha256'] != member['sha256']:
                            raise ValueError('Source cold snapshot changed')
                        local_proofs.append(item)
                recovery['runtime_import'] = {'producer': producer, 'source_runtime_identity': producer['runtime_identity'],
                    'consumer_runtime_identity': engine.identity, 'evidence': local_proofs}
                plan[new_case / old_phase.name / 'resume-source.json'] = recovery
                proofs.extend(local_proofs)
                counts['recovery_sources'] += 1
        review = old_case / 'review.json'
        if review.exists():
            plan[new_case / 'review.json'] = read_json(review)
            proofs.append(evidence(review)); counts['review_holds'] += 1
        if next_stage(reference, records) == 'candidate_collected':
            plan[new_case / 'result.json'] = {'status': 'candidate_collected', 'scientific_pass': False,
                'runtime_identity': engine.identity, 'original_producer_runtime_identity': origin,
                'strategy': 'verified_runtime_import'}
            counts['candidate_collected'] += 1
    # Validate the whole plan before writing anything. A crash during publication
    # has no completion marker and is safely replayable with identical content.
    for path, row in plan.items():
        if path.exists() and read_json(path) != row:
            raise ValueError('Refusing to overwrite destination evidence: ' + str(path))
    verify_evidence(proofs)
    for path, row in plan.items():
        if not path.exists():
            atomic_json(path, row)
    unique = {row['path']: row for row in proofs}
    summary = {'schema': 1, 'strategy': 'verified_runtime_import', 'runtime_identity': engine.identity,
        'source_binding': binding_proof, 'source_runtime_identity': origin, 'source': str(source),
        'source_evidence': list(unique.values()), 'counts': counts, 'native_capabilities_imported': False,
        'originals_modified': False, 'calculation_started': False, 'scientific_pass': False}
    atomic_json(summary_path, summary)
    engine._records_cache.clear()
    return summary
