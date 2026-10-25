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

"""Migration preserves original bytes and cannot import native capabilities."""
import copy
from pathlib import Path
import tempfile
import unittest

import resume as r
from migration import import_runtime, validate_source
from test_runtime import FakeData, FakeBackend, INPUT, FIELDS
from policy import next_stage
from fast_checkpoint import cold_snapshot, restore_snapshot


class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.data = FakeData(self.path, count=18)
        self.engine = r.Engine({'work_root': str(self.path / 'work')}, self.data, FakeBackend())
        self.engine.bind()
        self.old = self.engine.work / 'runtime-v2-fastpause'
        r.atomic(self.old / 'binding.json', {'version': '2.1-fastpause', 'runtime_identity': 'old',
            'parent_runner_identity': self.data.parent_identity, 'binaries': self.data.binaries,
            'reference_set_identity': self.data.manifest['reference_set_identity']})

    def stage(self):
        folder = self.old / 'jobs/OLD-001/opt-00'
        folder.mkdir(parents=True)
        row = {'status': 'collected', 'runtime_identity': 'old', 'name': 'opt-00',
            'reference_identity': self.data.candidates['OLD-001']['identity'], 'fields': FIELDS,
            'diagnostics': {'optimization_completed': True}, 'scientific_pass': False}
        for kind in ('checkpoint', 'fchk', 'log'):
            path = folder / kind
            path.write_text(kind)
            row[kind], row[kind + '_sha256'] = str(path), r.sha(path)
        r.atomic(folder / 'stage.json', row)
        return folder

    def test_import_preserves_results_review_and_original_bytes_idempotently(self):
        self.stage()
        r.atomic(self.old / 'jobs/OLD-018/review.json', {'status': 'needs_review', 'reason': 'timeout'})
        r.atomic(self.old / 'native-capabilities.json', {'passed': True})
        before = {str(p): r.sha(p) for p in self.old.rglob('*') if p.is_file()}
        result = import_runtime(self.engine, self.old)
        self.assertEqual(result['counts']['collected_parts'], 1)
        self.assertEqual(result['counts']['review_holds'], 1)
        self.assertEqual(r.read(self.engine.jobs / 'OLD-018/review.json')['reason'], 'timeout')
        row = self.engine.records('OLD-001')['opt-00']
        self.assertEqual(row['original_producer_runtime_identity'], 'old')
        self.assertFalse((self.engine.root / 'native-capabilities.json').exists())
        self.assertEqual(import_runtime(self.engine, self.old), result)
        self.assertEqual(before, {str(p): r.sha(p) for p in self.old.rglob('*') if p.is_file()})

    def test_tampering_prevents_all_destination_stage_writes(self):
        folder = self.stage()
        (folder / 'checkpoint').write_text('changed')
        with self.assertRaises(ValueError):
            import_runtime(self.engine, self.old)
        self.assertFalse(list(self.engine.jobs.glob('*/opt-00/stage.json')))
        self.assertFalse((self.engine.root / 'import-summary.json').exists())

    def test_different_binary_is_not_migratable(self):
        row = r.read(self.old / 'binding.json'); row['binaries'] = {'g16.exe': 'changed'}
        r.atomic(self.old / 'binding.json', row)
        with self.assertRaises(ValueError):
            import_runtime(self.engine, self.old)

    def test_interrupted_snapshot_keeps_old_producer_and_rejects_later_change(self):
        folder = self.old / 'jobs/OLD-001/opt-00/attempt-0001'
        folder.mkdir(parents=True)
        (folder / 'job.gjf').write_text(INPUT)
        (folder / 'job.chk').write_text('checkpoint')
        (folder / 'job.log').write_text('interrupted')
        phase = next_stage(self.data.candidates['OLD-001'], {})
        producer = dict(self.engine.producer('OLD-001', phase), runtime_identity='old')
        snapshot = cold_snapshot(folder, producer, writers_exited=True)
        r.atomic(folder / 'attempt.json', {'status': 'interrupted', 'runtime_identity': 'old',
            'reference_identity': producer['reference_identity'], 'binaries': self.data.binaries,
            'input_sha256': r.sha(folder / 'job.gjf'), 'checkpoint_snapshot': snapshot})
        import_runtime(self.engine, self.old)
        source = self.engine.interrupted_source('OLD-001', phase)
        self.assertEqual(validate_source(source, self.engine, 'OLD-001', phase), producer)
        target = self.path / 'independent'; target.mkdir()
        restore_snapshot(source['snapshot'], target, producer)
        self.assertEqual((target / 'job.chk').read_text(), 'checkpoint')
        (folder / 'job.gjf').write_text('changed')
        with self.assertRaises(ValueError):
            validate_source(source, self.engine, 'OLD-001', phase)


if __name__ == '__main__':
    unittest.main(verbosity=2)
