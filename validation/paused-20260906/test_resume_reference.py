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

import contextlib
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
import zipfile

import resume_reference as resume


class RecoveryTests(unittest.TestCase):
    def test_relocation_is_prefix_scoped_and_preserves_non_path_identity(self):
        value = {'path': resume.OLD_JOBS + r'\REF-001\OLD-013\initial\job.chk',
                 'identity': 'a' * 64, 'note': 'input mentioned ' + resume.OLD_JOBS}
        result = resume.relocated(value, [(resume.OLD_JOBS, Path('D:/continued/jobs'))])
        self.assertEqual(Path(result['path']), Path('D:/continued/jobs/REF-001/OLD-013/initial/job.chk'))
        self.assertEqual(result['identity'], value['identity'])
        self.assertEqual(result['note'], value['note'])
        self.assertEqual(value['path'], resume.OLD_JOBS + r'\REF-001\OLD-013\initial\job.chk')

    def test_archive_paths_cannot_escape(self):
        with tempfile.TemporaryDirectory() as temp:
            for name in ('../outside', '/outside', 'C:/outside', r'..\outside'):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    resume.contained(Path(temp), name)

    def archive_fixture(self, root):
        archive = root / 'part.zip'
        payload = b'original immutable FCHK evidence\n'
        with zipfile.ZipFile(archive, 'w') as bundle:
            bundle.writestr('campaign/input.fch', payload)
        member = {'path': 'campaign/input.fch', 'bytes': len(payload),
                  'sha256': hashlib.sha256(payload).hexdigest()}
        manifest = root / 'manifest.json'
        resume.atomic(manifest, {'archives': [{'name': archive.name, 'bytes': archive.stat().st_size,
            'sha256': resume.sha(archive), 'files': [member]}]})
        return SimpleNamespace(manifest=manifest, assets=root, data_root=root / 'data', extract=False)

    def test_verify_does_not_extract_or_launch(self):
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            args = self.archive_fixture(Path(temp))
            resume.verify_assets(args)
            self.assertFalse(args.data_root.exists())

    def test_extract_and_repeat_preserve_exact_bytes(self):
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            args = self.archive_fixture(Path(temp))
            args.extract = True
            resume.verify_assets(args)
            resume.verify_assets(args)
            self.assertEqual((args.data_root / 'campaign/input.fch').read_bytes(), b'original immutable FCHK evidence\n')

    def test_refuse_overwrite_different_evidence(self):
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            args = self.archive_fixture(Path(temp))
            args.extract = True
            target = args.data_root / 'campaign/input.fch'
            target.parent.mkdir(parents=True)
            target.write_bytes(b'different')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                resume.verify_assets(args)
            self.assertEqual(target.read_bytes(), b'different')

    def test_corrupt_archive_is_rejected_before_restore(self):
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            args = self.archive_fixture(Path(temp))
            (args.assets / 'part.zip').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                resume.verify_assets(args)
            self.assertFalse(args.data_root.exists())

    def test_preparation_rebases_only_working_copy(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data, work = root / 'data', root / 'work'
            resume.atomic(data / 'campaign/references/REF-001/reference-candidates.json', {'reference_set_identity': 'x'})
            state_path = data / 'jobs/REF-001/OLD-013/initial/stage.json'
            resume.atomic(state_path, {'status': 'running', 'checkpoint': resume.OLD_JOBS + r'\REF-001\OLD-013\job.chk'})
            before = state_path.read_bytes()
            resume.prepare_runtime(data, work)
            resume.prepare_runtime(data, work)
            self.assertEqual(state_path.read_bytes(), before)
            result = resume.read(work / 'jobs/REF-001/OLD-013/initial/stage.json')
            self.assertEqual(Path(result['checkpoint']), work / 'jobs/REF-001/OLD-013/job.chk')
            self.assertEqual(result['status'], 'running')
            self.assertFalse((work / 'jobs/REF-001/OLD-013/job.chk').exists())


if __name__ == '__main__':
    unittest.main()
