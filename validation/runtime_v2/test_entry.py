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

import tempfile
from pathlib import Path
import unittest
import os
from unittest.mock import patch
from entry import translate
from install_runtime import stub, activate, rollback
import resume
from test_runtime import FakeData, FakeBackend


class EntryTests(unittest.TestCase):
    def test_legacy_flags_preserve_selection_limit_and_worker_count(self):
        self.assertEqual(translate(['--legacy']), ['menu'])
        self.assertEqual(translate(['--legacy', '--execute', '--limit', '2', '--workers', '1', '--case', 'OLD-093']),
            ['run', '--limit', '2', '--workers', '1', '--case', 'OLD-093'])
        self.assertEqual(translate(['status']), ['status'])

    def test_modern_and_legacy_stubs_parse_with_space_and_unicode_paths(self):
        for legacy in (True, False):
            content = stub(Path('C:/COV 空格/Resume'), 'resume.py', legacy)
            compile(content, '<compatible entry>', 'exec')

    def test_activation_and_rollback_preserve_original_entries_and_results(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            home, fast = root / 'Resume', root / 'FastPause-2.1'
            data = FakeData(root)
            config = {'work_root': str(root / 'work')}
            engine = resume.Engine(config, data, FakeBackend()); engine.bind()
            resume.atomic(engine.root / 'import-summary.json', {'test_fixture': True})
            proof = root / 'proof.txt'; proof.write_text('fixture native evidence')
            resume.atomic(engine.root / 'native-capabilities.json', {'runtime_identity': engine.identity,
                'binaries': data.binaries, 'capabilities': {'ram_pause': {'passed': True,
                    'methods': ['RPBE1PBE', 'UPBE1PBE']}},
                'evidence': [{'path': str(proof), 'sha256': resume.sha(proof)}]})
            config_file = root / 'config.json'; resume.atomic(config_file, config)
            paths = [home / 'cov_resume_menu.py', home / 'README-先读我.md'] + [
                fast / 'validation/runtime_v2' / f for f in ('resume.py', 'config.json', 'README.md', 'FAST_PAUSE.md', 'WORK_PROMPT.md')]
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True); path.write_text('original ' + path.name)
            before = {p: p.read_bytes() for p in paths}
            prepared = {'bundle_identity': 'fixture', 'bundle': str(root / 'bundle'), 'config': str(config_file)}
            with patch.dict(os.environ, {'LOCALAPPDATA': str(root / 'local')}), \
                    patch.object(resume, 'FrozenData', return_value=data), \
                    patch('windows_job.existing_gaussian', return_value=[]):
                result = activate(home, fast, prepared)
                self.assertIn('active.json', (home / 'cov_resume_menu.py').read_text())
                receipt = Path(result['rollback_receipt'])
                (home / 'cov_resume_menu.py').write_text('later edit')
                with self.assertRaisesRegex(ValueError, 'Entry changed'):
                    rollback(home, receipt)
                (home / 'cov_resume_menu.py').write_bytes(stub(home, 'resume.py', True).encode('utf-8'))
                restored = rollback(home, receipt)
                self.assertEqual(restored['restored_files'], len(paths))
                self.assertEqual(before, {p: p.read_bytes() for p in paths})
                self.assertEqual(proof.read_text(), 'fixture native evidence')


if __name__ == '__main__':
    unittest.main(verbosity=2)
