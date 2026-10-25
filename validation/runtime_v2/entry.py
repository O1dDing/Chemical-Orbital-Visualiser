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

"""Shared launcher, including the original Resume helper command-line interface."""
from pathlib import Path
import argparse
import sys
import resume


def translate(argv):
    if '--legacy' not in argv:
        return argv or ['menu']
    argv = [arg for arg in argv if arg != '--legacy']
    parser = argparse.ArgumentParser(description='COV unified runtime (legacy Resume arguments supported)')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--limit', type=int, default=1)
    parser.add_argument('--workers', type=int, choices=(1, 2), default=1)
    parser.add_argument('--case', action='append', default=[])
    args = parser.parse_args(argv)
    if not args.execute:
        return ['menu']
    result = ['run', '--limit', str(args.limit), '--workers', str(args.workers)]
    for case in args.case:
        result.extend(['--case', case])
    return result


if __name__ == '__main__':
    # The compatibility stub supplies the shared configuration explicitly.
    config, raw = sys.argv[1], sys.argv[2:]
    sys.argv = [str(Path(resume.__file__)), '--config', config] + translate(raw)
    resume.main()
