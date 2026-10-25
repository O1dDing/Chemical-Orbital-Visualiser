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

"""Explicit isolated native acceptance and bounded seven-core OLD-108 smoke.

Never calls Engine.run, never submits the formal reference queue. Imports retain
original artifact identities; all real computations stay under a native-test root.
"""
from contextlib import ExitStack
from pathlib import Path
import argparse
import json
import os
import re
import sys
import time

import resume
import windows_job
from fast_checkpoint import persistent_input, cold_snapshot
from migration import import_runtime
from native_acceptance import run_acceptance
from policy import initial_part, route_signature, preferred_cores
from state_store import lease, atomic_json


def run(config_path, source, output, part):
    config = resume.read(config_path)
    data = resume.FrozenData(config)
    engine = resume.Engine(config, data, windows_job)
    host = windows_job.host_info()
    output.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        stack.enter_context(lease(Path(os.environ['LOCALAPPDATA']) / 'COV/gaussian-runtime-v2.lock'))
        stack.enter_context(lease(data.work / 'jobs/supervisor.lock'))
        if windows_job.existing_gaussian(data.gaussian):
            raise RuntimeError('Existing Gaussian trees block isolated smoke')
        windows_job.full_startup_affinity()
        engine.bind()
        if part == 'import':
            report = import_runtime(engine, source)
            engine.set_mode('pause')
            atomic_json(output / 'migration.json', report)
            print(json.dumps(report['counts']), flush=True)
            return
        if part == 'pause':
            print('Starting isolated water/NO native pause acceptance at adaptive tiers',flush=True)
            run_acceptance(engine, host, capability='ram_pause')
            for method in ('RPBE1PBE','UPBE1PBE'):
                engine.require_capability('ram_pause', method)
            atomic_json(output / 'native-pause.json',resume.read(engine.root / 'native-capabilities.json'))
            print('Native pause capability accepted for this exact runtime identity',flush=True)
            return
        if host['physical_cores'] != 16:
            raise RuntimeError('The seven-core host smoke requires 16 physical cores')
        case='OLD-108'; cores=preferred_cores(data.sizes[case],physical=16)
        assert cores==7 and data.sizes[case]==574
        directory=output / ('heavy-seven-core-' + str(time.time_ns()))
        directory.mkdir();(directory/'scratch').mkdir()
        frozen=data.references/'reference-inputs'/case/'initial.gjf'
        frozen_sha=resume.sha(frozen)
        original=resume.text(frozen)
        inp=persistent_input(initial_part(original,'opt',cores))
        assert route_signature(inp)==route_signature(persistent_input(initial_part(original,'opt',4)))
        (directory/'job.gjf').write_text(inp,encoding='ascii',newline='\n')
        started=[];samples=[];begin=time.monotonic()
        def start(value):
            started.append(value);atomic_json(directory/'started.json',value)
        def tick(value):
            samples.append(value)
            atomic_json(directory/'samples.json',samples)
        def stop():
            return time.monotonic()-begin>=60 and bool(re.search(r'NBasis\s*=\s*574',resume.text(directory/'job.log')))
        reservation=windows_job.new_reservation(cores,label='isolated-smoke/OLD-108',
            metadata={'runtime_identity':engine.identity,'reference_identity':data.candidates[case]['identity'],
                      'frozen_initial_sha256':frozen_sha,'actual_input_sha256':resume.sha(directory/'job.gjf')})
        try:
            with reservation.scope():
                result=windows_job.run_tree([data.gaussian/'g16.exe',directory/'job.gjf',directory/'job.log'],
                    directory,cores,32,120,cancel=stop,on_started=start,on_tick=tick)
        finally:reservation.release()
        logs=resume.text(directory/'job.log')
        receipt={'runtime_identity':engine.identity,'host':host,'case_id':case,'cores':cores,
            'nbasis':data.sizes[case],'frozen_sha256':frozen_sha,'input_sha256':resume.sha(directory/'job.gjf'),
            'directory':str(directory),'result':result,'started':started,'samples':samples,
            'normal_termination':'Normal termination of Gaussian' in logs,
            'gaussian_thread_banner':re.findall(r'.*(?:processors via shared memory|%nprocshared).*',logs,re.I),
            'scientific_pass':False,'formal_case_modified':False,
            'scope':'Bounded resource smoke; intentional stop, not completed optimization or a speedup benchmark'}
        receipt['passed']=bool(started and started[0]['cores']==7 and
            all(started[0]['child_environment'][k]=='7' for k in ('OMP_NUM_THREADS','NCPUS','OMP_THREAD_LIMIT','MKL_NUM_THREADS')) and
            result['active_processes_at_return']==0 and re.search(r'NBasis\s*=\s*574',logs) and
            re.search(r'(?:use up to\s+7\s+processors|%nprocshared=7)',logs,re.I) and
            resume.sha(frozen)==frozen_sha and 'Error termination' not in logs)
        receipt['checkpoint_snapshot']=cold_snapshot(directory,engine.producer(case,type('Phase',(),{'name':'opt-00'})()),
            writers_exited=True,reserve_bytes=2**30,include_scratch=True)
        atomic_json(directory/'smoke.json',receipt)
        atomic_json(output/'heavy-smoke.json',receipt)
        print(json.dumps({k:receipt[k] for k in ('passed','directory','cores','nbasis','normal_termination','scope')},indent=2),flush=True)
        if not receipt['passed']:raise RuntimeError('Heavy resource smoke did not satisfy its assertions')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('part',choices=('import','pause','heavy'))
    args=parser.parse_args();run(args.config,args.source,args.output,args.part)
