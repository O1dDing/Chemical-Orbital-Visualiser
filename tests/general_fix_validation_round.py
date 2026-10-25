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

"""Frozen native collection with bounded process trees and resumable case attempts.

Uses the existing real-control plan/evidence collector. All cases reach terminal
collection states before a separate scientific review is permitted. A collected
case is never labelled a scientific pass here. No Gaussian calculation is started.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import os
import queue
import shutil
import sys
import threading
import time

from validation_process import atomic_json, physical_core_masks, run_tree

LOCK=threading.RLock()
LOCAL=threading.local()

def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): value.update(block)
    return value.hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def checked_identity(manifest,root):
    for name,sha in manifest['programs'].items():
        if digest(root/'programs'/name)!=sha: raise RuntimeError('Frozen program differs: '+name)
    for name,sha in manifest['runner_sources'].items():
        if digest(Path(__file__).parent/name)!=sha: raise RuntimeError('Frozen runner differs: '+name)
    for case in manifest['cases']:
        if digest(case['input'])!=case['sha256']: raise RuntimeError('Frozen input differs: '+case['case_id'])
    for name,sha in manifest['input_artifacts'].items():
        if digest(root/name)!=sha: raise RuntimeError('Frozen input or sidecar differs: '+name)
    for name,sha in manifest['image_dependency_artifacts'].items():
        if digest(root/'image-deps'/name)!=sha: raise RuntimeError('Frozen image dependency differs: '+name)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--round',type=Path,required=True)
    parser.add_argument('--reference-deps',type=Path,required=True)
    args=parser.parse_args()
    root=args.round.resolve()
    manifest=read(root/'manifest.json')
    if manifest['case_count']!=len(manifest['cases']): raise RuntimeError('Incomplete frozen case manifest')
    sys.path.insert(0,str(args.reference_deps))
    sys.path.insert(0,str(root/'image-deps'))
    os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
                      NUMEXPR_NUM_THREADS='1')
    import native_collection_batch as native
    checked_identity(manifest,root)
    workers=manifest['resources']['workers']
    threads=manifest['resources']['threads_per_worker']
    quota=manifest['resources']
    offset=quota.get('core_slot_offset',0)
    if (not 1<=threads or not 1<=workers<=4 or not 0<=offset
            or not workers*threads<=quota['aggregate_cpu_cores']<=12-offset
            or not 0<quota['memory_gib_per_worker']*workers<=quota['aggregate_memory_gib']<=128):
        raise ValueError('Frozen worker quotas exceed the authorized global ceiling')
    masks=physical_core_masks(offset+workers*threads)[offset:]
    slots=queue.Queue()
    for i in range(workers): slots.put(sum(masks[i*threads:(i+1)*threads]))
    results={}
    running={}
    begin=time.perf_counter()
    started=time.time()
    # A local lock protects resumptions from launching concurrent supervisors
    # against the same round. A process exit releases this OS byte-range lock.
    import msvcrt
    lock_file=(root/'coordinator.lock').open('a+b')
    lock_file.seek(0)
    if lock_file.read(1)==b'': lock_file.write(b'0');lock_file.flush()
    lock_file.seek(0)
    msvcrt.locking(lock_file.fileno(),msvcrt.LK_NBLCK,1)

    def progress(*unused):
        with LOCK:
            rows=list(results.values())
            atomic_json(root/'progress.json',{'status':'collecting','started_epoch':started,
                'updated_epoch':time.time(),'expected_cases':manifest['case_count'],
                'terminal_cases':len(rows),'running':dict(running),
                'collection_complete':sum(x['collection_status']=='complete' for x in rows),
                'collection_gaps':sum(x['collection_status']=='completed_with_gaps' for x in rows),
                'collection_errors':sum(x['collection_status']=='error' for x in rows),
                'scientific_analysis':'pending full collection barrier','formal_case_passes':0})

    def run_process(options,unused_manifest,case_id,exe,arguments,destination,stdout_name,timeout):
        expected=manifest['programs'][Path(exe).name]
        if digest(exe)!=expected: raise RuntimeError('Program identity changed before launch')
        directory=destination/('process-'+Path(stdout_name).stem)
        directory.mkdir()
        env=dict(os.environ,COV_CPU_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',
                 OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        def on_started(record):
            atomic_json(directory/'started.json',record)
            with LOCK:
                running[case_id]={'stage':Path(exe).name,'pid':record['pid'],
                                 'directory':str(directory),'started_epoch':time.time()}
            progress()
        outcome=run_tree([exe,*arguments],directory,env,LOCAL.mask,
                         manifest['resources']['memory_gib_per_worker'],timeout,on_started=on_started)
        atomic_json(directory/'process.json',outcome)
        shutil.copy2(directory/'launcher.log',destination/stdout_name)
        if outcome['timed_out']: raise RuntimeError('Process-tree timeout; partial evidence retained')
        if digest(exe)!=expected: raise RuntimeError('Program identity changed during execution')
        return outcome

    native.run_process=run_process
    native.progress=progress
    native.STARTED=time.perf_counter()
    original_plan=native.make_plan

    def make_plan(production):
        plan,expected=original_plan(production)
        # Additional complete textures at two fixed resolutions for each
        # supplied spin's occupied/virtual frontier. Selection still traverses
        # the real browser controls; identities are derived from the input.
        orbitals=production['orbitals']
        frontiers=set()
        for spin in sorted({mo['spin'] for mo in orbitals}):
            occupied=[i for i,mo in enumerate(orbitals) if mo['spin']==spin and mo['occupation']>0]
            virtual=[i for i,mo in enumerate(orbitals) if mo['spin']==spin and mo['occupation']==0]
            if occupied: frontiers.add(max(occupied,key=lambda i:orbitals[i]['energy_hartree']))
            if virtual: frontiers.add(min(virtual,key=lambda i:orbitals[i]['energy_hartree']))
        extra=[]
        for resolution in (64,128):
            extra += [f'scene 0.92 0.0 0.0 2.2 0.03 {resolution}',
                      'seek "panel.browser"','click "browser.filter"','key "Home"','key "Down"','key "Enter"']
            for index in sorted(frontiers):
                extra += [native.command('text','browser.search',f'{index+1} '),
                          native.command('click',f'browser.mo.{index}'),
                          native.command('volume_full',f'frontier-{index+1:04}-{resolution}',index)]
        expected['full_grid_frontiers_zero_based']=sorted(frontiers)
        expected['full_grid_resolutions']=[64,128]
        return plan+'\n'.join(extra)+'\n',expected

    native.make_plan=make_plan

    def collect(case):
        case_id=case['case_id']
        slot=slots.get()
        LOCAL.mask=slot
        attempt_begin=time.perf_counter()
        case_dir=root/'cases'/case_id
        case_dir.mkdir(parents=True,exist_ok=True)
        existing=sorted(case_dir.glob('attempt-*'))
        attempt=case_dir/f'attempt-{len(existing)+1:03d}'
        attempt.mkdir()
        options=SimpleNamespace(output=attempt,build=root/'programs',threads=threads,
                                timeout=manifest['native_timeout_seconds'])
        try:
            if digest(case['input'])!=case['sha256']: raise RuntimeError('Input identity changed')
            for log in case['log_candidates']:
                if digest(log)!=manifest['input_artifacts'][str(Path(log).relative_to(root))]:
                    raise RuntimeError('Producer sidecar changed before collection')
            info=native.collect_one(options,manifest,case)
            evidence=attempt/'cases'/case_id
            # Ordinary and validation instrumentation call the same production
            # science. Collect both now; compare only after the full barrier.
            if (evidence/'production.json').exists():
                ordinary=run_process(options,manifest,case_id,root/'programs'/'cov_scientific_audit_dump_off.exe',
                    [Path(case['input']),'--interactions'],evidence,'production-off.json',900)
                info['ordinary_process']=ordinary
                if ordinary['exit_code']!=0:
                    info['collection_status']='completed_with_gaps'
                    info.setdefault('integrity_errors',[]).append('ordinary production data collection failed')
                elif (evidence/'production-off.json').exists():
                    read(evidence/'production-off.json')
            info['evidence_directory']=str(evidence.relative_to(root))
            info['source_sha256']=digest(case['input'])
            info['round_identity']=manifest['round_identity']
            if info['source_sha256']!=case['sha256']:
                info['collection_status']='error';info['error']='Input changed during collection'
            for log in case['log_candidates']:
                if digest(log)!=manifest['input_artifacts'][str(Path(log).relative_to(root))]:
                    info['collection_status']='error';info['error']='Producer sidecar changed during collection'
            # Capture exact raw evidence identities only after the process tree
            # and lossless screenshot conversion have finished.
            files={str(p.relative_to(evidence)):{'bytes':p.stat().st_size,'sha256':digest(p)}
                   for p in sorted(evidence.rglob('*')) if p.is_file()}
            atomic_json(attempt/'evidence-manifest.json',files)
            info['evidence_manifest']=str((attempt/'evidence-manifest.json').relative_to(root))
            info['evidence_manifest_sha256']=digest(attempt/'evidence-manifest.json')
        except Exception as error:
            info={'case_id':case_id,'collection_status':'error','error':str(error),
                  'round_identity':manifest['round_identity'],'evidence_directory':str(attempt.relative_to(root))}
        finally:
            slots.put(slot)
        info['wall_seconds_including_ordinary_and_packaging']=time.perf_counter()-attempt_begin
        atomic_json(case_dir/'terminal.json',info)
        with LOCK:
            results[case_id]=info;running.pop(case_id,None)
        progress()
        print(json.dumps({'case_id':case_id,'state':info['collection_status'],
                          'seconds':info.get('wall_seconds')},ensure_ascii=False),flush=True)

    try:
        pending=[]
        for case in manifest['cases']:
            terminal=root/'cases'/case['case_id']/'terminal.json'
            if terminal.exists():
                data=read(terminal)
                if data['round_identity']!=manifest['round_identity']: raise RuntimeError('Stale result identity')
                if data.get('evidence_manifest'):
                    evidence_list=root/data['evidence_manifest']
                    if digest(evidence_list)!=data['evidence_manifest_sha256']: raise RuntimeError('Prior evidence manifest changed')
                    directory=root/data['evidence_directory']
                    for name,record in read(evidence_list).items():
                        if digest(directory/name)!=record['sha256']: raise RuntimeError('Prior evidence changed: '+case['case_id']+'/'+name)
                # Terminal failures also remain terminal. A repair/retry needs a
                # newly frozen round, not silent replacement of failed evidence.
                results[case['case_id']]=data
            else: pending.append(case)
        progress()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(collect,case) for case in pending]
            for future in as_completed(futures): future.result()
        checked_identity(manifest,root)
        if len(results)!=manifest['case_count']: raise RuntimeError('Some cases lack terminal states')
        final={'status':'all_collection_terminal','round_identity':manifest['round_identity'],
               'all_terminal':True,'expected_cases':manifest['case_count'],'terminal_cases':len(results),
               'collection_counts':{key:sum(x['collection_status']==key for x in results.values())
                 for key in ('complete','completed_with_gaps','error')},
               'wall_seconds_this_session':time.perf_counter()-begin,'finished_epoch':time.time(),
               'scientific_analysis':'ready after this barrier; not yet performed','formal_case_passes':0,
               'cases':[results[key] for key in sorted(results)]}
        atomic_json(root/'collection-complete.json',final)
        atomic_json(root/'progress.json',{k:v for k,v in final.items() if k!='cases'})
        print(json.dumps({k:v for k,v in final.items() if k!='cases'},ensure_ascii=False),flush=True)
        return 0 if all(x['collection_status']=='complete' for x in results.values()) else 2
    finally:
        lock_file.seek(0);msvcrt.locking(lock_file.fileno(),msvcrt.LK_UNLCK,1);lock_file.close()

if __name__=='__main__': raise SystemExit(main())
