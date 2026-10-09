"""Read-only corpus evidence collection with isolated COV processes.

This stage checks evidence integrity only. It does not run scientific reference
comparisons, classify chemical errors, or repair the production implementation.
"""
import argparse
import concurrent.futures
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
LOCK = threading.RLock()
RUNNING = {}
RESULTS = {}
STARTED = time.perf_counter()


def now():
    return datetime.now().astimezone().isoformat()


def save(path, value):
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    temp.replace(path)


def make_manifest(args):
    old = {}
    registry = args.audit_dir/'existing_cases.jsonl'
    if registry.exists():
        for line in registry.read_text(encoding='utf-8-sig').splitlines():
            row = json.loads(line)
            if row.get('path') and row.get('case_id'):
                old[str(Path(row['path']).resolve()).casefold()] = row['case_id']
    logs = {}
    for path in args.logs.rglob('*'):
        if path.is_file() and path.suffix.casefold() in ('.out','.log'):
            logs.setdefault(path.stem.casefold(), []).append(str(path.resolve()))
    cases = []
    for i, path in enumerate(sorted(p for p in args.inputs.rglob('*') if p.is_file() and p.suffix.casefold() in ('.fch','.fchk'))):
        stat = path.stat()
        # Identity and producer records remain raw evidence, not an SCF verdict.
        case_id = old.get(str(path.resolve()).casefold(), f'FILE-{i+1:04}')
        cases.append(dict(case_id=case_id, input=str(path.resolve()), relative_path=str(path.relative_to(args.inputs)),
                          source_size=stat.st_size, source_mtime_ns=stat.st_mtime_ns,
                          log_candidates=sorted(logs.get(path.stem.casefold(), []))))
    if len({x['case_id'] for x in cases}) != len(cases):
        raise RuntimeError('Case IDs are not unique; do not start ambiguous evidence directories')
    if args.case_ids:
        selected = set(args.case_ids.split(','))
        cases = [c for c in cases if c['case_id'] in selected]
        if {c['case_id'] for c in cases} != selected:
            raise RuntimeError('Requested pilot case is absent from the input inventory')
    cases.sort(key=lambda c: c['case_id'])
    return dict(schema=1, created_at=now(), input_root=str(args.inputs.resolve()), recursive=True,
                extensions=['.fch','.fchk'], case_count=len(cases), cases=cases,
                git_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
                workers=args.workers, numerical_threads_per_process=args.threads,
                collection_only=True, scientific_analysis='deferred at user request',
                window_mode='background-hidden with actual production renderer and ImGui back-buffer readback',
                texture_scope='Every alpha/beta MO, 8192 deterministic indices before duplicate removal; not every voxel',
                coverage='Input/log snapshots, production data, all-MO textures, compact members and spin counterparts, structure, controls, actual exports')


def command(op, identity, value=None):
    return op+' '+json.dumps(str(identity),ensure_ascii=False)+((' '+json.dumps(str(value),ensure_ascii=False)) if value is not None else '')


def make_plan(production):
    import numpy as np
    coordinates = np.array(production['atoms'], dtype=float)[:,1:]
    if len(coordinates) > 1:
        _, _, axes = np.linalg.svd(coordinates-coordinates.mean(axis=0), full_matrices=True)
        normal=axes[-1]
    else:
        normal=np.array([0.,0.,1.])
    yaw=math.atan2(normal[2],normal[0]); pitch=math.asin(float(np.clip(normal[1],-1,1)))
    def scene(opacity, yaw_=yaw, pitch_=pitch, iso=0.03):
        return f'scene {opacity} {yaw_:.8f} {pitch_:.8f} 2.2 {iso} 128'
    primary=sorted({i for row in production['compact_rows'] for i in row['members'] if 0<=i<len(production['orbitals'])})
    counterparts=sorted({i for row in production['compact_rows'] for i in row['counterparts'] if 0<=i<len(production['orbitals'])}-set(primary))
    orbitals=production['orbitals']
    occupied=[i for i,m in enumerate(orbitals) if m['occupation']>0]
    selected=max(occupied,key=lambda i: orbitals[i]['energy_hartree']) if occupied else 0
    plan=['COV_VALIDATION 1', scene(0.02), 'capture "structure-face"',scene(0.02,yaw+0.8,pitch+0.35),
          'capture "structure-oblique"',scene(0.92),'seek "panel.browser"','click "browser.filter"',
          'key "Home"','key "Down"','key "Enter"']
    for i in range(len(orbitals)):
        plan += [command('text','browser.search',f'{i+1} '),command('click',f'browser.mo.{i}'),command('volume',f'mo-{i+1:04}',i)]

    def select_browser(i):
        return ['seek "panel.browser"','click "browser.filter"','key "Home"','key "Down"','key "Enter"',
                command('text','browser.search',f'{i+1} '),command('click',f'browser.mo.{i}'),
                'text "browser.search" ""','click "browser.filter"','key "Home"','key "Enter"','seek "panel.diagram"']

    export_names=[]
    def export_view(name):
        export_names.append(name)
        return [command('export-name',name),'click "diagram.export"',
                'hover "scene.viewport"',command('capture',name+'-state')]

    plan += select_browser(primary[0] if primary else selected)
    plan += export_view('export-compact-first')
    for i in primary:
        plan += [command('click',f'diagram.mo.{i}'),command('hover',f'diagram.mo.{i}'),
                 command('capture',f'member-{i+1:04}-tooltip'),'hover "scene.viewport"',command('capture',f'member-{i+1:04}-surface')]
    plan += export_view('export-compact-last')
    # The current product only exposes a merged beta member after browser
    # selection. Record that real path, without changing diagram hit targets.
    for i in counterparts:
        plan += select_browser(i)
        plan += [command('click',f'diagram.mo.{i}'),command('hover',f'diagram.mo.{i}'),
                 command('capture',f'counterpart-{i+1:04}-tooltip'),'hover "scene.viewport"',command('capture',f'counterpart-{i+1:04}-surface')]
    if counterparts: plan += export_view('export-beta-counterpart')
    plan += select_browser(selected)
    plan += ['hover "scene.viewport"','capture "compact"']
    if primary:
        # Restore a primary context for complementary top/bottom captures.
        plan += select_browser(primary[0])
        plan += [command('seek',f'diagram.mo.{primary[-1]}'),'hover "scene.viewport"','capture "compact-top"',
                 command('seek',f'diagram.mo.{primary[0]}'),'hover "scene.viewport"','capture "compact-bottom"']
    plan += ['click "diagram.compact"','hover "scene.viewport"','capture "expanded"',
             'click "diagram.compact"','hover "scene.viewport"','capture "restored"',
             'click "diagram.linear"','hover "scene.viewport"','capture "linear"','click "diagram.nonlinear"',
             'seek "panel.browser"','click "browser.unit"','key "Home"','key "Down"','key "Enter"',
             'seek "panel.diagram"','hover "scene.viewport"','capture "electron-volts"',
             'seek "panel.browser"','click "browser.unit"','key "Home"','key "Enter"',
             'seek "panel.diagram"','click "diagram.export"','hover "scene.viewport"','capture "export-state"']
    # Each state is still reached through the same controls. Preserve its
    # bundle before continuing, instead of overwriting the previous evidence.
    expanded=[]
    for command_line in plan:
        if command_line=='click "diagram.export"' and expanded and expanded[-1]=='seek "panel.diagram"':
            expanded.append(command('export-name','actual-export'))
            export_names.append('actual-export')
        expanded.append(command_line)
        for capture_name in ('expanded','restored','linear','electron-volts'):
            if command_line==command('capture',capture_name):
                expanded += export_view('export-'+capture_name)
                if capture_name=='linear':
                    for i in primary:
                        expanded += [command('click',f'diagram.mo.{i}'),
                            command('hover',f'diagram.mo.{i}'),
                            command('capture',f'linear-member-{i+1:04}')]
                    if primary: expanded += [command('click',f'diagram.mo.{primary[0]}')]
    plan=[]
    for line in expanded:
        plan.append(line)
        if line=='click "diagram.export"':
            plan += ['seek "diagram.export_options"','click "diagram.export_options"','seek "diagram.export_data"','click "diagram.export_data"','click "diagram.export_options"']
    details_mo=primary[0] if primary else selected
    if primary:
        # Navigate through actual controls and wheel events; preserve both ends
        # of the long window, closure, and the reopened selected-member state.
        plan += select_browser(details_mo)
        plan += ['click "diagram.details"','hover "diagram.details.close"','capture "orbital-details-top"',
                 'hover "diagram.details.close.bottom"','capture "orbital-details-bottom"',
                 'click "diagram.details.close.bottom"','hover "scene.viewport"','capture "orbital-details-closed"',
                 'click "diagram.details"','hover "diagram.details.close"','capture "orbital-details-reopened"',
                 'click "diagram.details.close"','hover "scene.viewport"']
    # Cover each new explanation mechanism through drawn controls. Numerical
    # rows only choose what to inspect; they are not an independent oracle.
    detail_requests={}
    def require_detail(mo,target,kind):
        if mo not in primary:return
        detail_requests.setdefault(mo,{})[target]=kind
    seen_scopes=set()
    for row in production['compact_rows']:
        scope=row.get('symmetry_explanation',{})
        signature=(scope.get('origin'),bool(scope.get('local_decomposition')),bool(scope.get('candidate_source')))
        if signature in seen_scopes or not row['members']:continue
        seen_scopes.add(signature);mo=row['members'][0]
        require_detail(mo,'details.symmetry.origin','symmetry-origin')
        require_detail(mo,'details.symmetry.members','symmetry-member-spins')
        require_detail(mo,'details.metal-ligand.end','metal-ligand-field-availability')
        if scope.get('local_decomposition'):
            require_detail(mo,'details.symmetry.represented-rank','local-projection-rank')
        if scope.get('local_assignment'):
            require_detail(mo,'details.symmetry.measure','local-projection-meaning')
    unrepresented_topologies=[]
    for family in production.get('pi_topology_assignments',[]):
        graph=family.get('topology_graph',{})
        if graph.get('source') in (None,'unavailable'):continue
        members=sorted(set(family['orbitals']) & set(primary))
        if not members:
            unrepresented_topologies.append(family['family_id']);continue
        mo=members[0];prefix='pi.rings.'+family['family_id']
        require_detail(mo,prefix+'.source','ring-connectivity-source')
        witnesses=graph.get('channel_ring_witnesses',[])
        if not witnesses and len(family['orientation_channels'])>1:
            require_detail(mo,prefix+'.scope','no-channel-ring-witness-scope')
        for witness in range(len(witnesses)):
            for ring in range(2):require_detail(mo,f'{prefix}.{witness}.path{ring}','ring-path')
            require_detail(mo,f'{prefix}.{witness}.additional-connection','ring-additional-connection')
    gaps=production.get('pi_interactions',[])+production.get('crystal_field_gaps',[])
    for index,gap in enumerate(gaps):
        row_index=gap['retained_level']
        if row_index>=len(production['compact_rows']):continue
        members=production['compact_rows'][row_index]['members']
        if members:require_detail(members[0],f'details.energy-gap.{index}.scope','typed-gap-scope')
    if gaps and primary:
        plan += select_browser(primary[0])
        for unit in range(1,6):
            plan += ['seek "panel.browser"','click "browser.unit"','key "Home"']+['key "Down"']*unit+['key "Enter"','seek "panel.diagram"']
            plan += export_view(f'export-gap-unit-{unit}')
        plan += ['seek "panel.browser"','click "browser.unit"','key "Home"','key "Enter"','seek "panel.diagram"']
    extra_details=[]
    for mo,targets in sorted(detail_requests.items()):
        plan += select_browser(mo)+['click "diagram.details"']
        for target,kind in targets.items():
            capture=f'common-details-{len(extra_details):03d}'
            plan += [command('hover',target),command('capture',capture)]
            extra_details.append(dict(capture=capture,mo=mo,target=target,kind=kind))
        plan += ['click "diagram.details.close"','hover "scene.viewport"']
    for language in range(1,4):
        plan += ['click "language"','key "Home"']+['key "Down"']*language+['key "Enter"','hover "scene.viewport"',command('capture',f'language-{language}')]
    plan += ['click "language"','key "Home"','key "Enter"','hover "scene.viewport"','capture "language-0"',
             scene(0.92,iso=0.02),'capture "iso-002"',scene(0.92,iso=0.04),'capture "iso-004"',scene(0.92),
             'hover "scene.viewport"','capture "final"']
    expected=dict(mos=len(orbitals),primary_members=primary,spin_counterparts=counterparts,
                  export_names=export_names,
                  captures=[json.loads(x[len('capture '):]) for x in plan if x.startswith('capture ')],
                  source_selection_for_controls=primary[0] if primary else selected)
    expected['orbital_details_member']=details_mo if primary else None
    expected['orbital_details_captures']=['orbital-details-top','orbital-details-bottom','orbital-details-closed','orbital-details-reopened'] if primary else []
    expected['common_root_details']=extra_details
    expected['topology_families_outside_compact_rows']=unrepresented_topologies
    return '\n'.join(plan)+'\n', expected


def progress(args, manifest, state='running'):
    with LOCK:
        values=list(RESULTS.values())
        save(args.output/'progress.json',dict(schema=1,state=state,updated_at=now(),total=manifest['case_count'],
             finished=len(values),complete=sum(v['collection_status']=='complete' for v in values),
             gaps=sum(v['collection_status']=='completed_with_gaps' for v in values),
             errors=sum(v['collection_status']=='error' for v in values),elapsed_seconds=time.perf_counter()-STARTED,
             running=dict(RUNNING),scientific_analysis='deferred at user request'))


def run_process(args, manifest, case_id, exe, arguments, destination, stdout_name, timeout):
    startup=subprocess.STARTUPINFO();startup.dwFlags|=subprocess.STARTF_USESHOWWINDOW;startup.wShowWindow=0
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
        env[name]=str(args.threads)
    with (destination/stdout_name).open('wb') as stdout, (destination/(stdout_name+'.stderr.txt')).open('wb') as stderr:
        proc=subprocess.Popen([str(exe),*map(str,arguments)],cwd=destination,stdout=stdout,stderr=stderr,env=env,startupinfo=startup)
        with LOCK:
            RUNNING[case_id]=dict(stage=exe.name,pid=proc.pid,started_at=now(),directory=str(destination))
            progress(args,manifest)
        try:
            code=proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill();proc.wait()
            raise RuntimeError(f'{exe.name} exceeded {timeout} seconds; partial evidence retained')
    return dict(pid=proc.pid,exit_code=code,executable=str(exe))


def collect_one(args, manifest, record):
    from PIL import Image,ImageChops
    case_id=record['case_id']; directory=args.output/'cases'/case_id
    directory.mkdir(parents=True)
    begin=time.perf_counter(); data=dict(case_id=case_id,input=record['input'],started_at=now(),collection_status='error',scientific_analysis='deferred')
    try:
        source=Path(record['input']); before=source.stat()
        if before.st_size!=record['source_size'] or before.st_mtime_ns!=record['source_mtime_ns']:
            raise RuntimeError('Input changed since manifest creation')
        shutil.copy2(source,directory/('source'+source.suffix.lower()))
        log_copies=[]
        for n, log in enumerate(record['log_candidates']):
            path=Path(log); dest=directory/f'source-log-{n+1}{path.suffix.lower()}'
            shutil.copy2(path,dest);log_copies.append(dict(original=str(path),copy=dest.name,size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
        save(directory/'source-identity.json',dict(**record,snapshot='source'+source.suffix.lower(),logs=log_copies))
        t=time.perf_counter()
        dump=run_process(args,manifest,case_id,args.build/'cov_scientific_audit_dump.exe',[source,'--interactions'],directory,'production.json',900)
        data['inspection_process']=dump;data['inspection_seconds']=time.perf_counter()-t
        if dump['exit_code']:
            raise RuntimeError('Production inspection failed; source and failure log retained')
        production=json.loads((directory/'production.json').read_text(encoding='utf-8'))
        plan,expected=make_plan(production)
        (directory/'case.plan').write_text(plan,encoding='utf-8');save(directory/'expected-evidence.json',expected)
        t=time.perf_counter()
        native=run_process(args,manifest,case_id,args.build/'cov_validation.exe',
                 [source,'--validation-plan',directory/'case.plan','--validation-output',directory/'native','--validation-background'],
                 directory,'native.log',args.timeout)
        data['native_process']=native;data['native_seconds']=time.perf_counter()-t
        native_dir=directory/'native'; session_path=native_dir/'session.json'
        session=json.loads(session_path.read_text(encoding='utf-8')) if session_path.exists() else None
        data['native_session']=session
        missing=[];invalid=[];mo_ids=[];points=[];full_grids=[]
        for path in sorted(native_dir.glob('*.volume.json')):
            try:
                volume=json.loads(path.read_text(encoding='utf-8'))
                if 'full_grid' in volume:
                    grid=volume['full_grid']
                    binary=native_dir/grid['file']
                    shape=[volume[k] for k in ('nx','ny','nz')]
                    if not binary.is_file() or binary.stat().st_size!=4*math.prod(shape) or grid['point_count']!=math.prod(shape):
                        invalid.append(path.name+': incomplete full texture')
                    full_grids.append(dict(mo=volume['rendered_mo'],shape=shape,metadata=path.name,binary=grid['file']))
                else:
                    mo_ids.append(volume['rendered_mo']);points.append(len(volume['samples']))
                if not volume['samples'] or any(not math.isfinite(v) for _,v in volume['samples']):invalid.append(path.name)
            except (KeyError,ValueError):invalid.append(path.name)
        missing += [f'MO {i+1} texture' for i in range(expected['mos']) if i not in mo_ids]
        if len(mo_ids)!=len(set(mo_ids)):invalid.append('duplicate MO texture identity')
        for i in expected.get('full_grid_frontiers_zero_based',[]):
            for resolution in expected.get('full_grid_resolutions',[]):
                found=[g for g in full_grids if g['mo']==i and g['shape']==[resolution]*3]
                if len(found)!=1:missing.append(f'MO {i+1} complete {resolution} texture')
        images=[]
        for path in sorted(native_dir.glob('*.bmp')):
            target=path.with_suffix('.png')
            with Image.open(path) as raw:
                rgb=raw.convert('RGB');rgb.save(target)
                with Image.open(target) as packaged:
                    if ImageChops.difference(rgb,packaged.convert('RGB')).getbbox() is not None:
                        raise RuntimeError('Lossless screenshot conversion verification failed')
                    if packaged.size!=(2100,1250):invalid.append(f'{target.name}: unexpected dimensions')
            # Only this run's verified lossless replacement removes its own BMP.
            path.unlink();images.append(target.name)
        selected=[]
        for capture in expected['captures']:
            path=native_dir/(capture+'.ui.json')
            if not path.exists() or not (native_dir/(capture+'.png')).exists():
                missing.append(capture);continue
            shot=json.loads(path.read_text(encoding='utf-8'))
            if capture.startswith(('member-','counterpart-')) and capture.endswith('-surface'):
                wanted=int(capture.split('-')[1])-1;actual=shot['state']['applied_mo']
                selected.append(dict(capture=capture,requested_mo=wanted,observed_mo=actual))
                if wanted!=actual:invalid.append(f'{capture}: selected {actual}, requested {wanted}')
        for export_name in expected.get('export_names',['actual-export']):
            for suffix in ('png','svg','json','csv'):
                path=native_dir/(export_name+'.mo.'+suffix)
                if not path.exists() or path.stat().st_size==0:missing.append(path.name)
        failed_actions=[]
        if (native_dir/'actions.jsonl').exists():
            failed_actions=[x for x in (json.loads(line) for line in (native_dir/'actions.jsonl').read_text(encoding='utf-8').splitlines()) if x['status']!='executed']
        if not session:missing.append('completed native session')
        if native['exit_code'] not in (0,2):invalid.append(f'native exit {native["exit_code"]}')
        identity=json.loads((native_dir/'identity.json').read_text(encoding='utf-8'))
        if identity['git_commit']!=manifest['git_commit']:invalid.append('build identity differs from frozen manifest')
        after=source.stat()
        if after.st_size!=before.st_size or after.st_mtime_ns!=before.st_mtime_ns:invalid.append('source changed during collection')
        data.update(collection_status='complete' if not missing and not invalid and not failed_actions else 'completed_with_gaps',
                    orbitals=expected['mos'],texture_count=len(mo_ids),sample_points_min=min(points,default=0),sample_points_max=max(points,default=0),
                    screenshots=len(images),export_png=(native_dir/'actual-export.mo.png').exists(),
                    complete_texture_grids=full_grids,
                    primary_members=len(expected['primary_members']),spin_counterparts=len(expected['spin_counterparts']),
                    selection_records=selected,missing=missing,integrity_errors=invalid,failed_actions=failed_actions,
                    note='Collection completeness only; scientific correctness and visual/chemical review deferred.')
        cards=''.join(f'<figure><a href="native/{html.escape(name)}"><img loading="lazy" src="native/{html.escape(name)}"></a><figcaption>{html.escape(name)}</figcaption></figure>' for name in images+(['actual-export.mo.png'] if data['export_png'] else []))
        page='<!doctype html><meta charset="utf-8"><title>'+case_id+' 取证</title><style>body{font:16px system-ui;margin:30px;background:#10141c;color:#edf3ff}a{color:#83bfff}.gallery{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}figure{margin:0}img{width:100%}pre{white-space:pre-wrap}</style>'
        page+=f'<h1>{case_id} 原始取证</h1><p>{html.escape(record["relative_path"])}</p><p>科学分析、图像和化学解释复核待下一阶段。当前状态仅表示证据采集完整性。</p><p><a href="collection.json">取证记录</a> · <a href="production.json">生产数据</a> · <a href="source-identity.json">原始文件身份</a></p><div class="gallery">{cards}</div>'
        (directory/'review.html').write_text(page,encoding='utf-8')
    except Exception as error:
        data['error']=f'{type(error).__name__}: {error}'
    data['finished_at']=now();data['wall_seconds']=time.perf_counter()-begin
    save(directory/'collection.json',data)
    with LOCK:
        RUNNING.pop(case_id,None);RESULTS[case_id]=data
        with (args.output/'cases.jsonl').open('a',encoding='utf-8') as stream:stream.write(json.dumps(data,ensure_ascii=False)+'\n')
        progress(args,manifest)
        print(json.dumps({k:data.get(k) for k in ('case_id','collection_status','orbitals','texture_count','screenshots','wall_seconds','error')},ensure_ascii=False),flush=True)
    return data


def main():
    ap=argparse.ArgumentParser()
    for name in ('inputs','logs','build','audit-dir','reference-deps','output'):
        ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--workers',type=int,default=4);ap.add_argument('--threads',type=int,default=3)
    ap.add_argument('--timeout',type=int,default=1800);ap.add_argument('--case-ids')
    args=ap.parse_args()
    args.output=args.output.resolve()
    if args.output.exists():raise RuntimeError('Use a fresh batch output directory')
    if not 1<=args.workers<=8:raise ValueError('Use between 1 and 8 independent processes')
    os.environ['OPENBLAS_NUM_THREADS']=str(args.threads)
    sys.path[:0]=[str(args.reference_deps),str(args.build.parent/'python-deps')]
    args.output.mkdir(parents=True)
    manifest=make_manifest(args);save(args.output/'manifest.json',manifest)
    progress(args,manifest)
    print(json.dumps(dict(batch=str(args.output),total=manifest['case_count'],workers=args.workers,started_at=now())),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers,thread_name_prefix='cov-evidence') as pool:
        pending=[pool.submit(collect_one,args,manifest,c) for c in manifest['cases']]
        for future in concurrent.futures.as_completed(pending):future.result()
    progress(args,manifest,'finished')
    values=sorted(RESULTS.values(),key=lambda x:x['case_id'])
    summary=dict(schema=1,state='finished',finished_at=now(),git_commit=manifest['git_commit'],files=len(values),
                 complete=sum(x['collection_status']=='complete' for x in values),
                 completed_with_gaps=sum(x['collection_status']=='completed_with_gaps' for x in values),
                 errors=sum(x['collection_status']=='error' for x in values),elapsed_seconds=time.perf_counter()-STARTED,
                 total_orbitals=sum(x.get('orbitals',0) for x in values),actual_textures=sum(x.get('texture_count',0) for x in values),
                 framebuffer_screenshots=sum(x.get('screenshots',0) for x in values),
                 actual_export_pngs=sum(bool(x.get('export_png')) for x in values),
                 scientific_analysis='deferred at user request',visual_review='deferred at user request',
                 attention=[x['case_id'] for x in values if x['collection_status']!='complete'])
    rows=''.join(f'<tr><td><a href="cases/{x["case_id"]}/review.html">{x["case_id"]}</a></td><td>{html.escape(Path(x["input"]).name)}</td><td>{x["collection_status"]}</td><td>{x.get("texture_count",0)}/{x.get("orbitals","?")}</td><td>{x.get("screenshots",0)}</td><td>{x["wall_seconds"]:.1f}s</td></tr>' for x in values)
    page='<!doctype html><meta charset="utf-8"><title>FCHK 全库并行取证</title><style>body{font:16px system-ui;margin:30px;background:#10141c;color:#edf3ff}a{color:#83bfff}td,th{padding:8px;border-bottom:1px solid #344255}table{border-collapse:collapse}</style>'
    page+=f'<h1>FCHK 全库并行取证</h1><p>本阶段只采集证据并检查文件完整性；科学判定、图像解释和纠错分析尚未进行。</p><p>{len(values)} 文件，{summary["complete"]} 完整，{summary["completed_with_gaps"]} 有缺口，{summary["errors"]} 运行错误。墙钟 {summary["elapsed_seconds"]:.1f} 秒。</p><p><a href="summary.json">汇总</a> · <a href="manifest.json">输入清单与构建身份</a></p><table><tr><th>案例</th><th>文件</th><th>取证状态</th><th>实际纹理/MO</th><th>截图</th><th>耗时</th></tr>{rows}</table>'
    (args.output/'index.html').write_text(page,encoding='utf-8')
    save(args.output/'summary.json',summary)
    save(args.output/'COMPLETED.json',dict(finished_at=now(),summary='summary.json',collection_only=True))
    print(json.dumps(summary,ensure_ascii=False),flush=True)
    return 0 if not summary['errors'] and not summary['completed_with_gaps'] else 2


if __name__=='__main__':
    sys.exit(main())
