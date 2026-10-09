"""One-molecule native collection + independent checking, with separate timings.

All molecular selections come from the input/reference or production inspection;
no molecule-name, fixed orbital-number or atom-order special cases are used.
Requires the retained audit helpers and independently installed IOData/GBasis.
The process returns 2 for scientific/UI failures, even when collection succeeds.
"""
import argparse
import csv
import html
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

START = time.perf_counter()
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '4'


def main():
    ap = argparse.ArgumentParser()
    for name in ('input', 'log', 'build', 'off-build', 'audit-dir', 'reference-deps', 'output'):
        ap.add_argument('--'+name, type=Path, required=True)
    ap.add_argument('--case-id', default='case')
    ap.add_argument('--group', required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise RuntimeError('Use a new output directory; existing evidence is never overwritten')
    args.output.mkdir(parents=True)
    out = args.output
    sys.path[:0] = [str(args.reference_deps), str(args.audit_dir), str(args.build.parent/'python-deps')]
    import numpy as np
    from iodata import load_one
    from gbasis.wrappers import from_iodata
    from gbasis.integrals.overlap import overlap_integral
    from gbasis.evals.eval import evaluate_basis
    from PIL import Image
    import audit_existing
    import audit_symmetry

    timing = {'runtime_imports_seconds': time.perf_counter()-START}
    checks = []
    def save(name, value):
        (out/name).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    def check(id, status, observed, tolerance=None):
        checks.append(dict(check_id=id, status=status, observed=observed, tolerance=tolerance))
    def process(exe, extra, name):
        p = subprocess.run([str(exe), *map(str, extra)], capture_output=True, encoding='utf-8', errors='replace', timeout=900)
        (out/name).write_text(p.stdout+p.stderr, encoding='utf-8')
        return p
    def elapsed(phase, begin):
        timing[phase] = time.perf_counter()-begin
        print(phase, round(timing[phase], 3), flush=True)

    t = time.perf_counter()
    audit_existing.BUILD = args.build
    audit_existing.COMMIT = subprocess.check_output(['git', '-C', str(Path(__file__).resolve().parents[1]), 'rev-parse', 'HEAD'], text=True).strip()
    identity = audit_existing.audit(args.input, args.case_id, {args.input.stem.casefold(): [args.log]})
    save('input.json', identity)
    _, _, fields = audit_existing.read_fchk(args.input)
    check('INPUT-01', 'pass' if not identity['issues'] else 'fail', identity['issues'])
    check('INPUT-02', 'pass' if identity['producer_quality']['identity_verified'] else 'insufficient', identity['producer_quality'])
    check('INPUT-03', 'insufficient' if not identity['producer_quality']['optimization_completed'] else 'review',
          'Input is assessed as the supplied calculation; a fixed-geometry SP does not establish an optimized minimum or stable ground state.')
    proc = process(args.build/'cov_scientific_audit_dump.exe', [args.input, '--interactions'], 'production.json')
    if proc.returncode:
        raise RuntimeError('Production inspection failed')
    production = json.loads(proc.stdout)
    off = process(args.off_build/'cov_scientific_audit_dump.exe', [args.input, '--interactions'], 'ordinary-build.json')
    check('NATIVE-01', 'pass' if proc.stdout == off.stdout and off.returncode == 0 else 'fail',
          'ON/OFF production parser, overlap, chemistry, connection graph and compact data compared byte-for-byte')
    raw_coeff = fields['Alpha MO coefficients'].reshape(-1, production['nbasis'])
    if 'Beta MO coefficients' in fields:
        raw_coeff = np.vstack([raw_coeff, fields['Beta MO coefficients'].reshape(-1, production['nbasis'])])
    cov_coeff = np.asarray([mo['coefficients'] for mo in production['orbitals']])
    coeff_error = float(np.max(np.abs(cov_coeff-raw_coeff)))
    energy_ref = fields['Alpha Orbital Energies']
    if 'Beta Orbital Energies' in fields:
        energy_ref = np.r_[energy_ref, fields['Beta Orbital Energies']]
    energy_error = float(np.max(np.abs(np.array([x['energy_hartree'] for x in production['orbitals']])-energy_ref)))
    check('INPUT-04', 'pass' if coeff_error <= 1e-5 and energy_error <= 1e-10 else 'fail',
          dict(coefficient_max_abs=coeff_error, energy_max_abs_hartree=energy_error), 'float coefficients 1e-5; energies 1e-10 Ha')
    elapsed('input_and_production_seconds', t)

    t = time.perf_counter()
    mol = load_one(str(args.input))
    basis = from_iodata(mol)
    overlap = overlap_integral(basis)
    eig = np.linalg.eigvalsh(overlap)
    blocks = [('alpha', mol.mo.coeffsa)]
    if mol.mo.kind == 'unrestricted':
        blocks.append(('beta', mol.mo.coeffsb))
    orthogonal = {}
    for spin, coeff in blocks:
        delta = coeff.T@overlap@coeff-np.eye(coeff.shape[1])
        diag = float(np.max(np.abs(np.diag(delta))))
        np.fill_diagonal(delta, 0)
        orthogonal[spin] = dict(nmo=coeff.shape[1], normalization_max_abs=diag, offdiagonal_max_abs=float(np.max(np.abs(delta))))
    density = (mol.mo.coeffs*mol.mo.occs)@mol.mo.coeffs.T
    trace = float(np.einsum('ij,ji->', density, overlap))
    check('NUM-01', 'pass' if all(max(x['normalization_max_abs'], x['offdiagonal_max_abs']) <= 2e-6 for x in orthogonal.values()) and abs(trace-sum(mol.mo.occs))<2e-5 else 'fail',
          dict(blocks=orthogonal, trace_ps=trace, condition_number=float(eig[-1]/eig[0])), '2e-6 orthogonality; 2e-5 electrons')
    # Gaussian Cartesian-g ordering differs from the COV internal order. Pure
    # spherical shells keep producer order here; their real-space phase is
    # checked independently below, never silently corrected in production.
    permutation, owners = [], []
    for shell, atom in zip(fields['Shell types'], fields['Shell to atom map']):
        n = 4 if shell == -1 else 2*abs(shell)+1 if shell < -1 else (shell+1)*(shell+2)//2
        order = [14,4,0,13,12,8,3,5,1,11,9,2,10,7,6] if shell == 4 else list(range(n))
        offset=len(permutation); permutation.extend(offset+i for i in order);owners.extend([int(atom)-1]*n)
    cov_s = np.asarray(production['overlap'])
    s_error = float(np.max(np.abs(cov_s.reshape(overlap.shape)-overlap[np.ix_(permutation,permutation)]))) if cov_s.size == overlap.size else None
    check('NUM-02', 'pass' if s_error is not None and s_error <= 2e-5 else 'fail', s_error, 'S maximum absolute difference <=2e-5')
    ps=density@overlap; products=ps*ps.T
    if mol.mo.kind == 'unrestricted':
        spin_density=(mol.mo.coeffsa*mol.mo.occsa)@mol.mo.coeffsa.T-(mol.mo.coeffsb*mol.mo.occsb)@mol.mo.coeffsb.T
        qs=spin_density@overlap;products+=qs*qs.T
    owners=np.array(owners); atom_ids=[np.flatnonzero(owners==i) for i in range(len(mol.atnums))]
    mayer=[]
    for a,b,value in production['bond_orders']:
        ref=float(products[np.ix_(atom_ids[int(a)],atom_ids[int(b)])].sum())
        mayer.append(dict(atoms=[int(a),int(b)], cov=value, reference=ref, error=abs(value-ref)))
    save('mayer.json', mayer)
    check('NUM-03', 'pass' if max((x['error'] for x in mayer),default=0)<=2e-4 else 'fail',
          dict(reported_pairs=len(mayer), max_abs=max((x['error'] for x in mayer), default=0)), '2e-4')
    invalid=[]
    for i,c in enumerate(production['chemistry']):
        for name in ('bonding','channel','manifold'):
            v=np.asarray(c[name]);status=c.get(name+'_status')
            if status==3 and np.all(v==0):continue
            if c['available'] and (not np.isfinite(v).all() or v.min() < -1e-8 or v.max()>1+1e-8 or abs(v.sum()-1)>2e-6):
                invalid.append(dict(source_mo=i+1, distribution=name, values=v.tolist()))
    save('distribution-errors.json',invalid)
    check('CHEM-01', 'fail' if invalid else 'pass',dict(invalid_distributions=len(invalid)), 'probabilities in [0,1], sum=1; NotApplicable zero allowed')
    sym=audit_symmetry.audit(dict(case_id=args.case_id,path=str(args.input),point_group=args.group))
    save('symmetry.json', sym)
    labels={i:g['label'].upper() for g in sym.get('groups',[]) if g['label'] for i in g['internal_indices']}
    compact=sorted({i for row in production['compact_rows'] for i in row['members']})
    sym_bad=[i for i in compact if labels.get(i)!=production['chemistry'][i]['symmetry'].upper()]
    check('SYM-01','pass' if not sym_bad and len(labels)>=len(compact) else 'fail',dict(labelled_mos=len(labels),compact_mismatches=sym_bad,operations=sym['operations']))
    # Pair chemistry keeps its physical scope explicit: local adjacent bonds
    # and remote pairs are retained separately; no occupied/virtual shortcut.
    pairs=[]
    for i in compact:
        c=mol.mo.coeffs[:,i]; row=[]
        for edge in production['render_edges']:
            a,b=edge['atoms'];ia,ib=atom_ids[a],atom_ids[b]
            row.append(dict(atoms=[a,b],kind=edge['kind'],overlap_per_unit_occupation=float(2*c[ia]@overlap[np.ix_(ia,ib)]@c[ib])))
        pairs.append(dict(internal_mo=i,source_mo=i+1,pairs=row))
    save('pi-pair-reference.json',pairs)
    elapsed('independent_integrals_symmetry_seconds', t)

    # Scene directions are derived from the geometry. These are scene settings,
    # while every MO selection/hover/filter/expand/export is normal UI input.
    _,_,axes=np.linalg.svd(mol.atcoords-mol.atcoords.mean(axis=0))
    normal=axes[-1]
    yaw=math.atan2(normal[2],normal[0]);pitch=math.asin(float(normal[1]))
    def scene(opacity,yaw_=yaw,pitch_=pitch,iso=0.03,res=128):
        return f'scene {opacity} {yaw_:.8f} {pitch_:.8f} 2.2 {iso} {res}'
    def command(op,id,value=None):
        return op+' '+json.dumps(str(id),ensure_ascii=False)+((' '+json.dumps(str(value),ensure_ascii=False)) if value is not None else '')
    plan=['COV_VALIDATION 1',scene(0.02),'capture "structure-face"',scene(0.02,yaw+0.8,pitch+0.35),
          'capture "structure-oblique"',scene(0.92), 'seek "panel.browser"', 'click "browser.filter"','key "Home"','key "Down"','key "Enter"']
    for i in range(len(production['orbitals'])):
        plan += [command('text','browser.search',f'{i+1} '), command('click',f'browser.mo.{i}'), command('volume',f'mo-{i+1:03}',i)]
    plan += ['click "browser.filter"','key "Home"','key "Enter"','text "browser.search" ""','seek "panel.diagram"']
    for i in compact:
        plan += [command('click',f'diagram.mo.{i}'),command('hover',f'diagram.mo.{i}'),command('capture',f'pi-{i+1:03}-tooltip'),
                 'hover "scene.viewport"',command('capture',f'pi-{i+1:03}-surface')]
    homo=int(np.flatnonzero(mol.mo.occs>0)[-1]);lumo=int(np.flatnonzero(mol.mo.occs==0)[0])
    plan += [command('click',f'diagram.mo.{homo}'),'hover "scene.viewport"','capture "compact"',
             command('seek',f'diagram.mo.{compact[-1]}'),'hover "scene.viewport"','capture "compact-top"',
             command('seek',f'diagram.mo.{compact[0]}'),'hover "scene.viewport"','capture "compact-bottom"',
             'click "diagram.compact"','hover "scene.viewport"','capture "expanded"','click "diagram.compact"','hover "scene.viewport"','capture "restored"',
             'click "diagram.linear"','hover "scene.viewport"','capture "linear"','click "diagram.nonlinear"',
             'seek "panel.browser"','click "browser.unit"','key "Home"','key "Down"','key "Enter"','seek "panel.diagram"','hover "scene.viewport"','capture "electron-volts"',
             'seek "panel.browser"','click "browser.unit"','key "Home"','key "Enter"',
             'seek "panel.diagram"','click "diagram.export"','seek "diagram.export_options"','click "diagram.export_options"','seek "diagram.export_data"','click "diagram.export_data"','click "diagram.export_options"','hover "scene.viewport"','capture "export-state"']
    for language in range(1,4):
        plan += ['click "language"','key "Home"']+['key "Down"']*language+['key "Enter"','hover "scene.viewport"',command('capture',f'language-{language}')]
    plan += ['click "language"','key "Home"','key "Enter"','hover "scene.viewport"','capture "language-0"',
             scene(0.92,iso=0.02),'capture "iso-002"',scene(0.92,iso=0.04),'capture "iso-004"',scene(0.92),
             'hover "scene.viewport"','capture "final"']
    plan_path=out/'case.plan';plan_path.write_text('\n'.join(plan)+'\n',encoding='utf-8')
    t=time.perf_counter()
    native=process(args.build/'cov_validation.exe',[args.input,'--validation-plan',plan_path,'--validation-output',out/'native'],'native.log')
    elapsed('native_collection_wall_seconds',t)
    if not (out/'native'/'session.json').exists():
        raise RuntimeError('Native run crashed or did not finish; see native.log')
    session=json.loads((out/'native'/'session.json').read_text(encoding='utf-8'))
    save('native-session.json',session)
    check('NATIVE-02','pass' if native.returncode==0 and session['failed_commands']==0 else 'fail',session)

    t=time.perf_counter();cache={};volumes=[]
    odd_m_phase=[]
    for shell in fields['Shell types']:
        if shell < -1:
            l=abs(int(shell));odd_m_phase += [1]+[(-1)**m for m in range(1,l+1) for _ in range(2)]
        else:
            n=4 if shell == -1 else (shell+1)*(shell+2)//2
            odd_m_phase += [1]*n
    odd_m_phase=np.array(odd_m_phase)
    for p in sorted((out/'native').glob('*.volume.json')):
        d=json.loads(p.read_text(encoding='utf-8'));mo=d['rendered_mo'];samples=np.array(d['samples'])
        ids=samples[:,0].astype(int);actual=samples[:,1]
        key=(tuple(d['grid_box_bohr']),d['nx'],d['ny'],d['nz'],tuple(ids))
        if key not in cache:
            ijk=np.column_stack((ids%d['nx'],(ids//d['nx'])%d['ny'],ids//(d['nx']*d['ny']))).astype(np.float32)
            lo=np.asarray(d['grid_box_bohr'][:3],dtype=np.float32);hi=np.asarray(d['grid_box_bohr'][3:],dtype=np.float32)
            # One final rounding models CUDA's fused multiply-add; the residual
            # coordinate ambiguity is far smaller than the documented tolerance.
            frac=ijk/np.array([d['nx']-1,d['ny']-1,d['nz']-1],dtype=np.float32)
            points=(lo.astype(float)+frac.astype(float)*(hi-lo).astype(float)).astype(np.float32).astype(float)
            aos=evaluate_basis(basis,points)
            cache[key]=(mol.mo.coeffs.T@aos,(mol.mo.coeffs*odd_m_phase[:,None]).T@aos)
        ref=cache[key][0][mo];diagnostic=cache[key][1][mo]
        diagnostic_nrms=float(np.linalg.norm(actual-diagnostic)/max(np.linalg.norm(diagnostic),1e-30))
        dot=float(actual@ref);sign=1 if dot>=0 else -1
        nrms=float(np.linalg.norm(actual-sign*ref)/max(np.linalg.norm(ref),1e-30))
        cosine=abs(dot)/max(float(np.linalg.norm(actual)*np.linalg.norm(ref)),1e-30)
        relative_max=float(np.max(np.abs(actual-sign*ref))/max(np.max(np.abs(ref)),1e-30))
        volumes.append(dict(source_mo=mo+1,points=len(ids),frame=d['frame'],generation=d['generation'],nrms=nrms,cosine=cosine,relative_max=relative_max,phase=sign,known_odd_m_phase_diagnostic_nrms=diagnostic_nrms,
                            passed=nrms<=1e-4 and cosine>=1-1e-7 and relative_max<=1e-3))
    save('actual-volume-comparison.json',volumes)
    check('NUM-05','diagnostic',dict(known_odd_m_phase_explained=sum(x['known_odd_m_phase_diagnostic_nrms']<=1e-4 for x in volumes),total=len(volumes),
          worst_diagnostic_nrms=max((x['known_odd_m_phase_diagnostic_nrms'] for x in volumes),default=None)),
          'Reference-only phase experiment. Production remains unchanged; diagnostic agreement never changes NUM-04 failure to pass.')
    check('NUM-04','pass' if len(volumes)==len(production['orbitals']) and all(x['passed'] for x in volumes) else 'fail',
          dict(expected=len(production['orbitals']),readbacks=len(volumes),passed=sum(x['passed'] for x in volumes),failed=sum(not x['passed'] for x in volumes),worst_nrms=max((x['nrms'] for x in volumes),default=None)),
          'All MOs; 8192 deterministic sample indices per actual texture, duplicates removed; NRMS 1e-4; abs cosine >=1-1e-7; relative max 1e-3. Not every voxel or all-space proof.')
    frame_rows=[json.loads(x) for x in (out/'native'/'frames.jsonl').read_text(encoding='utf-8').splitlines()]
    transition=[x for x in frame_rows if x['rendered_generation']!=x['volume_generation']]
    mismatches=[x for x in frame_rows if x['rendered_mo']!=x['drawn_ui_mo']]
    check('NATIVE-03','pass' if not mismatches else 'fail',dict(frames=len(frame_rows),deferred_update_frames=len(transition),scene_ui_mismatches=len(mismatches)),
          'Scene and UI were drawn with the same MO. A later applied selection may differ until the next frame; this is recorded, not hidden.')
    shots={p.stem.removesuffix('.ui'):json.loads(p.read_text(encoding='utf-8')) for p in (out/'native').glob('*.ui.json')}
    def members(name):return [x['data']['used_internal_mo'] for x in shots[name]['draw_trace'] if x['kind']=='draw.level']
    check('UI-01','pass' if members('compact')==compact and members('restored')==compact and shots['expanded']['state']['compact'] is False else 'fail',
          dict(compact=members('compact'),expanded=members('expanded'),restored=members('restored')))
    exp=json.loads((out/'native'/'actual-export.mo.json').read_text(encoding='utf-8'))
    exp_members=sorted({i for x in exp['diagram_rows'] for i in x['member_orbitals']})
    check('EXPORT-01','pass' if exp_members==compact else 'fail',dict(gui_members=compact,export_members=exp_members,export_mode=exp['mode'],export_rows=exp['diagram_row_count']))
    check('EXPORT-02','pass' if exp['electronic_state']['charge']==identity['charge'] and exp['electronic_state']['multiplicity']==identity['multiplicity'] else 'fail',exp['electronic_state'])
    # Values in JSON/CSV must still retain producer energy even if the wrong
    # builder selected their rows. This does not excuse EXPORT-01/02.
    energy_errors=[abs(x['energy_hartree']-energy_ref[x['index']]) for x in exp['orbitals']]
    check('EXPORT-03','pass' if max(energy_errors,default=0)<=1e-10 else 'fail',dict(max_energy_error_hartree=max(energy_errors,default=0)))
    csv_rows=list(csv.DictReader((out/'native'/'actual-export.mo.csv').open(encoding='utf-8')))
    csv_errors=[]
    for c,j in zip(csv_rows,exp['orbitals']):
        if int(c['index'])!=j['index'] or abs(float(c['energy_hartree'])-j['energy_hartree'])>1e-10 or float(c['occupation'])!=j['occupation'] or c['spin']!=j['spin'] or c['symmetry']!=j['symmetry'] or (c['included_in_diagram']=='1')!=j['included_in_diagram']:
            csv_errors.append(c['index'])
    import xml.etree.ElementTree as ET
    svg=ET.parse(out/'native'/'actual-export.mo.svg').getroot()
    check('EXPORT-04','pass' if not csv_errors and len(csv_rows)==len(exp['orbitals']) and int(svg.attrib['data-diagram-row-count'])==exp['diagram_row_count'] else 'fail',
          dict(csv_rows=len(csv_rows),csv_value_discrepancies=csv_errors,svg_rows=svg.attrib['data-diagram-row-count']),
          'Cross-format metadata consistency; PNG/SVG appearance still requires image review and cannot excuse wrong GUI/export membership.')
    controls={}
    for name in ('electron-volts','language-0','language-1','language-2','language-3','iso-002','iso-004','final'):
        d=shots[name]; languages=[x['data']['value'] for x in d['draw_trace'] if x['kind']=='draw.text' and x['data']['label']=='language']
        controls[name]=dict(selected=d['state']['applied_mo'],volume_generation=d['state']['volume_generation'],unit=d['state']['energy_unit'],language=languages)
    correct_languages=all(controls[f'language-{i}']['language']==[str(i)] for i in range(4))
    check('UI-02','pass' if correct_languages and controls['electron-volts']['unit']==1 and len({x['volume_generation'] for x in controls.values()})==1 and all(x['selected']==homo for x in controls.values()) else 'fail',controls)
    selection_errors=[i for i in compact if shots[f'pi-{i+1:03}-surface']['state']['applied_mo']!=i]
    check('UI-03','pass' if not selection_errors else 'fail',dict(individually_selected_pi_mos=len(compact),selection_errors=selection_errors))
    elapsed('independent_volume_and_ui_check_seconds',t)

    t=time.perf_counter()
    images=[]
    for p in sorted((out/'native').glob('*.bmp')):
        target=p.with_suffix('.png'); im=Image.open(p).convert('RGB');im.save(target)
        # These BMPs were created only by this run; the verified lossless PNG
        # replaces them. No original/calculation/older evidence is touched.
        if np.array_equal(np.asarray(im),np.asarray(Image.open(target))):p.unlink()
        images.append(target)
    elapsed('evidence_packaging_seconds',t)
    timing['total_automated_seconds']=time.perf_counter()-START
    summary=dict(case_id=args.case_id,input=str(args.input),timing=timing,
                 scientific_verdict='fail' if any(x['status']=='fail' for x in checks) else 'review',
                 visual_review='pending human/model image review',checks=checks,
                 scope='One supplied molecule. No new Gaussian SCF/optimization calculation. Native collection plus all-MO sampled texture comparison; additional scientific/visual limitations remain explicit.')
    save('result.json',summary)
    cards=''.join(f'<figure><a href="native/{p.name}" target="_blank"><img loading="lazy" src="native/{p.name}"></a><figcaption>{html.escape(p.stem)}</figcaption></figure>' for p in images)
    rows=''.join('<tr><td>'+html.escape(c['check_id'])+'</td><td>'+c['status']+'</td><td><pre>'+html.escape(json.dumps(c['observed'],ensure_ascii=False,indent=2))+'</pre></td></tr>' for c in checks)
    page='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>单分子验证审阅</title><style>
body{font:16px system-ui;margin:32px;background:#10141c;color:#edf3ff}a{color:#83bfff}h1{font-size:26px}table{border-collapse:collapse;width:100%}td,th{border:1px solid #344255;padding:12px;vertical-align:top}pre{white-space:pre-wrap;max-height:260px;overflow:auto;margin:0}.gallery{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}figure{margin:0}img{width:100%;border:1px solid #344255}figcaption{padding:6px}small{color:#b4c5db}</style>'''
    page+=f'<h1>{html.escape(args.case_id)} 单分子原生验证</h1><p>自动流程 {timing["total_automated_seconds"]:.3f} 秒；其中界面/GPU取证 {timing["native_collection_wall_seconds"]:.3f} 秒。科学判定：{summary["scientific_verdict"]}。图像待逐张复核。</p><p>点击图像可查看原尺寸。输入身份、参考值、逐帧状态和原始导出均保存在同目录。</p><table><tr><th>检查</th><th>状态</th><th>证据/数值</th></tr>{rows}</table><h2>实际帧缓冲图像</h2><div class="gallery">{cards}</div></html>'
    (out/'review.html').write_text(page,encoding='utf-8')
    print(json.dumps(dict(verdict=summary['scientific_verdict'],timing=timing,failed_checks=[x['check_id'] for x in checks if x['status']=='fail']),ensure_ascii=False),flush=True)
    return 2 if summary['scientific_verdict']=='fail' else 0


if __name__ == '__main__':
    sys.exit(main())
