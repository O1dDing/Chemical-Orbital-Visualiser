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

"""Independent numerical review after a complete native collection barrier.

IOData/GBasis define the reference basis and raw orbitals. COV data are measured
outputs only. This subset does not certify molecular physics or all UI semantics.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import sys
import time

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
sys.path.insert(0,r'F:\Dev\cov-validation-20260905\reference-deps')
import numpy as np
from iodata import load_one
from gbasis.wrappers import from_iodata
from gbasis.integrals.overlap import overlap_integral
from gbasis.evals.eval import evaluate_basis
from grid_reference_compare import grid_metrics
from validation_process import atomic_json
from review_density_matrices import review_density
from reference_cache import from_review_environment, array_identity

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))

def internal_basis_map(mol):
    # Internal component order from the documented representation contract;
    # source order/signs come independently from IOData's format conventions.
    cart={0:['1'],1:['x','y','z'],2:['xx','yy','zz','xy','xz','yz'],
          3:['xxx','yyy','zzz','xyy','xxy','xxz','xzz','yzz','yyz','xyz'],
          4:['xxxx','yyyy','zzzz','xxxy','xxxz','xyyy','yyyz','xzzz','yzzz','xxyy','xxzz','yyzz','xxyz','xyyz','xyzz']}
    indices=[];scales=[];offset=0
    for shell in mol.obasis.shells:
        for l,kind in zip(shell.angmoms,shell.kinds):
            source=mol.obasis.conventions[(l,kind)]
            if kind=='p':
                target=['c0']+[name for m in range(1,l+1) for name in (f'c{m}',f's{m}')]
                phases=[1]+[(-1)**m for m in range(1,l+1) for _ in range(2)]
            else:target=cart[l];phases=[1]*len(target)
            for name,phase in zip(target,phases):
                candidates=[i for i,label in enumerate(source) if label.lstrip('-')==name]
                if len(candidates)!=1:raise ValueError('Independent AO convention is ambiguous')
                index=candidates[0]
                indices.append(offset+index);scales.append(phase*(-1 if source[index].startswith('-') else 1))
            offset+=len(source)
    return np.asarray(indices,int),np.asarray(scales,float)

def points_at(meta,ids):
    ids=np.asarray(ids,dtype=np.int64)
    n=np.array([meta[k] for k in ('nx','ny','nz')],dtype=np.int64)
    ijk=np.column_stack((ids%n[0],(ids//n[0])%n[1],ids//(n[0]*n[1]))).astype(np.float32)
    lo=np.asarray(meta['grid_box_bohr'][:3],dtype=np.float32)
    hi=np.asarray(meta['grid_box_bohr'][3:],dtype=np.float32)
    frac=ijk/(n-1).astype(np.float32)
    return (lo.astype(float)+frac.astype(float)*(hi-lo).astype(float)).astype(np.float32).astype(float)

def review(root,case_id,out):
    start=time.perf_counter()
    manifest=load(root/'manifest.json')
    barrier=load(root/'collection-complete.json')
    if not barrier['all_terminal'] or barrier['terminal_cases']!=manifest['case_count']:
        raise RuntimeError('The full frozen collection must finish before numerical review')
    if barrier['round_identity']!=manifest['round_identity']:raise RuntimeError('Collection barrier identity differs')
    record=next(x for x in manifest['cases'] if x['case_id']==case_id)
    terminal=load(root/'cases'/case_id/'terminal.json')
    if terminal['round_identity']!=manifest['round_identity']:raise RuntimeError('Case/round identity differs')
    source=Path(record['input'])
    if sha(source)!=record['sha256']:raise RuntimeError('Input hash mismatch')
    directory=root/terminal['evidence_directory']
    artifact_list=root/terminal['evidence_manifest']
    if sha(artifact_list)!=terminal['evidence_manifest_sha256']:raise RuntimeError('Evidence manifest changed')
    for name,identity in load(artifact_list).items():
        if sha(directory/name)!=identity['sha256']:raise RuntimeError('Evidence changed: '+name)
    if out.exists():raise FileExistsError('Existing scientific review is never overwritten')
    out.mkdir(parents=True)
    production=load(directory/'production.json')
    limits=manifest['thresholds']
    checks=[]
    def check(name,passed,detail):checks.append({'check':name,'status':'pass' if passed else 'fail','observed':detail})
    ordinary=directory/'production-off.json'
    if ordinary.exists():
        check('BUILD-ON-OFF',load(ordinary)==production,{'scope':'Entire recorded production analysis and diagram data'})
    else:checks.append({'check':'BUILD-ON-OFF','status':'insufficient','observed':'Ordinary-build collection absent'})
    events=[json.loads(line) for line in (directory/'native'/'events.jsonl').read_text(encoding='utf-8').splitlines()]
    diagnostics=[event['data'] for event in events if event['kind']=='input.numerical_diagnostics']
    check('ACTUAL-EXE-METRIC-DATA',len(diagnostics)==1 and diagnostics[0]==production['numerical_diagnostics'],
          {'matching_actual_exe_input_events':len(diagnostics)})
    density_events=[event['data'] for event in events if event['kind']=='input.density_evidence']
    check('ACTUAL-EXE-DENSITY-DATA',len(density_events)==1 and density_events[0]==production['density_evidence'],
          {'matching_actual_exe_density_events':len(density_events)})
    mol=load_one(str(source));basis=tuple(from_iodata(mol));coeff=np.asarray(mol.mo.coeffs)
    reference_cache=from_review_environment(record['sha256'], Path(__file__))
    coefficient_identity=array_identity(coeff)
    source_s=(reference_cache.array('overlap', {'basis':'IOData source AO'},
              lambda: overlap_integral(basis)) if reference_cache else overlap_integral(basis))
    indices,scales=internal_basis_map(mol)
    expected_s=source_s[np.ix_(indices,indices)]*scales[:,None]*scales[None,:]
    expected_c=coeff[indices]*scales[:,None]
    actual_c=np.asarray([mo['coefficients'] for mo in production['orbitals']]).T
    raw_c=np.asarray([mo['gaussian_source_coefficients'] for mo in production['orbitals']]).T
    actual_s=np.asarray(production['overlap']).reshape(source_s.shape)
    np.savez_compressed(out/'independent-metric.npz',overlap=source_s,source_coefficients=coeff,
                        source_indices=indices,basis_scales=scales)
    overlap_error=float(np.max(np.abs(actual_s-expected_s)))
    check('AO-INTEGRAL',overlap_error<=limits['basis_overlap_max_absolute_error'],{'max_absolute_error':overlap_error})
    raw_error=float(np.max(np.abs(raw_c-coeff)))
    check('MO-RAW-PRECISION',raw_error<=limits['raw_mo_coefficients_max_absolute_error'],{'max_absolute_error':raw_error})
    mapped_error=float(np.max(np.abs(actual_c-expected_c)))
    check('MO-AO-TRANSFORM',mapped_error<=limits['mapped_mo_coefficients_max_absolute_error'],{'max_absolute_error':mapped_error})
    recorded=np.asarray(production['ao_transform'])
    mapping_ok=np.array_equal(recorded[:,0],np.arange(len(indices))) and np.array_equal(recorded[:,1],indices) and \
               np.array_equal(recorded[:,2],scales) and np.array_equal(recorded[:,3],scales)
    check('AO-TRANSFORM-PROVENANCE',mapping_ok,{'ao_count':len(indices)})
    checks.extend(review_density(source,mol,indices,scales,expected_s,actual_s,production['density_evidence'],
                                limits['electron_trace_error'],out))
    energy_error=float(np.max(np.abs(np.array([mo['energy_hartree'] for mo in production['orbitals']])-mol.mo.energies)))
    check('MO-ENERGIES',energy_error<=limits['energy_error_hartree'],{'max_absolute_error_hartree':energy_error})
    blocks=[]
    for spin in sorted({mo['spin'] for mo in production['orbitals']}):
        columns=[i for i,mo in enumerate(production['orbitals']) if mo['spin']==spin]
        actual=actual_c[:,columns].T@actual_s@actual_c[:,columns]
        reference=coeff[:,columns].T@source_s@coeff[:,columns]
        result={'spin':spin,'ao_count':len(indices),'mo_count':len(columns),
                'reference_gram_max_error':float(np.max(np.abs(reference-np.eye(len(columns))))),
                'production_gram_max_error':float(np.max(np.abs(actual-np.eye(len(columns))))),
                'production_reference_gram_max_error':float(np.max(np.abs(actual-reference)))}
        blocks.append(result)
        check('MO-GRAM-'+str(spin),max(result['reference_gram_max_error'],result['production_gram_max_error'])<=limits['mo_orthonormality_max_error'],result)
    occup=np.array([mo['occupation'] for mo in production['orbitals']])
    trace=float(np.sum(occup*np.einsum('ij,ij->j',actual_c,actual_s@actual_c)))
    check('MO-RECONSTRUCTED-ELECTRON-TRACE',abs(trace-float(mol.nelec))<=limits['electron_trace_error'],
          {'trace':trace,'expected_electrons':float(mol.nelec),'scope':'MO-reconstructed density, not all optional stored producer density fields'})
    distributions=[]
    for i,item in enumerate(production['chemistry']):
        if not item['available']:continue
        for name in ('bonding','channel','manifold'):
            values=np.asarray(item[name],float)
            if item.get(name+'_status')==3 and np.all(values==0):continue
            if not np.isfinite(values).all() or values.min()<-1e-8 or values.max()>1+1e-8 or abs(values.sum()-1)>2e-6:
                distributions.append({'mo':i,'distribution':name,'values':values.tolist()})
    check('DECOMPOSITION-CONSERVATION',not distributions,{'invalid_distributions':distributions,
          'available_mos':sum(item['available'] for item in production['chemistry']),
          'scope':'Every reported available distribution; unavailable chemistry is not certified'})
    atomic_json(out/'metric-checks.json',checks)
    thresholds={'nrms_max':limits['all_mo_nrms_max'],'abs_cosine_min':limits['all_mo_abs_cosine_min'],
                'relative_peak_error_max':limits['all_mo_relative_peak_error_max']}
    cache={};samples=[];full_groups={}
    frames={row['frame']:row for row in (json.loads(line) for line in (directory/'native'/'frames.jsonl').read_text().splitlines())}
    for path in sorted((directory/'native').glob('*.volume.json')):
        meta=load(path);mo=meta['rendered_mo'];values=np.asarray(meta['samples'])
        ids=values[:,0].astype(int)
        key=(tuple(meta['grid_box_bohr']),meta['nx'],meta['ny'],meta['nz'],tuple(ids))
        state=frames[meta['frame']]
        matched=state['rendered_mo']==mo and state['rendered_generation']==meta['generation'] and state['drawn_ui_mo']==mo
        if key not in cache:
            points=points_at(meta,ids)
            spec={'coordinates':array_identity(points),'coefficients':coefficient_identity}
            cache[key]=(reference_cache.array('sampled_mo',spec,
                        lambda: coeff.T@evaluate_basis(basis,points)) if reference_cache
                        else coeff.T@evaluate_basis(basis,points))
        metric=grid_metrics(values[:,1],cache[key][mo],thresholds)
        metric.update(mo=mo,metadata=path.name,frame_association=matched,full_grid='full_grid' in meta)
        if not matched:metric['pass']=False;metric['status']='fail'
        if 'full_grid' not in meta and len(ids)<limits['same_grid_minimum_unique_samples']:
            metric['pass']=False;metric['status']='insufficient'
        samples.append(metric)
        if 'full_grid' in meta:
            full_key=(tuple(meta['grid_box_bohr']),meta['nx'],meta['ny'],meta['nz'])
            full_groups.setdefault(full_key,[]).append((path,meta))
    atomic_json(out/'sampled-textures.json',samples)
    full_results=[]
    for group_index,items in enumerate(full_groups.values()):
        meta=items[0][1];count=meta['nx']*meta['ny']*meta['nz']
        mos=[m['rendered_mo'] for _,m in items]
        reference_path=out/f'full-reference-{group_index:02d}.npy'
        def write_full_reference(destination):
            reference=np.lib.format.open_memmap(destination,mode='w+',dtype='<f8',shape=(len(mos),count))
            for first in range(0,count,4096):
                last=min(first+4096,count)
                reference[:,first:last]=coeff[:,mos].T@evaluate_basis(basis,points_at(meta,np.arange(first,last)))
            reference.flush()
            del reference
        if reference_cache:
            spec={key:meta[key] for key in ('grid_box_bohr','nx','ny','nz')}
            spec.update(mos=mos,coefficients=coefficient_identity,block_size=4096,dtype='<f8')
            stored=reference_cache.materialize('full_mo',spec,write_full_reference)
            # Every review remains portable, independently hashed and auditable.
            shutil.copyfile(stored,reference_path)
        else:
            write_full_reference(reference_path)
        reference=np.load(reference_path,mmap_mode='r',allow_pickle=False)
        for row,(path,metadata) in enumerate(items):
            binary=path.parent/metadata['full_grid']['file']
            if binary.parent.resolve()!=path.parent.resolve() or binary.stat().st_size!=count*4:
                raise ValueError('Full texture path or byte count mismatch')
            actual=np.memmap(binary,dtype='<f4',mode='r',shape=(count,))
            wanted=np.asarray(metadata['samples'])
            if not np.array_equal(actual[wanted[:,0].astype(int)],wanted[:,1]):raise ValueError('Full/sparse texture values differ')
            result=grid_metrics(actual,reference[row],thresholds)
            prefix=path.name.removesuffix('.volume.json')
            ui=load(path.parent/(prefix+'.ui.json'))
            state=frames[metadata['frame']]
            result.update(mo=metadata['rendered_mo'],metadata=path.name,resolution=[meta[k] for k in ('nx','ny','nz')],
                          texture_sha256=sha(binary),reference_file=reference_path.name,reference_row=row,
                          same_frame_ui=ui['state']==state)
            if not result['same_frame_ui']:result['pass']=False;result['status']='fail'
            full_results.append(result)
        del reference
    atomic_json(out/'complete-textures.json',full_results)
    sampled=[x for x in samples if not x['full_grid']]
    expected=load(directory/'expected-evidence.json')
    full_expected=len(expected.get('full_grid_frontiers_zero_based',[]))*len(expected.get('full_grid_resolutions',[]))
    check('ALL-MO-SAMPLED-TEXTURES',len(sampled)==coeff.shape[1] and all(x['pass'] for x in sampled),
          {'expected_mos':coeff.shape[1],'measured':len(sampled),'passed':sum(x['pass'] for x in sampled),
           'maximum_nrms':max((x['nrms'] for x in sampled if 'nrms' in x),default=None)})
    check('FRONTIER-COMPLETE-TEXTURES',len(full_results)==full_expected and all(x['pass'] for x in full_results),
          {'expected':full_expected,'measured':len(full_results),'passed':sum(x['pass'] for x in full_results),
           'maximum_nrms':max((x['nrms'] for x in full_results if 'nrms' in x),default=None)})
    if sha(source)!=record['sha256']:raise RuntimeError('Input changed during review')
    atomic_json(out/'reference-cache.json',{'enabled':reference_cache is not None,
                'events':reference_cache.events if reference_cache else [],
                'cached_verdicts':False,'actual_evidence_recompared':True})
    result={'schema':1,'case_id':case_id,'round_identity':manifest['round_identity'],
            'status':'numeric_subset_pass' if all(x['status']=='pass' for x in checks) else 'numeric_subset_fail',
            'checks':checks,'wall_seconds':time.perf_counter()-start,'formal_case_pass':False,
            'coordinate_reference':'Recorded CUDA float grid interpolation with one final FMA rounding, as in the pre-existing independent checker',
            'remaining_scope':['reference physical conclusions paused','full independent stored-density checks',
                               'symmetry/topology/chemical adjudication','visual correctness and all UI/export semantics',
                               'all required invariants and ordinary-build behavior'],
            'source_sha256':record['sha256'],'reviewer_sha256':sha(Path(__file__)),
            'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in out.iterdir() if p.is_file()}}
    atomic_json(out/'review.json',result)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--round',type=Path,required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=review(args.round,args.case,args.output)
    print(json.dumps({k:result[k] for k in ('case_id','status','wall_seconds','formal_case_pass')}),flush=True)
    raise SystemExit(0 if result['status']=='numeric_subset_pass' else 2)
