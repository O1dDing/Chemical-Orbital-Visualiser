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

"""Independent same-input checks of the actual production P/Q matrices.

IOData supplies source coefficients, occupations, and spin layout. The literal
FCHK producer SCF density is parsed separately. This module never changes COV
outputs or borrows producer-MO decimal-print uncertainty for an algorithm test.
"""
from pathlib import Path
import math
import re
import numpy as np

U=np.finfo(np.float64).eps/2
TINY=np.nextafter(np.float64(0),np.float64(1))


def gamma(n):
    if not isinstance(n,int) or n<0 or n*U>=.1:raise ValueError('Invalid arithmetic operation bound')
    return n*U/(1-n*U)


def same_input_density(coefficients,weights):
    c=np.asarray(coefficients,dtype=np.float64);w=np.asarray(weights,dtype=np.float64)
    if c.ndim!=2 or not c.shape[0] or w.shape!=(c.shape[1],) or not np.isfinite(c).all() or not np.isfinite(w).all():
        raise ValueError('Finite AO-by-MO coefficients and one weight per MO required')
    # Zero weights may be retained. A known zero Q is an available matrix.
    k=c.shape[1];nops=12*k+64
    with np.errstate(over='raise',invalid='raise'):
        nominal=(c*w)@c.T
        magnitude=(np.abs(c)*np.abs(w))@np.abs(c).T
        radius=np.nextafter(gamma(nops)*magnitude+nops*TINY,np.inf)
    if not np.isfinite(nominal).all() or not np.isfinite(radius).all():raise ValueError('Overflowed reference')
    return nominal,radius,nops


def compare_same_input(observed,nominal,radius):
    a=np.asarray(observed,dtype=np.float64)
    if a.shape!=nominal.shape or not np.isfinite(a).all():raise ValueError('Invalid measured matrix')
    with np.errstate(over='raise',invalid='raise'):
        budget=np.nextafter(radius+gamma(8)*(np.abs(a)+np.abs(nominal)+radius)+8*TINY,np.inf)
        difference=np.abs(a-nominal);excess=np.maximum(difference-budget,0)
    if not np.isfinite(budget).all():raise ValueError('Overflowed comparison budget')
    return dict(pass_same_parsed_input=bool(np.all(excess==0)),violating_entries=int(np.count_nonzero(excess)),
        maximum_absolute_difference=float(np.max(difference)),maximum_arithmetic_budget=float(np.max(budget)),
        maximum_excess=float(np.max(excess)),source_print_uncertainty_used=False)


def unpack(packed,n):
    values=np.asarray(packed,dtype=np.float64)
    if values.shape!=(n*(n+1)//2,) or not np.isfinite(values).all():raise ValueError('Incomplete or nonfinite packed matrix')
    out=np.empty((n,n));lower=np.tril_indices(n);out[lower]=values;out[(lower[1],lower[0])]=values
    return out


def read_producer_scf_matrices(path,n):
    lines=Path(path).read_text(encoding='ascii').splitlines();result={}
    for name,channel in (('Total SCF Density','total'),('Spin SCF Density','spin')):
        positions=[i for i,line in enumerate(lines) if line[:40].strip()==name]
        if not positions:continue
        if len(positions)!=1:raise ValueError('Duplicate producer SCF density declaration')
        i=positions[0];match=re.fullmatch(r'R\s+N=\s*(\d+)\s*',lines[i][40:].strip())
        if match is None or int(match.group(1))!=n*(n+1)//2:raise ValueError('Invalid producer SCF density dimension')
        tokens=[]
        for line in lines[i+1:]:
            if len(tokens)>=int(match.group(1)):break
            tokens.extend(line.split())
        if len(tokens)!=int(match.group(1)):raise ValueError('Incomplete producer SCF density')
        result[channel]=unpack([float(v.replace('D','E')) for v in tokens],n)
    return result


def review_density(source,molecule,indices,scales,independent_s,production_s,evidence,trace_limit,out):
    checks=[];arrays={}
    def check(name,ok,detail):checks.append(dict(check=name,status='pass' if ok else 'fail',observed=detail))
    n=len(indices)
    check('DENSITY-EVIDENCE-CONTRACT',evidence['schema']==1 and evidence['basis_count']==n and
        evidence['source_format']=='FCHK' and evidence['scalar_type']=='IEEE754-float64' and
        evidence['packing']=='lower triangular: i*(i+1)/2+j, i>=j' and
        evidence['ao_representation']=='COV internal normalized real or Cartesian AO order',
        {'basis_count':n,'scope':'Actual production AO matrices and separate raw-MO reconstructions'})
    ca=np.asarray(molecule.mo.coeffsa);cb=np.asarray(molecule.mo.coeffsb)
    fa=np.asarray(molecule.mo.occsa);fb=np.asarray(molecule.mo.occsb)
    if molecule.mo.kind not in ('restricted','unrestricted'):raise ValueError('Unimplemented independent spin-layout reference')
    shared=molecule.mo.kind=='restricted'
    if shared:
        if not np.array_equal(ca,cb):raise ValueError('Inconsistent shared coefficient reference')
        c=ca;weights={'total':fa+fb,'spin':fa-fb};model='canonical-shared-source-identities'
    else:
        c=np.concatenate((ca,cb),axis=1);weights={'total':np.concatenate((fa,fb)),'spin':np.concatenate((fa,-fb))}
        model='explicit-spin-orbital-occupations'
    expected_count={'total':float(fa.sum()+fb.sum()),'spin':float(fa.sum()-fb.sum())}
    if not all(np.all((v==0)|(v==1)) for v in (fa,fb)):raise ValueError('FCHK canonical occupations are not complete integer spin counts')
    mapped=c[indices]*scales[:,None]
    stored=read_producer_scf_matrices(source,n)
    for channel in ('total','spin'):
        nominal,radius,nops=same_input_density(mapped,weights[channel])
        raw=evidence['production_mo_reconstruction'][channel]
        check('DENSITY-RECONSTRUCTION-MODEL-'+channel,raw['status']=='available' and raw['provenance']=='derived' and
            raw['occupation_model']==model and raw['model_provenance']=='producer' and raw['occupied_space_complete'] and
            raw['occupation_sum']==expected_count[channel] and raw['expected_electron_count']==expected_count[channel],
            {k:v for k,v in raw.items() if k!='packed'})
        reconstructed=unpack(raw['packed'],n)
        comparison=compare_same_input(reconstructed,nominal,radius)
        check('DENSITY-MO-SAME-INPUT-'+channel,comparison['pass_same_parsed_input'],{**comparison,
            'operation_bound':nops,'scope':'Only binary64 arithmetic; zero decimal-print uncertainty'})
        actual_record=evidence['actual'][channel];actual=unpack(actual_record['packed'],n)
        provenance='producer' if channel in stored else 'derived'
        check('DENSITY-ACTUAL-AVAILABILITY-'+channel,actual_record['status']=='available' and
            actual_record['provenance']==provenance and actual_record['occupied_space_complete'],
            {'expected_provenance':provenance,**{k:v for k,v in actual_record.items() if k!='packed'}})
        if channel in stored:
            expected=stored[channel][np.ix_(indices,indices)]*scales[:,None]*scales[None,:]
            check('DENSITY-PRODUCER-AO-TRANSFORM-'+channel,np.array_equal(actual,expected),
                {'maximum_absolute_error':float(np.max(np.abs(actual-expected))),
                 'scope':'Exact parsed producer values under the independently established signed AO permutation'})
        else:
            comparison=compare_same_input(actual,nominal,radius)
            check('DENSITY-ACTUAL-MO-FALLBACK-'+channel,comparison['pass_same_parsed_input'],comparison)
        for label,matrix,info in (('actual',actual,actual_record),('mo',reconstructed,raw)):
            trace=float(np.einsum('ij,ji',matrix,independent_s))
            check('DENSITY-INDEPENDENT-ELECTRON-TRACE-'+label+'-'+channel,
                abs(trace-expected_count[channel])<=trace_limit,
                {'trace':trace,'expected_electrons':expected_count[channel],'absolute_error':abs(trace-expected_count[channel]),
                 'tolerance':trace_limit,'scope':'Same-source density and independent AO metric; a failure may originate in the input'})
            implemented_trace=float(np.einsum('ij,ji',matrix,production_s))
            magnitude=float(np.sum(np.abs(matrix)*np.abs(production_s.T)))
            allowance=gamma(8*n*n+64)*magnitude+(8*n*n+64)*TINY
            stated=info['metric_trace']
            check('DENSITY-REPORTED-TRACE-'+label+'-'+channel,stated is not None and math.isfinite(stated) and
                abs(stated-implemented_trace)<=allowance,
                {'recorded':stated,'independently_recomputed_from_actual_matrix_and_COV_metric':implemented_trace,
                 'arithmetic_allowance':allowance})
        arrays.update({channel+'_actual':actual,channel+'_mo_production':reconstructed,
                       channel+'_mo_independent':nominal,channel+'_arithmetic_radius':radius})
    np.savez_compressed(Path(out)/'independent-density.npz',**arrays)
    return checks
