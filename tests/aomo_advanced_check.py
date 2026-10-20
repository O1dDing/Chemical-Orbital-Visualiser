"""Audit intended advanced selections against native records and raw NBO matrices.

This does not launch the application, evaluate GPU fields, or obtain expected
terms from production coefficients. Combine its state checks with GPU readback.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path
from aomo_numerics import np,LFN,reference_field_columns,read_w,tables,sha,save


def lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def flat(term):
    return dict(**term['orbital'],coefficient=term['coefficient']) if 'orbital' in term else term


def key(term):
    return tuple(term[k] for k in ('kind','spin','index'))


def check(package,evidence,expected,only=()):
    package,evidence=Path(package),Path(evidence)
    case=dict(id=package.name,fchk=str(package/'canonical.fchk'),report=str(package/'analysis.log'),
        archive=str(package/('FILE.47' if (package/'FILE.47').exists() else 'electronic.47')),
        matrices={k:str(package/f'FILE.{v}') for k,v in LFN.items()})
    reference=reference_field_columns(case);raw=tables(case['report'])
    rank=len(next(iter(raw['nao'].values())))
    matrices={k:read_w(case['matrices'][k],k,rank,rank) for k in ('NAONHO','NBONLMO')}
    actions=lines(evidence/'actions.jsonl');events=lines(evidence/'events.jsonl')
    selections=sorted((int(p.stem.split('-')[1]),p,json.loads(p.read_text())) for p in evidence.glob('selection-*.json'))
    volumes=[(p,json.loads(p.read_text())) for p in evidence.glob('*.volume.json')]
    results=[];other_checks=[]

    def terms_match(actual,wanted,errors,label,exact=True,tolerance=1e-6):
        actual=[flat(t) for t in actual];wanted=[flat(t) for t in wanted]
        a={key(t):t['coefficient'] for t in actual};w={key(t):t['coefficient'] for t in wanted}
        if len(a)!=len(actual):errors.append(label+': duplicate terms')
        if (set(a)!=set(w) if exact else not set(w)<=set(a)):errors.append(label+': term identities differ')
        err=max((abs(a[k]-v) for k,v in w.items() if k in a),default=0.)
        if err>tolerance:errors.append(f'{label}: signed coefficient error {err}')
        return err

    for intended in expected['checks']:
        trigger=intended['trigger']
        if only and trigger not in only:continue
        target=intended.get('popup_choice',trigger)
        matches=[(i,a) for i,a in enumerate(actions) if a.get('id')==target and a.get('op')=='click']
        # One expectation refers to one click, in chronological order. Repeated
        # context-setting clicks may intentionally have no field readback.
        occurrence=intended.get('occurrence',1)
        matches=matches[occurrence-1:occurrence]
        if intended.get('state_only'):
            state_errors=[]
            if not matches:state_errors.append('Expected state-only click missing')
            for i,a in matches:
                lower=actions[i-1]['frame'] if i else -1
                if a.get('status')!='executed':state_errors.append('State-only click failed')
                if any(lower<f<=a['frame'] for f,p,s in selections):state_errors.append('State-only click unexpectedly applied an orbital')
                subsequent=next((s for f,p,s in selections if f>a['frame']),None)
                if subsequent is None or len(subsequent['terms'])!=intended['selected_nao_term_count']:state_errors.append('Selected term count not demonstrated by following sum')
            results.append(dict(trigger=trigger,passed=not state_errors,errors=state_errors,scope='state-only click plus subsequent applied sum term count'))
            continue
        if 'mode_code' not in intended:
            other_checks.append(dict(trigger=trigger,status='not_covered_by_orbital_state_checker'))
            continue
        if not matches:
            results.append(dict(trigger=trigger,passed=False,errors=['Expected click missing']));continue
        for action_index,action in matches:
            errors=[];prior_frame=actions[action_index-1]['frame'] if action_index else -1
            found=[(f,p,s) for f,p,s in selections if prior_frame<f<=action['frame']]
            if action.get('status')!='executed':errors.append('Click did not execute')
            if len(found)!=1:
                results.append(dict(trigger=trigger,passed=False,errors=errors+[f'Expected one newly applied selection in click window, found {len(found)}']));continue
            frame,path,selected=found[0];terms=[flat(t) for t in selected['terms']]
            if not any(e.get('kind')=='aomo.selection' and e.get('file')==path.name for e in events):errors.append('Missing applied selection event')
            if not selected.get('available'):errors.append('Selection unavailable')
            if selected['mode']!=intended['mode_code']:errors.append('Selection mode differs')
            if selected.get('dataset_id')!=expected.get('integration_id'):errors.append('Dataset identity differs from fixed expectation')
            if 'target_canonical_index' in intended and selected.get('target_canonical_index')!=intended['target_canonical_index']:errors.append('Canonical target differs')
            if selected.get('normalize')!=intended.get('normalize',False):errors.append('Normalization differs')
            if 'term_count' in intended and len(terms)!=intended['term_count']:errors.append('Term count differs')
            if 'terms' in intended:terms_match(terms,intended['terms'],errors,'fixed expectation')
            if 'contains_terms' in intended:terms_match(terms,intended['contains_terms'],errors,'fixed subset',False)
            if 'excluded_parent_index' in intended and any(t['kind']=='NBO' and t['index']==intended['excluded_parent_index'] for t in terms):errors.append('Excluded parent remains in tail')
            for attr,termattr in [('terms_kind','kind'),('terms_spin','spin')]:
                if attr in intended and any(t[termattr]!=intended[attr] for t in terms):errors.append(attr+' differs')

            # Context comes from earlier intended native clicks, never from
            # observed coefficients, label text, or a fitted field sign.
            context={}
            for previous in actions[:action_index+1]:
                m=re.fullmatch(r'nbo\.typed\.(NHO|NLMO):(total|alpha|beta):(\d+)',previous.get('id',''))
                if m and previous.get('status')=='executed':context[m[1]]=(m[2],int(m[3]))
            derived=None;raw_error=0.;raw_parent=None;raw_angular={}
            if trigger.startswith('nbo.nho.') and trigger!='nbo.nho.sum.full':
                if 'NHO' not in context:errors.append('No preceding explicit NHO context')
                else:
                    spin,index=context['NHO'];column=matrices['NAONHO'][spin][:,index]
                    labels=raw['nao'][spin]
                    if 'term_angular' in intended:
                        members=[i for i in range(rank) if labels[i+1]['angular'].startswith(intended['term_angular'])]
                    elif 'matches_actual_nho' in intended:members=list(range(rank))
                    elif 'terms' in intended:members=[t['index'] for t in intended['terms']]
                    else:members=[int(trigger.rsplit('.',1)[1])] if trigger.startswith('nbo.nho.nao.') else []
                    derived=[dict(kind='NAO',spin=spin,index=i,coefficient=float(column[i])) for i in members]
                    raw_error=terms_match(terms,derived,errors,'raw NAONHO',tolerance=2e-8)
                    if 'source_nao_angular' in intended and any(labels[t['index']+1]['angular']!=intended['source_nao_angular'] for t in terms):errors.append('Literal source NAO angular label differs')
            if 'angular_nao_weights' in intended:
                spin,index=context.get('NHO',(None,None))
                if spin is None:errors.append('Missing NHO context for angular weights')
                else:
                    column=matrices['NAONHO'][spin][:,index]
                    for angular,w in intended['angular_nao_weights'].items():
                        raw_angular[angular]=float(sum(column[i]**2 for i in range(rank) if raw['nao'][spin][i+1]['angular'].startswith(angular)))
                        if abs(raw_angular[angular]-w)>1e-6:errors.append('Expected NHO weight disagrees with raw '+angular)
            if trigger.startswith('nbo.nlmo.') or 'verified_parent' in intended:
                if 'NLMO' not in context:errors.append('No preceding explicit NLMO context')
                else:
                    spin,index=context['NLMO'];row=raw['nlmo'][spin].get(index+1)
                    parents=[] if row is None else [i-1 for i,o in raw['nbo'][spin].items() if o['identity']==row['identity']]
                    if len(parents)!=1:errors.append('Printed NLMO parent is not unique')
                    else:
                        raw_parent=parents[0];column=matrices['NBONLMO'][spin][:,index]
                        if 'verified_parent' in intended and key(intended['verified_parent'])!=('NBO',spin,raw_parent):errors.append('Expected parent differs from literal report')
                        if 'printed_parent_percent' in intended and abs(row['percent']-intended['printed_parent_percent'])>1e-4:errors.append('Printed parent percent differs')
                        if trigger in ('nbo.nlmo.main','nbo.nlmo.tail'):
                            members=[raw_parent] if trigger.endswith('.main') else [i for i in range(rank) if i!=raw_parent]
                            derived=[dict(kind='NBO',spin=spin,index=i,coefficient=float(column[i])) for i in members]
                            raw_error=terms_match(terms,derived,errors,'raw NBONLMO',tolerance=2e-8)

            # Independent norm uses actual intended identities, but coefficients
            # must already pass fixed/raw expectations above. No JSON AO columns
            # enter this reference. Overlay norms remain separate per field.
            n=reference['overlap'].shape[0];components=[]
            for t in terms:
                spin=t['spin'];kind=t['kind'];index=t['index']
                if kind=='canonical':
                    if spin=='beta':index-=reference['canonical']['alpha'].shape[1]
                    column=reference['canonical'][spin][:,index]
                elif kind=='GaussianAO':column=np.eye(n)[:,index]
                else:column=reference['families'][kind][spin][:,index]
                components.append(t['coefficient']*column)
                target_mo=selected.get('target_canonical_index')
                if selected['mode'] in (1,2,3,4) and target_mo is not None and kind in reference['links']:
                    local=target_mo-reference['canonical']['alpha'].shape[1] if spin=='beta' else target_mo
                    if abs(t['coefficient']-reference['links'][kind][spin][index,local])>1e-6:errors.append('Component coefficient differs from independent MO link')
            fields=components if selected['mode']==4 else [np.sum(components,axis=0)]
            reconstruction_error=None
            if 'matches_actual_nho' in intended:
                kind,spin,index=intended['matches_actual_nho'].split(':')
                reconstruction_error=float(np.max(abs(fields[0]-reference['families'][kind][spin][:,int(index)])))
                if reconstruction_error>2e-6:errors.append('Complete NAO sum does not reconstruct raw NHO')
            norms=[float(c@reference['overlap']@c) for c in fields]
            if selected['normalize']:norms=[1. if v>1e-20 else 0. for v in norms]
            if len(selected['metric_norm2'])!=len(norms) or max((abs(a-b) for a,b in zip(selected['metric_norm2'],norms)),default=0.)>2e-6:errors.append('Reported norm differs from raw field metric')
            if 'metric_norm2_approx' in intended and (len(norms)!=1 or abs(norms[0]-intended['metric_norm2_approx'])>2e-6):errors.append('Metric norm differs from fixed expectation')
            field_count=intended.get('field_count',1)
            if len(selected.get('columns',[]))!=field_count or len(fields)!=field_count:errors.append('Actual field count differs')
            next_frame=min((f for f,p,s in selections if f>frame),default=float('inf'))
            linked=[(p,v) for p,v in volumes if frame<v['frame']<next_frame]
            if set(v.get('field_index',0) for p,v in linked)!=set(range(field_count)):errors.append('Missing or extra actual volume field indices')
            results.append(dict(trigger=trigger,occurrence=intended.get('occurrence'),selection_file=path.name,selection_frame=frame,action_frame=action['frame'],passed=not errors,errors=errors,raw_matrix_coefficient_max_error=raw_error,raw_parent_index=raw_parent,raw_angular_weights=raw_angular,independent_norm2=norms,full_nho_reconstruction_error=reconstruction_error,volume_files=[p.name for p,v in linked]))

    session=json.loads((evidence/'session.json').read_text());selection_errors=[e for e in events if e.get('kind')=='aomo.selection.error']
    return dict(schema='cov.aomo.advanced.state.audit.v1',passed=bool(results) and all(r['passed'] for r in results) and not selection_errors and session['failed_commands']==0,
        scope='orbital selection state and raw-matrix semantics; GPU signed fields and rendered UI need separate acceptance',
        requested_trigger_subset=list(only),planned_check_count=len(expected['checks']),executed_orbital_state_check_count=len(results),unverified_other_ui_check_count=len(other_checks),source_checks=reference['checks'],source_sha256={k:sha(v) for k,v in dict(NAONHO=case['matrices']['NAONHO'],NBONLMO=case['matrices']['NBONLMO'],fchk=case['fchk'],archive=case['archive'],report=case['report']).items()},
        selection_errors=selection_errors,failed_commands=session['failed_commands'],checks=results,other_ui_checks=other_checks,
        deferred_visual_assertions=expected.get('assertions',[]))


def main():
    p=argparse.ArgumentParser();p.add_argument('package',type=Path);p.add_argument('evidence',type=Path);p.add_argument('--expected',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--only-trigger',action='append',default=[]);a=p.parse_args()
    result=check(a.package,a.evidence,json.loads(a.expected.read_text()),a.only_trigger);save(a.output,result)
    print(json.dumps(dict(passed=result['passed'],checks=len(result['checks']),failed=[r for r in result['checks'] if not r['passed']]),ensure_ascii=False))
    return 0 if result['passed'] else 2


if __name__=='__main__':raise SystemExit(main())
