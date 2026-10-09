"""Check visible fallback state after replacing an already-selected NBO dataset."""
import argparse
import json
from pathlib import Path


def check(root):
    results=[]
    expected={
        'gaussian-only': ('canonical',False,{}),
        'report-only': ('canonical',False,{'report':'available','source_association':'missing','aomo':'missing'}),
        'mismatch': ('canonical',False,{'source_association':'rejected','aomo':'missing'}),
        'missing-nao': ('canonical',False,{'source_association':'available','nao':'missing','aomo':'missing','nbo':'available','nlmo':'available','interactions':'available'}),
        'bad-nho': ('canonical',True,{'nho':'rejected','aomo':'available','nbo':'available','nlmo':'available','charges':'available','wiberg':'available'}),
        'multiple': ('NAO',True,{}),
        'multi-file': ('NAO',True,{'source_association':'available','aomo':'available','nho':'available'}),
    }
    for name,(orbital,graph,capabilities) in expected.items():
        path=root/name
        view=json.loads((path/'view.ui.json').read_text(encoding='utf-8'))
        ui_path=path/'input-state.ui.json'
        ui=json.loads(ui_path.read_text(encoding='utf-8')) if ui_path.exists() else view
        state=view['state'];session=json.loads((path/'session.json').read_text(encoding='utf-8'))
        targets={x['id']:x for x in ui['targets']}
        fields={x['data']['label']:x['data']['value'] for x in ui['draw_trace'] if x['kind']=='draw.text'}
        integrations=sorted(path.glob('integration-*.json'),key=lambda p:int(p.stem.split('-')[1]))
        integration=json.loads(integrations[-1].read_text(encoding='utf-8')) if integrations else {}
        caps={x['key']:x['state'] for x in integration.get('capabilities',[])}
        assertions={
            'native_commands':session['failed_commands']==0,
            'rendered_kind':state['rendered_set']==orbital,
            'scene_matches_selection':state['scene_matches_applied'],
            'overall_graph':any(x['id'].startswith('aomo.node.') for x in view['targets'])==graph,
            'capability_isolation':all(caps.get(k)==v for k,v in capabilities.items()),
        }
        if orbital=='canonical':
            assertions['original_fchk_coefficients']=state['rendered_direct_fchk_coefficients'] is True
            assertions['no_stale_nao_label']='NAO 5' not in fields.get('scene.orbital_label','')
        if name=='multiple':
            assertions['explicit_candidate_choice']=all(f'nbo.candidate.{i}' in targets for i in (0,1))
            assertions['prior_view_retained']=state['rendered_source_index']==4
        if name=='gaussian-only':
            assertions['visible_missing_message']='missing' in fields.get('nbo.input.status','').lower()
        if name=='mismatch':
            assertions['visible_mismatch_message']='incompatible_dimensions' in fields.get('nbo.input.status','')
        results.append(dict(case=name,passed=all(assertions.values()),assertions=assertions,
                            visible_fields=fields,expected_capabilities=capabilities))
    report=dict(passed=all(x['passed'] for x in results),cases=results,
        scope='native production drop handler and drawn UI state; OS drag gesture is a separate ordinary-desktop check')
    (root/'semantic-check.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);a=p.parse_args()
    r=check(a.evidence);print(json.dumps(r,ensure_ascii=False));raise SystemExit(0 if r['passed'] else 2)
