from pathlib import Path
import argparse,json,os,sys,hashlib,shutil
ap=argparse.ArgumentParser(description='Replay NBO load failure, reattachment, selection and export through actual native controls.')
ap.add_argument('dest',type=Path)
ap.add_argument('--exe',type=Path,required=True)
ap.add_argument('--water-package',type=Path,required=True,help='Saved H2O package with analysis.log, canonical.fchk, electronic.47 and complete FILE matrices')
ap.add_argument('--o2-package',type=Path,required=True,help='Saved restricted-open-shell O2 package with canonical.fchk and drop.covnbopkg')
args=ap.parse_args()
from validation_process import run_tree
dest=args.dest.resolve();dest.mkdir(parents=True,exist_ok=False)
fixtures=dest/'fixtures';fixtures.mkdir()
water=args.water_package.resolve()
o2=args.o2_package.resolve()
bad=fixtures/'bad.fchk';bad.write_text('Broken fixture\nSP\nNumber of atoms                            I                2\nNumber of basis functions                  I                4\n',encoding='ascii')
unknown=fixtures/'unsupported.input';unknown.write_text('not a wavefunction\n',encoding='ascii')
ambiguous=fixtures/'ambiguous';ambiguous.mkdir()
for name in ['analysis.log','electronic.47','FILE.37','FILE.49','FILE.50','FILE.51','FILE.52','FILE.53']:
    shutil.copyfile(water/name,ambiguous/name)
(ambiguous/'other.log').write_bytes((water/'analysis.log').read_bytes()+b'\n# deliberate discovery duplicate with distinct payload\n')
partial=fixtures/'report-only';partial.mkdir();shutil.copyfile(water/'analysis.log',partial/'analysis.log')
prefix=['COV_VALIDATION 1','window 2100 1250','scene 0.26 0.65 0.35 1.7 0.05 96']
def load(p):return [f'drop "{p.as_posix()}"','wait "load"']
def capture(n):return ['hover "scene.viewport"',f'capture "{n}"']
def export(n,entry='nbo.export'):
    plan=[f'export-name "{n}"',f'seek "{entry}"',f'click "{entry}"']
    if entry=='diagram.export':
        plan += ['seek "diagram.export_options"','click "diagram.export_options"','seek "diagram.export_data"','click "diagram.export_data"','click "diagram.export_options"']
    return plan+capture(n)
def manual():
    p=[]
    for key,name in [('path','analysis.log'),('archive47','electronic.47'),('aonbo','FILE.37'),('nbomo','FILE.49'),('naomo','FILE.51'),('aonao','FILE.52'),('naonbo','FILE.53')]:
        p += [f'seek "nbo.{key}"',f'text "nbo.{key}" "{(water/name).as_posix()}"']
    return p+['seek "nbo.attach"','click "nbo.attach"','wait "attach"']
plans={}
p=prefix+load(o2/'drop.covnbopkg')
p+=['seek "nbo.orbital.pick"','click "nbo.orbital.pick"','click "nbo.typed.NAO:alpha:3"']+export('before')
p+=['seek "nbo.set.canonical"','click "nbo.set.canonical"']+capture('returned-canonical')
p+=load(bad)+export('bad-fchk','diagram.export')
p+=load(o2/'drop.covnbopkg')+export('restored')
p+=['seek "nbo.orbital.pick"','click "nbo.orbital.pick"','click "nbo.typed.NAO:alpha:3"']
p+=load(unknown)+export('bad-unknown','diagram.export')
p+=load(o2/'drop.covnbopkg')+export('restored-again')
plans['failure-recovery']=p
p=prefix+load(water)+load(ambiguous)+['seek "nbo.candidate.0"','click "nbo.candidate.0"']+export('candidate-selected')
p+=['seek "nbo.inputs.advanced"','click "nbo.inputs.advanced"']+manual()+export('manual-attached')
p+=['seek "nbo.overlay.wiberg"','click "nbo.overlay.wiberg"','seek "nbo.overlay.e2"','click "nbo.overlay.e2"']+capture('overlay-on')
p+=load(partial)+export('report-only')+load(water)+export('full-restored')
plans['manual-overlay']=p
checks=[]
def check(name,ok,observed=None):checks.append(dict(name=name,passed=bool(ok),observed=observed))
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def field(ui,key):return next((x['data']['value'] for x in ui['draw_trace'] if x['kind']=='draw.text' and x['data'].get('label')==key),None)
for name,plan in plans.items():
    p=dest/(name+'.plan');p.write_text('\n'.join(plan)+'\n',encoding='utf-8')
    out=dest/name
    result=run_tree([str(args.exe),str(o2/'canonical.fchk'),'--validation-plan',str(p),'--validation-output',str(out),'--validation-background'],dest,dict(os.environ,OMP_NUM_THREADS='3',MKL_NUM_THREADS='3',OPENBLAS_NUM_THREADS='1'),cpu_mask=7,memory_gib=20,timeout_seconds=180)
    check(name+':process',result['exit_code']==0,result)
    if (out/'session.json').exists():check(name+':native-actions',read(out/'session.json')['failed_commands']==0)
    if result['exit_code']:continue
    for path in out.glob('*.active-view.json'):
        label=path.name.removesuffix('.active-view.json');view=read(path);ui=read(out/(label+'.ui.json'));route=read(out/(label+'.analysis.json'));integration=read(out/(label+'.integration.json'))
        check(label+':view-equals-scene',view==json.loads(field(ui,'scene.active_view')))
        check(label+':scene-current',ui['state']['scene_matches_applied'])
        check(label+':route-identity',route['integration_id']==integration.get('id',''))
        if label.startswith('bad-'):
            check(label+':canonical',view['kind']=='canonical' and view['selection'] is None)
            check(label+':no-nbo-analysis',not route['integration_id'] and all(x['provider']=='legacy' for x in route['mo_composition']))
            check(label+':no-salc',read(out/(label+'.salc.json')).get('status')=='not_analysed')
        if label.startswith('restored'):
            check(label+':charge-zero',route['total_atomic_charge']['values']==[0,0])
            check(label+':success-status-refreshed',str(field(ui,'status.detail')).endswith('analysis.log'))
        if label=='candidate-selected':check(label+':two-candidates',field(ui,'nbo.input.candidate_count')=='2')
        if label=='manual-attached':
            check(label+':no-old-candidates',field(ui,'nbo.input.candidate_count')=='0')
            check(label+':actual-source',Path(integration['dataset']['source']['path']).resolve()==(water/'analysis.log').resolve())
        if label in ['report-only','full-restored']:
            for key in ['wiberg','e2']:check(label+':'+key+'-reset',field(ui,'nbo.overlay.'+key+'.checked')=='false')
    if name=='failure-recovery':
        ui=read(out/'returned-canonical.ui.json');active=json.loads(field(ui,'scene.active_view'))
        check('return:status-refreshed',active['kind']=='canonical' and active['selection'] is None and
              field(ui,'status.detail')==active['label'],dict(status=field(ui,'status.detail'),active=active))
    if name=='manual-overlay':
        ui=read(out/'overlay-on.ui.json')
        for key in ['wiberg','e2']:check('overlay-on:'+key,field(ui,'nbo.overlay.'+key+'.checked')=='true')
report=dict(passed=all(x['passed'] for x in checks),binary=str(args.exe),sha256=hashlib.sha256(args.exe.read_bytes()).hexdigest(),checks=checks)
(dest/'checks.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(dict(passed=report['passed'],checks=len(checks),failed=[x['name'] for x in checks if not x['passed']])))
raise SystemExit(0 if report['passed'] else 2)
