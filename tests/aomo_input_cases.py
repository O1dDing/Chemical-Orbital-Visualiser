"""Prepare and run visible input/fallback cases on copied vendor inputs."""
import argparse
import json
import shutil
import subprocess
from pathlib import Path


def prepare(root):
    source=root/'packages'/'NBO-01'
    target=root/'input-packages'
    target.mkdir(exist_ok=True)
    definitions={
        'gaussian-only': ['canonical.fchk'],
        'report-only': ['canonical.fchk','analysis.log'],
        'missing-nao': [p.name for p in source.iterdir() if p.name not in
                       ('FILE.51','FILE.52','FILE.53','FILE.57','FILE.58')],
        'bad-nho': [p.name for p in source.iterdir()],
        'mismatch': ['canonical.fchk','analysis.log','electronic.47'],
    }
    for name,files in definitions.items():
        folder=target/name
        folder.mkdir(exist_ok=False)
        for file in files: shutil.copy2(source/file,folder/file)
        if name=='mismatch':
            shutil.copy2(root/'packages'/'NBO-03'/'electronic.47',folder/'electronic.47')
        if name=='bad-nho':
            path=folder/'FILE.54';lines=path.read_text().splitlines()
            words=lines[3].split();words[0]='9.99999999E+09';lines[3]=' '.join(words)
            path.write_text('\n'.join(lines)+'\n')
    multi=target/'multiple';multi.mkdir(exist_ok=False)
    for name in ('NBO-01','NBO-03'):shutil.copytree(root/'packages'/name,multi/name)
    definitions['multiple']=[]
    plans=root/'input-plans';plans.mkdir(exist_ok=True)
    for name in definitions:
        lines=['COV_VALIDATION 1','scene 0.35 0.65 0.35 1.7 0.03 64',
               f'drop "{source.as_posix()}"','click "aomo.node.NAO:total:4"',
               f'drop "{(target/name).as_posix()}"','hover "scene.viewport"',
               'capture "view"','hover "panel.nbo"','capture "input-state"']
        (plans/f'{name}.plan').write_text('\n'.join(lines)+'\n')
    # Multiple simultaneous paths enter the same production GLFW drop handler.
    lines=['COV_VALIDATION 1','scene 0.35 0.65 0.35 1.7 0.03 64',
           'drop '+' '.join('"'+p.as_posix()+'"' for p in sorted(source.iterdir())),
           'click "aomo.node.NAO:total:4"','volume "nao" "0"',
           'hover "scene.viewport"','capture "view"']
    (plans/'multi-file.plan').write_text('\n'.join(lines)+'\n')


def run(root,exe,attempt):
    results=[];out=root/attempt;out.mkdir(exist_ok=False)
    for plan in sorted((root/'input-plans').glob('*.plan')):
        evidence=out/plan.stem
        result=subprocess.run([str(exe),str(root/'packages'/'NBO-01'/'canonical.fchk'),
             '--validation-plan',str(plan),'--validation-output',str(evidence),
             '--validation-background'],capture_output=True,text=True,timeout=120)
        (out/(plan.stem+'.log')).write_text(result.stdout+result.stderr)
        session=json.loads((evidence/'session.json').read_text()) if (evidence/'session.json').exists() else {}
        results.append(dict(case=plan.stem,returncode=result.returncode,session=session))
        print(plan.stem,result.returncode,session.get('failed_commands'),flush=True)
    (out/'summary.json').write_text(json.dumps(results,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path)
    p.add_argument('--prepare',action='store_true');p.add_argument('--exe',type=Path)
    p.add_argument('--attempt',default='inputs-r01');a=p.parse_args()
    if a.prepare:prepare(a.root.resolve())
    if a.exe:run(a.root.resolve(),a.exe.resolve(),a.attempt)
