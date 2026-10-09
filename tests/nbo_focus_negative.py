"""Real-source failure paths for NAO import; no chemistry or source mutation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
from nbo_integration_suite import run, digest, save


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--missing-beta-fchk', type=Path, required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    loaded=json.loads(args.cases.read_text(encoding='utf-8'))
    cases=loaded['cases'] if isinstance(loaded,dict) else loaded
    for case in cases:
        for key in ('fchk','report','archive47','aonbo','nbomo','naomo','aonao','naonbo'):
            if not Path(case[key]).is_absolute():
                case[key]=str((args.cases.resolve().parent/ case[key]).resolve())
    water=next(x for x in cases if x['case_id']=='NBO-01')
    radical=next(x for x in cases if x['case_id']=='NBO-07')
    keys=('fchk','report','archive47','aonbo','nbomo','naomo','aonao','naonbo')
    hashes={str(Path(c[k])):digest(Path(c[k])) for c in cases for k in keys}
    results=[]
    def probe(name, base=water, replacements=None, expected='reject'):
        folder=args.output/name;folder.mkdir()
        files={k:Path(base[k]) for k in keys}
        files.update(replacements or {})
        proc=run(args.build/'cov_nbo_probe.exe', [files[k] for k in keys[:5]]+
                 [folder/'probe']+[files[k] for k in keys[5:]],folder/'process.log')
        source=folder/'probe/dataset.json'
        data=json.loads(source.read_text(encoding='utf-8')) if source.exists() else None
        available=[m for m in (data or {}).get('mo_decompositions',[]) if m['available']]
        if expected=='reject':
            passed=proc['returncode']!=0 or (data is not None and not available)
        elif expected=='unavailable':
            passed=proc['returncode']==0 and data is not None and not available
        elif expected=='incomplete-reject':
            passed=proc['returncode']==2 and data is not None and not available and any(
                v['status']=='incomplete_nao_transforms' for v in data['nao_validation'])
        elif expected=='missing-beta':
            validation={v['spin']:v for v in (data or {}).get('nao_validation',[])}
            passed=proc['returncode']==0 and bool(available) and \
                'alpha' in validation and validation['alpha']['available'] and \
                'beta' in validation and not validation['beta']['direct_fchk_coefficients'] and \
                validation['beta']['available'] and validation['beta']['status']=='verified_archive_only' and \
                not validation['beta']['canonical_occupations'] and \
                not any(m['spin']=='beta' for m in available) and \
                all(m['occupation'] in (0,1) for m in available)
        elif expected=='phase-positive':
            validation={v['spin']:v for v in (data or {}).get('nao_validation',[])}
            passed=proc['returncode']==0 and set(validation)=={'alpha','beta'} and \
                all(v['available'] and v['direct_fchk_coefficients'] for v in validation.values())
        else:
            raise ValueError(expected)
        record=dict(id=name,expected=expected,status='pass' if passed else 'fail',
                    process=proc,available_decompositions=len(available),
                    nao_validation=(data or {}).get('nao_validation',[]))
        results.append(record);save(folder/'check.json',record)
    def variant(name, source, transform):
        path=args.output/name
        path.write_text(transform(Path(source).read_text(encoding='utf-8')),encoding='utf-8')
        return path
    probe('no-NAO-sidecars',replacements={k:'-' for k in ('naomo','aonao','naonbo')},expected='unavailable')
    probe('NAOMO-without-AONAO',replacements={'aonao':'-'},expected='incomplete-reject')
    truncated=variant('truncated-naomo.txt',water['naomo'],lambda s:'\n'.join(s.splitlines()[:10])+'\n')
    probe('truncated-NAOMO',replacements={'naomo':truncated})
    probe('wrong-matrix-kind',replacements={'naomo':Path(water['aonao'])})
    probe('wrong-NAONBO-kind',replacements={'naonbo':Path(water['aonbo'])})
    def perturb(s):
        lines=s.splitlines()
        for i,line in enumerate(lines[3:],3):
            match=re.search(r'[+-]?\d+\.\d+',line)
            if match:
                lines[i]=line[:match.start()]+format(float(match.group())+0.05,'.9f')+line[match.end():]
                break
        return '\n'.join(lines)+'\n'
    damaged=variant('perturbed-naomo.txt',water['naomo'],perturb)
    probe('perturbed-transform',replacements={'naomo':damaged})
    probe('foreign-spin-and-dimension',replacements={'naomo':Path(radical['naomo'])})
    repeated=variant('two-analysis-segments.txt',water['report'],lambda s:s+'\n'+s)
    probe('ambiguous-analysis',replacements={'report':repeated})
    probe('missing-beta-is-not-a-pass',base=radical,replacements={'fchk':args.missing_beta_fchk},expected='missing-beta')
    probe('actual-column-phase-positive',base=radical,expected='phase-positive')
    unchanged=all(digest(Path(p))==h for p,h in hashes.items())
    results.append(dict(id='original-producer-files-unchanged',status='pass' if unchanged else 'fail',hashes=hashes))
    save(args.output/'summary.json',results)
    print(json.dumps(dict(passed=sum(r['status']=='pass' for r in results),total=len(results))))
    return 0 if all(r['status']=='pass' for r in results) else 2


if __name__=='__main__':
    raise SystemExit(main())
