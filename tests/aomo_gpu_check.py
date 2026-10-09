"""Compare actual selected GPU fields with independent FCHK/vendor-matrix fields.

The selection record supplies identities and signed weights, never reference AO
coefficients. IOData and GBasis independently evaluate the original FCHK basis.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '2'
from nbo_reference import np, evaluate_basis
from aomo_numerics import reference_field_columns, LFN


def check(package: Path, evidence: Path, expected=None):
    case=dict(id=package.name,fchk=str(package/'canonical.fchk'),
        report=str(package/'analysis.log'),
        archive=str(package/('FILE.47' if (package/'FILE.47').exists() else 'electronic.47')),
        matrices={kind:str(package/f'FILE.{number}') for kind,number in LFN.items()
                  if (package/f'FILE.{number}').exists()})
    reference=reference_field_columns(case)
    mol=reference['molecule'];basis=reference['basis']
    n = mol.obasis.nbasis
    columns = {}
    for spin in reference['spins']:
        columns[spin]={'canonical':reference['canonical'][spin], 'GaussianAO':np.eye(n)}
        for kind,spins in reference['families'].items():columns[spin][kind]=spins[spin]
    alpha_count=next(iter(reference['canonical'].values())).shape[1]
    selections = sorted((int(p.stem.split('-')[1]), json.loads(p.read_text(encoding='utf-8')))
                        for p in evidence.glob('selection-*.json'))
    results = []
    grid_cache={}
    for path in sorted(evidence.glob('*.volume.json')):
        actual = json.loads(path.read_text(encoding='utf-8'))
        prior = [x for frame, x in selections if frame < actual['frame']]
        if not prior:
            results.append(dict(file=path.name, passed=False, reason='no applied typed selection'))
            continue
        selected = prior[-1]
        terms = selected['terms']
        intended=(expected or {}).get(path.name.removesuffix('.volume.json'))
        request_matches=intended is None or (len(terms)==1 and terms[0]['orbital']['kind']==intended['kind']
                                             and terms[0]['orbital']['index']==intended['index'])
        if selected['mode'] == 4:
            terms = [terms[actual.get('field_index', 0)]]
        reference_coeff = np.zeros(n)
        weight_errors = []
        for term in terms:
            ref = term['orbital']; spin = ref['spin']
            if spin not in columns and len(columns)==1:spin=next(iter(columns))
            index = ref['index']
            if ref['kind'] == 'canonical' and spin == 'beta':
                index -= alpha_count
            reference_coeff += term['coefficient']*columns[spin][ref['kind']][:, index]
            target = selected.get('target_canonical_index')
            if target is not None and selected['mode'] in (1, 2, 3, 4) and ref['kind'] in reference['links']:
                target_local = target-alpha_count if spin == 'beta' else target
                matrix = reference['links'][ref['kind']][spin]
                weight_errors.append(abs(term['coefficient']-float(matrix[index, target_local])))
        if selected['normalize']:
            metric=float(reference_coeff @ reference['overlap'] @ reference_coeff)
            if metric<=1e-20:raise ValueError('Cannot normalize a zero reference field')
            reference_coeff/=np.sqrt(metric)
        samples = np.asarray(actual['samples'])
        ids = samples[:, 0].astype(int)
        dims = np.array([actual['nx'], actual['ny'], actual['nz']], dtype=np.float32)
        ijk = np.column_stack((ids % actual['nx'], (ids//actual['nx']) % actual['ny'],
                               ids//(actual['nx']*actual['ny']))).astype(np.float32)
        lo = np.asarray(actual['grid_box_bohr'][:3], dtype=np.float32)
        hi = np.asarray(actual['grid_box_bohr'][3:], dtype=np.float32)
        points = (lo.astype(float)+(ijk/(dims-1)).astype(float)*(hi-lo).astype(float)).astype(np.float32).astype(float)
        grid_key=(tuple(dims),tuple(lo),tuple(hi),tuple(ids))
        if grid_key not in grid_cache:grid_cache[grid_key]=evaluate_basis(basis, points)
        ref_values = reference_coeff @ grid_cache[grid_key]
        difference = samples[:, 1]-ref_values
        nrms = float(np.linalg.norm(difference)/max(np.linalg.norm(ref_values), 1e-30))
        relative_max = float(np.max(np.abs(difference))/max(np.max(np.abs(ref_values)), 1e-30))
        results.append(dict(file=path.name, terms=terms, mode=selected['mode'], samples=len(ids),
                            signed_nrms=nrms, relative_max=relative_max,
                            raw_weight_max_error=max(weight_errors, default=0.0),
                            requested_identity=intended, request_matches=request_matches,
                            passed=nrms < 1e-4 and relative_max < 1e-3 and request_matches
                                and max(weight_errors, default=0.0)<1e-6))
    errors = [json.loads(line) for line in (evidence/'events.jsonl').read_text(encoding='utf-8').splitlines()
              if json.loads(line).get('kind') == 'aomo.selection.error']
    session = json.loads((evidence/'session.json').read_text(encoding='utf-8'))
    report = dict(passed=bool(results) and all(x['passed'] for x in results)
                    and not errors and session['failed_commands'] == 0,
                  numerical_reference='original FCHK canonical + raw AONAO/AONBO/AONHO/AONLMO/AOPNAO; independent LABEL/CONTRACT AO map; GBasis',
                  source_checks=reference['checks'],
                  phase_policy='signed comparison; no fitted global sign allowed',
                  selection_errors=errors, failed_commands=session['failed_commands'], fields=results)
    (evidence/'independent-gpu.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('package', type=Path)
    parser.add_argument('evidence', type=Path)
    parser.add_argument('--expected', type=Path)
    args = parser.parse_args()
    result = check(args.package, args.evidence,
                   json.loads(args.expected.read_text(encoding='utf-8')) if args.expected else None)
    print(json.dumps(result, ensure_ascii=False))
    raise SystemExit(0 if result['passed'] else 2)
