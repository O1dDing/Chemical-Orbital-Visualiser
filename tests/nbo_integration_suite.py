"""Scoped NBO acceptance through production parsing and real native UI input.

Every invocation retains a new evidence directory. Producer calculations are
never started here. GPU comparisons use independent IOData/GBasis AO values
and raw producer NBOMO matrices, not production render coefficients.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '2'
from nbo_reference import load_reference, np


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(executable, arguments, logfile, timeout=600):
    start = time.perf_counter()
    argv = [x.as_posix() if isinstance(x, Path) else str(x) for x in arguments]
    proc = subprocess.run([str(executable), *argv], capture_output=True,
                          encoding='utf-8', errors='replace', timeout=timeout)
    logfile.write_text(proc.stdout + proc.stderr, encoding='utf-8')
    return dict(returncode=proc.returncode, seconds=time.perf_counter()-start)


def command(op, name, value=None):
    return op+' '+json.dumps(str(name), ensure_ascii=False) + (
        ' '+json.dumps(str(value), ensure_ascii=False) if value is not None else '')


def selections(dataset):
    """Cover each spin and occupied/virtual and local orbital kinds from data."""
    chosen = []
    for spin in dict.fromkeys(o['spin'] for o in dataset['orbitals']):
        group = [(i, o) for i, o in enumerate(dataset['orbitals']) if o['spin'] == spin]
        for kind in dict.fromkeys(o['kind'] for _, o in group):
            candidates = [(i, o) for i, o in group if o['kind'] == kind]
            chosen.append(max(candidates, key=lambda x: x[1]['occupation'])[0])
        chosen.extend([min(group, key=lambda x: x[1]['occupation'])[0], group[-1][0]])
    return list(dict.fromkeys(chosen))


def case_check(case, build, out, native=True, languages=False):
    out.mkdir()
    checks = []
    def check(name, passed, observed, tolerance=None):
        checks.append(dict(id=name, status='pass' if passed else 'fail',
                           observed=observed, tolerance=tolerance))
    reference = load_reference(case)
    a = reference.analysis_dir
    evidence = reference.matrix_checks()
    save(out/'reference-matrices.json', evidence)
    # A missing beta FCHK column is an explicitly scoped unknown, not a pass.
    unknown = {k: v for k, v in evidence['checks'].items() if v['status'] == 'unknown'}
    permitted_unknown = (reference.archive_header['open_shell'] and
                         set(unknown) <= {'beta_fchk_vs_archive_mo'} and
                         evidence['checks'].get('fchk_spin_density', {}).get('status') == 'pass' and
                         evidence['checks'].get('fchk_scf_density', {}).get('status') == 'pass')
    check('independent_representation', evidence['archive_representation_status'] == 'pass'
          and (not unknown or permitted_unknown), evidence)
    producer = run(build/'cov_nbo_probe.exe', [reference.fchk_path, a/'analysis.log',
        a/'FILE.47', a/'FILE.37', a/'FILE.49', out/'probe'], out/'probe.log')
    check('production_attach', producer['returncode'] == 0, producer)
    if producer['returncode']:
        return checks
    dataset = json.loads((out/'probe/dataset.json').read_text(encoding='utf-8'))
    contract = json.loads((out/'probe/render-contract.json').read_text(encoding='utf-8'))
    expected_count = reference.nbasis * len(reference.spins)
    check('orbital_count_and_energy_identity', len(dataset['orbitals']) == expected_count
          and len(contract['orbitals']) == expected_count
          and all(o['canonical_energy'] is None for o in contract['orbitals']),
          dict(nbo=len(dataset['orbitals']), expected=expected_count,
               canonical_energies=[o['canonical_energy'] for o in contract['orbitals']]))
    for spin in reference.spins:
        tag = spin if len(reference.spins) == 2 else 'total'
        orbs = [o for o in dataset['orbitals'] if o['spin'] == tag]
        ne = float(reference.mo_occupations[spin].sum())
        check('occupation_trace_'+spin, abs(sum(o['occupation'] for o in orbs)-ne) <= 3e-4,
              dict(observed=sum(o['occupation'] for o in orbs), expected=ne), 3e-4)
        for kind, expected in [('AONBO', reference.aonbo[spin]), ('NBOMO', reference.nbomo[spin])]:
            matrices = [m for m in dataset['matrices'] if m['kind'] == kind and m['spin'] == tag]
            error = float(np.max(np.abs(np.asarray(matrices[0]['values']).reshape(expected.shape)-expected))) if len(matrices) == 1 else None
            check('raw_matrix_'+kind+'_'+spin, error is not None and error <= 1e-12, error, 1e-12)
    raw_lines = (a/'analysis.log').read_text(errors='replace').splitlines()
    provenance_errors = []
    for collection in ('populations', 'naos', 'orbitals', 'e2', 'wiberg'):
        for row in dataset[collection]:
            src = row['source']; line = src['line_begin']
            if not 0 < line <= len(raw_lines) or raw_lines[line-1] != src['raw']:
                provenance_errors.append([collection, line])
    check('literal_source_lines', not provenance_errors, provenance_errors)
    # Independent table parser is kept beside the reference numerical module.
    from nbo_reference import raw_table_checks
    table_result = raw_table_checks(a/'analysis.log', dataset)
    save(out/'reference-tables.json', table_result)
    check('raw_table_values', table_result['status'] == 'pass', table_result)
    if not native:
        return checks
    selected = selections(dataset)
    plan = ['COV_VALIDATION 1', 'scene 0.92 0.7 0.25 2.2 0.03 80',
            'capture "canonical-before"', 'volume "canonical-before" ""', 'seek "panel.nbo"']
    for target, file in [('path', 'analysis.log'), ('archive47', 'FILE.47'),
                          ('aonbo', 'FILE.37'), ('nbomo', 'FILE.49')]:
        plan.append(command('text', 'nbo.'+target, (a/file).as_posix()))
    plan += ['click "nbo.attach"', 'capture "attached"', 'click "nbo.set.nbo"',
             'seek "panel.nbo"']
    for i in selected:
        orbital = dataset['orbitals'][i]
        typed_id = f"nbo.typed.NBO:{orbital['spin']}:{orbital['id']-1}"
        plan += ['seek "nbo.orbital.pick"', 'click "nbo.orbital.pick"',
                 command('click', typed_id), command('volume', f'nbo-{i:03}', i)]
        if dataset['orbitals'][i]['kind'] in ('LP', '3C'):
            plan += ['hover "scene.viewport"', command('capture', f'local-{i:03}')]
    plan += ['hover "scene.viewport"', 'capture "nbo-selected"',
             'export-name "nbo-bundle"', 'click "nbo.export"', 'capture "exported"']
    if languages:
        for language in range(4):
            plan += ['click "language"', 'key "Home"'] + ['key "Down"']*language + ['key "Enter"',
                     'seek "panel.nbo"', 'seek "nbo.path"', 'hover "scene.viewport"',
                     command('capture', f'language-{language}-controls'),
                     'seek "nbo.set.canonical"', 'wheel "nbo.set.canonical" -4',
                     'hover "scene.viewport"', command('capture', f'language-{language}')]
    plan += ['click "nbo.set.canonical"', 'volume "canonical-after" ""',
             'capture "canonical-after"', 'seek "panel.nbo"', 'click "nbo.set.nbo"',
             'seek "panel.nbo"',
             'volume "nbo-restored" "'+str(selected[-1])+'"', 'capture "nbo-restored"']
    plan_path = out/'native.plan'
    plan_path.write_text('\n'.join(plan)+'\n', encoding='utf-8')
    process = run(build/'cov_validation.exe', [reference.fchk_path, '--validation-plan', plan_path,
                  '--validation-output', out/'native', '--validation-background'], out/'native.log')
    session_path = out/'native/session.json'
    session = json.loads(session_path.read_text()) if session_path.exists() else None
    check('native_execution', process['returncode'] == 0 and session is not None
          and session['failed_commands'] == 0, dict(process=process, session=session))
    if session is None:
        return checks
    volume_checks = []
    for path in sorted((out/'native').glob('*.volume.json')):
        d = json.loads(path.read_text()); samples = np.asarray(d['samples'])
        ids = samples[:, 0].astype(int); actual = samples[:, 1]
        dims = np.array([d['nx'], d['ny'], d['nz']], dtype=np.float32)
        ijk = np.column_stack((ids%d['nx'], (ids//d['nx'])%d['ny'], ids//(d['nx']*d['ny']))).astype(np.float32)
        lo = np.asarray(d['grid_box_bohr'][:3], dtype=np.float32)
        hi = np.asarray(d['grid_box_bohr'][3:], dtype=np.float32)
        frac = ijk/(dims-1)
        points = (lo.astype(float)+frac.astype(float)*(hi-lo).astype(float)).astype(np.float32).astype(float)
        spin = 'beta' if d['spin'] == 'beta' else 'alpha'
        ref = reference.evaluate(points, 'nbo' if d['orbital_set'] == 'nbo' else 'mo',
                                 spin, [d['source_index']+1])[0]
        dot = float(actual@ref); sign = 1 if dot >= 0 else -1
        nrms = float(np.linalg.norm(actual-sign*ref)/max(np.linalg.norm(ref), 1e-30))
        cosine = abs(dot)/max(float(np.linalg.norm(actual)*np.linalg.norm(ref)), 1e-30)
        relmax = float(np.max(np.abs(actual-sign*ref))/max(np.max(np.abs(ref)), 1e-30))
        identity_ok = (d['rendered_mo'] < len(dataset['orbitals']) and
                       dataset['orbitals'][d['rendered_mo']]['id'] == d['source_index']+1
                       and dataset['orbitals'][d['rendered_mo']]['spin'] in (d['spin'], 'total')) if d['orbital_set'] == 'nbo' else True
        volume_checks.append(dict(file=path.name, orbital_set=d['orbital_set'], spin=d['spin'],
            source_index=d['source_index'], sample_count=len(ids), nrms=nrms, abs_cosine=cosine,
            relative_max=relmax, global_phase=sign, identity_ok=identity_ok,
            passed=identity_ok and nrms <= 1e-4 and cosine >= 1-1e-7 and relmax <= 1e-3))
    save(out/'actual-volume-comparison.json', volume_checks)
    check('actual_gpu_texture_samples', len(volume_checks) == len(selected)+3 and
          all(x['passed'] for x in volume_checks), volume_checks,
          'NRMS <=1e-4; abs cosine >=1-1e-7; relative max <=1e-3; selected orbitals and <=8192 texture indices')
    before = json.loads((out/'native/canonical-before.volume.json').read_text())
    after = json.loads((out/'native/canonical-after.volume.json').read_text())
    check('canonical_restored', before['orbital_set'] == after['orbital_set'] == 'canonical'
          and before['source_index'] == after['source_index']
          and before['spin'] == after['spin'] and before['samples'] == after['samples'],
          dict(before_source=before['source_index'], after_source=after['source_index'],
               bitwise_sample_equality=before['samples'] == after['samples']))
    export = out/'native/nbo-bundle.nbo.json'
    if export.exists():
        exported = json.loads(export.read_text(encoding='utf-8'))
        check('json_export_matches_import', exported == dataset, dict(source=str(export)))
        for suffix, rows, names in [('npa', dataset['populations'], ('charge', 'core', 'valence', 'rydberg')),
                                    ('nao', dataset['naos'], ('occupation',)),
                                    ('wiberg', dataset['wiberg'], ('value',)),
                                    ('e2', dataset['e2'], ('value', 'energy_gap_hartree', 'fock_hartree'))]:
            file = out/f'native/nbo-bundle.{suffix}.csv'
            parsed = list(csv.DictReader(file.open(encoding='utf-8-sig', newline=''))) if file.exists() else []
            ok = len(parsed) == len(rows) and all(abs(float(p[k])-r[k]) <= 1e-12
                    for p, r in zip(parsed, rows) for k in names)
            for p, r in zip(parsed, rows):
                ok = ok and p['spin'] == r['spin']
                for key in ('atom', 'id', 'atom_a', 'atom_b', 'donor', 'acceptor'):
                    if key in r:
                        ok = ok and int(p[key]) == r[key]
                for key in ('diagonal_fock_hartree', 'spin_density', 'printing_threshold'):
                    if key in p:
                        ok = ok and ((p[key] == '' and r.get(key) is None) or
                            (p[key] != '' and r.get(key) is not None and abs(float(p[key])-r[key]) <= 1e-12))
            check('csv_'+suffix, ok, dict(exported_rows=len(parsed), expected_rows=len(rows)))
        nbo_rows = list(csv.DictReader((out/'native/nbo-bundle.nbo.csv').open(encoding='utf-8-sig', newline='')))
        expected_rows = sum(max(1, len(o['components'])) for o in dataset['orbitals'])
        nbo_csv_ok = len(nbo_rows) == expected_rows
        for p in nbo_rows:
            o = dataset['orbitals'][int(p['index'])]
            nbo_csv_ok = nbo_csv_ok and int(p['id']) == o['id'] and p['spin'] == o['spin'] and p['kind'] == o['kind'] and abs(float(p['occupation'])-o['occupation']) <= 1e-12
            e = o['diagonal_fock_hartree']
            nbo_csv_ok = nbo_csv_ok and ((e is None and p['diagonal_fock_hartree'] == '') or (e is not None and p['diagonal_fock_hartree'] != '' and abs(float(p['diagonal_fock_hartree'])-e) <= 1e-12))
        check('csv_nbo', nbo_csv_ok, dict(exported_rows=len(nbo_rows), expected_rows=expected_rows))
        view = json.loads((out/'native/nbo-bundle.view.json').read_text())
        canonical = view.get('canonical_mo') or {}
        cs = canonical.get('spin', 'alpha'); col = canonical.get('source_index')
        total_error = 0.0
        for item in view.get('contributions', []):
            expected = reference.nbomo[cs][item['matrix_row'], col]**2
            total_error = max(total_error, abs(item['weight_raw']-expected))
        check('nbomo_export', view.get('nbomo_status') == 'mapped' and
              len(view['contributions']) == reference.nbasis and total_error <= 1e-12 and
              view.get('selected_nbo_index') == selected[-1] and
              view.get('rendered_set') == 'nbo' and view.get('rendered_orbital_index') == selected[-1] and
              abs(view['shown_weight_raw']+view['remaining_weight_raw']-view['total_weight_raw']) <= 1e-12,
              dict(max_error=total_error, selected_nbo=view.get('selected_nbo_index'),
                   selection_expected=selected[-1]))
    else:
        check('exports_present', False, 'no bundle written')
    frames = [json.loads(line) for line in (out/'native/frames.jsonl').read_text().splitlines()]
    mismatches = [f['frame'] for f in frames if f['drawn_ui_mo'] != f['rendered_mo']]
    check('ui_scene_identity', not mismatches, dict(frames=len(frames), mismatches=mismatches))
    save(out/'selected-orbitals.json', [dict(index=i, **dataset['orbitals'][i]) for i in selected])
    return checks


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cases', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--case', action='append')
    p.add_argument('--parser-only', action='store_true')
    args = p.parse_args()
    args.output = args.output.resolve(); args.cases = args.cases.resolve(); args.build = args.build.resolve()
    args.output.mkdir(exist_ok=False)
    identity = {name: dict(path=str(args.build/name), sha256=digest(args.build/name))
                for name in ('cov_validation.exe', 'cov_nbo_probe.exe')}
    save(args.output/'build-identity.json', identity)
    reports = []
    for case in sorted(args.cases.glob('NBO-*')):
        if args.case and case.name not in args.case:
            continue
        print('START', case.name, flush=True)
        try:
            checks = case_check(case, args.build, args.output/case.name, not args.parser_only,
                                languages=case.name == 'NBO-02')
            report = dict(case=case.name, status='pass' if all(c['status'] == 'pass' for c in checks) else 'fail', checks=checks)
        except Exception as error:
            report = dict(case=case.name, status='error', reason=f'{type(error).__name__}: {error}')
            import traceback
            (args.output/case.name/'checker-error.txt').write_text(traceback.format_exc())
        save(args.output/case.name/'checks.json', report); reports.append(report)
        save(args.output/'summary.json', reports)
        print('END', case.name, report['status'], flush=True)
    return 0 if reports and all(r['status'] == 'pass' for r in reports) else 2


if __name__ == '__main__':
    raise SystemExit(main())
