"""Independent acceptance of full NAO transformations and the real focus UI.

This checker never launches chemical calculations. Each invocation keeps a new
evidence directory. Raw vendor matrices and report labels are the references;
production JSON is the observation, not the numerical oracle.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import traceback

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '2'
from nbo_reference import np, _read_w_square
from nbo_integration_suite import run, command


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def raw_labels(report):
    """The common NAO basis is described by the first (total) NAO table."""
    labels = {}
    inside = False
    pattern = re.compile(r'^\s*(\d+)\s+([A-Za-z]+)\s+(\d+)\s+(\S+)\s+(Cor|Val|Ryd)\(\s*(\d+)([spdfghik])\)\s+([\d.+-]+)')
    for line in Path(report).read_text(encoding='utf-8', errors='replace').splitlines():
        if 'NATURAL POPULATIONS:  Natural atomic orbital occupancies' in line:
            if labels:
                break
            inside = True
        if inside and 'Summary of Natural Population Analysis' in line:
            break
        match = pattern.match(line) if inside else None
        if match:
            identity, symbol, atom, angular, typ, principal, letter, occ = match.groups()
            labels[int(identity)] = dict(atom=int(atom), symbol=symbol, angular=angular,
                type=typ, principal_n=int(principal), angular_l='spdfghik'.index(letter),
                occupation=float(occ))
    if not labels:
        raise ValueError('No independently readable producer NAO table')
    return labels


def numerical_case(case, build, out):
    out.mkdir(parents=True)
    checks = []
    def check(name, condition, evidence):
        checks.append(dict(id=name, status='pass' if condition else 'fail', observed=evidence))
    keys = ('fchk', 'report', 'archive47', 'aonbo', 'nbomo', 'naomo', 'aonao', 'naonbo')
    files = {key: Path(case[key]) for key in keys}
    before = {key: digest(path) for key, path in files.items()}
    proc = run(build/'cov_nbo_probe.exe', [files[key] for key in keys[:5]] +
               [out/'probe'] + [files[key] for key in keys[5:]], out/'probe.log')
    check('production_attach', proc['returncode'] == 0, proc)
    dataset_path = out/'probe/dataset.json'
    if not dataset_path.exists():
        return checks
    data = json.loads(dataset_path.read_text(encoding='utf-8'))
    identity = json.loads((out/'probe/canonical-identity.json').read_text(encoding='utf-8'))
    check('canonical_identity_unchanged', identity['preserved'], identity)
    archive_text = files['archive47'].read_text(encoding='ascii')
    n = int(re.search(r'NBAS\s*=\s*(\d+)', archive_text[:1000]).group(1))
    header = archive_text.split('$END', 1)[0].upper()
    unrestricted = re.search(r'\bOPEN\b', header) is not None
    transforms = _read_w_square(files['naomo'], n, 'MOs in the NAO basis:', unrestricted)
    a = _read_w_square(files['aonao'], n, 'NAOs in the AO basis:', False)
    u = _read_w_square(files['naonbo'], n, 'NBOs in the NAO basis:', unrestricted)
    for name, reference in [('NAOMO', transforms), ('AONAO', a), ('NAONBO', u)]:
        errors = []
        for spin, matrix in reference.items():
            tag = spin if unrestricted and name != 'AONAO' else 'total'
            found = [m for m in data['matrices'] if m['kind'] == name and m['spin'] == tag]
            errors.append(float(np.max(np.abs(np.asarray(found[0]['values']).reshape(matrix.shape)-matrix)))
                          if len(found) == 1 else None)
        check('raw_matrix_'+name, all(x is not None and x <= 1e-12 for x in errors), errors)
    labels = raw_labels(files['report'])
    check('all_raw_NAO_labels', set(labels) == set(range(1, n+1)), len(labels))
    all_mo = data.get('mo_decompositions', [])
    check('decomposition_identity_coverage', len(all_mo) == len(identity['orbitals']) and
          {m['canonical_index'] for m in all_mo} == set(range(len(identity['orbitals']))),
          dict(decompositions=len(all_mo), canonical=len(identity['orbitals'])))
    errors, label_errors, occupation_errors = [], [], []
    populations = np.zeros(n)
    for mo in all_mo:
        if not mo['available']:
            continue
        spin = mo['spin'] if mo['spin'] != 'total' else 'alpha'
        column = transforms[spin][:, mo['source_orbital_index']]
        rows = mo['rows']
        if len(rows) != n:
            errors.append(dict(index=mo['canonical_index'], missing_rows=len(rows)))
            continue
        ids = np.array([r['nao_id']-1 for r in rows])
        weights = np.asarray([r['weight'] for r in rows])
        coefficients = np.asarray([r['coefficient'] for r in rows])
        reference = column[ids]
        phase_error = min(float(np.max(np.abs(coefficients-reference))),
                          float(np.max(np.abs(coefficients+reference))))
        weight_error = float(np.max(np.abs(weights-reference**2)))
        group_error = 0.0
        for group_type in ('atoms', 'shells'):
            expected = {}
            for nao_id, label in labels.items():
                key = (label['atom'],) if group_type == 'atoms' else (
                    label['atom'], label['principal_n'], label['angular_l'])
                expected[key] = expected.get(key, 0.0) + float(column[nao_id-1]**2)
            groups = mo[group_type]
            for group in groups:
                key = (group['atom'],) if group_type == 'atoms' else (
                    group['atom'], group['principal_n'], group['angular_l'])
                group_error = max(group_error, abs(group['weight']-expected.pop(key, float('inf'))))
            if expected:
                group_error = float('inf')
        if phase_error > 1e-12 or weight_error > 1e-12 or group_error > 2e-12:
            errors.append(dict(index=mo['canonical_index'], phase_error=phase_error,
                               weight_error=weight_error, group_error=str(group_error)))
        if abs(sum(weights)-1) > 1e-6 or abs(sum(weights)-mo['weight_sum']) > 2e-12:
            errors.append(dict(index=mo['canonical_index'], weight_sum=float(sum(weights))))
        for row in rows:
            raw = labels[row['nao_id']]
            for key in ('atom', 'symbol', 'angular', 'principal_n', 'angular_l'):
                if row[key] != raw[key]:
                    label_errors.append([mo['canonical_index'], row['nao_id'], key])
            expected_electrons = row['weight'] * mo['occupation']
            if row['electron_contribution'] is None or abs(row['electron_contribution']-expected_electrons) > 1e-12:
                occupation_errors.append([mo['canonical_index'], row['nao_id']])
            populations[row['nao_id']-1] += expected_electrons
    check('independent_weights_and_groups', not errors, errors)
    check('independent_NAO_labels', not label_errors, label_errors)
    check('occupation_weighted_contributions', not occupation_errors, occupation_errors)
    direct_spins = {v['spin'] for v in data.get('nao_validation', [])
                    if v['available'] and v['direct_fchk_coefficients']}
    complete = 'total' in direct_spins or direct_spins == {'alpha', 'beta'}
    if complete:
        population_error = max(abs(populations[k-1]-v['occupation']) for k,v in labels.items())
        check('independent_total_NAO_populations', population_error <= 1e-5,
              dict(max_error=population_error, printed_rounding=0.000005))
    else:
        check('partial_spin_scope_explicit', unrestricted and all(
            v['direct_fchk_coefficients'] or v['status']=='verified_archive_only'
            for v in data['nao_validation']) and not any(
                m['available'] and m['spin'] not in direct_spins for m in all_mo),
            data['nao_validation'])
    check('source_files_unchanged', before == {k:digest(v) for k,v in files.items()}, before)
    save(out/'checks.json', checks)
    return checks


def native_case(case, build, out):
    checks = []
    def check(name, condition, observed):
        checks.append(dict(id=name, status='pass' if condition else 'fail', observed=observed))
    plan = ['COV_VALIDATION 1', 'scene 0.92 0.7 0.25 2.2 0.03 80',
            'volume "canonical-before" ""', 'seek "panel.nbo"']
    for key, target in [('report', 'path'), ('archive47', 'archive47'),
                        ('aonbo', 'aonbo'), ('nbomo', 'nbomo'), ('naomo', 'naomo'),
                        ('aonao', 'aonao'), ('naonbo', 'naonbo')]:
        plan.append(command('text', 'nbo.'+target, Path(case[key]).as_posix()))
    plan += ['click "nbo.attach"', 'seek "nbo.focus"', 'capture "focus-attached"',
             'click "nbo.focus.mo"', 'key "Home"', 'key "Enter"',
             'click "nbo.focus.show_all_atoms"', 'click "nbo.focus.group.atom"',
             'seek "nbo.focus.graph"', 'wheel "nbo.focus.graph" -5', 'hover "scene.viewport"', 'capture "focus-atoms"',
             'seek "panel.nbo"', 'export-name "all-atoms"', 'click "nbo.export"',
             'seek "nbo.focus"', 'click "nbo.focus.hide_all_atoms"',
             'seek "panel.nbo"', 'export-name "hidden-atoms"', 'click "nbo.export"',
             'seek "nbo.focus"', 'click "nbo.focus.show_all_atoms"',
             'click "nbo.focus.group.shell"', 'seek "nbo.focus.graph"', 'wheel "nbo.focus.graph" -5', 'hover "scene.viewport"',
             'capture "focus-shells"', 'seek "panel.nbo"',
             'export-name "all-shells"', 'click "nbo.export"',
             'click "nbo.set.nbo"', 'seek "panel.nbo"', 'click "nbo.set.canonical"',
             'volume "canonical-after" ""', 'seek "nbo.focus"',
             'hover "scene.viewport"', 'capture "focus-restored"',
             'click "nbo.focus.inspect_3d"', 'volume "focused-inspection" ""',
             'seek "panel.nbo"', 'export-name "focused-inspection"', 'click "nbo.export"']
    ligand_groups = [[2,3],[4,5],[6,7],[8,9]] if case['case_id']=='NBO-F12' else []
    if ligand_groups:
        plan += ['seek "nbo.focus"', 'click "nbo.focus.group.clear"']
        for atom_ids in ligand_groups:
            plan += ['click "nbo.focus.group.clear_draft"']
            plan += [command('click','nbo.focus.ligand.'+str(atom)) for atom in atom_ids]
            plan += ['click "nbo.focus.group.add"']
        plan += ['click "nbo.focus.group.mode.ligand_l"', 'seek "nbo.focus.graph"', 'wheel "nbo.focus.graph" -5',
                 'hover "scene.viewport"', 'capture "four-ligands"',
                 'seek "panel.nbo"', 'export-name "four-ligands"', 'click "nbo.export"',
                 'seek "nbo.focus"', 'click "nbo.focus.include_rydberg"',
                 'seek "panel.nbo"', 'export-name "with-rydberg"', 'click "nbo.export"',
                 'seek "nbo.focus"', 'click "nbo.focus.include_rydberg"']
    if case['case_id'] == 'NBO-01':
        plan += ['click "language"', 'key "Home"', 'key "Down"', 'key "Enter"',
                 'seek "panel.nbo"', 'seek "nbo.focus"', 'hover "scene.viewport"', 'capture "focus-chinese-controls"',
                 'seek "nbo.focus.graph"', 'wheel "nbo.focus.graph" -5', 'hover "scene.viewport"', 'capture "focus-chinese-graph"']
    path = out/'native.plan'
    path.write_text('\n'.join(plan)+'\n', encoding='utf-8')
    proc = run(build/'cov_validation.exe', [Path(case['fchk']), '--validation-plan', path,
               '--validation-output', out/'native', '--validation-background'], out/'native.log')
    session_file = out/'native/session.json'
    session = json.loads(session_file.read_text()) if session_file.exists() else None
    check('native_real_UI_commands', proc['returncode'] == 0 and session is not None
          and session['failed_commands'] == 0, dict(process=proc, session=session))
    if not session:
        return checks
    exports = {}
    for name in ('all-atoms', 'hidden-atoms', 'all-shells'):
        base = out/'native'/name
        needed = [Path(str(base)+suffix) for suffix in
                  ('.focus.json', '.focus.csv', '.focus.svg', '.focus.png', '.nbo.json')]
        check('export_files_'+name, all(p.exists() and p.stat().st_size for p in needed),
              [str(p) for p in needed if not p.exists()])
        if not all(p.exists() for p in needed):
            continue
        view = json.loads(needed[0].read_text(encoding='utf-8'))
        rows = list(csv.DictReader(needed[1].open(encoding='utf-8-sig', newline='')))
        data = json.loads(needed[4].read_text(encoding='utf-8'))
        exports[name] = view
        expected_rows = sum(max(1, len(mo['rows'])) for mo in data['mo_decompositions'])
        check('full_numerical_export_'+name, len(rows) == expected_rows and
              len(view['mo_rows']) == len(data['mo_decompositions']),
              dict(csv_rows=len(rows), expected=expected_rows, mo_rows=len(view['mo_rows'])))
        known = {(mo['canonical_index'], row['nao_id']):row
                 for mo in data['mo_decompositions'] for row in mo['rows']}
        errors = []
        for row in rows:
            if not row['nao_id']:
                continue
            source = known.get((int(row['canonical_index']), int(row['nao_id'])))
            if source is None or any(abs(float(row[key])-source[key]) > 1e-12
                                     for key in ('coefficient', 'weight')):
                errors.append([row['canonical_index'], row['nao_id']])
        check('numerical_export_values_'+name, not errors, errors)
        visible_sum = sum(float(r['weight']) for r in rows if r['graph_visible']=='true')
        node_sum = sum(node['raw_weight'] for node in view['graph_nodes'])
        check('graph_selected_weight_'+name, abs(visible_sum-node_sum) < 1e-10,
              dict(rows=visible_sum, nodes=node_sum))
        check('graph_identity_'+name, view['status']=='available' and
              view['focused_mo_index'] in view['central_mo_indices'] and
              all(0 <= i < len(data['mo_decompositions']) for i in view['central_mo_indices']) and
              all(node['type']=='composition' and node['visibility_reason'] for node in view['graph_nodes']) and
              all(edge['type']=='composition' and edge['evidence'] for edge in view['graph_edges']) and
              len(view['graph_nodes']) == len(view['graph_edges']),
              dict(file=str(needed[0]), snapshot_id=view['snapshot_id'],
                   focused_mo_index=view['focused_mo_index'], central_mo_indices=view['central_mo_indices'],
                   nodes=len(view['graph_nodes']), edges=len(view['graph_edges'])))
    if len(exports) == 3:
        first, hidden, shells = [exports[x] for x in ('all-atoms','hidden-atoms','all-shells')]
        check('hiding_does_not_change_data_or_MO_set', first['mo_rows']==hidden['mo_rows']==shells['mo_rows']
              and first['central_mo_indices']==hidden['central_mo_indices']==shells['central_mo_indices']
              and first['snapshot_id']==hidden['snapshot_id']==shells['snapshot_id'],
              dict(snapshot_ids=[x['snapshot_id'] for x in (first,hidden,shells)],
                   central_indices=first['central_mo_indices']))
        check('explicit_graph_visibility', bool(first['graph_nodes']) and not hidden['graph_nodes']
              and bool(shells['graph_nodes']), [len(x['graph_nodes']) for x in (first,hidden,shells)])
    before = json.loads((out/'native/canonical-before.volume.json').read_text())
    after = json.loads((out/'native/canonical-after.volume.json').read_text())
    check('canonical_render_preserved', before['samples']==after['samples'] and
          before['source_index']==after['source_index'] and before['spin']==after['spin'],
          dict(before=before['source_index'], after=after['source_index']))
    inspected = json.loads((out/'native/focused-inspection.volume.json').read_text())
    inspected_view = json.loads((out/'native/focused-inspection.focus.json').read_text())
    focused_index = inspected_view['focused_mo_index']
    numerical = next(x for x in inspected_view['mo_rows'] if x['canonical_index']==focused_index)
    check('explicit_focus_to_3D', inspected['orbital_set']=='canonical' and
          inspected['source_index']==numerical['source_orbital_index'] and
          inspected_view['central_mo_indices']==exports['all-shells']['central_mo_indices'],
          dict(focused=focused_index, rendered=inspected['source_index'], spin=inspected['spin']))
    if ligand_groups:
        normal = json.loads((out/'native/four-ligands.focus.json').read_text())
        rydberg = json.loads((out/'native/with-rydberg.focus.json').read_text())
        errors=[]
        for label, view in [('valence',normal),('rydberg',rydberg)]:
            records = list(csv.DictReader((out/f'native/{"four-ligands" if label=="valence" else "with-rydberg"}.focus.csv').open(encoding='utf-8-sig', newline='')))
            mo = next(x for x in view['mo_rows'] if x['canonical_index']==view['focused_mo_index'])
            if [g['atom_ids'] for g in view['ligand_groups']] != ligand_groups:
                errors.append([label,'group-membership'])
            for group in view['ligand_groups']:
                for angular in group['angular_weights']:
                    members = [r for r in mo['rows'] if r['atom'] in group['atom_ids'] and r['angular_l']==angular['angular_l']]
                    expected_full=sum(r['weight'] for r in members)
                    expected_selected=sum(float(r['weight']) for r in records if r['graph_visible']=='true'
                        and int(r['atom']) in group['atom_ids'] and int(r['angular_l'])==angular['angular_l'])
                    if abs(angular['full_weight']-expected_full)>1e-12 or abs(angular['selected_shell_subtotal']-expected_selected)>1e-12:
                        errors.append([label,group['id'],angular['angular_l']])
            visible=sum(float(r['weight']) for r in records if r['graph_visible']=='true')
            node_sum=sum(n['raw_weight'] for n in view['graph_nodes'])
            if abs(visible-node_sum)>1e-10:
                errors.append([label,'double-counted-group',visible,node_sum])
        check('four_explicit_ligands_and_l_decomposition', not errors, errors)
        check('Rydberg_selection_preserves_data', normal['mo_rows']==rydberg['mo_rows'] and any(
            n.get('principal_n')==4 and n.get('angular_l')==1 and n['atom_ids']==[1]
            for n in rydberg['graph_nodes']), dict(nodes_before=len(normal['graph_nodes']),nodes_after=len(rydberg['graph_nodes'])))
    save(out/'native-checks.json', checks)
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', action='append')
    parser.add_argument('--native', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = json.loads(args.cases.read_text(encoding='utf-8'))
    if isinstance(cases, dict):
        cases = cases['cases']
    for case in cases:
        for key in ('fchk', 'report', 'archive47', 'aonbo', 'nbomo', 'naomo', 'aonao', 'naonbo', 'sao', 'reference'):
            if case.get(key) and not Path(case[key]).is_absolute():
                case[key] = str((args.cases.resolve().parent / case[key]).resolve())
    reports = []
    for case in cases:
        case_id = case.get('case_id', case.get('case'))
        if args.case and case_id not in args.case:
            continue
        print('START', case_id, flush=True)
        directory = args.output/case_id
        try:
            checks = numerical_case(case, args.build, directory)
            if args.native and all(x['status']=='pass' for x in checks):
                checks += native_case(case, args.build, directory)
            report = dict(case=case_id, status='pass' if checks and all(x['status']=='pass' for x in checks) else 'fail', checks=checks)
        except Exception as error:
            directory.mkdir(parents=True, exist_ok=True)
            (directory/'checker-error.txt').write_text(traceback.format_exc(), encoding='utf-8')
            report = dict(case=case_id, status='error', reason=str(error))
        reports.append(report)
        save(args.output/'summary.json', reports)
        print('END', case_id, report['status'], flush=True)
    return 0 if reports and all(r['status']=='pass' for r in reports) else 2


if __name__ == '__main__':
    raise SystemExit(main())
