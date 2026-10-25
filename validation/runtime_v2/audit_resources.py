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

"""Read-only REF-001 workload census. Never parses original small-basis NBasis."""
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

from policy import ELEMENTS, VERSION, basis_count, preferred_cores, historical_reference_digest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def census(home, output):
    config = json.loads((home / 'runtime/config.json').read_text())
    work = Path(config['work_root'])
    references = work / 'references/REF-001'
    manifest = json.loads((references / 'reference-candidates.json').read_text())
    campaign = json.loads((Path(config['data_root']) / 'campaign/campaign.json').read_text())
    originals = {r['case_id']: r for r in campaign['cases']}
    history = defaultdict(list)
    roots = [work / 'jobs/REF-001'] + [p / 'jobs' for p in sorted(work.glob('runtime*'))]
    for root in roots:
        for path in sorted(root.rglob('job.gjf')):
            match = re.search(r'OLD-\d{3}', str(path))
            if not match or any(part.startswith('snapshot-') for part in path.parts):
                continue
            contents = path.read_text(encoding='ascii', errors='replace')
            nproc = re.findall(r'%nprocshared\s*=\s*(\d+)', contents, re.I)
            log = path.with_name('job.log')
            if not log.exists():
                # The v1 runner used initial.log beside initial.gjf in some bundles.
                log = path.with_suffix('.log')
            logs = log.read_text(encoding='ascii', errors='replace') if log.exists() else ''
            nbasis = sorted(set(map(int, re.findall(r'NBasis\s*=\s*(\d+)', logs))))
            history[match[0]].append({'input': str(path), 'input_sha256': sha(path),
                'nprocshared': sorted(set(map(int, nproc))), 'log': str(log) if log.exists() else None,
                'log_sha256': sha(log) if log.exists() else None, 'nbasis': nbasis,
                'scf_convergence_failures': len(re.findall(r'Convergence failure|SCF has not converged', logs, re.I)),
                'normal_termination': 'Normal termination of Gaussian' in logs})
    rows = []
    for ref in manifest['candidates']:
        case = ref['case_id']; original = originals[case]; f = original['fchk_identity']
        if historical_reference_digest(ref) != ref['identity']:
            raise ValueError('Reference identity changed: ' + case)
        folder = references / 'reference-inputs' / case
        basis = folder / 'basis.gbs'; initial = folder / 'initial.gjf'
        if sha(basis) != ref['basis_sha256'] or sha(initial) != ref['initial_input_sha256']:
            raise ValueError('Frozen reference input changed: ' + case)
        contents = initial.read_text(encoding='ascii')
        strict_count = basis_count(basis.read_text(encoding='ascii'), f['Atomic numbers'], ref['shell_flags'] == '5D 7F')
        proofs = [r for r in history[case] if r['nbasis']]
        observed = {n for r in proofs for n in r['nbasis']}
        if observed and observed != {strict_count}:
            raise ValueError(f'{case} reference NBasis mismatch: {observed} vs {strict_count}')
        elements = Counter(ELEMENTS[z] for z in f['Atomic numbers'])
        symbols = (['C'] + (['H'] if 'H' in elements else []) + sorted(set(elements)-{'C','H'})) if 'C' in elements else sorted(elements)
        formula = ''.join(s + (str(elements[s]) if elements[s] != 1 else '') for s in symbols)
        route = re.search(r'(?ms)^#.*?(?=\n\s*\n)', contents)[0].replace('\n', ' ')
        tier16, tier8 = preferred_cores(strict_count, physical=16), preferred_cores(strict_count, physical=8)
        ecp = {ELEMENTS[int(z)]: n for z, n in ref['ecp_core_electrons_by_atomic_number'].items() if n}
        flags = ['hybrid exact exchange', 'VeryTight SCF', 'SuperFineGrid']
        for key in ('Opt=', 'CalcFC', 'Freq'):
            if key in route: flags.append(key.rstrip('='))
        if ref['basis_name'].endswith('D'): flags.append('diffuse functions')
        if ecp: flags.append('ECP')
        if f['Multiplicity'] != 1 or f.get('has_beta_coefficients'): flags.append('open shell / unrestricted')
        if max(f['Atomic numbers']) > 18: flags.append('heavy element')
        name = original.get('original_input', {}).get('title') or f.get('title') or formula
        rows.append({'case_id': case, 'name': name, 'name_source': 'frozen original title (label, not reidentified molecule)',
            'original_relative_path': original['original_relative_path'], 'formula': formula,
            'atom_count': len(f['Atomic numbers']), 'elements': dict(sorted(elements.items())),
            'charge': f['Charge'], 'multiplicity': f['Multiplicity'], 'basis': ref['basis_name'],
            'pure_cartesian': ref['shell_flags'], 'nbasis': strict_count,
            'nbasis_source': 'same REF-001 Gaussian log; strict shell count agrees' if proofs else 'strict frozen contracted-shell count',
            'nbasis_log_evidence': [{'path': r['log'], 'sha256': r['log_sha256']} for r in proofs],
            'basis_sha256': ref['basis_sha256'], 'reference_identity': ref['identity'],
            'frozen_input_sha256': ref['initial_input_sha256'], 'functional': ref['method'],
            'hybrid': True, 'exact_exchange_fraction': .25, 'opt': 'Opt=' in route,
            'calcfc': 'CalcFC' in route, 'freq': bool(re.search(r'\bFreq\b', route)),
            'scf': re.search(r'SCF=\([^)]*\)', route)[0], 'grid': 'SuperFineGrid',
            'ecp_core_electrons': ecp, 'route': route, 'frozen_nprocshared': 4,
            'historical_nprocshared': sorted({n for r in history[case] for n in r['nprocshared']}),
            'historical_execution': history[case], 'current_live_nprocshared': None,
            'convergence_failure_count': sum(r['scf_convergence_failures'] for r in history[case]),
            'cores_16': tier16, 'cores_8': tier8, 'policy_version': VERSION,
            'reason': f'NBasis={strict_count}; ' + ('>=450 heavy tier' if tier16 == 7 else '300-449 medium tier' if tier16 == 5 else '<300 small tier') + '; ' + ', '.join(flags),
            'collected_in_active_runtime': (work / config['runtime_directory'] / 'jobs' / case / 'result.json').exists()})
    assert [r['case_id'] for r in rows] == [f'OLD-{i:03d}' for i in range(1,274)]
    output.mkdir(parents=True, exist_ok=True)
    heavy = [r for r in rows if r['cores_16'] > 5]
    counts = {'16_physical': dict(Counter(r['cores_16'] for r in rows)), '8_physical': dict(Counter(r['cores_8'] for r in rows))}
    report = {'policy_version': VERSION, 'case_count': len(rows), 'counts': counts,
        'heavy_count': len(heavy), 'observed_reference_nbasis_cases': sum(bool(r['nbasis_log_evidence']) for r in rows),
        'current_live_nprocshared_note': 'Not inferred from historical files; use separate fresh OS process receipt.',
        'cases': rows}
    (output / 'all-273.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    for name, values in [('all-273.csv', rows), ('heavy-over-5.csv', heavy)]:
        with (output / name).open('w', encoding='utf-8-sig', newline='') as stream:
            fields = [k for k in rows[0] if k not in ('historical_execution', 'nbasis_log_evidence')]
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
            writer.writerows({k: json.dumps(v, ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items() if k in fields} for r in values)
    lines = ['# REF-001 adaptive scheduling census', '', f'273 cases; tier counts: `{counts}`.', '',
        'Names are frozen input labels. Formulas derive from atomic numbers; NBasis never derives from the old small-basis FCHK.',
        'The 7/5/4 profile is a conservative admission policy, not a measured performance optimum. All cases retain their frozen scientific settings.', '',
        '| Case | Frozen name | Formula | NBasis | 16-core tier | 8-core tier | Source |', '|---|---|---|---:|---:|---:|---|']
    lines += [f"| {r['case_id']} | {r['name'].replace('|','/')} | {r['formula']} | {r['nbasis']} | {r['cores_16']} | {r['cores_8']} | {'Gaussian log + count' if r['nbasis_log_evidence'] else 'strict count'} |" for r in rows]
    lines += ['', '## Complete heavy list (>5 cores on 16 physical cores)', '', f'{len(heavy)} cases. Completed cases are classified without scheduling a rerun.', '', '| Case | Formula | NBasis | Reason |', '|---|---|---:|---|']
    lines += [f"| {r['case_id']} | {r['formula']} | {r['nbasis']} | {r['reason']} |" for r in heavy]
    (output / 'all-273.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'cases'},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--home', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args(); census(args.home,args.output)
