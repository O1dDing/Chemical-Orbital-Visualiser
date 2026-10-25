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

"""Frozen original-input graph collection; chemical review follows the full barrier."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

import networkx as nx
import numpy as np
import rdkit
from rdkit import Chem

from structure_identity import exact_mapping, from_graph_record, proximity_records
from validation_process import atomic_json


def sha256(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def original_geometry(path, expected_atom_count):
    """Read the selected original FCHK fields without modifying their values."""
    wanted={b'Number of atoms':('I',None),b'Atomic numbers':('I',expected_atom_count),
            b'Current cartesian coordinates':('R',3*expected_atom_count)}
    found={}
    with path.open('rb') as stream:
        for line in stream:
            label=line[:43].strip()
            if label not in wanted:
                continue
            kind,count=wanted[label]
            fields=line[43:].decode('ascii').split()
            if fields[0]!=kind:
                raise ValueError(f'Unexpected FCHK type for {label!r}')
            if count is None:
                found[label.decode()]=int(fields[1])
            else:
                if fields[1]!='N=' or int(fields[2])!=count:
                    raise ValueError(f'Unexpected FCHK array count for {label!r}')
                values=[]
                while len(values)<count:
                    data=next(stream).decode('ascii').split()
                    values.extend(int(value) if kind=='I' else float(value.replace('D','E')) for value in data)
                if len(values)!=count:
                    raise ValueError(f'FCHK array overrun for {label!r}')
                found[label.decode()]=values
            if len(found)==len(wanted):
                break
    if len(found)!=len(wanted) or found['Number of atoms']!=expected_atom_count:
        raise ValueError('Missing or inconsistent original geometry fields')
    return found


def preflight(root):
    import unittest
    suite=unittest.defaultTestLoader.loadTestsFromName('test_structure_identity')
    started=time.perf_counter_ns()
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    value={'tests_run':result.testsRun,'expected_tests':9,'failures':len(result.failures),
           'errors':len(result.errors),'skipped':len(result.skipped),'expected_failures':len(result.expectedFailures),
           'wall_seconds':(time.perf_counter_ns()-started)/1e9,
           'passed':result.wasSuccessful() and result.testsRun==9 and not result.skipped and not result.expectedFailures,
           'formal_case_passes':0}
    atomic_json(root/'preflight/summary.json',value)
    return 0 if value['passed'] else 1


def collect(root):
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    for relative,expected in manifest['source_sha256'].items():
        if sha256(root/relative)!=expected:
            raise ValueError(f'Frozen source identity mismatch: {relative}')
    campaign_path=Path(manifest['campaign_path'])
    if sha256(campaign_path)!=manifest['campaign_sha256']:
        raise ValueError('Frozen original-case manifest changed')
    campaign=json.loads(campaign_path.read_text(encoding='utf-8'))
    composition_path=Path(manifest['composition_index_path'])
    if sha256(composition_path)!=manifest['composition_index_sha256']:
        raise ValueError('The frozen composition and purpose index changed')
    purposes={row['case_id']:row for row in json.loads(composition_path.read_text(encoding='utf-8'))['cases']}
    cases=campaign['cases']
    if len(cases)!=273 or [row['case_id'] for row in cases]!=[f'OLD-{i:03d}' for i in range(1,274)]:
        raise ValueError('The original 273-case set is not complete and ordered')
    if not json.loads((root/'preflight/summary.json').read_text())['passed']:
        raise ValueError('The complete graph preflight has not passed')
    versions={'networkx':nx.__version__,'numpy':np.__version__,'rdkit':rdkit.__version__}
    if versions!=manifest['expected_runtime_versions']:
        raise ValueError(f'Unexpected graph runtime versions: {versions}')
    periodic=Chem.GetPeriodicTable()
    numbers=sorted({z for row in cases for z in row['fchk_identity']['Atomic numbers']})
    radii={z:float(periodic.GetRcovalent(z)) for z in numbers}
    atomic_json(root/'radii.json',{'source':'RDKit GetPeriodicTable().GetRcovalent','versions':versions,
                                 'units':'angstrom','atomic_number_to_radius':radii,
                                 'radius_scale_policy':manifest['radius_scales'],'chemical_bond_truth':False})
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
    cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    rotation=np.eye(3)+np.sin(.72)*cross+(1-np.cos(.72))*(cross@cross)
    translation=np.array([2.1,-3.4,5.6])
    records=[]
    started=time.perf_counter_ns()
    for ordinal,entry in enumerate(cases,1):
        cid=entry['case_id']
        case_started=time.perf_counter_ns()
        output={'case_id':cid,'case_identity':entry['case_identity'],
                'original_fchk_sha256':entry['original_fchk_sha256'],
                'original_relative_path':entry['original_relative_path'],
                'reference_purpose':purposes[cid]['reference_purpose'],
                'original_input_title':purposes[cid]['original_input_title'],
                'charge':entry['fchk_identity']['Charge'],
                'specified_multiplicity':entry['fchk_identity']['Multiplicity'],
                'chemical_scope_reviewed':False,'confirmed_graphs':[],
                'chemical_graph_count_credit':False,'formal_case_pass':False}
        try:
            path=Path(entry['original_fchk'])
            if sha256(path)!=entry['original_fchk_sha256']:
                raise ValueError('Original FCHK identity mismatch')
            fields=original_geometry(path,entry['fchk_identity']['Number of atoms'])
            for key,value in fields.items():
                if value!=entry['fchk_identity'][key]:
                    raise ValueError(f'Original FCHK geometry cache mismatch: {key}')
            z=fields['Atomic numbers']
            xyz=np.array(fields['Current cartesian coordinates']).reshape((-1,3))*manifest['bohr_to_angstrom']
            original=proximity_records(z,xyz,radii,manifest['radius_scales'],manifest['boundary_tolerance_angstrom'])
            order=np.random.default_rng(manifest['permutation_seed']+ordinal).permutation(len(z))
            transformed=proximity_records(np.array(z)[order],(xyz@rotation+translation)[order],radii,
                                           manifest['radius_scales'],manifest['boundary_tolerance_angstrom'])
            checks=[]
            for first,second in zip(original,transformed):
                mapping=exact_mapping(from_graph_record(first['graph']),from_graph_record(second['graph']))
                checks.append({'scale':first['scale'],'element_and_edge_bijection':mapping is not None,
                               'mapping_zero_based':mapping,
                               'source_boundary_edge_count':len(first['boundary_edges']),
                               'transformed_boundary_edge_count':len(second['boundary_edges'])})
            output.update(status='proximity_evidence_collected',atomic_numbers=z,
                source_coordinates_angstrom=xyz.tolist(),proximity_graphs=original,
                rigid_transform_and_permutation={'rotation_rows':rotation.tolist(),'translation_angstrom':translation.tolist(),
                    'input_reorder_zero_based':order.tolist(),'checks':checks,
                    'all_nominal_graphs_isomorphic':all(row['element_and_edge_bijection'] for row in checks)},
                radii_scale_stable_edges=all(row['graph']['edges_zero_based']==original[0]['graph']['edges_zero_based'] for row in original))
            if not all(row['element_and_edge_bijection'] for row in checks):
                output['transformed_proximity_evidence']=transformed
        except Exception as error:
            output.update(status='collection_error',error=f'{type(error).__name__}: {error}',traceback=traceback.format_exc(),
                          atomic_numbers=entry['fchk_identity']['Atomic numbers'])
        output['collection_wall_seconds']=(time.perf_counter_ns()-case_started)/1e9
        path=root/'cases'/cid/'geometry-graphs.json'
        atomic_json(path,output)
        records.append({'case_id':cid,'status':output['status'],'evidence':str(path.relative_to(root)),
                        'sha256':sha256(path)})
        atomic_json(root/'progress.json',{'phase':'collecting_original_graph_evidence','terminal_cases':ordinal,
            'expected_cases':273,'all_terminal':False,'formal_case_passes':0,'external_molecules_added':0})
    errors=sum(row['status']=='collection_error' for row in records)
    batch={'round_identity':manifest['round_identity'],'all_terminal':True,'expected_cases':273,
           'terminal_cases':273,'collection_errors':errors,'records':records,
           'wall_seconds':(time.perf_counter_ns()-started)/1e9,'formal_case_passes':0,'external_molecules_added':0}
    atomic_json(root/'batch.json',batch)
    atomic_json(root/'progress.json',{key:value for key,value in batch.items() if key!='records'}|
                {'phase':'original_graph_collection_complete'})
    print(json.dumps({key:value for key,value in batch.items() if key!='records'}))
    return 0 if not errors else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=('preflight','collect'))
    parser.add_argument('--root',type=Path,required=True)
    arguments=parser.parse_args()
    raise SystemExit(preflight(arguments.root) if arguments.mode=='preflight' else collect(arguments.root))
