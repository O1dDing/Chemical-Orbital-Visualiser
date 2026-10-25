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

"""Conservative nuclear-connectivity evidence for validation-set deduplication.

This is independent validation infrastructure, not COV's bond perception. A
radius-based proximity graph is explicitly not a chemical bond assignment.
Graph comparison keeps every explicit nucleus and element, and ignores charge,
spin, bond order and stereochemistry for conservative new-molecule exclusion.
Those omitted distinctions remain metadata; equality here is not a claim of
complete chemical identity. Environment/primary-molecule scope needs review.
"""
from __future__ import annotations

from collections import Counter
import math

import networkx as nx
import numpy as np


def composition_key(atomic_numbers):
    numbers = tuple(atomic_numbers)
    if not numbers or any(isinstance(z, bool) or int(z) != z or not 1 <= z <= 118 for z in numbers):
        raise ValueError('One valid atomic number is required for every nucleus')
    return tuple(sorted(Counter(int(z) for z in numbers).items()))


def nuclear_graph(atomic_numbers, edges):
    numbers = tuple(atomic_numbers)
    composition_key(numbers)
    graph = nx.Graph()
    graph.add_nodes_from((index, {'atomic_number': int(z)}) for index, z in enumerate(numbers))
    for left, right in edges:
        if (isinstance(left, bool) or isinstance(right, bool) or int(left) != left or int(right) != right
                or not 0 <= left < len(numbers) or not 0 <= right < len(numbers) or left == right):
            raise ValueError('Edges require distinct in-range nucleus indices')
        graph.add_edge(int(left), int(right))
    return graph


def graph_record(graph):
    indices = sorted(graph.nodes)
    if indices != list(range(len(indices))):
        raise ValueError('Serialized graphs require contiguous zero-based nucleus indices')
    cyclic_blocks = sorted(sorted(part) for part in nx.biconnected_components(graph)
                           if len(part) >= 3 and graph.subgraph(part).number_of_edges() >= len(part))
    block_memberships = Counter(atom for block in cyclic_blocks for atom in block)
    return {'atomic_numbers': [graph.nodes[index]['atomic_number'] for index in indices],
            'edges_zero_based': sorted([min(a, b), max(a, b)] for a, b in graph.edges),
            'connected_components_zero_based': sorted(sorted(part) for part in nx.connected_components(graph)),
            'cycle_rank': graph.number_of_edges()-graph.number_of_nodes()+nx.number_connected_components(graph),
            'cyclic_blocks_zero_based': cyclic_blocks,
            'cycle_block_junction_atoms_zero_based': sorted(atom for atom, count in block_memberships.items() if count > 1),
            'hash_is_prefilter_only': nx.weisfeiler_lehman_graph_hash(graph, node_attr='atomic_number', iterations=3)}


def from_graph_record(record):
    return nuclear_graph(record['atomic_numbers'], record['edges_zero_based'])


def exact_mapping(first, second):
    """A checked, element-preserving whole-graph bijection or None; never hash equality."""
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        first, second, node_match=lambda a, b: a['atomic_number'] == b['atomic_number'])
    if not matcher.is_isomorphic():
        return None
    mapping = dict(matcher.mapping)
    if set(mapping) != set(first) or set(mapping.values()) != set(second):
        raise AssertionError('The graph matcher did not return a complete bijection')
    for left, target in mapping.items():
        if first.nodes[left]['atomic_number'] != second.nodes[target]['atomic_number']:
            raise AssertionError('The mapping changed an element')
    moved_edges = {frozenset((mapping[a], mapping[b])) for a, b in first.edges}
    if moved_edges != {frozenset(edge) for edge in second.edges}:
        raise AssertionError('The mapping did not preserve every edge and non-edge')
    return mapping


def proximity_records(atomic_numbers, coordinates_angstrom, radii_angstrom, scales,
                      boundary_tolerance_angstrom=1e-6):
    """Store nominal, definite and possible edges under explicitly declared cutoffs."""
    numbers = tuple(atomic_numbers)
    composition_key(numbers)
    xyz = np.asarray(coordinates_angstrom, dtype=float)
    if xyz.shape != (len(numbers), 3) or not np.isfinite(xyz).all():
        raise ValueError('Finite Cartesian coordinates in Angstrom are required')
    scales = tuple(float(scale) for scale in scales)
    if not scales or any(not math.isfinite(scale) or scale <= 0 for scale in scales):
        raise ValueError('Finite positive radius multipliers are required')
    if not math.isfinite(boundary_tolerance_angstrom) or boundary_tolerance_angstrom <= 0:
        raise ValueError('A finite positive boundary tolerance is required')
    radii = np.array([radii_angstrom[int(z)] for z in numbers], dtype=float)
    if not np.isfinite(radii).all() or (radii <= 0).any():
        raise ValueError('Finite positive radii are required for every element')
    pairs = []
    for left in range(len(numbers)):
        for right in range(left+1, len(numbers)):
            distance = float(np.linalg.norm(xyz[left]-xyz[right]))
            if distance <= boundary_tolerance_angstrom:
                raise ValueError(f'Coincident or numerically unresolved nuclei: {left}, {right}')
            radius_sum = float(radii[left]+radii[right])
            pairs.append((left, right, distance, radius_sum))
    records = []
    for scale in scales:
        nominal, definite, possible, boundary = [], [], [], []
        for left, right, distance, radius_sum in pairs:
            margin = distance-scale*radius_sum
            edge = [left, right]
            if margin <= 0:
                nominal.append(edge)
            if margin < -boundary_tolerance_angstrom:
                definite.append(edge)
            if margin <= boundary_tolerance_angstrom:
                possible.append(edge)
            if abs(margin) <= boundary_tolerance_angstrom:
                boundary.append({'edge_zero_based': edge, 'distance_angstrom': distance,
                                 'cutoff_angstrom': scale*radius_sum, 'margin_angstrom': margin})
        records.append({'scale': scale, 'graph': graph_record(nuclear_graph(numbers, nominal)),
                        'definite_edges_zero_based': definite, 'possible_edges_zero_based': possible,
                        'boundary_edges': boundary, 'is_chemical_bond_assignment': False})
    return records


def smiles_nuclear_graph(smiles):
    """An explicitly declared chemical control; no bond-order inference from coordinates."""
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError('The declared chemical structure could not be read')
    molecule = Chem.AddHs(molecule)
    return nuclear_graph([atom.GetAtomicNum() for atom in molecule.GetAtoms()],
                         [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in molecule.GetBonds()])


def assess_primary_molecule_novelty(candidate, baseline, *, primary_scope_reviewed):
    """Conservative exclusion only; this never certifies a quantum/COV validation pass.

    Each baseline item contains atomic_numbers, chemical_scope_reviewed and
    confirmed_graphs (explicitly reviewed graph records, not raw proximity
    guesses). A same-composition unresolved baseline blocks novelty credit.
    """
    if not primary_scope_reviewed:
        return {'status': 'needs_primary_molecule_scope_review', 'novelty_credit': False}
    key = composition_key([candidate.nodes[i]['atomic_number'] for i in candidate])
    unresolved = []
    for item in baseline:
        if composition_key(item['atomic_numbers']) != key:
            continue
        if not item['chemical_scope_reviewed'] or not item['confirmed_graphs']:
            unresolved.append(item['case_id'])
            continue
        for record in item['confirmed_graphs']:
            if composition_key(record['atomic_numbers']) != composition_key(item['atomic_numbers']):
                raise ValueError('A reviewed graph disagrees with its baseline nuclear composition')
            mapping = exact_mapping(candidate, from_graph_record(record))
            if mapping is not None:
                return {'status': 'previously_present_nuclear_connectivity', 'novelty_credit': False,
                        'matched_case_id': item['case_id'], 'mapping_zero_based': mapping}
    if unresolved:
        return {'status': 'needs_same_composition_baseline_review', 'novelty_credit': False,
                'unresolved_case_ids': unresolved}
    return {'status': 'new_reviewed_primary_nuclear_connectivity', 'novelty_credit': True,
            'quantum_or_cov_validation_pass': False}
