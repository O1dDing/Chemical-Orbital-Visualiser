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

"""Meaningful controls for conservative molecule-counting and exact graph matches."""
import unittest

import networkx as nx
import numpy as np

from structure_identity import (assess_primary_molecule_novelty, exact_mapping,
    graph_record, nuclear_graph, proximity_records, smiles_nuclear_graph)


class StructureIdentityTests(unittest.TestCase):
    def test_element_positions_matter_with_equal_composition_and_unlabeled_shape(self):
        a=nuclear_graph([8,6,6],[(0,1),(1,2)])
        b=nuclear_graph([6,8,6],[(0,1),(1,2)])
        self.assertTrue(nx.is_isomorphic(a,b))
        self.assertIsNone(exact_mapping(a,b))

    def test_wl_hash_collision_is_not_a_whole_graph_match(self):
        ring=nuclear_graph([6]*6,[(i,(i+1)%6) for i in range(6)])
        triangles=nuclear_graph([6]*6,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3)])
        self.assertEqual(graph_record(ring)['hash_is_prefilter_only'],graph_record(triangles)['hash_is_prefilter_only'])
        self.assertIsNone(exact_mapping(ring,triangles))

    def test_explicit_hydrogens_distinguish_allene_and_propyne(self):
        self.assertIsNone(exact_mapping(smiles_nuclear_graph('C=C=C'),smiles_nuclear_graph('CC#C')))

    def test_cycle_block_junction_distinguishes_shared_vertex_from_fused_rings(self):
        shared_vertex=nuclear_graph([6]*5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)])
        shared_edge=nuclear_graph([6]*4,[(0,1),(1,2),(2,0),(0,3),(3,1)])
        self.assertEqual(graph_record(shared_vertex)['cycle_block_junction_atoms_zero_based'],[0])
        self.assertEqual(graph_record(shared_edge)['cycle_block_junction_atoms_zero_based'],[])

    def test_rotation_translation_and_atom_permutation_preserve_graph_and_bijection(self):
        xyz=np.array([[0,0,0],[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],dtype=float)*.63
        numbers=np.array([6,1,1,1,1])
        radii={1:.31,6:.76}
        original=proximity_records(numbers,xyz,radii,[1.15,1.25,1.35])
        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
        cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        rotation=np.eye(3)+np.sin(.72)*cross+(1-np.cos(.72))*(cross@cross)
        order=np.array([3,1,4,0,2])
        moved=proximity_records(numbers[order],(xyz@rotation+[2.1,-3.4,5.6])[order],radii,[1.15,1.25,1.35])
        for a,b in zip(original,moved):
            self.assertFalse(a['boundary_edges'] or b['boundary_edges'])
            self.assertIsNotNone(exact_mapping(nuclear_graph(a['graph']['atomic_numbers'],a['graph']['edges_zero_based']),
                                             nuclear_graph(b['graph']['atomic_numbers'],b['graph']['edges_zero_based'])))

    def test_cutoff_boundary_is_retained_as_uncertainty(self):
        record=proximity_records([1,1],[[0,0,0],[.775,0,0]],{1:.31},[1.25])[0]
        self.assertEqual(record['definite_edges_zero_based'],[])
        self.assertEqual(record['possible_edges_zero_based'],[[0,1]])
        self.assertEqual(len(record['boundary_edges']),1)
        self.assertFalse(record['is_chemical_bond_assignment'])

    def test_nonfinite_and_coincident_geometries_are_not_silent_empty_graphs(self):
        for xyz in ([[0,0,0],[0,0,0]],[[0,0,0],[float('nan'),0,0]]):
            with self.assertRaises(ValueError):
                proximity_records([1,1],xyz,{1:.31},[1.25])

    def test_charge_state_does_not_create_novel_nuclear_connectivity(self):
        positive=smiles_nuclear_graph('[CH3+]')
        radical=smiles_nuclear_graph('[CH3]')
        baseline=[{'case_id':'CONTROL-OLD','atomic_numbers':[6,1,1,1],
                   'chemical_scope_reviewed':True,'confirmed_graphs':[graph_record(positive)]}]
        result=assess_primary_molecule_novelty(radical,baseline,primary_scope_reviewed=True)
        self.assertEqual(result['status'],'previously_present_nuclear_connectivity')
        self.assertFalse(result['novelty_credit'])

    def test_unreviewed_same_composition_and_primary_scope_block_novelty_credit(self):
        candidate=smiles_nuclear_graph('CC#C')
        baseline=[{'case_id':'CONTROL-UNRESOLVED','atomic_numbers':[6,6,6,1,1,1,1],
                   'chemical_scope_reviewed':False,'confirmed_graphs':[]}]
        self.assertFalse(assess_primary_molecule_novelty(candidate,baseline,primary_scope_reviewed=True)['novelty_credit'])
        self.assertFalse(assess_primary_molecule_novelty(candidate,[],primary_scope_reviewed=False)['novelty_credit'])
        reviewed=[dict(baseline[0],chemical_scope_reviewed=True,confirmed_graphs=[graph_record(smiles_nuclear_graph('C=C=C'))])]
        result=assess_primary_molecule_novelty(candidate,reviewed,primary_scope_reviewed=True)
        self.assertTrue(result['novelty_credit'])
        self.assertFalse(result['quantum_or_cov_validation_pass'])


if __name__=='__main__':
    unittest.main()
