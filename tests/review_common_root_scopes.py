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

"""Independent post-barrier checks of the declared local and graph scopes.

Reference S and source coefficients are reused from a completed, identified
IOData/GBasis review. COV values are observations, never reference values.
This is an implementation subset, not a molecular-physics acceptance report.
"""
from pathlib import Path
import argparse, hashlib, itertools, json, math, os, sys, time, traceback
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
sys.path.insert(0,r'F:\Dev\cov-validation-20260905\reference-deps')
import numpy as np
from iodata import load_one

CRITERIA={
    'projection_absolute_error':2e-8,
    'angular_function_absolute_error':1e-11,
    'rotation_absolute_error':2e-8,
    'energy_absolute_error_hartree':1e-11,
    'mayer_absolute_error':2e-8,
    'metric_gram_rank_multiplier':64,
    'metric_negative_eigenvalue_multiplier':64,
    'cycle_enumeration_state_limit':250000,
    'cycle_enumeration_count_limit':10000,
    'topology_mayer_cutoff':0.05,
    'topology_electronic_radius_factor':1.45,
    'topology_geometric_radius_factor':1.22,
    'angstrom_to_bohr':1.8897261254578281,
    'scope':'Unweighted S-metric projectors, source-member identities, Mayer/distance graph contract and exhaustive channel-associated cycle witnesses. No physical irrep or chemical-state acceptance inferred.',
    'reference_reuse':'Only identified independent S, source MO coefficients and validated same-source density arrays; COV output/pass decisions are not reference caches.',
    'projection_reference':'Independent symmetric-eigenvalue whitening of positive-semidefinite S followed by normalized-column SVD; negative eigenvalues within the declared arithmetic bound are zeroed, and larger negatives reject the reference. Angular functions use associated-Legendre recurrence and sphere quadrature.',
    'radii_scope':'Frozen declared structural-neighbour model constants; their universal chemical validity is not tested here.',
}
RADII=[.85,.31,.28,1.28,.96,.84,.76,.71,.66,.57,.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,
2.03,1.76,1.70,1.60,1.53,1.39,1.39,1.32,1.26,1.24,1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,
2.20,1.95,1.90,1.75,1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,
2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,
1.51,1.44,1.41,1.36,1.36,1.32,1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,2.06,2.00,1.96,1.90,1.87,1.80,1.69,
1.68,1.68,1.65,1.67,1.73,1.76,1.61,1.57,1.49,1.43,1.41,1.34,1.29,1.28,1.21,1.22,1.36,1.43,1.62,1.75,1.65,1.57]
ORDERS={0:[(0,0,0)],1:[(1,0,0),(0,1,0),(0,0,1)],
2:[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
3:[(3,0,0),(0,3,0),(0,0,3),(1,2,0),(2,1,0),(2,0,1),(1,0,2),(0,1,2),(0,2,1),(1,1,1)],
4:[(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),(0,3,1),(1,0,3),(0,1,3),(2,2,0),(2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)]}

def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def write(p,v):
    p=Path(p);t=p.with_suffix(p.suffix+'.tmp')
    t.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');t.replace(p)
def real_harmonics(l,unit):
    z=np.clip(unit[:,2],-1,1);phi=np.arctan2(unit[:,1],unit[:,0]);values=[]
    for m in range(l+1):
        pmm=np.ones_like(z)
        for q in range(1,m+1):pmm*=-(2*q-1)*np.sqrt(np.maximum(0,1-z*z))
        if l==m:p=pmm
        else:
            pm1=(2*m+1)*z*pmm
            if l==m+1:p=pm1
            else:
                previous,current=pmm,pm1
                for degree in range(m+2,l+1):
                    previous,current=current,((2*degree-1)*z*current-(degree+m-1)*previous)/(degree-m)
                p=current
        p=p*math.sqrt((2*l+1)/(4*math.pi)*math.factorial(l-m)/math.factorial(l+m))
        if m==0:values.append(p)
        else:values.extend((math.sqrt(2)*p*np.cos(m*phi),math.sqrt(2)*p*np.sin(m*phi)))
    return np.stack(values,axis=1)

z,zw=np.polynomial.legendre.leggauss(14);phi=np.arange(32)*2*math.pi/32
zz=np.repeat(z,32);pp=np.tile(phi,14)
POINTS=np.stack((np.sqrt(1-zz**2)*np.cos(pp),np.sqrt(1-zz**2)*np.sin(pp),zz),axis=1)
WEIGHTS=np.repeat(zw,32)*2*math.pi/32
def angular_reference(L,pure,l,rotation):
    if pure:source=real_harmonics(L,POINTS)
    else:
        source=np.stack([np.prod(POINTS**np.asarray(p),axis=1)/
            math.sqrt(math.prod(math.prod(range(2*k-1,0,-2)) for k in p)) for p in ORDERS[L]],axis=1)
    target=real_harmonics(l,POINTS@rotation)
    coefficients=np.linalg.lstsq(source*np.sqrt(WEIGHTS[:,None]),target*np.sqrt(WEIGHTS[:,None]),rcond=None)[0]
    error=float(np.max(np.abs(source@coefficients-target)))
    if error>CRITERIA['angular_function_absolute_error']:raise ValueError('Independent angular fit failed: '+str(error))
    return coefficients

def metric_root(S):
    eigenvalues,vectors=np.linalg.eigh((S+S.T)*.5)
    bound=CRITERIA['metric_negative_eigenvalue_multiplier']*len(S)*np.finfo(float).eps*max(float(np.max(np.abs(eigenvalues))),np.finfo(float).tiny)
    if eigenvalues[0]<-bound:raise ValueError('Independent AO metric is indefinite beyond arithmetic resolution')
    return np.sqrt(np.maximum(eigenvalues,0))[:,None]*vectors.T

class CaseReview:
    def __init__(self,mol,S,C,production,density):
        self.mol=mol;self.S=S;self.C=C;self.p=production;self.density=density
        self.checks=[];self.counts={'unique_projections':0,'scope_occurrences':0,'topology_assignments':0,'gap_records':0}
        self.maximum_projection_error=0.;self.workspace_cache={};self.projection_cache={};self.target_cache={}
        self.shells=[];offset=0;self.basis_atom=[]
        for shell in mol.obasis.shells:
            for L,kind in zip(shell.angmoms,shell.kinds):
                L=int(L);pure=kind=='p';count=2*L+1 if pure else (L+1)*(L+2)//2
                self.shells.append({'atom':int(shell.icenter),'L':L,'pure':pure,'offset':offset,'count':count})
                self.basis_atom.extend([int(shell.icenter)]*count);offset+=count
        assert offset==len(S) and len(RADII)==119
        self.basis_atom=np.asarray(self.basis_atom)
        self.root=metric_root(S)
        self.reference_spins=[0]*C.shape[1]
        if mol.mo.kind=='unrestricted':self.reference_spins[len(mol.mo.coeffsa[0]):]=[1]*len(mol.mo.coeffsb[0])
        self.check('SOURCE-SPIN-IDENTITY',[x['spin'] for x in production['orbitals']]==self.reference_spins)

    def check(self,name,ok,observed=None):
        self.checks.append({'check':name,'status':'pass' if ok else 'fail','observed':observed})
    def close(self,name,got,wanted,tol=None):
        tol=CRITERIA['projection_absolute_error'] if tol is None else tol
        error=abs(float(got)-float(wanted)) if got is not None and wanted is not None else None
        self.check(name,error is not None and math.isfinite(error) and error<=tol,{'recorded':got,'independent':wanted,'absolute_error':error,'tolerance':tol})
        if error is not None:self.maximum_projection_error=max(self.maximum_projection_error,error)

    def orthonormal(self,coeff):
        if not coeff.shape[1]:return np.empty((len(self.S),0))
        values=self.root@coeff
        norms=np.linalg.norm(values,axis=0);positive=norms>0
        values[:,positive]/=norms[positive];values[:,~positive]=0
        u,s,_=np.linalg.svd(values,full_matrices=False)
        threshold=CRITERIA['metric_gram_rank_multiplier']*max(values.shape)*np.finfo(float).eps*(s[0]**2 if len(s) else 0)
        return u[:,s*s>threshold]

    def workspace(self,atom,rotation):
        key=(atom,tuple(rotation.ravel()))
        if key in self.workspace_cache:return self.workspace_cache[key]
        n=len(self.S);indices=[];columns=[[] for _ in range(25)];shell_ids=[[] for _ in range(25)];source=[]
        for sid,shell in enumerate(self.shells):
            if shell['atom']!=atom:continue
            source.append(sid);L,pure,start,count=(shell[k] for k in ('L','pure','offset','count'))
            indices.extend(range(start,start+count))
            for l in ([L] if pure else range(L,-1,-2)):
                ref=angular_reference(L,pure,l,rotation)
                for m in range(2*l+1):
                    vec=np.zeros(n);vec[start:start+count]=ref[:,m]
                    columns[l*l+m].append(vec);shell_ids[l*l+m].append(sid)
        centre=self.orthonormal(np.eye(n)[:,indices])
        angular=[self.orthonormal(np.stack(c,axis=1) if c else np.empty((n,0))) for c in columns]
        result=(centre,angular,source,shell_ids);self.workspace_cache[key]=result;return result

    def projection(self,value,indices,path):
        identity=json.dumps([value,indices],sort_keys=True,separators=(',',':'))
        if identity in self.projection_cache:return self.projection_cache[identity]
        self.counts['unique_projections']+=1
        assert value['status']=='available','Complete frozen FCHK supplies local AO definitions and metric'
        atom=value['atom_index'];rotation=np.asarray(value['rotation_reference_to_input']).reshape(3,3)
        self.check(path+'/proper-rotation',np.max(np.abs(rotation.T@rotation-np.eye(3)))<=CRITERIA['rotation_absolute_error'] and abs(np.linalg.det(rotation)-1)<=CRITERIA['rotation_absolute_error'])
        centre,angular,sources,shell_ids=self.workspace(atom,rotation)
        self.check(path+'/source-shells',value['source_shell_indices']==sources)
        expected_spins={name:[i for i in indices if self.reference_spins[i]==spin] for spin,name in enumerate(('alpha','beta'))}
        expected_spins={k:v for k,v in expected_spins.items() if v}
        self.check(path+'/spin-blocks',{b['spin']:b['orbital_indices'] for b in value['spins']}==expected_spins)
        total=0.;component_total=np.zeros(25);rank=0
        for block in value['spins']:
            bi=block['orbital_indices'];key=tuple(bi)
            if key not in self.target_cache:self.target_cache[key]=self.orthonormal(self.C[:,bi])
            target=self.target_cache[key];tr=float(np.linalg.norm(centre.T@target)**2);r=target.shape[1]
            bp=path+'/'+block['spin'];metrics=block['centre']
            self.check(bp+'/ranks',metrics['reference']['numerical_rank']==centre.shape[1] and metrics['orbitals']['numerical_rank']==r)
            self.close(bp+'/centre-trace',metrics['subspace_overlap_trace'],tr)
            if r:self.close(bp+'/centre-mean',metrics['mean_orbital_fraction'],tr/r)
            computed=[]
            self.check(bp+'/component-count',len(block['components'])==25)
            for i,entry in enumerate(block['components']):
                l=int(math.isqrt(i));cp=bp+'/component-'+str(i);ar=angular[i]
                metric=entry['metric'];ct=float(np.linalg.norm(ar.T@target)**2);computed.append(ct)
                self.check(cp+'/scope',entry['angular_degree']==l and entry['component']==i-l*l and entry['source_shell_indices']==shell_ids[i])
                self.check(cp+'/ranks',metric['reference']['numerical_rank']==ar.shape[1] and metric['orbitals']['numerical_rank']==r)
                self.close(cp+'/trace',metric['subspace_overlap_trace'],ct)
                if r:self.close(cp+'/mean',metric['mean_orbital_fraction'],ct/r)
                if tr>0:self.close(cp+'/local-fraction',entry['fraction_of_local_projection'],ct/tr)
            self.close(bp+'/partition-residual',block['angular_partition_residual'],sum(computed)-tr)
            total+=tr;rank+=r;component_total+=computed
        self.check(path+'/represented-rank',value['represented_spin_orbital_rank']==rank)
        self.close(path+'/total-centre',value['centre_projection_trace'],total)
        if rank:self.close(path+'/total-mean',value['centre_mean_fraction'],total/rank)
        for i,v in enumerate(component_total):self.close(path+'/total-component-'+str(i),value['component_projection_traces'][i],float(v))
        self.close(path+'/total-partition',value['angular_partition_residual'],float(sum(component_total)-total))
        result=(total,rank,component_total);self.projection_cache[identity]=result;return result

    def scope(self,s,path,expected_members=None):
        self.counts['scope_occurrences']+=1
        indices=s['orbital_indices'];origin=s['origin'];decomp=s['local_decomposition'];assignment=s['local_assignment']
        self.check(path+'/indices',len(set(indices))==len(indices) and all(0<=i<self.C.shape[1] for i in indices))
        if expected_members is not None:self.check(path+'/target-members',set(indices)==set(expected_members))
        if s['candidate_source']:
            self.check(path+'/propagated-candidate',origin in ('pi-partner-candidate','spin-counterpart-candidate') and assignment is None and decomp is None and s['molecular_assignment'] is None)
            self.scope(s['candidate_source'],path+'/candidate-source')
        if origin=='producer':
            self.check(path+'/producer-no-local-axes',not s['axes_available'] and s['rotation_reference_to_input'] is None and assignment is None and decomp is None)
            self.check(path+'/producer-label-retained',all(self.p['chemistry'][i]['symmetry']==s['label'] for i in indices))
        if origin=='molecular-operations' and s['molecular_assignment']:
            self.check(path+'/molecular-member-scope',set(indices)<=set(s['molecular_assignment']['orbital_indices']) and decomp is None and assignment is None)
        if decomp:
            self.check(path+'/local-axis-scope',s['axes_available'] and s['point_group_basis']=='local-coordination-template' and s['atom_indices'][0]==decomp['atom_index'] and s['rotation_reference_to_input']==decomp['rotation_reference_to_input'])
            total,rank,components=self.projection(decomp,indices,path+'/projection')
            if assignment:
                self.check(path+'/metric-assignment',origin=='local-metric-projection' and assignment['source']=='metric-angular-projection' and assignment['projection']==decomp and assignment['label']==s['label'])
                shell=assignment['shell'];weight=float(sum(components[shell*shell:(shell+1)**2]))
                self.close(path+'/assigned-centre-fraction',assignment['centre_mean_fraction'],total/rank)
                self.close(path+'/assigned-angular-fraction',assignment['angular_fraction_within_centre'],weight/total)
                confidence=assignment['conditional_shell_purity']
                self.check(path+'/conditional-range',.55<=confidence<=1+CRITERIA['projection_absolute_error'])
                self.close(path+'/assigned-target-fraction',assignment['labelled_fraction_of_target'],confidence*weight/rank)
                self.check(path+'/measure-scope',assignment['measure']=='unweighted-S-metric-subspace-overlap-not-electron-population')
        self.check(path+'/no-dimension-fallback-on-complete-input',origin!='local-dimension-candidate')

    def molecular_graph(self):
        # Independent block contraction of same-source density with GBasis S.
        ps=self.density['total_actual']@self.S;qs=self.density['spin_actual']@self.S
        products=ps*ps.T+qs*qs.T
        selectors=np.equal(np.arange(len(self.mol.atnums))[:,None],self.basis_atom[None,:]).astype(float)
        mayer=selectors@products@selectors.T
        errors=[abs(row[2]-mayer[int(row[0]),int(row[1])]) for row in self.p['bond_orders']]
        self.check('GRAPH/independent-Mayer-values',max(errors,default=0)<=CRITERIA['mayer_absolute_error'],
            {'recorded_pairs':len(errors),'maximum_absolute_error':float(max(errors,default=0))})
        distances=np.linalg.norm(self.mol.atcoords[:,None,:]-self.mol.atcoords[None,:,:],axis=2)
        radii=np.array([RADII[int(z)] for z in self.mol.atnums])*CRITERIA['angstrom_to_bohr']
        limit=CRITERIA['topology_electronic_radius_factor']*(radii[:,None]+radii[None,:])
        n=len(radii)
        edges={(a,b) for a in range(n) for b in range(a+1,n)
               if mayer[a,b]>=CRITERIA['topology_mayer_cutoff'] and distances[a,b]<=limit[a,b]}
        return n,edges

    def topology(self):
        n,expected_edges=self.molecular_graph();cycle_cache={};stats=[]
        for index,assignment in enumerate(self.p['pi_topology_assignments']):
            self.counts['topology_assignments']+=1
            path='topology/'+str(index);g=assignment['topology_graph'];channels=assignment['orientation_channels']
            if g is None:
                self.check(path+'/graph-available',False);continue
            edges={tuple(sorted(e)) for e in g['edges']}
            self.check(path+'/graph-contract',g['source']=='mayer-and-covalent-distance-model' and g['atom_count']==n and
                g['atom_index_base']==0 and g['channel_index_base']==0 and g['closed_cycle_paths'] and
                g['minimum_mayer_order']==CRITERIA['topology_mayer_cutoff'] and
                g['maximum_covalent_radius_factor']==CRITERIA['topology_electronic_radius_factor'])
            self.check(path+'/independent-edges',edges==expected_edges and len(edges)==len(g['edges']),
                {'missing':sorted(expected_edges-edges),'extra':sorted(edges-expected_edges)})
            adjacency=[set() for _ in range(n)]
            for a,b in edges:
                assert 0<=a<n and 0<=b<n and a!=b
                adjacency[a].add(b);adjacency[b].add(a)
            csets=[set(c['atoms']) for c in channels]
            for i,c in enumerate(csets):self.check(path+'/channel-'+str(i),len(c)==len(channels[i]['atoms']) and all(0<=a<n for a in c))
            degree_hubs=[i for i in range(n) if len(adjacency[i])>=4] if len(channels)>=2 else []
            self.check(path+'/examined-hubs',g['examined_hubs']==degree_hubs)
            witnessed=set();plain_hubs=set();all_cycles=[]
            if degree_hubs:
                key=tuple(sorted(edges))
                if key not in cycle_cache:cycle_cache[key]=enumerate_cycles(adjacency)
                all_cycles=cycle_cache[key]
                for hub in degree_hubs:
                    containing=[set(c)-{hub} for c in all_cycles if hub in c]
                    if any(a.isdisjoint(b) for a,b in itertools.combinations(containing,2)):plain_hubs.add(hub)
                for a,b in itertools.combinations(range(len(csets)),2):
                    for hub in plain_hubs:
                        first=[set(c)-{hub} for c in all_cycles if hub in c and associated(set(c)-{hub},csets[a],csets[b],hub)]
                        second=[set(c)-{hub} for c in all_cycles if hub in c and associated(set(c)-{hub},csets[b],csets[a],hub)]
                        if any(x.isdisjoint(y) for x in first for y in second):witnessed.add((a,b,hub))
            self.check(path+'/unscoped-cycle-unions',set(g['unscoped_cycle_union_hubs'])==plain_hubs)
            recorded=set()
            for widx,w in enumerate(g['channel_ring_witnesses']):
                wp=path+'/witness-'+str(widx);a,b=w['channel_indices'];hub=w['hub'];cycles=w['cycles']
                recorded.add((a,b,hub))
                self.check(wp+'/indices',0<=a<b<len(csets) and hub in degree_hubs)
                sets=[]
                for j,cycle in enumerate(cycles):
                    valid=len(cycle)>=4 and cycle[0]==cycle[-1]==hub and len(set(cycle[:-1]))==len(cycle)-1
                    valid=valid and all(tuple(sorted((x,y))) in edges for x,y in zip(cycle,cycle[1:]))
                    self.check(wp+'/closed-simple-path-'+str(j),valid)
                    sets.append(set(cycle)-{hub})
                self.check(wp+'/one-shared-hub',len(sets)==2 and sets[0].isdisjoint(sets[1]))
                self.check(wp+'/channel-boundary',associated(sets[0],csets[a],csets[b],hub) and associated(sets[1],csets[b],csets[a],hub))
                self.check(wp+'/additional-connection',w['additional_connection_without_hub']==connected_without_hub(adjacency,sets[0],sets[1],hub))
            self.check(path+'/complete-association-search',g['channel_association_search_complete'] and recorded==witnessed and len(recorded)==len(g['channel_ring_witnesses']),
                {'missing':sorted(witnessed-recorded),'extra':sorted(recorded-witnessed),'cycles_independently_enumerated':len(all_cycles)})
            if assignment['topology']=='spiro':
                self.check(path+'/spiro-supported',bool(g['channel_ring_witnesses']))
            stats.append({'family_id':assignment['family_id'],'topology':assignment['topology'],'channels':len(channels),'witnesses':len(recorded),'graph_edges':len(edges)})
        return stats

    def gaps(self):
        energies=np.asarray(self.mol.mo.energies)
        for key,kind in (('pi_interactions','pi-partner'),('crystal_field_gaps','crystal-field')):
            for i,g in enumerate(self.p[key]):
                self.counts['gap_records']+=1;path=key+'/'+str(i)
                self.check(path+'/typed-scope',g['gap_kind']==kind and g['interpretation_scope']=='local-metal-ligand-subspace' and g['score_meaning']=='heuristic-support-not-probability')
                means=[]
                for side in ('lower','upper'):
                    members=g[side+'_orbitals']
                    self.check(path+'/'+side+'-members',bool(members) and len(members)==len(set(members)) and all(0<=j<len(energies) for j in members))
                    mean=float(np.mean(energies[members]));means.append(mean)
                    self.close(path+'/'+side+'-energy',g[side+'_energy_hartree'],mean,CRITERIA['energy_absolute_error_hartree'])
                    self.scope(g[side+'_symmetry_scope'],path+'/'+side+'-scope',members)
                self.close(path+'/splitting',g['splitting_hartree'],means[1]-means[0],CRITERIA['energy_absolute_error_hartree'])
                self.check(path+'/energy-order',means[1]>=means[0])
                self.check(path+'/evidence-type',(g['orbital_evidence'] is not None and g['crystal_field_evidence'] is None) if kind=='pi-partner' else
                    (g['orbital_evidence'] is None and g['crystal_field_evidence'] is not None))

def associated(cycle,own,other,hub):
    return len(cycle & (own-other-{hub}))>=2 and not (cycle & (other-{hub}))

def connected_without_hub(adjacency,first,second,hub):
    reached=set(first);todo=list(first)
    while todo:
        a=todo.pop()
        for b in adjacency[a]:
            if b==hub or b in reached:continue
            if b in second:return True
            reached.add(b);todo.append(b)
    return False

def enumerate_cycles(adjacency):
    # Enumerate by minimum vertex, retain one orientation. This operates on
    # the complete graph, independently of COV's channel-pair path search.
    cycles=[];states=0
    def visit(start,path,seen):
        nonlocal states
        states+=1
        if states>CRITERIA['cycle_enumeration_state_limit']:raise RuntimeError('Independent complete-cycle enumeration budget exhausted; absence is unproven')
        for nxt in sorted(adjacency[path[-1]]):
            if nxt==start and len(path)>=3:
                if path[1]<path[-1]:cycles.append(tuple(path))
                if len(cycles)>CRITERIA['cycle_enumeration_count_limit']:raise RuntimeError('Independent cycle-count budget exhausted')
            elif nxt>start and nxt not in seen:visit(start,path+[nxt],seen|{nxt})
    for start in range(len(adjacency)):visit(start,[start],{start})
    return cycles

def self_checks():
    checks={}
    dummy=object.__new__(CaseReview);dummy.S=np.array([[1.,.6],[.6,1.]])
    dummy.root=metric_root(dummy.S)
    first=dummy.orthonormal(np.array([[1.],[0.]]));second=dummy.orthonormal(np.array([[0.],[1.]]))
    checks['nonorthogonal-known-trace']=abs(float(np.linalg.norm(first.T@second)**2)-.36)<2e-14
    checks['duplicate-and-zero-columns-rank']=dummy.orthonormal(np.array([[1.,2.,0.],[0.,0.,0.]])).shape[1]==1
    basis=np.array([[1.,2.],[0.,1.]])
    full=dummy.orthonormal(basis)
    checks['full-space-basis-change']=np.max(np.abs(full@full.T-np.eye(2)))<2e-14
    singular=np.ones((2,2));root=metric_root(singular)
    checks['semidefinite-metric-supported']=np.max(np.abs(root.T@root-singular))<2e-14
    scale=math.sqrt(3/(4*math.pi));expected=np.array([[0,-scale,0],[0,0,-scale],[scale,0,0]])
    checks['cartesian-p-phase-and-order']=np.max(np.abs(angular_reference(1,False,1,np.eye(3))-expected))<2e-14
    def graph(n,edges):
        a=[set() for _ in range(n)]
        for i,j in edges:a[i].add(j);a[j].add(i)
        return a
    edges=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
    g=graph(5,edges);cycles=enumerate_cycles(g)
    checks['two-ring-positive']=set(map(frozenset,cycles))=={frozenset([0,1,2]),frozenset([0,3,4])}
    checks['exclusive-channel-positive']=associated({1,2},{1,2},{3,4},0)
    checks['other-channel-atom-rejected']=not associated({1,2,3},{1,2},{3,4},0)
    checks['hub-only-connection']=not connected_without_hub(g,{1,2},{3,4},0)
    bridged=graph(5,edges+[(2,4)])
    checks['additional-connection']=connected_without_hub(bridged,{1,2},{3,4},0)
    checks['acyclic-negative']=not enumerate_cycles(graph(5,[(0,1),(1,2),(2,3),(2,4)]))
    return {'checks':{k:bool(v) for k,v in checks.items()},'all_passed':all(checks.values())}

def review_case(root,numeric_root,case,out,manifest,review_manifest_sha):
    cid=case['case_id'];terminal=read(root/'cases'/cid/'terminal.json');nt=read(numeric_root/'cases'/cid/'terminal.json')
    assert terminal['round_identity']==nt['round_identity']==manifest['round_identity']
    assert nt['review_manifest_sha256']==review_manifest_sha
    directory=root/terminal['evidence_directory'];inventory=root/terminal['evidence_manifest']
    assert sha(inventory)==terminal['evidence_manifest_sha256']
    expected=read(inventory);assert sha(directory/'production.json')==expected['production.json']['sha256']
    source=Path(case['input']);assert sha(source)==case['sha256']
    attempt=numeric_root/nt['attempt'];result=attempt/'result'
    dependencies={}
    for name in ('independent-metric.npz','independent-density.npz','review.json'):
        p=result/name;entry=nt['artifacts'][str(Path('result')/name)]
        assert sha(p)==entry['sha256'],name;dependencies[str(p)]=entry['sha256']
    nr=read(result/'review.json');assert nr['case_id']==cid and nr['round_identity']==manifest['round_identity']
    reference=np.load(result/'independent-metric.npz',allow_pickle=False)
    density=np.load(result/'independent-density.npz',allow_pickle=False)
    indices=reference['source_indices'];scales=reference['basis_scales']
    S=reference['overlap'][np.ix_(indices,indices)]*scales[:,None]*scales[None,:]
    C=reference['source_coefficients'][indices]*scales[:,None]
    mol=load_one(str(source),fmt='fchk');production=read(directory/'production.json')
    review=CaseReview(mol,S,C,production,density)
    review.check('NUMERICAL-REFERENCE-REVIEW-PASSED',nt['status']=='numeric_subset_pass',nt['status'])
    for i,row in enumerate(production['compact_rows']):review.scope(row['symmetry_explanation'],'row/'+str(i),row['members'])
    topology=review.topology();review.gaps()
    if cid=='OLD-071':review.check('NAMED/allene-not-spiro',all(x['topology']!='spiro' for x in topology))
    if cid=='OLD-107':review.check('NAMED/spiropentadiene-positive-witness',any(x['topology']=='spiro' and x['witnesses']>0 for x in topology))
    failures=[c for c in review.checks if c['status']!='pass']
    return {'case_id':cid,'round_identity':manifest['round_identity'],'status':'scope_subset_fail' if failures else 'scope_subset_pass',
        'check_count':len(review.checks),'failure_count':len(failures),'counts':review.counts,'topology':topology,
        'maximum_projection_comparison_error':review.maximum_projection_error,'failures':failures,'checks':review.checks,
        'reused_independent_references':dependencies,'production_sha256':sha(directory/'production.json'),
        'formal_complete_case_pass':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--round',type=Path,required=True)
    parser.add_argument('--numeric-review',default='review-v1');parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.round;numeric=root/args.numeric_review;out=args.output
    manifest=read(root/'manifest.json');barrier=read(root/'collection-complete.json');nb=read(numeric/'review-complete.json')
    assert barrier['all_terminal'] and barrier['terminal_cases']==manifest['case_count'] and barrier['round_identity']==manifest['round_identity']
    assert nb['all_terminal'] and nb['terminal_cases']==manifest['case_count']
    review_manifest_sha=sha(numeric/'manifest.json')
    assert read(numeric/'manifest.json')['round_identity']==manifest['round_identity']
    out.mkdir(parents=True,exist_ok=False);start=time.time();results=[]
    write(out/'criteria.json',{**CRITERIA,'structural_radii_angstrom':RADII})
    identity={'round_identity':manifest['round_identity'],'collection_manifest_sha256':sha(root/'manifest.json'),
        'collection_barrier_sha256':sha(root/'collection-complete.json'),'numeric_barrier_sha256':sha(numeric/'review-complete.json'),
        'numeric_environment_manifest_sha256':review_manifest_sha,'reviewer_sha256':sha(__file__),
        'criteria_sha256':sha(out/'criteria.json'),'started_epoch':start,'formal_complete_case_passes':0}
    write(out/'identity.json',identity)
    controls=self_checks();write(out/'independent-oracle-controls.json',controls)
    if not controls['all_passed']:raise RuntimeError('Independent oracle controls failed')
    for case in manifest['cases']:
        t=time.perf_counter()
        try:record=review_case(root,numeric,case,out,manifest,review_manifest_sha)
        except Exception as error:
            record={'case_id':case['case_id'],'status':'review_error','error':type(error).__name__+': '+str(error),'traceback':traceback.format_exc(),'formal_complete_case_pass':False}
        record['wall_seconds']=time.perf_counter()-t
        dest=out/(case['case_id']+'.json');write(dest,record)
        results.append({'case_id':case['case_id'],'status':record['status'],'report':dest.name,'sha256':sha(dest),'wall_seconds':record['wall_seconds'],
            'failure_count':record.get('failure_count'), 'check_count':record.get('check_count')})
        write(out/'progress.json',{'terminal_cases':len(results),'expected_cases':manifest['case_count'],'counts':{s:sum(x['status']==s for x in results) for s in sorted({x['status'] for x in results})}})
        print(json.dumps(results[-1]),flush=True)
    counts={s:sum(x['status']==s for x in results) for s in sorted({x['status'] for x in results})}
    write(out/'summary.json',{'all_terminal':True,'all_passed':counts.get('scope_subset_pass')==manifest['case_count'],
        'round_identity':manifest['round_identity'],'identity_sha256':sha(out/'identity.json'),'criteria_sha256':sha(out/'criteria.json'),
        'case_count':len(results),'status_counts':counts,'cases':results,'wall_seconds':time.time()-start,
        'formal_complete_case_passes':0,'limitations':[CRITERIA['scope'],'Conditional label component fractions are checked; point-group catalogue physics, local-axis suitability and full molecular symmetry remain separate acceptance obligations.']})
    return 0 if counts.get('scope_subset_pass')==manifest['case_count'] else 2

if __name__=='__main__':raise SystemExit(main())
