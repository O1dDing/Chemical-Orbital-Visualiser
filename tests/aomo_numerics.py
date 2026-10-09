"""Independent AO/NBO numerical acceptance; never starts chemical producers.

Vendor W matrices, archive density/metric/canonical columns, literal report
tables and FCHK arrays are the oracle. Production JSON is an observation.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, shutil, subprocess, sys, time, traceback
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '2'
sys.path.insert(0, os.environ.get('COV_REFERENCE_DEPS', r'E:/Dev/cov-validation-20260905/reference-deps'))
import numpy as np

NUMBER = r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][+-]?\d+)?'
HEADINGS = {'AONAO':'NAOs in the AO basis:', 'AONBO':'NBOs in the AO basis:',
 'AONHO':'NHOs in the AO basis:', 'AONLMO':'NLMOs in the AO basis:', 'AOPNAO':'PNAOs in the AO basis:',
 'NAOMO':'MOs in the NAO basis:', 'NBOMO':'MOs in the NBO basis:', 'NLMOMO':'MOs in the NLMO basis:',
 'AOMO':'MOs in the AO basis:', 'NAONBO':'NBOs in the NAO basis:', 'NAONHO':'NHOs in the NAO basis:',
 'NAONLMO':'NLMOs in the NAO basis:', 'NHONBO':'NBOs in the NHO basis:', 'NBONLMO':'NLMOs in the NBO basis:'}
LFN = dict(zip(('AONBO','NBOMO','SAO','NAOMO','AONAO','NAONBO','AONHO','AONLMO','AOPNAO','NAONHO','NAONLMO','NHONBO','NBONLMO','NLMOMO','AOMO'),(37,49,50,51,52,53,54,55,56,57,58,59,60,61,62)))
def save(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
def values(text): return np.array([float(x.replace('D','E').replace('d','e')) for x in text.split()])
def error(a,b): return float(np.max(np.abs(a-b))) if a.size else 0.
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def fchk(path):
    lines=Path(path).read_text(errors='replace').splitlines(); out={};i=0
    while i<len(lines):
        match=re.match(r'^(.{40})\s+([IR])\s+(?:N=\s*(\d+)|(.+))$',lines[i]);i+=1
        if not match:continue
        name,typ,count,scalar=match.groups();name=name.strip()
        if count:
            words=[]
            while len(words)<int(count):words+=lines[i].split();i+=1
            out[name]=np.asarray([float(x.replace('D','E')) for x in words])
        else:out[name]=float(scalar)
    return out
def archive(path):
    text=Path(path).read_text();blocks={k.upper():v for k,v in re.findall(r'\$([A-Za-z0-9]+)\s*(.*?)\$END',text,re.S)}
    head=blocks['GENNBO'];n=int(re.search(r'NBAS\s*=\s*(\d+)',head).group(1));opened=bool(re.search(r'\bOPEN\b',head));spins=['alpha','beta'] if opened else ['total'];out={'n':n,'spins':spins,'blocks':blocks}
    for key in ('OVERLAP','DENSITY','FOCK','LCAOMO'):
        if key not in blocks:continue
        raw=values(blocks[key]); tags=['total'] if key=='OVERLAP' else spins;count=n*n if key=='LCAOMO' else n*(n+1)//2
        if raw.size!=count*len(tags):raise ValueError(f'{key} count {raw.size} != {count*len(tags)}')
        out[key]={}
        for i,spin in enumerate(tags):
            piece=raw[i*count:(i+1)*count]
            if key=='LCAOMO':m=piece.reshape((n,n),order='F')
            else:m=np.zeros((n,n));m[np.tril_indices(n)]=piece;m=m+np.tril(m,-1).T
            out[key][spin]=m
    coords=[]
    for line in blocks['COORD'].splitlines():
        row=line.split()
        if len(row)==5 and row[0].isdigit():coords.append([float(x) for x in row])
    out['atoms']=coords;return out
def identity(label):
    m=re.match(r'\s*([A-Za-z0-9]+\*?)\s*\(\s*(\d+)\s*\)\s*(.*)',label)
    return (m[1],int(m[2]),tuple((x,int(y)) for x,y in re.findall(r'([A-Z][a-z]?)\s+(\d+)',m[3]))) if m else None
def tables(path):
    result={k:{} for k in ('nao','npa','nbo','nlmo','wiberg')};spin='total';section='';cols=[]
    for line in Path(path).read_text(errors='replace').splitlines():
        if 'ALPHA SPIN ORBITALS' in line.upper():spin='alpha';section=''
        if 'BETA  SPIN ORBITALS' in line.upper() or 'BETA SPIN ORBITALS' in line.upper():spin='beta';section=''
        if 'NATURAL POPULATIONS:' in line:section='nao'
        if 'Summary of Natural Population Analysis:' in line:section='npa'
        if 'Wiberg bond index matrix' in line:section='wiberg'
        if 'Wiberg bond index, Totals' in line:section=''
        if '(Occupancy)' in line and 'Bond orbital' in line:section='nbo'
        if 'SECOND ORDER PERTURBATION' in line or 'NATURAL BOND ORBITALS (Summary)' in line:section=''
        if 'NLMO / Occupancy / Percent from Parent NBO' in line:section='nlmo'
        if section=='nao':
            m=re.match(rf'^\s*(\d+)\s+([A-Za-z]+)\s+(\d+)\s+(\S+)\s+([A-Za-z]+\([^)]*\))\s+({NUMBER})',line)
            if m:result['nao'].setdefault(spin,{})[int(m[1])]=dict(atom=int(m[3])-1,angular=m[4],type=m[5],occupation=float(m[6]))
        if section=='npa':
            m=re.match(rf'^\s*([A-Z][a-z]?)\s+(\d+)\s+((?:{NUMBER}\s+){{4}}{NUMBER})(?:\s+({NUMBER}))?\s*$',line)
            if m:result['npa'].setdefault(spin,{})[int(m[2])-1]=dict(zip(('charge','core','valence','rydberg','total'),values(m[3]).tolist()))
        if section=='nbo':
            m=re.match(rf'^\s*(\d+)\.\s*\(({NUMBER})\)\s*(.*)$',line)
            if m:result['nbo'].setdefault(spin,{})[int(m[1])]=dict(occupation=float(m[2]),identity=identity(m[3]))
        if section=='nlmo':
            m=re.match(rf'^\s*(\d+)\.\s*\(({NUMBER})\)\s*({NUMBER})%\s*(.*)$',line)
            if m:result['nlmo'].setdefault(spin,{})[int(m[1])]=dict(occupation=float(m[2]),percent=float(m[3]),identity=identity(m[4]))
        if section=='wiberg':
            if re.match(r'^\s*Atom\s+\d',line):cols=[int(x)-1 for x in re.findall(r'\d+',line)]
            m=re.match(r'^\s*(\d+)\.\s+[A-Z][a-z]?\s+(.*)$',line)
            if m:
                row=values(m[2]);assert len(row)==len(cols)
                for col,val in zip(cols,row):result['wiberg'].setdefault(spin,{})[(int(m[1])-1,col)]=float(val)
    return result
def read_w(path,kind,rows,cols):
    lines=Path(path).read_text().splitlines();assert lines[1].strip()==HEADINGS[kind],(path,lines[:3]);starts={}
    for i,line in enumerate(lines):
        if line.strip()=='ALPHA SPIN':starts['alpha']=i+1
        if line.strip() in ('BETA SPIN','BETA  SPIN'):starts['beta']=i+1
    if not starts:starts={'total':3}
    out={}
    for spin,first in starts.items():
        vals=[]
        for line in lines[first:]:
            words=line.split()
            if not words:continue
            if not all(re.fullmatch(NUMBER,x) and any(c in x for c in '.EeDd') for x in words):break
            vals.extend(float(x.replace('D','E')) for x in words)
            if len(vals)>=rows*cols:break
        if len(vals)!=rows*cols:raise ValueError(f'{kind}/{spin}: expected {rows}x{cols}, got {len(vals)} before footer')
        out[spin]=np.asarray(vals).reshape((rows,cols),order='F')
    return out
def reference_field_columns(case):
    """Return real local columns in IOData/GBasis FCHK basis order.

    The mapping comes solely from FILE47 LABEL angular polynomials and ordered
    primitive shells, using NBO7 manual B.7 pp. B-74/B-75 and IOData's FCHK
    convention table. No production coefficients/row map or observed field is
    used. This remains valid for rectangular local spaces and does not replace
    a local orbital by its projection onto the retained canonical space.

    Result: basis (GBasis), molecule (IOData), families[kind][spin] as AO x
    local-column arrays, canonical[spin], phase[spin], overlap in that same AO
    order, links[kind][spin] as local x canonical, and checks metadata. Spin is
    'total' for restricted inputs and 'alpha'/'beta' for unrestricted inputs.
    """
    from iodata import load_one
    from iodata.formats.fchk import CONVENTIONS
    from gbasis.wrappers import from_iodata
    if 'archive47' in case:case=dict(id=case.get('case_id'),fchk=case['fchk'],archive=case['archive47'],report=case['report'],matrices={k:case[k.lower()] for k in LFN if k.lower() in case})
    ar=archive(case['archive']);raw=fchk(case['fchk']);n=ar['n'];mol=load_one(case['fchk']);basis=from_iodata(mol)
    if 'CUBICF' in ar['blocks']['GENNBO'].upper():raise ValueError('Cubic spherical f convention requires its explicit polynomial mapping')
    def fields(block):
        matches=list(re.finditer(r'([A-Za-z][A-Za-z0-9]*)\s*=',block));return {m[1].upper():values(block[m.end():matches[i+1].start() if i+1<len(matches) else len(block)]) for i,m in enumerate(matches)}
    meta=fields(ar['blocks']['BASIS']);contract=fields(ar['blocks']['CONTRACT']);labels=meta['LABEL'].astype(int);centers=meta['CENTER'].astype(int);components=contract['NCOMP'].astype(int);mapping=[];offset=0;primitive_offset=0
    shell_types=raw['Shell types'].astype(int);shell_atoms=raw['Shell to atom map'].astype(int);primitives=raw['Number of primitives per shell'].astype(int)
    if len(shell_types)!=len(components):raise ValueError('Archive/FCHK primitive shell count mismatch')
    def angular_name(label):
        if label in (1,51):return '1'
        l,index=divmod(label,100)
        if l==1:return ['x','y','z'][(index-1)%50]
        if index<50:
            names=['x'*x+'y'*y+'z'*(l-x-y) for x in range(l,-1,-1) for y in range(l-x,-1,-1)]
            return names[index-1]
        if l==2:return {51:'s2',52:'c1',53:'s1',54:'c2',55:'c0'}[index]
        return (['c0']+[v for m in range(1,l+1) for v in (f'c{m}',f's{m}')])[index-51]
    for sh,(kind,count) in enumerate(zip(shell_types,components)):
        names=['1','x','y','z'] if kind==-1 else CONVENTIONS[(abs(int(kind)),'p' if kind<=-2 else 'c')]
        if len(names)!=count or np.any(centers[offset:offset+count]!=shell_atoms[sh]):raise ValueError('Archive/FCHK ordered atom/angular shell mismatch')
        nprim=int(primitives[sh]);start=int(contract['NPTR'][sh])-1
        if int(contract['NPRIM'][sh])!=nprim:raise ValueError('Primitive contraction length mismatch')
        if not np.allclose(contract['EXP'][start:start+nprim],raw['Primitive exponents'][primitive_offset:primitive_offset+nprim],rtol=1e-7,atol=1e-10):raise ValueError('Primitive exponents mismatch')
        for l in ([0,1] if kind==-1 else [abs(int(kind))]):
            coeff=raw['P(S=P) Contraction coefficients'] if kind==-1 and l==1 else raw['Contraction coefficients']
            if not np.allclose(contract[['CS','CP','CD','CF','CG'][l]][start:start+nprim],coeff[primitive_offset:primitive_offset+nprim],rtol=2e-7,atol=1e-10):raise ValueError('Primitive contraction coefficients mismatch')
        local=[names.index(angular_name(label)) for label in labels[offset:offset+count]]
        if len(set(local))!=count:raise ValueError('Angular labels are not a unique shell permutation')
        mapping.extend(offset+j for j in local);offset+=count;primitive_offset+=nprim
    if offset!=n or len(set(mapping))!=n:raise ValueError('Incomplete archive AO label map')
    fcoords=raw['Current cartesian coordinates'].reshape((-1,3));acoords=np.array(ar['atoms'])[:,2:];factor=1. if 'BOHR' in ar['blocks']['GENNBO'].upper() else 1.8897261254578281
    if error(fcoords,acoords*factor)>3e-5:raise ValueError('Archive/FCHK atom coordinates mismatch')
    if not np.array_equal(raw['Atomic numbers'].astype(int),np.array(ar['atoms'])[:,0].astype(int)):raise ValueError('Archive/FCHK element order mismatch')
    mapping=np.array(mapping);s=np.zeros((n,n));s[np.ix_(mapping,mapping)]=ar['OVERLAP']['total'];canonical={};phase={};checks={};families={};links={};tab=tables(case['report']);r=len(tab['nao'][ar['spins'][0]])
    for spin in ar['spins']:
        key=('Beta' if spin=='beta' else 'Alpha')+' MO coefficients'
        if key not in raw:raise ValueError(f'Missing direct FCHK {spin} canonical columns')
        canonical[spin]=raw[key].reshape((-1,n)).T;nc=canonical[spin].shape[1];archive_c=np.zeros((n,nc));archive_c[mapping,:]=ar['LCAOMO'][spin][:,:nc];dots=np.diag(canonical[spin].T@s@archive_c);phase[spin]=np.where(dots<0,-1.,1.);residual=error(canonical[spin],archive_c*phase[spin]);metric=error(canonical[spin].T@s@canonical[spin],np.eye(nc))
        if residual>3e-6 or metric>2e-5 or np.min(abs(dots))<.9999:raise ValueError(f'{spin}: metadata AO map failed direct coefficient/metric identity: {residual}/{metric}')
        checks[spin]=dict(phase_aligned_coefficient_error=residual,metric_error=metric,minimum_abs_column_overlap=float(min(abs(dots))))
    # IOData's coefficients must agree with the directly read FCHK ordering.
    io_coeff=np.asarray(mol.mo.coeffs);joined=np.hstack(list(canonical.values()))
    if io_coeff.shape!=joined.shape or error(io_coeff,joined)>1e-12:raise ValueError('IOData convention differs from its declared FCHK basis order')
    for kind in ('NAO','NBO','NHO','NLMO','PNAO'):
        blocks=read_w(case['matrices']['AO'+kind],'AO'+kind,n,r);families[kind]={};links[kind]={}
        for spin in ar['spins']:
            b=blocks.get(spin,blocks.get('total'));ordered=np.zeros_like(b);ordered[mapping,:]=b;families[kind][spin]=ordered
            if kind!='PNAO':links[kind][spin]=ordered.T@s@canonical[spin]
    return dict(basis=basis,molecule=mol,families=families,canonical=canonical,phase=phase,overlap=s,links=links,checks=checks,archive_to_fchk_ao=mapping.tolist(),spins=ar['spins'])
def collect(root,ids):
    root=Path(root);found=[]
    for case in sorted(root.iterdir()):
        if not case.is_dir() or (ids and case.name not in ids):continue
        choices=list(case.glob('*/manifest.json'));chosen=None
        for p in choices:
            m=json.loads(p.read_text(encoding='utf-8-sig'))
            if m.get('status')=='producer_complete_unvalidated' and m.get('paired_fchk') and m.get('matrix_archive47'):chosen=m
        stage=case/'aomo-from-two-stage-01'/'manifest.json'
        if stage.exists():
            candidate=json.loads(stage.read_text(encoding='utf-8-sig'))
            if candidate.get('status')=='producer_complete_unvalidated':
                chosen=candidate;chosen['paired_fchk']=candidate['canonical_fchk'];chosen['matrix_archive47']=str(stage.parent/'FILE.47')
        if chosen:found.append(dict(id=case.name,fchk=chosen['paired_fchk'],archive=chosen['matrix_archive47'],report=chosen['report'],matrices=chosen['matrices'],provenance=chosen.get('fchk_relationship',chosen.get('kind'))))
        elif ids:found.append(dict(id=case.name,unavailable='No complete paired manifest yet'))
    return found
def audit(case,out,probe,reuse=None):
    out.mkdir(parents=True);checks=[];start=time.perf_counter()
    def check(key,ok,observed,tolerance=None):checks.append(dict(id=key,status='pass' if ok else 'fail',observed=observed,tolerance=tolerance))
    if 'unavailable' in case:return dict(id=case['id'],status='unavailable',reason=case['unavailable'])
    inputs={case['fchk'],case['archive'],case['report'],*case['matrices'].values()};input_hashes={str(p):sha(p) for p in sorted(inputs)}
    a=archive(case['archive']);t=tables(case['report']);n=a['n'];s=a['OVERLAP']['total'];r=len(t['nao'][a['spins'][0]]);local={'NAO':r,'NHO':r,'NBO':r,'NLMO':r,'PNAO':r,'AO':n};active={spin:int(np.count_nonzero(np.max(abs(c),axis=0))) for spin,c in a['LCAOMO'].items()};local['MO']=max(active.values())
    mats={};errors={};density={};norms={};raw_summary={'ao_dimension':n,'local_dimension':r,'canonical_dimensions':active,'chains':{},'basis':{}}
    for kind,path in case['matrices'].items():
        if kind not in HEADINGS:continue
        origin=next(x for x in ('NLMO','NAO','NHO','NBO','AO') if kind.startswith(x));target=kind[len(origin):];mats[kind]=read_w(path,kind,local[origin],local[target])
    def get(kind,spin):return mats[kind].get(spin,mats[kind].get('total'))
    for spin in a['spins']:
        c=a['LCAOMO'][spin][:,:active[spin]];p=a['DENSITY'][spin];sp=s@p@s
        for family in ('NAO','NHO','NBO','NLMO'):
            b=get('AO'+family,spin);g=b.T@s@b;orth=error(g,np.eye(r));check(f'{spin}:{family}:metric',orth<2e-5,orth,2e-5);raw_summary['basis'][spin+':'+family]={'orthogonality':orth,'columns':r}
            if family in ('NAO','NBO','NLMO'):
                pop=np.diag(b.T@sp@b);table=t[family.lower()].get(spin,{})
                pe=max((abs(pop[i-1]-row['occupation']) for i,row in table.items()),default=0);check(f'{spin}:{family}:occupation',pe<4e-5,pe,4e-5)
            if family=='NAO':density[spin]=b.T@sp@b
        for kind in ('NAOMO','NBOMO','NLMOMO','NAONBO','NAONHO','NAONLMO','NHONBO','NBONLMO'):
            origin=next(x for x in ('NLMO','NAO','NHO','NBO') if kind.startswith(x));target=kind[len(origin):];left=get('AO'+origin,spin);x=get(kind,spin);right=c if target=='MO' else get('AO'+target,spin);projection=left.T@s@right;pr=error(projection,x);check(f'{spin}:{kind}:projection',pr<2e-5,pr,2e-5);delta=left@x-right;metric_error=float(np.sqrt(max(0,float(np.max(np.diag(delta.T@s@delta))))));raw_summary['chains'][spin+':'+kind]={'projection_error':pr,'reconstruction_metric_max':metric_error,'full_space':r>=active[spin] if target=='MO' else True};check(f'{spin}:{kind}:composition',metric_error<2e-5 or (target=='MO' and r<active[spin]),metric_error,2e-5)
        nh=get('NAONHO',spin);labels=t['nao'][spin];atom=np.array([labels[i+1]['atom'] for i in range(r)]);weights=np.stack([np.sum(nh[atom==j,:]**2,axis=0) for j in range(len(a['atoms']))]);loc=error(weights.max(axis=0),np.ones(r));check(spin+':NHO:locality',loc<2e-5,loc,2e-5)
        hybrid={letter:np.sum(nh[[labels[i+1]['angular'].startswith(letter) for i in range(r)],:]**2,axis=0).tolist() for letter in 'spdfg'};raw_summary['basis'][spin+':NHO']['angular_weights']=hybrid
        bl=get('NBONLMO',spin);parents={row['identity']:i for i,row in t['nbo'][spin].items()};parent_errors=[]
        for i,row in t['nlmo'].get(spin,{}).items():
            parent=parents.get(row['identity']);parent_errors.append(abs(100*bl[parent-1,i-1]**2-row['percent']) if parent else 1e9)
        pe=max(parent_errors,default=0);check(spin+':NLMO:printed_parent',pe<.0003,pe,.0003)
    for spin,rows in t['npa'].items():
        relevant=a['spins'] if spin=='total' else [spin];population_error=0
        for atom,row in rows.items():
            explicit=sum(sum(density[tag][i-1,i-1] for i,label in t['nao'][tag].items() if label['atom']==atom) for tag in relevant);z,zcore,*_=a['atoms'][atom];factor=1 if spin=='total' else .5;pop=explicit+(z-zcore)*factor;population_error=max(population_error,abs(pop-row['total']),abs(z*factor-pop-row['charge']))
        check(spin+':NPA:population_charge',population_error<4e-5,population_error,4e-5)
    for spin,rows in t['wiberg'].items():
        relevant=a['spins'] if spin=='total' else [spin];we=0
        for (i,j),value in rows.items():
            if i==j:continue
            predicted=0
            for tag in relevant:
                atoms=np.array([t['nao'][tag][k+1]['atom'] for k in range(r)]);predicted+=(1 if tag=='total' else 2)*np.sum(density[tag][np.ix_(atoms==i,atoms==j)]**2)
            we=max(we,abs(predicted-value))
        check(spin+':Wiberg:density',we<8e-5,we,8e-5)
    save(out/'raw-oracle.json',raw_summary)
    if probe or reuse:
        production_path=Path(reuse)/case['id']/'production.json' if reuse else out/'production.json'
        if reuse:code=0 if production_path.is_file() else 1
        else:
            proc=subprocess.run([str(probe),case['fchk'],str(Path(case['report']).parent),str(production_path)],capture_output=True,text=True,timeout=900);(out/'probe.log').write_text(proc.stdout+proc.stderr);code=proc.returncode
        check('production_exit',code==0,dict(returncode=code,reused=bool(reuse),evidence=str(production_path)))
        if code==0:
            observed=json.loads(production_path.read_text());d=observed['integration'];caps={x['key']:x for x in d['capabilities']};check('canonical_immutable',observed['canonical_preserved'],observed['canonical_preserved'])
            for kind,spinblocks in mats.items():
                for spin,ref in spinblocks.items():
                    matches=[m for m in d['dataset']['matrices'] if m['kind']==kind and m['spin']==spin];err=error(np.array(matches[0]['values']).reshape(ref.shape),ref) if len(matches)==1 and (matches[0]['rows'],matches[0]['columns'])==ref.shape else None;check('production_raw:'+kind+':'+spin,err is not None and err<1e-12,err,1e-12)
            for kind in ('nao','nho','nbo','nlmo'):check('production_capability:'+kind,caps[kind]['state']=='available',caps[kind])
            f=fchk(case['fchk']);fc_errors=[]
            for mo in observed['canonical']['orbitals']:
                key=('Beta' if mo['spin']=='beta' else 'Alpha')+' MO coefficients';coeff=f[key].reshape((-1,n))[mo['source_index']];fc_errors.append(error(coeff,np.asarray(mo['gaussian_coefficients'])))
            check('canonical_coefficients_vs_FCHK',max(fc_errors,default=0)<1e-12,max(fc_errors,default=0),1e-12)
            mapping=d['dataset']['association']['gaussian_row_zero_based'];scale=d['dataset']['association']['coefficient_scale'];trans=observed['canonical']['gaussian_transform'];inverse={x['source_index']:i for i,x in enumerate(trans)};coeff_errors=[]
            field_reference=reference_field_columns(case)
            check('independent_metadata_AO_mapping',mapping==field_reference['archive_to_fchk_ao'] and max(abs(np.asarray(scale)-1))<3e-6,dict(mapping_matches=mapping==field_reference['archive_to_fchk_ao'],scale_max_error=float(max(abs(np.asarray(scale)-1))),source_checks=field_reference['checks']),3e-6)
            canonical_raw={};phase_errors=[];link_errors=[]
            for spin in a['spins']:
                key=('Beta' if spin=='beta' else 'Alpha')+' MO coefficients';cgauss=f[key].reshape((-1,n)).T;cnbo=cgauss[np.asarray(mapping),:]/np.asarray(scale)[:,None];canonical_raw[spin]=cnbo;reference=a['LCAOMO'][spin][:,:cnbo.shape[1]];phase=np.where(np.sum(cnbo*reference,axis=0)<0,-1.,1.);phase_errors.append(error(cnbo,reference*phase));
            check('independent_archive_vs_FCHK_columns',max(phase_errors,default=1e9)<2e-5,max(phase_errors,default=1e9),2e-5)
            expected_links={}
            for spin,cnbo in canonical_raw.items():
                for kind in ('NAO','NHO','NBO','NLMO'):expected_links[(spin,kind)]=field_reference['links'][kind][spin]
            for link in d['links']:
                ref=link['orbital'];kind=ref['kind'];spin=ref['spin'];mo=observed['canonical']['orbitals'][link['canonical_index']]
                if (spin,kind) not in expected_links:continue
                expected=expected_links[(spin,kind)][ref['index'],mo['source_index']];link_errors.append(abs(link['coefficient']-expected));
                if link['weight'] is not None:link_errors.append(abs(link['weight']-expected**2))
            check('production_signed_links_vs_raw_projection',max(link_errors,default=0)<2e-5,max(link_errors,default=0),2e-5)
            observed_link_count=sum(1 for link in d['links'] if link['orbital']['kind'] in ('NAO','NHO','NBO','NLMO'))
            expected_link_count=4*r*sum(c.shape[1] for c in canonical_raw.values())
            check('production_complete_link_coverage',observed_link_count==expected_link_count,dict(observed=observed_link_count,expected=expected_link_count))
            for orbital in d['orbitals']:
                kind=orbital['ref']['kind'];spin=orbital['ref']['spin'];col=orbital['ref']['index']
                if 'AO'+kind not in mats:continue
                ref=get('AO'+kind,spin)[:,col];internal=np.zeros(n)
                for k,g in enumerate(mapping):idx=inverse[g];internal[idx]=ref[k]*scale[k]*trans[idx]['scale']
                coeff_errors.append(error(internal,np.array(orbital['coefficients'])))
            check('production_all_local_coefficients',max(coeff_errors,default=0)<1e-12,max(coeff_errors,default=0),1e-12)
            reconstruction=[row['error'] for row in observed['reconstruction'] if row['kind']=='GaussianAO' or r>=max(active.values())];check('production_typed_reconstruction',all(x is not None and x<2e-5 for x in reconstruction),max((x for x in reconstruction if x is not None),default=0),2e-5)
    check('input_bytes_preserved',all(sha(p)==value for p,value in input_hashes.items()),input_hashes)
    status='pass' if all(x['status']=='pass' for x in checks) else 'fail';result=dict(id=case['id'],status=status,seconds=time.perf_counter()-start,dimensions={'AO':n,'localized':r,'MO':active},provenance=case.get('provenance'),checks=checks);save(out/'result.json',result);return result
def negative_cases(case,foreign,out,probe):
    results=[]
    originals={p:sha(p) for p in [case['fchk'],case['archive'],case['report'],*case['matrices'].values()]}
    for name,kind in [('corrupt_nho','AONHO'),('corrupt_nlmo','AONLMO'),('corrupt_redundant_chain','NAONLMO'),('missing_nlmo','AONLMO'),('wrong_source','archive'),('ambiguous_steps',None)]:
        dest=out/name;inputs=dest/'input';inputs.mkdir(parents=True);shutil.copy2(case['archive'],inputs/'FILE.47');shutil.copy2(case['report'],inputs/'analysis.log')
        for matrix,path in case['matrices'].items():shutil.copy2(path,inputs/f'FILE.{LFN[matrix]}')
        if name.startswith('corrupt'):
            path=inputs/f'FILE.{LFN[kind]}';lines=path.read_text().splitlines()
            for i,line in enumerate(lines[3:],3):
                if not re.fullmatch(rf'\s*(?:{NUMBER}\s*)+',line):continue
                match=re.search(NUMBER,line);lines[i]=line[:match.start()]+f'{float(match[0])+.125:15.9f}'+line[match.end():];break
            path.write_text('\n'.join(lines)+'\n')
        if name=='missing_nlmo':(inputs/f'FILE.{LFN[kind]}').unlink()
        if name=='wrong_source':shutil.copy2(foreign['archive'],inputs/'FILE.47')
        if name=='ambiguous_steps':shutil.copytree(inputs,dest/'other-step')
        observed_path=dest/'production.json';proc=subprocess.run([str(probe),case['fchk'],str(dest if name=='ambiguous_steps' else inputs),str(observed_path)],capture_output=True,text=True,timeout=120);(dest/'probe.log').write_text(proc.stdout+proc.stderr)
        if name=='ambiguous_steps':ok=proc.returncode!=0 and 'candidate_count=2' in proc.stderr;detail={'returncode':proc.returncode,'stderr':proc.stderr}
        elif proc.returncode:ok=False;detail={'returncode':proc.returncode,'stderr':proc.stderr}
        else:
            observed=json.loads(observed_path.read_text());caps={x['key']:x['state'] for x in observed['integration']['capabilities']};ok=observed['canonical_preserved'] and caps['canonical']=='available'
            if name=='wrong_source':ok=ok and caps['source_association']=='rejected' and caps['nao']!='available'
            elif name=='corrupt_redundant_chain':ok=ok and caps['matrix:NAONLMO:total']=='rejected' and caps['nlmo']=='available' and caps['nho']=='available'
            else:blocked='nho' if kind=='AONHO' else 'nlmo';other='nlmo' if blocked=='nho' else 'nho';ok=ok and caps[blocked]!='available' and caps[other]=='available' and caps['nao']=='available' and caps['aomo']=='available'
            detail=caps
        row=dict(id=name,status='pass' if ok else 'fail',observed=detail);results.append(row);save(dest/'result.json',row);print(name,row['status'],flush=True)
    preserved=all(sha(p)==h for p,h in originals.items());results.append(dict(id='negative_source_bytes_preserved',status='pass' if preserved else 'fail',observed=originals));save(out/'summary.json',results);return 0 if all(x['status']=='pass' for x in results) else 1
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path);p.add_argument('--reuse',type=Path);p.add_argument('--negative',action='store_true');p.add_argument('--cases',nargs='*');p.add_argument('--manifest',type=Path);args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False);cases=json.loads(args.manifest.read_text()) if args.manifest else collect(args.root,args.cases)
    cases=[dict(id=c['case_id'],fchk=c['fchk'],archive=c['archive47'],report=c['report'],matrices={k:c[k.lower()] for k in LFN if k.lower() in c},provenance=c.get('fchk_relationship')) if 'case_id' in c else c for c in cases]
    if args.cases:cases=[c for c in cases if c['id'] in args.cases]
    if args.negative:return negative_cases(cases[0],cases[1],args.output,args.probe)
    save(args.output/'inputs.json',cases);results=[]
    for case in cases:
        try:result=audit(case,args.output/case['id'],args.probe,args.reuse)
        except Exception as e:result=dict(id=case['id'],status='error',error=str(e),traceback=traceback.format_exc());save(args.output/(case['id']+'-error.json'),result)
        results.append(result);print(case['id'],result['status'],flush=True);save(args.output/'summary.json',results)
    return 0 if results and all(x['status']=='pass' for x in results) else 1
if __name__=='__main__':raise SystemExit(main())
