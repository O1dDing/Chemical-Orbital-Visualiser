"""Independent real-space bond-axis channel check from raw NBO coefficients.

Uses IOData/GBasis values on cylindrical rings and angular Fourier power.
No production local-angular projector or production orbital coefficients are
used. Full-space ring power is not equated to per-centre S-metric population.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys
from pathlib import Path
from aomo_numerics import np,reference_field_columns,tables,LFN,save
from gbasis.evals.eval import evaluate_basis

def case_in_directory(folder,id):
    p=Path(folder)
    return dict(id=id,fchk=str(p/'canonical.fchk'),archive=str(p/'electronic.47'),report=str(p/'analysis.log'),matrices={k:str(p/f'FILE.{v}') for k,v in LFN.items()})

def ring_test(reference,atoms,coefficients,angle_offset=0.):
    coords=np.asarray(reference['molecule'].atcoords);a,b=coords[list(atoms)];z=(b-a)/np.linalg.norm(b-a);helper=np.eye(3)[np.argmin(abs(z))];x=np.cross(z,helper);x/=np.linalg.norm(x);y=np.cross(z,x)
    nangle=96;angles=np.arange(nangle)*2*np.pi/nangle+angle_offset
    # Two radial windows and multiple axial slices keep a nodal ring from
    # determining the conclusion. Quadrature is cylindrical-volume weighted.
    radii=np.array([.18,.35,.60,.9,1.3,1.8,2.4,3.2]);slices=np.linspace(-.35,1.35,19);points=[]
    for frac in slices:
        centre=a+frac*(b-a)
        for radius in radii:points.extend(centre+radius*(np.cos(angles)[:,None]*x+np.sin(angles)[:,None]*y))
    points=np.asarray(points);ao=evaluate_basis(reference['basis'],points);field=(coefficients.T@ao).reshape((-1,len(slices),len(radii),nangle));fft=np.fft.rfft(field,axis=-1)/nangle
    power=abs(fft)**2;power[...,1:-1]*=2
    weighted=np.sum(power*radii[None,None,:,None],axis=(1,2));weighted/=np.sum(weighted,axis=1)[:,None]
    norm=np.sum(field**2*radii[None,None,:,None],axis=(1,2,3));out=[]
    for k in range(field.shape[0]):
        dominant=int(np.argmax(weighted[k]));sign=(-1)**dominant
        halfturn=np.sqrt(np.sum((np.roll(field[k],nangle//2,axis=-1)-sign*field[k])**2*radii[None,:,None])/norm[k])
        quarter_sign=1 if dominant%4==0 else -1 if dominant%4==2 else None
        quarter=None if quarter_sign is None else float(np.sqrt(np.sum((np.roll(field[k],nangle//4,axis=-1)-quarter_sign*field[k])**2*radii[None,:,None])/norm[k]))
        radial=[]
        for inds in (slice(0,4),slice(4,8)):
            q=np.sum(power[k,:,inds,:]*radii[inds][None,:,None],axis=(0,1));q/=q.sum();radial.append(dict(dominant=int(np.argmax(q)),fraction=float(q[dominant])))
        out.append(dict(dominant_m=dominant,channel=('sigma','pi','delta','phi','gamma')[dominant] if dominant<5 else f'|m|={dominant}',dominant_fraction=float(weighted[k,dominant]),fractions=weighted[k,:9].tolist(),half_turn_relative_error=float(halfturn),quarter_turn_relative_error=quarter,radial_windows=radial))
    return out

def audit(case,output,observed=None):
    reference=reference_field_columns(case);raw=tables(case.get('report'));rows=[];observed_data={};source_binding={}
    if observed:
        data=json.loads(Path(observed).read_text());data=data.get('integration',data)
        observed_data={(x['ref']['spin'],x['ref']['index']):x for x in data['orbitals'] if x['ref']['kind']=='NBO'}
        # Product paths are provenance only; raw values remain the oracle.
        dataset=data['dataset'];pairs=[('report',case['report'],dataset['source']['path'])]
        matrices=case.get('matrices',case)
        aonbo_path=matrices.get('AONBO',matrices.get('aonbo'))
        for matrix in dataset['matrices']:
            if matrix['kind']=='AONBO':pairs.append(('AONBO:'+matrix['spin'],aonbo_path,matrix['source']['path']))
        for key,raw_path,product_path in pairs:
            expected=hashlib.sha256(Path(raw_path).read_bytes()).hexdigest();actual=hashlib.sha256(Path(product_path).read_bytes()).hexdigest()
            source_binding[key]=dict(raw_sha256=expected,product_sha256=actual,equal=expected==actual)
        if not any(key.startswith('AONBO:') for key in source_binding):raise ValueError('Observed product has no AONBO provenance')
    for spin,orbitals in raw['nbo'].items():
        # Explicit atom identities select the target pair; ordinal does not
        # participate in channel assignment or grouping.
        ids=[i for i,o in orbitals.items() if o['identity'] and o['identity'][0] in ('BD','BD*') and len(o['identity'][2])==2 and tuple(atom for symbol,atom in o['identity'][2])==(1,2)]
        if not ids:continue
        coeff=reference['families']['NBO'][spin][:,np.array(ids)-1]
        base=ring_test(reference,(0,1),coeff);rotated=ring_test(reference,(0,1),coeff,.319)
        for i,a,b in zip(ids,base,rotated):
            label=orbitals[i]['identity'];o=observed_data.get((spin,i-1));observed_channel=None
            if o:
                import re
                m=re.search(r'Derived local bond-axis channel: ([^\.]+)',o.get('detail',''));observed_channel=m[1] if m else None
            assigned=observed_channel and observed_channel!='mixed/unresolved'
            agreement=None if not assigned else observed_channel==a['channel']
            rows.append(dict(id=i,spin=spin,producer_identity=label,occupation=orbitals[i]['occupation'],reference=a,azimuth_rotated=b,azimuth_fraction_max_error=float(max(abs(np.array(a['fractions'])-np.array(b['fractions'])))),observed_channel=observed_channel,assigned_label_agrees=agreement,observed_detail=o.get('detail') if o else None))
    missing_labels=sum(r['observed_channel'] is None for r in rows)
    status='pending_observation' if not observed else 'fail' if not rows or missing_labels or any(not b['equal'] for b in source_binding.values()) or any(r['assigned_label_agrees'] is False for r in rows) else 'pass'
    result=dict(case_id=case.get('case_id',case.get('id')),status=status,source_checks=reference['checks'],method='Independent raw-column GBasis cylindrical-ring Fourier and rotation parity; not the production local S-metric projector',observed_path=str(observed) if observed else None,assigned_labels_checked=sum(r['assigned_label_agrees'] is not None for r in rows),missing_product_labels=missing_labels,unresolved_product_labels=sum(r['observed_channel']=='mixed/unresolved' for r in rows),rows=rows)
    result['observed_source_binding']=source_binding
    output.mkdir(parents=True,exist_ok=True);save(output/'result.json',result);print(result['case_id'],[(r['id'],r['reference']['channel'],round(r['reference']['dominant_fraction'],6),r['observed_channel'],r['assigned_label_agrees']) for r in rows],flush=True);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--observed',action='append',default=[]);args=p.parse_args();cases=json.loads(args.manifest.read_text());observed=dict(x.split('=',1) for x in args.observed);args.output.mkdir(parents=True,exist_ok=False);results=[]
    # A supplied manifest controls the exact calculation; evidence is never
    # transferred between similar molecules or different analysis steps.
    for case in cases:results.append(audit(case,args.output/case.get('case_id',case.get('id')),observed.get(case.get('case_id',case.get('id')))))
    save(args.output/'summary.json',results)
    return 0 if all(c['status']=='pass' for c in results) else 2 if any(c['status']=='pending_observation' for c in results) else 1
if __name__=='__main__':raise SystemExit(main())
