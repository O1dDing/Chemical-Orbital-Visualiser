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

"""Review preserved native export frames only after the whole collection ends.

Checks identity, members, raw values, selection, draw-call semantics and selected
PNG strokes. This is not a substitute for complete screenshot/chemistry review.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sys
import time
import xml.etree.ElementTree as ET

from validation_process import atomic_json


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))


def overlapping_rectangles(rectangles):
    result=[]
    for i,(identity,left,right,top,bottom) in enumerate(rectangles):
        for other,other_left,other_right,other_top,other_bottom in rectangles[i+1:]:
            if min(right,other_right)>max(left,other_left) and min(bottom,other_bottom)>max(top,other_top):
                result.append([identity,other])
    return result


def review_details_frame(data, selected_mo, visible_target, image_size):
    """Check actual control/capture evidence without claiming glyph correctness."""
    failures=[]
    state=data['state']
    if not state['scene_matches_applied'] or any(state[k]!=selected_mo for k in
            ('rendered_mo','applied_mo','drawn_ui_mo','requested_mo')):
        failures.append('selected-member-mismatch')
    details=[t for t in data['targets'] if t['id'].startswith('diagram.details.')]
    if visible_target is None:
        if details: failures.append('window-remains-open')
    else:
        targets=[t for t in details if t['id']==visible_target]
        if len(targets)!=1 or not targets[0]['visible']:
            failures.append('required-control-not-visible')
        elif not all(math.isfinite(v) for v in targets[0]['rect']):
            failures.append('nonfinite-control-bounds')
        else:
            x0,y0,x1,y1=targets[0]['rect'];width,height=image_size
            if not (0<=x0<x1<=width and 0<=y0<y1<=height):
                failures.append('required-control-outside-frame')
            clip=targets[0].get('clip_rect')
            if clip is None:
                failures.append('actual-clip-evidence-missing')
            elif not all(math.isfinite(v) for v in clip) or not (clip[0]<=x0<x1<=clip[2] and clip[1]<=y0<y1<=clip[3]):
                failures.append('required-control-partially-clipped')
    return failures


def review_case(root,record,Image):
    terminal=load(root/'cases'/record['case_id']/'terminal.json')
    evidence=root/terminal['evidence_directory']
    manifest_file=root/terminal['evidence_manifest']
    assert sha(manifest_file)==terminal['evidence_manifest_sha256']
    artifacts={Path(k).as_posix():v for k,v in load(manifest_file).items()}
    consumed={}
    def consume(relative):
        relative=Path(relative).as_posix()
        path=evidence/relative
        identity=sha(path)
        assert identity==artifacts[relative]['sha256'],relative+' changed'
        consumed[relative]=identity
        return path
    native=evidence/'native'
    expected=load(consume('expected-evidence.json'))
    production=load(consume('production.json'))
    events=[json.loads(line) for line in consume('native/events.jsonl').read_text(encoding='utf-8').splitlines()]
    traces={x['frame']:x['data'] for x in events if x['kind']=='export.frame-trace'}
    exports=[x for x in events if x['kind']=='export.actual']
    failures=[]; bundles=[]
    def check(ok,code,detail):
        if not ok: failures.append({'check':code,'detail':detail})
    actions=[json.loads(line) for line in consume('native/actions.jsonl').read_text(encoding='utf-8').splitlines()]
    failed_actions=[row for row in actions if row['status']!='executed']
    check(not failed_actions,'VIEW-ACTION-FAILURES',failed_actions)
    names=expected.get('export_names',[])
    check(bool(names),'VIEW-COVERAGE','No frozen export-state coverage specified')
    check(len(names)==len(set(names))==len(exports),'VIEW-EXPORT-COUNT',
          {'expected':len(names),'measured':len(exports)})
    details_frames=[]
    if expected.get('orbital_details_captures'):
        controls={'orbital-details-top':'diagram.details.close',
                  'orbital-details-bottom':'diagram.details.scope',
                  'orbital-details-closed':None,
                  'orbital-details-reopened':'diagram.details.close'}
        check(set(expected['orbital_details_captures'])==set(controls),
              'VIEW-DETAILS-COVERAGE',expected['orbital_details_captures'])
        for capture,target in controls.items():
            data=load(consume('native/'+capture+'.ui.json'))
            with Image.open(consume('native/'+capture+'.png')) as image:
                image.load();size=image.size
            issues=review_details_frame(data,expected['orbital_details_member'],target,size)
            if capture=='orbital-details-bottom':
                issues += review_details_frame(data,expected['orbital_details_member'],'diagram.details.close.bottom',size)
            check(not issues,'VIEW-DETAILS-CONTROL',{'capture':capture,'failures':issues})
            details_frames.append({'capture':capture,'selected_mo':expected['orbital_details_member'],
                                   'required_visible_control':target,'image_size':list(size),'pass':not issues})
    common_details=[]
    for request in expected.get('common_root_details',[]):
        capture=request['capture'];ui=load(consume('native/'+capture+'.ui.json'))
        with Image.open(consume('native/'+capture+'.png')) as image:size=image.size
        state=ui['state'];targets=[t for t in ui['targets'] if t['id']==request['target']]
        visible=len(targets)==1 and targets[0]['visible']
        if visible:
            x0,y0,x1,y1=targets[0]['rect'];clip=targets[0]['clip_rect']
            visible=(0<=x0<x1<=size[0] and 0<=y0<y1<=size[1] and
                     clip[0]<=x0<x1<=clip[2] and clip[1]<=y0<y1<=clip[3])
        passed=visible and state['scene_matches_applied'] and all(state[k]==request['mo'] for k in
            ('rendered_mo','applied_mo','drawn_ui_mo','requested_mo'))
        check(passed,'VIEW-COMMON-DETAILS',request)
        common_details.append({**request,'pass':passed})
    linear_clicks=[]
    for mo in expected.get('primary_members',[]):
        capture=f'linear-member-{mo+1:04}'
        state=load(consume('native/'+capture+'.ui.json'))['state']
        consume('native/'+capture+'.png')
        passed=state['axis_mode']==0 and state['scene_matches_applied'] and all(
            state[k]==mo for k in ('rendered_mo','applied_mo','drawn_ui_mo','requested_mo'))
        check(passed,'VIEW-LINEAR-MEMBER-CLICK',{'capture':capture,'expected_mo':mo,'actual':state})
        linear_clicks.append({'mo':mo,'pass':passed,'capture':capture})
    for name in names:
        begin=len(failures)
        matches=[x for x in exports if Path(x['data']['base']).name==name]
        check(len(matches)==1,'VIEW-EXPORT-EVENT',name)
        if len(matches)!=1: continue
        event=matches[0]; capture_trace=traces.get(event['frame'],[])
        snapshot_records=[x['data'] for x in capture_trace if x['kind']=='diagram.snapshot']
        check(len(snapshot_records)==1,'VIEW-FRAME-SNAPSHOT',name)
        if len(snapshot_records)!=1: continue
        drawn=snapshot_records[0]
        paths={ext:consume('native/'+name+'.mo.'+ext) for ext in ('json','csv','svg','png')}
        data=load(paths['json']); context=data.get('view_snapshot') or {}
        # Scoped evidence can exceed the CSV module's default 128 KiB field
        # limit. A decoded field cannot be larger than its verified input file.
        csv.field_size_limit(max(csv.field_size_limit(),paths['csv'].stat().st_size))
        with paths['csv'].open(encoding='utf-8',newline='') as csv_stream:
            csv_rows=list(csv.DictReader(csv_stream))
        svg=ET.parse(paths['svg']).getroot()
        with Image.open(paths['png']) as image:
            png_id=image.info.get('cov.view.snapshot')
            image.load(); pixels=image.convert('RGB')
        check(int(svg.attrib['width'])==pixels.width and int(svg.attrib['height'])==pixels.height,
              'VIEW-EXPORT-DIMENSIONS',name)
        check(event['data']['success'] and context.get('origin')=='interactive-canvas' and
              drawn['id']==event['data'].get('snapshot_id')==context.get('id')==
              svg.attrib.get('data-view-snapshot')==png_id,'VIEW-IDENTITY',name)
        selected=context.get('inspected_orbital_index')
        check(selected==drawn['inspected_orbital_index']==event['data'].get('selected_index'),
              'VIEW-INSPECTION',name)
        metadata={m['index']:m for m in data['orbitals']}
        csv_metadata={int(m['index']):m for m in csv_rows}
        check(len(metadata)==len(data['orbitals']) and set(metadata)==set(csv_metadata),
              'VIEW-RAW-METADATA-COVERAGE',name)
        for i,mo in metadata.items():
            source=production['orbitals'][i]; item=csv_metadata.get(i,{})
            check(abs(mo['energy_hartree']-source['energy_hartree'])<=1e-12 and
                  abs(mo['occupation']-source['occupation'])<=1e-7,
                  'VIEW-RAW-VALUES',{'bundle':name,'mo':i})
            check(mo['selected']==(i==selected) and item.get('selected')==str(int(i==selected)) and
                  item.get('view_snapshot_id')==context.get('id') and
                  item.get('energy_unit')==data['energy_unit'] and
                  math.isclose(float(item.get('energy_hartree','nan')),mo['energy_hartree'],abs_tol=1e-12),
                  'VIEW-CSV-JSON-SELECTION',{'bundle':name,'mo':i})
            check((mo['spin'].lower().startswith('b'))==(source['spin']==1),
                  'VIEW-RAW-SPIN',{'bundle':name,'mo':i})
        rows=data['diagram_rows']
        # CODATA 2022; independent conversion values, same tolerance as the
        # dedicated six-unit control. Record fields keep the physical type.
        joule=4.3597447222060e-18*6.02214076e23
        factors={'Ha':1.0,'eV':27.211386245981,'J/mol':joule,'kJ/mol':joule/1000,
                 'cal/mol':joule/4.184,'kcal/mol':joule/4184}
        for array_name,kind in [('pi_interactions','pi-partner'),('crystal_field_gaps','crystal-field')]:
            records=data.get(array_name,[])
            svg_records=[node for node in svg.iter() if node.get('data-gap-kind')==kind]
            check(len(records)==len(svg_records),'VIEW-TYPED-GAP-COUNT',{'bundle':name,'kind':kind})
            for gap,node in zip(records,svg_records):
                lower=[production['orbitals'][i]['energy_hartree'] for i in gap['lower_orbitals']]
                upper=[production['orbitals'][i]['energy_hartree'] for i in gap['upper_orbitals']]
                valid=bool(lower and upper) and gap['gap_kind']==kind
                if valid:
                    lower_mean=sum(lower)/len(lower);upper_mean=sum(upper)/len(upper)
                    split=upper_mean-lower_mean
                    valid=(math.isclose(lower_mean,gap['lower_energy_hartree'],abs_tol=1e-12) and
                           math.isclose(upper_mean,gap['upper_energy_hartree'],abs_tol=1e-12) and
                           math.isclose(split,gap['splitting_hartree'],abs_tol=1e-12) and
                           math.isclose(split*factors[gap['display_unit']],gap['splitting_display'],rel_tol=1e-10,abs_tol=1e-12))
                check(valid,'VIEW-TYPED-GAP-VALUES',{'bundle':name,'kind':kind,'lower':gap['lower_orbitals'],'upper':gap['upper_orbitals']})
                check(json.loads(node.attrib['data-evidence'])==gap,'VIEW-TYPED-GAP-SVG',{'bundle':name,'kind':kind})
            column='pi_partner_evidence_json' if kind=='pi-partner' else 'crystal_field_gaps_json'
            check(all(json.loads(row[column])==records for row in csv_rows),'VIEW-TYPED-GAP-CSV',{'bundle':name,'kind':kind})
        svg_rows=[node for node in svg if node.attrib.get('class')=='mo-level']
        check(len(rows)==data['diagram_row_count']==drawn['row_count']==len(svg_rows),
              'VIEW-ROW-COUNT',name)
        native_members=[x['data'] for x in capture_trace if x['kind']=='draw.level']
        native_rects=[]
        for member in native_members:
            if not all(k in member for k in ('x','y','hit_half_width','hit_half_height')):
                check(False,'VIEW-HIT-GEOMETRY-MISSING',{'bundle':name,'member':member.get('used_internal_mo')})
                continue
            x,y=member['x'],member['y']; w,h=member['hit_half_width'],member['hit_half_height']
            native_rects.append((member['used_internal_mo'],x-w,x+w,y-h,y+h))
        collisions=overlapping_rectangles(native_rects)
        check(not collisions,'VIEW-HIT-REGIONS-OVERLAP',{'bundle':name,'pairs':collisions})
        expected_native=[]
        selected_png=[]
        svg_rects=[]
        for row_index,row in enumerate(rows):
            members=row.get('members',[])
            check(bool(members),'VIEW-MEMBERS-MISSING',{'bundle':name,'row':row_index})
            marks=[] if row_index>=len(svg_rows) else list(svg_rows[row_index].iter())
            marks=[m for m in marks if m.attrib.get('class')=='mo-member']
            check(len(marks)==len(members),'VIEW-SVG-MEMBER-COUNT',{'bundle':name,'row':row_index})
            energies=[]
            for member_index,member in enumerate(members):
                indices=[member['orbital_index']]+([member['spin_counterpart']]
                    if member['spin_counterpart'] is not None else [])
                expected_selected=selected in indices
                expected_used=selected if expected_selected else indices[0]
                check(member['selected']==expected_selected and member['inspected_orbital_index']==expected_used,
                      'VIEW-MEMBER-SELECTION',{'bundle':name,'row':row_index,'member':member_index})
                for i in indices:
                    check(i in metadata and metadata[i]['diagram_row_index']==row_index and
                          metadata[i]['included_in_diagram'],'VIEW-MEMBER-RAW-LINK',{'bundle':name,'mo':i})
                    energies.append(production['orbitals'][i]['energy_hartree'])
                expected_native.append((expected_used,row['representative_metadata_index'],
                    row['layout_energy_hartree'],member['alpha_arrows'],member['beta_arrows']))
                if member_index>=len(marks): continue
                attrs=marks[member_index].attrib
                half=float(attrs['stroke-width'])/2
                svg_rects.append((member['orbital_index'],float(attrs['x1']),float(attrs['x2']),
                                  float(attrs['y1'])-half,float(attrs['y1'])+half))
                check(0<=float(attrs['x1'])<float(attrs['x2'])<pixels.width,
                      'VIEW-EXPORT-MEMBER-BOUNDS',{'bundle':name,'member':member['orbital_index']})
                check(int(attrs['data-orbital-index'])==indices[0] and
                      int(attrs['data-inspected-orbital-index'])==expected_used and
                      (attrs['data-selected']=='true')==expected_selected,
                      'VIEW-SVG-SELECTION',{'bundle':name,'row':row_index,'member':member_index})
                if expected_selected:
                    # Inspect actual decoded PNG pixels at the selected stroke,
                    # away from the two electron-arrow shafts.
                    x0,x1,y=float(attrs['x1']),float(attrs['x2']),float(attrs['y1'])
                    probes=[(round(x),round(y)+dy) for x in (x0+2,x1-2) for dy in (-1,0,1)]
                    visible=any(0<=x<pixels.width and 0<=yy<pixels.height and
                                pixels.getpixel((x,yy))==(40,109,224) for x,yy in probes)
                    check(visible,'VIEW-PNG-SELECTED-STROKE',{'bundle':name,'row':row_index,'probes':probes})
                    selected_png.append(visible)
            check(bool(energies) and math.isclose(min(energies),row.get('all_members_energy_min_hartree',math.nan),abs_tol=1e-12)
                  and math.isclose(max(energies),row.get('all_members_energy_max_hartree',math.nan),abs_tol=1e-12),
                  'VIEW-ROW-ENERGY-RANGE',{'bundle':name,'row':row_index})
        collisions=overlapping_rectangles(svg_rects)
        check(not collisions,'VIEW-EXPORT-STROKES-OVERLAP',{'bundle':name,'pairs':collisions})
        check(len(expected_native)==len(native_members),'VIEW-FRAME-MEMBER-COUNT',name)
        for wanted,actual in zip(expected_native,native_members):
            check(wanted[0]==actual['used_internal_mo'] and wanted[1]==actual['level_metadata_internal_mo'] and
                  math.isclose(wanted[2],actual['layout_energy_hartree'],abs_tol=1e-12) and
                  wanted[3]==actual['alpha_arrows'] and wanted[4]==actual['beta_arrows'],
                  'VIEW-DRAW-CALL-CONSISTENCY',{'bundle':name,'wanted':wanted,'actual':actual})
        bundles.append({'name':name,'status':'pass' if len(failures)==begin else 'fail',
            'snapshot_id':context.get('id'),'export_frame':event['frame'],'selected_mo':selected,
            'rows':len(rows),'members':len(expected_native),'selected_png_checks':len(selected_png)})
    return {'case_id':record['case_id'],'status':'view_subset_pass' if not failures else 'view_subset_fail',
        'linear_member_clicks':linear_clicks,'details_control_frames':details_frames,'common_root_details':common_details,
        'bundles':bundles,'failures':failures,'consumed_artifacts':consumed,'formal_case_pass':False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--round',type=Path,required=True)
    parser.add_argument('--review-id',default='view-review-v1')
    args=parser.parse_args(); root=args.round.resolve()
    if Path(args.review_id).name!=args.review_id: raise ValueError('review-id must be a directory name')
    manifest=load(root/'manifest.json'); barrier=load(root/'collection-complete.json')
    assert barrier['all_terminal'] and barrier['terminal_cases']==manifest['case_count']
    assert barrier['round_identity']==manifest['round_identity']
    out=root/args.review_id
    if out.exists(): raise FileExistsError('Existing review evidence is never overwritten')
    for relative,identity in manifest['image_dependency_artifacts'].items():
        assert sha(root/'image-deps'/relative)==identity
    sys.path.insert(0,str(root/'image-deps'))
    from PIL import Image
    out.mkdir(); start=time.perf_counter(); results=[]
    for record in manifest['cases']:
        try: result=review_case(root,record,Image)
        except Exception as error:
            result={'case_id':record['case_id'],'status':'view_subset_error',
                    'error':repr(error),'formal_case_pass':False}
        atomic_json(out/(record['case_id']+'.json'),result); results.append(result)
    result={'round_identity':manifest['round_identity'],'case_count':len(results),
        'all_reviewed':len(results)==manifest['case_count'],
        'status_counts':{s:sum(x['status']==s for x in results) for s in sorted({x['status'] for x in results})},
        'bundle_count':sum(len(x.get('bundles',[])) for x in results),
        'details_control_frame_count':sum(len(x.get('details_control_frames',[])) for x in results),
        'failure_count':sum(len(x.get('failures',[])) for x in results),
        'wall_seconds':time.perf_counter()-start,'reviewer_sha256':sha(Path(__file__)),
        'collection_barrier_sha256':sha(root/'collection-complete.json'),'formal_case_passes':0,
        'limitations':['Draw-call semantics do not prove all visible text/glyph positions or clipping.',
                       'PNG stroke probes do not certify all output pixels.',
                       'Chemical/symmetry interpretation and independently managed physical references remain separate.'],
        'case_reports':{r['case_id']+'.json':sha(out/(r['case_id']+'.json')) for r in results}}
    atomic_json(out/'summary.json',result)
    print(json.dumps(result),flush=True)
    return 0 if result['status_counts']=={'view_subset_pass':manifest['case_count']} else 2


if __name__=='__main__': raise SystemExit(main())
