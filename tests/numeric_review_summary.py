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

"""Summarize a terminal full numerical review without relabelling formal cases."""
from collections import Counter
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import statistics
import time


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def distribution(values):
    values = sorted(values)
    if not values:
        return None
    def percentile(fraction):
        position = fraction * (len(values) - 1)
        lo, hi = math.floor(position), math.ceil(position)
        return values[lo] + (values[hi] - values[lo]) * (position - lo)
    return {'n': len(values), 'sum_seconds': sum(values), 'minimum_seconds': values[0],
            'median_seconds': statistics.median(values), 'p95_seconds': percentile(.95),
            'maximum_seconds': values[-1]}


parser = argparse.ArgumentParser()
parser.add_argument('--round', type=Path, required=True)
parser.add_argument('--review-id', default='review-v1')
args = parser.parse_args()
root = args.round.resolve()
review_root = root / args.review_id
manifest = read(root / 'manifest.json')
collection = read(root / 'collection-complete.json')
reviews = read(review_root / 'review-complete.json')
expected = {case['case_id'] for case in manifest['cases']}
if (not collection['all_terminal'] or not reviews['all_terminal']
        or {case['case_id'] for case in collection['cases']} != expected
        or {case['case_id'] for case in reviews['cases']} != expected
        or collection['round_identity'] != manifest['round_identity']
        or reviews['round_identity'] != manifest['round_identity']):
    raise RuntimeError('Matching complete collection/review barriers are required')
output = review_root / 'unified-review'
if output.exists():
    raise FileExistsError('An existing numerical summary is never overwritten')
output.mkdir()
checks = {}
failures = []
texture_failures = []
rows = []
review_times = []
sampled_count = complete_count = 0
sources = {case['case_id']: case for case in manifest['cases']}
collected = {case['case_id']: case for case in collection['cases']}
for terminal in reviews['cases']:
    case_id = terminal['case_id']
    row = {'case_id': case_id, 'input': sources[case_id]['relative_path'],
           'collection_status': collected[case_id]['collection_status'],
           'numerical_status': terminal['status'], 'formal_case_pass': False,
           'collection_seconds': collected[case_id].get('wall_seconds_including_ordinary_and_packaging'),
           'numerical_review_seconds': None, 'failed_checks': []}
    if terminal['status'] == 'review_error':
        failures.append({'case_id': case_id, 'check': 'REVIEW-EXECUTION',
                         'status': 'error', 'detail': terminal['error'], 'evidence': terminal['attempt']})
        row['failed_checks'] = ['REVIEW-EXECUTION']
    else:
        report_path = review_root / terminal['review']
        for consumed in (report_path, report_path.parent / 'sampled-textures.json',
                         report_path.parent / 'complete-textures.json'):
            key = str(consumed.relative_to(review_root / terminal['attempt']))
            if sha(consumed) != terminal['artifacts'][key]['sha256']:
                raise RuntimeError('Review artifact changed before aggregation: ' + str(consumed))
        result = read(report_path)
        if result['round_identity'] != manifest['round_identity'] or result['case_id'] != case_id:
            raise RuntimeError('Result identity mismatch')
        for check in result['checks']:
            checks.setdefault(check['check'], Counter())[check['status']] += 1
            if check['status'] != 'pass':
                failures.append({'case_id': case_id, **check, 'evidence': terminal['review']})
                row['failed_checks'].append(check['check'])
        row['numerical_review_seconds'] = result['wall_seconds']
        review_times.append(result['wall_seconds'])
        for name, kind in (('sampled-textures.json', 'sampled'), ('complete-textures.json', 'complete')):
            for item in read(report_path.parent / name):
                if name == 'sampled-textures.json' and item['full_grid']:
                    continue
                if kind == 'sampled':
                    sampled_count += 1
                else:
                    complete_count += 1
                if not item['pass']:
                    texture_failures.append({'case_id': case_id, 'kind': kind, **item})
    rows.append(row)
counts = dict(Counter(row['numerical_status'] for row in rows))
result = {'schema': 1, 'round_identity': manifest['round_identity'],
          'collection_sha256': sha(root / 'collection-complete.json'),
          'review_barrier_sha256': sha(review_root / 'review-complete.json'),
          'summary_script_sha256': sha(Path(__file__)), 'created_epoch': time.time(),
          'case_count': len(expected), 'collection_counts': collection['collection_counts'],
          'numeric_subset_counts': counts, 'check_counts': checks,
          'sampled_mos': sampled_count, 'complete_frontier_textures': complete_count,
          'failed_check_count': len(failures), 'failed_texture_count': len(texture_failures),
          'failures': failures, 'texture_failures': texture_failures, 'cases': rows,
          'timing': {'collection_wall_seconds': collection['wall_seconds_this_session'],
                     'independent_review_wall_seconds': reviews['wall_seconds_this_session'],
                     'collection_per_case': distribution([row['collection_seconds'] for row in rows
                                                          if row['collection_seconds'] is not None]),
                     'independent_review_per_case': distribution(review_times)},
          'formal_case_passes': 0,
          'remaining': ['paused physical reference calculations and adjudication',
                        'full independent stored-density checks',
                        'symmetry/topology/scoped chemistry/UI/export issue closure',
                        'all required invariants and ordinary-build interaction consistency']}
(output / 'summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
with (output / 'case-status.csv').open('w', newline='', encoding='utf-8-sig') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    for row in rows:
        writer.writerow({**row, 'failed_checks': ';'.join(row['failed_checks'])})
text = [f'# {manifest["name"]} 独立数值复核', '',
        f'全部 {len(expected)} 个案例完成取证并进入数值复核终态。程序提交 `{manifest["git_commit"]}`。', '',
        f'数值子集状态：{counts}。正式完整通过数仍为 **0**。', '',
        f'共检查 {sampled_count} 个 MO 的实际纹理采样及 {complete_count} 份完整前线纹理；',
        f'发现 {len(failures)} 项未通过的检查，涉及 {len({row["case_id"] for row in failures})} 个案例。', '',
        f'取证墙钟 {collection["wall_seconds_this_session"]:.3f} 秒；独立复核墙钟 {reviews["wall_seconds_this_session"]:.3f} 秒。', '',
        '这些耗时没有包含仍暂停的 Gaussian 参考计算、未完成的图像/化学裁决和后续修复。', '',
        '逐案例状态见 `case-status.csv`；逐检查数值、失败纹理及证据位置见 `summary.json`。', '',
        '| 检查 | 通过 | 失败 | 待补证据 |', '|---|---:|---:|---:|']
for name, count in sorted(checks.items()):
    text.append(f'| {name} | {count.get("pass",0)} | {count.get("fail",0)} | {count.get("insufficient",0)} |')
text += ['', '数值通过仅覆盖本轮明列的项目。未完成项没有被改成“不适用”，失败案例也没有移出正式集合。', '']
(output / '数值复核结果.md').write_text('\n'.join(text), encoding='utf-8')
print(json.dumps({key: result[key] for key in ('case_count', 'numeric_subset_counts', 'sampled_mos',
                                               'failed_check_count', 'failed_texture_count', 'timing')}, ensure_ascii=False))
