"""One-entry hygienic lowering of the saved Q4 moment substitutions."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_moment_lowering_v4'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_moment_lowering_v4.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_moment_lowering_fixture_v4.py'
RUNNER = 'scripts/run_original_q4_moment_lowering_v4_once.py'
P = 'experiments/root_st073/'
Q4 = P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
JOIN = P + 'original_q4_saved_moment_join_v1_receipt.json.gz'
FAILED = P + 'original_q4_moment_lowering_v3_unresolved_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, Q4, JOIN, FAILED)
FIXED = {Q4: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
         JOIN: '93802554dde9b7b212709b3f85545d53d158a342ea0577d48800101f11eef151',
         FAILED: '4991441ebc3f06b0d7be414bb9c074b3a0073842858cac5542651add209869f8'}

def sha(b): return hashlib.sha256(b).hexdigest()
def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f: f.write(raw)
    return sha(raw)
def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None: raise RuntimeError('new module unavailable')
    obj = importlib.util.module_from_spec(spec); sys.modules[name] = obj; spec.loader.exec_module(obj)
    return obj
def decode(raw): return json.loads(gzip.decompress(raw))

def main():
    if ENTRY.exists() or RESULT.exists(): raise RuntimeError('one-entry component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == COMPONENT and brief['phase'] == COMPONENT + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
    rows = manifest['files']
    assert len(rows) == len(REQUIRED) and {r['path'] for r in rows} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in rows:
        path = row['path']; raw[path] = (ROOT / path).read_bytes(); hashes[path] = sha(raw[path])
        assert hashes[path] == row['sha256']
        if path in FIXED: assert hashes[path] == FIXED[path]
    entry = save(ENTRY, {'component': COMPONENT, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'input_hashes': hashes, 'raw_inputs_authenticated_before_decode': True,
                         'ancestor_called': False, 'quadrature_or_recurrence_recomputed': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'input_hashes': hashes, 'raw_inputs_authenticated_before_decode': True,
           'ancestor_called': False, 'quadrature_or_recurrence_recomputed': False,
           'fixture': None, 'lowerings': []}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.build); out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('hygienic-lowering synthetic fixture failed')
        failed = decode(raw[FAILED])['result']
        assert failed['status'] == 'returned' and len(failed['lowerings']) == 13
        assert all(item['lowering']['status'] == 'unresolved' for item in failed['lowerings'])
        prior = decode(raw[JOIN])['result']
        assert prior['status'] == 'returned' and prior['input_hashes'][Q4] == hashes[Q4]
        joined = prior['archive_result']; assert joined['status'] == 'matched'
        assert joined['scope']['capture_safe_substitution_proved'] is False
        targets = joined['evidence']['moment_rows']; assert len(targets) == 13
        graph = decode(raw[Q4]); assert graph['source_family'] == joined['evidence']['source_family']
        for row in targets:
            inlet = row['q4']['inlet_node']
            lowering = leaf.build(graph, inlet, graph['source_family'])
            if not isinstance(lowering, dict) or lowering.get('status') not in ('lowered', 'unresolved'):
                raise ValueError('invalid lowering result')
            out['lowerings'].append({'moment': row['name'], 'ordinary_Z_order': row['ordinary_Z_order'],
                                     'inlet_node': inlet, 'lowering': lowering})
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    counts = {}
    for item in out['lowerings']:
        key = item['lowering']['status']; counts[key] = counts.get(key, 0) + 1
    print(json.dumps({'status': out['status'], 'lowering_counts': counts,
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
