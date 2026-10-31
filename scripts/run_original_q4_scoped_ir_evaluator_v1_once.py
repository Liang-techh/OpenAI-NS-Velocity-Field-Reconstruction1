"""One frozen evaluator entry, with synthetic tests and no real-source oracle."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_scoped_ir_evaluator_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_scoped_ir_evaluator_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_scoped_ir_evaluator_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_scoped_ir_evaluator_v1_once.py'
PRIOR = 'experiments/root_st073/original_q4_moment_lowering_v4_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR)
PRIOR_SHA = '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f:
        f.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new component unavailable')
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == NAME and brief['phase'] == NAME + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in manifest['files']:
        p = row['path']
        raw[p] = (ROOT / p).read_bytes()
        hashes[p] = sha(raw[p])
        assert hashes[p] == row['sha256']
    assert hashes[PRIOR] == PRIOR_SHA
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_authenticated_before_decode': True,
                         'ancestor_called': False, 'real_source_quadrature_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'real_source_quadrature_called': False, 'fixture': None, 'actual_IR_gate': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        # Convert Fractions in synthetic case records to strings only when serializing.
        got = fixture.run_fixture(leaf.evaluate)
        out['fixture'] = json.loads(json.dumps(got, default=str))
        assert got['status'] == 'passed' and got['failed'] == 0 and got['passed'] > 0
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert prior['status'] == 'returned' and len(prior['lowerings']) == 13
        first = prior['lowerings'][0]
        assert first['lowering']['status'] == 'lowered'
        gate = leaf.evaluate(first['lowering'], oracle=None, free_variables={}, backend=None)
        assert gate['status'] == 'rejected' and gate['value'] is None
        out['actual_IR_gate'] = {'saved_moment': first['moment'], 'ordinary_Z_order': first['ordinary_Z_order'],
                                 'result': gate, 'real_source_oracle_supplied': False,
                                 'actual_source_value_computed': False, 'saved_lowering_count': 13}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture': (out['fixture'] or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
