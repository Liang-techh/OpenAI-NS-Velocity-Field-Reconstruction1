"""Freeze-first, single-entry Q4 source-port instance construction."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_sc_flat_instance_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_sc_flat_instance_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_sc_flat_instance_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_sc_flat_instance_v1_once.py'
INSTANCE = 'experiments/root_st073/original_q4_sc_instance_v1_receipt.json.gz'
FLAT = 'experiments/root_st073/original_q4_slow_B_flat_v1_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, INSTANCE, FLAT)
FIXED = {
    INSTANCE: '2438f3016b3442d369a2f199424910d820a6f77eb3bd4f83f6d83d7907319a37',
    FLAT: 'eaee3662dc501d65b282517e59c1ddb55cc91bcbc7bbda9703cc7faa2d27e92f',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new module unavailable')
    out = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('one-entry component already entered')
    brief_raw, manifest_raw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(brief_raw), json.loads(manifest_raw)
    assert brief['component'] == COMPONENT and brief['phase'] == COMPONENT + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(manifest_raw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
    rows = manifest['files']
    assert len(rows) == len(REQUIRED) and {r['path'] for r in rows} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in rows:
        path = row['path']
        raw[path] = (ROOT / path).read_bytes()
        hashes[path] = sha(raw[path])
        assert hashes[path] == row['sha256']
        if path in FIXED:
            assert hashes[path] == FIXED[path]
    entry = save(ENTRY, {'component': COMPONENT, 'source_sha256': sha(manifest_raw),
                         'brief_sha256': sha(brief_raw), 'input_hashes': hashes,
                         'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry,
           'source_sha256': sha(manifest_raw), 'input_hashes': hashes,
           'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False,
           'fixture': None, 'port_fixture': None, 'instance_result': None}
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        for key, fun, argument in (('fixture', fixture.run_fixture, leaf.join),
                                   ('port_fixture', fixture.run_port_fixture, leaf.evaluate_port)):
            got = fun(argument)
            out[key] = got
            if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
                raise ValueError(key + ' failed')
        instance = json.loads(gzip.decompress(raw[INSTANCE]))
        flat = json.loads(gzip.decompress(raw[FLAT]))
        out['instance_result'] = leaf.join(instance, flat)
        if out['instance_result'].get('status') not in ('verified', 'mismatch', 'unresolved'):
            raise ValueError('invalid instance result')
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'instance_status': (out['instance_result'] or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
