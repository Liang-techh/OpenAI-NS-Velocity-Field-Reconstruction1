"""One frozen entry for the old/current endpoint log-radius source packet."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_logradius_source_binding_v2'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
P = 'experiments/root_st073/'
LEAF = 'src/openai_ns_reconstruction/original_logradius_source_binding_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_logradius_source_binding_fixture_v2.py'
RUNNER = 'scripts/run_original_logradius_source_binding_v2_once.py'
PRIOR = P + 'original_f_main_gap_signed_seam_v1_receipt.json.gz'
OWN = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.json.gz'
CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields_check.json'
FAILURE = P + 'original_logradius_source_binding_v1_failure_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, OWN, CHECK, FAILURE)
PIN = {PRIOR: 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889',
       OWN: '2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b',
       CHECK: '0bc8229a1716871d896161caac907c572f7cd0a58540900aaa00869d3a920fb1',
       FAILURE: '1b52dab01a5a43f7f0a76f5c48c49734cdf9736b46788973aca6ce3b953d05c6'}


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
        if p in PIN:
            assert hashes[p] == PIN[p]
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_authenticated_before_decode': True,
                         'ancestor_called': False, 'protected_graph_refs_indexed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'protected_graph_refs_indexed': False, 'fixture': None, 'archive_result': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf.build)
        assert out['fixture']['passed'] is True
        assert all(v is True for v in out['fixture']['checks'].values())
        saved = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert saved['status'] == 'returned'
        previous = saved['archive_result']
        assert previous['status'] == 'constructed'
        own = json.loads(gzip.decompress(raw[OWN]))
        check = json.loads(raw[CHECK])
        assert check['all_passed'] is True
        assert check['input_hashes'][Path(OWN).name] == hashes[OWN]
        assert own['source_family'] == check['source_family'] == previous['source_family']['own']['fields']
        out['archive_result'] = leaf.build(own, previous)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture_passed': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
