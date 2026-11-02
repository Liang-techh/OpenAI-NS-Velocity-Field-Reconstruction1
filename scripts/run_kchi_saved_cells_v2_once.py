"""Enter the new Kchi saved-cell adapter once, preserving the complete result."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_saved_cells_v2'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_kchi_first_saved_cells_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_kchi_first_saved_cells_fixture_v2.py'
RUNNER = 'scripts/run_kchi_saved_cells_v2_once.py'
PRIOR = 'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz'
FAILURE = 'experiments/root_st073/kchi_saved_cells_v1_failure_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, FAILURE)
PIN = {PRIOR: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       FAILURE: '172cf7a071cc0a27177eaf8cc201daf61715ed3c2264f4760eb3f5a86129fd7a'}


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
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == NAME and brief['phase'] == NAME + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert sha(Path(brief['source_checkpoint']['preimage_path']).read_bytes()) == brief['source_checkpoint']['preimage_sha256']
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert len(manifest['files']) == len(REQUIRED) and {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in manifest['files']:
        p = row['path']
        raw[p] = (ROOT / p).read_bytes()
        hashes[p] = sha(raw[p])
        assert hashes[p] == row['sha256']
        if p in PIN:
            assert hashes[p] == PIN[p]
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_gzip_authenticated_before_decompression': True,
                         'ancestor_called': False, 'source_value_or_integral_computed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'source_value_or_integral_computed': False, 'fixture': None, 'saved_cells': None}
    try:
        leaf = module(LEAF, NAME + '_leaf')
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert previous['status'] == 'returned'
        out['saved_cells'] = leaf.bind(previous['output'])
        assert len(out['saved_cells']['cells']) == 8
        assert out['saved_cells']['source_identity_verified'] is False
        assert out['saved_cells']['tail_is_error_only'] is True
        assert out['saved_cells']['appended_to_cumulative_m15'] is False
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
