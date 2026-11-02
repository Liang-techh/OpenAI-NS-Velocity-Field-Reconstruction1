"""One entry of new target-cell physical-compensation tail propagation."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_target_tail_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_kchi_compensation_target_tail_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_kchi_compensation_target_tail_fixture_v1.py'
RUNNER = 'scripts/run_kchi_target_tail_v1_once.py'
LOW = 'experiments/root_st073/n1_core_low_word_v3_receipt.json.gz'
KCHI = 'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz'
GEOMETRY = 'experiments/root_st073/kchi_geometry_bridge_v2_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, LOW, KCHI, GEOMETRY)
PIN = {LOW: '5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',
       KCHI: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       GEOMETRY: '88e5aeecce42f34e758a429862744c096b47cd00051876495ffe1b8b11d74d84'}


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


def encode_intervals(value):
    if hasattr(value, '_mpi_'):
        return {'lower_exact_mpf_tuple': list(value._mpi_[0]),
                'upper_exact_mpf_tuple': list(value._mpi_[1]), 'directed': True}
    if isinstance(value, dict):
        return {str(k): encode_intervals(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode_intervals(v) for v in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError('cannot serialize result type ' + type(value).__name__)


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
                         'ancestor_called': False, 'fixed_exit_or_m15_recomputed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'fixed_exit_or_m15_recomputed': False, 'fixture': None, 'target_tail': None}
    try:
        leaf = module(LEAF, NAME + '_leaf')
        fixture = module(FIXTURE, NAME + '_fixture')
        c = MPIntervalContext()
        c.prec = 1800
        out['fixture'] = encode_intervals(fixture.run_fixture(leaf.build, c))
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        low = json.loads(gzip.decompress(raw[LOW]))['result']
        kchi = json.loads(gzip.decompress(raw[KCHI]))['result']
        geometry = json.loads(gzip.decompress(raw[GEOMETRY]))['result']
        assert low['status'] == kchi['status'] == geometry['status'] == 'returned'
        assert geometry['geometry']['coordinate_correspondence_verified'] is True
        result = leaf.build(c, low['output'], kchi['output'])
        out['target_tail'] = encode_intervals(result)
        assert result['all_coefficients_symmetric_error_only'] is True
        assert result['source_owner_identity_proved'] is False and result['appended_to_ledger'] is False
        assert result['fixed_x2_exit_computed'] is False and result['complete_W1'] is False
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
