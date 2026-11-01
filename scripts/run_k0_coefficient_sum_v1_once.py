"""One entry of new directed coefficient grouping over authenticated saved rows."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_coefficient_sum_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
DYAD = 'src/openai_ns_reconstruction/original_bounded_dyadic_sum_v1.py'
DYAD_FIX = 'src/openai_ns_reconstruction/original_bounded_dyadic_sum_fixture_v1.py'
LEAF = 'src/openai_ns_reconstruction/original_k0_grouped_dyadic_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_k0_grouped_dyadic_fixture_v1.py'
RUNNER = 'scripts/run_k0_coefficient_sum_v1_once.py'
PRIOR = 'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz'
LOG = 'experiments/root_st073/k0_log_dyadic_v1_receipt.json.gz'
REQUIRED = (DYAD, DYAD_FIX, LEAF, FIXTURE, RUNNER, PRIOR, LOG)
PIN = {PRIOR: '288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       LOG: 'b60d747b1bbd41442d31d11707e4a54c29a79161c64d5f16b2d59b9c259ba631'}


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
           'source_value_or_integral_computed': False, 'dyadic_fixture': None,
           'group_fixture': None, 'grouped_coefficients': None}
    try:
        dyad = module(DYAD, NAME + '_dyad')
        df = module(DYAD_FIX, NAME + '_dyad_fixture')
        out['dyadic_fixture'] = df.run_fixture(dyad)
        assert out['dyadic_fixture']['passed'] is True and all(out['dyadic_fixture']['checks'].values())
        leaf = module(LEAF, NAME + '_group')
        fixture = module(FIXTURE, NAME + '_group_fixture')
        out['group_fixture'] = fixture.run_fixture(leaf, dyad)
        assert out['group_fixture']['passed'] is True and all(out['group_fixture']['checks'].values())
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert previous['status'] == 'returned'
        symbolic = previous['real_source_expressions']
        assert symbolic['status'] == 'symbolic_unintegrated'
        assert symbolic['source_identity_verified'] is False and symbolic['unknown_terms_are_not_zero'] is True
        assert symbolic['term_counts'] == {'axial': {'Z0': 540, 'Z1': 1620},
                                           'swirl': {'Z0': 396, 'Z1': 1584}}
        grouped = leaf.group(symbolic, dyad, precision=2048)
        out['grouped_coefficients'] = grouped
        assert grouped['member_count'] == grouped['input_term_count'] == 4140
        assert grouped['input_term_counts'] == symbolic['term_counts']
        assert grouped['members_consumed_exactly_once'] is True
        assert grouped['source_identity_verified'] is False and grouped['unknown_terms_are_not_zero'] is True
        assert grouped['appended_to_212_ledger'] is False and grouped['integration_performed'] is False
        assert grouped['sum_intervals_receipt'] == dyad.ROUNDING_RECEIPT
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'dyadic_fixture': (out['dyadic_fixture'] or {}).get('passed'),
                      'group_fixture': (out['group_fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
