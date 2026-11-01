"""Test new dyadic products and conditionally bind a saved log connection."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_log_dyadic_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
DYAD = 'src/openai_ns_reconstruction/original_unbounded_dyadic_product_v1.py'
DYAD_FIX = 'src/openai_ns_reconstruction/original_unbounded_dyadic_product_fixture_v1.py'
LEAF = 'src/openai_ns_reconstruction/original_k0_log_connection_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_k0_log_connection_fixture_v1.py'
RUNNER = 'scripts/run_k0_log_dyadic_v1_once.py'
PRIOR = 'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz'
SWIRL = 'experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz'
CORRECTION = 'experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json'
SOURCE_F = 'src/openai_ns_reconstruction/n1_leading_F_connection_v1.py'
CORE = 'experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json'
REQUIRED = (DYAD, DYAD_FIX, LEAF, FIXTURE, RUNNER, PRIOR, SWIRL, CORRECTION, SOURCE_F, CORE)
PIN = {PRIOR:'288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',
       CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',
       SOURCE_F:'b3a79ea05d5b2786b4530f66b97d31558fe75c8b97f8031b42ec7364ac9198b0',
       CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096'}


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
                         'ancestor_called': False, 'original_source_value_or_integral_computed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'original_source_value_or_integral_computed': False, 'dyadic_fixture': None,
           'connection_fixture': None, 'saved_log_connection': None}
    try:
        dyad = module(DYAD, NAME + '_dyad')
        df = module(DYAD_FIX, NAME + '_dyad_fixture')
        out['dyadic_fixture'] = df.run_fixture(dyad)
        assert out['dyadic_fixture']['passed'] is True and all(out['dyadic_fixture']['checks'].values())
        leaf = module(LEAF, NAME + '_connection')
        fixture = module(FIXTURE, NAME + '_connection_fixture')
        out['connection_fixture'] = fixture.run_fixture(leaf)
        assert out['connection_fixture']['passed'] is True and all(out['connection_fixture']['checks'].values())
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        swirl_receipt = json.loads(gzip.decompress(raw[SWIRL]))
        correction = json.loads(raw[CORRECTION])
        core = json.loads(raw[CORE])
        assert previous['status'] == swirl_receipt['result']['status'] == 'returned'
        assert correction['receipt_sha256'] == hashes[SWIRL] and correction['numerical_result_changed'] is False
        assert correction['input_id'] == 'F0_connection'
        pins = {row['path']: row['sha256'] for row in swirl_receipt['source']['files']}
        assert pins[SOURCE_F] == hashes[SOURCE_F] and pins[CORE] == hashes[CORE]
        symbolic = previous['real_source_expressions']
        swirl = swirl_receipt['result']['output']
        assert symbolic['source_family'] == swirl['source_family']
        assert symbolic['source_logs'] == swirl['source_logs']
        assert all(core[k] == symbolic['source_family'][k] for k in
                   ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256'))
        connection = {'ell_over_lambda': swirl['F0_connection']['ell_over_lambda'],
                      'source_family': swirl['source_family'], 'source_logs': swirl['source_logs'],
                      'authenticated_raw_receipt_hash': hashes[SWIRL],
                      'provenance_correction': {'declared_by_runner': True, 'receipt_sha256': hashes[SWIRL],
                                                'correction_sha256': hashes[CORRECTION], 'record': correction},
                      'saved_core_Z_domain': core['whole_core']['Z'],
                      'same_live_F0_owner_proved': False, 'source_identity_proved': False}
        out['saved_log_connection'] = leaf.bind(symbolic, connection)
        expected = sum(1 for cell in symbolic['cells'] for term in cell['terms']['swirl']['Z1']
                       if term['formal_source_powers'][5] == 1)
        assert expected == 396 and out['saved_log_connection']['term_counts']['swirl']['Z1'] == expected
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'dyadic_fixture': (out['dyadic_fixture'] or {}).get('passed'),
                      'connection_fixture': (out['connection_fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
