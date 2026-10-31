"""Enter the new guarded integral fixture once; never invoke an ancestor."""
import gzip
import hashlib
import dataclasses
import importlib.util
import json
import sys
import traceback
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_integral_coverage_guard_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_integral_coverage_guard_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_integral_coverage_guard_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_integral_coverage_guard_v1_once.py'
PRIOR = 'experiments/root_st073/original_q4_interval_ir_runtime_v1_failure_receipt.json.gz'
RUNTIME = 'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_v2.py'
RUNTIME_RECEIPT = 'experiments/root_st073/original_q4_interval_ir_runtime_v2_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, RUNTIME, RUNTIME_RECEIPT)
PIN = {PRIOR: '66e9461771d3cf08334417846a8acca662c48b21be4c197adc6bfba21d554093',
       RUNTIME_RECEIPT: '1c55c1cb4699eb8fea3e874209777455fb52467ec02d40ca8633b48dcf581c56'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    def encode_unknown(value):
        if isinstance(value, Fraction):
            return {'type': 'Fraction', 'numerator': value.numerator, 'denominator': value.denominator}
        if dataclasses.is_dataclass(value) and not isinstance(value, type):
            return {'type': type(value).__name__, **dataclasses.asdict(value)}
        raise TypeError('unsupported result value: ' + type(value).__name__)
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False,
                      default=encode_unknown) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new component unavailable')
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def rational(obj):
    assert obj['type'] == 'Fraction'
    n, d = obj['numerator'], obj['denominator']
    assert type(n) is int and type(d) is int and d > 0
    return Fraction(n, d)


def endpoints(obj):
    assert obj['type'] == 'FreshInterval' and obj['kind'] == 'enclosure'
    lo, hi = rational(obj['lo']), rational(obj['hi'])
    assert lo <= hi
    return lo, hi


def saved_fixture_geometry(raw):
    saved = json.loads(gzip.decompress(raw))['result']
    assert saved['status'] == 'failed' and saved['fixture']['status'] == 'passed'
    records = saved['fixture']['successful_result']['integral_receipts']
    assert len(records) == 9
    counts = []
    for rec in records:
        lo, lo2 = endpoints(rec['lower'])
        hi, hi2 = endpoints(rec['upper'])
        assert lo == lo2 and hi == hi2 and lo < hi
        cells = rec['coverage_cells']
        assert rec['coverage_count'] == len(cells) and cells
        cursor = lo
        for cell in cells:
            left, right = endpoints(cell)
            assert left == cursor and left < right
            cursor = right
        assert cursor == hi
        counts.append(len(cells))
    return {'saved_records': len(records), 'cell_counts': counts,
            'stored_geometry_contiguous': True, 'callback_trace_available': False,
            'sum_verified': False, 'source_enclosures_verified': False}


def main():
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
    assert len(manifest['files']) == len(REQUIRED)
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
                         'ancestor_called': False, 'original_source_or_quadrature_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'original_source_or_quadrature_called': False, 'fixture': None, 'saved_fixture_geometry': None}
    try:
        previous = json.loads(gzip.decompress(raw[RUNTIME_RECEIPT]))
        pins = [r['sha256'] for r in previous['source']['files'] if r['path'] == RUNTIME]
        assert pins == [hashes[RUNTIME]]
        out['runtime_source_matches_published_receipt'] = True
        out['saved_fixture_geometry'] = saved_fixture_geometry(raw[PRIOR])
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
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
