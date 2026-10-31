"""One frozen entry for interval IR and strict directed-endpoint conversion."""
import dataclasses
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import sys
import traceback
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_interval_ir_runtime_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
P = 'experiments/root_st073/'
LEAF = 'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_fixture_v1.py'
CONVERTER = 'src/openai_ns_reconstruction/original_mpmath_directed_endpoints_v1.py'
CONVERTER_FIXTURE = 'src/openai_ns_reconstruction/original_mpmath_directed_endpoints_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_interval_ir_runtime_v1_once.py'
PRIOR = P + 'original_q4_moment_lowering_v4_receipt.json.gz'
KERNEL_RECEIPT = P + 'original_q4_fraction_interval_v2_receipt.json.gz'
ADAPTER_RECEIPT = P + 'original_q4_core_enclosure_adapter_v2_receipt.json.gz'
KERNEL = 'src/openai_ns_reconstruction/original_q4_fraction_interval_v2.py'
ADAPTER = 'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_v2.py'
REQUIRED = (LEAF, FIXTURE, CONVERTER, CONVERTER_FIXTURE, RUNNER, PRIOR,
            KERNEL_RECEIPT, ADAPTER_RECEIPT, KERNEL, ADAPTER)
PIN = {PRIOR: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
       KERNEL_RECEIPT: 'ce463dfa4e2bee2a81eebb4ad27da224b78d50dffa5162314c6fdeb9c444120d',
       ADAPTER_RECEIPT: '8afcf0f531fbbf458d3283ff11acc0d9160fb76abe2c6f26ad74fca73136b9ce'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode_unknown(value):
    if isinstance(value, Fraction):
        return {'type': 'Fraction', 'numerator': value.numerator, 'denominator': value.denominator}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {'type': type(value).__name__, **dataclasses.asdict(value)}
    raise TypeError('unsupported result value: ' + type(value).__name__)


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False,
                      default=encode_unknown) + '\n').encode()
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


def source_pin(receipt, path):
    rows = receipt['source']['files']
    matches = [row['sha256'] for row in rows if row['path'] == path]
    if len(matches) != 1:
        raise ValueError('missing prior source pin')
    return matches[0]


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
                         'ancestor_called': False, 'original_source_or_quadrature_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'original_source_or_quadrature_called': False, 'fixture': None,
           'conversion_fixture': None, 'actual_IR_gate': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        converter = module(CONVERTER, NAME + '_converter')
        converter_fixture = module(CONVERTER_FIXTURE, NAME + '_converter_fixture')
        out['fixture'] = fixture.run_fixture(leaf.evaluate_interval)
        assert out['fixture']['status'] == 'passed' and all(out['fixture']['checks'].values())
        out['conversion_fixture'] = converter_fixture.run_fixture(converter)
        assert out['conversion_fixture']['passed'] is True and all(out['conversion_fixture']['checks'].values())
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        kr = json.loads(gzip.decompress(raw[KERNEL_RECEIPT]))
        ar = json.loads(gzip.decompress(raw[ADAPTER_RECEIPT]))
        assert prior['status'] == kr['result']['status'] == ar['result']['status'] == 'returned'
        assert source_pin(kr, KERNEL) == hashes[KERNEL] and source_pin(ar, ADAPTER) == hashes[ADAPTER]
        assert len(prior['lowerings']) == 13
        first = prior['lowerings'][0]
        assert first['lowering']['status'] == 'lowered'
        gate = leaf.evaluate_interval(first['lowering'], oracle=None, backend=None, interval_type=None,
                                      extract_component=None, endpoint_converter=None, partitions=1)
        assert gate['status'] == 'rejected' and gate['value'] is None and gate['interval'] is None
        out['actual_IR_gate'] = {'moment': first['moment'], 'ordinary_Z_order': first['ordinary_Z_order'],
                                 'result': gate, 'saved_lowering_count': 13,
                                 'kernel_source_pin_matches': True, 'adapter_source_pin_matches': True,
                                 'old_kernel_or_adapter_imported': False, 'actual_source_value_computed': False,
                                 'mpmath_distribution_version_metadata': importlib.metadata.version('mpmath'),
                                 'actual_interval_object_inspected': False}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture': (out['fixture'] or {}).get('status'),
                      'conversion_fixture': (out['conversion_fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
