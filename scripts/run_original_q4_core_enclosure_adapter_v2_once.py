"""One frozen entry: synthetic adapter tests and saved source-record preflight."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_core_enclosure_adapter_v2'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
P = 'experiments/root_st073/'
LEAF = 'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_fixture_v2.py'
RUNNER = 'scripts/run_original_q4_core_enclosure_adapter_v2_once.py'
PRIOR = P + 'original_q4_moment_lowering_v4_receipt.json.gz'
PROVIDER = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json'
FAILURE = P + 'original_q4_core_enclosure_adapter_v1_failure_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, PROVIDER, CHECK, FAILURE)
PIN = {PRIOR: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
       PROVIDER: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       CHECK: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead',
       FAILURE: '228e921e3a34d6d147e13e7c3fb35369f2ed45a1efb148a069dea1e40d111176'}


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
                         'ancestor_called': False, 'real_source_provider_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'real_source_provider_called': False, 'fixture': None, 'saved_source_preflight': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf.make_adapter)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        check = json.loads(raw[CHECK])
        assert prior['status'] == 'returned' and check['all_passed'] is True
        lowering = prior['lowerings'][0]['lowering']
        assert lowering['status'] == 'lowered' and check['source_family'] == lowering['source_family']
        assert check['input_hashes'][Path(PROVIDER).name] == hashes[PROVIDER]
        components = [a for a in lowering['source_atoms']
                      if a['op'] == 'current_original_unique_analytic_core_component']
        assert len(components) == 1
        record = components[0]['row']
        descriptor = record['source_provider']
        adapter = leaf.make_adapter(provider=None, source_family=lowering['source_family'],
                                    provider_descriptor=descriptor, degree=6, depth=3,
                                    interval_validator=lambda _: True)
        request = adapter._read_request(record)
        out['saved_source_preflight'] = {
            'quantity': request['quantity'], 'radial_order': request['radial_order'],
            'axial_order': request['axial_order'], 'provider_descriptor': descriptor,
            'provider_descriptor_matches_source_sha': isinstance(descriptor, dict)
                and descriptor.get('module_sha256') == hashes[PROVIDER],
            'source_family_matches_check_receipt': True,
            'provider_invoked': False, 'numeric_source_value_computed': False,
            'actual_caller_binding_certified': False,
        }
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
