"""Enter only the new m15 density fixture, preserving source-scale evidence."""
import dataclasses
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'scaled_core_m15_density_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_core_scaled_m15_density_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_core_scaled_m15_density_fixture_v1.py'
RUNNER = 'scripts/run_scaled_core_m15_density_v1_once.py'
ADAPTER = 'src/openai_ns_reconstruction/n1_core_w1_m15_adapter_v2.py'
CORE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
PRIOR = 'experiments/root_st073/original_q4_core_saved_receipt_link_v1_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, ADAPTER, CORE, PRIOR)
PIN = {ADAPTER:'647dfd6e0280c9de1859d3e96c116bbe3ab92b6a83aab5d8336b64e76d1d5949',
       CORE:'d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       PRIOR:'d169fa93bde9ad8738abf69d3f601dab3321057b45cbee39febb8bcfb6ab957d'}


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
                         'ancestor_called': False, 'original_source_values_computed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'original_source_values_computed': False, 'fixture': None, 'saved_scope': None}
    try:
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert previous['status'] == 'returned'
        evidence = previous['saved_evidence']
        historical = evidence['historical_actual_oracle']
        out['saved_scope'] = {'source_family': evidence['source_family'],
                              'saved_sample_rho': historical['rho'], 'saved_sample_Z': historical['Z'],
                              'saved_point_used_as_full_cell': False,
                              'prior_component_graph_binding_proved': evidence['component_graph_binding_proved'],
                              'source_and_datum_certified_by_this_entry': False}
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
