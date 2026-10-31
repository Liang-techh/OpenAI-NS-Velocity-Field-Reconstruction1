"""Single-entry archived Q4 slow B flat source consumer.

Raw inputs are hashed before gzip decoding.  No original provider is imported,
constructed or invoked.  Synthetic fixture result and real archive outcome are
saved before they are reported.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_slow_B_flat_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE = DIRECTORY / 'source.json'
ENTRY = DIRECTORY / 'entry.json'
RESULT = DIRECTORY / 'result.json'
LEAF = 'src/openai_ns_reconstruction/original_q4_slow_B_flat_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_slow_B_flat_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_slow_B_flat_v1_once.py'
Q4_ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP_ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
Q4_SOURCE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
LOOP_SOURCE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
PRIOR_SOURCE = 'src/openai_ns_reconstruction/original_native_lower_physical_flux_identity.py'
PRIOR_REPORT = 'experiments/root_st073/original_native_lower_physical_flux_identity_execution.json'
FLAT_RECEIPT = 'experiments/root_st073/original_q4_archived_flat_collar_v1_receipt.json.gz'
PHYSICAL_RECEIPT = 'experiments/root_st073/original_q4_archived_physical_rows_v1_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, Q4_ARCHIVE, LOOP_ARCHIVE, Q4_SOURCE, LOOP_SOURCE, PRIOR_SOURCE, PRIOR_REPORT, FLAT_RECEIPT, PHYSICAL_RECEIPT)
FIXED = {
    Q4_ARCHIVE: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP_ARCHIVE: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    PRIOR_REPORT: 'af9ad418486f96462e38289ffbae2afff4b426735418f978064c8e9ccbec58c2',
    FLAT_RECEIPT: 'bc519a2825141521b3e805d72ca276f697191521641556ed8e2c04e48f7b3a87',
    PHYSICAL_RECEIPT: '083bf7fb02519b498023da86806b3db7bbbd5940d14fcbeca43b45edd5afba99',
}
DECODED = {
    Q4_ARCHIVE: '0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75',
    LOOP_ARCHIVE: 'f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, data):
    raw = (json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new module cannot be loaded')
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('one-entry component already entered')
    brief_raw = BRIEF.read_bytes()
    manifest_raw = SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    manifest = json.loads(manifest_raw)
    assert brief['component'] == COMPONENT and brief['phase'] == COMPONENT + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['one_entry_only'] is True and brief['do_not_reenter'] is True
    assert manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert brief['source_checkpoint']['sha256'] == sha(manifest_raw)
    rows = manifest['files']
    assert len(rows) == len(REQUIRED) and {r['path'] for r in rows} == set(REQUIRED)
    hashes = {}
    raw_inputs = {}
    for row in rows:
        path = row['path']
        raw = (ROOT / path).read_bytes()
        assert sha(raw) == row['sha256']
        if path in FIXED:
            assert sha(raw) == FIXED[path]
        hashes[path] = row['sha256']
        raw_inputs[path] = raw
    entry_sha = save(ENTRY, {
        'component': COMPONENT, 'entered': True, 'source_sha256': sha(manifest_raw),
        'brief_sha256': sha(brief_raw), 'input_hashes': hashes,
        'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False,
    })
    outcome = {
        'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry_sha,
        'source_sha256': sha(manifest_raw), 'input_hashes': hashes,
        'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False,
        'fixture': None, 'archive_result': None, 'decoded_sha256': {},
    }
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.verify)
        outcome['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True:
            raise ValueError('synthetic fixture did not pass')
        checks = got.get('checks')
        if not isinstance(checks, dict) or not checks or not all(v is True for v in checks.values()):
            raise ValueError('synthetic fixture has missing or failed checks')
        decoded = {}
        for path in (Q4_ARCHIVE, LOOP_ARCHIVE):
            source_bytes = gzip.decompress(raw_inputs[path])
            outcome['decoded_sha256'][path] = sha(source_bytes)
            if sha(source_bytes) != DECODED[path]:
                raise ValueError('decoded archive SHA mismatch: ' + path)
            decoded[path] = json.loads(source_bytes)
        prior_flat = json.loads(gzip.decompress(raw_inputs[FLAT_RECEIPT]))
        prior_physical = json.loads(gzip.decompress(raw_inputs[PHYSICAL_RECEIPT]))
        actual = leaf.verify(decoded[Q4_ARCHIVE], decoded[LOOP_ARCHIVE], prior_flat, prior_physical)
        outcome['archive_result'] = actual
        if not isinstance(actual, dict) or actual.get('status') not in ('verified', 'mismatch', 'unresolved'):
            raise ValueError('archive consumer returned invalid status')
        outcome['status'] = 'returned'
    except BaseException as exc:
        outcome['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, outcome)
    print(json.dumps({'status': outcome['status'], 'archive_status': (outcome['archive_result'] or {}).get('status'),
                      'entry_sha256': entry_sha, 'result_sha256': result_sha}))
    return 0 if outcome['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())


