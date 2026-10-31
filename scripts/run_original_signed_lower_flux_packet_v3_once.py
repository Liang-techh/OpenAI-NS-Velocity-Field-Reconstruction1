"""One-entry synthetic validation of the signed lower-packet adapter.

No original source provider is imported or run by this script.  Source inputs
are raw-hash pinned to make the proposed integration reviewable.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_signed_lower_flux_packet_v3'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE = DIRECTORY / 'source.json'
ENTRY = DIRECTORY / 'entry.json'
RESULT = DIRECTORY / 'result.json'
LEAF = 'src/openai_ns_reconstruction/original_signed_lower_flux_packet_v3.py'
FIXTURE = 'src/openai_ns_reconstruction/original_signed_lower_flux_packet_fixture_v3.py'
RUNNER = 'scripts/run_original_signed_lower_flux_packet_v3_once.py'
VIEW = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields.py'
QUERY = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_signed_functions.py'
RECOVERY = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
SCALED = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_native_signed_input_enclosures.py'
SEAM = 'src/openai_ns_reconstruction/original_n1_first_micro_seam.py'
PHYSICAL_SOURCE = 'docs/CURRENT_ORIGINAL_N1_FIRST_MICRO_SEAM_2026_10_17.md'
OLD_LEAF = 'src/openai_ns_reconstruction/original_signed_lower_flux_packet_v1.py'
OLD_RECEIPT = 'experiments/root_st073/original_signed_lower_flux_packet_v1_receipt.json.gz'
FAILED_RECEIPT = 'experiments/root_st073/original_signed_lower_flux_packet_v2_failure.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, VIEW, QUERY, RECOVERY, SCALED, SEAM, PHYSICAL_SOURCE, OLD_LEAF, OLD_RECEIPT, FAILED_RECEIPT)


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
        raise RuntimeError('this one-entry component was already entered')
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
    for row in rows:
        path = row['path']
        raw = (ROOT / path).read_bytes()
        assert sha(raw) == row['sha256']
        hashes[path] = row['sha256']
    entry_sha = save(ENTRY, {
        'component': COMPONENT, 'entered': True, 'source_sha256': sha(manifest_raw),
        'brief_sha256': sha(brief_raw), 'input_hashes': hashes,
        'actual_source_provider_called': False,
    })
    outcome = {
        'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry_sha,
        'source_sha256': sha(manifest_raw), 'input_hashes': hashes,
        'actual_source_provider_called': False, 'real_source_validation_claimed': False,
        'fixture': None,
    }
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.evaluate)
        outcome['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True:
            raise ValueError('synthetic fixture did not pass')
        checks = got.get('checks')
        if not isinstance(checks, dict) or not checks or not all(v is True for v in checks.values()):
            raise ValueError('synthetic fixture has missing or failed checks')
        outcome['status'] = 'returned'
    except BaseException as exc:
        outcome['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, outcome)
    print(json.dumps({'status': outcome['status'], 'entry_sha256': entry_sha, 'result_sha256': result_sha}))
    return 0 if outcome['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())


