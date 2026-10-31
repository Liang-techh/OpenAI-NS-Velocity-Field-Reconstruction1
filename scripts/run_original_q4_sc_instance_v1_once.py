"""Single-entry new explicit Q4 sc instance; no ancestor code is invoked."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_sc_instance_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_sc_instance_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_sc_instance_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_sc_instance_v1_once.py'
SC_RECEIPT = 'experiments/root_st073/original_q4_archived_sc_endpoint_v1_receipt.json.gz'
FLAT_RECEIPT = 'experiments/root_st073/original_q4_slow_B_flat_v1_receipt.json.gz'
INNER = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_inner_exit_strict_collar.json'
INNER_CHECK = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_inner_exit_strict_collar_check.json'
REQUIRED = (LEAF, FIXTURE, RUNNER, SC_RECEIPT, FLAT_RECEIPT, INNER, INNER_CHECK)
FIXED = {
    SC_RECEIPT: '43ffdb58909074230587f721003a77ff0022cd433f514900347daeccbb838f6e',
    FLAT_RECEIPT: 'eaee3662dc501d65b282517e59c1ddb55cc91bcbc7bbda9703cc7faa2d27e92f',
    INNER: '53c619ff8f4e4bd31d6fa88f21b31f70020b278c3403f005daea91c482cc8163',
    INNER_CHECK: '043ad00271ae167e81046f89c82bb3ec07d3a64458d3d470980fdc6b14159d8d',
}
SC_HASH = 'a50d8fdc6728d9cd6ce0394b703025be77c68746608a20b6fa05cc4ff5be353a'
Q4 = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def save(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)

def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new module unavailable')
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('one-entry component already entered')
    brief_raw, manifest_raw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(brief_raw), json.loads(manifest_raw)
    assert brief['component'] == COMPONENT and brief['phase'] == COMPONENT + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(manifest_raw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
    rows = manifest['files']
    assert len(rows) == len(REQUIRED) and {row['path'] for row in rows} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in rows:
        path = row['path']
        raw[path] = (ROOT / path).read_bytes()
        hashes[path] = sha(raw[path])
        assert hashes[path] == row['sha256']
        if path in FIXED:
            assert hashes[path] == FIXED[path]
    entry = save(ENTRY, {'component': COMPONENT, 'source_sha256': sha(manifest_raw),
                         'brief_sha256': sha(brief_raw), 'input_hashes': hashes,
                         'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry,
           'source_sha256': sha(manifest_raw), 'input_hashes': hashes,
           'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False,
           'fixture': None, 'instance_result': None}
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.bind)
        out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('fresh synthetic fixture failed')
        sc = json.loads(gzip.decompress(raw[SC_RECEIPT]))
        flat = json.loads(gzip.decompress(raw[FLAT_RECEIPT]))
        source = json.loads(raw[INNER])
        check = json.loads(raw[INNER_CHECK])
        assert check['all_passed'] is True
        assert check['input_hashes'][Path(INNER).name] == hashes[INNER]
        assert source['actual_five_defect_family_sha256'] == check['actual_five_defect_family_sha256']
        assert source['cross_graph_attachment_sha256'] == check['cross_graph_attachment_sha256']
        prior = sc['result']; latest = flat['result']
        for path in (Q4, LOOP):
            assert prior['input_hashes'][path] == latest['input_hashes'][path]
        assert latest['archive_result']['status'] == 'verified'
        assert latest['archive_result']['scope']['Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive'] is True
        selected = source['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint']
        out['instance_result'] = leaf.bind(sc, selected, 10, SC_HASH)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'instance_status': (out['instance_result'] or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__':
    raise SystemExit(main())
