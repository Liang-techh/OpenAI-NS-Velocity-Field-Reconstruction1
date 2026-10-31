"""Single-entry lexical analysis of previously saved Q4 substitution binders."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_binder_reachability_v3'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_binder_reachability_v3.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_binder_reachability_fixture_v3.py'
RUNNER = 'scripts/run_original_q4_binder_reachability_v3_once.py'
PRIOR = 'experiments/root_st073/original_q4_lexical_closure_v1_receipt.json.gz'
ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
FAILURE = 'experiments/root_st073/original_q4_binder_reachability_v1_failure.json.gz'
FAILURE_V2 = 'experiments/root_st073/original_q4_binder_reachability_v2_failure.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, ARCHIVE, FAILURE, FAILURE_V2)
FIXED = {
    FAILURE_V2: '84a5c6afabdada756a2d449898a6527c689817d3c35c6d1be2eec9b3d1d8aa46',
    FAILURE: 'e868d9ec7c3d6d1d11d8a7b44361e3dbe9a98856796284f769b0ab9b6c9d0a82',
    PRIOR: 'dc2738f64c4bf2d52546ce87778f5369e7527ce459d82c7d127a8366f257be4e',
    ARCHIVE: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new module unavailable')
    out = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out


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
    assert len(rows) == len(REQUIRED) and {r['path'] for r in rows} == set(REQUIRED)
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
           'fixture': None, 'archive_result': None}
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.analyze)
        out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('synthetic binder fixture failed')
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert prior['status'] == 'returned' and prior['input_hashes'][ARCHIVE] == hashes[ARCHIVE]
        assert prior['archive_result']['source_archive_sha256'] == hashes[ARCHIVE]
        previous_rows = prior['archive_result']['rows']
        assert len(previous_rows) == 3 and [r['port']['ordinary_Z_order'] for r in previous_rows] == [0, 1, 2]
        row_targets = []
        for row in previous_rows:
            assert row['status'] == 'candidate_closures' and row['scope']['old_errors_removed'] is False
            assert row['port']['archive_sha256'] == hashes[ARCHIVE]
            closures = row['unique_closures']
            assert closures and all(c['status'] == 'closure_created' for c in closures)
            refs = [c['node'] for c in closures]
            assert all(type(ref) is int and ref >= 0 for ref in refs) and len(refs) == len(set(refs))
            row_targets.append(refs)
        graph = json.loads(gzip.decompress(raw[ARCHIVE]))
        results = []
        for ref in sorted(set().union(*map(set, row_targets))):
            result = leaf.analyze(graph, ref)
            if not isinstance(result, dict) or result.get('status') not in ('complete', 'unresolved', 'mismatch'):
                raise ValueError('invalid binder analyzer result')
            results.append({'substitution_node': ref, 'analysis': result})
        out['archive_result'] = {'row_targets': row_targets, 'unique_analyses': results,
                                 'source_archive_sha256': hashes[ARCHIVE], 'prior_result_sha256': hashes[PRIOR],
                                 'old_families_reentered': False, 'numeric_field_value_computed': False}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    actual = out['archive_result'] or {}
    counts = {}
    for item in actual.get('unique_analyses', []):
        key = item['analysis']['status']; counts[key] = counts.get(key, 0) + 1
    print(json.dumps({'status': out['status'], 'analysis_statuses': counts,
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())


