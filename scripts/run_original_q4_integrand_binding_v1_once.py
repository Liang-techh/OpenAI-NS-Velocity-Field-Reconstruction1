"""Single-entry inspection of saved Q4 core-inlet integral frontiers."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_integrand_binding_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_integrand_binding_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_integrand_binding_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_integrand_binding_v1_once.py'
PRIOR = 'experiments/root_st073/original_q4_binder_reachability_v3_receipt.json.gz'
ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, ARCHIVE)
FIXED = {
    PRIOR: '94f05f4f580e99f7e630ed186113fa852d13b17f5e8bc72bfdc0c31886a02ceb',
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
           'fixture': None, 'core_fixture': None, 'archive_result': None}
    try:
        leaf = module(LEAF, COMPONENT)
        fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.inspect)
        out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('synthetic integrand fixture failed')
        core_got = fixture.run_core_fixture(leaf.inspect_core)
        out['core_fixture'] = core_got
        if not isinstance(core_got, dict) or core_got.get('passed') is not True or not core_got.get('checks') or not all(v is True for v in core_got['checks'].values()):
            raise ValueError('synthetic core source fixture failed')
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert prior['status'] == 'returned' and prior['input_hashes'][ARCHIVE] == hashes[ARCHIVE]
        assert prior['archive_result']['source_archive_sha256'] == hashes[ARCHIVE]
        frontiers, cores = {}, {}
        for item in prior['archive_result']['unique_analyses']:
            assert item['analysis']['scope']['complete_identity_proved'] is False
            for frontier in item['analysis']['frontier']:
                if frontier.get('operation') == 'current_original_unique_analytic_core_component':
                    core = dict(frontier)
                    core['expected_binder_node'] = item['analysis']['target']['binder_node_id']
                    ref = core['node_id']
                    if ref in cores:
                        assert cores[ref] == core
                    cores[ref] = core
                if frontier.get('operation') != 'definite_integral':
                    continue
                ref = frontier['node_id']
                assert type(ref) is int and ref >= 0
                if ref in frontiers:
                    assert frontiers[ref] == frontier
                frontiers[ref] = frontier
        assert frontiers and cores
        graph = json.loads(gzip.decompress(raw[ARCHIVE]))
        results = []
        for ref in sorted(frontiers):
            result = leaf.inspect(graph, frontiers[ref])
            if not isinstance(result, dict) or result.get('status') not in ('observed', 'mismatch'):
                raise ValueError('invalid integral inspection status')
            results.append({'frontier': frontiers[ref], 'inspection': result})
        core_results = []
        for ref in sorted(cores):
            result = leaf.inspect_core(graph, cores[ref])
            if not isinstance(result, dict) or result.get('status') not in ('source_ref', 'mismatch'):
                raise ValueError('invalid core source status')
            core_results.append({'frontier': cores[ref], 'inspection': result})
        out['archive_result'] = {'inspections': results, 'source_archive_sha256': hashes[ARCHIVE],
                                 'core_sources': core_results,
                                 'prior_result_sha256': hashes[PRIOR], 'old_families_reentered': False,
                                 'numeric_field_value_computed': False}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    actual = out['archive_result'] or {}
    counts = {}
    for item in actual.get('inspections', []):
        key = item['inspection']['status']; counts[key] = counts.get(key, 0) + 1
    core_counts = {}
    for item in actual.get('core_sources', []):
        key = item['inspection']['status']; core_counts[key] = core_counts.get(key, 0) + 1
    print(json.dumps({'status': out['status'], 'inspection_statuses': counts, 'core_statuses': core_counts,
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
