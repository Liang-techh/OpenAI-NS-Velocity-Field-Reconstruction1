"""Single-entry expansion of previously opaque authenticated Q4 source rows."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_endpoint_expand_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_endpoint_expand_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_endpoint_expand_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_endpoint_expand_v1_once.py'
PRIOR = 'experiments/root_st073/original_q4_endpoint_resolver_v1_receipt.json.gz'
PORTS = 'experiments/root_st073/original_q4_sc_flat_instance_v1_receipt.json.gz'
ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, PORTS, ARCHIVE)
FIXED = {
    PRIOR: 'ee8ac7e6a410ab565f89f3e85b2d7040389ffceb2b97072328f10455918e6f17',
    PORTS: '1e4a3855f3894044303cfa36794f5bbd2e990cb5dd55bf34baa5b7024f0eb7a3',
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
        got = fixture.run_fixture(leaf.expand)
        out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('cases') or not all(x.get('passed') is True for x in got['cases']):
            raise ValueError('synthetic fixture failed')
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        port_result = json.loads(gzip.decompress(raw[PORTS]))['result']['instance_result']
        assert prior['status'] == 'returned' and prior['input_hashes'][ARCHIVE] == hashes[ARCHIVE]
        assert prior['input_hashes'][PORTS] == hashes[PORTS]
        assert prior['archive_result']['archive_sha256'] == hashes[ARCHIVE]
        assert port_result['status'] == 'verified' and port_result['scope']['old_errors_removed'] is False
        ports = port_result['source_ports']['U_y']
        saved = prior['archive_result']['results']
        assert len(ports) == len(saved) == 3 and [p['ordinary_Z_order'] for p in ports] == [0, 1, 2]
        assert all(p['archive_sha256'] == hashes[ARCHIVE] for p in ports)
        graph = json.loads(gzip.decompress(raw[ARCHIVE]))
        results = []
        for port, old in zip(ports, saved):
            result = leaf.expand(graph, old)
            results.append({'port': port, 'expansion': result})
        out['archive_result'] = {'rows': results, 'source_archive_sha256': hashes[ARCHIVE],
                                 'prior_result_sha256': hashes[PRIOR],
                                 'old_resolver_reentered': False, 'numerical_field_value_computed': False}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    actual = out['archive_result'] or {}
    print(json.dumps({'status': out['status'], 'row_statuses': [x['expansion']['status'] for x in actual.get('rows', [])],
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
