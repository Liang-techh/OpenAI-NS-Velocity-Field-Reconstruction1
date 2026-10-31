"""Single-entry, source-only AST witness and bounded selected-row comparison."""
import gzip
import hashlib
import importlib.util
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_moment_ast_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_moment_ast_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_moment_ast_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_moment_ast_v1_once.py'
PREFIX = 'experiments/root_st073/'
GRAPH = PREFIX + 'lei_ren_part1_paper_compliant_current_native_Rc_function_transport.py'
CONTROLS = PREFIX + 'lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py'
CORE_SOURCE = PREFIX + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py'
CORE = PREFIX + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.json.gz'
Q4 = PREFIX + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
PRIOR = PREFIX + 'original_q4_integrand_binding_v1_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, GRAPH, CONTROLS, CORE_SOURCE, CORE, Q4, PRIOR)
FIXED = {
 GRAPH: '57863c9c3e5746f770ece127cb233193d0e022735dd69994b8ca06c9a584d9d5',
 CONTROLS: 'cb6f11d4a62d59a45f05d0961236d5c6236bd544dca8b18499f68b2d8d98d009',
 CORE_SOURCE: '29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890',
 CORE: 'f930fb87cf9278bd37e6f9b11a6c2bf6431b531426c3eab43540733d61a5089f',
 Q4: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
 PRIOR: '2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2'}

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f: f.write(raw)
    return sha(raw)

def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None: raise RuntimeError('new module unavailable')
    out = importlib.util.module_from_spec(spec); spec.loader.exec_module(out)
    return out

def main():
    if ENTRY.exists() or RESULT.exists(): raise RuntimeError('one-entry component already entered')
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
        path = row['path']; raw[path] = (ROOT / path).read_bytes(); hashes[path] = sha(raw[path])
        assert hashes[path] == row['sha256']
        if path in FIXED: assert hashes[path] == FIXED[path]
    entry = save(ENTRY, {'component': COMPONENT, 'source_sha256': sha(manifest_raw),
                         'brief_sha256': sha(brief_raw), 'input_hashes': hashes,
                         'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry,
           'source_sha256': sha(manifest_raw), 'input_hashes': hashes,
           'raw_inputs_authenticated_before_decode': True, 'ancestor_called': False,
           'source_fixture': None, 'graph_fixture': None, 'source_result': None, 'archive_result': None}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        for key, got in (('source_fixture', fixture.run_fixture(leaf.check_source)),
                         ('graph_fixture', fixture.run_graph_fixture(leaf.bind))):
            out[key] = got
            if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
                raise ValueError(key + ' failed')
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert prior['status'] == 'returned' and prior['input_hashes'][Q4] == hashes[Q4]
        assert prior['archive_result']['source_archive_sha256'] == hashes[Q4]
        inspections = prior['archive_result']['inspections']
        assert len(inspections) == 13
        source = leaf.check_source({'graph': raw[GRAPH].decode('utf-8'), 'controls': raw[CONTROLS].decode('utf-8'),
                                    'core_inlet': raw[CORE_SOURCE].decode('utf-8')})
        out['source_result'] = source
        q4, core = json.loads(gzip.decompress(raw[Q4])), json.loads(gzip.decompress(raw[CORE]))
        comparison = leaf.bind(q4, core, inspections, source)
        links = {}
        for archive_name, archive in (('q4', q4), ('core', core)):
            links[archive_name] = {}
            for path in (GRAPH, CONTROLS, CORE_SOURCE, CORE):
                name = Path(path).name; recorded = archive.get('input_hashes', {}).get(name)
                links[archive_name][name] = {'recorded': recorded, 'current_sha256': hashes[path],
                                             'matches': recorded == hashes[path] if recorded is not None else None}
        out['archive_result'] = {'comparison': comparison, 'source_hash_links': links,
                                 'inspections_count': len(inspections), 'old_families_reentered': False}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'source_status': (out['source_result'] or {}).get('status'),
                      'comparison_status': ((out['archive_result'] or {}).get('comparison') or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': result_sha}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
