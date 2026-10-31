"""One-entry bridge from selected Q4 source rows to saved inlet ranges."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_q4_saved_moment_join_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_saved_moment_join_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_saved_moment_join_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_saved_moment_join_v1_once.py'
P = 'experiments/root_st073/'
CORE = P + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.json.gz'
Q4 = P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
INTERIOR = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments.json.gz'
CHECK = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments_check.json'
PRIOR = P + 'original_q4_integrand_binding_v1_receipt.json.gz'
AST = P + 'original_q4_moment_ast_v1_receipt.json.gz'
PHYSICAL = P + 'original_core_physical_moment_binary_v2_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, CORE, Q4, INTERIOR, CHECK, PRIOR, AST, PHYSICAL)
FIXED = {
 CORE: 'f930fb87cf9278bd37e6f9b11a6c2bf6431b531426c3eab43540733d61a5089f',
 Q4: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
 INTERIOR: 'de8c97d25db958a692454b597116049739e03c2a97cab3800181fc8dd780be72',
 CHECK: 'cf1297a115b12cee236c97a8b837eac4d57e1e880d87c104170056e61cad70b9',
 PRIOR: '2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2',
 AST: '39c262869cb857a24afc66398f24d6d7f5af27e15f81be289150366d0ac3cae8',
 PHYSICAL: 'e068e2d66c61aa7b7c41e210e97e91c2943253f949c4eb462721a7f4f2a505f8'}

def sha(b): return hashlib.sha256(b).hexdigest()
def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f: f.write(raw)
    return sha(raw)
def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None: raise RuntimeError('new module unavailable')
    obj = importlib.util.module_from_spec(spec); sys.modules[name] = obj; spec.loader.exec_module(obj)
    return obj
def decode(raw): return json.loads(gzip.decompress(raw))

def main():
    if ENTRY.exists() or RESULT.exists(): raise RuntimeError('one-entry component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == COMPONENT and brief['phase'] == COMPONENT + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
    rows = manifest['files']
    assert len(rows) == len(REQUIRED) and {r['path'] for r in rows} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in rows:
        path = row['path']; raw[path] = (ROOT / path).read_bytes(); hashes[path] = sha(raw[path])
        assert hashes[path] == row['sha256']
        if path in FIXED: assert hashes[path] == FIXED[path]
    entry = save(ENTRY, {'component': COMPONENT, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'input_hashes': hashes, 'raw_inputs_authenticated_before_decode': True,
                         'ancestor_called': False, 'quadrature_or_recurrence_recomputed': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'input_hashes': hashes, 'raw_inputs_authenticated_before_decode': True,
           'ancestor_called': False, 'quadrature_or_recurrence_recomputed': False,
           'fixture': None, 'archive_result': None}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.bridge); out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('new synthetic moment join fixture failed')
        prior, ast, physical = decode(raw[PRIOR])['result'], decode(raw[AST])['result'], decode(raw[PHYSICAL])['result']
        assert prior['status'] == 'returned' and prior['input_hashes'][Q4] == hashes[Q4]
        assert ast['status'] == 'returned' and ast['source_result']['status'] == 'verified'
        assert ast['archive_result']['comparison']['status'] == 'matched'
        assert ast['input_hashes'][Q4] == hashes[Q4] and ast['input_hashes'][CORE] == hashes[CORE]
        assert physical['status'] == 'returned' and physical['input_hashes'][INTERIOR] == hashes[INTERIOR]
        check = json.loads(raw[CHECK]); assert check['all_passed'] is True
        assert check['input_hashes'][Path(INTERIOR).name] == hashes[INTERIOR]
        core, q4, interior = decode(raw[CORE]), decode(raw[Q4]), decode(raw[INTERIOR])
        out['archive_result'] = leaf.bridge(core, q4, prior['archive_result']['inspections'], interior, physical)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'bridge_status': (out['archive_result'] or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
