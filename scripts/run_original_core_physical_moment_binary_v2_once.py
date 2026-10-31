"""One-entry dyadic scaling of previously integrated inlet moments."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_core_physical_moment_binary_v2'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_core_physical_moment_binary_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_core_physical_moment_binary_v2_fixture_v2.py'
ADAPTER = 'src/openai_ns_reconstruction/original_core_saved_inlet_binary_v2.py'
RUNNER = 'scripts/run_original_core_physical_moment_binary_v2_once.py'
P = 'experiments/root_st073/'
ARCHIVE = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments.json.gz'
CHECK = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments_check.json'
FAILED = P + 'original_core_physical_moment_scale_v1_failed_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, ADAPTER, RUNNER, ARCHIVE, CHECK, FAILED)
FIXED = {
 ARCHIVE: 'de8c97d25db958a692454b597116049739e03c2a97cab3800181fc8dd780be72',
 CHECK: 'cf1297a115b12cee236c97a8b837eac4d57e1e880d87c104170056e61cad70b9',
 FAILED: '6fb9214c37f0d82ac996fb2a37cbc41047b46a43aecb3926a0ecbd42e35ee9fe'}

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
           'fixture': None, 'family_evidence': None, 'views': {}}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        adapter = module(ADAPTER, COMPONENT + '_adapter')
        got = fixture.run_fixture(leaf.scale); out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('new dyadic arithmetic fixture failed')
        prior = json.loads(gzip.decompress(raw[FAILED]))['result']
        assert prior['status'] == 'failed' and prior['error']['message'] == 'endpoint exceeds conversion limit'
        assert prior['fixture']['passed'] and prior['rounding_fixture']['passed']
        check = json.loads(raw[CHECK])
        assert check['input_hashes'][Path(ARCHIVE).name] == hashes[ARCHIVE]
        data = json.loads(gzip.decompress(raw[ARCHIVE]))
        family = {k: data[k] for k in adapter.FAMILY_FIELDS}
        evidence = prior['family_evidence']
        assert evidence['saved_core_family'] == family
        out['family_evidence'] = evidence
        for name in ('inlet', 'fresh_inlet'):
            packet = adapter.saved_inlet(data, check, name)
            physical = leaf.scale(packet['moments'], packet['log_derivative'], source_family=packet['family'])
            out['views'][name] = {'Z': packet['Z'], 'rho': packet['rho'],
                                  'input_moments': packet['moments'], 'F0_relative_Z1': packet['log_derivative'],
                                  'physical': physical}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'views': list(out['views']), 'entry_sha256': entry,
                      'result_sha256': digest, 'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
