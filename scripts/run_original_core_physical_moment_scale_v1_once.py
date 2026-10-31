"""Single-entry reuse of saved inlet moments with symbolic physical scales."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_core_physical_moment_scale_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_core_physical_moment_scale_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_core_physical_moment_scale_fixture_v1.py'
ADAPTER = 'src/openai_ns_reconstruction/original_core_inlet_saved_scale_adapter_v1.py'
RUNNER = 'scripts/run_original_core_physical_moment_scale_v1_once.py'
P = 'experiments/root_st073/'
ARCHIVE = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments.json.gz'
CHECK = P + 'lei_ren_part1_paper_compliant_current_core_interior_moments_check.json'
Q4 = P + 'original_q4_integrand_binding_v1_receipt.json.gz'
FIXED_CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json'
REQUIRED = (LEAF, FIXTURE, ADAPTER, RUNNER, ARCHIVE, CHECK, Q4, FIXED_CHECK)
FIXED = {
 ARCHIVE: 'de8c97d25db958a692454b597116049739e03c2a97cab3800181fc8dd780be72',
 CHECK: 'cf1297a115b12cee236c97a8b837eac4d57e1e880d87c104170056e61cad70b9',
 Q4: '2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2',
 FIXED_CHECK: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead'}

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

def wire(value, interval_type):
    if isinstance(value, Fraction): return str(value)
    if isinstance(value, interval_type): return [str(value.lower), str(value.upper)]
    if isinstance(value, dict): return {k: wire(v, interval_type) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [wire(v, interval_type) for v in value]
    return value

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
           'fixture': None, 'rounding_fixture': None, 'family_evidence': None, 'views': {}}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        adapter = module(ADAPTER, COMPONENT + '_adapter')
        got = fixture.run_fixture(leaf.scale); out['fixture'] = wire(got, leaf.RationalInterval)
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('synthetic physical scale fixture failed')
        unit = Fraction(1, 1 << adapter.BITS)
        tests = {'positive_extreme': adapter.outward_mpf_tuple([0,1,-10**18,1]) == (0,unit),
                 'negative_extreme': adapter.outward_mpf_tuple([1,1,-10**18,1]) == (-unit,0),
                 'positive_half_grid': adapter.outward_mpf_tuple([0,3,-adapter.BITS-1,2]) == (unit,2*unit),
                 'negative_half_grid': adapter.outward_mpf_tuple([1,3,-adapter.BITS-1,2]) == (-2*unit,-unit),
                 'exact_integer': adapter.outward_mpf_tuple([0,1,2,1]) == (4,4)}
        out['rounding_fixture'] = {'passed': all(tests.values()), 'checks': tests}
        if not all(tests.values()): raise ValueError('dyadic outward-rounding fixture failed')
        check = json.loads(raw[CHECK]); fixed = json.loads(raw[FIXED_CHECK])
        assert check['input_hashes'][Path(ARCHIVE).name] == hashes[ARCHIVE]
        assert fixed['all_passed'] is True
        data = json.loads(gzip.decompress(raw[ARCHIVE]))
        prior = json.loads(gzip.decompress(raw[Q4]))['result']
        assert prior['status'] == 'returned'
        family = {k: data[k] for k in adapter.FAMILY_FIELDS}
        core_rows = prior['archive_result']['core_sources']
        assert len(core_rows) == 3 and all(x['inspection']['status'] == 'source_ref' for x in core_rows)
        q4_families = [x['inspection']['evidence']['source_row']['source_family'] for x in core_rows]
        out['family_evidence'] = {'saved_core_family': family, 'fixed_point_family': fixed['source_family'],
             'all_three_Q4_core_family_records_match': all(x == family for x in q4_families),
             'fixed_point_family_matches': fixed['source_family'] == family,
             'source_family_match_is_not_capture_safe_substitution_proof': True}
        for name in ('inlet', 'fresh_inlet'):
            packet = adapter.saved_inlet(data, check, name)
            physical = leaf.scale(packet['moments'], packet['log_derivative'], source_family=packet['family'])
            out['views'][name] = {'Z': packet['Z'], 'rho': packet['rho'],
                                  'input_moments': packet['moments'], 'F0_relative_Z1': packet['log_derivative'],
                                  'physical': wire(physical, leaf.RationalInterval)}
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'views': list(out['views']), 'entry_sha256': entry,
                      'result_sha256': digest, 'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
