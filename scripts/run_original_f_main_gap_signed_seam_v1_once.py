"""Single-entry source-bound construction of the current F main/gap seam."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_f_main_gap_signed_seam_v1'
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_f_main_gap_signed_seam_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_f_main_gap_signed_seam_fixture_v1.py'
DUAL = 'src/openai_ns_reconstruction/original_f_main_gap_signed_seam_dual_fixture_v1.py'
RUNNER = 'scripts/run_original_f_main_gap_signed_seam_v1_once.py'
P = 'experiments/root_st073/'
F = P + 'lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz'
FCHECK = P + 'lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json'
OWN = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.json.gz'
OCHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields_check.json'
REQUIRED = (LEAF, FIXTURE, DUAL, RUNNER, F, FCHECK, OWN, OCHECK)
FIXED = {F: '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
         FCHECK: '320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15',
         OWN: '2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b',
         OCHECK: '0bc8229a1716871d896161caac907c572f7cd0a58540900aaa00869d3a920fb1'}

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
                         'ancestor_called': False, 'protected_ref_indexed': False})
    out = {'component': COMPONENT, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'input_hashes': hashes, 'raw_inputs_authenticated_before_decode': True,
           'ancestor_called': False, 'protected_ref_indexed': False, 'fixture': None, 'dual_fixture': None, 'archive_result': None}
    try:
        leaf = module(LEAF, COMPONENT); fixture = module(FIXTURE, COMPONENT + '_fixture')
        got = fixture.run_fixture(leaf.build); out['fixture'] = got
        if not isinstance(got, dict) or got.get('passed') is not True or not got.get('checks') or not all(v is True for v in got['checks'].values()):
            raise ValueError('synthetic source-bound seam fixture failed')
        dual = module(DUAL, COMPONENT + '_dual_fixture')
        dual_got = dual.run_fixture(leaf.build, fixture._archives); out['dual_fixture'] = dual_got
        if not isinstance(dual_got, dict) or dual_got.get('passed') is not True or not dual_got.get('checks') or not all(v is True for v in dual_got['checks'].values()):
            raise ValueError('independent synthetic dual-number fixture failed')
        fcheck, ocheck = json.loads(raw[FCHECK]), json.loads(raw[OCHECK])
        assert fcheck['all_passed'] is True and ocheck['all_passed'] is True
        assert fcheck['input_hashes'][Path(F).name] == hashes[F]
        assert fcheck['input_hashes'][Path(OWN).name] == hashes[OWN]
        assert fcheck['input_hashes'][Path(OCHECK).name] == hashes[OCHECK]
        assert ocheck['input_hashes'][Path(OWN).name] == hashes[OWN]
        f, own = json.loads(gzip.decompress(raw[F])), json.loads(gzip.decompress(raw[OWN]))
        assert f['source_family'] == own['source_family'] == fcheck['source_family'] == ocheck['source_family']
        assert f['input_hashes'][Path(OWN).name] == hashes[OWN]
        assert f['input_hashes'][Path(OCHECK).name] == hashes[OCHECK]
        out['archive_result'] = leaf.build(f, own)
        assert isinstance(out['archive_result'], dict) and out['archive_result'].get('status') in ('constructed', 'unresolved')
        assert all(out['archive_result']['scope'][k] is False for k in
                   ('source_field_radius_owner_identity', 'capture_safe_endpoint', 'numerical_interval_result',
                    'seam_cancellation', 'ordinary_Z_commutation', 'm15_closure', 'old_errors_removed'))
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'seam_status': (out['archive_result'] or {}).get('status'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1

if __name__ == '__main__': raise SystemExit(main())
