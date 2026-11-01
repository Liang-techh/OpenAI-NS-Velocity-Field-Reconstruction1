"""Read saved K0 W1 cells once, without rerunning the producer."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_w1_saved_cells_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_saved_k0_w1_cells_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_saved_k0_w1_cells_fixture_v1.py'
RUNNER = 'scripts/run_k0_w1_saved_cells_v1_once.py'
K0 = 'experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz'
CORE = 'experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json'
PRIOR = 'experiments/root_st073/scaled_core_m15_density_v1_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, K0, CORE, PRIOR)
PIN = {K0:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',
       CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',
       PRIOR:'eb0b3f0551ff2b860cfe26057741f53da41d4d925b018954435d9ece286e8a48'}
FAMILY_KEYS = ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f:
        f.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new component unavailable')
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == NAME and brief['phase'] == NAME + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert sha(Path(brief['source_checkpoint']['preimage_path']).read_bytes()) == brief['source_checkpoint']['preimage_sha256']
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert len(manifest['files']) == len(REQUIRED) and {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in manifest['files']:
        p = row['path']
        raw[p] = (ROOT / p).read_bytes()
        hashes[p] = sha(raw[p])
        assert hashes[p] == row['sha256']
        if p in PIN:
            assert hashes[p] == PIN[p]
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_gzip_authenticated_before_decompression': True,
                         'ancestor_called': False, 'new_source_values_or_integrals_computed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'new_source_values_or_integrals_computed': False, 'fixture': None, 'saved_binding': None,
           'saved_background': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert previous['status'] == 'returned'
        prior_family = previous['saved_scope']['source_family']
        core = json.loads(raw[CORE])
        assert all(core[k] == prior_family[k] for k in FAMILY_KEYS)
        archived = json.loads(gzip.decompress(raw[K0]))
        assert archived['result']['status'] == 'returned'
        rows = archived['source']['files']
        pins = [r['sha256'] for r in rows if r['path'] == CORE]
        assert pins == [hashes[CORE]]
        output = archived['result']['output']
        assert output['source_family'] == prior_family
        out['saved_binding'] = leaf.bind(output)
        whole = core['whole_core']
        out['saved_background'] = {
            'source_family': prior_family,
            'rho_domain': whole['rho'], 'Z_domain': whole['Z'],
            'U_lead_Z0': whole['Uz_derivatives']['rho0_Z0'],
            'U_lead_Z1': whole['Uz_derivatives']['rho0_Z1'],
            'Phi_Z0': whole['Phi_derivatives']['rho0_Z0'],
            'Phi_Z1': whole['Phi_derivatives']['rho0_Z1'],
            'model_tail_record': whole['explicit_model_factorial_tail'],
            'nonlinear_correction_record': whole['nonlinear_correction_bounds'],
            'core_hash_matched_original_profile_receipt': True,
            'Delta_U_packet_available': False,
            'F0_local_and_anchor_identity_proved': False,
            'new_enclosure_certification_performed': False,
        }
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
