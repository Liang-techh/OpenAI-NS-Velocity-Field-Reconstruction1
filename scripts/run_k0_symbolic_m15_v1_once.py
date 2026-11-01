"""Bind authenticated saved K0 cells to symbolic m15 density, once."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_symbolic_m15_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_k0_symbolic_m15_cells_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_k0_symbolic_m15_cells_fixture_v1.py'
RUNNER = 'scripts/run_k0_symbolic_m15_v1_once.py'
PRIOR = 'experiments/root_st073/k0_w1_saved_cells_v1_receipt.json.gz'
DOC = 'docs/CURRENT_N1_CORE_W0_M15_2026_10_30.md'
ADAPTER = 'src/openai_ns_reconstruction/n1_core_w1_m15_adapter_v2.py'
CUTOFF = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py'
CORE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, DOC, ADAPTER, CUTOFF, CORE)
PIN = {PRIOR:'a212575248d4eab13a11447f578d6c45205f70cde55a5c2d47f62fe8af50d752',
       DOC:'7d16e3c3b152d10301c63bb095b1e2d79e690fe7e563a104f3a6a165eae7102b',
       ADAPTER:'647dfd6e0280c9de1859d3e96c116bbe3ab92b6a83aab5d8336b64e76d1d5949',
       CUTOFF:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',
       CORE:'d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59'}


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
           'new_source_values_or_integrals_computed': False, 'fixture': None, 'real_source_expressions': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        previous = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert previous['status'] == 'returned'
        binding, background = previous['saved_binding'], previous['saved_background']
        assert binding['source_family'] == background['source_family']
        # Attach the precise source path to each previously saved interval;
        # keep the original record intact and do not reinterpret its endpoints.
        tagged_background = dict(background)
        for name, source_field in (('U_lead_Z0','Uz_derivatives/rho0_Z0'),
                                   ('U_lead_Z1','Uz_derivatives/rho0_Z1'),
                                   ('Phi_Z0','Phi_derivatives/rho0_Z0'),
                                   ('Phi_Z1','Phi_derivatives/rho0_Z1')):
            tagged_background[name] = {
                'interval': background[name],
                'source_id': 'core_physical_field:whole_core/' + source_field,
                'source_path': 'experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json',
                'source_receipt_sha256': hashes[PRIOR],
            }
        out['real_source_expressions'] = leaf.build(binding, tagged_background)
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
