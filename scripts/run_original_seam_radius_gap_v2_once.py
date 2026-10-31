"""One frozen entry: conditional seam differences tied to saved evidence."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_seam_radius_gap_v2'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_seam_radius_gap_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_seam_radius_gap_fixture_v2.py'
RUNNER = 'scripts/run_original_seam_radius_gap_v2_once.py'
RECEIPT = 'experiments/root_st073/original_f_main_gap_signed_seam_v1_receipt.json.gz'
FAILURE = 'experiments/root_st073/original_seam_radius_gap_v1_failure_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, FAILURE, RECEIPT)
EXPECTED = 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as out:
        out.write(raw)
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
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw = {}
    for row in manifest['files']:
        raw[row['path']] = (ROOT / row['path']).read_bytes()
        assert sha(raw[row['path']]) == row['sha256']
    assert sha(raw[RECEIPT]) == EXPECTED
    assert sha(raw[FAILURE]) == '73e5a8e06601fde2cc01fe5beb329c94829cdf8a8598ec9b4197303208cc3b9f'
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw),
                         'brief_sha256': sha(braw), 'raw_authenticated_before_decode': True,
                         'ancestor_called': False, 'protected_graph_refs_indexed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry,
           'source_sha256': sha(mraw), 'prior_receipt_sha256': EXPECTED,
           'ancestor_called': False, 'protected_graph_refs_indexed': False,
           'fixture': None, 'sides': {}}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf.build)
        assert out['fixture']['passed'] is True
        assert all(v is True for v in out['fixture']['checks'].values())
        saved = json.loads(gzip.decompress(raw[RECEIPT]))['result']
        assert saved['status'] == 'returned'
        prior = saved['archive_result']
        assert prior['status'] == 'constructed'
        for side, chart in (('main.upper', 'main'), ('gap.lower', 'gap')):
            ep = prior['endpoints'][side]
            old = prior['candidate_old_C3_endpoints']['selected_endpoint_refs'][chart]
            # Ports carry labels, not values. The mapping below is only a proposed
            # provenance bridge: no equality or endpoint transport is inferred.
            ports = {key: {'role': side + '.' + key} for key in leaf.PORT_NAMES}
            packet = leaf.build(ports)
            out['sides'][side] = {
                'conditional_gap': packet,
                'prior_current_F_endpoint_candidate': {k: ep[k] for k in
                    ('R', 'R_Z', 'U', 'U_Z1', 'Q', 'Q_Z1', 'Q_y', 'Q_y_Z1')},
                'prior_old_C3_composite_endpoint_refs': old,
                'port_binding_proposals': {
                    'Rf_Rf_Z': 'prior current F R and formal R_Z; derivative semantics unresolved',
                    'Hf_Hf_Z': 'candidate from current F Q,U,delta,Z and their ordinary derivatives',
                    'qyf_qyf_Z': 'candidate from current F Q_y and its saved Z1 row',
                    'Ho_Ho_Z_qyo_qyo_Z': 'old source-field and endpoint-product transport unresolved',
                    'D': 'log(R_old/Rf), including old/current coordinate and reference-scale differences',
                    'D_Z': 'ordinary derivative of the full D, not assumed zero',
                },
                'binding_proposals_are_proofs': False,
                'old_to_F_endpoint_identity_proved': False,
            }
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture_passed': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
