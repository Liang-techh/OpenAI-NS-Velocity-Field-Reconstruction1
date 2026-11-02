"""One entry of higher formal Kchi collar words on a saved core inlet."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'formal_kchi_collar_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_formal_axis_kchi_volterra_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_formal_axis_kchi_volterra_fixture_v1.py'
RUNNER = 'scripts/run_formal_kchi_collar_v1_once.py'
BASE = 'experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz'
KCHI = 'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz'
BIND = 'experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz'
FIRST = 'experiments/root_st073/n1_collar_core_inlet_potential_v5_receipt.json.gz'
OLDTAIL = 'experiments/root_st073/n1_collar_inlet_j2_tail_v3_receipt.json.gz'
SELECT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
REQUIRED = (LEAF, FIXTURE, RUNNER, BASE, KCHI, BIND, FIRST, OLDTAIL, SELECT)
PIN = {BASE: '59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d',
       KCHI: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       BIND: '5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',
       FIRST: '1c644fb0f41a8ae407c09f56b00e7aee1720f91d08163a0c7483b4b1b6c61e7a',
       OLDTAIL: 'dd28d5be7e090afc720f7b71d48ddbeb0e542e8e4d5097e86a86d78ab8c36891',
       SELECT: '0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}


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


def encode_intervals(value):
    if hasattr(value, '_mpi_'):
        return {'lower_exact_mpf_tuple': list(value._mpi_[0]),
                'upper_exact_mpf_tuple': list(value._mpi_[1]), 'directed': True}
    if isinstance(value, dict):
        return {str(k): encode_intervals(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode_intervals(v) for v in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError('cannot serialize result type ' + type(value).__name__)


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
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
                         'ancestor_called': False, 'previous_words_not_reentered': True})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'fixture': None, 'higher_collar_word': None}
    try:
        leaf = module(LEAF, NAME + '_leaf')
        fixture = module(FIXTURE, NAME + '_fixture')
        c = MPIntervalContext()
        c.prec = 1800
        out['fixture'] = encode_intervals(fixture.run_fixture(leaf.build, c))
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        decoded = {p: json.loads(gzip.decompress(raw[p])) for p in (BASE, KCHI, BIND, FIRST, OLDTAIL)}
        assert all(v['result']['status'] == 'returned' for v in decoded.values())
        base = decoded[BASE]['result']['output']
        kchi = decoded[KCHI]['result']['output']
        binding = decoded[BIND]['result']['output']
        selected = json.loads(raw[SELECT])
        assert base['source_family'] == kchi['source_family'] == binding['source_family'] == selected['source_family']
        assert base['source_logs'] == kchi['source_logs'] == binding['source_logs']
        assert binding['source_checks']['one_query_jet_shared_in_integration'] is True
        contract = selected['actual_cutoff_complete_support_measure']
        assert contract['passed'] is True and contract['measure_is_Z_independent'] is True
        assert '0<hb<1/2' in contract['inequality_proof'] and contract['source_family'] == base['source_family']
        previous_pins = {v['path']: v['sha256'] for v in decoded[FIRST]['source']['files']}
        assert previous_pins[BASE] == hashes[BASE] and previous_pins[KCHI] == hashes[KCHI]
        sc_bits = tuple(contract['source_sc_exact_mpf_tuple'])
        assert len(sc_bits) == 4
        sc = c.make_mpf((sc_bits, sc_bits))
        saved = base['baseline_G0g_fixed_exit']['G0g']['aggregate']
        inlet = {'F': {}, 'Fx': {}}
        index = {}
        for field, component in (('F1', 'F'), ('d_x_F1', 'Fx')):
            for order in ('Z0', 'Z1'):
                rows = []
                for i, item in enumerate(saved[order][field]):
                    path = ['baseline_G0g_fixed_exit', 'G0g', 'aggregate', order, field, i]
                    sid = 'saved_core_G0g:' + field + ':' + order + ':' + str(i)
                    origin = item.get('origin') if isinstance(item.get('origin'), dict) else {}
                    parent_ids = list(item.get('source_ids', [])) + list(origin.get('source_ids', []))
                    classes = list(item.get('error_classes', [])) + list(origin.get('error_classes', []))
                    powers = item.get('powers', [0, item.get('F0_power', 0), 0, item.get('Lambda_power', 0)])
                    rows.append({'coefficient': item['coefficient'], 'powers': powers,
                                 'source_ids': parent_ids + [sid], 'source_path': path,
                                 'error_classes': classes, 'uncertainty': True,
                                 'error_only': False, 'original_term': item})
                    index[sid] = {'receipt': BASE, 'path': path, 'parent_source_ids': parent_ids}
                assert rows
                inlet[component][order] = rows
        axis = []
        for order in (0, 1):
            axis.append({'coefficient': kchi['chi']['chi'][order], 'powers': [0, 0, 0, 0],
                         'source_ids': ['axis_chi:whole_Z', 'axis_chi:Z' + str(order)],
                         'source_path': ['chi', 'chi', order],
                         'error_classes': ['whole_Z_axis_potential_interval'], 'uncertainty': True,
                         'scope': 'same source axis Kchi, whole Z'})
        forcing = {'Z0': [], 'Z1': []}
        result = leaf.build(c, inlet, forcing, axis[0], axis[1], sc, max_order=2, min_order=2)
        out['higher_collar_word'] = encode_intervals(result)
        out['input_index'] = index
        out['saved_input_scope'] = {'family': base['source_family'], 'source_logs': base['source_logs'],
                                    'inlet_receipt': BASE, 'axis_receipt': KCHI,
                                    'axis_owner_evidence': BIND, 'previous_j1_receipt': FIRST,
                                    'previous_j2_tail_receipt': OLDTAIL,
                                    'forcing_supplied': False, 'same_live_owner_proved': False,
                                    'old_and_new_tail_combined': False}
        assert result['min_order'] == result['max_order'] == 2 and result['tail_starts_at_j'] == 3
        assert result['alternative_envelopes']['combined'] is False
        assert result['scope']['ledger_written'] is False
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
