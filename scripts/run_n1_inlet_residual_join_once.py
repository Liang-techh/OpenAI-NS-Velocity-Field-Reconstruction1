"""One-entry join of saved raw W1 and native first-inlet operator sectors.

No ancestor module is imported or executed. Raw input hashes are checked
before the selected JSON and compressed operator archive are decoded.
"""
from __future__ import annotations

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
COMPONENT = 'n1_inlet_residual_join_v1'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / x for x in ('source.json', 'entered.json', 'result.json'))
BASE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_'
REPORT = BASE + 'selected_words.json'
OPERATOR = BASE + 'inner_operator.json.gz'
CHECK = BASE + 'inner_operator_check.json'
OP_SOURCE = BASE + 'inner_operator.py'
COUPLED = BASE + 'coupled_first_step.py'
VOLTERRA = BASE + 'inner_volterra_functions.py'
NUMERIC = 'src/openai_ns_reconstruction/original_n1_selected_word_numeric.py'
IDENTITY = 'src/openai_ns_reconstruction/original_n1_first_volterra_remainder.py'
LEAF = 'src/openai_ns_reconstruction/n1_inlet_word_residual.py'
RUNNER = 'scripts/run_n1_inlet_residual_join_once.py'
REQUIRED = (REPORT, OPERATOR, CHECK, OP_SOURCE, COUPLED, VOLTERRA, NUMERIC, IDENTITY, LEAF, RUNNER)
FIXED = {
    REPORT: '0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',
    OPERATOR: '63b181010f1715d3b20589b06c1c2139ef63bcc73d7b0f3c5b7c43d9c5e03f52',
    CHECK: '846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4',
    COUPLED: '97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7',
    VOLTERRA: 'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',
    OP_SOURCE: '268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',
}
FIELDS = ('F1', 'Uz1', 'K1', 'P1', 'd_x_F1', 'd_x_Uz1')


def sha(raw): return hashlib.sha256(raw).hexdigest()
def require(condition, message):
    if not condition: raise ValueError(message)


def exclusive(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()
    with path.open('xb') as handle: handle.write(raw)
    return sha(raw)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(), 'component already entered')
    brief_raw = BRIEF.read_bytes(); brief = json.loads(brief_raw)
    require(brief.get('component') == COMPONENT and brief.get('phase') == COMPONENT + '_pending'
            and brief.get('goal_active') is True and brief.get('full_goal_complete') is False
            and brief.get('one_entry_only') is True and brief.get('do_not_reenter') is True, 'brief guard')
    checkpoint = brief['source_checkpoint']; source_raw = SOURCE.read_bytes()
    require(Path(checkpoint['path']).resolve() == SOURCE.resolve() and sha(source_raw) == checkpoint['sha256'], 'source checkpoint')
    source = json.loads(source_raw)
    require(source['component'] == COMPONENT and source['one_entry_only'] is True, 'source manifest scope')
    hashes = {}
    for row in source['files']:
        path = row['path']; require(path in REQUIRED and path not in hashes, 'unexpected input')
        digest = sha((ROOT / path).read_bytes())
        require(digest == row['sha256'] and (path not in FIXED or digest == FIXED[path]), 'input hash changed: '+path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), 'missing input')
    return hashes, sha(source_raw), sha(brief_raw)


def record(c, value):
    lo = tuple(value['lower_exact_mpf_tuple']); hi = tuple(value['upper_exact_mpf_tuple'])
    require(len(lo) == len(hi) == 4, 'four-entry exact interval endpoints required')
    low, high = c.make_mpf((lo, lo)), c.make_mpf((hi, hi))
    require(low <= high, 'reversed exact interval')
    return {'lower': low, 'upper': high}


def serialize(value):
    if hasattr(value, '_mpi_'):
        return {'lower_exact_mpf_tuple': list(value._mpi_[0]), 'upper_exact_mpf_tuple': list(value._mpi_[1]), 'directed': True}
    if isinstance(value, dict):
        return {(','.join(str(v) for v in k) if isinstance(k, tuple) else str(k)): serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [serialize(v) for v in value]
    if isinstance(value, (str, int, bool)) or value is None: return value
    raise TypeError('unsupported result type '+type(value).__name__)


def consume(hashes):
    raw = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes(); require(sha(data) == hashes[path], 'input changed after entry: '+path); raw[path] = data
    # The report and operator are decoded only after the raw-byte SHA checks.
    report = json.loads(raw[REPORT]); archive = json.loads(gzip.decompress(raw[OPERATOR])); check = json.loads(raw[CHECK])
    require(report['checks']['passed'] is True and check['all_passed'] is True, 'accepted input checks')
    require(check['input_hashes'][Path(OPERATOR).name] == hashes[OPERATOR], 'operator check does not bind archive')
    enclosure = report['actual_selected_first_word_enclosure']; sample = archive['samples']['first_inlet']
    require(enclosure['passed'] is True and enclosure['ordinary_Z_orders'] == [0, 1]
            and enclosure['core_quadrature_and_source_interval_errors_retained'] is True
            and enclosure['actual_directed_sector_coefficients_and_local_radii_evaluated'] is True, 'selected finite-word receipt')
    require(report['source_family'] == enclosure['source_family'] == archive['source_family'] == sample['source_family'], 'source family mismatch')
    require(sample['Z'] == enclosure['selected_Z_interval'], 'selected Z tube mismatch')
    require(sample['source_log_basis_names'] == enclosure['source_log_basis_names'] == ['hb', 'F0_anchor', 'Pstar_squared', 'Lambda']
            and sample['source_log_bases'] == enclosure['source_log_bases'], 'source anchor mismatch')
    require(sample['amplitude_anchor_fixed_during_jet_algebra'] is True
            and sample['ordinary_Z_amplitude_derivatives_in_coefficients'] is True
            and sample['outputs_are_enclosures_of_original_functions_not_selected_point_values'] is True, 'operator jet semantics')
    require(sample['geometry']['chart'] == 'first', 'not native first chart')
    require(enclosure['source_pins'][COUPLED]['expected_sha256'] == hashes[COUPLED]
            and enclosure['source_pins'][COUPLED]['passed'] is True, 'coupled source pin')
    source_text = {key: raw[key].decode('utf-8') for key in (NUMERIC, COUPLED, VOLTERRA, IDENTITY)}
    require(all(fragment in source_text[NUMERIC] for fragment in (
        'sample = archive["samples"]["core_exit"]', 'sample["computed_preconditioned_indices"] == [0, 1]',
        'inlets = {order: dict(W0=core_fields', '("B0_feedback", "B1_derivative_feedback")',
        'w1 = transport(c, inlets[order]["W1"], y_upper)')), 'word source binding')
    require(all(fragment in source_text[COUPLED] for fragment in (
        'Wbase=(I+J)Gg; H=(I+J)G(B0-Akin+B1*d_Z)',
        'B0_feedback={name:op.export(value,1)', 'B1_derivative_feedback={name:op.export(value,1)')), 'coupled word definition')
    require('L(sum_0^N Wk)-g=-(C0*WN+B1*dZ(WN))' in source_text[IDENTITY]
            and "C0={key:q for key,q in matrices['B0'].items() if key not in ('1,5','2,6','3,6')}" in source_text[VOLTERRA], 'canonical residual identity')
    q0 = [row for row in enclosure['selected_coordinate_rows'] if row['q'] == '0']
    require(len(q0) == 1, 'unique q=0 saved row required'); q0 = q0[0]
    c = MPIntervalContext(); c.prec = 1800
    for key in ('phase', 'rho_enclosure_only'):
        r = record(c, sample['geometry'][key]); target = c.mpf(0 if key == 'phase' else 4)
        require(r['lower'] == target and r['upper'] == target, 'inlet coordinate mismatch: '+key)
    y = record(c, q0['y_upper']); require(y['lower'] == 0 and y['upper'] == 0, 'q0 transport not exact seam')
    W1 = {'Z0': {}, 'Z1': {}}
    for order in (0, 1):
        row = q0['ordinary_Z'][str(order)]
        require(row['allword_tail_included'] is False, 'saved finite word contains omitted tail')
        for col, name in enumerate(FIELDS, 1):
            error = record(c, row['local_raw_error_radii']['W1'][name])
            require(error['lower'] == 0 and error['upper'] == 0, 'nonzero q0 local radius')
            W1['Z'+str(order)][col] = [dict(source_exponents=t['source_exponents'],
                combined_source_log_interval=record(c, t['combined_source_log_interval']),
                coefficient_interval=record(c, t['coefficient_interval'])) for t in row['W1'][name]]
    matrices = {}
    for matrix in ('B0', 'B1'):
        matrices[matrix] = {'Z0': {}, 'Z1': {}}
        for key, row in sample['regular_system'][matrix].items():
            pair = tuple(int(x) for x in key.split(',')); require(len(pair) == 2, 'matrix index')
            for order in (0, 1):
                matrices[matrix]['Z'+str(order)][pair] = [dict(source_exponents=t['source_exponents'],
                    combined_source_log_interval=record(c, t['combined_source_log']),
                    coefficient_interval=record(c, t['ordinary_Z_derivatives'][order])) for t in row['terms']]
    spec = importlib.util.spec_from_file_location('n1_inlet_residual_leaf', ROOT / LEAF)
    require(spec is not None and spec.loader is not None, 'leaf import'); module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module; spec.loader.exec_module(module)
    output = module.compute_inlet_word_residual(B0=matrices['B0'], B1=matrices['B1'], W1=W1)
    require(output['ordinary_Z']['Z1_fully_evaluated'] is False, 'missing W1_ZZ misrepresented')
    return {'source_family': report['source_family'], 'selected_Z_interval': enclosure['selected_Z_interval'],
        'input_binding': {'same_source_family': True, 'same_Z_tube': True, 'same_four_log_anchors': True,
          'rho_is_4_and_phase_zero': True, 'q0_local_W1_radius_zero': True,
          'core_cell_quadrature_intervals_retained': True, 'W1_Z2_unavailable': True},
        'residual_enclosure': serialize(output), 'all_Z_functional_coverage': False,
        'actual_allword_or_NS_residual': False, 'signed_m15_integral_computed': False}


def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {'component': COMPONENT, 'source_sha256': source_sha,
        'brief_sha256': brief_sha, 'input_hashes': hashes, 'entered': True})
    result = {'component': COMPONENT, 'source_sha256': source_sha, 'entry_sha256': entry_sha, 'status': 'failed'}
    try:
        result['output'] = consume(hashes); result['status'] = 'returned'; result['error'] = None
    except BaseException as exc:
        result['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = exclusive(RESULT, result)
    print(json.dumps({'component': COMPONENT, 'status': result['status'], 'result_sha256': digest}))
    return 0 if result['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
