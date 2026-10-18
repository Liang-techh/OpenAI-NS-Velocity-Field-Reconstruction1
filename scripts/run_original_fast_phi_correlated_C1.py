"""New correlated fast-row calculation; saved inspection is nonnumerical."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
import time

from mpmath.libmp import mpf_cmp


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'experiments/root_st073'
RUNNER = 'scripts/run_original_fast_phi_correlated_C1.py'
PARENT = 'lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json'
PARENT_SHA = '78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb'
PARENT_BYTES = 3159660
NAME = 'lei_ren_part1_paper_compliant_current_original_fast_phi_correlated_C1.json'
REPORT = BASE / NAME
CHECKPOINT = BASE / (NAME + '.new_math_checkpoint')
MODULES = (
    'src/openai_ns_reconstruction/original_post_restore_fast_phi_identity.py',
    'src/openai_ns_reconstruction/original_post_restore_fast_phi_bound.py',
    'src/openai_ns_reconstruction/original_post_restore_closed_phase.py',
)
FIXTURES = (
    'original_post_restore_fast_phi_identity_fixture.json',
    'original_post_restore_fast_phi_bound_fixture.json',
    'original_post_restore_closed_phase_fixture.json',
)
REVIEWED = {
    MODULES[0]: 'b1d77a1839edd86713263a781068da3d583cfa4c2cbdad2cf81d8ebf6850692a',
    MODULES[1]: '9cd734153944d11c1725ad22154a1dce61aaddc022803795eb88a8a0a54d75f8',
    MODULES[2]: '25b9a21c20ad8d0c005ad5ff42c5500cc56e2463cd0535e26ae7f03349a4163f',
}
FALSE_SCOPE = (
    'actual_signed_native_p2_or_phase_values_installed', 'signed_pressure_join_completed',
    'independent_reserved_full_pulse_m15_matched', 'global_five_moment_closure',
    'exact_heat_exterior', 'admissible_stress_and_flat_remainder', 'true_n_dependent_recursion',
    'realizable_oscillatory_correction', 'independent_Cartesian_NS_validation',
    'legacy_packets_recertified', 'inherited_global_guard_admitted', 'full_goal_complete',
)
C5 = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
C5_SHA = '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d'
C5_AST = 'fffc175d358e7cce71d3046ec10303e22d46b5ebb3145f9a44c0d1fe88796f8c'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def sha(path):
    return digest(Path(path).read_bytes())


def input_path(name):
    path = Path(name)
    require(not path.is_absolute(), 'Absolute source path')
    if '..' in path.parts:
        require(name == '../../src/openai_ns_reconstruction/schedule_pressure.py', 'Source traversal')
        return (BASE / path).resolve()
    return ROOT / path if path.parts[0] == 'experiments' else BASE / path.name


def preflight():
    for name in MODULES + (RUNNER,):
        raw = (ROOT / name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'),
                'New Python must use LF and one EOF newline: ' + name)
        ast.parse(raw, filename=name)
    for name, expected in REVIEWED.items():
        require(expected != 'REVIEW_PENDING' and sha(ROOT / name) == expected, 'Review pending/changed: ' + name)


def parent():
    raw = (BASE / PARENT).read_bytes()
    require(len(raw) == PARENT_BYTES and digest(raw) == PARENT_SHA, 'Parent changed before decoding')
    parsed = json.loads(raw)
    for name, expected in parsed['input_hashes'].items():
        require(sha(input_path(name)) == expected, 'Input changed: ' + name)
    for name, expected in parsed['implementation_hashes'].items():
        path = Path(name)
        require(not path.is_absolute() and '..' not in path.parts, 'Unsafe implementation path')
        require(sha(ROOT / name) == expected, 'Implementation changed: ' + name)
    return raw, parsed


def endpoint(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4 and all(type(x) is int for x in value),
            'Malformed exact MPF tuple')
    sign, mantissa, exponent, bits = value
    require(sign in (0, 1) and mantissa >= 0 and bits == mantissa.bit_length(), 'Nonfinite endpoint')
    require((mantissa == 0 and tuple(value) == (0, 0, 0, 0)) or (mantissa > 0 and mantissa & 1),
            'Noncanonical endpoint')
    return tuple(value)


def record(value):
    require(isinstance(value, dict) and value.get('directed') is True, 'Directed record missing')
    require(mpf_cmp(endpoint(value['lower_exact_mpf_tuple']), endpoint(value['upper_exact_mpf_tuple'])) <= 0,
            'Reversed directed interval')
    return 1


def records(tree):
    if isinstance(tree, dict):
        if 'lower_exact_mpf_tuple' in tree or 'upper_exact_mpf_tuple' in tree:
            return record(tree)
        return sum(records(value) for value in tree.values())
    if isinstance(tree, list):
        return sum(records(value) for value in tree)
    return 0


def cap(row, order):
    require(row['ordinary_Z_order'] == order and type(row['exact_zero']) is bool, 'Cap order/marker changed')
    value = row['log_absolute_upper_bound']
    if row['exact_zero']:
        require(value is None, 'Exact zero has a cap')
    else:
        record(value)


def below(row, threshold):
    return row['exact_zero'] or mpf_cmp(endpoint(row['log_absolute_upper_bound']['upper_exact_mpf_tuple']),
        endpoint(threshold['lower_exact_mpf_tuple'])) < 0


def source_guard():
    raw = (BASE / C5).read_bytes()
    require(digest(raw) == C5_SHA, 'C5 source changed')
    tree = ast.parse(raw, filename=C5)
    node = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'loop_functions')
    actual = digest(ast.dump(node, include_attributes=False).encode())
    require(actual == C5_AST, 'C5 loop AST changed')
    return dict(path='experiments/root_st073/' + C5, source_sha256=C5_SHA, function='loop_functions',
                AST_sha256=actual, source_executed=False,
                source_ports=['B_over_Pstar_phi', 'direction_at_inverse', 'phase_lambda_at_inverse', 'Phi'],
                post_conditions=['a=4/5', 'b=t0=0', 'q,K constant in Z', 'Pstar,E0 frozen in Z'])


def fixture_packets():
    result = {}
    for name, module in zip(FIXTURES, MODULES):
        raw = (BASE / name).read_bytes()
        receipt = json.loads(raw)
        require(receipt['passed'] is True and receipt['module_sha256'] == REVIEWED[module],
                'Fixture failed or source unbound: ' + name)
        require(receipt.get('checks') and all(value is True for value in receipt['checks'].values()),
                'Fresh check missing or failed: ' + name)
        result[name] = dict(bytes=len(raw), sha256=digest(raw), receipt=receipt)
    identity = result[FIXTURES[0]]['receipt']
    require(identity['checks']['exact_chainrule_contraction'] is True
            and identity['exact_algebra_proof']['residual'] == '0'
            and identity['ratio_certified'] is False
            and identity['source_domain_premise']['actual_signed_native_parameters_bound'] is False,
            'Exact/generic source proof scope changed')
    proof = result[FIXTURES[2]]['receipt']['proof']
    require(all(proof[key] == '0' for key in ('I_psi_residual_mod_unit_circle', 'I_r_residual',
                't_u_residual_mod_unit_circle', 't_psi_residual')), 'Closed primitive exact proof changed')
    return result


def pin_maps(old, fixtures):
    inputs = dict(old['input_hashes'])
    inputs[PARENT] = PARENT_SHA
    inputs.update({name: value['sha256'] for name, value in fixtures.items()})
    implementations = dict(old['implementation_hashes'])
    implementations.update({name: sha(ROOT / name) for name in MODULES + (RUNNER,)})
    return inputs, implementations


def produce(raw, old):
    require(not CHECKPOINT.exists() and not REPORT.exists(), 'Saved new math exists; replay refused')
    sys.path.insert(0, str(ROOT / 'src'))
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_post_restore_fast_phi_bound as provider
    c = MPIntervalContext()
    c.prec = 1600
    started = time.monotonic()
    fast = provider.compute(c, raw)
    checkpoint = dict(fast=fast, parent_report_sha256=PARENT_SHA,
        math_source_hashes={name: sha(ROOT / name) for name in MODULES},
        computation_runner_sha256=sha(ROOT / RUNNER), interval_precision_bits=c.prec,
        accepted_ancestors_replayed=False, execution_seconds=time.monotonic() - started)
    CHECKPOINT.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                          encoding='utf-8', newline='\n')
    started = time.monotonic()
    checkpoint['endpoint_C1'] = provider.endpoint_C1(c, old, fast['new_fast_Z1'])
    checkpoint['execution_seconds'] += time.monotonic() - started
    CHECKPOINT.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                          encoding='utf-8', newline='\n')


def validate(old, report):
    fixtures = fixture_packets()
    inputs, implementations = pin_maps(old, fixtures)
    require(report['input_hashes'] == inputs and report['implementation_hashes'] == implementations, 'Pin maps changed')
    require(report['fresh_fixtures'] == fixtures and report['source_guard'] == source_guard(), 'Proof receipt changed')
    require(all(report['scope'][key] is False for key in FALSE_SCOPE), 'Downstream scope overclaimed')
    packets = old['actual_directed_Rm_controls']
    rm, loop = packets['Rm_endpoint'], packets['loop_mixed']
    fast, new = report['actual_fast_phi_correlated_C1']['fast'], report['actual_fast_phi_correlated_C1']['endpoint_C1']
    require(fast['passed'] is True and new['passed'] is True, 'New component failed')
    require(fast['input_report_sha256'] == PARENT_SHA and fast['chart'] == new['chart'] == 'post_restore',
            'Source/chart binding changed')
    require(fast['old_B_phi_over_Pstar_over_frozen_E0'] == loop['B_phi_over_Pstar_over_frozen_E0'], 'Old fast caps changed')
    require(fast['old_value_row_0_0'] == loop['B_phi_over_Pstar_over_frozen_E0'][0], 'Value row changed')
    for index, row in enumerate(fast['new_B_phi_over_Pstar_over_frozen_E0']):
        if index != 1:
            require(row == loop['B_phi_over_Pstar_over_frozen_E0'][index], 'Unrelated mixed row changed')
    cap(fast['new_fast_Z1'], 1)
    cap(fast['physical_deltaUy_Z1'], 1)
    for key in ('new_deltaUy_Z1', 'new_deltaQy_Z1', 'new_deltaQy_Z1_over_Rm', 'new_combined_Qy_over_Rm_Z1'):
        cap(new[key], 1)
    require(new['old_deltaUy_Z0_preserved'] == rm['finite_deltaUy']['log_absolute_upper_bounds'][0], 'deltaUy0 replayed/changed')
    require(new['old_deltaMy_Z0_Z2_preserved'] == rm['finite_deltaMy']['log_absolute_upper_bounds'][:3], 'Memory reset/changed')
    require(new['old_deltaQy_Z0_preserved'] == rm['finite_deltaQy']['log_absolute_upper_bounds'][0], 'deltaQy0 changed')
    require(new['native_Qy_over_Rm_signed_log_scale'] == rm['normalization']['native_Qy_over_Rm_signed_log_scale'],
            'Signed native endpoint changed')
    require(new['external_pressure_boundary'] == rm['pressure_boundary']
            and new['external_pressure_boundary']['coefficient'] == -1, 'External endpoint removed')
    require(new['caps_used_as_signed_values'] is False and new['old_caps_retained'] is True, 'Cap/value or error contract changed')
    threshold = rm['normalization']['threshold']
    expected = below(rm['normalization']['combined_Qy_over_Rm_log_absolute_bounds'][0], threshold)
    expected = expected and below(new['new_combined_Qy_over_Rm_Z1'], threshold)
    require(report['scope']['Rm_normalized_C1_endpoint_target_met'] is expected, 'Threshold claim differs from exact tuples')
    enclosure = report['full_normalized_endpoint_C1_enclosure']
    if expected:
        require(enclosure['normalized_by'] == 'Wpost=Rm^2'
                and enclosure['enclosure_not_actual_zero_or_selected_value'] is True, 'Full endpoint scope changed')
        for k, row in enumerate(enclosure['rows']):
            require(row['ordinary_Z_order'] == k, 'Full endpoint row order changed')
            record(row['value'])
            require(endpoint(row['value']['lower_exact_mpf_tuple']) == (1, 1, -1024, 1)
                    and endpoint(row['value']['upper_exact_mpf_tuple']) == (0, 1, -1024, 1), 'Wrong endpoint radius')
    else:
        require(enclosure is None, 'Unsupported small endpoint enclosure')
    require(report['execution_audit']['accepted_ancestors_replayed'] is False, 'Old mathematics replay')
    return records(report['actual_fast_phi_correlated_C1'])


def finalize(old):
    require(not REPORT.exists() and CHECKPOINT.exists(), 'No checkpoint or saved report exists')
    raw = CHECKPOINT.read_bytes()
    cached = json.loads(raw)
    require(cached['math_source_hashes'] == REVIEWED and cached['parent_report_sha256'] == PARENT_SHA,
            'Cached math/source binding changed')
    require('endpoint_C1' in cached, 'New endpoint component not complete; do not repeat completed fast component')
    fixtures = fixture_packets()
    inputs, implementations = pin_maps(old, fixtures)
    norm = old['actual_directed_Rm_controls']['Rm_endpoint']['normalization']
    passed = below(norm['combined_Qy_over_Rm_log_absolute_bounds'][0], norm['threshold'])
    passed = passed and below(cached['endpoint_C1']['new_combined_Qy_over_Rm_Z1'], norm['threshold'])
    enclosure = None
    if passed:
        value = dict(lower_exact_mpf_tuple=[1, 1, -1024, 1], upper_exact_mpf_tuple=[0, 1, -1024, 1],
            lower_display='-2^-1024', upper_display='2^-1024', directed=True)
        enclosure = dict(normalized_by='Wpost=Rm^2', source_quantity='Rm*Qy/Wpost=Qy/Rm',
            rows=[dict(ordinary_Z_order=k, value=value) for k in (0, 1)],
            enclosure_not_actual_zero_or_selected_value=True, external_pressure_coefficient=-1,
            native_source_and_finite_error_kept_separate=True)
    report = dict(stage='current_original_fast_phi_correlated_C1', parent_report=PARENT,
        parent_report_sha256=PARENT_SHA, parent_report_bytes=PARENT_BYTES,
        input_hashes=inputs, implementation_hashes=implementations, source_guard=source_guard(),
        fresh_fixtures=fixtures, actual_fast_phi_correlated_C1=dict(fast=cached['fast'], endpoint_C1=cached['endpoint_C1']),
        full_normalized_endpoint_C1_enclosure=enclosure,
        scope=dict(new_signed_source_derivative_interface=True, new_finite_input_closed_phase_function=True,
            new_source_correlated_fast_Z1_bound=True, Rm_normalized_C1_endpoint_target_met=passed,
            inherited_cap_retained_as_premise=True, full_inherited_cap_revalidation=False,
            **{name: False for name in FALSE_SCOPE}),
        execution_audit=dict(accepted_ancestors_replayed=False, only_new_correlated_row_and_operator_computed=True,
            checkpoint_sha256=digest(raw), checkpoint_bytes=len(raw),
            computation_runner_sha256=cached['computation_runner_sha256'], actual_signed_native_inputs_available=False),
        interval_precision_bits=cached['interval_precision_bits'], execution_seconds=cached['execution_seconds'])
    validate(old, report)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                      encoding='utf-8', newline='\n')
    return report


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write-report', action='store_true')
    mode.add_argument('--finalize-checkpoint', action='store_true')
    mode.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    preflight()
    fixture_packets()
    raw, old = parent()
    if args.write_report:
        produce(raw, old)
        report = finalize(old)
    elif args.finalize_checkpoint:
        report = finalize(old)
    else:
        report = json.loads(REPORT.read_bytes())
    count = validate(old, report)
    print(json.dumps(dict(passed=True, report=str(REPORT), sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        canonical_records=count, input_pins=len(report['input_hashes']), implementation_pins=len(report['implementation_hashes']),
        execution_seconds=report['execution_seconds'], scope=report['scope'])))


if __name__ == '__main__':
    main()
