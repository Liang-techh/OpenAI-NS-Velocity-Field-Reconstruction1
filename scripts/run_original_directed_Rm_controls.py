"""New post-restore mixed/Rm calculation; inspection reads saved data only."""
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
RUNNER = 'scripts/run_original_directed_Rm_controls.py'
MODULES = (
    'src/openai_ns_reconstruction/original_directed_post_restore_mixed_native.py',
    'src/openai_ns_reconstruction/original_directed_post_restore_mixed_loop.py',
    'src/openai_ns_reconstruction/original_directed_Rm_endpoint_kernel.py',
)
PARENTS = {
    'native': ('lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json',
               42177499, '10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a'),
    'loop': ('lei_ren_part1_paper_compliant_current_original_directed_loop_controls.json',
             1752702, 'adefdae708754a2572eae15a90662274bf17e88a0d72d5db9befbc9e2858b252'),
    'finite_N': ('lei_ren_part1_paper_compliant_current_original_directed_finite_N_controls.json',
                 1368720, '152d211cd137e7e0498fd157a0ac1467541dc6214ae59c13af62844112003d9a'),
}
NAME = 'lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json'
REPORT = BASE / NAME
CHECKPOINT = BASE / (NAME + '.new_math_checkpoint')
FIXTURES = (
    'original_directed_post_restore_mixed_native_fixture.json',
    'original_directed_post_restore_mixed_loop_fixture.json',
    'original_directed_Rm_endpoint_fixture.json',
)
REVIEWED = {
    MODULES[0]: 'cd524586754ce41a70fcb9fb7ac4434a8b01a75dc110d68712aecce636797c95',
    MODULES[1]: '770c81fc856e60dfaafbba48c06c3dbd002c59c2d5e7da2d1246ac7ffc6b92c3',
    MODULES[2]: '6bf8bfa91a8c3a54a12bf43199fc5f52b37662f136ef83b73118110b6000dbb9',
}
FALSE_SCOPE = (
    'global_five_moment_closure', 'signed_pressure_join_completed',
    'independent_reserved_full_pulse_m15_matched', 'exact_heat_exterior',
    'admissible_stress_and_flat_remainder', 'true_n_dependent_recursion',
    'realizable_oscillatory_correction', 'independent_Cartesian_NS_validation',
    'legacy_packets_recertified', 'inherited_global_guard_admitted', 'full_goal_complete',
)
MIXED = [(i, k) for i in (0, 1) for k in range(4)]


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
                'New Python must use LF and exactly one EOF newline: ' + name)
        ast.parse(raw, filename=name)
    for name, expected in REVIEWED.items():
        require(expected != 'REVIEW_PENDING' and sha(ROOT / name) == expected,
                'Review pending or source changed: ' + name)


def parents():
    raw, parsed = {}, {}
    for key, (name, size, expected) in PARENTS.items():
        value = (BASE / name).read_bytes()
        require(len(value) == size and digest(value) == expected, 'Parent changed before decode: ' + key)
        raw[key], parsed[key] = value, json.loads(value)
    finite = parsed['finite_N']
    require(finite['parent_report_sha256'] == PARENTS['loop'][2]
            and finite['native_report_sha256'] == PARENTS['native'][2], 'Parent chain changed')
    require(finite['scope']['current_N_I1_and_B_target_met'] is True, 'Finite-N premise missing')
    for name, expected in finite['input_hashes'].items():
        require(sha(input_path(name)) == expected, 'Input changed: ' + name)
    for name, expected in finite['implementation_hashes'].items():
        path = Path(name)
        require(not path.is_absolute() and '..' not in path.parts, 'Unsafe implementation path')
        require(sha(ROOT / path) == expected, 'Implementation changed: ' + name)
    return raw, parsed


def endpoint(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4
            and all(type(x) is int for x in value), 'Malformed endpoint tuple')
    sign, mantissa, exponent, bits = value
    require(sign in (0, 1) and mantissa >= 0 and bits == mantissa.bit_length(), 'Nonfinite endpoint')
    require((mantissa == 0 and tuple(value) == (0, 0, 0, 0))
            or (mantissa > 0 and mantissa & 1), 'Noncanonical endpoint')
    return tuple(value)


def record(value):
    require(isinstance(value, dict) and value.get('directed') is True, 'Directed record missing')
    require(mpf_cmp(endpoint(value['lower_exact_mpf_tuple']),
                    endpoint(value['upper_exact_mpf_tuple'])) <= 0, 'Reversed record')
    return 1


def records(tree):
    if isinstance(tree, dict):
        if 'lower_exact_mpf_tuple' in tree or 'upper_exact_mpf_tuple' in tree:
            return record(tree)
        return sum(records(value) for value in tree.values())
    if isinstance(tree, list):
        return sum(records(value) for value in tree)
    return 0


def cap_row(row):
    require(type(row.get('exact_zero')) is bool, 'Exact zero marker missing')
    value = row.get('log_absolute_upper_bound')
    if row['exact_zero']:
        require(value is None, 'Exact zero carries a bound')
    else:
        record(value)


def mixed(rows):
    require(isinstance(rows, list) and len(rows) == 8, 'Eight mixed rows required')
    for index, row in zip(MIXED, rows):
        require((row['physical_y_order'], row['ordinary_Z_order']) == index, 'Mixed row order changed')
        cap_row(row)


def same_cap(left, right):
    if left is None or right is None:
        return left is right
    return all(endpoint(left[key]) == endpoint(right[key]) for key in
               ('lower_exact_mpf_tuple', 'upper_exact_mpf_tuple'))


def z_rows(rows, count):
    require(isinstance(rows, list) and len(rows) == count, 'Z row count changed')
    for k, row in enumerate(rows):
        require(row['ordinary_Z_order'] == k, 'Z order changed')
        if 'value' in row:
            record(row['value'])
        else:
            cap_row(row)


def source_guards():
    names = {
        'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py': ('run',),
        'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py': ('densities', 'histories'),
        'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py': ('loop_functions',),
        'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': ('total_radial_functions', 'reduced_pressure_S'),
    }
    result = {}
    for name, selected in names.items():
        path = BASE / name
        raw = path.read_bytes()
        tree = ast.parse(raw, filename=name)
        functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
        require(all(item in functions for item in selected), 'Source function changed: ' + name)
        result[name] = dict(source_sha256=digest(raw), source_executed=False,
            function_AST_sha256={item: digest(ast.dump(functions[item], include_attributes=False).encode())
                                 for item in selected})
    return result


def validate_packets(parent, packets):
    native, loop, rm = (packets[name] for name in ('native_mixed', 'loop_mixed', 'Rm_endpoint'))
    require(native['passed'] is True and loop['passed'] is True and rm['passed'] is True, 'New packet failed')
    require(native['source_native_report_sha256'] == PARENTS['native'][2], 'Native source binding')
    require(native['chart'] == loop['chart'] == 'post_restore', 'Wrong mixed chart')
    for packet in (native, loop, rm):
        require(packet['accepted_ancestors_replayed'] is False, 'Ancestor mathematics replay')
    saved = parent['native']['actual_native_directed_controls_v2']['raw_interval_tree']['restoration']['charts']['post_restore']
    require(native['conditioning'] == saved['conditioning'], 'Saved positive conditioning changed')
    require(native['input_record_provenance']['native_caps'] == saved['native_caps']
            and native['input_record_provenance']['root_caps'] == saved['root_caps'],
            'Exact saved source receipts changed')
    for key, rows in native['root_mixed_caps'].items():
        mixed(rows)
        if key in ('p2', 'q', 'u', 'h', 'r', 'alpha', 'K'):
            require(all(same_cap(left['log_absolute_upper_bound'], right['log_absolute_upper_bound'])
                        for left, right in zip(rows[:4], saved['root_caps'][key])),
                    'Accepted y0 root rows changed: ' + key)
    for key, rows in native['native_caps'].items():
        mixed(rows)
        require(all(same_cap(left['log_absolute_upper_bound'], right['log_absolute_upper_bound'])
                    for left, right in zip(rows[:4], saved['native_caps'][key][:4])),
                'Accepted native y0 rows changed: ' + key)
    require(all(same_cap(left['log_absolute_upper_bound'], right['log_absolute_upper_bound'])
                for left, right in zip(native['root_mixed_caps']['E_over_frozen_E0'][:4], saved['E_relative'][:4])),
            'Frozen-E normalization y0 rows changed')
    for key in ('A', 'B_over_Pstar_over_frozen_E0', 'B_phi_over_Pstar_over_frozen_E0', 'inverse_Psi'):
        mixed(loop[key])
    phase_indices = [(i, k, j) for total in range(5) for i in (0, 1)
                     for k in range(4) for j in range(5) if i + k + j == total]
    require(loop['phase_intermediate_row_count'] == len(phase_indices) == 24,
            'Mixed phase grid incomplete')
    for rows in loop['phase_intermediates'].values():
        require(len(rows) == 24, 'Phase row count changed')
        for index, row in zip(phase_indices, rows):
            require((row['physical_y_order'], row['ordinary_Z_order'], row['free_angle_order']) == index,
                    'Free-angle mixed order changed')
            cap_row(row)
    require(loop['fast_phi_separate'] is True and loop['B_phi_has_no_N_factor'] is True, 'Slow/fast source mixed')
    require(loop['B_phi_value_control']['derivative_bound_not_inferred_from_value_bound'] is True,
            'Value cap was differentiated')
    require(loop['shared_y0_algebra_recalculated_for_new_grid'] is True
            and loop['old_loop_builder_called'] is False, 'New-grid dependency receipt incorrect')
    require(native['global_five_moment_closure'] is False and loop['full_goal_complete'] is False,
            'Downstream scope overclaimed')
    corrected = parent['native']['actual_restoration_post_jet_repair']['corrected_background']
    require(rm['native_endpoint']['source_rows']
            == corrected['charts']['post_restore']['endpoint_fields']['Qy_end'][:2]
            and rm['native_endpoint']['endpoint_recomputed'] is False, 'Signed native endpoint replaced')
    z_rows(rm['native_endpoint']['signed_Qy'], 2)
    for actual, original in zip(rm['native_endpoint']['signed_Qy'], rm['native_endpoint']['source_rows']):
        require(same_cap(actual['value'], original['value']), 'Signed native endpoints changed')
    finite_saved = parent['finite_N']['actual_directed_finite_N_controls']['charts']['post_restore']['physical_error_log_bounds']
    for key, count in (('delta_U', 4), ('delta_M', 4), ('delta_Q', 3)):
        z_rows(rm['saved_finite_error_bounds'][key], count)
        require(all(same_cap(a['log_absolute_upper_bound'], b['log_absolute_upper_bound'])
                    for a, b in zip(rm['saved_finite_error_bounds'][key], finite_saved[key])),
                'Inherited finite data changed: ' + key)
    for owner in ('finite_deltaUy', 'finite_deltaMy'):
        z_rows(rm[owner]['log_absolute_upper_bounds'], 4)
    z_rows(rm['finite_deltaQy']['log_absolute_upper_bounds'], 2)
    require(all(rm['finite_deltaQy']['operator_receipt'][key] is True for key in
                ('positive_L_lower_used', 'q1_denominator_derivative_correction_included',
                 'delta_positive_and_less_than_one_checked')),
            'Pressure reciprocal omitted lower floor')
    norm = rm['normalization']
    physical_scale = parent['native']['actual_restoration_background_I1_numeric']['physical_scale']['log_Rm']
    require(same_cap(norm['log_Rm'], physical_scale)
            and same_cap(norm['log_Rm'], corrected['charts']['post_restore']['log_R_upper'])
            and norm['Wpost'] == 'Rm^2' and norm['RQy_over_Wpost'] == 'Qy/Rm', 'Physical scale changed')
    z_rows(norm['finite_delta_Qy_over_Rm_log_absolute_bounds'], 2)
    z_rows(norm['combined_Qy_over_Rm_log_absolute_bounds'], 2)
    threshold = norm['threshold']
    record(threshold)
    passed = all(row['exact_zero'] or mpf_cmp(
        endpoint(row['log_absolute_upper_bound']['upper_exact_mpf_tuple']),
        endpoint(threshold['lower_exact_mpf_tuple'])) < 0
        for row in norm['combined_Qy_over_Rm_log_absolute_bounds'])
    require(norm['threshold_passed'] is passed
            and rm['target_status']['normalized_2_power_minus_1024_target_passed'] is passed,
            'Threshold claim does not match exact tuples')
    if passed:
        box = norm['honest_symmetric_enclosure']
        record(box)
        require(endpoint(box['lower_exact_mpf_tuple']) == (1, 1, -1024, 1)
                and endpoint(box['upper_exact_mpf_tuple']) == (0, 1, -1024, 1), 'Wrong enclosure radius')
    else:
        require(norm['honest_symmetric_enclosure'] is None, 'Unsupported small enclosure')
    boundary = rm['pressure_boundary']
    require(boundary['coefficient'] == -1 and boundary['external_source_handles'] == [183696, 183697]
            and boundary['signed_coefficients_by_Z'] == {'0': -1, '1': -1}, 'External endpoint removed')
    require(all(rm['contracts'][key] is False for key in
                ('native_Qy_recomputed', 'caps_used_as_signed_values', 'symmetric_interval_used_as_zero',
                 'N_materialized', 'Rm_materialized')), 'Endpoint contract overclaimed')
    return records(packets)


def fixtures():
    result = {}
    for name, module in zip(FIXTURES, MODULES):
        raw = (BASE / name).read_bytes()
        value = json.loads(raw)
        require(value['passed'] is True and value['module_sha256'] == REVIEWED[module],
                'Fixture failed or source unbound: ' + name)
        require(value.get('checks') and all(item is True for item in value['checks'].values()),
                'Fixture check missing or failed: ' + name)
        result[name] = dict(sha256=digest(raw), bytes=len(raw), receipt=value)
    return result


def pin_maps(parent, fixture_packets):
    inputs = dict(parent['finite_N']['input_hashes'])
    inputs.update({name: expected for name, size, expected in PARENTS.values()})
    inputs.update({name: value['sha256'] for name, value in fixture_packets.items()})
    implementations = dict(parent['finite_N']['implementation_hashes'])
    implementations.update({name: sha(ROOT / name) for name in MODULES + (RUNNER,)})
    return inputs, implementations


def produce(raw, resume=False):
    require(not REPORT.exists(), 'Saved report exists; replay refused')
    if resume:
        require(CHECKPOINT.exists(), 'No new-math checkpoint to resume')
        checkpoint = json.loads(CHECKPOINT.read_bytes())
        require(checkpoint['parent_hashes'] == {key: item[2] for key, item in PARENTS.items()}, 'Checkpoint parent changed')
        for key, name in zip(('native_mixed', 'loop_mixed', 'Rm_endpoint'), MODULES):
            if key in checkpoint['packets']:
                require(checkpoint['math_source_hashes'][name] == REVIEWED[name], 'Completed component changed')
    else:
        require(not CHECKPOINT.exists(), 'New-math checkpoint exists; replay refused')
        checkpoint = dict(packets={}, math_source_hashes={},
            parent_hashes={key: item[2] for key, item in PARENTS.items()},
            execution_seconds=0, interval_precision_bits=1600,
            computation_runner_sha256=sha(ROOT / RUNNER), accepted_ancestors_replayed=False,
            component_execution_records=[], precheckpoint_failure=dict(
                new_native_attempts=1, completed_packets=0, accepted_ancestors_replayed=False,
                source_sha256='71744f50cabed6b5ad4b7ebb51d36ad215fb5d664cc693c77f2c5a4328f0e642',
                exception='ValueError: Square-root log bounds are reversed',
                resolution='Use the certified upper endpoint for the sqrt upper-log interface; add tiny-u regression.'))
    def save():
        CHECKPOINT.write_text(json.dumps(checkpoint, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                              encoding='utf-8', newline='\n')
    save()
    sys.path.insert(0, str(ROOT / 'src'))
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_directed_post_restore_mixed_native as native_provider
    from openai_ns_reconstruction import original_directed_post_restore_mixed_loop as loop_provider
    from openai_ns_reconstruction import original_directed_Rm_endpoint_kernel as endpoint_provider
    c = MPIntervalContext()
    c.prec = 1600
    calls = (
        ('native_mixed', MODULES[0], lambda: native_provider.compute(c, raw['native'])),
        ('loop_mixed', MODULES[1], lambda: loop_provider.compute(c, checkpoint['packets']['native_mixed'], raw['loop'])),
        ('Rm_endpoint', MODULES[2], lambda: endpoint_provider.compute(c, checkpoint['packets']['native_mixed'],
            checkpoint['packets']['loop_mixed'], raw['finite_N'], raw['native'])),
    )
    for key, name, call in calls:
        if key in checkpoint['packets']:
            continue
        started = time.monotonic()
        try:
            packet = call()
        except Exception as error:
            elapsed = time.monotonic() - started
            checkpoint['component_execution_records'].append(dict(component=key, passed=False,
                source_sha256=sha(ROOT / name), seconds=elapsed, exception=type(error).__name__ + ': ' + str(error)))
            checkpoint['execution_seconds'] += elapsed
            save()
            raise
        elapsed = time.monotonic() - started
        checkpoint['packets'][key] = packet
        checkpoint['math_source_hashes'][name] = sha(ROOT / name)
        checkpoint['component_execution_records'].append(dict(component=key, passed=True,
            source_sha256=sha(ROOT / name), seconds=elapsed))
        checkpoint['execution_seconds'] += elapsed
        save()


def validate_report(parent, report):
    fixture_packets = fixtures()
    inputs, implementations = pin_maps(parent, fixture_packets)
    require(report['input_hashes'] == inputs and report['implementation_hashes'] == implementations, 'Pin maps changed')
    require(report['fresh_fixtures'] == fixture_packets, 'Fixture receipts changed')
    require(report['source_guards'] == source_guards(), 'Source guard changed')
    require(all(report['scope'][name] is False for name in FALSE_SCOPE), 'Downstream scope overclaimed')
    require(report['execution_audit']['accepted_ancestors_replayed'] is False, 'Replay flag')
    require(report['historical_payload_retained_by_parent_hash'] is True, 'Historical preservation missing')
    require(report['scope']['Rm_normalized_C1_endpoint_target_met']
            is report['actual_directed_Rm_controls']['Rm_endpoint']['normalization']['threshold_passed'],
            'Endpoint scope mismatch')
    return validate_packets(parent, report['actual_directed_Rm_controls'])


def finalize(parent):
    require(not REPORT.exists() and CHECKPOINT.exists(), 'No new checkpoint or report already exists')
    raw = CHECKPOINT.read_bytes()
    cached = json.loads(raw)
    require(cached['math_source_hashes'] == REVIEWED, 'Cached math source changed')
    require(cached['parent_hashes'] == {key: item[2] for key, item in PARENTS.items()}, 'Cached input changed')
    fixture_packets = fixtures()
    inputs, implementations = pin_maps(parent, fixture_packets)
    report = dict(stage='current_original_directed_Rm_controls',
        parents={key: dict(report=item[0], bytes=item[1], sha256=item[2]) for key, item in PARENTS.items()},
        input_hashes=inputs, implementation_hashes=implementations,
        source_guards=source_guards(), fresh_fixtures=fixture_packets,
        actual_directed_Rm_controls=cached['packets'], historical_payload_retained_by_parent_hash=True,
        scope=dict(new_post_restore_mixed_controls_computed=True,
                   Rm_normalized_C1_endpoint_target_met=cached['packets']['Rm_endpoint']['normalization']['threshold_passed'],
                   inherited_cap_retained_as_premise=True,
                   full_inherited_cap_revalidation=False, **{name: False for name in FALSE_SCOPE}),
        execution_audit=dict(accepted_ancestors_replayed=False, new_inputs_only=True,
            computation_runner_sha256=cached['computation_runner_sha256'],
            checkpoint_sha256=digest(raw), checkpoint_bytes=len(raw),
            precheckpoint_failure=cached['precheckpoint_failure'],
            component_execution_records=cached['component_execution_records']),
        execution_seconds=cached['execution_seconds'], interval_precision_bits=cached['interval_precision_bits'])
    validate_report(parent, report)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                      encoding='utf-8', newline='\n')
    return report


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write-report', action='store_true')
    mode.add_argument('--finalize-checkpoint', action='store_true')
    mode.add_argument('--resume-new-math', action='store_true')
    mode.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    preflight()
    fixtures()
    raw, parent = parents()
    if args.write_report:
        produce(raw)
        report = finalize(parent)
    elif args.resume_new_math:
        produce(raw, resume=True)
        report = finalize(parent)
    elif args.finalize_checkpoint:
        report = finalize(parent)
    else:
        report = json.loads(REPORT.read_bytes())
    count = validate_report(parent, report)
    print(json.dumps(dict(passed=True, report=str(REPORT), sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        canonical_records=count, input_pins=len(report['input_hashes']),
        implementation_pins=len(report['implementation_hashes']), execution_seconds=report['execution_seconds'],
        scope=report['scope'], Rm_target_met=report['scope']['Rm_normalized_C1_endpoint_target_met'])))


if __name__ == '__main__':
    main()
