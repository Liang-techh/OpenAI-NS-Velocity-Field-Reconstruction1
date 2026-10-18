"""Evaluate new signed native sources once; inspect saved bytes without math."""
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
BASE = ROOT/'experiments/root_st073'
RUNNER = 'scripts/run_original_signed_native_sources.py'
PARENT = 'lei_ren_part1_paper_compliant_current_original_fast_phi_correlated_C1.json'
PARENT_SHA = '8b77d380fcfd4837c4bcfff10a1e4fe512b37d47b30340c5093552bb3b1ca108'
NATIVE = 'lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json'
NATIVE_SHA = '10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a'
NAME = 'lei_ren_part1_paper_compliant_current_original_signed_native_sources.json'
REPORT = BASE/NAME
CHECKPOINT_DIR = ROOT/'tmp/original_signed_native_sources'
MODULES = tuple('src/openai_ns_reconstruction/'+name for name in (
    'original_signed_native_history_kernel.py', 'original_restoration_gamma2.py',
    'original_signed_native_p2_kernel.py', 'original_signed_native_history_data.py'))
REVIEWED = dict(zip(MODULES, (
    'd81af96ac03a486b28bc0b0a6f106fb5b7e19af4324d9ee5928db35bf30a225e',
    'c6663728251083fd41f47f45bffa486d91b3a0231e9993799d3f139c64d3af1c',
    '17b3fda3a8e4bcdf732557cc5ea391deff98ad3d456e5707abe2504c823b5a84',
    '5ce14c5ada9f75a36f4c7bf4fa3020b69f04bd695a305e7043f063d22a58a323')))
FIXTURES = ('original_signed_native_history_fixture.json', 'original_restoration_gamma2_fixture.json',
            'original_signed_native_p2_fixture.json')
FALSE_SCOPE = ('signed_pressure_join_completed', 'independent_reserved_full_pulse_m15_matched',
    'global_five_moment_closure', 'exact_heat_exterior', 'admissible_stress_and_flat_remainder',
    'true_n_dependent_recursion', 'realizable_oscillatory_correction',
    'independent_Cartesian_NS_validation', 'legacy_packets_recertified',
    'full_phase_or_inverse_evaluated', 'full_goal_complete')


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def sha(path):
    return digest(Path(path).read_bytes())


def write(path, value):
    require(not path.exists(), 'Saved new calculation exists; replay/overwrite refused: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def checkpoint(name):
    return CHECKPOINT_DIR/(name+'.json')


def input_path(name):
    p = Path(name)
    require(not p.is_absolute(), 'Absolute input path')
    if '..' in p.parts:
        require(name == '../../src/openai_ns_reconstruction/schedule_pressure.py', 'Input traversal')
        return (BASE/p).resolve()
    return ROOT/p if p.parts[0] == 'experiments' else BASE/p.name


def parent():
    raw = (BASE/PARENT).read_bytes()
    require(len(raw) == 497555 and digest(raw) == PARENT_SHA, 'Parent changed before decode')
    old = json.loads(raw)
    for name, expected in old['input_hashes'].items():
        require(sha(input_path(name)) == expected, 'Input changed: '+name)
    for name, expected in old['implementation_hashes'].items():
        require(sha(ROOT/name) == expected, 'Implementation changed: '+name)
    return old


def preflight():
    for name in MODULES+(RUNNER,):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF: '+name)
        ast.parse(raw, filename=name)
    for name, expected in REVIEWED.items():
        require(sha(ROOT/name) == expected, 'Reviewed source changed: '+name)


def fixtures():
    result = {}
    for name, module in zip(FIXTURES, MODULES[:3]):
        raw = (BASE/name).read_bytes()
        value = json.loads(raw)
        require(value['passed'] is True and value['module_sha256'] == REVIEWED[module], 'Fixture/source mismatch: '+name)
        require(all(value.get('checks', {'passed': True}).values()), 'Fresh fixture failed: '+name)
        result[name] = dict(bytes=len(raw), sha256=digest(raw), receipt=value)
    return result


def source_guard():
    prefix = 'lei_ren_part1_paper_compliant_current_original_'
    functions = ((prefix+'prescribed_R100_Rm_C6_functions.py', 'phi_from_Utheta'),
                 (prefix+'prescribed_R100_Rm_C6_functions.py', '__init__'),
                 (prefix+'prescribed_bridge_C6_functions.py', 'generic_fields'),
                 (prefix+'core_inlet_axial6_generic_source.py', 'generic_functions'))
    result = []
    for name, function in functions:
        raw = (BASE/name).read_bytes()
        nodes = [v for v in ast.walk(ast.parse(raw)) if isinstance(v, ast.FunctionDef) and v.name == function]
        require(len(nodes) == 1, 'Source function is ambiguous: '+name+':'+function)
        result.append(dict(path='experiments/root_st073/'+name, function=function,
                           source_sha256=digest(raw), AST_sha256=digest(ast.dump(nodes[0], include_attributes=False).encode()),
                           source_executed=False))
    return result


def produce(component):
    path = checkpoint(component)
    require(not path.exists() and not REPORT.exists(), 'Accepted/new saved math exists; replay refused')
    preflight()
    parent()
    fixtures()
    sys.path.insert(0, str(ROOT/'src'))
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_signed_native_history_data as provider
    native_raw = (BASE/NATIVE).read_bytes()
    require(len(native_raw) == 42177499 and digest(native_raw) == NATIVE_SHA, 'Native packet changed')
    c = MPIntervalContext()
    c.prec = 1600
    started = time.monotonic()
    if component == 'gamma2':
        from openai_ns_reconstruction import original_restoration_gamma2 as gamma
        value = gamma.validated_gamma2(tolerance='1e-8', dps=100, root=ROOT).as_dict()
        value['total_directed_record'] = dict(lower_exact_mpf_tuple=value['lower_mpf_tuple'],
            upper_exact_mpf_tuple=value['upper_mpf_tuple'], directed=True)
    elif component == 'histories':
        gamma = json.loads(checkpoint('gamma2').read_bytes())
        value = provider.compute(c, native_raw, gamma['payload']['total_directed_record'])
    else:
        histories = json.loads(checkpoint('histories').read_bytes())
        value = provider.compute_p2(c, native_raw, histories['payload'])
    # Persist each new calculation before any report inspection/packaging.
    result = dict(component=component, payload=value, execution_seconds=time.monotonic()-started,
        interval_precision_bits=c.prec if component != 'gamma2' else None,
        math_source_hashes=REVIEWED, computation_runner_sha256=sha(ROOT/RUNNER),
        native_report_sha256=NATIVE_SHA, accepted_ancestors_replayed=False)
    write(path, result)
    print(json.dumps(dict(component=component, checkpoint=str(path), sha256=sha(path), seconds=result['execution_seconds'])))


def endpoint(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4 and all(type(v) is int for v in value), 'Malformed MPF')
    s, m, e, b = value
    require(s in (0, 1) and m >= 0 and b == m.bit_length(), 'Nonfinite MPF')
    require((m == 0 and tuple(value) == (0, 0, 0, 0)) or (m > 0 and m & 1), 'Noncanonical MPF')
    return tuple(value)


def records(tree):
    count = 0
    if isinstance(tree, dict):
        for lower, upper in (('lower_exact_mpf_tuple', 'upper_exact_mpf_tuple'), ('lower_mpf_tuple', 'upper_mpf_tuple')):
            if lower in tree and upper in tree:
                require(mpf_cmp(endpoint(tree[lower]), endpoint(tree[upper])) <= 0, 'Reversed interval')
                count += 1
        count += sum(records(v) for v in tree.values() if isinstance(v, (dict, list)))
    elif isinstance(tree, list):
        count += sum(records(v) for v in tree if isinstance(v, (dict, list)))
    return count


def package():
    require(not REPORT.exists(), 'Saved report exists; overwrite refused')
    preflight()
    old, fx = parent(), fixtures()
    calculations = {name: json.loads(checkpoint(name).read_bytes()) for name in ('gamma2', 'histories', 'p2')}
    for value in calculations.values():
        require(value['math_source_hashes'] == REVIEWED and not value['accepted_ancestors_replayed'], 'Calculation/source mismatch')
    inputs, implementations = dict(old['input_hashes']), dict(old['implementation_hashes'])
    inputs[PARENT], inputs[NATIVE] = PARENT_SHA, NATIVE_SHA
    inputs.update({name: value['sha256'] for name, value in fx.items()})
    implementations.update({name: sha(ROOT/name) for name in MODULES+(RUNNER,)})
    report = dict(schema='original_signed_native_sources_v1', parent_report=PARENT, parent_report_sha256=PARENT_SHA,
        native_report=NATIVE, native_report_sha256=NATIVE_SHA, input_hashes=inputs, implementation_hashes=implementations,
        fresh_fixtures=fx, source_guard=source_guard(), actual_signed_native_sources=calculations,
        execution_seconds=sum(v['execution_seconds'] for v in calculations.values()),
        execution_audit=dict(accepted_ancestors_replayed=False, each_new_component_evaluated_once=True,
            checkpoint_hashes={name: sha(checkpoint(name)) for name in calculations}),
        scope={**dict.fromkeys(FALSE_SCOPE, False), 'actual_signed_scaled_native_p2_C1_installed': True,
               'actual_signed_post_restore_A_b_c_Z0_Z2_installed': True,
               'Rm_normalized_C1_endpoint_target_met': old['scope']['Rm_normalized_C1_endpoint_target_met']},
        preserved=dict(cumulative_terms=212, signed_core_lower_terms=6, reserved_n1_own_swirl_registration_count=1,
                       all_inherited_errors=True, external_Rm_handles=[183696, 183697], external_coefficient=-1,
                       P1_is_not_m15=True, pressure_Q_is_not_loop_root_q=True),
        limits=['Signed interval enclosures of source functions; no selected midpoint or cap value.',
                'Old whole-chart bounds remain inherited premises; no blanket recertification.',
                'Physical p2=positive R*Utheta/2 times signed K; physical radius not materialized.',
                'No phase/inverse evaluation or pressure endpoint join performed in this packet.'])
    write(REPORT, report)
    print(json.dumps(dict(report=str(REPORT), bytes=REPORT.stat().st_size, sha256=sha(REPORT))))


def inspect():
    preflight()
    old, fx = parent(), fixtures()
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    expected_inputs = {**old['input_hashes'], PARENT: PARENT_SHA, NATIVE: NATIVE_SHA,
                       **{name: value['sha256'] for name, value in fx.items()}}
    expected_implementations = {**old['implementation_hashes'], **{name: sha(ROOT/name) for name in MODULES+(RUNNER,)}}
    require(report['input_hashes'] == expected_inputs and report['implementation_hashes'] == expected_implementations, 'Source pins changed')
    require(report['fresh_fixtures'] == fx and report['source_guard'] == source_guard(), 'Fixture/AST source changed')
    require(all(report['scope'][k] is False for k in FALSE_SCOPE), 'Downstream scope overclaimed')
    packet = report['actual_signed_native_sources']
    gamma, histories, p2 = (packet[k]['payload'] for k in ('gamma2', 'histories', 'p2'))
    require(mpf_cmp(endpoint(gamma['width_upper_mpf_tuple']), endpoint(gamma['tolerance_mpf_tuple'])) <= 0, 'Gamma2 width exceeds tolerance')
    require(histories['passed'] and histories['incoming_A_B_C_retained'] and not histories['caps_promoted_to_values'], 'History scope failed')
    require(p2['passed'] and p2['full_native_p2_signed_scaled_enclosure'] and not p2['midpoint_or_logcap_promoted'], 'Signed p2 scope failed')
    for name in ('post_restore', 'Rm'):
        chart = p2['charts'][name]
        require(chart['p2_value_sign'] in ('positive', 'negative') and len(chart['K']) == 2, 'Signed p2 data unavailable')
        for order, row in enumerate(chart['K']):
            require(row['ordinary_Z_order'] == order, 'Ordinary-Z order changed')
        value = chart['K'][0]['value']
        zero = (0, 0, 0, 0)
        require(mpf_cmp(endpoint(value['upper_exact_mpf_tuple']), zero) < 0 if chart['p2_value_sign'] == 'negative' else
                mpf_cmp(endpoint(value['lower_exact_mpf_tuple']), zero) > 0, 'p2 sign claim is unsupported')
    receipt = dict(passed=True, report=NAME, sha256=digest(raw), bytes=len(raw), canonical_records=records(packet),
        input_pins=len(expected_inputs), implementation_pins=len(expected_implementations),
        inspection_math_providers_imported=False, accepted_ancestors_replayed=False)
    print(json.dumps(receipt))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--compute', choices=('gamma2', 'histories', 'p2'))
    actions.add_argument('--write-report', action='store_true')
    actions.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    if args.compute:
        produce(args.compute)
    elif args.write_report:
        package()
    else:
        inspect()
