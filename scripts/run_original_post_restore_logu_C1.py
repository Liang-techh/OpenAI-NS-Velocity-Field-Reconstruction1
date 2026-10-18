"""New logarithmic-u contraction calculation with saved-only inspection."""
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
RUNNER = 'scripts/run_original_post_restore_logu_C1.py'
PARENT = 'lei_ren_part1_paper_compliant_current_original_signed_native_sources.json'
PARENT_SHA = '9887a353aeab1b24f301a5b6c240f83f58a5d54d001eaac6d16f932b87979561'
NAME = 'lei_ren_part1_paper_compliant_current_original_post_restore_logu_C1.json'
REPORT = BASE/NAME
CHECKPOINT = ROOT/'tmp/original_post_restore_logu_C1/new_math.json'
MODULES = tuple('src/openai_ns_reconstruction/'+name for name in (
    'original_post_restore_logu_contraction.py', 'original_post_restore_logu_endpoint.py'))
REVIEWED = dict(zip(MODULES, (
    '691b0179a603577596dc27a65480406d42063e5c3fb23b9692d714a64af9d2e3',
    'e6252fcae41c982460d7db1a72fd17f437d26c799ab65a0aa88c47475113c465')))
FIXTURES = ('original_post_restore_logu_contraction_fixture.json',
            'original_post_restore_logu_endpoint_fixture.json')
FALSE_SCOPE = ('signed_pressure_join_completed', 'independent_reserved_full_pulse_m15_matched',
    'global_five_moment_closure', 'exact_heat_exterior', 'admissible_stress_and_flat_remainder',
    'true_n_dependent_recursion', 'realizable_oscillatory_correction',
    'independent_Cartesian_NS_validation', 'legacy_packets_recertified',
    'actual_phase_or_inverse_values_installed', 'full_goal_complete')


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def sha(path):
    return digest(Path(path).read_bytes())


def input_path(name):
    p = Path(name)
    require(not p.is_absolute(), 'Absolute input path')
    if '..' in p.parts:
        require(name == '../../src/openai_ns_reconstruction/schedule_pressure.py', 'Input traversal')
        return (BASE/p).resolve()
    return ROOT/p if p.parts[0] == 'experiments' else BASE/p.name


def parent():
    raw = (BASE/PARENT).read_bytes()
    require(len(raw) == 566530 and digest(raw) == PARENT_SHA, 'Parent changed before decode')
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
        require(expected != 'REVIEW_PENDING' and sha(ROOT/name) == expected, 'Reviewed source changed: '+name)


def fixtures():
    result = {}
    for name, module in zip(FIXTURES, MODULES):
        raw = (BASE/name).read_bytes()
        value = json.loads(raw)
        require(value['passed'] is True and value['module_sha256'] == REVIEWED[module], 'Fixture/source mismatch: '+name)
        require(value.get('checks') and all(value['checks'].values()), 'Fresh fixture failed: '+name)
        result[name] = dict(bytes=len(raw), sha256=digest(raw), receipt=value)
    return result


def write(path, value):
    require(not path.exists(), 'Saved new calculation exists; replay/overwrite refused: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def produce():
    require(not CHECKPOINT.exists() and not REPORT.exists(), 'Saved new math exists; replay refused')
    preflight()
    parent()
    fixtures()
    sys.path.insert(0, str(ROOT/'src'))
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_post_restore_logu_endpoint as provider
    c = MPIntervalContext()
    c.prec = 1600
    started = time.monotonic()
    payload = provider.compute(c, *((ROOT/name).read_bytes() for name in (
        provider.SIGNED_REPORT_PATH, provider.FAST_REPORT_PATH,
        provider.RM_REPORT_PATH, provider.NATIVE_REPORT_PATH)))
    endpoint_packet = payload['endpoint_C1']
    if endpoint_packet['threshold_passed']:
        radius = c.mpf(2)**-1024
        symmetric = c.mpf([-radius.b, radius.b])
        payload['honest_normalized_C1_enclosure'] = dict(
            lower_exact_mpf_tuple=list(symmetric._mpi_[0]), upper_exact_mpf_tuple=list(symmetric._mpi_[1]),
            directed=True, ordinary_Z_orders=[0, 1], enclosure_not_actual_zero_or_value=True)
    else:
        payload['honest_normalized_C1_enclosure'] = None
    # Persist the only new actual calculation before packaging or inspection.
    write(CHECKPOINT, dict(payload=payload, execution_seconds=time.monotonic()-started,
        interval_precision_bits=c.prec, math_source_hashes=REVIEWED,
        computation_runner_sha256=sha(ROOT/RUNNER), parent_report_sha256=PARENT_SHA,
        accepted_ancestors_replayed=False))
    print(json.dumps(dict(checkpoint=str(CHECKPOINT), sha256=sha(CHECKPOINT))))


def package():
    require(not REPORT.exists(), 'Report overwrite refused')
    preflight()
    old, fx = parent(), fixtures()
    calc = json.loads(CHECKPOINT.read_bytes())
    require(calc['math_source_hashes'] == REVIEWED and not calc['accepted_ancestors_replayed'], 'Checkpoint/source mismatch')
    inputs, implementations = dict(old['input_hashes']), dict(old['implementation_hashes'])
    inputs[PARENT] = PARENT_SHA
    inputs.update({name: value['sha256'] for name, value in fx.items()})
    implementations.update({name: sha(ROOT/name) for name in MODULES+(RUNNER,)})
    report = dict(schema='original_post_restore_logu_C1_v1', parent_report=PARENT, parent_report_sha256=PARENT_SHA,
        input_hashes=inputs, implementation_hashes=implementations, fresh_fixtures=fx,
        actual_post_restore_logu_C1=calc['payload'], execution_seconds=calc['execution_seconds'],
        interval_precision_bits=calc['interval_precision_bits'],
        execution_audit=dict(accepted_ancestors_replayed=False, new_saved_data_production_attempts=1,
            checkpoint_sha256=sha(CHECKPOINT), computation_runner_sha256=calc['computation_runner_sha256']),
        scope={**dict.fromkeys(FALSE_SCOPE, False), 'actual_signed_scaled_native_p2_C1_installed': True,
               'Rm_normalized_C1_endpoint_target_met': calc['payload']['endpoint_C1']['threshold_passed']},
        preserved=old['preserved'],
        limits=['Uniform analytic fixed-phase contraction bound, not selected phase values.',
                'Saved p2 signed function enclosures enter u_Z/u; absolute caps remain errors.',
                'An endpoint smallness certificate does not prove a terminal functional moment identity.'])
    write(REPORT, report)
    print(json.dumps(dict(report=str(REPORT), bytes=REPORT.stat().st_size, sha256=sha(REPORT))))


def endpoint(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4 and all(type(v) is int for v in value), 'Malformed MPF')
    s, m, e, b = value
    require(s in (0, 1) and m >= 0 and b == m.bit_length(), 'Nonfinite MPF')
    require((m == 0 and tuple(value) == (0, 0, 0, 0)) or (m > 0 and m & 1), 'Noncanonical MPF')
    return tuple(value)


def records(tree):
    if isinstance(tree, dict):
        count = 0
        if 'lower_exact_mpf_tuple' in tree:
            require(mpf_cmp(endpoint(tree['lower_exact_mpf_tuple']), endpoint(tree['upper_exact_mpf_tuple'])) <= 0,
                    'Reversed directed interval')
            count += 1
        return count+sum(records(v) for v in tree.values() if isinstance(v, (dict, list)))
    if isinstance(tree, list):
        return sum(records(v) for v in tree if isinstance(v, (dict, list)))
    return 0


def inspect():
    preflight()
    old, fx = parent(), fixtures()
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    expected_inputs = {**old['input_hashes'], PARENT: PARENT_SHA, **{k:v['sha256'] for k,v in fx.items()}}
    expected_implementations = {**old['implementation_hashes'], **{name:sha(ROOT/name) for name in MODULES+(RUNNER,)}}
    require(report['input_hashes'] == expected_inputs and report['implementation_hashes'] == expected_implementations,
            'Source pins changed')
    require(report['fresh_fixtures'] == fx and report['preserved'] == old['preserved'], 'Saved evidence changed')
    require(all(report['scope'][k] is False for k in FALSE_SCOPE), 'Downstream scope overclaimed')
    packet = report['actual_post_restore_logu_C1']
    require(packet['passed'] is True, 'New adapter failed')
    require(packet['uniform_derivative_bound']['W_upper'] == 600, 'Uniform contraction bound changed')
    proof = fx[FIXTURES[0]]['receipt']
    require(proof['bound']['proved_bound'] == [1061, 2] and proof['bound']['safe_bound'] == [600, 1],
            'Uniform exact proof changed')
    domain = packet['domain_receipt']
    require(domain['certified_abs_u_ge_one'] is True and domain['q_interval'] == ['1/2', '3/2']
            and domain['stored_logs_not_exponentiated'] is True, 'Actual source domain missing')
    require(mpf_cmp(endpoint(domain['log_u_lower']['lower_exact_mpf_tuple']), (0,0,0,0)) > 0,
            'Actual log-u lower bound does not prove domain')
    new = packet['endpoint_C1']
    threshold = endpoint(new['threshold']['lower_exact_mpf_tuple'])
    def below(row):
        return row['exact_zero'] or mpf_cmp(endpoint(row['log_absolute_upper_bound']['upper_exact_mpf_tuple']), threshold) < 0
    target = below(new['both_normalized_rows']['old_Z0']) and below(new['both_normalized_rows']['new_Z1'])
    require(target == new['threshold_passed'] == report['scope']['Rm_normalized_C1_endpoint_target_met'],
            'Endpoint target comparison changed')
    require((packet['honest_normalized_C1_enclosure'] is not None) == target, 'Conditional enclosure changed')
    require(packet['external_Rm_handles'] == [183696,183697] and packet['external_Rm_coefficient'] == -1,
            'External endpoint changed')
    signed_history = old['actual_signed_native_sources']['histories']['payload']
    require(packet['history_error_receipts'] == signed_history, 'Inherited history errors changed')
    rm_raw = (BASE/'lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json').read_bytes()
    require(digest(rm_raw) == '78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb', 'Rm parent changed')
    rm = json.loads(rm_raw)['actual_directed_Rm_controls']['Rm_endpoint']
    require(new['old_deltaUy_Z0_preserved'] == rm['finite_deltaUy']['log_absolute_upper_bounds'][0]
        and new['old_deltaMy_Z0_Z2_preserved'] == rm['finite_deltaMy']['log_absolute_upper_bounds'][:3]
        and new['old_deltaQy_Z0_preserved'] == rm['finite_deltaQy']['log_absolute_upper_bounds'][0]
        and new['native_Qy_over_Rm_signed_log_scale'] == rm['normalization']['native_Qy_over_Rm_signed_log_scale']
        and new['external_pressure_boundary'] == rm['pressure_boundary'], 'Saved source rows changed')
    require(new['arithmetic']['saved_denominator_correction'] == rm['finite_deltaQy']['operator_receipt']['q1_denominator_correction']
        and new['arithmetic']['full_log_Rm_subtracted_once'] == rm['normalization']['log_Rm'],
        'Denominator correction/radius source changed')
    print(json.dumps(dict(passed=True, report=NAME, sha256=digest(raw), bytes=len(raw),
        canonical_records=records(packet), input_pins=len(expected_inputs), implementation_pins=len(expected_implementations),
        Rm_normalized_C1_endpoint_target_met=target,
        inspection_math_providers_imported=False, accepted_ancestors_replayed=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--compute-new', action='store_true')
    actions.add_argument('--write-report', action='store_true')
    actions.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    if args.compute_new:
        produce()
    elif args.write_report:
        package()
    else:
        inspect()
