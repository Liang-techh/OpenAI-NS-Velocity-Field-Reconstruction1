"""Compute new loop controls from saved corrected native data, or inspect JSON.

Inspection performs only tuple, hash, AST and structural checks. Production
checkpoints the new data before validation so a validator repair cannot require
repeating an accepted mathematical computation.
"""
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
PARENT_NAME = 'lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json'
PARENT_SHA = '10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a'
PARENT_BYTES = 42177499
NAME = 'lei_ren_part1_paper_compliant_current_original_directed_loop_controls.json'
REPORT = BASE / NAME
PENDING = BASE / (NAME + '.pending')
MODULE = 'src/openai_ns_reconstruction/original_directed_loop_controls.py'
RUNNER = 'scripts/run_original_directed_loop_controls.py'
C5 = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
FIXTURE_NAME = 'original_directed_loop_controls_fresh_fixture.json'
FIXTURE_SHA = '4670cf41b2547fc0a6895ff7bd9229955ad48ca8c0b174061d270bc6f8c4cb6e'
REVIEWED = {
    MODULE: '9f2cd71e4de45468f89bb2a3bbc6bdcb9e7932090b6c943944c843d7084dcec8',
    C5: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    'src/openai_ns_reconstruction/original_log_Z_phase_bounds.py': '21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59',
    'src/openai_ns_reconstruction/original_log_mixed_jets.py': 'ceeebf88a02207b02fedfdff727bc2ae1d09f4a52ebf27b8a84ff8876719361a',
}
CHARTS = ('long_reshape', 'reference', 'restoration', 'post_restore')
ROWS = ('A', 'B_over_E', 'B_over_Pstar_over_frozen_E0', 'A_phi',
        'B_over_Pstar_over_frozen_E0_phi', 'inverse_Psi')
FALSE_SCOPE = ('full_C5_Z5_bounds', 'physical_y_bounds', 'current_N_revalidated',
               'full_pressure_revalidated', 'global_five_moment_closure', 'full_goal_complete',
               'legacy_loop_packets_recertified', 'inherited_log_caps_revalidated')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(name, *, input_pin=False):
    path = Path(name)
    require(not path.is_absolute(), 'Absolute source pin')
    if '..' in path.parts:
        require(input_pin and path.as_posix() == '../../src/openai_ns_reconstruction/schedule_pressure.py',
                'Unregistered source traversal')
        return (BASE / path).resolve()
    if input_pin and not path.parts[0] == 'experiments':
        return BASE / path.name
    return ROOT / path


def parent_data():
    path = BASE / PARENT_NAME
    raw = path.read_bytes()
    require(len(raw) == PARENT_BYTES and hashlib.sha256(raw).hexdigest() == PARENT_SHA,
            'Saved corrected native parent changed')
    parent = json.loads(raw)
    for group in ('input_hashes', 'implementation_hashes'):
        for name, expected in parent[group].items():
            require(sha(source_path(name, input_pin=group == 'input_hashes')) == expected,
                    'Parent source pin changed: ' + name)
    packet = parent['actual_native_directed_controls_v2']
    require(packet['scope']['corrected_native_controls_verified'] is True,
            'Corrected native premise missing')
    require(packet['scope']['inherited_log_caps_revalidated'] is False,
            'Legacy caps must remain unrevalidated')
    raw_tree = packet['raw_interval_tree']
    cut_rows = packet['cut_numeric_receipt']['rows']
    require({row['chart'] for row in cut_rows} == set(CHARTS) and len(cut_rows) == 4,
            'Cut evidence does not cover four charts')
    for row in cut_rows:
        native = raw_tree[row['section']]['charts'][row['chart']]
        require(row['eta_source_log'] == raw_tree['shared']['logeta']
                and row['kappa_log_cap'] == native['root_caps']['kappa'][0]['log_absolute_upper_bound'],
                'Cut evidence is not bound to the current source')
        for key in ('eta_source_log', 'kappa_log_cap', 'numeric_lower_bound'):
            record(row[key])
        require(mpf_cmp(endpoint(row['eta_source_log']['upper_exact_mpf_tuple']), (0,0,0,0)) <= 0
                and mpf_cmp(endpoint(row['kappa_log_cap']['upper_exact_mpf_tuple']), (0,0,0,0)) < 0
                and mpf_cmp(endpoint(row['numeric_lower_bound']['lower_exact_mpf_tuple']), (0,1,1,1)) >= 0,
                'Strict active flat-cut inequality failed')
        for key in ('positive_a_lower', 'positive_gamma_lower'):
            record(native['conditioning'][key])
            require(mpf_cmp(endpoint(native['conditioning'][key]['lower_exact_mpf_tuple']), (0,0,0,0)) > 0,
                    'Positive root premise failed: ' + key)
        for key in ('log_D_lower','log_lambda_lower'):
            record(native['conditioning'][key])
    return parent


def endpoint(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4
            and all(type(part) is int for part in value), 'Malformed MPF tuple')
    sign, mantissa, exponent, bits = value
    require(sign in (0, 1) and mantissa >= 0 and bits == mantissa.bit_length(),
            'Nonfinite or malformed MPF tuple')
    require((mantissa == 0 and tuple(value) == (0, 0, 0, 0))
            or (mantissa > 0 and mantissa & 1), 'Noncanonical MPF tuple')
    return tuple(value)


def record(value):
    require(isinstance(value, dict) and value.get('directed') is True,
            'Directed interval required')
    low = endpoint(value['lower_exact_mpf_tuple'])
    high = endpoint(value['upper_exact_mpf_tuple'])
    require(mpf_cmp(low, high) <= 0, 'Reversed interval')
    return 1


def record_count(tree):
    if isinstance(tree, dict):
        if 'lower_exact_mpf_tuple' in tree or 'upper_exact_mpf_tuple' in tree:
            return record(tree)
        return sum(record_count(value) for value in tree.values())
    if isinstance(tree, list):
        return sum(record_count(value) for value in tree)
    return 0


def rows(values, expected_orders):
    require(isinstance(values, list) and len(values) == len(expected_orders), 'Incomplete rows')
    for index, row in zip(expected_orders, values):
        require(type(row.get('exact_zero')) is bool, 'Exact-zero marker missing')
        if isinstance(index, tuple):
            require((row.get('Z_order'), row.get('free_angle_order')) == index,
                    'Unordered phase rows')
        else:
            require(row.get('ordinary_Z_order') == index, 'Unordered ordinary rows')
        bound = row.get('log_absolute_upper_bound')
        if row['exact_zero']:
            require(bound is None, 'Exact-zero row carries a bound')
        else:
            record(bound)


def source_guard():
    raw = (ROOT / C5).read_bytes()
    tree = ast.parse(raw, filename=C5)
    node = next(node for node in tree.body
                if isinstance(node, ast.FunctionDef) and node.name == 'loop_functions')
    dumped = ast.dump(node, include_attributes=False)
    return dict(path=C5, source_sha256=hashlib.sha256(raw).hexdigest(),
                function='loop_functions', ast_sha256=hashlib.sha256(dumped.encode()).hexdigest(),
                source_executed=False)


def validate(parent, report):
    require(report['parent_report'] == PARENT_NAME and report['parent_report_sha256'] == PARENT_SHA,
            'Parent binding changed')
    require(report['input_hashes'] == dict(parent['input_hashes'], **{PARENT_NAME:PARENT_SHA,FIXTURE_NAME:FIXTURE_SHA}),
            'Input pins changed')
    require(sha(BASE/FIXTURE_NAME) == FIXTURE_SHA, 'Fresh fixture receipt changed')
    fixture = json.loads((BASE/FIXTURE_NAME).read_bytes())
    require(fixture['passed'] is True and all(value is True for value in fixture['checks'].values())
            and report['fresh_fixture'] == fixture, 'Fresh fixture evidence changed')
    implementations = dict(parent['implementation_hashes'])
    implementations.update({MODULE:sha(ROOT/MODULE), RUNNER:sha(ROOT/RUNNER)})
    require(report['implementation_hashes'] == implementations, 'Implementation pins changed')
    for path, expected in REVIEWED.items():
        require(expected != 'REVIEW_PENDING' and sha(ROOT/path) == expected,
                'Reviewed source changed: ' + path)
    require(report['source_guard'] == source_guard(), 'C5 loop source guard changed')
    require(report['historical_payload_retained_by_parent_hash'] is True,
            'Historical data preservation link missing')
    scope = report['scope']
    require(scope['fixed_target_phase_Z0_through_Z3_verified'] is True,
            'New loop scope missing')
    require(all(scope[key] is False for key in FALSE_SCOPE), 'Downstream scope overclaimed')
    packet = report['actual_directed_loop_controls']
    require(packet['passed'] is True
            and packet['native_parent_report_sha256'] == PARENT_SHA,
            'Component source binding missing')
    require(packet['scope']['fixed_target_phase_Z0_through_Z3_verified'] is True
            and all(packet['scope'][key] is False for key in FALSE_SCOPE),
            'Component scope missing or overclaimed')
    require(set(packet['charts']) == set(CHARTS), 'Wrong chart set')
    raw = parent['actual_native_directed_controls_v2']['raw_interval_tree']
    phase_orders = tuple((k, j) for total in range(4) for k in range(4)
                         for j in range(4) if k+j == total)
    for name, chart in packet['charts'].items():
        section = 'long_reference' if name in CHARTS[:2] else 'restoration'
        native = raw[section]['charts'][name]
        require(chart['input_conditioning'] == native['conditioning'],
                'Conditioning was recomputed or replaced: ' + name)
        require(chart['input_root_caps'] == native['root_caps'],
                'Root inputs were recomputed or replaced: ' + name)
        require(chart['input_E_relative'] == native['E_relative'],
                'Amplitude input changed: ' + name)
        require(chart['conditioning'] == chart['input_conditioning']
                and chart['E_relative'] == native['E_relative'][:4],
                'Public conditioning or relative amplitude changed: ' + name)
        for field, values in chart['root_caps'].items():
            require(values == native['root_caps'][field],
                    'Public root cap differs from its input: ' + name + '.' + field)
        for field in ROWS:
            rows(chart[field], tuple(range(4)))
        for field, values in chart['phase_caps'].items():
            rows(values, phase_orders)
        require(chart['fixed_target_phase'] is True
                and chart['signed_value_selected'] is False, 'Selected phase or cap value')
    audit = report['execution_audit']
    require(audit['accepted_ancestors_replayed'] is False
            and audit['native_builder_replayed'] is False
            and audit['old_loop_builders_called'] is False, 'Accepted stage replay')
    return record_count(packet)


def write_json(path, value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',
                          encoding='utf-8',newline='\n')


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write-report', action='store_true')
    mode.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    parent = parent_data()
    if args.write_report:
        require(not REPORT.exists() and not PENDING.exists(),
                'Saved or pending data exists; computation replay refused')
        for path, expected in REVIEWED.items():
            require(expected != 'REVIEW_PENDING' and sha(ROOT/path) == expected,
                    'Source review pending or changed: ' + path)
        sys.path.insert(0,str(ROOT/'src'))
        from mpmath.ctx_iv import MPIntervalContext
        from openai_ns_reconstruction import original_directed_loop_controls as provider
        c = MPIntervalContext()
        c.prec = 1600
        started = time.monotonic()
        packet = provider.compute(c,(BASE/PARENT_NAME).read_bytes())
        inputs = dict(parent['input_hashes'])
        inputs[PARENT_NAME] = PARENT_SHA
        require(sha(BASE/FIXTURE_NAME) == FIXTURE_SHA, 'Fresh fixture receipt changed')
        inputs[FIXTURE_NAME] = FIXTURE_SHA
        implementations = dict(parent['implementation_hashes'])
        implementations.update({MODULE:sha(ROOT/MODULE),RUNNER:sha(ROOT/RUNNER)})
        report = dict(stage='current_original_directed_loop_controls',
            parent_report=PARENT_NAME,parent_report_sha256=PARENT_SHA,
            input_hashes=inputs,implementation_hashes=implementations,
            source_guard=source_guard(),historical_payload_retained_by_parent_hash=True,
            fresh_fixture=json.loads((BASE/FIXTURE_NAME).read_bytes()),
            actual_directed_loop_controls=packet,
            scope=dict(fixed_target_phase_Z0_through_Z3_verified=True,
                       **{key:False for key in FALSE_SCOPE}),
            execution_audit=dict(accepted_ancestors_replayed=False,native_builder_replayed=False,
                old_loop_builders_called=False,only_new_loop_inputs_propagated=True),
            execution_seconds=time.monotonic()-started)
        write_json(PENDING, report)
        count = validate(parent,report)
        PENDING.rename(REPORT)
    else:
        report = json.loads(REPORT.read_bytes())
        count = validate(parent,report)
    print(json.dumps(dict(passed=True,report=str(REPORT),sha256=sha(REPORT),
        bytes=REPORT.stat().st_size,canonical_records=count,
        input_pins=len(report['input_hashes']),implementation_pins=len(report['implementation_hashes']),
        scope=report['scope'],execution_seconds=report['execution_seconds'])))


if __name__ == '__main__':
    main()
