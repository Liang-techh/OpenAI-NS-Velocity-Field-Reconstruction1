"""Checkpoint R100, both switches and the post-switch power source projection."""
from __future__ import annotations

import argparse
import ast
import base64
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
from inspect import signature
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'experiments/root_st073'
PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX+'Rm_physical_Qy_Z_calculus.json'
PARENT_SHA = '8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746'
POSITIVE = PREFIX+'macro_evolution_source_projection.json'
POSITIVE_SHA = '880abe2721144f86b0a9775f62735dccc8c5650a5fd8faab27613997daa4ae16'
LOADER = 'scripts/diagnose_original_Rm_physical_qy0_access.py'
LOADER_SHA = '765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb'
RUNNER = 'scripts/run_original_R100_switch_C2_C6_source_projection_v2.py'
KERNELS = dict(switch_projection='src/openai_ns_reconstruction/original_R100_switch_C2_C6_source_projection_v2.py')
GRAPH_SHA = 'a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3'
CK = ROOT/'tmp/original_R100_switch_C2_C6_source_projection_v2'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def write_new(path, value):
    require(not path.exists(), 'artifact exists: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write('\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def paths(key):
    return tuple(CK/('new_'+key+'_'+kind+'.json') for kind in ('source', 'entered', 'result'))


def record(name, expected):
    path = BASE/name
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, 'source record changed: '+name)
    # Raw archive authentication is deliberately before every decode.
    decoded = gzip.decompress(raw) if name.endswith('.gz') else raw
    value = dict(path='experiments/root_st073/'+name, sha256=expected, source_text=decoded.decode('utf-8'))
    if name.endswith('.gz'):
        value.update(encoding='gzip', raw_bytes_base64=base64.b64encode(raw).decode('ascii'),
            decoded_sha256=hashlib.sha256(decoded).hexdigest())
    return value



def frozen_prior_attempt():
    code_hashes = {
        'scripts/run_original_R100_switch_C2_C6_source_projection.py': '2c233660f3d77fa905599c90f1657dba849f64371b6f888b318cc51fa42faaa6',
        'src/openai_ns_reconstruction/original_R100_switch_C2_C6_source_projection.py': '0dc182f8cbe050cf4a0c1d957001713654fa2388a333818a1c9bf5bee9b00754',
    }
    old_base = 'tmp/original_R100_switch_C2_C6_source_projection/'
    checkpoint_hashes = {
        old_base+'new_switch_projection_source.json': '0c00a72d64068698ff74d6a7832a739a14f83e9f300344f70fe22b46a7f3902f',
        old_base+'new_switch_projection_entered.json': '57cb171236ca002a7811cd927d5c0a39829a91dc0930a7a7fe590c28fc69c991',
        old_base+'new_switch_projection_result.json': '8633696809197e44d14f039a4294ccf19ab52c9b2fc30e27eadd506d76cd4a75',
    }
    for name, expected in {**code_hashes, **checkpoint_hashes}.items():
        require(sha(ROOT/name) == expected, 'immutable prior attempt changed: '+name)
    source = json.loads((ROOT/(old_base+'new_switch_projection_source.json')).read_bytes())
    entry = json.loads((ROOT/(old_base+'new_switch_projection_entered.json')).read_bytes())
    result = json.loads((ROOT/(old_base+'new_switch_projection_result.json')).read_bytes())
    require(result.get('passed') is False and result.get('error') is None
        and not result.get('failed_attempt_preserved')
        and result.get('unresolved_gates') == ['ast_witness_missing:second_sigma_prefix_and_phase_derivative_identity']
        and result.get('numeric_values_computed') is False
        and result.get('old_operators_replayed') is False, 'prior normal-negative scope changed')
    require(source['implementation_hashes'][list(code_hashes)[0]] == code_hashes[list(code_hashes)[0]]
        and source['implementation_hashes'][list(code_hashes)[1]] == code_hashes[list(code_hashes)[1]]
        and entry['source_sha256'] == result['source_checkpoint_sha256'] == checkpoint_hashes[old_base+'new_switch_projection_source.json']
        and result['entry_checkpoint_sha256'] == checkpoint_hashes[old_base+'new_switch_projection_entered.json'], 'prior attempt checkpoint link changed')
    return dict(implementation_hashes=code_hashes, checkpoint_hashes=checkpoint_hashes,
        result_sha256=checkpoint_hashes[old_base+'new_switch_projection_result.json'],
        passed=False, unresolved_gates=result['unresolved_gates'],
        source_component_calls=1, accepted_ancestor_operator_calls=0,
        numeric_quadrature_recurrence_or_norm_calls=0,
        correction='Use structural AST equality for unary-negative expressions; source equations and acceptance conditions are unchanged.')


def setup(key):
    require(key in KERNELS, 'unknown component')
    prior_attempt = frozen_prior_attempt()
    require(sha(BASE/PARENT) == PARENT_SHA and sha(BASE/POSITIVE) == POSITIVE_SHA
        and sha(ROOT/LOADER) == LOADER_SHA, 'saved source parents changed')
    parent = json.loads((BASE/PARENT).read_bytes())
    positive = json.loads((BASE/POSITIVE).read_bytes())
    require(isinstance(positive, dict) and isinstance(positive.get('scope'), dict)
        and isinstance(positive.get('component_result'), dict), 'accepted macro stage mapping required')
    scope, macro = positive['scope'], positive['component_result']
    ports = macro.get('exact_macro_inlet_function_ports', {})
    require(macro.get('passed') is True
        and macro.get('full_macro_background_evolution_projection_proved') is True
        and macro.get('full_macro_own_history_projection_proved') is True
        and macro.get('nodes_digest') == GRAPH_SHA
        and ports.get('graph_ports_valid') is True and ports.get('C2_source_projection_admitted') is True
        and scope.get('ordinary_Z_orders') == [0, 1, 2]
        and scope.get('same_density_certificate') is False
        and scope.get('finite_N_history_identity') is False
        and scope.get('physical_Q_y_Z1') is False
        and scope.get('global_five_moment_closure') is False
        and scope.get('full_goal_complete') is False,
        'accepted macro source stage scope changed')
    for name in (RUNNER, KERNELS[key]):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    loader = module('Rm_R100_switch_saved_graph_loader', ROOT/LOADER)
    nodes, _, _ = loader.load_inputs()
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_SHA, 'actual source graph changed')
    inherited = {}
    for pins in (parent['input_hashes'], parent['implementation_hashes'],
                 positive['input_hashes'], positive['implementation_hashes'], {PARENT:PARENT_SHA, POSITIVE:POSITIVE_SHA}):
        require(isinstance(pins, dict), 'source pins must be a mapping')
        for name, expected in pins.items():
            require(isinstance(name, str) and isinstance(expected, str)
                and (name not in inherited or inherited[name] == expected), 'source pin conflict: '+str(name))
            inherited[name] = expected
    for name, expected in inherited.items():
        require(sha(loader.path_for(name)) == expected, 'inherited source pin changed: '+name)
    for name, expected in positive['implementation_hashes'].items():
        require(sha(ROOT/name) == expected, 'accepted source code changed: '+name)
    kernel = module('new_R100_switch_source_declarations_'+key, ROOT/KERNELS[key])
    required = kernel.RECORD_HASHES
    require(isinstance(required, dict) and required, 'fixed source records missing')
    for name, expected in required.items():
        require(inherited.get(name) == expected, 'source outside accepted parent: '+name)
    records = {name: record(name, expected) for name, expected in required.items()}
    values = dict(source_records=records, expected_nodes_digest=GRAPH_SHA)
    packet = dict(schema='R100_switch_C2_C6_source_inputs_v2', component=key, prior_attempt=prior_attempt,
        parent_report=PARENT, parent_sha256=PARENT_SHA, accepted_source_stage_report=POSITIVE,
        accepted_source_stage_sha256=POSITIVE_SHA, graph_digest=GRAPH_SHA, graph_node_count=len(nodes),
        input_digest=digest(values), record_hashes=required,
        implementation_hashes={name: sha(ROOT/name) for name in (RUNNER, LOADER, KERNELS[key])},
        decoded_source_digests={name: hashlib.sha256(value['source_text'].encode()).hexdigest() for name, value in records.items()},
        schema_expected=kernel.SCHEMA, accepted_mathematical_operators_replayed=False,
        lower_bounds_not_function_values=True, numeric_field_evaluated=False, full_goal_complete=False)
    return parent, nodes, values, packet


def prepare(key):
    *_, packet = setup(key)
    source, entered, result = paths(key)
    require(not any(path.exists() for path in (source, entered, result)), 'attempt already prepared')
    write_new(source, packet)
    return dict(passed=True, component=key, source_sha256=sha(source), input_digest=packet['input_digest'])


def validate_result(result, packet):
    require(isinstance(result, dict) and result.get('schema') == packet['schema_expected']
        and type(result.get('passed')) is bool, 'source bridge schema changed')
    require(result.get('source_only') is True and result.get('nodes_digest') == GRAPH_SHA, 'actual source graph changed')
    for flag in ('numeric_values_computed', 'old_operators_replayed', 'full_goal_complete'):
        require(result.get(flag) is False, 'source scope changed: '+flag)
    if result['passed']:
        require(not result.get('unresolved_gates'), 'source bridge passed with unresolved gates')


def run_component(key):
    _, nodes, values, packet = setup(key)
    source, entered, target = paths(key)
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists(), 'attempt entered or changed')
    kernel = module('actual_new_R100_switch_source_'+key, ROOT/KERNELS[key])
    require(isinstance(values, dict) and set(values) == {'source_records', 'expected_nodes_digest'}
        and isinstance(values['source_records'], dict), 'explicit source argument mapping required')
    signature(kernel.derive).bind(nodes, **values)
    write_new(entered, dict(source_sha256=sha(source), component=key, input_digest=packet['input_digest'],
        kernel_sha256=packet['implementation_hashes'][KERNELS[key]], entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result = time.perf_counter(), None
    try:
        result = kernel.derive(nodes, **values)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered), runner_full_goal_complete=False)
    write_new(target, result)
    # Persist the returned result before validating or inspecting its content.
    require(digest(nodes) == GRAPH_SHA and digest(values) == packet['input_digest'], 'source input mutated')
    if not result.get('failed_attempt_preserved'):
        validate_result(result, packet)
    return dict(component=key, passed=result['passed'], sha256=sha(target), seconds=result['execution_seconds'],
        error=result.get('error'), unresolved_gates=result.get('unresolved_gates', []))


def saved_component(key):
    parent, _, _, packet = setup(key)
    source, entered, target = paths(key)
    result, entry = json.loads(target.read_bytes()), json.loads(entered.read_bytes())
    require(json.loads(source.read_bytes()) == packet and entry['component'] == key
        and result['source_checkpoint_sha256'] == entry['source_sha256'] == sha(source)
        and result['entry_checkpoint_sha256'] == sha(entered)
        and entry['input_digest'] == packet['input_digest']
        and entry['kernel_sha256'] == packet['implementation_hashes'][KERNELS[key]], 'saved source/entry/result changed')
    if not result.get('failed_attempt_preserved'):
        validate_result(result, packet)
    else:
        require(result['passed'] is False and isinstance(result.get('traceback'), str), 'malformed failure')
    return parent, packet, result


def inspect(key):
    _, _, result = saved_component(key)
    return dict(passed=True, saved_only=True, component=key, source_bridge_passed=result['passed'],
        result_sha256=sha(paths(key)[2]), full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=tuple(action for key in KERNELS for action in ('prepare_'+key, key, 'inspect_'+key)))
    action = parser.parse_args().action
    print(json.dumps(prepare(action.removeprefix('prepare_')) if action.startswith('prepare_')
        else inspect(action.removeprefix('inspect_')) if action.startswith('inspect_') else run_component(action)))
