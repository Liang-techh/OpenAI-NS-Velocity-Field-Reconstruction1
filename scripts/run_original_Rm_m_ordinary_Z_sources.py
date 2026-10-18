"""Checkpoint new leading/history ordinary-Z source bindings only."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import importlib.util
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
LOADER = 'scripts/diagnose_original_Rm_physical_qy0_access.py'
LOADER_SHA = '765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb'
RUNNER = 'scripts/run_original_Rm_m_ordinary_Z_sources.py'
KERNELS = dict(leading='src/openai_ns_reconstruction/original_Rm_leading_m_ordinary_Z_source.py',
    history='src/openai_ns_reconstruction/original_Rm_history_m_ordinary_Z_source.py')
SCHEMAS = dict(leading='original_Rm_leading_m_ordinary_Z_source_v1', history='original_Rm_history_m_ordinary_Z_source_v1')
GRAPH_SHA = 'a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3'
CK = ROOT/'tmp/original_Rm_m_ordinary_Z_sources'
SOURCES = {key:CK/('new_'+key+'_source.json') for key in KERNELS}
RESULTS = {key:CK/('new_'+key+'_result.json') for key in KERNELS}
ENTRIES = {key:CK/('new_'+key+'_entered.json') for key in KERNELS}


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


def setup(component=None):
    require(sha(BASE/PARENT) == PARENT_SHA and sha(ROOT/LOADER) == LOADER_SHA, 'parent or data loader changed')
    parent = json.loads((BASE/PARENT).read_bytes())
    require(parent['actual_Rm_physical_Qy_Z_derivative_provenance']['outer_calculus_proved'] is True
        and parent['scope']['physical_Q_y_Z1_value_seam_proved'] is False, 'saved baseline scope changed')
    selected = KERNELS if component is None else {component:KERNELS[component]}
    for name in (RUNNER, *selected.values()):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    loader = module('Rm_m_Z_saved_data_loader', ROOT/LOADER)
    nodes, _, _ = loader.load_inputs()
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_SHA, 'actual graph changed')
    inherited = {**parent['input_hashes'], **parent['implementation_hashes']}
    for name, expected in inherited.items():
        require(sha(loader.path_for(name)) == expected, 'inherited source pin changed: '+name)
    inputs, record_hashes = {}, {}
    for key, path in selected.items():
        kernel = module('new_m_Z_source_declarations_'+key, ROOT/path)
        required = kernel.RECORD_HASHES
        require(isinstance(required, dict) and required, 'fixed source records missing')
        records = {}
        for name, expected in required.items():
            require(inherited.get(name) == expected, 'source record is outside accepted parent: '+name)
            raw = (BASE/name).read_bytes()
            require(hashlib.sha256(raw).hexdigest() == expected, 'source record changed: '+name)
            records[name] = dict(path='experiments/root_st073/'+name, sha256=expected, source_text=raw.decode('utf-8'))
        inputs[key] = dict(source_records=records, expected_nodes_digest=GRAPH_SHA)
        record_hashes[key] = dict(required)
    packet = dict(schema='Rm_m_ordinary_Z_source_inputs_v1', parent_report=PARENT, parent_sha256=PARENT_SHA,
        graph_digest=GRAPH_SHA, graph_node_count=len(nodes), input_digests={key:digest(value) for key,value in inputs.items()},
        record_hashes=record_hashes, implementation_hashes={name:sha(ROOT/name) for name in (RUNNER, LOADER, *selected.values())},
        accepted_mathematical_operators_replayed=False, numeric_field_evaluated=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def prepare(component):
    *_, packet = setup(component)
    require(not any(p.exists() for p in (SOURCES[component], RESULTS[component], ENTRIES[component])), 'attempt already prepared')
    write_new(SOURCES[component], packet)
    return dict(passed=True, component=component, source_sha256=sha(SOURCES[component]), input_digests=packet['input_digests'])


def validate_result(key, result):
    require(isinstance(result, dict) and result.get('schema') == SCHEMAS[key] and type(result.get('passed')) is bool,
        'source result schema changed')
    require(result.get('source_only') is True and result.get('nodes_digest') == GRAPH_SHA, 'source result graph changed')
    rows = result.get('rows', [])
    expected = dict(leading={'postrestore':[96436,96437,96438], 'patch':[195000,195001,195002]},
        history={'postrestore':[121602,121603,121604], 'patch':[198466,198467,198468]})[key]
    require(len(rows) == 4, 'four source edges required')
    require({(r.get('side'),r.get('ordinary_order'),r.get('source_handle'),r.get('derivative_handle')) for r in rows}
        == {(side,order,ports[order],ports[order+1]) for side,ports in expected.items() for order in (0,1)}, 'actual source edge changed')
    if result['passed']:
        require(all(r.get('source_derivative_proved') is True and r.get('regularity_proved') is True for r in rows)
            and not result.get('unresolved_gates'), 'open source edge promoted')
    for flag in ('numeric_values_computed', 'old_operators_replayed', 'full_goal_complete'):
        require(result.get(flag) is False, 'source execution/scope flag changed: '+flag)


def run_component(key):
    _, nodes, inputs, packet = setup(key)
    source = SOURCES[key]
    require(json.loads(source.read_bytes()) == packet, 'frozen foundation changed')
    require(not ENTRIES[key].exists() and not RESULTS[key].exists(), 'new source already entered')
    values = inputs[key]
    kernel = module('actual_new_m_Z_source_'+key, ROOT/KERNELS[key])
    write_new(ENTRIES[key], dict(source_sha256=sha(source), component=key, input_digest=digest(values),
        kernel_sha256=packet['implementation_hashes'][KERNELS[key]], entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result = time.perf_counter(), None
    try:
        result = kernel.derive(nodes, **values)
        require(digest(nodes) == GRAPH_SHA and digest(values) == packet['input_digests'][key], 'source input mutated')
        validate_result(key, result)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(ENTRIES[key]), runner_full_goal_complete=False)
    write_new(RESULTS[key], result)
    return dict(component=key, passed=result['passed'], sha256=sha(RESULTS[key]), seconds=result['execution_seconds'],
        error=result.get('error'), unresolved_gates=result.get('unresolved_gates', []))


def saved():
    results, packets = {}, {}
    for key, path in RESULTS.items():
        parent, _, inputs, packet = setup(key)
        require(json.loads(SOURCES[key].read_bytes()) == packet, 'saved foundation changed')
        result, entry = json.loads(path.read_bytes()), json.loads(ENTRIES[key].read_bytes())
        require(result['source_checkpoint_sha256'] == entry['source_sha256'] == sha(SOURCES[key])
            and result['entry_checkpoint_sha256'] == sha(ENTRIES[key]) and entry['component'] == key
            and entry['input_digest'] == packet['input_digests'][key]
            and entry['kernel_sha256'] == packet['implementation_hashes'][KERNELS[key]], 'saved source/entry/result changed')
        if not result.get('failed_attempt_preserved'):
            validate_result(key, result)
        else:
            require(result['passed'] is False and isinstance(result.get('traceback'), str), 'malformed failure')
        results[key] = result
        packets[key] = packet
    return parent, dict(parent_sha256=PARENT_SHA, components=packets,
        implementation_hashes={name:value for packet in packets.values() for name,value in packet['implementation_hashes'].items()}), results


def inspect():
    _, packet, results = saved()
    return dict(passed=True, saved_only=True, parent_sha256=packet['parent_sha256'],
        components={key:dict(passed=result['passed'], sha256=sha(RESULTS[key])) for key,result in results.items()},
        full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_leading', 'prepare_history', 'leading', 'history', 'inspect'))
    action = parser.parse_args().action
    print(json.dumps(prepare(action.removeprefix('prepare_')) if action.startswith('prepare_')
        else run_component(action) if action in KERNELS else inspect()))
