"""Checkpoint new linear-m density and conditional transport source bridges."""
from __future__ import annotations

import argparse
import ast
import base64
from datetime import datetime, timezone
import gzip
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
POSITIVE = PREFIX+'Rm_history_a_positivity_sources.json'
POSITIVE_SHA = '4c3b730e57966c0c4a8c2481564a6bf6ec19b1c1136ab01f7e7e2ade1246f675'
LOADER = 'scripts/diagnose_original_Rm_physical_qy0_access.py'
LOADER_SHA = '765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb'
RUNNER = 'scripts/run_original_Rm_m_C2_C5_source_bridge.py'
KERNELS = dict(density='src/openai_ns_reconstruction/original_Rm_m_density_C2_C5_source_bridge.py',
    transport='src/openai_ns_reconstruction/original_Rm_m_C2_C5_transport_bridge.py',
    local='src/openai_ns_reconstruction/original_Rm_local_m_C2_regularity.py')
GRAPH_SHA = 'a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3'
CK = ROOT/'tmp/original_Rm_m_C2_C5_source_bridge'


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


def setup(key):
    require(key in KERNELS, 'unknown component')
    require(sha(BASE/PARENT) == PARENT_SHA and sha(BASE/POSITIVE) == POSITIVE_SHA
        and sha(ROOT/LOADER) == LOADER_SHA, 'saved source parents changed')
    parent = json.loads((BASE/PARENT).read_bytes())
    positive = json.loads((BASE/POSITIVE).read_bytes())
    require(positive['scope']['actual_whole_domain_a_positive_sources_proved'] == 2
        and not positive['scope']['finite_N_history_C2_source_proved']
        and not positive['scope']['physical_Q_y_Z1_value_seam_proved'], 'accepted scope changed')
    for name in (RUNNER, KERNELS[key]):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    loader = module('Rm_m_bridge_saved_graph_loader', ROOT/LOADER)
    nodes, _, _ = loader.load_inputs()
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_SHA, 'actual source graph changed')
    inherited = {**parent['input_hashes'], **parent['implementation_hashes']}
    for name, expected in inherited.items():
        require(sha(loader.path_for(name)) == expected, 'inherited source pin changed: '+name)
    for name, expected in positive['implementation_hashes'].items():
        require(sha(ROOT/name) == expected, 'accepted positivity code changed: '+name)
    kernel = module('new_m_bridge_source_declarations_'+key, ROOT/KERNELS[key])
    required = kernel.RECORD_HASHES
    require(isinstance(required, dict) and required, 'fixed source records missing')
    for name, expected in required.items():
        require(inherited.get(name) == expected, 'source outside accepted parent: '+name)
    records = {name: record(name, expected) for name, expected in required.items()}
    values = dict(source_records=records, expected_nodes_digest=GRAPH_SHA)
    if key == 'local':
        values['positivity_packet'] = positive
    packet = dict(schema='Rm_m_C2_C5_source_bridge_inputs_v1', component=key,
        parent_report=PARENT, parent_sha256=PARENT_SHA, accepted_a_positive_report=POSITIVE,
        accepted_a_positive_sha256=POSITIVE_SHA, graph_digest=GRAPH_SHA, graph_node_count=len(nodes),
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
    kernel = module('actual_new_m_source_bridge_'+key, ROOT/KERNELS[key])
    write_new(entered, dict(source_sha256=sha(source), component=key, input_digest=packet['input_digest'],
        kernel_sha256=packet['implementation_hashes'][KERNELS[key]], entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result = time.perf_counter(), None
    try:
        result = kernel.derive(nodes, **values)
        require(digest(nodes) == GRAPH_SHA and digest(values) == packet['input_digest'], 'source input mutated')
        validate_result(result, packet)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered), runner_full_goal_complete=False)
    write_new(target, result)
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
