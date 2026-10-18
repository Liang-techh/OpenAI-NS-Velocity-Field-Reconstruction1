"""Checkpoint new whole-domain positivity source bindings for history roots."""
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
LEADING = PREFIX+'Rm_leading_m_ordinary_Z_sources.json'
LEADING_SHA = '02cb39ce0272143655566f2a876bfdc1a083ea90c7aaf869506d401ed10ffcf8'
LOADER = 'scripts/diagnose_original_Rm_physical_qy0_access.py'
LOADER_SHA = '765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb'
RUNNER = 'scripts/run_original_Rm_history_a_positivity.py'
KERNELS = dict(postrestore='src/openai_ns_reconstruction/original_Rm_postrestore_a_positive_source.py',
    patch='src/openai_ns_reconstruction/original_Rm_patch_a_positive_source.py')
SCHEMAS = dict(postrestore='original_Rm_postrestore_a_positive_source_v1', patch='original_Rm_patch_a_positive_source_v1')
TARGETS = dict(postrestore=119435, patch=196071)
GRAPH_SHA = 'a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3'
CK = ROOT/'tmp/original_Rm_history_a_positivity'
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


def setup(key):
    require(sha(BASE/PARENT) == PARENT_SHA and sha(BASE/LEADING) == LEADING_SHA
        and sha(ROOT/LOADER) == LOADER_SHA, 'saved parent/source loader changed')
    parent = json.loads((BASE/PARENT).read_bytes())
    leading = json.loads((BASE/LEADING).read_bytes())
    require(leading['actual_Rm_leading_m_ordinary_Z_source']['passed'] is True
        and leading['scope']['physical_Q_y_Z1_value_seam_proved'] is False, 'saved leading scope changed')
    for name in (RUNNER, KERNELS[key]):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    loader = module('Rm_history_a_saved_data_loader', ROOT/LOADER)
    nodes, _, _ = loader.load_inputs()
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_SHA, 'actual graph changed')
    inherited = {**parent['input_hashes'], **parent['implementation_hashes']}
    for name, expected in inherited.items():
        require(sha(loader.path_for(name)) == expected, 'inherited source pin changed: '+name)
    for name, expected in leading['implementation_hashes'].items():
        require(sha(ROOT/name) == expected, 'accepted leading implementation changed: '+name)
    kernel = module('new_history_a_source_declarations_'+key, ROOT/KERNELS[key])
    required = kernel.RECORD_HASHES
    require(isinstance(required, dict) and required, 'fixed source records missing')
    records = {}
    for name, expected in required.items():
        require(inherited.get(name) == expected, 'source record is outside accepted parent: '+name)
        raw = (BASE/name).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == expected, 'source record changed: '+name)
        records[name] = dict(path='experiments/root_st073/'+name, sha256=expected, source_text=raw.decode('utf-8'))
    values = dict(source_records=records, expected_nodes_digest=GRAPH_SHA)
    packet = dict(schema='Rm_history_a_positive_source_inputs_v1',component=key,
        parent_report=PARENT,parent_sha256=PARENT_SHA,accepted_leading_report=LEADING,accepted_leading_sha256=LEADING_SHA,
        graph_digest=GRAPH_SHA,graph_node_count=len(nodes),input_digest=digest(values),record_hashes=required,
        implementation_hashes={name:sha(ROOT/name) for name in (RUNNER,LOADER,KERNELS[key])},
        accepted_mathematical_operators_replayed=False,lower_bound_is_not_function_value=True,
        numeric_field_evaluated=False,full_goal_complete=False)
    return parent, nodes, values, packet


def prepare(key):
    *_, packet = setup(key)
    require(not any(p.exists() for p in (SOURCES[key],RESULTS[key],ENTRIES[key])), 'attempt already prepared')
    write_new(SOURCES[key],packet)
    return dict(passed=True,component=key,source_sha256=sha(SOURCES[key]),input_digest=packet['input_digest'])


def validate_result(key,result):
    require(isinstance(result,dict) and result.get('schema') == SCHEMAS[key] and type(result.get('passed')) is bool,
        'source positivity schema changed')
    require(result.get('source_only') is True and result.get('nodes_digest') == GRAPH_SHA
        and result.get('actual_target_handle') == TARGETS[key], 'actual source target changed')
    require(result.get('whole_native_domain') == (['0','1'] if key == 'postrestore' else ['1','exp(1)'])
        and result.get('whole_Z_domain') == ['-1','1'], 'whole source domain changed')
    if result['passed']:
        require(result.get('strict_positive_proved') is True and not result.get('unresolved_gates')
            and result.get('lower_bound_exact') in (('4/5',) if key == 'postrestore' else ('7/10','1/10')),
            'source lower bound promoted without proof')
    for flag in ('numeric_values_computed','old_operators_replayed','full_goal_complete'):
        require(result.get(flag) is False,'source scope changed: '+flag)


def run_component(key):
    _,nodes,values,packet = setup(key)
    require(json.loads(SOURCES[key].read_bytes()) == packet and not ENTRIES[key].exists()
        and not RESULTS[key].exists(), 'attempt entered or changed')
    kernel = module('actual_new_history_a_positive_'+key,ROOT/KERNELS[key])
    write_new(ENTRIES[key],dict(source_sha256=sha(SOURCES[key]),component=key,input_digest=packet['input_digest'],
        kernel_sha256=packet['implementation_hashes'][KERNELS[key]],entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started,result = time.perf_counter(),None
    try:
        result = kernel.derive(nodes,**values)
        require(digest(nodes) == GRAPH_SHA and digest(values) == packet['input_digest'],'source input mutated')
        validate_result(key,result)
    except Exception as error:
        result = dict(passed=False,failed_attempt_preserved=True,error_type=type(error).__name__,error=str(error),
            traceback=''.join(traceback.format_exception(error)),rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started,source_checkpoint_sha256=sha(SOURCES[key]),
        entry_checkpoint_sha256=sha(ENTRIES[key]),runner_full_goal_complete=False)
    write_new(RESULTS[key],result)
    return dict(component=key,passed=result['passed'],sha256=sha(RESULTS[key]),seconds=result['execution_seconds'],
        lower_bound_exact=result.get('lower_bound_exact'),error=result.get('error'),unresolved_gates=result.get('unresolved_gates',[]))


def saved_component(key):
    parent,_,values,packet = setup(key)
    result,entry = json.loads(RESULTS[key].read_bytes()),json.loads(ENTRIES[key].read_bytes())
    require(json.loads(SOURCES[key].read_bytes()) == packet
        and result['source_checkpoint_sha256'] == entry['source_sha256'] == sha(SOURCES[key])
        and result['entry_checkpoint_sha256'] == sha(ENTRIES[key]) and entry['component'] == key
        and entry['input_digest'] == packet['input_digest']
        and entry['kernel_sha256'] == packet['implementation_hashes'][KERNELS[key]],'saved source/entry/result changed')
    if not result.get('failed_attempt_preserved'):
        validate_result(key,result)
    else:
        require(result['passed'] is False and isinstance(result.get('traceback'),str),'malformed failure')
    return parent,packet,result


def inspect(key):
    _,_,result = saved_component(key)
    return dict(passed=True,saved_only=True,component=key,positive_source_passed=result['passed'],
        result_sha256=sha(RESULTS[key]),full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare_postrestore','prepare_patch','postrestore','patch','inspect_postrestore','inspect_patch'))
    action = parser.parse_args().action
    print(json.dumps(prepare(action.removeprefix('prepare_')) if action.startswith('prepare_')
        else inspect(action.removeprefix('inspect_')) if action.startswith('inspect_') else run_component(action)))
