"""Freeze and enter one J3 successor; persist stage diagnostics before validation."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUNNER = 'scripts/run_original_native_J3_endpoint_transport_v2.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_J3_endpoint_transport_v2.py'
PARENT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport.json'
PARENT_SHA = 'cc276ae3ab77c8cd6730e4f44b3ea60e95ea4af4b1e50d2153915781b26b16d0'
CK = ROOT / 'tmp/original_native_J3_endpoint_transport_v2'
COMPONENT = 'endpoint_transport'
INHERITED_SOURCE_COUNT = 1857
FAILED_ATTEMPT_HASHES = {'scripts/run_original_native_J2_U_endpoint_transport.py': '5b10fa16db583ff5ccbed09c1c57c1bf840feb85c244a5f3263975e6c0d0e076', 'src/openai_ns_reconstruction/original_native_J2_U_endpoint_transport.py': '582f6e95125027b979c91c87c685b640ec2306f3080a10d70b5b2a8e7bdfaed0', 'tmp/original_native_J2_U_endpoint_transport/new_J2_U_endpoint_transport_source.json': '9b3a523c78a63005ca67f793fcfc556a4cfbec0fabda9cb9ffe72f6ea634c105', 'tmp/original_native_J2_U_endpoint_transport/new_J2_U_endpoint_transport_entered.json': 'b849260ac12d818b224853e48bb99544d179d3bd7eeef7f73b7bdbe020c17ac1', 'tmp/original_native_J2_U_endpoint_transport/new_J2_U_endpoint_transport_result.json': '72fd536f035c36366a0031801db114cb19b06250ca868a3a13358f702694dbc3', 'scripts/run_original_native_micro_macro_B_jump_source.py': '583d26e94a4ec25913094f497eb442df6ab46c7201bd1a9521cbf47098f909f0', 'src/openai_ns_reconstruction/original_native_micro_macro_B_jump_source.py': 'c934743d502afe0d247862acf302f637a6f81f5686807a8dc45efacdbc268d3d', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_source.json': 'e0d7782d6bfc37d11c75ae08e05f7fb30588783bc3a7df6e9ceac6721e6c9f58', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_entered.json': '6b2815660d2b0d55717cddb4d576b9e23ef781ad901a713127d42a1c47596a8f', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_result.json': '78f7d54fdd38ad4a4856b64fb232c6d5fc7e994dda49aade74e84140947d332e', 'scripts/run_original_native_J2_mean_endpoint_composition.py': '7581781d037b9d6b59a0cd86934d3521712e3588c17685d93bfbc2affb337400', 'src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py': 'b10b2600f1f1a3b4d0f387e66e412da668b666dcd8996a71305caaa071e3f65f', 'tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_source.json': '85421129e9d7ac7a3f9a6288f8bb48c184e738cc9101d7c4dad3d3aa274d45c0', 'tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_entered.json': '439bbacb222c1da88b84858147b6a6df4012d3a3c0ee4b9473dc1ced85514981', 'tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_result.json': '483ec74f99290a5f77fd79de0ce896f8a0fecea70dcc8cf8f493dbdc8f9fa1f2', 'scripts/run_original_native_J3_endpoint_transport.py': 'ebbe29d91389b4a531e548ee6f18cc4cd705406a801bbeadf5c89087acf4aebf', 'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': '65bad0c995cd6e508cb76a7986269dc5312fd0f3761ccaad14515a8863d24956', 'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_source.json': '5c7a83cafb2cca70ad6625520f1ab54884a5e566e1d215e61c409ba4618f6267', 'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_entered.json': '4a7fa9f818bf54b80b65dd3155be7cf1c386a9d591f1b587bad87897d4a425f3', 'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_result.json': '393cbc84869592cfb0d412225e81ad85054532e410799e463d98e106cf466c0d'}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def literal_configuration():
    tree = ast.parse((ROOT / KERNEL).read_text(encoding='utf-8'))
    assignments = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name):
                    assignments[target.id] = statement.value

    def resolve(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return resolve(assignments[node.id])
        if isinstance(node, (ast.List, ast.Tuple)):
            return [resolve(item) for item in node.elts]
        if isinstance(node, ast.Dict):
            return {resolve(key): resolve(value) for key, value in zip(node.keys, node.values)}
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return resolve(node.left) + resolve(node.right)
        raise ValueError('Nonliteral configuration: ' + ast.dump(node))

    return {name: resolve(assignments[name]) for name in
            ('SCHEMA', 'RECORD_NAMES', 'RECORD_HASHES', 'STAGE_FLAGS', 'FALSE_FLAGS')}


def source_path(name):
    candidate = ROOT / name
    return candidate if candidate.exists() else ROOT / 'experiments/root_st073' / name


def read_record(name, expected):
    path = source_path(name)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == expected, 'Changed raw source: ' + name)
    decoded = gzip.decompress(raw) if name.endswith('.gz') else raw
    return dict(source_text=decoded.decode('utf-8'), raw_sha256=actual,
        decoded_sha256=hashlib.sha256(decoded).hexdigest(), raw_authenticated_before_decode=True,
        is_gzip=name.endswith('.gz'), source_path=path.relative_to(ROOT).as_posix())


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def validate(result, packet):
    require(result.get('schema') == packet['schema_expected']
        and result.get('component') == COMPONENT
        and result.get('source_only') is True and type(result.get('passed')) is bool,
        'New source consumer schema/scope changed')
    require(result.get('ordinary_Z_orders') == [0, 1]
        and result.get('ordinary_Z_not_temporal_n') is True, 'Ordinary Z scope changed')
    for flag in packet['false_flags']:
        require(result.get(flag) is False, 'Unadmitted scope changed: ' + flag)
    require(all(type(result.get(flag)) is bool for flag in packet['stage_flags']), 'Missing local flags')
    require(type(result.get('unresolved_gates')) is list, 'Missing local prerequisite gates')
    require(type(result.get('remaining_construction_gates')) is list
        and bool(result['remaining_construction_gates']), 'Missing open construction gates')
    if result['passed']:
        require(all(result[flag] for flag in packet['stage_flags'])
            and not result['unresolved_gates'], 'Passed without local prerequisites')
    else:
        require(bool(result['unresolved_gates']), 'Negative source result without local gates')
    require(result['source_inspection_audit']['new_source_consumer_calls'] == 1
        and result['source_inspection_audit']['accepted_ancestor_operator_calls'] == 0,
        'Source consumer call audit changed')


def paths():
    return tuple(CK / ('new_endpoint_transport_v2_' + kind + '.json')
                 for kind in ('source', 'entered', 'result'))


def diagnostic_path():
    return CK / 'new_endpoint_transport_v2_diagnostics.jsonl'


def setup():
    config = literal_configuration()
    parent_record = read_record(PARENT, PARENT_SHA)
    parent = json.loads(parent_record['source_text'])
    require(parent['schema'] == 'original_native_J3_endpoint_transport_failure_stage_v1'
        and parent['component_result']['passed'] is False
        and parent['accepted_identity'] is False and parent['full_goal_complete'] is False,
        'Frozen failure parent scope changed')
    require(parent['remaining_B_identity_nodes'] == list(range(523, 536))
        and parent['native_IN_reduced_source_identity'] ==
        'I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)',
        'Inherited accepted B or I_N scope changed')
    inherited = parent['source_pins']
    require(len(inherited) == INHERITED_SOURCE_COUNT, 'Parent source inventory changed')
    require(parent['preserved_failed_artifact_hashes_union'] == FAILED_ATTEMPT_HASHES
        and len(FAILED_ATTEMPT_HASHES) == 20, 'Preserved failed-artifact union changed')
    for name, expected in inherited.items():
        require(sha(source_path(name)) == expected, 'Changed inherited source: ' + name)
    require(PARENT in config['RECORD_NAMES'], 'Failure parent missing from new kernel')
    expected_hashes = config['RECORD_HASHES']
    require(set(expected_hashes) == set(config['RECORD_NAMES'])
        and expected_hashes[PARENT] == PARENT_SHA, 'Changed new record inventory')
    for name, expected in expected_hashes.items():
        if name in inherited:
            require(inherited[name] == expected, 'Conflicting inherited record: ' + name)
        else:
            relative = source_path(name).relative_to(ROOT).as_posix()
            committed = subprocess.check_output(['git', 'cat-file', 'blob', 'HEAD:' + relative], cwd=ROOT)
            require(hashlib.sha256(committed).hexdigest() == expected,
                'Uncommitted source dependency: ' + name)
    records = {name: parent_record if name == PARENT else read_record(name, expected)
               for name, expected in expected_hashes.items()}
    for name, expected in FAILED_ATTEMPT_HASHES.items():
        require(sha(ROOT / name) == expected, 'Preserved failed attempt changed: ' + name)
    input_digest = digest(records)
    packet = dict(schema='native_J3_endpoint_transport_inputs_v2', component=COMPONENT,
        input_digest=input_digest, record_hashes=expected_hashes,
        decoded_record_hashes={name: record['decoded_sha256'] for name, record in records.items()},
        inherited_source_hashes=inherited,
        preserved_failed_attempt_hashes=FAILED_ATTEMPT_HASHES,
        implementation_hashes={name: sha(ROOT / name) for name in (RUNNER, KERNEL)},
        schema_expected=config['SCHEMA'], stage_flags=config['STAGE_FLAGS'],
        false_flags=config['FALSE_FLAGS'], ancestor_operators_replayed=False,
        kernel_imported_during_setup=False, numeric_values_computed=False, full_goal_complete=False,
        diagnostics_policy='append-only JSONL flush and fsync before each stage callback returns',
        diagnostics_checkpoint_additional_to_source_entered_result=True)
    return records, packet


def prepare():
    source, entered, target = paths()
    require(not any(path.exists() for path in (source, entered, target, diagnostic_path())),
        'Successor attempt already exists')
    _, packet = setup()
    write_new(source, packet)
    return dict(passed=True, prepared_only=True, source_sha256=sha(source),
        input_digest=packet['input_digest'], mathematical_consumer_calls=0)


def run():
    records, packet = setup()
    source, entered, target = paths()
    journal = diagnostic_path()
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists() and not journal.exists(),
        'Frozen input changed or successor attempt already entered')
    write_new(entered, dict(component=COMPONENT, source_sha256=sha(source),
        input_digest=packet['input_digest'], implementation_hashes=packet['implementation_hashes'],
        kernel_sha256=packet['implementation_hashes'][KERNEL],
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    diagnostics = []

    def persist_diagnostic(stage, payload):
        require(isinstance(stage, str) and bool(stage), 'Invalid diagnostic stage')
        # JSON conversion snapshots the payload before any later stage mutation.
        frozen = json.loads(json.dumps(payload, sort_keys=True, ensure_ascii=True, allow_nan=False))
        event = dict(schema='native_J3_v2_stage_diagnostic_v1', index=len(diagnostics),
            stage=stage, input_digest=packet['input_digest'], payload=frozen)
        raw = (json.dumps(event, sort_keys=True, ensure_ascii=True, allow_nan=False) + '\n').encode('utf-8')
        with journal.open('xb' if not diagnostics else 'ab') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        diagnostics.append(event)

    started, caught = time.perf_counter(), False
    try:
        spec = importlib.util.spec_from_file_location('new_native_J3_endpoint_transport_v2', ROOT / KERNEL)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = module.derive(records, expected_input_digest=packet['input_digest'],
            diagnostic_callback=persist_diagnostic)
    except Exception as error:
        caught = True
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback=''.join(traceback.format_exception(error)))
    result.update(execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(source), entry_checkpoint_sha256=sha(entered),
        runner_caught_exception=caught, runner_full_goal_complete=False,
        stage_diagnostics=diagnostics,
        diagnostic_checkpoint_path=journal.relative_to(ROOT).as_posix(),
        diagnostic_checkpoint_sha256=sha(journal) if journal.exists() else None,
        diagnostic_events_saved_before_validation=len(diagnostics))
    write_new(target, result)
    require(digest(records) == packet['input_digest'], 'Input mutated during successor consumer')
    if not caught:
        validate(result, packet)
    return dict(passed=result['passed'], result_sha256=sha(target), seconds=result['execution_seconds'],
        error=result.get('error'), unresolved_gates=result.get('unresolved_gates', []),
        diagnostic_events=len(diagnostics))


def saved():
    records, packet = setup()
    source, entered, target = paths()
    result, entry = json.loads(target.read_bytes()), json.loads(entered.read_bytes())
    require(json.loads(source.read_bytes()) == packet
        and result['source_checkpoint_sha256'] == sha(source)
        and entry['source_sha256'] == sha(source)
        and result['entry_checkpoint_sha256'] == sha(entered)
        and entry['implementation_hashes'] == packet['implementation_hashes'], 'Saved checkpoint changed')
    events = result['stage_diagnostics']
    journal = diagnostic_path()
    if events:
        require(result['diagnostic_checkpoint_sha256'] == sha(journal), 'Diagnostic journal changed')
        require([json.loads(line) for line in journal.read_text(encoding='utf-8').splitlines()] == events,
            'Saved diagnostic payloads changed')
    else:
        require(not journal.exists() and result['diagnostic_checkpoint_sha256'] is None,
            'Unexpected diagnostic journal')
    if result['runner_caught_exception']:
        require(result['passed'] is False and result['failed_attempt_preserved'] is True,
            'Malformed saved successor failure')
    else:
        validate(result, packet)
    return dict(passed=True, saved_only=True, component_passed=result['passed'], full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_endpoint_transport', 'endpoint_transport', 'inspect_endpoint_transport'))
    action = parser.parse_args().action
    print(json.dumps(prepare() if action == 'prepare_endpoint_transport' else
                     run() if action == 'endpoint_transport' else saved()))
