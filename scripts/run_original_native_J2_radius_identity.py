"""Freeze, enter once, and save the new actual J2 radius C1 real exp/log identity."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import gzip
import hashlib
import importlib.util
import json
import subprocess
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUNNER = 'scripts/run_original_native_J2_radius_identity.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_J2_radius_identity.py'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_remaining_B_seam_sources.json'
PARENT_SHA = 'f765b451fcb962cbf59824db2d0f7e994b3de72c30c3ee4bf46ab9f01c0cc1e6'
CK = ROOT / 'tmp/original_native_J2_radius_identity'
COMPONENT = 'J2_radius_identity'
FAILED_ATTEMPT_HASHES = {'scripts/run_original_native_micro_macro_B_jump_source.py': '583d26e94a4ec25913094f497eb442df6ab46c7201bd1a9521cbf47098f909f0', 'src/openai_ns_reconstruction/original_native_micro_macro_B_jump_source.py': 'c934743d502afe0d247862acf302f637a6f81f5686807a8dc45efacdbc268d3d', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_source.json': 'e0d7782d6bfc37d11c75ae08e05f7fb30588783bc3a7df6e9ceac6721e6c9f58', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_entered.json': '6b2815660d2b0d55717cddb4d576b9e23ef781ad901a713127d42a1c47596a8f', 'tmp/original_native_micro_macro_B_jump_source/new_micro_macro_B_jump_result.json': '78f7d54fdd38ad4a4856b64fb232c6d5fc7e994dda49aade74e84140947d332e'}


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


def paths():
    return tuple(CK / ('new_J2_radius_identity_' + kind + '.json')
                 for kind in ('source', 'entered', 'result'))


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def setup():
    config = literal_configuration()
    parent_record = read_record(PARENT, PARENT_SHA)
    parent = json.loads(parent_record['source_text'])
    require(parent['component_result']['passed'] is True
        and parent['full_goal_complete'] is False, 'Unaccepted parent')
    inherited = parent['source_pins']
    require(len(inherited) == 1830, 'Parent source inventory changed')
    for name, expected in inherited.items():
        require(sha(source_path(name)) == expected, 'Changed inherited source: ' + name)
    require(PARENT in config['RECORD_NAMES'], 'Parent missing from new kernel records')
    expected_hashes = config['RECORD_HASHES']
    require(set(expected_hashes) == set(config['RECORD_NAMES'])
        and expected_hashes[PARENT] == PARENT_SHA, 'Changed new record inventory')
    for name, expected in expected_hashes.items():
        if name in inherited:
            require(inherited[name] == expected, 'Conflicting inherited record: ' + name)
        elif name != 'tmp/original_native_remaining_B_seam_sources/native_J2_saved_radius_DAG_gap_2026_10_21.json':
            relative = source_path(name).relative_to(ROOT).as_posix()
            committed = subprocess.check_output(['git', 'cat-file', 'blob', 'HEAD:' + relative], cwd=ROOT)
            require(hashlib.sha256(committed).hexdigest() == expected,
                'Uncommitted legacy source dependency: ' + name)
        else:
            import base64
            archived = base64.b64decode(parent['saved_gap_evidence']['raw_base64'], validate=True)
            require(hashlib.sha256(archived).hexdigest() == expected
                and archived == source_path(name).read_bytes(), 'Changed saved gap evidence')
    records = {name: parent_record if name == PARENT else read_record(name, expected)
               for name, expected in expected_hashes.items()}
    for name, expected in FAILED_ATTEMPT_HASHES.items():
        require(sha(ROOT / name) == expected, 'Preserved failed attempt changed: ' + name)
    input_digest = digest(records)
    packet = dict(schema='native_J2_radius_identity_inputs_v1', component=COMPONENT,
        input_digest=input_digest, record_hashes=expected_hashes,
        decoded_record_hashes={name: record['decoded_sha256'] for name, record in records.items()},
        inherited_source_hashes=inherited,
        preserved_failed_attempt_hashes=FAILED_ATTEMPT_HASHES,
        implementation_hashes={name: sha(ROOT / name) for name in (RUNNER, KERNEL)},
        schema_expected=config['SCHEMA'], stage_flags=config['STAGE_FLAGS'],
        false_flags=config['FALSE_FLAGS'], ancestor_operators_replayed=False,
        kernel_imported_during_setup=False, numeric_values_computed=False, full_goal_complete=False)
    return records, packet


def validate(result, packet):
    require(result.get('schema') == packet['schema_expected']
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


def prepare():
    source, entered, target = paths()
    require(not any(path.exists() for path in (source, entered, target)), 'Attempt already exists')
    _, packet = setup()
    write_new(source, packet)
    return dict(passed=True, prepared_only=True, source_sha256=sha(source),
        input_digest=packet['input_digest'], mathematical_consumer_calls=0)


def run():
    records, packet = setup()
    source, entered, target = paths()
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists(), 'Frozen input changed or attempt already entered')
    write_new(entered, dict(component=COMPONENT, source_sha256=sha(source),
        input_digest=packet['input_digest'], implementation_hashes=packet['implementation_hashes'],
        kernel_sha256=packet['implementation_hashes'][KERNEL],
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, caught = time.perf_counter(), False
    try:
        spec = importlib.util.spec_from_file_location('new_native_J2_radius_identity', ROOT / KERNEL)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = module.derive(records, expected_input_digest=packet['input_digest'])
    except Exception as error:
        caught = True
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback=''.join(traceback.format_exception(error)))
    result.update(execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(source), entry_checkpoint_sha256=sha(entered),
        runner_caught_exception=caught, runner_full_goal_complete=False)
    write_new(target, result)
    require(digest(records) == packet['input_digest'], 'Input mutated during new consumer')
    if not caught:
        validate(result, packet)
    return dict(passed=result['passed'], result_sha256=sha(target), seconds=result['execution_seconds'],
        error=result.get('error'), unresolved_gates=result.get('unresolved_gates', []))


def saved():
    records, packet = setup()
    source, entered, target = paths()
    result, entry = json.loads(target.read_bytes()), json.loads(entered.read_bytes())
    require(json.loads(source.read_bytes()) == packet
        and result['source_checkpoint_sha256'] == sha(source)
        and entry['source_sha256'] == sha(source)
        and result['entry_checkpoint_sha256'] == sha(entered)
        and entry['implementation_hashes'] == packet['implementation_hashes'], 'Saved checkpoint changed')
    if result['runner_caught_exception']:
        require(result['passed'] is False and result['failed_attempt_preserved'] is True,
            'Malformed saved failure')
    else:
        validate(result, packet)
    return dict(passed=True, saved_only=True, component_passed=result['passed'], full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_J2_radius_identity', 'J2_radius_identity', 'inspect_J2_radius_identity'))
    action = parser.parse_args().action
    print(json.dumps(prepare() if action == 'prepare_J2_radius_identity' else
                     run() if action == 'J2_radius_identity' else saved()))
