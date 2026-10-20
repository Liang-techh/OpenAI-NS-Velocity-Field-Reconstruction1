"""Once-only native Rc-to-Rp signed pressure source consumer; no ancestor execution."""
from pathlib import Path
from datetime import datetime, timezone
from inspect import signature
import argparse
import gzip
import hashlib
import importlib.util
import json
import time
import subprocess
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUNNER = 'scripts/run_original_native_Rc_Rp_signed_pressure_source.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_Rc_Rp_signed_pressure_source.py'
CK = ROOT / 'tmp/original_native_Rc_Rp_signed_pressure_source'
PARENT_REPORT = 'lei_ren_part1_paper_compliant_current_original_native_Rm_m_Z_source_lift.json'
PARENT_SHA = '2e2654a8be84f2f92ae149f3ed09ad3ef85679180024626138aa5f613211e254'
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'natural_physical_join_defects_zero_proved',
)
STAGE_FLAGS = (
    'native_Rc_incoming_source_bound', 'native_compact_quiet_source_continuation_proved',
    'signed_pressure_cocycle_proved', 'same_source_internal_boundary_cancellation_proved',
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def source_path(name):
    path = ROOT / name
    return path if path.exists() else ROOT / 'experiments/root_st073' / name


def paths():
    return tuple(CK / ('new_signed_pressure_' + kind + '.json') for kind in ('source', 'entered', 'result'))


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def load_kernel():
    spec = importlib.util.spec_from_file_location('new_native_Rc_Rp_signed_pressure_declarations', ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_record(name, expected):
    path = source_path(name)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == expected, 'changed raw source: ' + name)
    # The fixed raw digest is authenticated before text or gzip decoding.
    decoded = gzip.decompress(raw) if name.endswith('.gz') else raw
    return dict(source_text=decoded.decode('utf-8'), raw_sha256=actual,
        decoded_sha256=hashlib.sha256(decoded).hexdigest(), raw_authenticated_before_decode=True,
        is_gzip=name.endswith('.gz'), source_path=path.relative_to(ROOT).as_posix())


def setup():
    kernel = load_kernel()
    require(tuple(kernel.STAGE_FLAGS) == STAGE_FLAGS and tuple(kernel.FALSE_FLAGS) == FALSE_FLAGS, 'kernel scope mismatch')
    require(kernel.RECORD_HASHES.get(PARENT_REPORT) == PARENT_SHA, 'current source parent pin missing')
    records = {name: read_record(name, expected) for name, expected in kernel.RECORD_HASHES.items()}
    for name, expected in kernel.DECODED_RECORD_HASHES.items():
        require(records[name]['decoded_sha256'] == expected, 'decoded source changed: ' + name)
    parent = json.loads(records[PARENT_REPORT]['source_text'])
    require(parent['component_result']['passed'] is True
        and parent['component_result']['Rm_physical_Qy_Z1_function_identity_proved'] is True,
        'current native source parent not accepted')
    inherited = parent['source_pins']
    require(isinstance(inherited, dict), 'inherited source pins missing')
    for name, expected in inherited.items():
        require(sha(source_path(name)) == expected, 'changed inherited dependency: ' + name)
    for name, expected in kernel.RECORD_HASHES.items():
        if name != PARENT_REPORT:
            if name in inherited:
                require(inherited[name] == expected, 'changed inherited source owner: ' + name)
            else:
                committed = subprocess.check_output(['git', 'cat-file', 'blob',
                    'HEAD:' + source_path(name).relative_to(ROOT).as_posix()], cwd=ROOT)
                require(hashlib.sha256(committed).hexdigest() == expected,
                    'historical source report is not committed: ' + name)
    input_digest = digest(records)
    values = dict(source_records=records, expected_input_digest=input_digest)
    signature(kernel.derive).bind(**values)
    packet = dict(schema='native_Rc_Rp_signed_pressure_source_inputs_v1', component='signed_pressure',
        input_digest=input_digest, record_hashes=dict(kernel.RECORD_HASHES),
        decoded_record_hashes=dict(kernel.DECODED_RECORD_HASHES), inherited_source_hashes=inherited,
        implementation_hashes={name: sha(ROOT / name) for name in (RUNNER, KERNEL)},
        schema_expected=kernel.SCHEMA, accepted_mathematical_operators_replayed=False,
        old_definition_graphs_rebuilt=False, numeric_values_computed=False, full_goal_complete=False)
    return kernel, values, packet


def prepare():
    _, _, packet = setup()
    source, entered, target = paths()
    require(not any(path.exists() for path in (source, entered, target)), 'attempt already prepared')
    write_new(source, packet)
    return dict(passed=True, source_sha256=sha(source), input_digest=packet['input_digest'])


def validate_result(result, packet):
    require(isinstance(result, dict) and result.get('schema') == packet['schema_expected']
        and type(result.get('passed')) is bool and result.get('source_only') is True,
        'Rc-to-Rp signed pressure source schema/scope changed')
    require(result.get('ordinary_Z_orders') == [0, 1, 2, 3], 'pressure ordinary Z scope changed')
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False, 'unadmitted global scope: ' + flag)
    require(all(type(result.get(flag)) is bool for flag in STAGE_FLAGS), 'non-boolean source status')
    gates = result.get('unresolved_gates')
    require(isinstance(gates, list) and all(isinstance(item, str) for item in gates), 'invalid gate list')
    if result['passed']:
        require(all(result.get(flag) is True for flag in STAGE_FLAGS) and gates == [],
            'Rc-to-Rp signed pressure source passed without exact premises')
    else:
        require(result.get('signed_pressure_cocycle_proved') is False and bool(gates),
            'negative Rc-to-Rp signed pressure source needs failed identity and gates')


def run():
    kernel, values, packet = setup()
    source, entered, target = paths()
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists(), 'source changed or attempt already entered')
    signature(kernel.derive).bind(**values)
    write_new(entered, dict(component='signed_pressure', source_sha256=sha(source),
        input_digest=packet['input_digest'], kernel_sha256=packet['implementation_hashes'][KERNEL],
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result, caught = time.perf_counter(), None, False
    try:
        result = kernel.derive(**values)
    except Exception as error:
        caught = True
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered), runner_full_goal_complete=False, runner_caught_exception=caught)
    write_new(target, result)
    require(digest(values['source_records']) == packet['input_digest'], 'source input mutated')
    if caught:
        require(result.get('failed_attempt_preserved') is True and isinstance(result.get('traceback'), str),
            'malformed runner exception envelope')
    else:
        require(result.get('failed_attempt_preserved') is not True, 'consumer cannot self-authorize validation bypass')
        validate_result(result, packet)
    return dict(passed=result['passed'], result_sha256=sha(target), seconds=result['execution_seconds'],
        error=result.get('error'), unresolved_gates=result.get('unresolved_gates', []))


def saved():
    _, _, packet = setup()
    source, entered, target = paths()
    entry, result = json.loads(entered.read_bytes()), json.loads(target.read_bytes())
    require(json.loads(source.read_bytes()) == packet and entry['source_sha256'] == sha(source)
        and result['source_checkpoint_sha256'] == sha(source)
        and result['entry_checkpoint_sha256'] == sha(entered)
        and entry['input_digest'] == packet['input_digest']
        and entry['kernel_sha256'] == packet['implementation_hashes'][KERNEL], 'saved attempt changed')
    if result.get('runner_caught_exception') is True:
        require(result.get('failed_attempt_preserved') is True and result.get('passed') is False
            and isinstance(result.get('error_type'), str) and isinstance(result.get('traceback'), str),
            'malformed saved runner exception envelope')
    else:
        require(result.get('runner_caught_exception') is False
            and result.get('failed_attempt_preserved') is not True, 'invalid normal-return origin')
        validate_result(result, packet)
    return packet, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_signed_pressure', 'signed_pressure', 'inspect_signed_pressure'))
    action = parser.parse_args().action
    if action == 'prepare_signed_pressure':
        output = prepare()
    elif action == 'signed_pressure':
        output = run()
    else:
        _, result = saved()
        output = dict(passed=True, saved_only=True, component_passed=result['passed'], full_goal_complete=False)
    print(json.dumps(output))
