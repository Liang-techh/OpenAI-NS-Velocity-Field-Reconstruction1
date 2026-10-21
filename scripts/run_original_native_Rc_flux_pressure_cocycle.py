"""Once-only actual Rc flux pressure cocycle; no ancestor execution."""
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
RUNNER = 'scripts/run_original_native_Rc_flux_pressure_cocycle.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_Rc_flux_pressure_cocycle.py'
CK = ROOT / 'tmp/original_native_Rc_flux_pressure_cocycle'
PARENT_REPORT = 'lei_ren_part1_paper_compliant_current_original_native_lower_physical_flux_identity.json'
PARENT_SHA = 'f1f98828579a76fc2e71f812b7da0f356620f73e5c707bb91edb1bf789792035'
FALSE_FLAGS = (
    'native_compact_m15_source_identity_proved',
    'native_compact_V_m_source_identity_proved',
    'native_active_m15_source_identity_proved',
    'terminal_difference_zero_proved', 'old_error_terms_removed',
    'same_density_certificate',
    'finite_N_history_identity',
    'physical_Q_y_Z1',
    'global_five_moment_closure',
    'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved',
    'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed',
    'numeric_values_computed',
    'old_operators_replayed',
    'full_goal_complete',
    'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
)
STAGE_FLAGS = (
    'actual_original_Rc_active_compact_flux_C1_identity_proved',
    'actual_native_Rc_cumulative_pressure_cocycle_bound',
    'actual_corrected_active_pressure_cocycle_bound',
    'internal_Rc_and_lower_fluxes_cancelled_by_source_identities',
    'signed_terminal_reduced_to_actual_density_integral_difference',
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
    return tuple(CK / ('new_Rc_cocycle_' + kind + '.json')
                 for kind in ('source', 'entered', 'result'))


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def load_kernel():
    spec = importlib.util.spec_from_file_location(
        'new_native_Rc_flux_pressure_cocycle_declarations', ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_record(name, expected):
    path = source_path(name)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == expected, 'changed raw source: ' + name)
    # Authenticate the fixed raw digest before any text or gzip decoding.
    decoded = gzip.decompress(raw) if name.endswith('.gz') else raw
    return dict(source_text=decoded.decode('utf-8'), raw_sha256=actual,
        decoded_sha256=hashlib.sha256(decoded).hexdigest(), raw_authenticated_before_decode=True,
        is_gzip=name.endswith('.gz'), source_path=path.relative_to(ROOT).as_posix())


def setup():
    kernel = load_kernel()
    require(tuple(kernel.STAGE_FLAGS) == STAGE_FLAGS
        and tuple(kernel.FALSE_FLAGS) == FALSE_FLAGS, 'kernel scope mismatch')
    require(kernel.RECORD_HASHES.get(PARENT_REPORT) == PARENT_SHA,
        'current lower flux source parent pin missing')
    records = {name: read_record(name, expected)
               for name, expected in kernel.RECORD_HASHES.items()}
    for name, expected in kernel.DECODED_RECORD_HASHES.items():
        require(records[name]['decoded_sha256'] == expected, 'decoded source changed: ' + name)
    parent = json.loads(records[PARENT_REPORT]['source_text'])
    require(parent['component_result']['passed'] is True
        and parent['component_result']['lower_physical_RQy_source_C1_identity_proved_on_owned_Z_cell'] is True
        and parent['component_result']['signed_terminal_Rc_flux_correlation_preserved'] is True,
        'current lower flux source parent not accepted')
    inherited = parent['source_pins']
    require(isinstance(inherited, dict), 'inherited source pins missing')
    for name, expected in inherited.items():
        require(sha(source_path(name)) == expected, 'changed inherited dependency: ' + name)
    for name, expected in kernel.RECORD_HASHES.items():
        if name != PARENT_REPORT:
            if name in inherited:
                require(inherited[name] == expected, 'changed inherited source owner: ' + name)
            else:
                committed = subprocess.check_output([
                    'git', 'cat-file', 'blob',
                    'HEAD:' + source_path(name).relative_to(ROOT).as_posix()], cwd=ROOT)
                require(hashlib.sha256(committed).hexdigest() == expected,
                    'historical source report is not committed: ' + name)
    input_digest = digest(records)
    values = dict(source_records=records, expected_input_digest=input_digest)
    signature(kernel.derive).bind(**values)
    packet = dict(schema='native_Rc_flux_pressure_cocycle_inputs_v1', component='Rc_cocycle',
        input_digest=input_digest, record_hashes=dict(kernel.RECORD_HASHES),
        decoded_record_hashes=dict(kernel.DECODED_RECORD_HASHES), inherited_source_hashes=inherited,
        implementation_hashes={name: sha(ROOT / name) for name in (RUNNER, KERNEL)},
        schema_expected=kernel.SCHEMA, accepted_mathematical_operators_replayed=False,
        old_definition_graphs_rebuilt=False, numeric_values_computed=False, full_goal_complete=False)
    return kernel, values, packet


def prepare():
    _, _, packet = setup()
    source, entered, target = paths()
    require(not any(path.exists() for path in (source, entered, target)),
        'attempt already prepared')
    write_new(source, packet)
    return dict(passed=True, source_sha256=sha(source), input_digest=packet['input_digest'])


def validate_result(result, packet):
    require(isinstance(result, dict) and result.get('schema') == packet['schema_expected']
        and type(result.get('passed')) is bool and result.get('source_only') is True,
        'Rc flux pressure cocycle source identity schema/scope changed')
    require(result.get('ordinary_Z_orders') == [0, 1],
        'Rc flux pressure cocycle ordinary Z scope changed')
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False, 'unadmitted global scope: ' + flag)
    require(all(type(result.get(flag)) is bool for flag in STAGE_FLAGS),
        'non-boolean Rc flux pressure cocycle status')
    gates = result.get('unresolved_gates')
    require(isinstance(gates, list) and all(isinstance(item, str) for item in gates),
        'invalid gate list')
    if result['passed']:
        require(all(result.get(flag) is True for flag in STAGE_FLAGS) and gates == [],
            'Rc flux pressure cocycle source identity passed without exact premises')
    else:
        require(result.get('internal_Rc_and_lower_fluxes_cancelled_by_source_identities') is False
            and bool(gates),
            'negative Rc flux pressure cocycle source identity needs failed identity and gates')


def run():
    kernel, values, packet = setup()
    source, entered, target = paths()
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists(), 'source changed or attempt already entered')
    signature(kernel.derive).bind(**values)
    write_new(entered, dict(component='Rc_cocycle', source_sha256=sha(source),
        input_digest=packet['input_digest'], kernel_sha256=packet['implementation_hashes'][KERNEL],
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result, caught = time.perf_counter(), None, False
    try:
        result = kernel.derive(**values)
    except Exception as error:
        caught = True
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback=''.join(traceback.format_exception(error)),
            rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started,
        source_checkpoint_sha256=sha(source), entry_checkpoint_sha256=sha(entered),
        runner_full_goal_complete=False, runner_caught_exception=caught)
    write_new(target, result)
    require(digest(values['source_records']) == packet['input_digest'], 'source input mutated')
    if caught:
        require(result.get('failed_attempt_preserved') is True
            and isinstance(result.get('traceback'), str), 'malformed runner exception envelope')
    else:
        require(result.get('failed_attempt_preserved') is not True,
            'consumer cannot self-authorize validation bypass')
        validate_result(result, packet)
    return dict(passed=result['passed'], result_sha256=sha(target),
        seconds=result['execution_seconds'], error=result.get('error'),
        unresolved_gates=result.get('unresolved_gates', []))


def saved():
    _, _, packet = setup()
    source, entered, target = paths()
    entry, result = json.loads(entered.read_bytes()), json.loads(target.read_bytes())
    require(json.loads(source.read_bytes()) == packet and entry['source_sha256'] == sha(source)
        and result['source_checkpoint_sha256'] == sha(source)
        and result['entry_checkpoint_sha256'] == sha(entered)
        and entry['input_digest'] == packet['input_digest']
        and entry['kernel_sha256'] == packet['implementation_hashes'][KERNEL],
        'saved attempt changed')
    if result.get('runner_caught_exception') is True:
        require(result.get('failed_attempt_preserved') is True and result.get('passed') is False
            and isinstance(result.get('error_type'), str)
            and isinstance(result.get('traceback'), str),
            'malformed saved runner exception envelope')
    else:
        require(result.get('runner_caught_exception') is False
            and result.get('failed_attempt_preserved') is not True,
            'invalid normal-return origin')
        validate_result(result, packet)
    return packet, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=(
        'prepare_Rc_cocycle', 'Rc_cocycle', 'inspect_Rc_cocycle'))
    action = parser.parse_args().action
    if action == 'prepare_Rc_cocycle':
        output = prepare()
    elif action == 'Rc_cocycle':
        output = run()
    else:
        _, result = saved()
        output = dict(passed=True, saved_only=True,
                      component_passed=result['passed'], full_goal_complete=False)
    print(json.dumps(output))
