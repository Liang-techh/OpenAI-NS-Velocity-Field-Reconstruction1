"""Once-only admission overlay for the frozen native finite-N candidate DAG."""
from pathlib import Path
from datetime import datetime, timezone
from inspect import signature
import argparse
import hashlib
import importlib.util
import json
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUNNER = 'scripts/run_original_native_C2_patch_finite_N_admission.py'
KERNEL = 'src/openai_ns_reconstruction/original_native_C2_patch_finite_N_admission.py'
CK = ROOT / 'tmp/original_native_C2_patch_finite_N_admission'
INPUTS = {
    'saved_attempt': ('lei_ren_part1_paper_compliant_current_original_native_C2_patch_finite_N_source.json',
                      '5194f631aba2eebd4946003e1af52181ac8c050b92d636c36b6bdd2058245b5a'),
    'latest_report': ('lei_ren_part1_paper_compliant_current_original_native_C2_patch_composition.json',
                      'f426f396b3b9f55aea1bf1b6c07cd7678addae7644093b97babb110a47619136'),
    'provider_report': ('lei_ren_part1_paper_compliant_current_original_finite_N_provider_source_binding.json',
                        '1f23b6188abffb6d99cf00fb1274a58a2fba1100a326334a1a798d0b920f6c41'),
}
FROZEN = {
    'src/openai_ns_reconstruction/original_native_C2_patch_finite_N_source.py':
        '0f85613e17a59bcb1d5825acfafb720a8dc1f4580993cb06d88c444dfcc62bd2',
    'scripts/run_original_native_C2_patch_finite_N_source.py':
        '7669dd9a78a6eeb3c9a72f5cf372f52a64cad46411f5df767389e4e5ee3c8429',
    'tmp/original_native_C2_patch_finite_N_source/new_finite_N_source_source.json':
        '681681cc1e00a1043dd38c40a43be5667ed6d7cb18b1e089ac1341f05880186f',
    'tmp/original_native_C2_patch_finite_N_source/new_finite_N_source_entered.json':
        '6dc59bb9013397868b19397b89924c952dc7c81e9aa24516036198190448f6cc',
    'tmp/original_native_C2_patch_finite_N_source/new_finite_N_source_result.json':
        'ac8b4e100124d47e7c2696d21945edae05726052a4d7bba1a47563f666950379',
}
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete',
)
PATCH_FLAGS = (
    'native_C2_patch_finite_N_fields_defined', 'native_C2_C5_patch_provider_identity_proved',
    'native_C2_C5_patch_density_identity_proved', 'native_C2_patch_correction_inlet_identity_proved',
    'native_C2_C5_patch_history_identity_proved', 'patch_finite_N_source_identity',
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def paths():
    return tuple(CK / ('new_admission_' + kind + '.json') for kind in ('source', 'entered', 'result'))


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def load_kernel():
    spec = importlib.util.spec_from_file_location('new_native_finite_N_admission_declarations', ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def setup():
    values = {}
    for label, (name, expected) in INPUTS.items():
        path = ROOT / 'experiments/root_st073' / name
        require(sha(path) == expected, 'changed saved report: ' + name)
        values[label] = json.loads(path.read_bytes())
    for name, expected in FROZEN.items():
        require(sha(ROOT / name) == expected, 'changed frozen attempt: ' + name)
    require(values['saved_attempt']['component_result'] == json.loads((ROOT / next(
        name for name in FROZEN if name.endswith('_result.json'))).read_bytes()),
        'published negative and frozen result differ')
    for label in tuple(INPUTS):
        prefix = {'saved_attempt': 'attempt', 'latest_report': 'latest', 'provider_report': 'provider'}[label]
        values['expected_' + prefix + '_digest'] = digest(values[label])
    kernel = load_kernel()
    signature(kernel.derive).bind(**values)
    packet = dict(schema='native_C2_patch_finite_N_admission_inputs_v1',
        component='admission', input_reports={name: expected for name, expected in INPUTS.values()},
        frozen_attempt_hashes=FROZEN, input_digest=digest(values),
        expected_canonical_digests={name: value for name, value in values.items() if name.startswith('expected_')},
        implementation_hashes={name: sha(ROOT / name) for name in (RUNNER, KERNEL)},
        schema_expected=kernel.SCHEMA, old_definition_graphs_rebuilt=False,
        accepted_mathematical_operators_replayed=False, numeric_values_computed=False,
        full_goal_complete=False)
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
        'admission schema/scope changed')
    require(result.get('ordinary_Z_orders') == [0, 1, 2], 'ordinary Z scope changed')
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False, 'unadmitted global scope: ' + flag)
    require(all(type(result.get(flag)) is bool for flag in PATCH_FLAGS), 'non-boolean patch status')
    gates = result.get('unresolved_gates')
    require(isinstance(gates, list) and all(isinstance(item, str) for item in gates), 'invalid gate list')
    if result['passed']:
        require(all(result.get(flag) is True for flag in PATCH_FLAGS)
            and gates == [], 'admission passed without exact patch premises')
    else:
        require(result.get('patch_finite_N_source_identity') is False and bool(gates),
            'negative result needs a failed identity and explicit gates')


def run():
    kernel, values, packet = setup()
    source, entered, target = paths()
    require(json.loads(source.read_bytes()) == packet and not entered.exists()
        and not target.exists(), 'source changed or entered attempt')
    signature(kernel.derive).bind(**values)
    write_new(entered, dict(component='admission', source_sha256=sha(source), input_digest=packet['input_digest'],
        kernel_sha256=packet['implementation_hashes'][KERNEL], entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result = time.perf_counter(), None
    try:
        result = kernel.derive(**values)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered), runner_full_goal_complete=False)
    write_new(target, result)
    require(digest(values) == packet['input_digest'], 'saved metadata input mutated')
    if not result.get('failed_attempt_preserved'):
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
    if not result.get('failed_attempt_preserved'):
        validate_result(result, packet)
    return packet, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare_admission', 'admission', 'inspect_admission'))
    action = parser.parse_args().action
    if action == 'prepare_admission':
        output = prepare()
    elif action == 'admission':
        output = run()
    else:
        _, result = saved()
        output = dict(passed=True, saved_only=True, component_passed=result['passed'], full_goal_complete=False)
    print(json.dumps(output))
