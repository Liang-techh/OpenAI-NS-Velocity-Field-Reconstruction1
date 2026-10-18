"""Enter only the new Qy outer calculus; consume frozen primitive data."""
from __future__ import annotations

import argparse
import ast
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
PRIOR_RUNNER = 'scripts/run_original_Rm_physical_Qy_Z_source_ports.py'
PRIOR_RUNNER_SHA = '8a1a01aed6ecb765b1bc34d35645a5547f0134211774da27d378a22deef08074'
PARENT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_source_ports.json'
PARENT_SHA = 'f69ddf50c1f67c61199e447fe1dad9efaae8f0b7fb0e19babf4e7b460b60150e'
RUNNER = 'scripts/run_original_Rm_physical_Qy_Z_calculus.py'
KERNEL = 'src/openai_ns_reconstruction/original_Rm_physical_qy_Z_calculus.py'
KERNEL_SHA = '3475d8e213094bfce0a90461e421ce0352a938e8ca420217b27fdc1a89b2f82d'
CHECKPOINTS = ROOT/'tmp/original_Rm_physical_Qy_Z_calculus'
SOURCE = CHECKPOINTS/'new_source.json'
ENTRY = CHECKPOINTS/'new_operator_entered.json'
RESULT = CHECKPOINTS/'new_outer_result.json'
REPORT = ROOT/'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def write_new(path, value):
    require(not path.exists(), 'artifact already exists: '+str(path))
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


def setup():
    require(sha(ROOT/PRIOR_RUNNER) == PRIOR_RUNNER_SHA and sha(ROOT/PARENT) == PARENT_SHA, 'frozen parent changed')
    require(sha(ROOT/KERNEL) == KERNEL_SHA, 'new calculus source changed')
    for name in (RUNNER, KERNEL):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    old = module('frozen_Qy_Z_source_data', ROOT/PRIOR_RUNNER)
    _, nodes, inputs, prior_packet = old.setup()
    require(json.loads(old.SOURCE.read_bytes()) == prior_packet, 'frozen prior foundation changed')
    values, dependencies = old.actual_inputs('outer', inputs)
    parent = json.loads((ROOT/PARENT).read_bytes())
    require(parent['actual_Rm_physical_primitive_Z_jets'] == values['primitive_Z_source_receipt'], 'parent primitive differs')
    require(parent['actual_Rm_physical_Qy_Z_derivative_provenance']['outer_calculus_proved'] is False, 'old outer status changed')
    checkpoints = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(old.CHECKPOINTS.glob('*.json'))}
    require(len(checkpoints) == 7, 'frozen prior checkpoint set changed')
    require(parent['physical_Qy_Z_provenance_checkpoint_hashes'] == {Path(p).name: h for p, h in checkpoints.items()},
        'parent/checkpoint binding changed')
    packet = dict(schema='actual_Rm_physical_Qy_Z_calculus_inputs_v1', parent_report=PARENT, parent_sha256=PARENT_SHA,
        graph_digest=old.GRAPH_HASH, graph_node_count=len(nodes), input_digest=digest(values),
        implementation_hashes={RUNNER:sha(ROOT/RUNNER), KERNEL:KERNEL_SHA, PRIOR_RUNNER:PRIOR_RUNNER_SHA},
        frozen_checkpoint_hashes=checkpoints, primitive_dependency_hashes=dependencies,
        primitive_operator_replayed=False, accepted_mathematical_operator_replayed=False, full_goal_complete=False)
    return old, parent, nodes, values, packet


def prepare():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), 'new calculus attempt already prepared')
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE))


def run():
    old, _, nodes, values, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet, 'new foundation changed')
    require(not ENTRY.exists() and not RESULT.exists(), 'new outer already entered')
    kernel = module('new_Qy_Z_actual_calculus', ROOT/KERNEL)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), input_digest=packet['input_digest'], kernel_sha256=KERNEL_SHA,
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started = time.perf_counter()
    result = None
    try:
        result = kernel.derive(nodes, **values)
        require(digest(nodes) == old.GRAPH_HASH and digest(values) == packet['input_digest'], 'input mutated')
        old.validate_result('outer', result, values)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY), runner_full_goal_complete=False)
    write_new(RESULT, result)
    return dict(passed=result['passed'], outer_calculus_proved=result.get('outer_calculus_proved'),
        sha256=sha(RESULT), seconds=result['execution_seconds'], error=result.get('error'))


def saved():
    old, parent, _, values, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet, 'saved source changed')
    result, entry = json.loads(RESULT.read_bytes()), json.loads(ENTRY.read_bytes())
    require(result['source_checkpoint_sha256'] == entry['source_sha256'] == sha(SOURCE)
        and result['entry_checkpoint_sha256'] == sha(ENTRY) and entry['kernel_sha256'] == KERNEL_SHA
        and entry['input_digest'] == packet['input_digest'], 'saved entry/result changed')
    if not result.get('failed_attempt_preserved'):
        old.validate_result('outer', result, values)
    else:
        require(result['passed'] is False and isinstance(result.get('traceback'), str), 'malformed failure')
    return parent, packet, result


def build_report(parent, packet, result):
    report = copy.deepcopy(parent)
    report.update(stage='actual_Rm_physical_Qy_ordinary_Z_calculus_after_Z_atom_fix', parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA, actual_Rm_physical_Qy_Z_calculus_inputs=packet,
        preserved_first_Qy_Z_outer_result=parent['actual_Rm_physical_Qy_Z_derivative_provenance'],
        actual_Rm_physical_Qy_Z_derivative_provenance=result,
        physical_Qy_Z_calculus_checkpoint_hashes={p.name:sha(p) for p in sorted(CHECKPOINTS.glob('*.json'))})
    report['input_hashes'][Path(PARENT).name] = PARENT_SHA
    report['implementation_hashes'].update(packet['implementation_hashes'])
    report['scope'].update(physical_Q_y_Z1_value_seam_proved=result['passed'], current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report['execution_seconds']['new_outer_calculus_after_Z_atom_fix'] = result['execution_seconds']
    report['current_turn_physical_Qy_Z_execution_audit'].update(preserved_failed_primitive_calls=1,
        preserved_unresolved_outer_calls=1, new_outer_after_Z_atom_fix_calls=1, primitive_operator_replayed=False)
    return report


def write_report():
    parent, packet, result = saved()
    write_new(REPORT, build_report(parent, packet, result))
    return dict(saved=True, sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        outer_calculus_proved=result.get('outer_calculus_proved'), Qy_Z1_passed=result['passed'])


def inspect():
    parent, packet, result = saved()
    require(json.loads(REPORT.read_bytes()) == build_report(parent, packet, result), 'serialized report differs')
    return dict(passed=True, serialized_only=True, sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        outer_calculus_proved=result.get('outer_calculus_proved'), Qy_Z1_passed=result['passed'], full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'write_report', 'inspect'))
    print(json.dumps(globals()[parser.parse_args().action]()))
