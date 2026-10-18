"""Enter only the corrected leading-m source bridge; preserve first failure."""
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
OLD_RUNNER = 'scripts/run_original_Rm_m_ordinary_Z_sources.py'
OLD_RUNNER_SHA = '9fd9efc0f70ce149ae8126aaa7924e5e20e8d94f6f701bfe6954160f92fa40d1'
RUNNER = 'scripts/run_original_Rm_leading_m_Z_source_access.py'
KERNEL = 'src/openai_ns_reconstruction/original_Rm_leading_m_ordinary_Z_source_access.py'
KERNEL_SHA = 'afe78796fa3fed60918b62ec1d358eedeef7f9762b7dc60e06005b8a41539ed4'
FROZEN = {
    'src/openai_ns_reconstruction/original_Rm_leading_m_ordinary_Z_source.py':'f64309b284b3e26d4c66f51b3593c087377c68770b55d387c7314517cec589eb',
    'tmp/original_Rm_m_ordinary_Z_sources/new_leading_source.json':'746a4fbbf443026e0f81aa3edb07c156233b40d386788239d41b04da080d6799',
    'tmp/original_Rm_m_ordinary_Z_sources/new_leading_entered.json':'39a794ad369f4763c9f0202d8914765a08d2028efd739550912034d79b75d03d',
    'tmp/original_Rm_m_ordinary_Z_sources/new_leading_result.json':'d753f18d84f17f5ef83a4bf20fe654e3e7a8a6cd9a20e8d757f4c0473b6cdbd5',
}
CK = ROOT/'tmp/original_Rm_leading_m_Z_source_access'
SOURCE, ENTRY, RESULT = (CK/name for name in ('new_source.json','new_operator_entered.json','new_result.json'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def write_new(path, value):
    require(not path.exists(), 'artifact exists: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write('\n')


def setup():
    require(sha(ROOT/OLD_RUNNER) == OLD_RUNNER_SHA and sha(ROOT/KERNEL) == KERNEL_SHA, 'source changed')
    for name, expected in FROZEN.items():
        require(sha(ROOT/name) == expected, 'frozen failed source changed: '+name)
    for name in (RUNNER, KERNEL):
        raw = (ROOT/name).read_bytes()
        require(b'\r' not in raw and raw.endswith(b'\n') and not raw.endswith(b'\n\n'), 'LF/EOF changed')
        ast.parse(raw, filename=name)
    old = module('Rm_leading_Z_saved_source_inputs', ROOT/OLD_RUNNER)
    parent, nodes, inputs, prior = old.setup('leading')
    require(json.loads(old.SOURCES['leading'].read_bytes()) == prior, 'original failed foundation changed')
    values = inputs['leading']
    packet = dict(schema='Rm_leading_m_Z_source_access_inputs_v1', parent_sha256=old.PARENT_SHA,
        graph_digest=old.GRAPH_SHA, input_digest=digest(values), frozen_failure_hashes=FROZEN,
        inherited_source_record_hashes=prior['record_hashes']['leading'],
        implementation_hashes={RUNNER:sha(ROOT/RUNNER), KERNEL:KERNEL_SHA, OLD_RUNNER:OLD_RUNNER_SHA},
        accepted_mathematical_operators_replayed=False, full_goal_complete=False)
    return old, parent, nodes, values, packet


def prepare():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), 'new attempt already prepared')
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE))


def run():
    old, _, nodes, values, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet and not ENTRY.exists() and not RESULT.exists(), 'attempt already entered or changed')
    kernel = module('actual_Rm_leading_Z_access_source', ROOT/KERNEL)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_SHA,
        input_digest=digest(values), entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started, result = time.perf_counter(), None
    try:
        result = kernel.derive(nodes, **values)
        require(digest(nodes) == old.GRAPH_SHA and digest(values) == packet['input_digest'], 'source input mutated')
        old.validate_result('leading', result)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback=''.join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY), runner_full_goal_complete=False)
    write_new(RESULT, result)
    return dict(passed=result['passed'], sha256=sha(RESULT), seconds=result['execution_seconds'], error=result.get('error'),
        unresolved_gates=result.get('unresolved_gates', []))


def saved():
    old, parent, _, values, packet = setup()
    result, entry = json.loads(RESULT.read_bytes()), json.loads(ENTRY.read_bytes())
    require(json.loads(SOURCE.read_bytes()) == packet and result['source_checkpoint_sha256'] == entry['source_sha256'] == sha(SOURCE)
        and result['entry_checkpoint_sha256'] == sha(ENTRY) and entry['input_digest'] == packet['input_digest']
        and entry['kernel_sha256'] == KERNEL_SHA, 'saved source entry result changed')
    if not result.get('failed_attempt_preserved'):
        old.validate_result('leading', result)
    else:
        require(result['passed'] is False and isinstance(result.get('traceback'), str), 'malformed failed result')
    return parent, packet, result


def inspect():
    _, _, result = saved()
    return dict(passed=True, saved_only=True, leading_source_passed=result['passed'], result_sha256=sha(RESULT), full_goal_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','run','inspect'))
    print(json.dumps(globals()[parser.parse_args().action]()))
