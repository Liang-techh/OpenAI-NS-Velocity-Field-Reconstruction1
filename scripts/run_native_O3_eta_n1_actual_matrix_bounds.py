"""Single-use source-bound interval matrix consumer; never replay ancestors.

Root must freeze the code and source capsule, persist source.json and the
pending resume brief, then invoke this runner once. The entry checkpoint is
written exclusively before importing the new kernel. The complete returned
object, or the complete exception, is saved before any result review.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

R = Path(__file__).resolve().parents[1]
LEAF = R / 'tmp/native_O3_eta_n1_actual_matrix_bounds_v1'
SOURCE = LEAF / 'source.json'
CAPSULE = LEAF / 'capsule.json'
ENTRY = LEAF / 'entered.json'
RESULT = LEAF / 'result.json'
BRIEF = R.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
CODE = R / 'src/openai_ns_reconstruction/native_O3_eta_n1_actual_matrix_bounds.py'
COMPONENT = 'native_O3_eta_n1_actual_matrix_bounds'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, value):
    raw = (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':')) + '\n').encode('utf-8')
    with path.open('xb') as out:
        out.write(raw)
    return hashlib.sha256(raw).hexdigest()


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('this source consumer has already been entered')
    source = json.loads(SOURCE.read_bytes())
    brief = json.loads(BRIEF.read_bytes())
    if (brief.get('component') != COMPONENT or brief.get('phase') != 'prepared_unexecuted_actual_matrix_bounds'
            or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False
            or brief.get('do_not_reenter') is not True):
        raise RuntimeError('pending brief lifecycle mismatch')
    if brief.get('actual_matrix_source_checkpoint', {}).get('sha256') != sha(SOURCE):
        raise RuntimeError('pending brief/source hash mismatch')
    for item in source['files']:
        path = Path(item['path'])
        if sha(path) != item['sha256']:
            raise RuntimeError('frozen source changed: ' + str(path))
    if source['capsule_sha256'] != sha(CAPSULE):
        raise RuntimeError('capsule hash mismatch')
    entry = {'component': COMPONENT, 'source_sha256': sha(SOURCE),
             'capsule_sha256': sha(CAPSULE), 'entered': True, 'consumer_calls_started': 0,
             'ancestor_math_calls': 0, 'runtime_model_metadata': None}
    entry_sha = write_new(ENTRY, entry)
    try:
        spec = importlib.util.spec_from_file_location('_native_actual_matrix_bounds', CODE)
        if spec is None or spec.loader is None:
            raise RuntimeError('new kernel cannot be imported')
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        capsule = json.loads(CAPSULE.read_bytes())
        capsule['entered'] = True
        result = module.consume_frozen_source_bundle(capsule, interval_dps=280)
        complete = {'component': COMPONENT, 'status': 'returned', 'source_sha256': sha(SOURCE),
                    'entry_sha256': entry_sha, 'capsule_sha256': sha(CAPSULE),
                    'consumer_calls': 1, 'ancestor_math_calls': 0,
                    'runtime_model_metadata': None, 'result': result}
    except Exception as exc:
        complete = {'component': COMPONENT, 'status': 'failed', 'source_sha256': sha(SOURCE),
                    'entry_sha256': entry_sha, 'capsule_sha256': sha(CAPSULE),
                    'consumer_calls': 1, 'ancestor_math_calls': 0,
                    'runtime_model_metadata': None, 'exception': type(exc).__name__ + ': ' + str(exc),
                    'traceback': traceback.format_exc()}
    result_sha = write_new(RESULT, complete)
    print(json.dumps({'component': COMPONENT, 'status': complete['status'], 'result_sha256': result_sha,
                      'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None}))
    return 0 if complete['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
