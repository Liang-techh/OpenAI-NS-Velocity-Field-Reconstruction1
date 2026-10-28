"""Enter the frozen actual gap source-expression builder exactly once."""
from pathlib import Path
import hashlib, importlib.util, json, sys, traceback
R = Path(__file__).resolve().parents[1]
B = R.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
L = R / 'tmp/pulse_gap_actual_source_program_v1'
SOURCE, CAPSULE, ENTRY, RESULT = (L / n for n in ('source.json','capsule.json','entered.json','result.json'))
CODE = R / 'src/openai_ns_reconstruction/pulse_gap_actual_source_program.py'
COMPONENT = 'pulse_gap_actual_source_program_v1'
def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()
def save_new(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False) + '\n').encode('utf-8')
    with path.open('xb') as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()
def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('already entered; never reenter')
    source = json.loads(SOURCE.read_bytes())
    brief = json.loads(B.read_bytes())
    if brief.get('component') != COMPONENT or brief.get('phase') != 'prepared_pulse_gap_actual_source' or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False or brief.get('do_not_reenter') is not True:
        raise RuntimeError('pending brief mismatch')
    if brief['pulse_gap_actual_source_checkpoint']['sha256'] != sha(SOURCE):
        raise RuntimeError('source checkpoint mismatch')
    for item in source['files']:
        if sha(Path(item['path'])) != item['sha256']:
            raise RuntimeError('frozen file changed: ' + item['path'])
    receipt = save_new(ENTRY, {'component': COMPONENT, 'source_sha256': sha(SOURCE), 'capsule_sha256': sha(CAPSULE), 'entered': True, 'consumer_calls_started': 0, 'ancestor_math_calls': 0, 'runtime_model_metadata': None})
    try:
        spec = importlib.util.spec_from_file_location('_pulse_gap_actual_source_program', CODE)
        if spec is None or spec.loader is None:
            raise RuntimeError('source program load failed')
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        value = mod.build_actual_gap_source_program(R)
        out = {'component': COMPONENT, 'status': 'returned', 'source_sha256': sha(SOURCE), 'entry_sha256': receipt, 'capsule_sha256': sha(CAPSULE), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'result': value}
    except Exception as exc:
        out = {'component': COMPONENT, 'status': 'failed', 'source_sha256': sha(SOURCE), 'entry_sha256': receipt, 'capsule_sha256': sha(CAPSULE), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'exception': type(exc).__name__ + ': ' + str(exc), 'traceback': traceback.format_exc()}
    digest = save_new(RESULT, out)
    print(json.dumps({'component': COMPONENT, 'status': out['status'], 'result_sha256': digest, 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None}))
    return 0 if out['status'] == 'returned' else 1
if __name__ == '__main__':
    raise SystemExit(main())