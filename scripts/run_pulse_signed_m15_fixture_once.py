"""Enter the frozen, new synthetic signed-pulse fixture exactly once."""
from pathlib import Path
import hashlib, importlib.util, json, sys, traceback
R = Path(__file__).resolve().parents[1]
B = R.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
L = R / 'tmp/pulse_signed_m15_fixture_v1'
SOURCE, CAPSULE, ENTRY, RESULT = (L / n for n in ('source.json','capsule.json','entered.json','result.json'))
FIXTURE = R / 'scripts/check_pulse_signed_m15_correction_fixture.py'
COMPONENT = 'pulse_signed_m15_fixture_v1'
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
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
    if brief.get('component') != COMPONENT or brief.get('phase') != 'prepared_pulse_signed_m15_fixture' or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False or brief.get('do_not_reenter') is not True:
        raise RuntimeError('pending brief mismatch')
    if brief['pulse_signed_m15_source_checkpoint']['sha256'] != sha(SOURCE):
        raise RuntimeError('source checkpoint mismatch')
    for item in source['files']:
        if sha(Path(item['path'])) != item['sha256']:
            raise RuntimeError('frozen file changed: ' + item['path'])
    entry = save_new(ENTRY, {'component': COMPONENT, 'source_sha256': sha(SOURCE), 'capsule_sha256': sha(CAPSULE), 'entered': True, 'consumer_calls_started': 0, 'ancestor_math_calls': 0, 'runtime_model_metadata': None})
    try:
        sys.path.insert(0, str(R / 'src'))
        spec = importlib.util.spec_from_file_location('_pulse_signed_m15_fixture', FIXTURE)
        if spec is None or spec.loader is None:
            raise RuntimeError('fixture load failed')
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        result = mod.run_fixture()
        out = {'component': COMPONENT, 'status': 'returned', 'source_sha256': sha(SOURCE), 'entry_sha256': entry, 'capsule_sha256': sha(CAPSULE), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'result': result}
    except Exception as exc:
        out = {'component': COMPONENT, 'status': 'failed', 'source_sha256': sha(SOURCE), 'entry_sha256': entry, 'capsule_sha256': sha(CAPSULE), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'exception': type(exc).__name__ + ': ' + str(exc), 'traceback': traceback.format_exc()}
    digest = save_new(RESULT, out)
    print(json.dumps({'component': COMPONENT, 'status': out['status'], 'result_sha256': digest, 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None}))
    return 0 if out['status'] == 'returned' else 1
if __name__ == '__main__':
    raise SystemExit(main())