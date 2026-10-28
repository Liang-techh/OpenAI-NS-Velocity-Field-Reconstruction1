"""Enter each frozen new pulse mathematical leaf at most once."""
from pathlib import Path
import hashlib, importlib.util, json, sys, traceback
R = Path(__file__).resolve().parents[1]
B = R.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
LEAVES = {
    'chart': ('pulse_exact_chart_fixture_v1', 'scripts/check_pulse_source_chart_exact_fixture.py'),
    'density': ('pulse_correlated_density_fixture_v1', 'scripts/check_pulse_correlated_density_delta_fixture.py'),
}
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def save_new(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False) + '\n').encode('utf-8')
    with path.open('xb') as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()
def main():
    if len(sys.argv) != 2 or sys.argv[1] not in LEAVES:
        raise RuntimeError('select exactly one frozen leaf')
    component, relative_fixture = LEAVES[sys.argv[1]]
    leaf = R / 'tmp' / component
    source, capsule, entry, result = (leaf / n for n in ('source.json','capsule.json','entered.json','result.json'))
    if entry.exists() or result.exists():
        raise RuntimeError('already entered; never reenter')
    frozen = json.loads(source.read_bytes())
    brief = json.loads(B.read_bytes())
    if brief.get('component') != component or brief.get('phase') != 'prepared_pulse_exact_leaf' or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False or brief.get('do_not_reenter') is not True:
        raise RuntimeError('pending brief mismatch')
    if brief['pulse_exact_leaf_source_checkpoint']['sha256'] != sha(source):
        raise RuntimeError('source checkpoint mismatch')
    for item in frozen['files']:
        if sha(Path(item['path'])) != item['sha256']:
            raise RuntimeError('frozen file changed: ' + item['path'])
    receipt = save_new(entry, {'component': component, 'source_sha256': sha(source), 'capsule_sha256': sha(capsule), 'entered': True, 'consumer_calls_started': 0, 'ancestor_math_calls': 0, 'runtime_model_metadata': None})
    try:
        sys.path.insert(0, str(R / 'src'))
        spec = importlib.util.spec_from_file_location('_' + component, R / relative_fixture)
        if spec is None or spec.loader is None:
            raise RuntimeError('fixture load failed')
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        value = mod.run_fixture()
        out = {'component': component, 'status': 'returned', 'source_sha256': sha(source), 'entry_sha256': receipt, 'capsule_sha256': sha(capsule), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'result': value}
    except Exception as exc:
        out = {'component': component, 'status': 'failed', 'source_sha256': sha(source), 'entry_sha256': receipt, 'capsule_sha256': sha(capsule), 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None, 'exception': type(exc).__name__ + ': ' + str(exc), 'traceback': traceback.format_exc()}
    digest = save_new(result, out)
    print(json.dumps({'component': component, 'status': out['status'], 'result_sha256': digest, 'consumer_calls': 1, 'ancestor_math_calls': 0, 'runtime_model_metadata': None}))
    return 0 if out['status'] == 'returned' else 1
if __name__ == '__main__':
    raise SystemExit(main())