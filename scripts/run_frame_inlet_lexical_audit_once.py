"""Run one frozen frame-inlet lexical audit and persist its complete result."""
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / 'tmp/frame_inlet_lexical_audit_v1'
BRIEF = Path(r'C:\Users\z5242\Documents\Codex\2026-09-13\before-doing-any-project-work-verify\RESUME_BRIEF_2026_10_16.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


state = json.loads(BRIEF.read_bytes())
if state.get('phase') != 'frame_inlet_lexical_audit_pending' or state.get('component') != 'frame_inlet_lexical_audit_v1':
    raise RuntimeError('wrong pending state')
source = json.loads((DIRECTORY / 'source.json').read_bytes())
for name, expected in source['sha256'].items():
    if sha(ROOT / name) != expected:
        raise RuntimeError('frozen input changed: ' + name)
with (DIRECTORY / 'entered.json').open('x', encoding='utf-8') as stream:
    json.dump({'component': state['component'], 'entered': True,
               'source_sha256': sha(DIRECTORY / 'source.json')}, stream, sort_keys=True)
try:
    path = ROOT / 'src/openai_ns_reconstruction/frame_selected_inlet_scope.py'
    spec = importlib.util.spec_from_file_location('frame_inlet_lexical_audit_new_leaf', path)
    leaf = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = leaf
    spec.loader.exec_module(leaf)
    result = {'status': 'returned', 'result': leaf.recover(ROOT),
              'consumer_calls': 1, 'runtime_model_metadata': None}
except Exception:
    result = {'status': 'failed', 'traceback': traceback.format_exc(),
              'consumer_calls': 1, 'runtime_model_metadata': None}
with (DIRECTORY / 'result.json').open('x', encoding='utf-8') as stream:
    json.dump(result, stream, sort_keys=True)
print(json.dumps({'status': result['status'], 'audit_status': result.get('result', {}).get('status')}))
