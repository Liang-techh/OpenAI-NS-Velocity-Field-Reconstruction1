"""Single-entry runner for the bounded inlet replacement free-name leaf."""
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

root = Path(__file__).resolve().parents[1]
directory = root / 'tmp/frame_replacement_free_names_v2'
brief = Path(r'C:\Users\z5242\Documents\Codex\2026-09-13\before-doing-any-project-work-verify\RESUME_BRIEF_2026_10_16.json')
state = json.loads(brief.read_bytes())
if state.get('phase') != 'frame_replacement_free_names_v2_pending' or state.get('component') != 'frame_replacement_free_names_v2':
    raise RuntimeError('wrong pending state')
source = json.loads((directory / 'source.json').read_bytes())
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
for name, expected in source['sha256'].items():
    if sha(root / name) != expected:
        raise RuntimeError('frozen input changed: ' + name)
with (directory / 'entered.json').open('x', encoding='utf-8') as f:
    json.dump({'component': state['component'], 'entered': True,
               'source_sha256': sha(directory / 'source.json')}, f, sort_keys=True)
try:
    path = root / 'src/openai_ns_reconstruction/frame_replacement_free_names_v2.py'
    spec = importlib.util.spec_from_file_location('frame_replacement_free_names_new_leaf', path)
    leaf = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = leaf
    spec.loader.exec_module(leaf)
    result = {'status': 'returned', 'result': leaf.recover(root),
              'consumer_calls': 1, 'runtime_model_metadata': None}
except Exception:
    result = {'status': 'failed', 'traceback': traceback.format_exc(),
              'consumer_calls': 1, 'runtime_model_metadata': None}
with (directory / 'result.json').open('x', encoding='utf-8') as f:
    json.dump(result, f, sort_keys=True)
print(json.dumps({'status': result['status'], 'leaf_status': result.get('result', {}).get('status')}))
