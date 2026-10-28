"""Single-entry runner for historical frame gap saved roles."""
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

root=Path(__file__).resolve().parents[1]
dir=root/'tmp/historical_frame_gap_deferred_roles_v1'
brief=Path(r'C:\Users\z5242\Documents\Codex\2026-09-13\before-doing-any-project-work-verify\RESUME_BRIEF_2026_10_16.json')
state=json.loads(brief.read_bytes())
if state.get('phase')!='historical_frame_gap_deferred_roles_pending' or state.get('component')!='historical_frame_gap_deferred_roles_v1':
    raise RuntimeError('continuation brief is not in the frozen pending phase')
source=json.loads((dir/'source.json').read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name,digest in source['sha256'].items():
    if sha(root/name)!=digest:raise RuntimeError('frozen input changed: '+name)
with (dir/'entered.json').open('x',encoding='utf-8') as f:
    json.dump({'component':state['component'],'entered':True,'source_sha256':sha(dir/'source.json')},f,sort_keys=True)
try:
    path=root/'src/openai_ns_reconstruction/historical_frame_gap_deferred_roles.py'
    spec=importlib.util.spec_from_file_location('historical_frame_gap_deferred_roles_leaf',path)
    leaf=importlib.util.module_from_spec(spec);sys.modules[spec.name]=leaf;spec.loader.exec_module(leaf)
    result={'status':'returned','result':leaf.recover(root),'consumer_calls':1,'runtime_model_metadata':None}
except Exception:
    result={'status':'failed','traceback':traceback.format_exc(),'consumer_calls':1,'runtime_model_metadata':None}
with (dir/'result.json').open('x',encoding='utf-8') as f:json.dump(result,f,sort_keys=True)
print(json.dumps({'status':result['status'],'result_status':result.get('result',{}).get('status')}))
