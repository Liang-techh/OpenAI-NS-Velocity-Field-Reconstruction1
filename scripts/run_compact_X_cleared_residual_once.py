"""Single-entry runner for the compact X cleared-residual leaf."""
import hashlib,importlib.util,json,sys,traceback
from pathlib import Path
root=Path(__file__).resolve().parents[1];directory=root/'tmp/compact_X_cleared_residual_v1'
brief=Path(r'C:\Users\z5242\Documents\Codex\2026-09-13\before-doing-any-project-work-verify\RESUME_BRIEF_2026_10_16.json')
state=json.loads(brief.read_bytes())
if state.get('phase')!='compact_X_cleared_residual_pending' or state.get('component')!='compact_X_cleared_residual_v1':raise RuntimeError('wrong pending state')
source=json.loads((directory/'source.json').read_bytes())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
for name,expected in source['sha256'].items():
    if sha(root/name)!=expected:raise RuntimeError('frozen input changed: '+name)
with (directory/'entered.json').open('x',encoding='utf-8') as f:json.dump({'component':state['component'],'entered':True,'source_sha256':sha(directory/'source.json')},f,sort_keys=True)
try:
    path=root/'src/openai_ns_reconstruction/compact_X_cleared_residual.py';spec=importlib.util.spec_from_file_location('compact_X_cleared_residual_new_leaf',path)
    leaf=importlib.util.module_from_spec(spec);sys.modules[spec.name]=leaf;spec.loader.exec_module(leaf)
    result={'status':'returned','result':leaf.recover(root),'consumer_calls':1,'runtime_model_metadata':None}
except Exception:result={'status':'failed','traceback':traceback.format_exc(),'consumer_calls':1,'runtime_model_metadata':None}
with (directory/'result.json').open('x',encoding='utf-8') as f:json.dump(result,f,sort_keys=True)
print(json.dumps({'status':result['status'],'leaf_status':result.get('result',{}).get('status'),'checks_passed':result.get('result',{}).get('checks_passed')}))
