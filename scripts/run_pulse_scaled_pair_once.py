"""Entry-frozen runner for the single new synthetic scaled-pair fixture."""
import hashlib
import json
import traceback
from pathlib import Path

root=Path(__file__).resolve().parents[1]
dir=root/'tmp/pulse_scaled_pair_v1'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
brief=Path(r'C:\Users\z5242\Documents\Codex\2026-09-13\before-doing-any-project-work-verify\RESUME_BRIEF_2026_10_16.json')
pending=json.loads(brief.read_text(encoding='utf-8'))
if pending.get('phase')!='pulse_scaled_pair_pending' or pending.get('component')!='pulse_scaled_pair_v1':
    raise RuntimeError('brief is not in the frozen pending phase')
source=json.loads((dir/'source.json').read_text(encoding='utf-8'))
for name,digest in source['sha256'].items():
    if sha(root/name)!=digest: raise RuntimeError('frozen input hash changed: '+name)
entry={'component':'pulse_scaled_pair_v1','entered':True,'source_sha256':sha(dir/'source.json')}
with (dir/'entered.json').open('x',encoding='utf-8') as f: json.dump(entry,f,sort_keys=True)
try:
    import fixture_pulse_scaled_pair_once as fixture
    value=fixture.run()
    result={'status':'returned','result':value,'consumer_invocations':1,'ancestor_math_calls':0,'runtime_model_metadata':None}
except Exception:
    result={'status':'failed','traceback':traceback.format_exc(),'consumer_invocations':1,'ancestor_math_calls':0,'runtime_model_metadata':None}
with (dir/'result.json').open('x',encoding='utf-8') as f: json.dump(result,f,sort_keys=True)
print(json.dumps({'status':result['status'],'all_passed':result.get('result',{}).get('all_passed')}))
