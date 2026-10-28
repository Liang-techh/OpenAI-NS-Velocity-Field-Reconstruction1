"""Execute the frozen local log-Rp source bridge once; never replay ancestors."""
from pathlib import Path
import hashlib, importlib.util, json, sys, traceback
R=Path(__file__).resolve().parents[1]
B=R.parent.parent/'RESUME_BRIEF_2026_10_16.json'
L=R/'tmp/original_native_logRp_source_bridge_v1'
SOURCE,CAPSULE,ENTRY,RESULT=(L/n for n in ('source.json','capsule.json','entered.json','result.json'))
CODE=R/'src/openai_ns_reconstruction/original_native_logRp_source_bridge.py'
COMPONENT='original_native_logRp_source_bridge_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save_new(p,v):
    raw=(json.dumps(v,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode('utf-8')
    with p.open('xb') as f:f.write(raw)
    return hashlib.sha256(raw).hexdigest()
def main():
    if ENTRY.exists() or RESULT.exists():raise RuntimeError('entered or saved: no reentry')
    source=json.loads(SOURCE.read_bytes());brief=json.loads(B.read_bytes())
    if brief.get('component')!=COMPONENT or brief.get('phase')!='prepared_logRp_bridge' or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False or brief.get('do_not_reenter') is not True:raise RuntimeError('pending brief mismatch')
    if brief['logRp_bridge_source_checkpoint']['sha256']!=sha(SOURCE):raise RuntimeError('checkpoint mismatch')
    for f in source['files']:
        if sha(Path(f['path']))!=f['sha256']:raise RuntimeError('frozen file mismatch '+f['path'])
    e=save_new(ENTRY,{'component':COMPONENT,'source_sha256':sha(SOURCE),'capsule_sha256':sha(CAPSULE),'entered':True,'consumer_calls_started':0,'ancestor_math_calls':0,'runtime_model_metadata':None})
    try:
        spec=importlib.util.spec_from_file_location('_logRp_source_bridge',CODE)
        if spec is None or spec.loader is None:raise RuntimeError('cannot load new kernel')
        mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
        result=mod.load_logrp_bridge_once(R)
        out={'component':COMPONENT,'status':'returned','source_sha256':sha(SOURCE),'entry_sha256':e,'capsule_sha256':sha(CAPSULE),'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None,'result':result}
    except Exception as exc:
        out={'component':COMPONENT,'status':'failed','source_sha256':sha(SOURCE),'entry_sha256':e,'capsule_sha256':sha(CAPSULE),'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None,'exception':type(exc).__name__+': '+str(exc),'traceback':traceback.format_exc()}
    h=save_new(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':h,'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
