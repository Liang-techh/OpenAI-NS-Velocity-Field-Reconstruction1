"""Single-use literal-geometry bridge; source/entry/result are append-only."""
from pathlib import Path
import hashlib, importlib.util, json, sys, traceback
R = Path(__file__).resolve().parents[1]
BRIEF = R.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
LEAF = R / 'tmp/original_native_pulse_geometry_literal_bridge_v1'
SOURCE, CAPSULE, ENTRY, RESULT = (LEAF / n for n in ('source.json','capsule.json','entered.json','result.json'))
CODE = R / 'src/openai_ns_reconstruction/original_native_pulse_geometry_literal_bridge.py'
COMPONENT = 'original_native_pulse_geometry_literal_bridge_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_new(p, x):
    raw=(json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':'))+'\n').encode('utf-8')
    with p.open('xb') as out: out.write(raw)
    return hashlib.sha256(raw).hexdigest()
def main():
    if ENTRY.exists() or RESULT.exists(): raise RuntimeError('already entered; no reentry')
    source=json.loads(SOURCE.read_bytes()); brief=json.loads(BRIEF.read_bytes())
    if brief.get('component')!=COMPONENT or brief.get('phase')!='prepared_geometry_bridge' or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False or brief.get('do_not_reenter') is not True: raise RuntimeError('pending brief mismatch')
    if brief['geometry_bridge_source_checkpoint']['sha256']!=sha(SOURCE): raise RuntimeError('source checkpoint mismatch')
    for row in source['files']:
        if sha(Path(row['path']))!=row['sha256']: raise RuntimeError('frozen file changed: '+row['path'])
    entry_sha=write_new(ENTRY,{'component':COMPONENT,'source_sha256':sha(SOURCE),'capsule_sha256':sha(CAPSULE),'entered':True,'consumer_calls_started':0,'ancestor_math_calls':0,'runtime_model_metadata':None})
    try:
        spec=importlib.util.spec_from_file_location('_geometry_literal_bridge',CODE)
        if spec is None or spec.loader is None: raise RuntimeError('cannot load new kernel')
        module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
        result=module.load_geometry_bridge_once(R)
        complete={'component':COMPONENT,'status':'returned','source_sha256':sha(SOURCE),'entry_sha256':entry_sha,'capsule_sha256':sha(CAPSULE),'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None,'result':result}
    except Exception as exc:
        complete={'component':COMPONENT,'status':'failed','source_sha256':sha(SOURCE),'entry_sha256':entry_sha,'capsule_sha256':sha(CAPSULE),'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None,'exception':type(exc).__name__+': '+str(exc),'traceback':traceback.format_exc()}
    result_sha=write_new(RESULT,complete)
    print(json.dumps({'component':COMPONENT,'status':complete['status'],'result_sha256':result_sha,'consumer_calls':1,'ancestor_math_calls':0,'runtime_model_metadata':None}))
    return 0 if complete['status']=='returned' else 1
if __name__=='__main__': raise SystemExit(main())
