"""One-entry source-AST and symbolic collar pressure identity check."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
COMPONENT = 'collar_pressure_regular_singular_v6'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / x for x in ('source.json','entered.json','result.json'))
FORWARD = 'src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
CANONICAL = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py'
LEAF = 'src/openai_ns_reconstruction/collar_pressure_regular_singular_proof_v6.py'
RUNNER = 'scripts/run_collar_pressure_regular_singular_v6_once.py'
REQUIRED = frozenset((FORWARD,CANONICAL,LEAF,RUNNER))

def sha(b): return hashlib.sha256(b).hexdigest()
def exclusive(p, value):
    raw=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with p.open('xb') as f:f.write(raw)
    return sha(raw)
def require(value, message):
    if not value:raise ValueError(message)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    braw=BRIEF.read_bytes();b=json.loads(braw)
    require(b.get('component')==COMPONENT and b.get('phase')==COMPONENT+'_pending'
        and b.get('goal_active') is True and b.get('full_goal_complete') is False
        and b.get('do_not_reenter') is True and b.get('one_entry_only') is True,'pending guard')
    cp=b['source_checkpoint'];require(Path(cp['path']).resolve()==SOURCE.resolve(),'checkpoint path')
    sraw=SOURCE.read_bytes();require(sha(sraw)==cp['sha256'],'checkpoint hash')
    source=json.loads(sraw)
    require(source.get('component')==COMPONENT and source.get('one_entry_only') is True,'source scope')
    hashes={}
    for item in source['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected input')
        actual=sha((ROOT/p).read_bytes());require(actual==item['sha256'],'changed input '+p);hashes[p]=actual
    require(set(hashes)==REQUIRED,'missing input')
    return hashes,sha(sraw),sha(braw)
def consume(hashes):
    raw={}
    for name in REQUIRED:
        data=(ROOT/name).read_bytes();require(sha(data)==hashes[name],'changed after entry '+name);raw[name]=data
    spec=importlib.util.spec_from_file_location('collar_pressure_regular_singular_v6_leaf',ROOT/LEAF)
    require(spec is not None and spec.loader is not None,'leaf import')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module.prove(raw[FORWARD].decode('utf-8'),raw[CANONICAL].decode('utf-8'))
def main():
    hashes,source_sha,brief_sha=preflight()
    entry_sha=exclusive(ENTRY,{'component':COMPONENT,'entered':True,'source_sha256':source_sha,
        'brief_sha256':brief_sha,'input_hashes':hashes,'consumer_calls':0})
    out={'component':COMPONENT,'source_sha256':source_sha,'entry_sha256':entry_sha,
        'status':'failed','consumer_calls':0}
    try:
        out['consumer_calls']=1;out['proof']=consume(hashes);out['status']='returned';out['error']=None
    except BaseException as exc:
        out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    result_sha=exclusive(RESULT,out)
    print(json.dumps({'component':COMPONENT,'status':out['status'],'source_sha256':source_sha,
        'entry_sha256':entry_sha,'result_sha256':result_sha}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
