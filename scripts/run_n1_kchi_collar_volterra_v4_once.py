"""Freeze and execute only the independent variable-Kchi collar synthetic fixture."""
from __future__ import annotations
import hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_kchi_collar_volterra_v4';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LEAF='src/openai_ns_reconstruction/n1_kchi_collar_volterra_v4.py';FIXTURE='src/openai_ns_reconstruction/n1_kchi_collar_volterra_fixture_v4.py';RUNNER='scripts/run_n1_kchi_collar_volterra_v4_once.py'
FAILURE='experiments/root_st073/n1_kchi_collar_volterra_v2_failure.json.gz'
PRIOR='experiments/root_st073/n1_kchi_collar_volterra_v3_receipt.json.gz'
REQUIRED=(LEAF,FIXTURE,RUNNER,FAILURE,PRIOR)
FIXED={FAILURE:'612517aaa5069577a48a274be74c297e5252a964fe697719c23d6814020b55d2',PRIOR:'4ec675c128283440bf41ea80070d00b77e1a72f5a5bb123bdf486c88ebf14c9a'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def consume(hashes):
    for path in REQUIRED:require(sha((ROOT/path).read_bytes())==hashes[path],'changed after entry '+path)
    c=MPIntervalContext();c.prec=200;leaf=module(LEAF,'n1_kchi_collar_volterra_v4');fixture=module(FIXTURE,'n1_kchi_collar_volterra_fixture_v4');proof=fixture.run_fixture(leaf,c)
    require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    return {'fixture':proof,'synthetic_only':True,'actual_source_data_consumed':False,'actual_collar_seam_identity_proved':False,'feedback_word_tail_included':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
