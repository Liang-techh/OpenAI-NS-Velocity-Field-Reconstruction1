"""One-entry independent synthetic checks of source-factored collar leaves."""
from __future__ import annotations
import hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_collar_factored_source_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
REQUIRED=('src/openai_ns_reconstruction/n1_first_collar_forcing_jets.py',
 'src/openai_ns_reconstruction/n1_first_collar_factored_green.py',
 'src/openai_ns_reconstruction/n1_first_collar_forcing_fixture.py',
 'src/openai_ns_reconstruction/n1_first_collar_factored_fixture.py',
 'scripts/run_n1_collar_factored_source_once.py',
 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py',
 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py',
 'docs/CURRENT_N1_FIRST_COLLAR_SIGNED_GREEN_2026_10_30.md')
def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True
      and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected source')
        h=sha((ROOT/p).read_bytes());require(h==item['sha256'],'source changed: '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'source missing');return hashes,sha(src),sha(b)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);obj=importlib.util.module_from_spec(spec)
    sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:
        for p,h in hashes.items():require(sha((ROOT/p).read_bytes())==h,'post-entry input changed')
        force=module(REQUIRED[0],'collar_forcing_leaf');green=module(REQUIRED[1],'collar_factored_green_leaf')
        ff=module(REQUIRED[2],'collar_forcing_fixture');gf=module(REQUIRED[3],'collar_factored_fixture')
        c=MPIntervalContext();c.prec=300
        first=ff.run(force,c);require(first['passed'] is True,'forcing checks failed: '+str([k for k,v in first['checks'].items() if not v]))
        second=gf.run(green,c);require(second['passed'] is True,'Green checks failed: '+str([k for k,v in second['checks'].items() if not v]))
        out.update(status='returned',output={'forcing_fixture':first,'factored_green_fixture':second,
          'actual_source_called':False,'same_source_binding_proved':False,'hb_inside_forcing_grids_preserved':False,
          'actual_signed_collar_values_computed':False,'allword_tail_added':False,'m15_computed':False},error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
