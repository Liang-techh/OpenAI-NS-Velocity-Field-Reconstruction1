"""One-entry actual same-source collar own-mean enclosure."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_collar_mean_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
PRIOR='experiments/root_st073/n1_actual_collar_phase_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_actual_collar_mean_enclosure.py';FIXTURE='src/openai_ns_reconstruction/n1_actual_collar_mean_fixture.py';RUNNER='scripts/run_n1_actual_collar_mean_once.py'
REQUIRED=(PRIOR,LEAF,FIXTURE,RUNNER)
FIXED={PRIOR:'47c2caaf3b8c7006984445681b39a1f01cf6f1bd4788384602bfa36ef3dd2b11'}
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
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected input')
        h=sha((ROOT/p).read_bytes());require(h==item['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'source missing');return hashes,sha(src),sha(b)
def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']);hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'directed exact endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    prior=json.loads(gzip.decompress(raw[PRIOR]));require(prior['result']['status']=='returned' and prior['result_sha256']=='593c9ac5a58c490828ae35139d2647552e412a3c426436951841ebaf5f03296d','accepted actual phase receipt')
    saved=prior['result']['output'];require(saved['same_source_and_analytic_bound_receipts_verified'] is True
      and saved['actual_collar_F_and_U_enclosed_under_admitted_analytic_assumptions'] is True
      and saved['separate_cauchy_error_retained'] is True and saved['actual_signed_point_value_computed'] is False
      and saved['fixture']['passed'] is True,'actual source and separate error premises')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,saved['sc']);phase=iv(c,saved['enclosure']['phase_interval'])
    family=saved['source_family'];atoms=saved['source_logs'];require(saved['enclosure']['source_family']==family and saved['enclosure']['source_atoms']==atoms,'same source atom binding')
    core={k:iv(c,saved['own_mean_core_inlet'][str(k)]) for k in range(6)}
    U={}
    for k in range(6):
        terms=saved['enclosure']['fields']['U'][str(k)];U[k]=[{**v,'coefficient':iv(c,v['coefficient'])} for v in terms]
        require(U[k][0]['powers']==[0,0,0,0] and U[k][0]['semantics']=='actual_core_interval_not_selected_value'
          and all(t['semantics']=='symmetric_certified_error_not_field_value' and t['powers'][0]>=1 for t in U[k][1:]),'actual U source/error separation')
    certificate={k:True for k in ('same_source_checked','core_mean_inlet_actual','U_whole_phase_enclosure','hb_Z_independent','rho4_inlet_history')}
    leaf=module(LEAF,'n1_actual_collar_mean_leaf');fixture=module(FIXTURE,'n1_actual_collar_mean_fixture')
    proof=fixture.run(leaf,c);require(proof['passed'] is True,'independent fixture failed: '+str([k for k,v in proof['checks'].items() if not v]))
    out=leaf.build(c,sc=sc,phase=phase,core_mean_Z=core,U_Z=U,source_family=family,source_atoms=atoms,admitted_source_contract=certificate)
    return {'source_family':family,'source_logs':atoms,'sc':enc(sc),'phase':enc(phase),'fixture':proof,'mean':enc(out),
      'same_source_and_inlet_receipt_verified':True,'actual_mean_enclosed_under_admitted_analytic_assumptions':True,
      'input_phase_uncertainty_propagated_separately':True,'prior_result_sha256':prior['result_sha256'],
      'actual_signed_point_value_computed':False,'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(hashes),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
