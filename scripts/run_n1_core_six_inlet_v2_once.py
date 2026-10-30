"""One-entry, source-bound fixed-x=2 six-field baseline and packet assembly."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_six_inlet_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';Z2='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz';EXIT='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_g0_six_exit_v2.py';FIXTURE='src/openai_ns_reconstruction/n1_g0_six_exit_fixture_v2.py';RUNNER='scripts/run_n1_core_six_inlet_v2_once.py'
FAILURE='experiments/root_st073/n1_core_six_inlet_v1_failure.json.gz'
REQUIRED=(LOW,Z2,EXIT,LEAF,FIXTURE,RUNNER,FAILURE)
FIXED={FAILURE:'f46f37d28faa444c66002b77ad3bdb1304383fa63d4ae78b4a2665a084acfa11',LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',Z2:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',EXIT:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921'}
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
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def conv(c,x):
    if isinstance(x,dict):
        if 'lower_exact_mpf_tuple' in x and 'upper_exact_mpf_tuple' in x:return iv(c,x)
        return {k:conv(c,v) for k,v in x.items()}
    if isinstance(x,list):return [conv(c,v) for v in x]
    return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError(type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def consume(hashes):
    raw={}
    for path in REQUIRED:
        b=(ROOT/path).read_bytes();require(sha(b)==hashes[path],'changed after entry '+path);raw[path]=b
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_g0_six_exit_v2');fixture=module(FIXTURE,'n1_g0_six_exit_fixture_v2')
    proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    low=json.loads(gzip.decompress(raw[LOW]));z2=json.loads(gzip.decompress(raw[Z2]));exit=json.loads(gzip.decompress(raw[EXIT]))
    require(all(v['result']['status']=='returned' for v in (low,z2,exit)),'accepted receipt status')
    lp={v['path']:v['sha256'] for v in low['source']['files']};ep={v['path']:v['sha256'] for v in exit['source']['files']};require(lp[Z2]==hashes[Z2] and ep[LOW]==hashes[LOW],'authenticated receipt chain')
    a=low['result']['output'];z=z2['result']['output'];e=exit['result']['output'];require(a['source_family']==z['source_family']==e['source_family'],'same source family')
    require(iv(c,z['source_logs']['F0_of_Z'])._mpi_==iv(c,a['source_logs']['F0_anchor'])._mpi_ and iv(c,z['source_logs']['Lambda'])._mpi_==iv(c,a['source_logs']['Lambda'])._mpi_ and a['source_logs']==e['source_logs'],'same formal source logs')
    require(z['whole_core_intervals_not_point_values'] is True and z['positive_source_factors_not_materialized'] is True,'source enclosure semantics')
    sources=conv(c,z['global_low_source']);index={}
    for order in ('0','1'):
        for row in ('4','5','6'):
            terms=sources[order][row];expected={'4':[(0,0),(0,-1)],'5':[(1,-1)],'6':[(0,-1),(0,-2)]}[row]
            require([(v['F0_power'],v['Lambda_power']) for v in terms]==expected,'formal source powers')
            for i,v in enumerate(terms):
                sid=f'global_g:Z{order}:row{row}:term{i}';index[sid]={'receipt':Z2,'path':['global_low_source',order,row,i]};v['source_id']=sid;v['source_ids']=[sid];v['origin']={'source_ids':[sid],'error_classes':['shared_whole_core_source_enclosure'],'uncertainty':True};v['error_classes']=['shared_whole_core_source_enclosure'];v['uncertainty']=True;v['powers']=[0,v['F0_power'],0,v['Lambda_power']]
    outer=a['result']['outer_cells'];require(len(outer)==8,'eight saved source cells');cells=[{'t':iv(c,cell['t'])} for cell in outer]
    baseline=leaf.build_g0g(c,cells,low_source=sources)
    require(e['feedback_word_tail_included'] is False and e['full_six_component_inlet_computed'] is False,'saved exit limitations')
    packets={'E00_saved_enclosure':{'receipt':LOW,'path':['result','output','result','W1_Z0_core_exit'],'ordinary_Z1_path':['result','output','result','W1_Z1_core_exit']},
      'Dg_saved_fixed_exit':{'receipt':EXIT,'path':['result','output','saved_W0_first_correction_fixed_exit']},
      'raw_Q_finite':{'receipt':EXIT,'path':['result','output','raw_Q_fixed_exit']},
      'P_finite':{'receipt':EXIT,'path':['result','output','saved_potential_compensation_fixed_exit']},
      'inner_kernel_tail':{'receipt':EXIT,'path':['result','output','inner_W0_kernel_tail_fixed_exit']},
      'outer_kernel_tail':{'receipt':EXIT,'path':['result','output','outer_kernel_tail_fixed_exit']},
      'potential_kernel_tail':{'receipt':EXIT,'path':['result','output','potential_compensation_kernel_tail_fixed_exit']}}
    return {'source_family':a['source_family'],'source_logs':a['source_logs'],'source_index':index,'fixture':proof,'baseline_G0g_fixed_exit':enc(baseline),'existing_disjoint_packets':packets,'input_receipt_sha256':{p:hashes[p] for p in (LOW,Z2,EXIT)},'shared_global_source_reused_across_cells':True,'point_values_computed':False,'existing_packets_added_or_cancelled':False,'feedback_word_tail_included':False,'complete_six_component_inlet':False,'full_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
