"""One-entry no-reset first-collar join of saved core and local low words."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_collar_low_join_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
CORE='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';COLLAR='experiments/root_st073/n1_actual_collar_low_word_v1_receipt.json.gz';Z2='experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz'
CAUCHY='experiments/root_st073/n1_word1_z2_first_collar_v1_receipt.json.gz';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_core_collar_low_join_v2.py';FIXTURE='src/openai_ns_reconstruction/n1_core_collar_low_join_fixture_v2.py';RUNNER='scripts/run_n1_core_collar_low_join_v2_once.py'
FAILURE='experiments/root_st073/n1_core_collar_low_join_v1_failure.json.gz'
REQUIRED=(CORE,COLLAR,Z2,CAUCHY,TAIL,FAILURE,LEAF,FIXTURE,RUNNER)
FIXED={FAILURE:'2d5a78ebe5f8936dc712a69561b7edbe7e22e17e8df0b439a418eab9845b4f02',CORE:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',COLLAR:'f4ed142fd3c8ec79b919af92ad12f3e160b06a6c8b13a932bb4fc17df27df234',Z2:'fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5',CAUCHY:'89e4266ceaa89aa69c8ac41ac680b3f8fc92c6e6ee13e7b10f616b825b15cce9',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
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
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');hashes={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes());require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(src),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
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
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    core=json.loads(gzip.decompress(raw[CORE]));collar=json.loads(gzip.decompress(raw[COLLAR]));z2=json.loads(gzip.decompress(raw[Z2]));cauchy=json.loads(gzip.decompress(raw[CAUCHY]));selected=json.loads(raw[TAIL])
    for r in (core,collar,z2,cauchy):require(r['result']['status']=='returned','input not returned')
    a=core['result']['output'];b=collar['result']['output'];z=z2['result']['output'];old=cauchy['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(b['source_family']==z['source_family']==old['source_family']==selected['source_family']==family,'source family mismatch')
    require(b['source_logs']==z['source_logs']==atoms and b['sc']==z['sc'] and a['input_receipts'][COLLAR]==hashes[COLLAR],'same atom, phase, and parent receipt')
    require(a['saved_W0_outer_cells_reused'] and a['Z2_not_recomputed'] and a['allword_tail_added'] is False and b['core_W1_Z0_Z1_inlet_missing'] and b['result']['core_W1_inlet_included'] is False and b['allword_tail_added'] is False,'separate core/local inputs')
    require(a['result']['coordinate']=='rho=4, core exit, real Z in [-1,1]' and a['result']['ordinary_Z_orders']==[0,1] and a['result']['input_W0_ordinary_Z_orders']==[0,1,2] and a['result']['first_collar_added'] is False and b['domain']=='rho=4 exp(hb*s), 0<=s<=sc/4, real Z in [-1,1]' and b['result']['ordinary_Z_fixed_phase_and_R'] and b['result']['local_word_orders']==[0,1] and b['result']['source_power_basis']==['hb','F0_anchor','Pstar_squared','Lambda'],'same coordinate and ordinary-Z convention')
    require(b['prior_Z2_green_receipt']==Z2 and b['alternative_cauchy_receipt']==CAUCHY and b['allword_tail_receipt']==TAIL and z['allword_tail_added'] is False,'previous alternatives and tail')
    measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'],'geometry contract')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,b['sc']);phase=iv(c,b['phase']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and phase.a==0 and phase.b==(sc/4).b,'exact same phase')
    index={};coredata={};localdata={}
    for order in (0,1):
        fields=a['result'][f'W1_Z{order}_core_exit'];coredata[order]={};localdata[order]={}
        for field,terms in fields.items():
            group=[]
            for i,t in enumerate(terms):
                sid=f'join_core:{order}:{field}:{i}';index[sid]={'receipt':CORE,'path':['result','W1_Z'+str(order)+'_core_exit',field,i]}
                group.append({**t,'coefficient':iv(c,t['coefficient']),'input_id':sid})
            coredata[order][field]=group
        for field,terms in b['result']['local_signed_sectors'][str(order)].items():
            group=[]
            for i,t in enumerate(terms):
                sid=f'join_local:{order}:{field}:{i}';index[sid]={'receipt':COLLAR,'path':['result','local_signed_sectors',str(order),field,i],'parent_source_ids':t['origin']['source_ids']}
                group.append({**t,'coefficient':iv(c,t['coefficient']),'input_id':sid})
            localdata[order][field]=group
    leaf=module(LEAF,'n1_core_collar_low_join_v2');fixture=module(FIXTURE,'n1_core_collar_low_join_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    geom={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True}
    result=leaf.build(c,sc=sc,phase=phase,core=coredata,local=localdata,source_family=family,source_atoms=atoms,geometry_certificate=geom)
    for order in (0,1):
        for field,terms in localdata[order].items():
            saved=result['local_contribution'][order][field];require(len(saved)==len(terms),'local duplication or omission')
            for original,joined in zip(terms,saved):
                require(original['powers']==joined['powers'] and original['coefficient']._mpi_==joined['coefficient']._mpi_ and all(q in joined['origin']['source_ids'] for q in original['origin']['source_ids']) and original['input_id'] in joined['origin']['source_ids'] and original['origin']['error_classes']==joined['origin']['error_classes'] and original['origin']['uncertainty']==joined['origin']['uncertainty'],'local input changed')
            require(result['finite_sum'][order][field]==result['core_transport'][order][field]+saved,'finite sum is not exact concatenation')
    seen=set()
    for group in ('core_transport','local_contribution'):
        for fields in result[group].values():
            for terms in fields.values():
                for term in terms:seen.update(s for s in term['origin']['source_ids'] if s.startswith('join_'))
    require(seen==set(index),'lost or unexpected input provenance')
    return {'source_family':family,'source_logs':atoms,'sc':b['sc'],'phase':b['phase'],'domain':b['domain'],'fixture':proof,'input_index':index,'result':enc(result),
      'core_receipt':CORE,'local_receipt':COLLAR,'prior_Z2_green_receipt':Z2,'alternative_cauchy_receipt':CAUCHY,'allword_tail_receipt':TAIL,
      'whole_collar_uniform_enclosure_not_point_values':True,'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
