"""One-entry binding of saved actual core source sectors to collar inlet."""
from __future__ import annotations
import gzip, hashlib, importlib.util, json, sys, traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]; BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_core_inlet_v1'; DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
G3='experiments/root_st073/n1_core_g3_word0_v1_receipt.json.gz'
CORE='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz'
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz'
PHASE='experiments/root_st073/n1_actual_collar_phase_v1_receipt.json.gz'
CAUCHY='experiments/root_st073/n1_word1_z2_first_collar_v1_receipt.json.gz'
TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_actual_core_inlet_bind.py'
FIXTURE='src/openai_ns_reconstruction/n1_actual_core_inlet_fixture.py'
RUNNER='scripts/run_n1_actual_core_inlet_v1_once.py'
REQUIRED=(G3,CORE,FORCING,PHASE,CAUCHY,TAIL,LEAF,FIXTURE,RUNNER)
FIXED={G3:'a109cd2b88f3f6740ed4569cb9407d527e98a43287bc56609b41ff08c1052fc8',CORE:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',PHASE:'47c2caaf3b8c7006984445681b39a1f01cf6f1bd4788384602bfa36ef3dd2b11',CAUCHY:'89e4266ceaa89aa69c8ac41ac680b3f8fc92c6e6ee13e7b10f616b825b15cce9',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
def sha(b): return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v: raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f: f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes(); brief=json.loads(b); raw=SOURCE.read_bytes(); ck=brief['source_checkpoint']
    require(brief['component']==COMPONENT and brief['phase']==COMPONENT+'_pending' and brief['goal_active'] is True and brief['full_goal_complete'] is False and brief['one_entry_only'] is True and brief['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(raw),'source checkpoint')
    manifest=json.loads(raw);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for row in manifest['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(raw),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints')
    x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
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
    g3=json.loads(gzip.decompress(raw[G3]));core=json.loads(gzip.decompress(raw[CORE]));forcing=json.loads(gzip.decompress(raw[FORCING]));phase=json.loads(gzip.decompress(raw[PHASE]));cauchy=json.loads(gzip.decompress(raw[CAUCHY]))
    for receipt in (g3,core,forcing,phase,cauchy):require(receipt['result']['status']=='returned','input not returned')
    a=g3['result']['output'];b=core['result']['output'];f=forcing['result']['output'];p=phase['result']['output'];old=cauchy['result']['output'];family=b['source_family']
    require(a['source_family']==f['source_family']==p['source_family']==old['source_family']==family,'source family mismatch')
    require(b['prior_g3_result_sha256']==g3['frozen_byte_sha256']['result'],'g3 result binding')
    require(f['same_source_radial_and_core_receipts_verified'] is True and f['source_logs']==p['source_logs'],'same collar source')
    c=MPIntervalContext();c.prec=1800
    for key,other in (('F0_of_Z','F0_anchor'),('Lambda','Lambda')):
        require(iv(c,a['source_logs'][key])._mpi_==iv(c,b['source_logs'][key])._mpi_==iv(c,p['source_logs'][other])._mpi_,'different source atom '+key)
    g=a['W0_Z3_green_enclosure'];require(g['coordinate']=='core endpoint rho=4, X=2' and g['ordinary_Z_order']==3 and g['cell_count']==8,'g3 endpoint')
    require(b['result']['coordinate']=='rho=4, core exit, real Z in [-1,1]' and b['independent_source_matrix_fixture']['passed'] and b['independent_word1_green_fixture']['passed'],'core history')
    low={int(k):{row:[{**t,'coefficient':iv(c,t['coefficient'])} for t in terms] for row,terms in rows.items()} for k,rows in b['global_low_source'].items()}
    gfields={field:[{**t,'coefficient_interval':iv(c,t['coefficient_interval'])} for t in terms] for field,terms in g['fields'].items()}
    wfields={field:[{**t,'coefficient':iv(c,t['coefficient'])} for t in terms] for field,terms in b['result']['W1_Z2_core_exit'].items()}
    leaf=module(LEAF,'n1_actual_core_inlet_bind');fixture=module(FIXTURE,'n1_actual_core_inlet_fixture')
    proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str(proof))
    data=leaf.build(c,low,gfields,wfields)
    return {'source_family':family,'source_logs':p['source_logs'],'sc':p['sc'],'domain':p['domain'],'fixture':proof,'inlet':enc(data),
      'input_receipts':{x:hashes[x] for x in (G3,CORE,FORCING,PHASE,CAUCHY,TAIL)},
      'core_g3_is_endpoint_enclosure_not_last_cell_point_value':True,'low_orders_are_full_core_source_enclosures':True,
      'core_W1_history_integrated_once':True,'source_factors_formal':True,'cauchy_alternative_receipt':CAUCHY,
      'allword_tail_reference':TAIL,'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
