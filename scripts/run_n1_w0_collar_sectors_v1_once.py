"""One-entry actual first-collar W0 low-order four-cell enclosures."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_w0_collar_sectors_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
GREEN='experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz';FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';AXIAL='experiments/root_st073/n1_first_word_axial_m15_v1_receipt.json.gz';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_w0_collar_sectors_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_w0_collar_sectors_fixture_v1.py';RUNNER='scripts/run_n1_w0_collar_sectors_v1_once.py'
REQUIRED=(GREEN,FORCING,SWIRL,AXIAL,TAIL,LEAF,FIXTURE,RUNNER)
FIXED={GREEN:'fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',AXIAL:'02acdcb22d9b0f03b4201f502388814d54d49937124ebb2288e66b9706770ab5',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
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
    green=json.loads(gzip.decompress(raw[GREEN]));forcing=json.loads(gzip.decompress(raw[FORCING]));swirl=json.loads(gzip.decompress(raw[SWIRL]));axial=json.loads(gzip.decompress(raw[AXIAL]));selected=json.loads(raw[TAIL])
    for receipt in (green,forcing,swirl,axial):require(receipt['result']['status']=='returned','accepted receipt')
    a=green['result']['output'];f=forcing['result']['output'];w=swirl['result']['output'];u=axial['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(f['source_family']==w['source_family']==u['source_family']==selected['source_family']==family and f['source_logs']==w['source_logs']==u['source_logs']==atoms,'same source family and logs')
    require(a['sc']==w['sc']==u['sc'] and f['same_source_radial_and_core_receipts_verified'] and a['allword_tail_added'] is False and w['appended_to_cumulative_m15'] is False and u['appended_to_cumulative_m15'] is False,'same geometry and independent prior contributions')
    measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'],'geometry source')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,a['sc']);phase=iv(c,a['phase']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and phase.a==0 and phase.b==(sc/4).b,'same exact phase')
    index={};core={};g={}
    for order in (0,1):
        core[order]={};g[order]={}
        for field,terms in a['result']['core_W0_transport'][str(order)].items():
            core[order][field]=[]
            for i,t in enumerate(terms):
                sid=f'w0_core:{order}:{field}:{i}';index[sid]={'receipt':GREEN,'path':['result','core_W0_transport',str(order),field,i],'parent_source_ids':t['origin']['source_ids']}
                core[order][field].append({**t,'coefficient':iv(c,t['coefficient']),'input_id':sid})
        for row,terms in f['g'][str(order)].items():
            g[order][row]=[]
            for i,t in enumerate(terms):
                sid=f'w0_g:{order}:{row}:{i}';index[sid]={'receipt':FORCING,'path':['g',str(order),row,i],'parent_source_ids':t['origin']['source_ids']}
                g[order][row].append({**t,'coefficient':iv(c,t['coefficient']),'input_id':sid})
    leaf=module(LEAF,'n1_w0_collar_sectors_v1');fixture=module(FIXTURE,'n1_w0_collar_sectors_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True};result=leaf.build(c,sc=sc,core=core,g=g,source_family=family,source_atoms=atoms,certificate=cert)
    seen=set();classes=set()
    for fields in result['core_transport'].values():
        for terms in fields.values():
            for term in terms:seen.update(q for q in term['origin']['source_ids'] if q.startswith('w0_'));classes.update(term['origin']['error_classes'])
    for cells in result['local'].values():
        for fields in cells.values():
            for terms in fields.values():
                for term in terms:seen.update(q for q in term['origin']['source_ids'] if q.startswith('w0_'));classes.update(term['origin']['error_classes'])
    require(seen==set(index) and classes,'lost input provenance')
    return {'source_family':family,'source_logs':atoms,'sc':a['sc'],'domain':a['domain'],'fixture':proof,'input_index':index,'result':enc(result),'output_error_classes':sorted(classes),
      'saved_core_transport_receipt':GREEN,'source_forcing_receipt':FORCING,'prior_axial_receipt':AXIAL,'prior_swirl_receipt':SWIRL,'allword_tail_receipt':TAIL,
      'W0_green_is_new_local_cell_enclosure':True,'old_W1_Z2_not_recomputed':True,'allword_tail_added':False,'m15_computed':False,'appended_to_cumulative_m15':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
