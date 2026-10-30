"""One-entry actual-source first-collar finite Green word enclosure."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_collar_green_v3';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
INLET='experiments/root_st073/n1_actual_core_inlet_v1_receipt.json.gz'
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz'
PHASE='experiments/root_st073/n1_actual_collar_phase_v1_receipt.json.gz'
CAUCHY='experiments/root_st073/n1_word1_z2_first_collar_v1_receipt.json.gz'
TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_actual_collar_green_compact_v3.py'
FIXTURE='src/openai_ns_reconstruction/n1_actual_collar_green_fixture_v3.py'
RUNNER='scripts/run_n1_actual_collar_green_v3_once.py'
FAILURE='experiments/root_st073/n1_actual_collar_green_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_actual_collar_green_v2_failure.json.gz'
REQUIRED=(INLET,FORCING,PHASE,CAUCHY,TAIL,FAILURE,FAILURE2,LEAF,FIXTURE,RUNNER)
FIXED={FAILURE2:'3877a42aaf085f190eb0e44f148c56b6bd24fd861c3e99945e26a68b2957ca88',FAILURE:'eacf3a66da3466ad761061eba39acd1d8e5da004ee072abcd133a2131dc73aa0',INLET:'722c60ce802d5b1adcf8e4735406a791c5c58896ee6d27f0da040dc0e4c42d34',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',PHASE:'47c2caaf3b8c7006984445681b39a1f01cf6f1bd4788384602bfa36ef3dd2b11',CAUCHY:'89e4266ceaa89aa69c8ac41ac680b3f8fc92c6e6ee13e7b10f616b825b15cce9',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
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
    for row in manifest['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(src),sha(b)
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
    inlet=json.loads(gzip.decompress(raw[INLET]));forcing=json.loads(gzip.decompress(raw[FORCING]));phase=json.loads(gzip.decompress(raw[PHASE]));cauchy=json.loads(gzip.decompress(raw[CAUCHY]));selected=json.loads(raw[TAIL])
    for r in (inlet,forcing,phase,cauchy):require(r['result']['status']=='returned','input not returned')
    a=inlet['result']['output'];f=forcing['result']['output'];p=phase['result']['output'];old=cauchy['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(f['source_family']==p['source_family']==old['source_family']==selected['source_family']==family,'source family mismatch')
    require(atoms==p['source_logs']==f['source_logs'] and a['sc']==p['sc'],'source atoms and sc')
    require(a['core_W1_history_integrated_once'] and a['core_g3_is_endpoint_enclosure_not_last_cell_point_value'] and f['same_source_radial_and_core_receipts_verified'] and p['same_source_and_analytic_bound_receipts_verified'] and p['separate_cauchy_error_retained'],'bound actual inputs')
    require(a['allword_tail_added'] is False and f['allword_tail_added'] is False and p['allword_tail_added'] is False,'tail kept separate')
    require(a['input_receipts'][FORCING]==hashes[FORCING] and a['input_receipts'][PHASE]==hashes[PHASE] and a['input_receipts'][CAUCHY]==hashes[CAUCHY] and a['input_receipts'][TAIL]==hashes[TAIL],'same receipt chain')
    measure=selected['actual_cutoff_complete_support_measure'];words=selected['actual_selected_first_word_enclosure'];tail=selected['actual_selected_raw_word_tail']
    require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'] and words['passed'] and tail['passed'],'geometry/tail contract')
    require(old['same_source_binding']['actual_first_collar_uniform_norm'] and old['same_source_binding']['Cauchy_gap_verified'],'prior Cauchy bound')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,a['sc']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'exact sc identity')
    phase_iv=c.make_mpf((c.mpf(0)._mpi_[0],(sc/4).b._mpi_[1]))
    source_index={};core=a['inlet'];core_index=core['source_index']
    def coreterms(terms):
        out=[]
        for t in terms:
            require(t['origin']['source_ids'] and all(sid in core_index for sid in t['origin']['source_ids']) and (not t['origin']['uncertainty'] or t['origin']['error_classes']),'missing core provenance')
            out.append({**t,'coefficient':iv(c,t['coefficient'])})
        return out
    core0={int(k):{field:coreterms(terms) for field,terms in rows.items()} for k,rows in core['core_W0_Z'].items()}
    core1={field:coreterms(terms) for field,terms in core['core_W1_Z2'].items()}
    def localterms(terms,path):
        out=[]
        for i,t in enumerate(terms):
            sid='local:'+':'.join(map(str,path))+':'+str(i);origin=t['origin']
            require(all(q in f['source_index'] for q in origin['source_ids']), 'unmapped forcing parent IDs '+str(path)+':'+str(i))
            require(not origin['uncertainty'] or origin['error_classes'], 'unclassified forcing uncertainty '+str(path)+':'+str(i))
            source_index[sid]={'receipt':FORCING,'path':path+[i],'parent_source_ids':origin['source_ids']}
            out.append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':origin['uncertainty'],'error_classes':origin['error_classes'],'path_count':origin.get('path_count',1),'source_count':origin.get('source_count',1)}})
        return out
    g={int(k):{row:localterms(terms,['g',k,row]) for row,terms in rows.items()} for k,rows in f['g'].items()}
    matrix={name:{int(k):{key:localterms(terms,['matrix',name,k,key]) for key,terms in rows.items()} for k,rows in orders.items()} for name,orders in f['matrix'].items()}
    leaf=module(LEAF,'n1_actual_collar_green_compact_v1');fixture=module(FIXTURE,'n1_actual_collar_green_fixture')
    proof=fixture.run(leaf,c);require(proof['passed'],'independent fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    geom={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True}
    result=leaf.enclose(c,sc=sc,phase=phase_iv,core_W0_Z=core0,core_W1_Z2=core1,forcing_Z=g,matrix_Z=matrix,inner_cells=1,outer_cells=1,source_family=family,source_atoms=atoms,geometry_certificate=geom)
    allowed=set(core_index)|set(source_index);seen=set();classes=set()
    for group in ('core_transport','core_W0_transport','local_signed_sectors','finite_W1_Z2'):
        collection=result[group]
        if group=='core_W0_transport':collections=[terms for fields in collection.values() for terms in fields.values()]
        else:collections=list(collection.values())
        for terms in collections:
            for term in terms:
                origin=term['origin'];require(all(s in allowed for s in origin['source_ids']),'unmapped output source')
                seen.update(origin['source_ids']);classes.update(origin['error_classes'])
    require(seen and classes,'missing source or errors')
    return {'source_family':family,'source_logs':atoms,'domain':a['domain'],'sc':a['sc'],'phase':enc(phase_iv),
      'fixture':proof,'source_index_receipt':INLET,'core_source_index_in_parent_receipt':True,'local_source_index':source_index,
      'prior_forcing_source_index_in_parent_receipt':True,'forcing_groups_without_radial_parent_ids':sum(not x['parent_source_ids'] for x in source_index.values()),'distinct_referenced_input_ids':len(seen),'output_error_classes':sorted(classes),
      'result':enc(result),'alternative_cauchy_receipt':CAUCHY,'allword_tail_receipt':TAIL,'allword_tail_added':False,
      'whole_collar_uniform_input_enclosures_not_point_values':True,'interval_grouping_not_identity':True,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
