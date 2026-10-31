"""One-entry physical m15 projection of saved local-g5 potential words."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_local_g5_potential_physical_m15_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
FIRST='experiments/root_st073/n1_collar_kchi_source_first_v3_receipt.json.gz';TAIL='experiments/root_st073/n1_collar_local_g5_j2_tail_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';PREVIOUS='experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz';W0='experiments/root_st073/n1_w0_physical_m15_v1_receipt.json.gz';E00='experiments/root_st073/n1_e00_physical_m15_cells_v2_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_local_g5_potential_physical_m15_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_local_g5_potential_physical_m15_fixture_v1.py';RUNNER='scripts/run_n1_local_g5_potential_physical_m15_v1_once.py'
REQUIRED=(FIRST,TAIL,PRESS,SWIRL,CORRECTION,SELECT,PREVIOUS,W0,E00,LEAF,FIXTURE,RUNNER)
FIXED={FIRST:'bdfd705e43db41badd69b737a014f9ff288eaa58a18060fe736aebb196f18262',TAIL:'67894bac58f4fbbcd5437bf970dadafeb918fe2d9638425f8e5165d21af4530b',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',PREVIOUS:'48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2',W0:'f827579eb137c936eb0a2c00fbba873cccf6e8204d00e03d4c478531d0b770ab',E00:'209218dbe33f3f17c4e78d8e31c54d6398e0e13841a5993c48b9c099d64bd2c1'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint'];require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard');require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
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
    for path in REQUIRED:
        v=(ROOT/path).read_bytes();require(sha(v)==hashes[path],'changed after entry '+path);raw[path]=v
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_local_g5_potential_physical_m15_v1');fixture=module(FIXTURE,'n1_local_g5_potential_physical_m15_fixture_v1');proof=fixture.run_fixture(leaf.build,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={p:json.loads(gzip.decompress(raw[p])) for p in (FIRST,TAIL,PRESS,SWIRL,PREVIOUS)};require(all(r['result']['status']=='returned' for r in receipts.values()),'receipt status');a,b,p,f,prev=(receipts[k]['result']['output'] for k in (FIRST,TAIL,PRESS,SWIRL,PREVIOUS));selected=json.loads(raw[SELECT]);correction=json.loads(raw[CORRECTION]);family=a['source_family'];atoms=a['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (b,p,f,prev)) and selected['source_family']==family,'common source family and logs');pa={r['path']:r['sha256'] for r in receipts[FIRST]['source']['files']};pb={r['path']:r['sha256'] for r in receipts[TAIL]['source']['files']};pp={r['path']:r['sha256'] for r in receipts[PREVIOUS]['source']['files']};require(pb[FIRST]==hashes[FIRST] and all(pa[k]==pb[k] for k in ('experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz','experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz',PRESS,SELECT)) and all(pp[k]==hashes[k] for k in (SWIRL,CORRECTION,SELECT)),'authenticated local and leading chains');require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and correction['input_id']=='F0_connection','amplitude correction')
    aj=a['actual_first_correction'];bj=b['actual_local_j2_tail'];require(aj['whole_domain_source_box'] is True and aj['kernel_mass_definition']=='safe positive whole-cell interval geometry at q*=1' and iv(c,aj['target_q']).a==iv(c,aj['target_q']).b==1,'whole phase source and first word mass');require(b['first_potential_reference']['sha256']==hashes[FIRST] and b['local_j2_tail_only'] and b['m15_computed'] is False,'j2 source scope');require(a['source_index']==b['source_index'],'same raw local source index');contract=selected['actual_cutoff_complete_support_measure'];sc=iv(c,p['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0 and f['sc']==p['sc'] and contract['passed'] and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'physical geometry');require(aj['geometry']['sc_interval']==p['sc'] and bj['geometry']['sc_interval']==p['sc'],'same sc in local words')
    index={};packets={'J1':{},'J2':{}};input_counts={}
    for family_name,receipt,obj in (('J1',FIRST,aj),('J2',TAIL,bj)):
        for order in (0,1):
            z='Z'+str(order);items=obj['response']['F1' if family_name=='J1' else 'F'][z];out=[]
            for i,t in enumerate(items):
                sid=f'local_{family_name}:{z}:{i}';parents=t.get('source_ids',t['origin'].get('source_ids',[]));require(parents and any(str(x).startswith('collar_g5:') for x in parents),'local forcing provenance');index[sid]={'receipt':receipt,'path':['result','output','actual_first_correction' if family_name=='J1' else 'actual_local_j2_tail','response','F1' if family_name=='J1' else 'F',z,i],'parent_source_ids':parents,'chi_order':t.get('chi_order'),'source_order':t.get('source_order'),'error_classes':t.get('error_classes',[])};out.append({'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':[sid],'source_paths':['source_index:'+sid],'error_classes':t.get('error_classes',[]),'uncertainty':t.get('uncertainty',False),'error_only':t.get('error_only',False),'q_endpoint':c.mpf(1) if family_name=='J1' else None})
            require(len(out)==(16 if order==0 else 32),'local source term count');packets[family_name][z]=out;input_counts[family_name+':'+z]=len(out)
    lead={}
    for order in (0,1):
        z='Z'+str(order);terms=[]
        for i,t in enumerate(f['F0_connection']['leading'][str(order)]):
            sid=f'lead_F:{order}:{i}';index[sid]={'receipt':SWIRL,'path':['result','output','F0_connection','leading',str(order),i],'parent_source_ids':t['origin']['source_ids'],'connection_path_correction':CORRECTION};require(prev['source_index'][sid]==index[sid],'same leading field source path');terms.append({'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':[sid],'source_paths':['source_index:'+sid],'error_classes':t['origin']['error_classes'],'uncertainty':t['origin']['uncertainty'],'error_only':False})
        require(terms,'whole phase leading source');lead[z]=terms
    cert={'hb':c.mpf([0,'0.5']),'strict_hb_admitted':True,'source_chi_g_full_q_enclosure':True,'unified_exponent_bound':True,'geometry_Z_independent':True,'original_cutoff_mass_5_over_8':True};result=leaf.build(c,sc,packets['J1'],packets['J2'],lead,cert);require(result['schema']=='n1_local_g5_potential_physical_m15_v1' and result['ledger_212_written'] is False and result['chi_reapplied'] is False,'scope');counts={};seen=set()
    for group in ('J1','J2'):
        for z in ('Z0','Z1'):
            rows=result[group][z];counts[group+':'+z]=len(rows)
            for t in rows:
                require(all(sid in index for sid in t['source_ids']),'unindexed projection');require(group!='J2' or (t['error_only'] is True and t['symmetric_error'] is True),'tail must stay symmetric error');seen.update(t['source_ids'])
    return {'source_family':family,'source_logs':atoms,'source_index':index,'parent_local_source_index':a['source_index'],'input_counts':input_counts,'fixture':proof,'counts':counts,'used_source_id_count':len(seen),'projection':enc({'geometry':result['geometry'],'J1':result['J1'],'J2':result['J2'],'J1_contract':result['J1_contract'],'J2_contract':result['J2_contract']}),'reference_only':{W0:hashes[W0],E00:hashes[E00],PRESS:hashes[PRESS]},'G0g_recomputed':False,'pressure_recomputed':False,'old_prefix_disjointness_proved':False,'old_errors_removed':False,'appended_to_ledger':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
