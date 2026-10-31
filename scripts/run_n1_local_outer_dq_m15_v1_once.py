"""One-entry outer potential applied to the saved local raw-Q word."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_local_outer_dq_m15_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_actual_collar_low_word_v1_receipt.json.gz';INLET='experiments/root_st073/n1_actual_core_inlet_v1_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';PREVIOUS='experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz';LOCAL='experiments/root_st073/n1_local_g5_potential_physical_m15_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_local_outer_dq_m15_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_local_outer_dq_m15_fixture_v1.py';RUNNER='scripts/run_n1_local_outer_dq_m15_v1_once.py'
REQUIRED=(LOW,INLET,KCHI,BIND,PRESS,SWIRL,CORRECTION,SELECT,PREVIOUS,LOCAL,LEAF,FIXTURE,RUNNER)
FIXED={LOW:'f4ed142fd3c8ec79b919af92ad12f3e160b06a6c8b13a932bb4fc17df27df234',INLET:'722c60ce802d5b1adcf8e4735406a791c5c58896ee6d27f0da040dc0e4c42d34',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',PREVIOUS:'48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2',LOCAL:'55f25fa95c455f03f3e2af99645da3f29590040b1732f0b266275efb35dc10ba'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_local_outer_dq_m15_v1');fixture=module(FIXTURE,'n1_local_outer_dq_m15_fixture_v1');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={p:json.loads(gzip.decompress(raw[p])) for p in (LOW,INLET,KCHI,BIND,PRESS,SWIRL,PREVIOUS)};require(all(r['result']['status']=='returned' for r in receipts.values()),'receipt status');low,inlet,kchi,bind,p,f,prev=(receipts[k]['result']['output'] for k in (LOW,INLET,KCHI,BIND,PRESS,SWIRL,PREVIOUS));selected=json.loads(raw[SELECT]);correction=json.loads(raw[CORRECTION]);family=low['source_family'];atoms=low['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (inlet,kchi,bind,p,f,prev)) and selected['source_family']==family,'common source family and logs');lp={r['path']:r['sha256'] for r in receipts[LOW]['source']['files']};bp={r['path']:r['sha256'] for r in receipts[BIND]['source']['files']};pp={r['path']:r['sha256'] for r in receipts[PREVIOUS]['source']['files']};require(lp[INLET]==hashes[INLET] and bp[KCHI]==hashes[KCHI] and all(pp[k]==hashes[k] for k in (SWIRL,CORRECTION,SELECT)),'authenticated source chains');require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and correction['input_id']=='F0_connection','leading amplitude correction')
    require(low['whole_collar_uniform_input_enclosures_not_point_values'] and low['core_W1_Z0_Z1_inlet_missing'] and low['allword_tail_added'] is False and low['result']['core_W1_inlet_included'] is False and low['source_index_receipt']==INLET,'local-only whole-phase word');require(bind['source_checks']['one_query_jet_shared_in_integration'] is True and bind['phase_shear_cutoff_identified_with_axis_potential'] is False,'native collar axis chi');contract=selected['actual_cutoff_complete_support_measure'];sc=iv(c,p['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0 and low['sc']==p['sc'] and f['sc']==p['sc'] and contract['passed'] and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'physical geometry');source_index={};core_index=inlet['inlet']['source_index'];local_index=low['local_source_index'];W={};chi={};lead={};counts={}
    for order in (0,1):
        z='Z'+str(order);W[z]=[]
        for i,t in enumerate(low['result']['local_signed_sectors'][str(order)]['F1']):
            sid=f'local_W1_F:{order}:{i}';origin=t['origin'];parents=origin['source_ids'];require(parents and all(x in core_index or x in local_index for x in parents),'unindexed local raw-Q parent');require(not origin['uncertainty'] or origin['error_classes'],'lost inherited class');source_index[sid]={'receipt':LOW,'path':['result','output','result','local_signed_sectors',str(order),'F1',i],'parent_source_ids':parents,'parent_origin':origin};W[z].append({'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':[sid],'source_paths':['source_index:'+sid],'error_classes':origin['error_classes'],'uncertainty':origin['uncertainty']})
        counts[z]=len(W[z]);require(counts[z]==137,'actual local word F sector count')
        sid=f'chi:Z{order}';jet=bind['saved_chi_Z0_Z1'][str(order)]['coefficient'];require(jet==kchi['chi']['chi'][order],'same captured axis chi jet');source_index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',order],'binding_receipt':BIND};chi[z]=[{'coefficient':iv(c,jet),'powers':[0,0,0,0],'source_ids':[sid],'source_paths':['source_index:'+sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}]
        terms=[]
        for i,t in enumerate(f['F0_connection']['leading'][str(order)]):
            sid=f'lead_F:{order}:{i}';source_index[sid]={'receipt':SWIRL,'path':['result','output','F0_connection','leading',str(order),i],'parent_source_ids':t['origin']['source_ids'],'connection_path_correction':CORRECTION};require(prev['source_index'][sid]==source_index[sid],'same leading source path');terms.append({'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':[sid],'source_paths':['source_index:'+sid],'error_classes':t['origin']['error_classes'],'uncertainty':t['origin']['uncertainty']})
        require(terms,'missing leading F');lead[z]=terms
    cert={'hb':c.mpf([0,'0.5']),'strict_hb_admitted':True,'geometry_Z_independent':True,'q_all_phase':True,'axis_chi_is_not_phase_cutoff':True,'original_cutoff_mass_5_over_8':True};result=leaf.build(c,sc=sc,low_word=W,axis_chi=chi,leading_F=lead,source_family=family,source_atoms=atoms,certificate=cert);require(result['schema']=='n1_local_outer_dq_m15_v1','schema');seen=set();out={};out_counts={}
    for family_name,key in (('field','field_intermediate'),('m15','m15')):
        rows=result[key]['F1'] if family_name=='field' else result[key];out[family_name]={}
        for z in ('Z0','Z1'):
            terms=rows[z];out_counts[family_name+':'+z]=len(terms)
            for t in terms:
                require(all(sid in source_index for sid in t['source_ids']),'unindexed output');seen.update(t['source_ids'])
            out[family_name][z]=terms
    require(result['ledger_written'] is False and result['same_source_cancellation_claimed'] is False,'ledger/cancellation scope');return {'source_family':family,'source_logs':atoms,'source_index':source_index,'core_parent_index_reference':{'receipt':INLET,'path':['result','output','inlet','source_index']},'local_parent_index':local_index,'input_counts':counts,'fixture':proof,'counts':out_counts,'used_source_id_count':len(seen),'result':enc({'kernel':result['kernel'],'geometry':result['geometry'],'field':out['field'],'m15':out['m15']}),'local_Dg_reference_only':{'receipt':LOCAL,'sha256':hashes[LOCAL]},'P_compensation_included':False,'core_W1_inlet_included':False,'higher_outer_tail_included':False,'old_prefix_disjointness_proved':False,'old_errors_removed':False,'appended_to_ledger':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
