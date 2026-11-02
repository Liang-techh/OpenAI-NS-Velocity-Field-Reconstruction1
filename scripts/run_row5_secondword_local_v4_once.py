"""One-entry source-bound raw-K0 row-5 local second-feedback consumer."""
from __future__ import annotations
import base64,gc,gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';NAME='row5_secondword_local_v4';DIR=ROOT/'tmp'/NAME;SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LEAF='src/openai_ns_reconstruction/row5_secondword_local_v4.py';FIXTURE='src/openai_ns_reconstruction/row5_secondword_local_fixture_v4.py';RUNNER='scripts/run_row5_secondword_local_v4_once.py'
FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';BIND='experiments/root_st073/row5_g5_binding_v1_receipt.json.gz';BASE='experiments/root_st073/row5_baseline_four_cell_v4_receipt.json.gz';Z2='experiments/root_st073/w0_z2_four_cell_v2_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
FAILURE='experiments/root_st073/row5_secondword_local_v1_failure.json.gz'
FAILURE2='experiments/root_st073/row5_secondword_local_v2_failure.json.gz'
FAILURE3='experiments/root_st073/row5_secondword_local_v3_failure.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,FIRST,BIND,BASE,Z2,SELECT,FAILURE,FAILURE2,FAILURE3)
PIN={FAILURE3:'70cca25faaea2b9ad6ff7cf3c30db8776d99d76542ef687542a0c83911fe6005',FAILURE2:'083b6cbe55af387de99177087b99400894197c575aa7036e2e5a2598bd3b7db8',FAILURE:'c8446208821b7ab6055e4538e7a344113ee834755b6a625ae5cdac2410818144',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',BIND:'8b8378de335a230bccca1268d634ce36d15d6dcdf35b35b388e8b300d3ca8309',BASE:'866fae8b0cc0471eceef78c727d7be6d52d413ffefb2679881b3100291aea8ca',Z2:'751569810106ef47969f331da12f4d7d369e6fee8e30908e58c5bbc9b7f229dc',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
def sha(x):return hashlib.sha256(x).hexdigest()
def require(x,m):
 if not x:raise ValueError(m)
def save(path,value):
 raw=(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
 with path.open('xb') as f:f.write(raw)
 return sha(raw)
def encode(x):
 if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
 if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [encode(v) for v in x]
 if isinstance(x,(str,int,float,bool)) or x is None:return x
 raise TypeError('unsupported result '+type(x).__name__)
def interval(c,r):
 lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4 and all(type(n) is int for n in lo+hi),'exact interval');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def module(path,name):
 require(path in (LEAF,FIXTURE),'unexpected import');spec=importlib.util.spec_from_file_location(name,ROOT/path);require(spec is not None and spec.loader is not None,'module spec');obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def preflight():
 require(not ENTRY.exists() and not RESULT.exists(),'already entered');br=BRIEF.read_bytes();b=json.loads(br);sr=SOURCE.read_bytes();m=json.loads(sr)
 require(b['component']==NAME and b['phase']==NAME+'_pending' and b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True,'brief guard');ck=b['source_checkpoint'];require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(sr),'checkpoint');require(sha(Path(ck['preimage_path']).read_bytes())==ck['preimage_sha256'],'preimage');require(m['component']==NAME and m['one_entry_only'] is True,'manifest');h={}
 for row in m['files']:
  p=row['path'];require(p in FILES and p not in h,'unexpected source');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in PIN or v==PIN[p]),'hash mismatch '+p);h[p]=v
 require(set(h)==set(FILES),'missing source');return h,sha(sr),sha(br)
def unwrap_archive(raw,path):
 archive=json.loads(gzip.decompress(raw[path]));sr=base64.b64decode(archive['records_base64']['source.json'],validate=True);rr=base64.b64decode(archive['records_base64']['result.json'],validate=True);m=json.loads(sr);result=json.loads(rr);require(sha(sr)==result['source_sha256'] and result['status']=='returned','embedded archive status/source '+path);return m,result

def consume(h,out):
 raw={}
 for p in FILES:
  data=(ROOT/p).read_bytes();require(sha(data)==h[p],'changed after entry '+p);raw[p]=data
 c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,NAME+'_leaf');fixture=module(FIXTURE,NAME+'_fixture')
 proof=fixture.run_fixture(leaf.build,c);out['fixture']=encode(proof);require(proof.get('passed') is True,'synthetic fixture failed')
 first=json.loads(gzip.decompress(raw[FIRST]));bind_manifest,bind_result=unwrap_archive(raw,BIND);z2_manifest,z2_result=unwrap_archive(raw,Z2)
 # Baseline is decoded only to authenticate its embedded source manifest, not to recompute its result.
 archive=json.loads(gzip.decompress(raw[BASE]));base_manifest=json.loads(base64.b64decode(archive['records_base64']['source.json'],validate=True));del archive;gc.collect()
 selected=json.loads(raw[SELECT]);require(first['result']['status']=='returned','firstword receipt status')
 bp={x['path']:x['sha256'] for x in bind_manifest['files']};bap={x['path']:x['sha256'] for x in base_manifest['files']};zp={x['path']:x['sha256'] for x in z2_manifest['files']}
 require(bp[FIRST]==h[FIRST] and bp[BASE]==h[BASE] and bap[FIRST]==h[FIRST] and bap[Z2]==h[Z2] and zp[SELECT]==h[SELECT],'source chain pins')
 fo=first['result']['output'];bo=bind_result['output'];zo=z2_result['output'];family=fo['source_family'];logs=fo['source_logs']
 require(bo['source_family']==family and zo['source_family']==family and selected['source_family']==family,'family mismatch');require(bo['source_logs']==logs and zo['source_logs']==logs,'source logs mismatch')
 require(bo['input_hashes'][FIRST]==h[FIRST] and bo['result']['mapping_count']==32,'g5 binding missing')
 require(fo['domain']==zo['domain'] and fo['same_k0_firstword_construction']['full_collar_domain_equal'] is True and fo['same_k0_firstword_construction']['phase_equal'] is True,'firstword/geometry domain')
 require(fo['scope']['physical_dx_applied'] is False and fo['scope']['physical_h5_materialized'] is False and fo['result']['unknown_residual']['preserved'] is True,'firstword scope')
 require(fo['result']['status']=='conditional_firstword_substituted_expression' and fo['result']['formal_atoms_materialized'] is False and fo['result']['full_h5'] is False and fo['result']['interval_sum_performed'] is False,'firstword product semantics')
 measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] is True and measure['source_family']==family and measure['measure_is_Z_independent'] is True and '0<hb<1/2' in measure['inequality_proof'],'geometry premise')
 sc=interval(c,zo['sc']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'sc mismatch')
 saved_cells=zo['result']['q_cells'];require(len(saved_cells)==4,'saved q cell count');q_cells=[]
 for index,cell in enumerate(saved_cells):
  lo=interval(c,cell['q_lower']);hi=interval(c,cell['q_upper'])
  require(lo._mpi_==(c.mpf(index)/4)._mpi_ and hi._mpi_==(c.mpf(index+1)/4)._mpi_,'saved q partition')
  q_cells.append(c.make_mpf((lo._mpi_[0],hi._mpi_[1])))
 # Actual leaf API is finalized after static review; only this one call is receipt-backed.
 out['actual_build_calls']=1
 result=leaf.build(c, fo['result'], sc=sc, hb_admissibility=True, q_cells=q_cells)
 require(isinstance(result,dict) and result.get('schema')==NAME,'leaf result schema')
 require(len(result['cell_factors'])==4 and result['scope']['raw_K0_C0_B1_feedback_only'] is True,'local geometry/scope')
 require(result['scope']['source_only_g5_skipped'] is True and result['scope']['axis_plus_2chi_skipped'] is True and result['scope']['axis_skip_is_algebraic_cancellation'] is False,'source split')
 require(result['scope']['formal_hb_materialized'] is False and result['scope']['no_error_cancellation'] is True and result['scope']['complete_Rhi_bound'] is False,'formal/remainder scope')
 require(result['core_feedback_inlet']['known'] is False and result['core_feedback_inlet']['reset_to_zero'] is False,'nonzero core not discarded')
 require(all(result['product_group_ids_by_row'][row] for row in ('h5','h5_Z1')),'missing feedback order')
 return {'source_family':family,'source_logs':logs,'input_hashes':{p:h[p] for p in (FIRST,BIND,BASE,Z2,SELECT)},'sc':zo['sc'],'result':encode(result),'same_live_owner_proved':False,'core_secondword_inlet_bound':False,'higher_feedback_bound':False,'m15_ledger_changed':False}
def main():
 h,s,b=preflight();e=save(ENTRY,{'component':NAME,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':NAME,'source_sha256':s,'entry_sha256':e,'status':'failed'}
 try:out.update(status='returned',output=consume(h,out),error=None)
 except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
 d=save(RESULT,out);print(json.dumps({'component':NAME,'status':out['status'],'result_sha256':d}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
