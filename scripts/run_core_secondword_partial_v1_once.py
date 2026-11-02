"""One-entry core second-feedback finite-prefix consumer."""
from __future__ import annotations
import gc,gzip,hashlib,importlib.util,io,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';NAME='core_secondword_partial_v1';DIR=ROOT/'tmp'/NAME;SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LEAF='src/openai_ns_reconstruction/core_secondword_partial_v1.py';FIXTURE='src/openai_ns_reconstruction/core_secondword_partial_fixture_v1.py';RUNNER='scripts/run_core_secondword_partial_v1_once.py'
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';PROFILE='experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz';JOIN='experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz';FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';LOCAL='experiments/root_st073/row5_secondword_local_v4_compact_receipt.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,LOW,PROFILE,JOIN,FIRST,LOCAL)
PIN={LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PROFILE:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',JOIN:'a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',LOCAL:'70ce750f17202fc59a8f39d61eab7da864fc447be6a8ea9c79a3e625ea34ea55'}
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
def decode(c,x):
 if isinstance(x,dict):
  if 'lower_exact_mpf_tuple' in x or 'upper_exact_mpf_tuple' in x:
   lo=tuple(x['lower_exact_mpf_tuple']);hi=tuple(x['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4 and all(type(n) is int for n in lo+hi),'invalid endpoints');v=c.make_mpf((lo,hi));require(v.a<=v.b,'reversed endpoints');return v
  return {k:decode(c,v) for k,v in x.items()}
 if isinstance(x,list):return [decode(c,v) for v in x]
 return x
def module(path,name):
 require(path in (LEAF,FIXTURE),'unexpected import');spec=importlib.util.spec_from_file_location(name,ROOT/path);require(spec is not None and spec.loader is not None,'module spec');obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def preflight():
 require(not ENTRY.exists() and not RESULT.exists(),'already entered');br=BRIEF.read_bytes();b=json.loads(br);sr=SOURCE.read_bytes();m=json.loads(sr)
 require(b['component']==NAME and b['phase']==NAME+'_pending' and b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True,'brief guard');ck=b['source_checkpoint'];require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(sr),'checkpoint');require(sha(Path(ck['preimage_path']).read_bytes())==ck['preimage_sha256'],'preimage');require(m['component']==NAME and m['one_entry_only'] is True,'manifest');h={}
 for row in m['files']:
  p=row['path'];require(p in FILES and p not in h,'unexpected source');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in PIN or v==PIN[p]),'hash mismatch '+p);h[p]=v
 require(set(h)==set(FILES),'missing source');return h,sha(sr),sha(br)
def receipt(raw,path):
 obj=json.loads(gzip.decompress(raw[path]));require(obj['result']['status']=='returned','input receipt '+path);return obj
def local_manifest(data):
 with gzip.GzipFile(fileobj=io.BytesIO(data),mode='rb') as z:
  line=z.readline()
 prefix=b'{"component":"row5_secondword_local_v4","source":';require(line.startswith(prefix),'local compact source envelope');return json.loads(line[len(prefix):])
def consume(h,out):
 raw={}
 for p in FILES:
  data=(ROOT/p).read_bytes();require(sha(data)==h[p],'changed after entry '+p);raw[p]=data
 c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,NAME+'_leaf');fixture=module(FIXTURE,NAME+'_fixture')
 proof=fixture.run_fixture(leaf.build,c);out['fixture']=encode(proof);require(proof.get('passed') is True,'synthetic fixture failed')
 local_m=local_manifest(raw[LOCAL]);lm={x['path']:x['sha256'] for x in local_m['files']};require(lm[FIRST]==h[FIRST],'local firstword parent')
 first=receipt(raw,FIRST);fm={x['path']:x['sha256'] for x in first['source']['files']};require(fm[JOIN]==h[JOIN],'firstword join parent');first_family=first['result']['output']['source_family'];first_logs=first['result']['output']['source_logs'];del first;gc.collect()
 join=receipt(raw,JOIN);jm={x['path']:x['sha256'] for x in join['source']['files']};require(jm[LOW]==h[LOW],'join core parent');jo=join['result']['output'];require(jo['source_family']==first_family and jo['source_logs']==first_logs,'join family/logs');del join,jo;gc.collect()
 low=receipt(raw,LOW);profile=receipt(raw,PROFILE);pm={x['path']:x['sha256'] for x in profile['source']['files']};require(pm[LOW]==h[LOW],'profile core parent');lo=low['result']['output'];po=profile['result']['output']
 require(lo['source_family']==po['source_family']==first_family and lo['source_logs']==po['source_logs']==first_logs,'same saved family and logs')
 require(lo['saved_W0_outer_cells_reused'] is True and lo['allword_tail_added'] is False and po['allword_tail_added'] is False and po['appended_to_cumulative_m15'] is False,'saved scope')
 require(po['source_kernel_scope']=='K0-only nested finite contribution; not a full Kchi W1 enclosure' and po['Kchi_minus_K0_included'] is False and po['complete_m15_computed'] is False,'profile kernel scope')
 pcells=po['profile']['target_cells'];mcells=lo['matrix_cells'];require(len(pcells)==len(mcells)==8,'eight cells')
 fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');targets=[];matrices=[]
 for i in range(8):
  cell=pcells[i];require(cell['index']==i,'profile target cell index');by={}
  for order in ('Z0','Z1'):
   by[order]={}
   for field in fields:
    terms=[]
    for j,raw_term in enumerate(cell['W1_by_order'][order][field]):
     term=decode(c,raw_term);require(term.get('target_cell',i)==i,'profile term target cell');term['input_id']=f'profile:cell{i}:{order}:{field}:{j}';term['receipt_input_path']=[PROFILE,'result','output','profile','target_cells',i,'W1_by_order',order,field,j];term.setdefault('source_path',term['receipt_input_path']);terms.append(term)
    by[order][field]=terms
  targets.append({'index':i,'t':decode(c,cell['t']),'rho':decode(c,cell['rho']),'W1_by_order':by})
  mcell=mcells[i];assembled={'t':decode(c,mcell['t']),'rho':decode(c,mcell['rho'])}
  for family,keys in (('C0',('5,1','5,2','5,3','5,5')),('B1',('5,1','5,2','5,3'))):
   assembled[family]={}
   for order in ('0','1'):
    assembled[family][order]={}
    for key in keys:
     terms=[]
     for j,raw_term in enumerate(mcell[family][order][key]):
      term=decode(c,raw_term);term['input_id']=f'matrix:cell{i}:{family}:Z{order}:{key}:{j}';term['receipt_input_path']=[LOW,'result','output','matrix_cells',i,family,order,key,j];term.setdefault('source_path',term['receipt_input_path']);terms.append(term)
     assembled[family][order][key]=terms
  matrices.append(assembled)
 # A new entry on saved inputs, not a call to any previous Green builder.
 out['actual_build_calls']=1
 result=leaf.build(c,targets,matrices)
 require(isinstance(result,dict) and result['schema']==NAME,'result schema')
 require(result['ordinary_Z_orders']==[0,1] and result['complete_secondword'] is False and result['no_cancellation'] is True and result['cancellation_performed'] is False,'result scope')
 require(result['raw_matrix_columns']=={'C0':[1,2,3,5],'B1':[1,2,3]} and result['missing_product_rule']=='B1[0]*W[2]' and bool(result['unknown_ids_by_kind']['missing_B1_W1_Z2']),'result branches')
 require(all(result['terms_by_order'][order][field] for order in ('Z0','Z1') for field in ('F1','d_x_F1')),'empty product terms')
 return {'source_family':lo['source_family'],'source_logs':lo['source_logs'],'input_hashes':{p:h[p] for p in (LOW,PROFILE,JOIN,FIRST,LOCAL)},'result':encode(result),'actual_build_calls':1,'physical_owner_proved':False,'Kchi_compensation_included':False,'allword_tail_added':False,'m15_ledger_changed':False}
def main():
 h,s,b=preflight();e=save(ENTRY,{'component':NAME,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':NAME,'source_sha256':s,'entry_sha256':e,'status':'failed'}
 try:out.update(status='returned',output=consume(h,out),error=None)
 except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
 d=save(RESULT,out);print(json.dumps({'component':NAME,'status':out['status'],'result_sha256':d}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
