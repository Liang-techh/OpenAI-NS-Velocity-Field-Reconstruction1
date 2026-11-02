"""One-entry direct binding of saved g5 source occurrences across row-5 branches."""
from __future__ import annotations
import base64,gc,gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
NAME='row5_g5_binding_v1';DIR=ROOT/'tmp'/NAME;SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LEAF='src/openai_ns_reconstruction/row5_g5_binding_v1.py';FIXTURE='src/openai_ns_reconstruction/row5_g5_binding_fixture_v1.py';RUNNER='scripts/run_row5_g5_binding_v1_once.py'
AFFINE='experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz';FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';BASE='experiments/root_st073/row5_baseline_four_cell_v4_receipt.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,AFFINE,FIRST,BASE)
PIN={AFFINE:'b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',BASE:'866fae8b0cc0471eceef78c727d7be6d52d413ffefb2679881b3100291aea8ca'}
def sha(x):return hashlib.sha256(x).hexdigest()
def require(x,m):
 if not x:raise ValueError(m)
def save(path,value):
 raw=(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
 with path.open('xb') as f:f.write(raw)
 return sha(raw)
def module(path,name):
 require(path in (LEAF,FIXTURE),'unexpected module');spec=importlib.util.spec_from_file_location(name,ROOT/path);require(spec is not None and spec.loader is not None,'missing module spec');obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def preflight():
 require(not ENTRY.exists() and not RESULT.exists(),'already entered');braw=BRIEF.read_bytes();b=json.loads(braw);sraw=SOURCE.read_bytes();m=json.loads(sraw)
 require(b['component']==NAME and b['phase']==NAME+'_pending' and b['goal_active'] is True and b['full_goal_complete'] is False,'brief guard');require(b['one_entry_only'] is True and b['do_not_reenter'] is True,'entry guard')
 ck=b['source_checkpoint'];require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(sraw),'checkpoint');require(sha(Path(ck['preimage_path']).read_bytes())==ck['preimage_sha256'],'preimage')
 require(m['component']==NAME and m['one_entry_only'] is True,'manifest');h={}
 for row in m['files']:
  p=row['path'];require(p in FILES and p not in h,'unexpected source');digest=sha((ROOT/p).read_bytes());require(digest==row['sha256'] and (p not in PIN or digest==PIN[p]),'hash '+p);h[p]=digest
 require(set(h)==set(FILES),'missing source');return h,sha(sraw),sha(braw)
def consume(h,out):
 raw={}
 for p in FILES:
  data=(ROOT/p).read_bytes();require(sha(data)==h[p],'changed after entry '+p);raw[p]=data
 leaf=module(LEAF,NAME+'_leaf');fixture=module(FIXTURE,NAME+'_fixture')
 proof=fixture.run_fixture(leaf.build);out['fixture']=proof;require(isinstance(proof,dict) and proof.get('passed') is True,'synthetic fixture failed')
 affine=json.loads(gzip.decompress(raw[AFFINE]));first=json.loads(gzip.decompress(raw[FIRST]));archive=json.loads(gzip.decompress(raw[BASE]))
 base_source_raw=base64.b64decode(archive['records_base64']['source.json'],validate=True)
 base_result_raw=base64.b64decode(archive['records_base64']['result.json'],validate=True)
 base_manifest=json.loads(base_source_raw);base_result=json.loads(base_result_raw);del archive,base_result_raw;gc.collect()
 require(affine['result']['status']=='returned' and first['result']['status']=='returned' and base_result['status']=='returned','saved receipt status')
 require(sha(base_source_raw)==base_result['source_sha256'],'baseline embedded source SHA')
 first_pins={row['path']:row['sha256'] for row in first['source']['files']};base_pins={row['path']:row['sha256'] for row in base_manifest['files']}
 require(first_pins[AFFINE]==h[AFFINE] and base_pins[AFFINE]==h[AFFINE] and base_pins[FIRST]==h[FIRST],'common parent receipt hashes')
 a=affine['result']['output'];f=first['result']['output'];b=base_result['output']
 require(a['source_family']==f['source_family']==b['source_family'],'source family mismatch')
 require(a['source_logs']==f['source_logs']==b['source_logs'],'source logs mismatch')
 require(f['input_receipt_sha256'][AFFINE]==h[AFFINE] and b['input_hashes'][AFFINE]==h[AFFINE] and b['input_hashes'][FIRST]==h[FIRST],'saved input hashes mismatch')
 first_product=f['result'];base_product=b['products'];expression=a['expression']
 aw={'receipt':{'path':AFFINE,'sha256':h[AFFINE],'json_path_prefix':['result','output','expression'],'encoding':'gzip-json'},'expression':expression}
 fw={'receipt':{'path':FIRST,'sha256':h[FIRST],'json_path_prefix':['result','output','result'],'encoding':'gzip-json'},'result':first_product}
 bw={'receipt':{'path':BASE,'sha256':h[BASE],'json_path_prefix':['records_base64','result.json','output','products'],'encoding':'gzip-json then base64-json'},'products':base_product}
 out['actual_binding_calls']=1
 result=leaf.build(aw,fw,bw)
 require(result['schema']==NAME and result['mapping_count']==32,'binding count or schema')
 require(result['exact_once']['affine_occurrence_to_firstword_source_only_output'] is True and result['exact_once']['baseline_skipped_has_no_product'] is True,'g5 exact-once check')
 require(len(result['cell_combinations'])==4,'four cell references')
 return {'source_family':a['source_family'],'source_logs':a['source_logs'],'input_hashes':{p:h[p] for p in (AFFINE,FIRST,BASE)},'result':result,'actual_binding_calls':1,'ancestor_called':False,'same_live_owner_proved':False,'physical_h5_proved':False,'m15_ledger_changed':False}
def main():
 h,s,b=preflight();e=save(ENTRY,{'component':NAME,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':NAME,'source_sha256':s,'entry_sha256':e,'status':'failed'}
 try:out.update(status='returned',output=consume(h,out),error=None)
 except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
 d=save(RESULT,out);print(json.dumps({'component':NAME,'status':out['status'],'result_sha256':d}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
