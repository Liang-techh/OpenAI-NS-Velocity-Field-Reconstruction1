"""One-entry saved-W0/row5 baseline consumer; no ancestor code is called."""
from __future__ import annotations
import base64,gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
NAME='row5_baseline_four_cell_v4';DIR=ROOT/'tmp'/NAME
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LEAF='src/openai_ns_reconstruction/w0_joined_jets_v4.py';FIXTURE='src/openai_ns_reconstruction/w0_joined_jets_fixture_v4.py'
PRODUCT='src/openai_ns_reconstruction/row5_baseline_product_v4.py';PRODUCT_FIXTURE='src/openai_ns_reconstruction/row5_baseline_product_fixture_v4.py'
RUNNER='scripts/run_row5_baseline_four_cell_v4_once.py'
OLD='experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz';Z2='experiments/root_st073/w0_z2_four_cell_v2_receipt.json.gz'
AFFINE='experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz';FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz'
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz'
FAILURE='experiments/root_st073/row5_baseline_four_cell_v1_failure.json.gz'
FAILURE2='experiments/root_st073/row5_baseline_four_cell_v2_failure.json.gz'
FAILURE3='experiments/root_st073/row5_baseline_four_cell_v3_failure.json.gz'
FILES=(LEAF,FIXTURE,PRODUCT,PRODUCT_FIXTURE,RUNNER,OLD,Z2,AFFINE,FIRST,FORCING,FAILURE,FAILURE2,FAILURE3)
PIN={FAILURE3:'34dccd5503648d2ab525bb5ea00594fca6bed56b20b3082fd602adcc2a407edd',FAILURE2:'0615f5a14f8d3690c8cf08d1b98fed91656c6fa1e4892b16ebabbae1c3840069',FAILURE:'01359e0b7b2b9bd298ccff1bb34a5b568bc0fd88be9f2f4dc8da4ccf449721bf',OLD:'4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06',Z2:'751569810106ef47969f331da12f4d7d369e6fee8e30908e58c5bbc9b7f229dc',AFFINE:'b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2'}
ATOMS=('hb','F0_anchor','Pstar_squared','Lambda')
def sha(x):return hashlib.sha256(x).hexdigest()
def require(x,m):
 if not x:raise ValueError(m)
def save(path,value):
 raw=(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
 with path.open('xb') as f:f.write(raw)
 return sha(raw)
def encode(x):
 if hasattr(x,'_mpi_'):
  return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
 if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [encode(v) for v in x]
 if isinstance(x,(str,int,float,bool)) or x is None:return x
 raise TypeError('unserializable '+type(x).__name__)
def decode(c,x):
 if isinstance(x,dict):
  if 'lower_exact_mpf_tuple' in x or 'upper_exact_mpf_tuple' in x:
   lo=tuple(x['lower_exact_mpf_tuple']);hi=tuple(x['upper_exact_mpf_tuple'])
   require(len(lo)==len(hi)==4 and all(type(n) is int for n in lo+hi),'invalid interval endpoints')
   v=c.make_mpf((lo,hi));require(v.a<=v.b,'reversed interval');return v
  return {k:decode(c,v) for k,v in x.items()}
 if isinstance(x,list):return [decode(c,v) for v in x]
 return x
def module(path,name):
 require(path in (LEAF,FIXTURE,PRODUCT,PRODUCT_FIXTURE),'unexpected module')
 spec=importlib.util.spec_from_file_location(name,ROOT/path);require(spec is not None and spec.loader is not None,'module spec missing')
 obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def preflight():
 require(not ENTRY.exists() and not RESULT.exists(),'entry already exists')
 brief_raw=BRIEF.read_bytes();brief=json.loads(brief_raw);source_raw=SOURCE.read_bytes();manifest=json.loads(source_raw)
 require(brief['component']==NAME and brief['phase']==NAME+'_pending' and brief['goal_active'] is True and brief['full_goal_complete'] is False,'brief guard')
 require(brief['one_entry_only'] is True and brief['do_not_reenter'] is True,'entry guard')
 ck=brief['source_checkpoint'];require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(source_raw),'checkpoint')
 require(sha(Path(ck['preimage_path']).read_bytes())==ck['preimage_sha256'],'brief preimage')
 require(manifest['component']==NAME and manifest['one_entry_only'] is True,'manifest')
 hashes={}
 for row in manifest['files']:
  path=row['path'];require(path in FILES and path not in hashes,'unexpected source path')
  h=sha((ROOT/path).read_bytes());require(h==row['sha256'] and (path not in PIN or h==PIN[path]),'hash mismatch '+path);hashes[path]=h
 require(set(hashes)==set(FILES),'missing source path')
 return hashes,sha(source_raw),sha(brief_raw)
def receipt(raw,path):
 x=json.loads(gzip.decompress(raw[path]));require(x['result']['status']=='returned','not returned '+path);return x
def consume(hashes,outcome):
 raw={}
 for path in FILES:
  b=(ROOT/path).read_bytes();require(sha(b)==hashes[path],'post-entry hash '+path);raw[path]=b
 c=MPIntervalContext();c.prec=1800
 join=module(LEAF,NAME+'_join');join_fixture=module(FIXTURE,NAME+'_join_fixture')
 product=module(PRODUCT,NAME+'_product');product_fixture=module(PRODUCT_FIXTURE,NAME+'_product_fixture')
 jp=join_fixture.run_fixture(join.build,c);outcome['join_fixture']=encode(jp);require(jp.get('passed') is True,'join fixture failed')
 pp=product_fixture.run_fixture(product.build,c);outcome['product_fixture']=encode(pp);require(pp.get('passed') is True,'product fixture failed')
 old=receipt(raw,OLD);affine=receipt(raw,AFFINE)
 new_archive=json.loads(gzip.decompress(raw[Z2]));new_result=json.loads(base64.b64decode(new_archive['records_base64']['result.json'],validate=True))
 new_manifest=json.loads(base64.b64decode(new_archive['records_base64']['source.json'],validate=True))
 require(new_result['status']=='returned' and new_result['component']=='w0_z2_four_cell_v2','Z2 result')
 require(sha(base64.b64decode(new_archive['records_base64']['source.json'],validate=True))==new_result['source_sha256'],'Z2 manifest digest')
 new_pins={row['path']:row['sha256'] for row in new_manifest['files']}
 old_pins={row['path']:row['sha256'] for row in old['source']['files']}
 affine_pins={row['path']:row['sha256'] for row in affine['source']['files']}
 require(new_pins[OLD]==hashes[OLD] and new_pins[FORCING]==hashes[FORCING],'Z2 input provenance')
 require(old_pins[FORCING]==hashes[FORCING] and affine_pins[FORCING]==hashes[FORCING],'shared forcing provenance')
 old_o=old['result']['output'];new_o=new_result['output'];aff_o=affine['result']['output']
 family=old_o['source_family'];logs=old_o['source_logs']
 require(new_o['source_family']==family and aff_o['source_family']==family,'family mismatch')
 require(new_o['source_logs']==logs and aff_o['source_logs']==logs,'source logs mismatch')
 require(old_o['domain']==new_o['domain'],'domain mismatch')
 require(old_o['allword_tail_added'] is False and old_o['appended_to_cumulative_m15'] is False,'old scope')
 require(new_o['allword_tail_added'] is False and new_o['ledger_appended'] is False,'new scope')
 old_payload=decode(c,old_o['result']);new_payload=decode(c,new_o['result'])
 # Atom logs are opaque source records, not interval values to materialize.
 old_payload['source_atoms']=old_o['result']['source_atoms']
 packets=[]
 for role,payload in (('z0_z1',old_payload),('z2',new_payload)):
  packets.append({'source_verified':True,'packet_role':role,'source_family':family,'source_atoms':logs,'source_power_basis':list(ATOMS),'payload':payload})
 outcome['actual_join_calls']=1
 joined=join.build(*packets,source_family=family,source_atoms=logs)
 require(joined['core_retained_once'] is True and joined['g5_added'] is False and joined['m15_ledger_changed'] is False,'join scope')
 state=joined['state']
 # Add receipt locations without changing original origin, errors or formal powers.
 for cell,entries in state.items():
  for (component,order),sectors in entries.items():
   receipt_path=OLD if order in (0,1) else Z2
   field={1:'F1',2:'Uz1',3:'K1',5:'d_x_F1'}[component]
   for sector in ('core','local'):
    for index,term in enumerate(sectors[sector]):
     term.setdefault('input_id',f'W0:{order}:{component}:{sector}:{index}' + (f':cell{cell}' if sector=='local' else ''))
     if order in (0,1):
      parts=['result','output','result', 'core_transport',str(order),field,str(index)] if sector=='core' else ['result','output','result','local',str(order),str(cell),field,str(index)]
     else:
      parts=['records_base64','result.json','output','result','core_transport',field,str(index)] if sector=='core' else ['records_base64','result.json','output','result','local',str(cell),field,str(index)]
     pointer=(receipt_path+'#base64-json(records_base64/result.json)/'+'/'.join(parts[2:])) if order==2 else (receipt_path+'#/'+'/'.join(parts))
     term['binding_receipt_location']={'receipt':receipt_path,'json_path':parts,'embedded_result_decoding':'base64 JSON' if order==2 else None}
     term['binding_receipt_path']=pointer
     term.setdefault('input_path',pointer)
 expression=decode(c,aff_o['expression'])
 outcome['actual_product_calls']=1
 products=product.build(c,expression,state)
 require(products['schema']=='row5_baseline_product_v4' and products['g5_repeated_here'] is False,'product scope')
 require(products['core_local_separate'] is True and products['formal_atoms_materialized'] is False,'formal product scope')
 require(len(products['skipped_g5'])==32,'g5 occurrence count')
 return {'source_family':family,'source_logs':logs,'domain':old_o['domain'],'joined':encode(joined),'products':encode(products),
  'input_hashes':{p:hashes[p] for p in (OLD,Z2,AFFINE,FIRST,FORCING)},'firstword_receipt_raw_only':True,
  'firstword_g5_provenance':'published firstword v4 receipt pin, not re-entered',
  'fixed_normalized_x_identity_proved':False,'live_owner_proved':False,'gauge_proved':False,
  'higher_feedback_remainder_bound':False,'endpoint_cancellation_proved':False,'m15_ledger_changed':False}
def main():
 hashes,source_sha,brief_sha=preflight();entry_sha=save(ENTRY,{'component':NAME,'source_sha256':source_sha,'brief_sha256':brief_sha,'input_hashes':hashes,'entered':True})
 out={'component':NAME,'source_sha256':source_sha,'entry_sha256':entry_sha,'status':'failed'}
 try:out.update(status='returned',output=consume(hashes,out),error=None)
 except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
 result_sha=save(RESULT,out);print(json.dumps({'component':NAME,'status':out['status'],'result_sha256':result_sha}))
 return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
