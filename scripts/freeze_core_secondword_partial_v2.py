"""Freeze the raw K0 second-feedback core-exit source consumer."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';NAME='core_secondword_partial_v2';DIR=ROOT/'tmp'/NAME;PRE=Path('C:/Users/z5242/Documents/Codex/core_secondword_v2_pre.json')
LEAF='src/openai_ns_reconstruction/core_secondword_partial_v2.py';FIXTURE='src/openai_ns_reconstruction/core_secondword_partial_fixture_v2.py';RUNNER='scripts/run_core_secondword_partial_v2_once.py'
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';PROFILE='experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz';JOIN='experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz';FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';LOCAL='experiments/root_st073/row5_secondword_local_v4_compact_receipt.json.gz'
FAILURE='experiments/root_st073/core_secondword_partial_v1_failure.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,LOW,PROFILE,JOIN,FIRST,LOCAL,FAILURE)
PIN={FAILURE:'a229abba26f655ad15bf5e02e8bad2974c858bf9f6b9a67273b33a73fbfea709',LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PROFILE:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',JOIN:'a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',LOCAL:'70ce750f17202fc59a8f39d61eab7da864fc447be6a8ea9c79a3e625ea34ea55'}
def sha(x):return hashlib.sha256(x).hexdigest()
def main():
 raw=BRIEF.read_bytes();b=json.loads(raw);assert b['component']=='core_secondword_partial_v1' and b['phase']=='core_secondword_partial_v1_pending';assert b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True;assert not DIR.exists() and not PRE.exists()
 rows=[]
 for path in FILES:
  h=sha((ROOT/path).read_bytes());assert path not in PIN or h==PIN[path],path
  rows.append({'path':path,'sha256':h,'raw_sha256_before_decode':path in PIN})
 m={'component':NAME,'one_entry_only':True,'previous_source_checkpoint':b['source_checkpoint'],'files':rows,'policy':'No ancestor calls; fixture first; use saved 8-cell firstword profile and matrix; retain missing B1 W_ZZ, Kchi, whole tail, all errors and history.'};src=(json.dumps(m,ensure_ascii=False,separators=(',',':'))+'\n').encode()
 b['component']=NAME;b['phase']=NAME+'_pending';b['source_checkpoint']={'path':str(DIR/'source.json'),'sha256':sha(src),'preimage_path':str(PRE),'preimage_sha256':sha(raw)};b['next_task']='Run frozen second-feedback core source consumer once.';pending=(json.dumps(b,ensure_ascii=False,separators=(',',':'))+'\n').encode();assert len(pending)<=97000
 DIR.mkdir(parents=True)
 with PRE.open('xb') as f:f.write(raw)
 with (DIR/'source.json').open('xb') as f:f.write(src)
 BRIEF.write_bytes(pending);print(json.dumps({'source_sha256':sha(src),'brief_bytes':len(pending)}))
if __name__=='__main__':main()
