"""Freeze saved firstword feedback into a new K0 local secondword consumer."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';NAME='row5_secondword_local_v4';DIR=ROOT/'tmp'/NAME;PRE=Path('C:/Users/z5242/Documents/Codex/row5_w2_local_v4_pre.json')
LEAF='src/openai_ns_reconstruction/row5_secondword_local_v4.py';FIXTURE='src/openai_ns_reconstruction/row5_secondword_local_fixture_v4.py';RUNNER='scripts/run_row5_secondword_local_v4_once.py'
FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';BIND='experiments/root_st073/row5_g5_binding_v1_receipt.json.gz';BASE='experiments/root_st073/row5_baseline_four_cell_v4_receipt.json.gz';Z2='experiments/root_st073/w0_z2_four_cell_v2_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
FAILURE='experiments/root_st073/row5_secondword_local_v1_failure.json.gz'
FAILURE2='experiments/root_st073/row5_secondword_local_v2_failure.json.gz'
FAILURE3='experiments/root_st073/row5_secondword_local_v3_failure.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,FIRST,BIND,BASE,Z2,SELECT,FAILURE,FAILURE2,FAILURE3)
PIN={FAILURE3:'70cca25faaea2b9ad6ff7cf3c30db8776d99d76542ef687542a0c83911fe6005',FAILURE2:'083b6cbe55af387de99177087b99400894197c575aa7036e2e5a2598bd3b7db8',FAILURE:'c8446208821b7ab6055e4538e7a344113ee834755b6a625ae5cdac2410818144',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',BIND:'8b8378de335a230bccca1268d634ce36d15d6dcdf35b35b388e8b300d3ca8309',BASE:'866fae8b0cc0471eceef78c727d7be6d52d413ffefb2679881b3100291aea8ca',Z2:'751569810106ef47969f331da12f4d7d369e6fee8e30908e58c5bbc9b7f229dc',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
def sha(x):return hashlib.sha256(x).hexdigest()
def main():
 raw=BRIEF.read_bytes();b=json.loads(raw);assert b['component']=='row5_secondword_local_v3' and b['phase']=='row5_secondword_local_v3_pending';assert b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True;assert not DIR.exists() and not PRE.exists()
 rows=[]
 for path in FILES:
  h=sha((ROOT/path).read_bytes());assert path not in PIN or h==PIN[path],path
  rows.append({'path':path,'sha256':h,'authenticate_raw_before_decode':path in PIN})
 m={'component':NAME,'one_entry_only':True,'previous_source_checkpoint':b['source_checkpoint'],'files':rows,'policy':'No ancestor calls; fixture first; source-only g5 and explicit axis correction excluded by raw K0 operator ownership; retain core inlet unknown and all prior errors.'};src=(json.dumps(m,ensure_ascii=False,separators=(',',':'))+'\n').encode()
 b['component']=NAME;b['phase']=NAME+'_pending';b['source_checkpoint']={'path':str(DIR/'source.json'),'sha256':sha(src),'preimage_path':str(PRE),'preimage_sha256':sha(raw)};b['next_task']='Run frozen local K0 secondword consumer once.';pending=(json.dumps(b,ensure_ascii=False,separators=(',',':'))+'\n').encode();assert len(pending)<=97000
 DIR.mkdir(parents=True)
 with PRE.open('xb') as f:f.write(raw)
 with (DIR/'source.json').open('xb') as f:f.write(src)
 BRIEF.write_bytes(pending);print(json.dumps({'source_sha256':sha(src),'brief_bytes':len(pending)}))
if __name__=='__main__':main()
