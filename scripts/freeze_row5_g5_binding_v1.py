"""Freeze one new row-5 g5 occurrence binding consumer."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
NAME='row5_g5_binding_v1';DIR=ROOT/'tmp'/NAME;PRE=Path('C:/Users/z5242/Documents/Codex/row5_g5_bind_v1_pre.json')
LEAF='src/openai_ns_reconstruction/row5_g5_binding_v1.py';FIXTURE='src/openai_ns_reconstruction/row5_g5_binding_fixture_v1.py';RUNNER='scripts/run_row5_g5_binding_v1_once.py'
AFFINE='experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz';FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz';BASE='experiments/root_st073/row5_baseline_four_cell_v4_receipt.json.gz'
FILES=(LEAF,FIXTURE,RUNNER,AFFINE,FIRST,BASE)
PIN={AFFINE:'b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',BASE:'866fae8b0cc0471eceef78c727d7be6d52d413ffefb2679881b3100291aea8ca'}
def sha(x):return hashlib.sha256(x).hexdigest()
def main():
 raw=BRIEF.read_bytes();b=json.loads(raw);assert b['component']=='row5_baseline_four_cell_v4' and b['phase']=='component_published';assert b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True;assert not DIR.exists() and not PRE.exists()
 rows=[]
 for path in FILES:
  h=sha((ROOT/path).read_bytes());assert path not in PIN or h==PIN[path],path
  rows.append({'path':path,'sha256':h,'raw_sha256_before_decode':path in PIN})
 m={'component':NAME,'one_entry_only':True,'previous_source_checkpoint':b['source_checkpoint'],'files':rows,'policy':'No ancestor calls; authenticate raw gzip before decode; fixture before actual sources; compare individual affine/firstword/baseline occurrences; preserve unknown feedback.'}
 src=(json.dumps(m,ensure_ascii=False,separators=(',',':'))+'\n').encode();b['component']=NAME;b['phase']=NAME+'_pending';b['source_checkpoint']={'path':str(DIR/'source.json'),'sha256':sha(src),'preimage_path':str(PRE),'preimage_sha256':sha(raw)};b['next_task']='Run frozen row5 g5 source-occurrence binding once.';pending=(json.dumps(b,ensure_ascii=False,separators=(',',':'))+'\n').encode();assert len(pending)<=97000
 DIR.mkdir(parents=True)
 with PRE.open('xb') as f:f.write(raw)
 with (DIR/'source.json').open('xb') as f:f.write(src)
 BRIEF.write_bytes(pending);print(json.dumps({'source_sha256':sha(src),'brief_bytes':len(pending)}))
if __name__=='__main__':main()
