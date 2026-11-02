"""Freeze a new baseline W0 row-5 source consumer, without entering it."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
NAME='row5_baseline_four_cell_v4'
DIR=ROOT/'tmp'/NAME
PRE=Path('C:/Users/z5242/Documents/Codex/row5_w0_v4_pre.json')
LEAF='src/openai_ns_reconstruction/w0_joined_jets_v4.py'
FIXTURE='src/openai_ns_reconstruction/w0_joined_jets_fixture_v4.py'
PRODUCT='src/openai_ns_reconstruction/row5_baseline_product_v4.py'
PRODUCT_FIXTURE='src/openai_ns_reconstruction/row5_baseline_product_fixture_v4.py'
RUNNER='scripts/run_row5_baseline_four_cell_v4_once.py'
OLD='experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz'
Z2='experiments/root_st073/w0_z2_four_cell_v2_receipt.json.gz'
AFFINE='experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz'
FIRST='experiments/root_st073/row5_firstword_product_v4_receipt.json.gz'
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz'
FAILURE='experiments/root_st073/row5_baseline_four_cell_v1_failure.json.gz'
FAILURE2='experiments/root_st073/row5_baseline_four_cell_v2_failure.json.gz'
FAILURE3='experiments/root_st073/row5_baseline_four_cell_v3_failure.json.gz'
FILES=(LEAF,FIXTURE,PRODUCT,PRODUCT_FIXTURE,RUNNER,OLD,Z2,AFFINE,FIRST,FORCING,FAILURE,FAILURE2,FAILURE3)
PIN={FAILURE3:'34dccd5503648d2ab525bb5ea00594fca6bed56b20b3082fd602adcc2a407edd',FAILURE2:'0615f5a14f8d3690c8cf08d1b98fed91656c6fa1e4892b16ebabbae1c3840069',FAILURE:'01359e0b7b2b9bd298ccff1bb34a5b568bc0fd88be9f2f4dc8da4ccf449721bf',OLD:'4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06',Z2:'751569810106ef47969f331da12f4d7d369e6fee8e30908e58c5bbc9b7f229dc',AFFINE:'b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515',FIRST:'1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2'}
def sha(x):return hashlib.sha256(x).hexdigest()
def main():
 raw=BRIEF.read_bytes();b=json.loads(raw)
 assert b['component']=='row5_baseline_four_cell_v3' and b['phase']=='row5_baseline_four_cell_v3_pending'
 assert b['goal_active'] is True and b['full_goal_complete'] is False and b['one_entry_only'] is True and b['do_not_reenter'] is True
 assert not DIR.exists() and not PRE.exists()
 rows=[]
 for path in FILES:
  h=sha((ROOT/path).read_bytes())
  if path in PIN:assert h==PIN[path],path
  rows.append({'path':path,'sha256':h,'raw_sha256_only':path in (FIRST,FORCING,FAILURE,FAILURE2,FAILURE3)})
 manifest={'component':NAME,'one_entry_only':True,'previous_component':b['component'],'previous_source_checkpoint':b['source_checkpoint'],'files':rows,'policy':'No ancestor imports or calls. Authenticate raw gzip before decode; firstword receipt raw-only. Synthetic fixtures before receipt decode. Preserve g5 once and unknown remainder.'}
 source=(json.dumps(manifest,ensure_ascii=False,separators=(',',':'))+'\n').encode()
 b['component']=NAME;b['phase']=NAME+'_pending';b['source_checkpoint']={'path':str(DIR/'source.json'),'sha256':sha(source),'preimage_path':str(PRE),'preimage_sha256':sha(raw)};b['next_task']='Run frozen baseline row5 four-cell consumer once.'
 pending=(json.dumps(b,ensure_ascii=False,separators=(',',':'))+'\n').encode();assert len(pending)<=97000
 DIR.mkdir(parents=True)
 with PRE.open('xb') as f:f.write(raw)
 with (DIR/'source.json').open('xb') as f:f.write(source)
 BRIEF.write_bytes(pending)
 print(json.dumps({'source_sha256':sha(source),'brief_bytes':len(pending)}))
if __name__=='__main__':main()
