"""Freeze the new scoped-lowering consumer before its single entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / 'original_q4_moment_lowering_v4'
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_moment_lowering_v4_preimage.json')
P = 'experiments/root_st073/'
FILES = ('src/openai_ns_reconstruction/original_q4_moment_lowering_v4.py',
         'src/openai_ns_reconstruction/original_q4_moment_lowering_fixture_v4.py',
         'scripts/run_original_q4_moment_lowering_v4_once.py',
         P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz',
         P + 'original_q4_saved_moment_join_v1_receipt.json.gz',
         P + 'original_q4_moment_lowering_v3_unresolved_receipt.json.gz')

def sha(x): return hashlib.sha256(x).hexdigest()
def encoded(x): return (json.dumps(x, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')
def main():
    raw = BRIEF.read_bytes(); brief = json.loads(raw)
    assert brief['component'] == 'original_q4_moment_lowering_v3' and brief['phase'] == 'original_q4_moment_lowering_v3_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    failed = (ROOT / 'tmp/original_q4_moment_lowering_v3/result.json').read_bytes()
    assert sha(failed) == 'a3212693b4de21888311327d6c8785cb331acc46401103eec34eb62c6c800d3d'
    prior = json.loads(failed)
    assert prior['status'] == 'returned' and len(prior['lowerings']) == 13
    assert all(item['lowering']['status'] == 'unresolved' for item in prior['lowerings'])
    assert not DIR.exists() and not PREIMAGE.exists()
    manifest = {'component': 'original_q4_moment_lowering_v4', 'one_entry_only': True,
                'previous_failure_result_sha256': sha(failed),
                'files': [{'path': x, 'sha256': sha((ROOT / x).read_bytes())} for x in FILES]}
    mraw = encoded(manifest)
    brief['component'] = 'original_q4_moment_lowering_v4'; brief['phase'] = 'original_q4_moment_lowering_v4_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(mraw),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute v4 Q4 scoped lowering once; preserve prior attempts and all inherited errors.'
    pending = encoded(brief)
    assert len(pending) <= 97000
    DIR.mkdir(parents=True)
    with PREIMAGE.open('xb') as f: f.write(raw)
    with (DIR / 'source.json').open('xb') as f: f.write(mraw)
    BRIEF.write_bytes(pending)
    print(json.dumps({'manifest_sha256': sha(mraw), 'brief_bytes': len(pending), 'preimage_sha256': sha(raw)}))

if __name__ == '__main__': main()
