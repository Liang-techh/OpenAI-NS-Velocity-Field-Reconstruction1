"""Freeze the replacement consumer and preserve the failed v1 entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / 'original_core_physical_moment_binary_v2'
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/core_binary_v2_brief_preimage.json')
P = 'experiments/root_st073/'
FILES = ('src/openai_ns_reconstruction/original_core_physical_moment_binary_v2.py',
         'src/openai_ns_reconstruction/original_core_physical_moment_binary_v2_fixture_v2.py',
         'src/openai_ns_reconstruction/original_core_saved_inlet_binary_v2.py',
         'scripts/run_original_core_physical_moment_binary_v2_once.py',
         P + 'lei_ren_part1_paper_compliant_current_core_interior_moments.json.gz',
         P + 'lei_ren_part1_paper_compliant_current_core_interior_moments_check.json',
         P + 'original_core_physical_moment_scale_v1_failed_receipt.json.gz')

def sha(x): return hashlib.sha256(x).hexdigest()
def encoded(x): return (json.dumps(x, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')
def main():
    raw = BRIEF.read_bytes(); brief = json.loads(raw)
    assert brief['component'] == 'original_core_physical_moment_scale_v1'
    assert brief['phase'] == 'original_core_physical_moment_scale_v1_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    failed = (ROOT / 'tmp/original_core_physical_moment_scale_v1/result.json').read_bytes()
    assert sha(failed) == 'fd6fed257cc09ab28b7d0b4d57d971080c5928c1596e53072723f4e5be58c078'
    assert json.loads(failed)['status'] == 'failed'
    assert not DIR.exists() and not PREIMAGE.exists()
    manifest = {'component': 'original_core_physical_moment_binary_v2', 'one_entry_only': True,
                'previous_failure_result_sha256': sha(failed),
                'files': [{'path': x, 'sha256': sha((ROOT / x).read_bytes())} for x in FILES]}
    mraw = encoded(manifest)
    brief['component'] = 'original_core_physical_moment_binary_v2'
    brief['phase'] = 'original_core_physical_moment_binary_v2_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(mraw),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute new dyadic saved-moment scale consumer once; preserve v1 failure and all errors.'
    pending = encoded(brief)
    assert len(pending) <= 97000
    DIR.mkdir(parents=True)
    with PREIMAGE.open('xb') as f: f.write(raw)
    with (DIR / 'source.json').open('xb') as f: f.write(mraw)
    BRIEF.write_bytes(pending)
    print(json.dumps({'manifest_sha256': sha(mraw), 'brief_bytes': len(pending), 'preimage_sha256': sha(raw)}))

if __name__ == '__main__': main()
