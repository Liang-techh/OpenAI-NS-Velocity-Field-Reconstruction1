"""Freeze the new Q4 saved-moment bridge and its input archive hashes."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIR = ROOT / 'tmp' / 'original_q4_saved_moment_join_v1'
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_saved_moment_join_preimage.json')
P = 'experiments/root_st073/'
FILES = ('src/openai_ns_reconstruction/original_q4_saved_moment_join_v1.py',
         'src/openai_ns_reconstruction/original_q4_saved_moment_join_fixture_v1.py',
         'scripts/run_original_q4_saved_moment_join_v1_once.py',
         P + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.json.gz',
         P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz',
         P + 'lei_ren_part1_paper_compliant_current_core_interior_moments.json.gz',
         P + 'lei_ren_part1_paper_compliant_current_core_interior_moments_check.json',
         P + 'original_q4_integrand_binding_v1_receipt.json.gz',
         P + 'original_q4_moment_ast_v1_receipt.json.gz',
         P + 'original_core_physical_moment_binary_v2_receipt.json.gz')

def sha(x): return hashlib.sha256(x).hexdigest()
def encoded(x): return (json.dumps(x, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')
def main():
    raw = BRIEF.read_bytes(); brief = json.loads(raw)
    assert brief['component'] == 'original_core_physical_moment_binary_v2' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    manifest = {'component': 'original_q4_saved_moment_join_v1', 'one_entry_only': True,
                'files': [{'path': x, 'sha256': sha((ROOT / x).read_bytes())} for x in FILES]}
    mraw = encoded(manifest)
    brief['component'] = 'original_q4_saved_moment_join_v1'; brief['phase'] = 'original_q4_saved_moment_join_v1_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(mraw),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute the frozen Q4 saved-moment join once; retain all inherited errors.'
    pending = encoded(brief)
    assert len(pending) <= 97000
    DIR.mkdir(parents=True)
    with PREIMAGE.open('xb') as f: f.write(raw)
    with (DIR / 'source.json').open('xb') as f: f.write(mraw)
    BRIEF.write_bytes(pending)
    print(json.dumps({'manifest_sha256': sha(mraw), 'brief_bytes': len(pending), 'preimage_sha256': sha(raw)}))

if __name__ == '__main__': main()
