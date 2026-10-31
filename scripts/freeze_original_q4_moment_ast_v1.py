"""Freeze this new consumer and its inputs before its single permitted entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
DIRECTORY = ROOT / 'tmp' / 'original_q4_moment_ast_v1'
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_moment_ast_v1_brief_preimage.json')
P = 'experiments/root_st073/'
FILES = ('src/openai_ns_reconstruction/original_q4_moment_ast_v1.py',
         'src/openai_ns_reconstruction/original_q4_moment_ast_fixture_v1.py',
         'scripts/run_original_q4_moment_ast_v1_once.py',
         P + 'lei_ren_part1_paper_compliant_current_native_Rc_function_transport.py',
         P + 'lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py',
         P + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py',
         P + 'lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.json.gz',
         P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz',
         P + 'original_q4_integrand_binding_v1_receipt.json.gz')

def sha(data): return hashlib.sha256(data).hexdigest()
def encoded(obj): return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf-8')

def main():
    original = BRIEF.read_bytes(); brief = json.loads(original)
    assert brief['component'] == 'original_q4_integrand_binding_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIRECTORY.exists() and not PREIMAGE.exists()
    manifest = {'component': 'original_q4_moment_ast_v1', 'one_entry_only': True,
                'files': [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]}
    manifest_raw = encoded(manifest)
    brief['component'] = 'original_q4_moment_ast_v1'; brief['phase'] = 'original_q4_moment_ast_v1_pending'
    brief['source_checkpoint'] = {'path': str(DIRECTORY / 'source.json'), 'sha256': sha(manifest_raw),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(original)}
    brief['next_task'] = 'Execute frozen Q4 source AST and selected archive-row witness once; preserve all errors.'
    pending = encoded(brief)
    assert len(pending) <= 97000
    DIRECTORY.mkdir(parents=True)
    with PREIMAGE.open('xb') as f: f.write(original)
    with (DIRECTORY / 'source.json').open('xb') as f: f.write(manifest_raw)
    BRIEF.write_bytes(pending)
    print(json.dumps({'manifest_sha256': sha(manifest_raw), 'brief_bytes': len(pending), 'preimage_sha256': sha(original)}))

if __name__ == '__main__': main()
