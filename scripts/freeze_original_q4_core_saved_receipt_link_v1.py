"""Freeze a single source/receipt identity audit without rerunning providers."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_core_saved_receipt_link_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_core_receipt_link_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'scripts/run_original_q4_core_saved_receipt_link_v1_once.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    P + 'lei_ren_part1_paper_compliant_core_physical_field.py',
    P + 'lei_ren_part1_paper_compliant_current_core_common_fixed_point.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json',
    P + 'lei_ren_part1_paper_compliant_core_transfer.json',
    P + 'original_q4_frontier_closure_audit_v1_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py',
    P + 'lei_ren_part1_paper_compliant_current_core_scaled_swirl_source.py',
    P + 'lei_ren_part1_paper_compliant_core_coefficient_rebuild.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check.json',
    P + 'lei_ren_part1_paper_shared_bump_constants.py',
    P + 'lei_ren_part1_paper_shared_bump_constants.json',
    P + 'lei_ren_part1_paper_shared_analytic_tube.py',
    P + 'lei_ren_part1_paper_shared_analytic_tube.json',
    P + 'lei_ren_part1_paper_compliant_core_transfer.py',
)
PIN = {
    FILES[1]: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
    FILES[2]: '19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',
    FILES[3]: 'dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f',
    FILES[4]: 'bd9df45c31ba0224527bf21ed4fb0a62d093896ddc4662ff9e18f7f71df9d9d7',
    FILES[5]: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead',
    FILES[6]: '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',
    FILES[7]: '7ac5a829fc68db860b79013700e9b5538828a083d1f971ae3c7b0941fef36294',
    FILES[8]: 'd09f1bd0e1a55b14fd8a4006fd38f0991ea7787b88bc55ebb136c576736abf9f',
    FILES[9]: '367dc790a0b151ef5e19bc4e8db754a48b3db6e08eaaf8e7882fd48881a732a4',
    FILES[10]: 'd37bfa6fe89b8955e17b38a00e1c8ae53b6725bb5a73c4648067737ce3be990f',
    FILES[11]: '40e1a83f20642465a25936601b08e1a0237bdb7b18bbd85a31a1eca0c6971392',
    FILES[12]: '1c986f72e6967db685d8a464377d4b7b8c679a6eea2bf16aa138702f065f5ced',
    FILES[13]: '9dd940dab760c8d1976e6b543bd2170e6374400c142b8d1e2d09abdf90abf77f',
    FILES[14]: '38bde24df951cc3db6e33bae4e58c345cb6d2ca589444cf0357be1e13ba14b7b',
    FILES[15]: '1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6',
    FILES[16]: '2ad468af4815afe650eb3ebedb29aa392062c4b8a9c7b29e663540a1a2d199d6',
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_frontier_closure_audit_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen core j source owner and saved oracle receipt audit once.'
    pending = encode(brief)
    assert len(pending) <= 97000
    DIR.mkdir(parents=True)
    with PREIMAGE.open('xb') as f:
        f.write(raw)
    with (DIR / 'source.json').open('xb') as f:
        f.write(source)
    BRIEF.write_bytes(pending)
    print(json.dumps({'source_sha256': sha(source), 'brief_bytes': len(pending)}))


if __name__ == '__main__':
    main()
