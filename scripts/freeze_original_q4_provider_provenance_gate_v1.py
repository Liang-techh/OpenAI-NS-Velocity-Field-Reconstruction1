"""Freeze the new static provider provenance gate before its only entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_provider_provenance_gate_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_provider_gate_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_q4_provider_provenance_gate_v1.py',
    'src/openai_ns_reconstruction/original_q4_provider_provenance_gate_fixture_v1.py',
    'scripts/run_original_q4_provider_provenance_gate_v1_once.py',
    P + 'original_q4_moment_lowering_v4_receipt.json.gz',
    P + 'original_q4_core_enclosure_adapter_v2_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json',
    P + 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py',
    P + 'lei_ren_part1_paper_compliant_switch_signed_integrals.py',
    P + 'lei_ren_part1_paper_compliant_current_core_common_fixed_point.py',
)
PIN = {FILES[3]: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
       FILES[4]: '8afcf0f531fbbf458d3283ff11acc0d9160fb76abe2c6f26ad74fca73136b9ce',
       FILES[5]: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       FILES[6]: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead',
       FILES[7]: 'aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d',
       FILES[8]: '713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a',
       FILES[9]: 'dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_core_enclosure_adapter_v2' and brief['phase'] == 'component_published'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen callable provenance gate once; retain all errors.'
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
