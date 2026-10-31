"""Freeze one scale-preserving m15 density implementation entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'scaled_core_m15_density_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/m15_scaled_v1_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_core_scaled_m15_density_v1.py',
    'src/openai_ns_reconstruction/original_core_scaled_m15_density_fixture_v1.py',
    'scripts/run_scaled_core_m15_density_v1_once.py',
    'src/openai_ns_reconstruction/n1_core_w1_m15_adapter_v2.py',
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    'experiments/root_st073/original_q4_core_saved_receipt_link_v1_receipt.json.gz',
)
PIN = {FILES[3]: '647dfd6e0280c9de1859d3e96c116bbe3ab92b6a83aab5d8336b64e76d1d5949',
       FILES[4]: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       FILES[5]: 'd169fa93bde9ad8738abf69d3f601dab3321057b45cbee39febb8bcfb6ab957d'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_core_saved_receipt_link_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen scaled m15 density fixture once.'
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
