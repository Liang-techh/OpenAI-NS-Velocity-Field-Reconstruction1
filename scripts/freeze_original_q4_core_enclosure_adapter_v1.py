"""Freeze the core-enclosure adapter before its only entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_core_enclosure_adapter_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_core_enclosure_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_v1.py',
    'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_fixture_v1.py',
    'scripts/run_original_q4_core_enclosure_adapter_v1_once.py',
    P + 'original_q4_moment_lowering_v4_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json',
)
FIXED = {FILES[3]: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
         FILES[4]: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
         FILES[5]: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_scoped_ir_evaluator_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == FIXED[row['path']] for row in rows if row['path'] in FIXED)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen enclosure adapter once; no original provider or recurrence.'
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
