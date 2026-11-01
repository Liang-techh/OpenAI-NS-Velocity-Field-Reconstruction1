"""Freeze one source-backed symbolic K0 m15 cell bridge."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_symbolic_m15_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/k0_sym_m15_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_k0_symbolic_m15_cells_v1.py',
    'src/openai_ns_reconstruction/original_k0_symbolic_m15_cells_fixture_v1.py',
    'scripts/run_k0_symbolic_m15_v1_once.py',
    'experiments/root_st073/k0_w1_saved_cells_v1_receipt.json.gz',
    'docs/CURRENT_N1_CORE_W0_M15_2026_10_30.md',
    'src/openai_ns_reconstruction/n1_core_w1_m15_adapter_v2.py',
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py',
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
)
PIN = {FILES[3]:'a212575248d4eab13a11447f578d6c45205f70cde55a5c2d47f62fe8af50d752',
       FILES[4]:'7d16e3c3b152d10301c63bb095b1e2d79e690fe7e563a104f3a6a165eae7102b',
       FILES[5]:'647dfd6e0280c9de1859d3e96c116bbe3ab92b6a83aab5d8336b64e76d1d5949',
       FILES[6]:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',
       FILES[7]:'d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'k0_w1_saved_cells_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Bind saved K0 symbolic density cells once.'
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
