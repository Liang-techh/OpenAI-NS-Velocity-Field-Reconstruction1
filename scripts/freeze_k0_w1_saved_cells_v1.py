"""Freeze a one-entry read of saved K0 finite W1 target cells."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_w1_saved_cells_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/k0_w1_cells_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_saved_k0_w1_cells_v1.py',
    'src/openai_ns_reconstruction/original_saved_k0_w1_cells_fixture_v1.py',
    'scripts/run_k0_w1_saved_cells_v1_once.py',
    'experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz',
    'experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json',
    'experiments/root_st073/scaled_core_m15_density_v1_receipt.json.gz',
)
PIN = {FILES[3]: 'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',
       FILES[4]: 'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',
       FILES[5]: 'eb0b3f0551ff2b860cfe26057741f53da41d4d925b018954435d9ece286e8a48'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'scaled_core_m15_density_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Bind the saved K0 W1 cells once.'
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
