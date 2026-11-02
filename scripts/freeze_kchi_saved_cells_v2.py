"""Freeze one entry for the already saved first Kchi cells and tail."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_saved_cells_v2'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/kchi_cells_v2_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_kchi_first_saved_cells_v2.py',
    'src/openai_ns_reconstruction/original_kchi_first_saved_cells_fixture_v2.py',
    'scripts/run_kchi_saved_cells_v2_once.py',
    'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz',
    'experiments/root_st073/kchi_saved_cells_v1_failure_receipt.json.gz',
)
PIN = {FILES[3]: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       FILES[4]: '172cf7a071cc0a27177eaf8cc201daf61715ed3c2264f4760eb3f5a86129fd7a'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'kchi_saved_cells_v1' and brief['phase'] == 'kchi_saved_cells_v1_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    failure = ROOT / 'experiments/root_st073/kchi_saved_cells_v1_failure_receipt.json.gz'
    assert sha(failure.read_bytes()) == '172cf7a071cc0a27177eaf8cc201daf61715ed3c2264f4760eb3f5a86129fd7a'
    assert json.loads((ROOT / 'tmp/kchi_saved_cells_v1/result.json').read_bytes())['status'] == 'failed'
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen Kchi saved-cell adapter once.'
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
