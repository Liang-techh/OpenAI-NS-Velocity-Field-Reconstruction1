"""Freeze new target-cell compensation-tail propagation."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_target_tail_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/kchi_ttail_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_kchi_compensation_target_tail_v1.py',
    'src/openai_ns_reconstruction/original_kchi_compensation_target_tail_fixture_v1.py',
    'scripts/run_kchi_target_tail_v1_once.py',
    'experiments/root_st073/n1_core_low_word_v3_receipt.json.gz',
    'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz',
    'experiments/root_st073/kchi_geometry_bridge_v2_receipt.json.gz',
)
PIN = {FILES[3]: '5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',
       FILES[4]: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       FILES[5]: '88e5aeecce42f34e758a429862744c096b47cd00051876495ffe1b8b11d74d84'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'kchi_geometry_bridge_v2' and brief['phase'] == 'kchi_geometry_bridge_v2_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert json.loads((ROOT / 'tmp/kchi_geometry_bridge_v2/result.json').read_bytes())['status'] == 'returned'
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen P target-tail propagation once.'
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
