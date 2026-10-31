"""Freeze a new conditional endpoint-gap consumer before its sole entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_seam_radius_gap_v2'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/seam_radius_gap_v2_preimage.json')
FILES = (
    'src/openai_ns_reconstruction/original_seam_radius_gap_v2.py',
    'src/openai_ns_reconstruction/original_seam_radius_gap_fixture_v2.py',
    'scripts/run_original_seam_radius_gap_v2_once.py',
    'experiments/root_st073/original_seam_radius_gap_v1_failure_receipt.json.gz',
    'experiments/root_st073/original_f_main_gap_signed_seam_v1_receipt.json.gz',
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_seam_radius_gap_v1'
    assert brief['phase'] == 'original_seam_radius_gap_v1_pending'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert rows[-2]['sha256'] == '73e5a8e06601fde2cc01fe5beb329c94829cdf8a8598ec9b4197303208cc3b9f'
    assert rows[-1]['sha256'] == 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889'
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute the frozen conditional source-seam gap packet once; retain all errors.'
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
