"""Freeze a single, new conditional coefficient-grouping entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_coefficient_sum_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/k0_sum_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_bounded_dyadic_sum_v1.py',
    'src/openai_ns_reconstruction/original_bounded_dyadic_sum_fixture_v1.py',
    'src/openai_ns_reconstruction/original_k0_grouped_dyadic_v1.py',
    'src/openai_ns_reconstruction/original_k0_grouped_dyadic_fixture_v1.py',
    'scripts/run_k0_coefficient_sum_v1_once.py',
    'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz',
    'experiments/root_st073/k0_log_dyadic_v1_receipt.json.gz',
)
PIN = {FILES[5]: '288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       FILES[6]: 'b60d747b1bbd41442d31d11707e4a54c29a79161c64d5f16b2d59b9c259ba631'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'k0_log_dyadic_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen K0 coefficient sum once.'
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
