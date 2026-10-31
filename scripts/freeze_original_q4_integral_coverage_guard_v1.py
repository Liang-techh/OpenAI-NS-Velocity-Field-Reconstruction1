"""Freeze a single fresh guarded-integration entry and its provenance."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_integral_coverage_guard_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_coverage_guard_v1_preimage.json')
FILES = (
    'src/openai_ns_reconstruction/original_q4_integral_coverage_guard_v1.py',
    'src/openai_ns_reconstruction/original_q4_integral_coverage_guard_fixture_v1.py',
    'scripts/run_original_q4_integral_coverage_guard_v1_once.py',
    'experiments/root_st073/original_q4_interval_ir_runtime_v1_failure_receipt.json.gz',
    'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_v2.py',
    'experiments/root_st073/original_q4_interval_ir_runtime_v2_receipt.json.gz',
)
PIN = {FILES[3]: '66e9461771d3cf08334417846a8acca662c48b21be4c197adc6bfba21d554093',
       FILES[5]: '1c55c1cb4699eb8fea3e874209777455fb52467ec02d40ca8633b48dcf581c56'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_interval_ir_runtime_v2' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen guarded-integration fixture and saved geometry audit once.'
    pending = encode(brief)
    assert len(pending) <= 97000
    DIR.mkdir(parents=True)
    with PREIMAGE.open('xb') as stream:
        stream.write(raw)
    with (DIR / 'source.json').open('xb') as stream:
        stream.write(source)
    BRIEF.write_bytes(pending)
    print(json.dumps({'source_sha256': sha(source), 'brief_bytes': len(pending)}))


if __name__ == '__main__':
    main()
