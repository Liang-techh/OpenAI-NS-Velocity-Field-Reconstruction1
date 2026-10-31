"""Freeze corrected endpoint fixture and real-IR fail-closed gate."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_interval_ir_runtime_v2'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_interval_ir_v2_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_v1.py',
    'src/openai_ns_reconstruction/original_q4_interval_ir_runtime_v2.py',
    'src/openai_ns_reconstruction/original_mpmath_directed_endpoints_v1.py',
    'src/openai_ns_reconstruction/original_mpmath_directed_endpoints_v2.py',
    'src/openai_ns_reconstruction/original_mpmath_directed_endpoints_fixture_v2.py',
    'scripts/run_original_q4_interval_ir_runtime_v2_once.py',
    P + 'original_q4_interval_ir_runtime_v1_failure_receipt.json.gz',
    P + 'original_q4_moment_lowering_v4_receipt.json.gz',
    P + 'original_q4_fraction_interval_v2_receipt.json.gz',
    P + 'original_q4_core_enclosure_adapter_v2_receipt.json.gz',
    'src/openai_ns_reconstruction/original_q4_fraction_interval_v2.py',
    'src/openai_ns_reconstruction/original_q4_core_enclosure_adapter_v2.py',
)
PIN = {FILES[6]: '66e9461771d3cf08334417846a8acca662c48b21be4c197adc6bfba21d554093',
       FILES[7]: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
       FILES[8]: 'ce463dfa4e2bee2a81eebb4ad27da224b78d50dffa5162314c6fdeb9c444120d',
       FILES[9]: '8afcf0f531fbbf458d3283ff11acc0d9160fb76abe2c6f26ad74fca73136b9ce'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_interval_ir_runtime_v1' and brief['phase'] == 'original_q4_interval_ir_runtime_v1_pending'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute corrected frozen endpoint fixture and fail-closed Q4 gate once.'
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
