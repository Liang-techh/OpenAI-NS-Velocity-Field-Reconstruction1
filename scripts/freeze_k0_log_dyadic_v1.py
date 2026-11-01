"""Freeze a single saved log-connection and new dyadic-kernel entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'k0_log_dyadic_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/k0_log_dy_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_unbounded_dyadic_product_v1.py',
    'src/openai_ns_reconstruction/original_unbounded_dyadic_product_fixture_v1.py',
    'src/openai_ns_reconstruction/original_k0_log_connection_v1.py',
    'src/openai_ns_reconstruction/original_k0_log_connection_fixture_v1.py',
    'scripts/run_k0_log_dyadic_v1_once.py',
    'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz',
    'experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz',
    'experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json',
    'src/openai_ns_reconstruction/n1_leading_F_connection_v1.py',
    'experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json',
)
PIN = {FILES[5]:'288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       FILES[6]:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',
       FILES[7]:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',
       FILES[8]:'b3a79ea05d5b2786b4530f66b97d31558fe75c8b97f8031b42ec7364ac9198b0',
       FILES[9]:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'k0_symbolic_m15_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen log binding and dyadic product once.'
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
