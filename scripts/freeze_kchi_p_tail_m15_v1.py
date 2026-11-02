"""Freeze one conditional propagation of saved P-tail errors to m15."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_p_tail_m15_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/p_tail_m15_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_kchi_p_tail_m15_v1.py',
    'src/openai_ns_reconstruction/original_kchi_p_tail_m15_fixture_v1.py',
    'scripts/run_kchi_p_tail_m15_v1_once.py',
    'experiments/root_st073/kchi_target_tail_v1_receipt.json.gz',
    'experiments/root_st073/kchi_geometry_bridge_v2_receipt.json.gz',
    'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz',
    'experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz',
    'experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json',
)
PIN = {FILES[3]: '7e50ace01f0618a39a89528b487ee1bd959d817ad5eddc0ad18a4c6032d7c9fa',
       FILES[4]: '88e5aeecce42f34e758a429862744c096b47cd00051876495ffe1b8b11d74d84',
       FILES[5]: '288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       FILES[6]: 'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',
       FILES[7]: 'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'kchi_target_tail_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen P-tail m15 error propagation once.'
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
