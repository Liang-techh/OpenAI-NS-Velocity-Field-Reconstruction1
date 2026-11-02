"""Freeze one geometry bridge over two already saved cell packets."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'kchi_geometry_bridge_v2'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/kchi_geo_v2_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_kchi_k0_geometry_bridge_v2.py',
    'src/openai_ns_reconstruction/original_kchi_k0_geometry_bridge_fixture_v2.py',
    'scripts/run_kchi_geometry_bridge_v2_once.py',
    'experiments/root_st073/kchi_saved_cells_v2_receipt.json.gz',
    'experiments/root_st073/k0_symbolic_m15_v1_receipt.json.gz',
    'experiments/root_st073/kchi_geometry_bridge_v1_failure_receipt.json.gz',
)
PIN = {FILES[3]: '49695518f6a550540a256f56403bedda002677437b31c6eeaa559c071ae70fc3',
       FILES[4]: '288896681eae3f86559d27b1e01d8527c54621a34fa1d00da8cf2d3314eabcbf',
       FILES[5]: '6e7c00af611821ceba06d5115e3165925b1c40af00a1f05cbdfd0c59aba0a8b1'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'kchi_geometry_bridge_v1' and brief['phase'] == 'kchi_geometry_bridge_v1_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    failure = ROOT / 'experiments/root_st073/kchi_geometry_bridge_v1_failure_receipt.json.gz'
    assert sha(failure.read_bytes()) == '6e7c00af611821ceba06d5115e3165925b1c40af00a1f05cbdfd0c59aba0a8b1'
    assert json.loads((ROOT / 'tmp/kchi_geometry_bridge_v1/result.json').read_bytes())['status'] == 'failed'
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen Kchi geometry bridge once.'
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
