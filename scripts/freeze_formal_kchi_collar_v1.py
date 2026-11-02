"""Freeze a new higher-word formal collar consumer and its saved inputs."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'formal_kchi_collar_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/formal_kchi_pre.json')
FILES = (
    'src/openai_ns_reconstruction/original_formal_axis_kchi_volterra_v1.py',
    'src/openai_ns_reconstruction/original_formal_axis_kchi_volterra_fixture_v1.py',
    'scripts/run_formal_kchi_collar_v1_once.py',
    'experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz',
    'experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz',
    'experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz',
    'experiments/root_st073/n1_collar_core_inlet_potential_v5_receipt.json.gz',
    'experiments/root_st073/n1_collar_inlet_j2_tail_v3_receipt.json.gz',
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json',
)
PIN = {FILES[3]: '59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d',
       FILES[4]: 'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',
       FILES[5]: '5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',
       FILES[6]: '1c644fb0f41a8ae407c09f56b00e7aee1720f91d08163a0c7483b4b1b6c61e7a',
       FILES[7]: 'dd28d5be7e090afc720f7b71d48ddbeb0e542e8e4d5097e86a86d78ab8c36891',
       FILES[8]: '0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    if not __debug__:
        raise RuntimeError('optimized mode disables required preflight checks')
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'kchi_p_tail_m15_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Run frozen formal Kchi collar higher word once.'
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
