"""Freeze the source-reference radius bridge before its only entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_logradius_source_binding_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/logradius_binding_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_logradius_source_binding_v1.py',
    'src/openai_ns_reconstruction/original_logradius_source_binding_fixture_v1.py',
    'scripts/run_original_logradius_source_binding_v1_once.py',
    P + 'original_f_main_gap_signed_seam_v1_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields_check.json',
)
FIXED = {FILES[3]: 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889',
         FILES[4]: '2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b',
         FILES[5]: '0bc8229a1716871d896161caac907c572f7cd0a58540900aaa00869d3a920fb1'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_c3_endpoint_z1_transport_v2' and brief['phase'] == 'component_published'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == FIXED[row['path']] for row in rows if row['path'] in FIXED)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen source-bound log-radius packet once; retain all errors.'
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
