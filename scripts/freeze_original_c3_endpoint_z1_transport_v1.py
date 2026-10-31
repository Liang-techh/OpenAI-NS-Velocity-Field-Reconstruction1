"""Freeze the new source-reference endpoint transport before its only entry."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_c3_endpoint_z1_transport_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/c3_endpoint_z1_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_c3_endpoint_z1_transport_v1.py',
    'src/openai_ns_reconstruction/original_c3_endpoint_z1_transport_fixture_v1.py',
    'scripts/run_original_c3_endpoint_z1_transport_v1_once.py',
    P + 'original_f_main_gap_signed_seam_v1_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check.json',
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_seam_radius_gap_v2' and brief['phase'] == 'component_published'
    assert brief['goal_active'] and not brief['full_goal_complete']
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert rows[3]['sha256'] == 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889'
    assert rows[4]['sha256'] == '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen C3 endpoint Z1 transport once; retain all errors.'
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
