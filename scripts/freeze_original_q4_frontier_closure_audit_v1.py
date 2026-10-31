"""Freeze one new narrow original Q4 frontier source audit."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_frontier_closure_audit_v1'
DIR = ROOT / 'tmp' / NAME
PREIMAGE = Path('C:/Users/z5242/Documents/Codex/q4_frontier_closure_v1_preimage.json')
P = 'experiments/root_st073/'
FILES = (
    'src/openai_ns_reconstruction/original_q4_frontier_closure_audit_v1.py',
    'src/openai_ns_reconstruction/original_q4_frontier_closure_audit_fixture_v1.py',
    'scripts/run_original_q4_frontier_closure_audit_v1_once.py',
    P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz',
    P + 'original_q4_integrand_binding_v1_receipt.json.gz',
    P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py',
    P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    P + 'lei_ren_part1_paper_compliant_core_physical_field.py',
)
PIN = {FILES[3]: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
       FILES[4]: '2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2',
       FILES[5]: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
       FILES[6]: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       FILES[7]: '19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf8')


def main():
    raw = BRIEF.read_bytes()
    brief = json.loads(raw)
    assert brief['component'] == 'original_q4_integral_coverage_guard_v1' and brief['phase'] == 'component_published'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert not DIR.exists() and not PREIMAGE.exists()
    rows = [{'path': p, 'sha256': sha((ROOT / p).read_bytes())} for p in FILES]
    assert all(row['sha256'] == PIN[row['path']] for row in rows if row['path'] in PIN)
    source = encode({'component': NAME, 'one_entry_only': True, 'files': rows})
    brief['component'] = NAME
    brief['phase'] = NAME + '_pending'
    brief['source_checkpoint'] = {'path': str(DIR / 'source.json'), 'sha256': sha(source),
                                  'preimage_path': str(PREIMAGE), 'preimage_sha256': sha(raw)}
    brief['next_task'] = 'Execute frozen narrow Q4 frontier source audit once.'
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
