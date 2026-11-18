"""Read the pinned regular-collar sources once and preserve the full result."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_finite_collar_inlet_bound'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
ORIGINAL = PREFIX + 'compliant_current_original_'
SRC = 'src/openai_ns_reconstruction/'
INPUTS = {
    'pressure': (PREFIX + 'compliant_pressure_source.json',
                 '936a30c1d5dd8df266cdf4464e5fe8cd542dbec240f593caae58bb993eea4f9e'),
    'core': (PREFIX + 'compliant_core_transfer.json',
             '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497'),
    'tube': (PREFIX + 'shared_analytic_tube.json',
             '1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6'),
    'q2': (ORIGINAL + 'n1_correlated_q2_inlet.json',
           '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13'),
    'whole': (BASE + 'original_collar_regular_whole_z_v2_receipt.json.gz',
              '85dae8766190ffd39c4728d35c25d4d406143e1d794ed789feecc763e04b9861'),
    'finite': (BASE + 'original_collar_finite_defect_receipt.json.gz',
               'b425948f2301a86cd481b43cd9a771f89b625431a8862f47ef24eeb2a9a64e74'),
    'step_code': (PREFIX + 'logarithmic_core_step.py',
                  '52210fcdde9a538ca28a7010264daa6802b1cd02052741e05001162617d34fb6'),
    'initial_code': (PREFIX + 'functional_core_step.py',
                     '03aab21f3e6e665dd2ff0e9a3d2e6187e90813c300ccbb2d35e23a6d99098998'),
    'pressure_code': (PREFIX + 'candidate_pressure_function.py',
                      '4cc4c87488e5f468aa9a459d1e8c9f754e85f12dee811ca9d8e066fd4d4986dc'),
}
MODULES = (
    'scripts/run_original_finite_collar_inlet_bound_once.py',
    'src/openai_ns_reconstruction/original_finite_collar_inlet_bound.py',
    'src/openai_ns_reconstruction/original_collar_whole_z.py',
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
        stream.write('\n')


def main():
    manifest_path = Path(sys.argv[1]).resolve()
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    aliases = {role: path for role, (path, _) in INPUTS.items()}
    pins = {path: digest for path,digest in INPUTS.values()}
    files = manifest['files']
    if (manifest['component'] != COMPONENT or manifest['one_entry_only'] is not True
            or manifest['inputs'] != aliases or set(files) != set(pins) | set(MODULES)
            or any(files[path] != digest for path,digest in pins.items())):
        raise ValueError('wrong frozen source or file set')
    raw = {path: (ROOT/path).read_bytes() for path in files}
    if any(files[path] != sha(payload) for path,payload in raw.items()):
        raise ValueError('source SHA-256 changed')

    brief_path = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
    brief_bytes = brief_path.read_bytes()
    brief = json.loads(brief_bytes)
    manifest_sha = sha(manifest_bytes)
    if (brief['phase'] != COMPONENT+'_pending'
            or brief['source_checkpoint']['sha256'] != manifest_sha):
        raise ValueError('missing pending brief for this manifest')

    entry = dict(component=COMPONENT, source_sha256=manifest_sha, brief_sha256=sha(brief_bytes),
                 frozen_inputs=files, raw_archives_verified_before_decode=True)
    save(manifest_path.parent/'entry.json', entry)
    receipt = dict(entry)
    sys.path.insert(0, str(ROOT/'src'))
    try:
        from openai_ns_reconstruction.original_finite_collar_inlet_bound import evaluate

        inputs = {role: raw[path].decode('utf-8') if role.endswith('_code') else
                  json.loads(gzip.decompress(raw[path]) if path.endswith('.gz') else raw[path])
                  for role, path in aliases.items()}
        receipt.update(evaluate(inputs=inputs,pins={role:files[path] for role,path in aliases.items()}))
    except Exception as error:
        receipt.update(passed=False, error=str(error), traceback=traceback.format_exc())
    result = manifest_path.parent/'result.json'
    save(result, receipt)
    destination = ROOT/BASE/(COMPONENT+'_receipt.json.gz')
    with destination.open('xb') as output:
        output.write(gzip.compress(result.read_bytes(), mtime=0))
    return 0 if receipt.get('passed') else 1


if __name__ == '__main__':
    sys.exit(main())
