"""Read the pinned regular-collar sources once and preserve the full result."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_finite_defect'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
ORIGINAL = PREFIX + 'compliant_current_original_'
SRC = 'src/openai_ns_reconstruction/'
INPUTS = {
    'pressure': (PREFIX + 'compliant_pressure_source.json',
                 '936a30c1d5dd8df266cdf4464e5fe8cd542dbec240f593caae58bb993eea4f9e'),
    'q2': (ORIGINAL + 'n1_correlated_q2_inlet.json',
           '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13'),
    'whole': (BASE + 'original_collar_regular_whole_z_v2_receipt.json.gz',
              '85dae8766190ffd39c4728d35c25d4d406143e1d794ed789feecc763e04b9861'),
    'canonical': (BASE + 'original_collar_regular_pressure_v2_receipt.json.gz',
                  'c5ea49224dbc5f1b40e69d539e2011138ccea9b2f09551c4b8c628b48b72ab25'),
    'inner': (ORIGINAL + 'n1_inner_m15_numeric.json',
              '0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95'),
    'jets': (ORIGINAL + 'first_collar_jets.json',
             'a1a2818d18bf734c59460701d9ec4b342d6a8392ab97f581915766e9e6701a2f'),
    'core_code': (SRC + 'original_core_numeric.py',
                  '9c152f5ff50bd760ee1a860a11e48c1b3a2fd36027681ea05a3e7a0920181603'),
    'inner_code': (SRC + 'original_n1_inner_m15_numeric.py',
                   '004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b'),
    'jets_code': (SRC + 'original_first_collar_jets.py',
                  '71c424d11b6bd06f26e86c01fca2640ecdf13fcacec69253b297dc7a31972e5c'),
    'transport_code': (SRC + 'original_first_collar_transport.py',
                       '10e4552fa09c9392e4e42a1fc1a76930905121f8c6bc5ac9d96e4d9418bd4305'),
}
MODULES = (
    'scripts/run_original_collar_finite_defect_once.py',
    'src/openai_ns_reconstruction/original_collar_finite_defect.py',
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
        from openai_ns_reconstruction.original_collar_finite_defect import evaluate

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
