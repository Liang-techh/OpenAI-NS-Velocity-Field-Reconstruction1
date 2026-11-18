"""One frozen run of the canonical ST073 signed collar over [-1,1]."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_whole_z'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
PATHS = {
    'pressure': PREFIX + 'compliant_pressure_source.json',
    'pressure_check': PREFIX + 'compliant_current_pressure_terminal_closure_check.json',
    'core': PREFIX + 'compliant_core_transfer.json',
    'tube': PREFIX + 'shared_analytic_tube.json',
    'domain': PREFIX + 'compliant_current_original_n1_inner_analytic_domain.json.gz',
    'q2': PREFIX + 'compliant_current_original_n1_correlated_q2_inlet.json',
    'point': BASE + 'original_collar_signed_average_receipt.json',
}
PINS = {
    'pressure': '936a30c1d5dd8df266cdf4464e5fe8cd542dbec240f593caae58bb993eea4f9e',
    'pressure_check': '4c90ed7bfa7baf3808950781caae250957890d1dde05d0832e678e12ffb62806',
    'core': '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',
    'tube': '1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6',
    'domain': 'c9905f7a9d13d4f651d1f55036932ea6c69db4fc5fa890e67a9839e611af4db7',
    'q2': '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13',
    'point': '93a13f6779681539398a221770b1d49c27fbd3356c48b4675f6e9064b1515b02',
}
LEAF = 'src/openai_ns_reconstruction/original_collar_whole_z.py'
RUNNER = 'scripts/run_original_collar_whole_z_once.py'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
        stream.write('\n')


def main():
    path = Path(sys.argv[1]).resolve()
    source_bytes = path.read_bytes()
    manifest = json.loads(source_bytes)
    files = manifest['files']
    if (manifest['component'] != COMPONENT or manifest['one_entry_only'] is not True
            or manifest['inputs'] != PATHS or set(files) != {*PATHS.values(), LEAF, RUNNER}
            or any(files[PATHS[key]] != value for key, value in PINS.items())):
        raise ValueError('invalid frozen inputs')
    raw = {name: (ROOT / name).read_bytes() for name in files}
    if any(sha(blob) != files[name] for name, blob in raw.items()):
        raise ValueError('raw source digest mismatch')
    brief_bytes = Path(manifest['brief_path']).read_bytes()
    brief = json.loads(brief_bytes)
    if brief['phase'] != COMPONENT + '_pending' or brief['source_checkpoint']['sha256'] != sha(source_bytes):
        raise ValueError('missing matching pending brief')
    entry = {'component': COMPONENT, 'source_sha256': sha(source_bytes), 'brief_sha256': sha(brief_bytes)}
    save(path.parent / 'entry.json', entry)
    sys.path.insert(0, str(ROOT / 'src'))
    receipt = {**entry, 'frozen_inputs': files, 'raw_archives_verified_before_decode': True}
    try:
        from openai_ns_reconstruction.original_collar_whole_z import evaluate

        inputs = {role: json.loads(gzip.decompress(raw[name]) if name.endswith('.gz') else raw[name])
                  for role, name in PATHS.items()}
        receipt.update(evaluate(**inputs))
    except Exception as error:
        receipt.update(passed=False, error=str(error), traceback=traceback.format_exc())
    save(path.parent / 'result.json', receipt)
    return 0 if receipt['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
