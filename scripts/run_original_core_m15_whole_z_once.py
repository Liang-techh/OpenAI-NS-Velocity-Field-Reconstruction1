"""Read the pinned regular-collar sources once and preserve the full result."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_core_m15_whole_z'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
ORIGINAL = PREFIX + 'compliant_current_original_'
INPUTS = {
    'q2': (ORIGINAL + 'n1_correlated_q2_inlet.json',
           '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13'),
    'core': (PREFIX + 'compliant_core_transfer.json',
             '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497'),
    'tube': (PREFIX + 'shared_analytic_tube.json',
             '1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6'),
    'row3': (BASE + 'original_collar_cut_axial_whole_z_v2_receipt.json.gz',
             'c65abf05dbe4210c79fda560e210f0b4da127f75c091dbde44499122442307cf'),
    'row4': (BASE + 'original_collar_regular_whole_z_v2_receipt.json.gz',
             '85dae8766190ffd39c4728d35c25d4d406143e1d794ed789feecc763e04b9861'),
    'row6': (BASE + 'original_collar_whole_z_receipt.json.gz',
             '28268c3f2ef89941c02d792faf859cfbe3cbdf3b4816b7a336298fe9b20c884a'),
    'regular': (BASE + 'original_collar_regular_pressure_v2_receipt.json.gz',
                'c5ea49224dbc5f1b40e69d539e2011138ccea9b2f09551c4b8c628b48b72ab25'),
    'lower': (ORIGINAL + 'n1_lower_collar_functions.json.gz',
              'b2f983ab934875ae9aeb24b292277911dd5647e3a23943ddcf2b3b6e82595cb5'),
    'lower_code': (ORIGINAL + 'n1_lower_collar_functions.py',
                   'd00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa'),
    'cutoff': (ORIGINAL + 'n1_inner_cutoff_functions.json.gz',
               '600f4c6e69cff5b95e43312c1a9da6e77ba5e204d9873eabf27d72ba09100078'),
    'quiet': (ORIGINAL + 'compact_quiet_m15_numeric.json',
              '7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071'),
    'inner': (ORIGINAL + 'n1_inner_m15_numeric.json',
              '0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95'),
}
MODULES = (
    'scripts/run_original_core_m15_whole_z_once.py',
    'src/openai_ns_reconstruction/original_core_m15_whole_z.py',
    'src/openai_ns_reconstruction/original_collar_whole_z.py',
    'src/openai_ns_reconstruction/original_collar_cut_axial_whole_z_v2.py',
    'src/openai_ns_reconstruction/original_signed_density_expression.py',
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
        from openai_ns_reconstruction.original_core_m15_whole_z import evaluate

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
