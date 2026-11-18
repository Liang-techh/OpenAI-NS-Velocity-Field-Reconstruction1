"""Read the pinned regular-collar sources once and preserve the full result."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_regular_whole_z_v2'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
ORIGINAL = PREFIX + 'compliant_current_original_'
INPUTS = {
    'pressure': (PREFIX + 'compliant_pressure_source.json',
                 '936a30c1d5dd8df266cdf4464e5fe8cd542dbec240f593caae58bb993eea4f9e'),
    'pressure_check': (PREFIX + 'compliant_current_pressure_terminal_closure_check.json',
                       '4c90ed7bfa7baf3808950781caae250957890d1dde05d0832e678e12ffb62806'),
    'core': (PREFIX + 'compliant_core_transfer.json',
             '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497'),
    'tube': (PREFIX + 'shared_analytic_tube.json',
             '1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6'),
    'domain': (ORIGINAL + 'n1_inner_analytic_domain.json.gz',
               'c9905f7a9d13d4f651d1f55036932ea6c69db4fc5fa890e67a9839e611af4db7'),
    'domain_check': (ORIGINAL + 'n1_inner_analytic_domain_check.json',
                     '01a4a5f497a224d3dd6ebcb7b44809df5b3702e8d32028170a764903238e0ae2'),
    'q2': (ORIGINAL + 'n1_correlated_q2_inlet.json',
           '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13'),
    'point': (BASE + 'original_collar_signed_average_receipt.json',
              '93a13f6779681539398a221770b1d49c27fbd3356c48b4675f6e9064b1515b02'),
    'whole': (BASE + 'original_collar_whole_z_receipt.json.gz',
              '28268c3f2ef89941c02d792faf859cfbe3cbdf3b4816b7a336298fe9b20c884a'),
    'regular': (BASE + 'original_collar_regular_pressure_v2_receipt.json.gz',
                'c5ea49224dbc5f1b40e69d539e2011138ccea9b2f09551c4b8c628b48b72ab25'),
    'singular': (BASE + 'original_singular_forward_compare_v2_receipt.json.gz',
                 '66510cffb63e748591ed61c4152f8628c15ee2a14e10cf527dfd5a4130d07144'),
    'inner': (ORIGINAL + 'n1_inner_m15_numeric.json',
              '0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95'),
    'quiet': (ORIGINAL + 'compact_quiet_m15_numeric.json',
              '7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071'),
}
MODULES = (
    'scripts/run_original_collar_regular_whole_z_v2_once.py',
    'src/openai_ns_reconstruction/original_collar_regular_whole_z_v2.py',
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
        from openai_ns_reconstruction.original_collar_regular_whole_z_v2 import evaluate

        inputs = {role: json.loads(gzip.decompress(raw[path]) if path.endswith('.gz') else raw[path])
                  for role, path in aliases.items()}
        receipt.update(evaluate(**inputs))
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
