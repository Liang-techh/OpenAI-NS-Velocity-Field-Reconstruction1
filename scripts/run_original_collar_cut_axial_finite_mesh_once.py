"""Read the pinned regular-collar sources once and preserve the full result."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_cut_axial_finite_mesh'
BASE = 'experiments/root_st073/'
PREFIX = BASE + 'lei_ren_part1_paper_'
ORIGINAL = PREFIX + 'compliant_current_original_'
SRC = 'src/openai_ns_reconstruction/'
INPUTS = {
    'q2': (ORIGINAL + 'n1_correlated_q2_inlet.json',
           '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13'),
    'whole': (BASE + 'original_collar_cut_axial_whole_z_v2_receipt.json.gz',
              'c65abf05dbe4210c79fda560e210f0b4da127f75c091dbde44499122442307cf'),
    'earlier_cut': (BASE + 'original_collar_cut_axial_binding_receipt.json.gz',
                    '37373941750517804af4e099ea70d49ecd596492289d39ee7bf1d76e255f59a1'),
    'forward_code': (SRC + 'original_n1_inner_m15_numeric.py',
                     '004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b'),
    'cutoff_code': (SRC + 'original_n1_cutoff_numeric.py',
                    '7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731'),
}
MODULES = (
    'scripts/run_original_collar_cut_axial_finite_mesh_once.py',
    'src/openai_ns_reconstruction/original_collar_cut_axial_finite_mesh.py',
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
        from openai_ns_reconstruction.original_collar_cut_axial_finite_mesh import evaluate

        inputs = {role: raw[path].decode('utf-8') if role.endswith('_code') else
                  json.loads(gzip.decompress(raw[path]) if path.endswith('.gz') else raw[path])
                  for role, path in aliases.items()}
        pins = {role:files[path] for role,path in aliases.items()}
        pins.update(forward_path=aliases['forward_code'],cutoff_path=aliases['cutoff_code'])
        receipt.update(evaluate(inputs=inputs,pins=pins))
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
