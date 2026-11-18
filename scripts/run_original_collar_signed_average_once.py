"""Run the frozen ST073 signed collar average once."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_signed_average'
FTC = 'experiments/root_st073/original_collar_signed_ftc_v2_receipt.json.gz'
Q2 = ('experiments/root_st073/'
      'lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json')
PINS = {
    FTC: '88c36dde29a4cc3eb38cb0800686dd6731159f20f127ff8fc2402c397331c26f',
    Q2: '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13',
}
LEAF = 'src/openai_ns_reconstruction/original_collar_signed_average.py'
RUNNER = 'scripts/run_original_collar_signed_average_once.py'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def new_json(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
        stream.write('\n')


def main():
    manifest = Path(sys.argv[1]).resolve()
    raw_manifest = manifest.read_bytes()
    source = json.loads(raw_manifest)
    files = source['files']
    if (source['component'] != COMPONENT or source['one_entry_only'] is not True
            or set(files) != {*PINS, LEAF, RUNNER}
            or any(files[name] != digest for name, digest in PINS.items())):
        raise ValueError('invalid frozen inputs')
    raw = {}
    for name, digest in files.items():
        raw[name] = (ROOT / name).read_bytes()
        if sha(raw[name]) != digest:
            raise ValueError(f'changed frozen input {name}')
    raw_brief = Path(source['brief_path']).read_bytes()
    brief = json.loads(raw_brief)
    if brief['phase'] != COMPONENT + '_pending' or brief['source_checkpoint']['sha256'] != sha(raw_manifest):
        raise ValueError('missing matching pending brief')
    out = manifest.parent
    entry = {'component': COMPONENT, 'source_sha256': sha(raw_manifest), 'brief_sha256': sha(raw_brief)}
    new_json(out / 'entry.json', entry)
    sys.path.insert(0, str(ROOT / 'src'))
    receipt = {**entry, 'frozen_inputs': files, 'raw_archives_verified_before_decode': True}
    try:
        from openai_ns_reconstruction.original_collar_signed_average import evaluate

        receipt.update(evaluate(q2_report=json.loads(raw[Q2]), ftc_receipt=json.loads(gzip.decompress(raw[FTC]))))
    except Exception as error:
        receipt.update(passed=False, error=str(error), traceback=traceback.format_exc())
    new_json(out / 'result.json', receipt)
    return 0 if receipt['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
