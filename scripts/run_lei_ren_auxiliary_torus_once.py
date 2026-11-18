"""Run the frozen Part II auxiliary fixture once and retain its full receipt."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "lei_ren_auxiliary_torus"
FILES = {
    "src/openai_ns_reconstruction/lei_ren_auxiliary_torus.py",
    "src/openai_ns_reconstruction/charts.py",
    "scripts/check_lei_ren_auxiliary_torus_fixture.py",
    "scripts/run_lei_ren_auxiliary_torus_once.py",
}
PDF_SHA = "979e65cc4fd0b236bab04cef0284fa5808000962fa8ffa83923f3345c3519e5a"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_new(path, payload):
    with path.open('x', encoding='utf-8', newline='\n') as file:
        json.dump(payload, file, indent=2, sort_keys=True, allow_nan=False)
        file.write('\n')


def main():
    manifest_path = Path(sys.argv[1]).resolve()
    raw = manifest_path.read_bytes()
    manifest = json.loads(raw)
    expected = manifest['files']
    if (manifest['component'] != COMPONENT or manifest['one_entry_only'] is not True
            or manifest['paper_sha256'] != PDF_SHA or set(expected) != FILES):
        raise ValueError('invalid frozen Part II manifest')
    for filename, wanted in expected.items():
        if digest((ROOT / filename).read_bytes()) != wanted:
            raise ValueError(f'changed input {filename}')
    brief = Path(manifest['brief_path'])
    brief_raw = brief.read_bytes()
    state = json.loads(brief_raw)
    if state['phase'] != COMPONENT + '_pending' or state['source_checkpoint']['sha256'] != digest(raw):
        raise ValueError('pending continuation brief does not name the frozen source')
    work = manifest_path.parent
    entry = {'component': COMPONENT, 'manifest_sha256': digest(raw), 'brief_sha256': digest(brief_raw)}
    write_new(work / 'entry.json', entry)
    sys.path.insert(0, str(ROOT / 'src'))
    try:
        from check_lei_ren_auxiliary_torus_fixture import run

        result = {'passed': True, 'fixture': run()}
    except Exception as error:
        result = {'passed': False, 'error': str(error), 'traceback': traceback.format_exc()}
    receipt = {**entry, **result, 'frozen_files': expected, 'paper_sha256': PDF_SHA,
               'historical_calculations_replayed': False, 'physical_ns_validation': False}
    write_new(work / 'result.json', receipt)
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
