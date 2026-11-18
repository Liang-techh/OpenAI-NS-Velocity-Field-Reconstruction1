"""One frozen run of the canonical regular-pressure source consumer."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_regular_pressure_v2'
LEAF = 'src/openai_ns_reconstruction/original_collar_regular_pressure_v2.py'
RUNNER = 'scripts/run_original_collar_regular_pressure_v2_once.py'
POINT = 'src/openai_ns_reconstruction/original_collar_regular_pressure_point.py'
GRAPH = 'src/openai_ns_reconstruction/original_signed_density_expression.py'
ROLES = frozenset({'lower', 'cutoff', 'cut_check', 'ftc', 'q2', 'forward', 'native', 'forward_code'})


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _repo_path(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError(f'{label} must be a non-empty path string')
    candidate = Path(value)
    resolved = (candidate if candidate.is_absolute() else ROOT / candidate).resolve()
    try:
        resolved.relative_to(ROOT)
    except ValueError as error:
        raise ValueError(f'{label} is outside the repository: {value}') from error
    return resolved


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
        stream.write('\n')


def _decode(role, name, raw):
    payload = gzip.decompress(raw) if name.endswith('.gz') else raw
    if role.endswith('_code'):
        return payload.decode('utf-8')
    return json.loads(payload.decode('utf-8'))


def main():
    if len(sys.argv) != 2:
        raise ValueError('expected one source manifest path')
    manifest_path = Path(sys.argv[1]).resolve()
    _repo_path(str(manifest_path), 'source manifest')
    manifest_raw = manifest_path.read_bytes()
    manifest = json.loads(manifest_raw)
    if (manifest.get('component') != COMPONENT
            or manifest.get('one_entry_only') is not True):
        raise ValueError('invalid frozen manifest identity')

    if manifest.get('component') != COMPONENT or manifest.get('one_entry_only') is not True:
        raise ValueError('wrong frozen component')
    inputs = manifest.get('inputs')
    if not isinstance(inputs, dict) or set(inputs) != ROLES:
        raise ValueError('invalid frozen input roles')
    for role, name in inputs.items():
        _repo_path(name, f'input {role}')
    input_paths = set(inputs.values())

    files = manifest.get('files')
    if not isinstance(files, dict) or set(files) != input_paths | {LEAF, RUNNER, POINT, GRAPH}:
        raise ValueError('invalid frozen file set')
    for name, expected in files.items():
        _repo_path(name, f'frozen file {name}')
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError(f'invalid SHA-256 for {name}')

    # Authenticate every listed file as raw bytes before writing the exclusive
    # entry receipt, importing the consumer, or decoding any input archive.
    raw = {name: _repo_path(name, f'frozen file {name}').read_bytes() for name in files}
    for name, blob in raw.items():
        if sha(blob) != files[name]:
            raise ValueError(f'raw source digest mismatch: {name}')

    brief_path = Path(manifest['brief_path']).resolve()
    if brief_path != ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json':
        raise ValueError('wrong authoritative pending brief')
    brief_bytes = brief_path.read_bytes()
    brief = json.loads(brief_bytes)
    if (brief.get('phase') != COMPONENT + '_pending'
            or brief.get('source_checkpoint', {}).get('sha256') != sha(manifest_raw)):
        raise ValueError('missing matching pending brief')

    entry = {
        'component': COMPONENT,
        'source_sha256': sha(manifest_raw),
        'manifest_sha256': sha(manifest_raw),
        'brief_sha256': sha(brief_bytes),
        'phase': COMPONENT + '_pending',
        'frozen_inputs': files,
        'raw_archives_verified_before_decode': True,
    }
    save(manifest_path.parent / 'entry.json', entry)

    sys.path.insert(0, str(ROOT / 'src'))
    receipt = dict(entry)
    try:
        from openai_ns_reconstruction.original_collar_regular_pressure_v2 import evaluate

        decoded = {role: _decode(role, name, raw[name]) for role, name in inputs.items()}
        receipt.update(evaluate(inputs=decoded, pins={role: files[name] for role,name in inputs.items()}))
    except Exception as error:
        receipt.update(passed=False, error=str(error), traceback=traceback.format_exc())
    save(manifest_path.parent / 'result.json', receipt)
    return 0 if receipt.get('passed') else 1


if __name__ == '__main__':
    sys.exit(main())
