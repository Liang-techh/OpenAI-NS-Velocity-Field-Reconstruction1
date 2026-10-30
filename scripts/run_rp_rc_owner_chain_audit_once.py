"""One-entry, hash-pinned static audit; no historical producer is called."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
COMPONENT = 'rp_rc_owner_chain_audit_v1'
DIRECTORY = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / name for name in ('source.json', 'entered.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/rp_rc_owner_chain_audit.py'
RUNNER = 'scripts/run_rp_rc_owner_chain_audit_once.py'
PREFIX = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_'
SOURCE_PATHS = tuple(PREFIX + name for name in (
    'prescribed_R100_Rm_C6_functions.py', 'prescribed_loop_C5_functions.py',
    'modified_Q4_pressure_C3_functions.py', 'true_Rm_patch_C6_functions.py',
    'true_Rm_patch_finite_N_C5_functions.py', 'true_outer_C6_C5_functions.py',
    'true_Rm_patch_C6_functions_check.py',
))
REQUIRED = frozenset(SOURCE_PATHS + (LEAF, RUNNER))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def exclusive(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(raw)
    return sha(raw)


def preflight():
    if ENTRY.exists() or RESULT.exists():
        raise ValueError('component already entered')
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    if (brief.get('component') != COMPONENT or brief.get('phase') != COMPONENT + '_pending'
            or brief.get('goal_active') is not True or brief.get('full_goal_complete') is not False
            or brief.get('do_not_reenter') is not True or brief.get('one_entry_only') is not True):
        raise ValueError('pending brief guard failed')
    checkpoint = brief['source_checkpoint']
    if Path(checkpoint['path']).resolve() != SOURCE.resolve():
        raise ValueError('source checkpoint path mismatch')
    source_raw = SOURCE.read_bytes()
    if sha(source_raw) != checkpoint['sha256']:
        raise ValueError('source checkpoint hash mismatch')
    source = json.loads(source_raw)
    if (source.get('component') != COMPONENT or source.get('schema') != COMPONENT + '_checkpoint_v1'
            or source.get('one_entry_only') is not True):
        raise ValueError('source checkpoint scope mismatch')
    hashes = {}
    for item in source['files']:
        path, expected = item['path'], item['sha256']
        if path not in REQUIRED or path in hashes:
            raise ValueError('unexpected or duplicate input')
        actual = sha((ROOT / path).read_bytes())
        if actual != expected:
            raise ValueError('frozen raw hash mismatch: ' + path)
        hashes[path] = actual
    if set(hashes) != REQUIRED:
        raise ValueError('missing frozen input')
    return hashes, sha(source_raw), sha(brief_raw)


def consume(hashes):
    texts = {}
    for path in SOURCE_PATHS:
        raw = (ROOT / path).read_bytes()
        if sha(raw) != hashes[path]:
            raise ValueError('source changed after entry: ' + path)
        texts[path] = raw.decode('utf-8')
    leaf = ROOT / LEAF
    if sha(leaf.read_bytes()) != hashes[LEAF]:
        raise ValueError('leaf changed after entry')
    spec = importlib.util.spec_from_file_location('rp_rc_owner_chain_audit_v1_leaf', leaf)
    if spec is None or spec.loader is None:
        raise ValueError('leaf cannot be loaded')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if tuple(module.REQUIRED_SOURCE_PATHS) != SOURCE_PATHS:
        raise ValueError('leaf source manifest mismatch')
    return module.prove(texts)


def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {'component': COMPONENT, 'entered': True,
        'source_checkpoint_sha256': source_sha, 'brief_sha256': brief_sha,
        'input_hashes': hashes, 'consumer_calls': 0})
    result = {'component': COMPONENT, 'entry_sha256': entry_sha,
        'source_checkpoint_sha256': source_sha, 'consumer_calls': 0, 'status': 'failed'}
    try:
        result['consumer_calls'] = 1
        result['proof'] = consume(hashes)
        result['status'] = 'returned'
        result['error'] = None
    except BaseException as exc:
        result['error'] = {'type': type(exc).__name__, 'message': str(exc),
                           'traceback': traceback.format_exc()}
    result_sha = exclusive(RESULT, result)
    print(json.dumps({'component': COMPONENT, 'status': result['status'],
        'source_sha256': source_sha, 'entry_sha256': entry_sha, 'result_sha256': result_sha}))
    return 0 if result['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
