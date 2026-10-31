"""One frozen entry: exact interval kernel and saved IR operator coverage."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_fraction_interval_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_fraction_interval_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_fraction_interval_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_fraction_interval_v1_once.py'
PRIOR = 'experiments/root_st073/original_q4_moment_lowering_v4_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR)
PRIOR_SHA = '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f:
        f.write(raw)
    return sha(raw)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError('new component unavailable')
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def ir_operator_coverage(roots):
    """Inspect only IR expression operands, never source rows or graph IDs."""
    counts, unary = Counter(), Counter()
    pending = list(roots)
    seen_count = 0
    while pending:
        node = pending.pop()
        if not isinstance(node, dict):
            raise ValueError('IR expression node must be a mapping')
        seen_count += 1
        if seen_count > 100000:
            raise ValueError('IR expression budget exceeded')
        if node.get('kind') == 'fresh_binder':
            counts['fresh_binder'] += 1
            continue
        op = node.get('op')
        if not isinstance(op, str):
            raise ValueError('IR expression has no operator')
        counts[op] += 1
        if op == 'sum':
            pending.extend(node['terms'])
        elif op == 'product':
            pending.extend(node['factors'])
        elif op == 'negative':
            pending.append(node['argument'])
        elif op == 'positive_quotient':
            pending.extend((node['numerator'], node['denominator']))
        elif op == 'analytic_unary':
            unary[node['function']] += 1
            pending.append(node['argument'])
        elif op == 'definite_integral':
            pending.extend((node['lower'], node['upper'], node['integrand'], node['binder']))
        elif op == 'current_original_unique_analytic_core_component':
            pending.extend((node['arguments']['rho'], node['arguments']['Z_variable']))
        elif op not in ('exact_rational', 'bound_variable', 'current_original_source_parameter'):
            raise ValueError('unsupported IR operator: ' + op)
    return {'expression_nodes': seen_count, 'operator_counts': dict(counts), 'unary_counts': dict(unary),
            'unsupported_unary_names': sorted(name for name in unary if name != 'abs'),
            'source_rows_traversed': False, 'graph_nodes_indexed': False,
            'source_functions_evaluated': False, 'identity_or_error_certified': False}


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == NAME and brief['phase'] == NAME + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in manifest['files']:
        p = row['path']
        raw[p] = (ROOT / p).read_bytes()
        hashes[p] = sha(raw[p])
        assert hashes[p] == row['sha256']
    assert hashes[PRIOR] == PRIOR_SHA
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_authenticated_before_decode': True,
                         'ancestor_called': False, 'original_quadrature_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'original_quadrature_called': False, 'fixture': None, 'actual_IR_coverage': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert prior['status'] == 'returned' and len(prior['lowerings']) == 13
        roots = []
        for item in prior['lowerings']:
            assert item['lowering']['status'] == 'lowered'
            roots.append(item['lowering']['ir'])
        out['actual_IR_coverage'] = ir_operator_coverage(roots)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
