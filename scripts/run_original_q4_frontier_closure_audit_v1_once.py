"""Single frozen read of ten specified Q4 graph rows, no traversal."""
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_frontier_closure_audit_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
LEAF = 'src/openai_ns_reconstruction/original_q4_frontier_closure_audit_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_frontier_closure_audit_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_frontier_closure_audit_v1_once.py'
P = 'experiments/root_st073/'
ARCHIVE = P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
PRIOR = P + 'original_q4_integrand_binding_v1_receipt.json.gz'
PRODUCER = P + 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
CORE = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
PHYSICAL = P + 'lei_ren_part1_paper_compliant_core_physical_field.py'
REQUIRED = (LEAF, FIXTURE, RUNNER, ARCHIVE, PRIOR, PRODUCER, CORE, PHYSICAL)
PIN = {ARCHIVE: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
       PRIOR: '2d269c56062722fcb3e4bf4c5bf3c6968048e2cfe060ba7ca1e18207a7b9f7c2',
       PRODUCER: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
       CORE: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       PHYSICAL: '19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54'}
TARGETS = (35559, 35674, 35675, 35676, 35677, 35678, 35679, 35680, 35558, 20)
PROTECTED = frozenset((0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982))


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


def source_row_view(ref, row):
    """Expose only directly declared component handles; no child is read.

    The schema-specific fields are evidence of declarations, not a complete
    child list or a proof of lexical binding.
    """
    edges = []
    if row.get('operation') == 'current_original_unique_analytic_core_component':
        for role in ('rho', 'Z_variable'):
            value = row.get(role)
            if type(value) is not int:
                raise ValueError('component handle is not an integer: ' + role)
            edges.append({'node_id': value, 'role': role})
    return {'node_id': ref, 'metadata': row, 'declared_edges': edges,
            'declared_edges_complete': False, 'scope_aware': False}


def main():
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError('component already entered')
    braw, mraw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief, manifest = json.loads(braw), json.loads(mraw)
    assert brief['component'] == NAME and brief['phase'] == NAME + '_pending'
    assert brief['goal_active'] is True and brief['full_goal_complete'] is False
    assert brief['source_checkpoint']['sha256'] == sha(mraw)
    assert Path(brief['source_checkpoint']['path']).resolve() == SOURCE.resolve()
    assert sha(Path(brief['source_checkpoint']['preimage_path']).read_bytes()) == brief['source_checkpoint']['preimage_sha256']
    assert manifest['component'] == NAME and manifest['one_entry_only'] is True
    assert len(manifest['files']) == len(REQUIRED) and {r['path'] for r in manifest['files']} == set(REQUIRED)
    raw, hashes = {}, {}
    for row in manifest['files']:
        p = row['path']
        raw[p] = (ROOT / p).read_bytes()
        hashes[p] = sha(raw[p])
        assert hashes[p] == row['sha256']
        if p in PIN:
            assert hashes[p] == PIN[p]
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_gzip_authenticated_before_decompression': True,
                         'ancestor_called': False, 'source_values_or_quadrature_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'source_values_or_quadrature_called': False, 'fixture': None, 'source_rows': None, 'audit': None}
    try:
        assert len(set(TARGETS)) == 10 and not PROTECTED.intersection(TARGETS)
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        archive = json.loads(gzip.decompress(raw[ARCHIVE]))
        nodes = archive['exact_graph_nodes']
        assert isinstance(nodes, list)
        selected = {}
        for ref in TARGETS:
            assert type(ref) is int and ref not in PROTECTED and 0 <= ref < len(nodes)
            selected[ref] = nodes[ref]
            assert isinstance(selected[ref], dict)
        out['source_rows'] = selected
        out['audit'] = leaf.audit(selected, adapter=source_row_view)
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
