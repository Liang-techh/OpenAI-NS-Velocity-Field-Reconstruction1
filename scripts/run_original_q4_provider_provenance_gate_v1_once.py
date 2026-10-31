"""One frozen entry: static source-provider provenance gate."""
import ast
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_provider_provenance_gate_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
P = 'experiments/root_st073/'
LEAF = 'src/openai_ns_reconstruction/original_q4_provider_provenance_gate_v1.py'
FIXTURE = 'src/openai_ns_reconstruction/original_q4_provider_provenance_gate_fixture_v1.py'
RUNNER = 'scripts/run_original_q4_provider_provenance_gate_v1_once.py'
PRIOR = P + 'original_q4_moment_lowering_v4_receipt.json.gz'
ADAPTER = P + 'original_q4_core_enclosure_adapter_v2_receipt.json.gz'
PROVIDER = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json'
BINDING = P + 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py'
SWITCH = P + 'lei_ren_part1_paper_compliant_switch_signed_integrals.py'
COMMON = P + 'lei_ren_part1_paper_compliant_current_core_common_fixed_point.py'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, ADAPTER, PROVIDER, CHECK, BINDING, SWITCH, COMMON)
PIN = {PRIOR: '94698d280f97d006a8bd6c72a267c1adc18035866adb1fc046c8059545932b45',
       ADAPTER: '8afcf0f531fbbf458d3283ff11acc0d9160fb76abe2c6f26ad74fca73136b9ce',
       PROVIDER: 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       CHECK: '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead',
       BINDING: 'aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d',
       SWITCH: '713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a',
       COMMON: 'dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f'}


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
        if p in PIN:
            assert hashes[p] == PIN[p]
    entry = save(ENTRY, {'component': NAME, 'source_sha256': sha(mraw), 'brief_sha256': sha(braw),
                         'hashes': hashes, 'raw_authenticated_before_decode': True,
                         'ancestor_called': False, 'source_provider_called': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'source_provider_called': False, 'fixture': None, 'saved_descriptors': None}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf.verify, leaf.verify_wrapped)
        assert out['fixture']['passed'] is True and all(out['fixture']['checks'].values())
        prior = json.loads(gzip.decompress(raw[PRIOR]))['result']
        earlier = json.loads(gzip.decompress(raw[ADAPTER]))['result']
        check = json.loads(raw[CHECK])
        assert prior['status'] == earlier['status'] == 'returned' and check['all_passed'] is True
        assert check['input_hashes'][Path(PROVIDER).name] == hashes[PROVIDER]
        assert len(prior['lowerings']) == 13
        family = check['source_family']
        descriptors = []
        occurrences = 0
        for item in prior['lowerings']:
            lower = item['lowering']
            assert lower['status'] == 'lowered' and lower['source_family'] == family
            for atom in lower['source_atoms']:
                if atom['op'] != 'current_original_unique_analytic_core_component':
                    continue
                descriptor = atom['row']['source_provider']
                occurrences += 1
                if descriptor not in descriptors:
                    descriptors.append(descriptor)
        assert occurrences > 0 and descriptors
        expected = 'CurrentOriginalN1CoreFixedPointFunctions.component_enclosure'
        gates = [leaf.verify(d, raw[PROVIDER], Path(PROVIDER).name, expected_qualname=expected)
                 for d in descriptors]
        wrapped = [leaf.verify_wrapped(d, raw[PROVIDER], Path(PROVIDER).name,
                                      raw[SWITCH], Path(SWITCH).name, expected_qualname=expected)
                   for d in descriptors]
        binder_tree = ast.parse(raw[BINDING].decode('utf8'))
        binder = [x for x in binder_tree.body if isinstance(x, ast.FunctionDef) and x.name == 'ast_binding']
        assert len(binder) == 1
        binder_text = ast.unparse(binder[0])
        binding_checks = {
            'uses_inspect_getfile': 'inspect.getfile(fn)' in binder_text,
            'uses_source_for_ast': 'inspect.getsource(fn)' in binder_text,
            'hashes_named_module': 'sha(name)' in binder_text,
            'hashes_parsed_ast': 'ast.dump(ast.parse(text))' in binder_text,
            'common_imports_wrapper': any(isinstance(x, ast.ImportFrom)
                and x.module == Path(SWITCH).stem
                and any(alias.name == 'source_precision' and alias.asname is None for alias in x.names)
                for x in ast.parse(raw[COMMON].decode('utf8')).body),
            'provider_imports_common_alias': any(isinstance(x, ast.Import)
                and any(alias.name == Path(COMMON).stem and alias.asname == 'common' for alias in x.names)
                for x in ast.parse(raw[PROVIDER].decode('utf8')).body),
            'wrapper_imports_functools_wraps': any(isinstance(x, ast.ImportFrom)
                and x.module == 'functools'
                and any(alias.name == 'wraps' and alias.asname is None for alias in x.names)
                for x in ast.parse(raw[SWITCH].decode('utf8')).body),
        }
        assert all(binding_checks.values())
        out['saved_descriptors'] = {'source_atom_count': occurrences, 'unique_count': len(descriptors),
                                    'descriptors': descriptors, 'direct_gate_results': gates,
                                    'wrapped_gate_results': wrapped,
                                    'binding_function_ast_checks': binding_checks,
                                    'previous_adapter_descriptor': earlier['saved_source_preflight']['provider_descriptor'],
                                    'switch_source_sha256': hashes[SWITCH],
                                    'actual_callable_not_loaded_or_executed': True}
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
