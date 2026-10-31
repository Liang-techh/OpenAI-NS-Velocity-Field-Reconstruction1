"""One frozen entry for saved C3 endpoint Z1 source-reference transport."""
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
NAME = 'original_c3_endpoint_z1_transport_v2'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
P = 'experiments/root_st073/'
LEAF = 'src/openai_ns_reconstruction/original_c3_endpoint_z1_transport_v2.py'
FIXTURE = 'src/openai_ns_reconstruction/original_c3_endpoint_z1_transport_fixture_v2.py'
RUNNER = 'scripts/run_original_c3_endpoint_z1_transport_v2_once.py'
PRIOR = P + 'original_f_main_gap_signed_seam_v1_receipt.json.gz'
F = P + 'lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz'
C1 = P + 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py'
PULSE = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py'
C3 = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py'
C3CHECKER = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check.py'
C3CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check.json'
FAILURE = P + 'original_c3_endpoint_z1_transport_v1_failure_receipt.json.gz'
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, F, C1, PULSE, C3, C3CHECKER, C3CHECK, FAILURE)
PIN = {PRIOR: 'ebeb28cbda30149c01e4f41e6be9b1454f88b9ab8d04c58953addcb2ef131889',
       F: '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
       C3CHECK: 'cbac4685fc9031a2ed5c5a7ef624b438e1a4ce3f32b85b569512aecb1d38feb8',
       FAILURE: 'b22f5cc9ed3ac345df16290e2c4820f1ec64f8b2311f6d0f868c0e0b90d367aa'}


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


def source_function(raw, name, owner=None):
    tree = ast.parse(raw.decode('utf8'))
    scope = tree.body
    if owner:
        classes = [x for x in scope if isinstance(x, ast.ClassDef) and x.name == owner]
        if len(classes) != 1:
            raise ValueError('expected one class ' + owner)
        scope = classes[0].body
    matches = [x for x in scope if isinstance(x, ast.FunctionDef) and x.name == name]
    if len(matches) != 1:
        raise ValueError('expected one function ' + name)
    return ast.unparse(matches[0])


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
                         'ancestor_called': False, 'protected_graph_refs_indexed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_authenticated_before_decode': True, 'ancestor_called': False,
           'protected_graph_refs_indexed': False, 'fixture': None, 'source_audit': None, 'sides': {}}
    try:
        leaf = module(LEAF, NAME)
        fixture = module(FIXTURE, NAME + '_fixture')
        out['fixture'] = fixture.run_fixture(leaf.build)
        assert out['fixture']['passed'] is True
        assert all(v is True for v in out['fixture']['checks'].values())
        saved = json.loads(gzip.decompress(raw[PRIOR]))['result']
        assert saved['status'] == 'returned'
        previous = saved['archive_result']
        assert previous['status'] == 'constructed'
        archive = json.loads(gzip.decompress(raw[F]))
        assert archive['source_family'] == previous['source_family']['F']['fields']
        pins = archive['input_hashes']
        historical = json.loads(raw[C3CHECK])
        assert historical['all_passed'] is True and historical['source_family'] == archive['source_family']
        assert historical['input_hashes'][Path(C3).name] == hashes[C3]
        assert historical['input_hashes'][Path(C3CHECKER).name] == hashes[C3CHECKER]
        old_radius_check = historical['checks']['pulse']['actual_endpoint_coordinates_and_radii_Z_independent']
        assert old_radius_check is True
        historical_link = {
            'old_endpoint_coordinate_radius_Z_independent_recorded': old_radius_check,
            'old_C3_archive_hash_in_old_receipt': historical['input_hashes'].get(Path(C3).stem + '.json.gz'),
            'old_C3_archive_hash_in_F_inputs': pins.get(Path(C3).stem + '.json.gz'),
            'old_receipt_hash_in_F_inputs': pins.get(Path(C3CHECK).name),
            'old_receipt_actual_hash': hashes[C3CHECK],
            'check_scope': 'historical graph differentiation; not current F identity or physical scale transport',
        }
        source_pins = {Path(p).name: {'archived': pins.get(Path(p).name), 'actual': hashes[p],
                                    'matches': pins.get(Path(p).name) == hashes[p]} for p in (C1, PULSE, C3)}
        fixed = source_function(raw[C1], 'fixed', 'C1Algebra')
        endpoint = source_function(raw[C1], 'endpoint_functions')
        pressure = source_function(raw[C1], 'pressure_functions')
        pulse = source_function(raw[PULSE], 'pulse_pressure_functions')
        extend = source_function(raw[C3], 'extend', 'HigherZDerivative')
        c3init = source_function(raw[C3], '__init__', 'CurrentOriginalN1PulsePressureC3Functions')
        checks = {
            'fixed_sets_z_to_graph_zero': 'return source.C1Function(ref, self.g.zero)' in fixed,
            'endpoint_uses_fixed_radius': 'R = a.fixed(radius)' in endpoint,
            'pressure_uses_fixed_delta': 'c = a.fixed(delta)' in pressure,
            'pulse_endpoint_radius_fixed_before_at': 'radius = a.at(a.fixed(R), x, edge).value' in pulse,
            'pulse_endpoint_uses_saved_fields': "parent.endpoint_functions(g, at(f['U']), at(f['Q']), at(f['Q_y']), radius)" in pulse,
            'c3_extend_keeps_existing_prefix': 'existing = list(rows(value))' in extend and 'existing.append(self.derivative(existing[-1]))' in extend,
            'c3_extends_per_chart_endpoints': "self.pulse.bindings['per_chart_finite_endpoints']" in c3init and 'd.extend(value, 3)' in c3init,
        }
        out['source_audit'] = {'ast_checks': checks, 'source_file_pins': source_pins,
                               'historical_receipt': historical_link,
                               'graph_rows_indexed': False, 'source_code_executed': False,
                               'ast_shape_is_not_source_function_identity': True}
        assert all(checks.values())
        for side, chart in (('main.upper', 'main'), ('gap.lower', 'gap')):
            old = previous['candidate_old_C3_endpoints']['selected_endpoint_refs'][chart]
            refs = {}
            for letter, role in (('A', 'R_squared_q'), ('B', 'R_squared_q_squared'),
                                 ('C', 'R_squared_U_squared'), ('D', 'R_q_y')):
                refs[letter + '0'], refs[letter + '1'] = old[role][0], old[role][1]
            for role in ('Z', 'delta', 'rho_Z', 'delta_Z'):
                refs[role] = {'kind': 'unresolved_formal_port', 'role': side + '.' + role,
                              'source_value_selected': False}
            out['sides'][side] = leaf.build(refs)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'], 'fixture_passed': (out['fixture'] or {}).get('passed'),
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
