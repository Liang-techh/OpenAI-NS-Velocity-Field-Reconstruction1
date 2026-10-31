"""Single frozen link of saved core-oracle receipt and Q4 declaration."""
import gzip
import ast
import hashlib
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
NAME = 'original_q4_core_saved_receipt_link_v1'
DIR = ROOT / 'tmp' / NAME
SOURCE, ENTRY, RESULT = (DIR / p for p in ('source.json', 'entry.json', 'result.json'))
RUNNER = 'scripts/run_original_q4_core_saved_receipt_link_v1_once.py'
P = 'experiments/root_st073/'
CORE = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py'
PHYSICAL = P + 'lei_ren_part1_paper_compliant_core_physical_field.py'
COMMON = P + 'lei_ren_part1_paper_compliant_current_core_common_fixed_point.py'
CHECK = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.py'
RECEIPT = P + 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json'
TRANSFER = P + 'lei_ren_part1_paper_compliant_core_transfer.json'
PRIOR = P + 'original_q4_frontier_closure_audit_v1_receipt.json.gz'
OPERATOR = P + 'lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py'
SCALED = P + 'lei_ren_part1_paper_compliant_current_core_scaled_swirl_source.py'
REBUILD = P + 'lei_ren_part1_paper_compliant_core_coefficient_rebuild.py'
BRIDGE = P + 'lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check.json'
BUMP_PY = P + 'lei_ren_part1_paper_shared_bump_constants.py'
BUMP = P + 'lei_ren_part1_paper_shared_bump_constants.json'
TUBE_PY = P + 'lei_ren_part1_paper_shared_analytic_tube.py'
TUBE = P + 'lei_ren_part1_paper_shared_analytic_tube.json'
TRANSFER_PY = P + 'lei_ren_part1_paper_compliant_core_transfer.py'
REQUIRED = (RUNNER, CORE, PHYSICAL, COMMON, CHECK, RECEIPT, TRANSFER, PRIOR, OPERATOR, SCALED, REBUILD, BRIDGE, BUMP_PY, BUMP, TUBE_PY, TUBE, TRANSFER_PY)
PIN = {CORE:'d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59',
       PHYSICAL:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',
       COMMON:'dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f',
       CHECK:'bd9df45c31ba0224527bf21ed4fb0a62d093896ddc4662ff9e18f7f71df9d9d7',
       RECEIPT:'45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead',
       TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',
       PRIOR:'7ac5a829fc68db860b79013700e9b5538828a083d1f971ae3c7b0941fef36294',
       OPERATOR:'d09f1bd0e1a55b14fd8a4006fd38f0991ea7787b88bc55ebb136c576736abf9f',
       SCALED:'367dc790a0b151ef5e19bc4e8db754a48b3db6e08eaaf8e7882fd48881a732a4',
       REBUILD:'d37bfa6fe89b8955e17b38a00e1c8ae53b6725bb5a73c4648067737ce3be990f',
       BRIDGE:'40e1a83f20642465a25936601b08e1a0237bdb7b18bbd85a31a1eca0c6971392',
       BUMP_PY:'1c986f72e6967db685d8a464377d4b7b8c679a6eea2bf16aa138702f065f5ced',
       BUMP:'9dd940dab760c8d1976e6b543bd2170e6374400c142b8d1e2d09abdf90abf77f',
       TUBE_PY:'38bde24df951cc3db6e33bae4e58c345cb6d2ca589444cf0357be1e13ba14b7b',
       TUBE:'1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6',
       TRANSFER_PY:'2ad468af4815afe650eb3ebedb29aa392062c4b8a9c7b29e663540a1a2d199d6'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    raw = (json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
    with path.open('xb') as f:
        f.write(raw)
    return sha(raw)


def saved_evidence(raw, hashes):
    prior = json.loads(gzip.decompress(raw[PRIOR]))
    receipt = json.loads(raw[RECEIPT])
    transfer = json.loads(raw[TRANSFER])
    bridge = json.loads(raw[BRIDGE])
    bump = json.loads(raw[BUMP])
    tube = json.loads(raw[TUBE])
    rows = prior['result']['source_rows']
    assert prior['result']['status'] == 'returned' and len(rows) == 10
    family = rows['35559']['source_family']
    assert receipt['source_family'] == family
    assert receipt['input_hashes'][Path(CORE).name] == hashes[CORE]
    assert receipt['input_hashes'][Path(PHYSICAL).name] == hashes[PHYSICAL]
    assert receipt['input_hashes'][Path(COMMON).name] == hashes[COMMON]
    for path in (OPERATOR, SCALED, REBUILD):
        assert receipt['input_hashes'][Path(path).name] == hashes[path]
    assert receipt['input_hashes'][Path(CHECK).name] == hashes[CHECK]
    assert receipt['input_hashes'][Path(TRANSFER).name] == hashes[TRANSFER]
    assert transfer['input_hashes'][Path(TUBE).name] == hashes[TUBE]
    assert transfer['input_hashes'][Path(BUMP).name] == hashes[BUMP]
    assert tube['input_hashes'][Path(BUMP).name] == hashes[BUMP]
    for path in (BUMP_PY, TUBE_PY, TRANSFER_PY):
        assert receipt['input_hashes'][Path(path).name] == hashes[path]
    assert transfer['implicit_source_sha256'] == family['implicit_source_sha256']
    assert transfer['datum_enclosure_sha256'] == family['datum_enclosure_sha256']
    assert bridge['source_family'] == family
    assert bridge['input_hashes'][Path(RECEIPT).name] == hashes[RECEIPT]
    jrow = rows['35559']
    assert jrow['source_path'] == 'common.core.j' and jrow['source_provider']['module_sha256'] == hashes[CORE]
    components = [rows[str(i)] for i in range(35674,35681)]
    assert all(row['source_family'] == family and row['rho'] == 35558 and row['Z_variable'] == 20 for row in components)
    assert all(row['source_provider'] == components[0]['source_provider'] for row in components)
    actual = receipt['checks']['actual_oracle']
    identity = actual['exact_parameter_and_native_P0_identity']
    scalar = identity['actual_source_scalar_and_unique_anchor_identity']
    return {
        'source_family': family,
        'saved_check_source_hashes_match': True,
        'saved_transfer_source_and_datum_match': True,
        'q4_parameter_provider_core_hash_matches': True,
        'q4_components_declare_same_family_and_provider_descriptor': True,
        'historical_check_all_passed': receipt['all_passed'],
        'historical_actual_oracle': actual,
        'historical_identity_passed': identity['passed'],
        'historical_scalar_identity': scalar,
        'required_j_record': transfer['required_j'],
        'required_j_record_is_enclosure_not_selected_value': True,
        'saved_j_propagation': {
            'bump_required_j': bump['required_j'], 'tube_j': tube['j'],
            'transfer_j': transfer['required_j'],
            'bump_and_tube_exact_endpoint_tuples_match': all(
                bump['required_j'][k] == tube['j'][k] for k in ('lower_exact_mpf_tuple','upper_exact_mpf_tuple')),
            'tube_and_transfer_exact_endpoint_tuples_match': all(
                tube['j'][k] == transfer['required_j'][k] for k in ('lower_exact_mpf_tuple','upper_exact_mpf_tuple')),
            'bump_constants': {k:bump[k] for k in ('KN','CS','CA','CQ')},
            'source_branch_candidate_denominator': 16000000*bump['KN']*bump['CS']*bump['CA'],
            'candidate_branch_proved_by_this_audit': False,
            'interval_endpoints_not_selected_as_value': True,
        },
        'historical_inner_bridge': {
            'core_receipt_hash_matches': True,
            'all_passed': bridge['all_passed'],
            'original_numeric_n1_solution_oracle_executed': bridge['original_numeric_n1_solution_oracle_executed'],
            'checks': bridge['checks'],
            'query_results': bridge['query_results'],
            'source_identity': bridge['source_identity'],
            'positive_order_five_moment_repair_closed': bridge['positive_order_five_moment_repair_closed'],
            'full_corrected_Navier_Stokes_solution': bridge['full_corrected_Navier_Stokes_solution'],
        },
        'component_graph_binding_proved': False,
        'same_live_object_now_proved': False,
        'new_original_source_values_computed': False,
        'source_occurrence_complete_proved': False,
        'sample_scope': 'historical single rho/Z point only',
    }


def owner_chain_source_ast(raw):
    """Inspect fixed class constructor assignments, without importing them."""
    requirements = (
        (COMMON, 'CurrentCoreCommonFixedPoint', ('self.source=self.operator.source', 'self.core=self.source.core', 'self.rebuild=self.source.rebuild')),
        (OPERATOR, 'CurrentCoreNonlinearOperator', ('self.source=source if source is not None else CurrentCoreScaledSwirlSource()',)),
        (SCALED, 'CurrentCoreScaledSwirlSource', ('self.atoms=self.bridge.upstream.comparison.atoms', 'self.rebuild=self.atoms.rebuild', 'self.core=self.rebuild.core')),
        (REBUILD, 'CompliantCoreCoefficientRebuild', ('self.core=CompliantCorePhysicalField()',)),
    )
    found = []
    for path, class_name, statements in requirements:
        tree = ast.parse(raw[path].decode('utf8'))
        classes = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name]
        assert len(classes) == 1
        methods = [n for n in classes[0].body if isinstance(n, ast.FunctionDef) and n.name == '__init__']
        assert len(methods) == 1
        assignments = [n for n in methods[0].body if isinstance(n, ast.Assign)]
        for statement in statements:
            expected = ast.dump(ast.parse(statement).body[0], include_attributes=False)
            matches = [n for n in assignments if ast.dump(n, include_attributes=False) == expected]
            assert len(matches) == 1
            found.append({'file': path, 'class': class_name, 'statement': statement, 'line': matches[0].lineno})
    return {'assignments': found, 'all_expected_ast_assignments_observed': True,
            'dynamic_owner_identity_proved': False, 'source_code_executed': False}


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
                         'ancestor_called': False, 'new_original_source_values_computed': False})
    out = {'component': NAME, 'status': 'failed', 'entry_sha256': entry, 'source_sha256': sha(mraw),
           'raw_gzip_authenticated_before_decompression': True, 'ancestor_called': False,
           'new_original_source_values_computed': False, 'owner_chain_ast': None, 'saved_evidence': None}
    try:
        out['owner_chain_ast'] = owner_chain_source_ast(raw)
        out['saved_evidence'] = saved_evidence(raw, hashes)
        out['status'] = 'returned'
    except BaseException as exc:
        out['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    digest = save(RESULT, out)
    print(json.dumps({'status': out['status'],
                      'entry_sha256': entry, 'result_sha256': digest,
                      'error': out.get('error', {}).get('message')}))
    return 0 if out['status'] == 'returned' else 1


if __name__ == '__main__':
    raise SystemExit(main())
