"""Source-only initial inlet projection for the original first micro chart.

This module authenticates already accepted source and receipt records and
checks their AST operands.  It deliberately never imports a producer, reads
the workspace, evaluates a node, or replays an interval/operator.  Its one
mathematical conclusion is the source projection at ``s=0`` (``rho=4``): the
actual bridge packet, C2 initial consumer, C6 inlet/bridge consumers, and the
common fixed-point six-integral recipes all consume the same analytic core
functions.  Tail terms remain enclosures throughout the projection.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping


SCHEMA = "original_first_micro_initial_source_identity_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
POSITIVE_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_native_transport_pressure_sources.json"
POSITIVE_SHA256 = "7c88052348a78ab7e26728f492296fb57858cead1794b9d8cb6e57ad873a7b9b"


# Every record is a basename in experiments/root_st073.  These are all
# inherited pins from the two accepted reports; no new archive is admitted.
RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    POSITIVE_REPORT: POSITIVE_SHA256,
    "lei_ren_part1_paper_compliant_actual_bridge_integrals.py": "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3",
    "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": "575fc6b5c1c78957bf9e0d178582bfcff26c167d4de037cd7aea41fc69f52145",
    "lei_ren_part1_paper_compliant_comparison_point_integrals.py": "a7103469587709ff80b087b99fcf5490c8445afcd25c2794a4097b73963d98b7",
    "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": "3b2cd05bd8732acaf3b53b4f97422c55721a4ab9db9c760758b5cbdd0e31171b",
    "lei_ren_part1_paper_compliant_core_integral_atoms.py": "1c24860d3d4f6b070731b251f9da77a31d3096a7d917a4107c3140ed127f6261",
    "lei_ren_part1_paper_compliant_core_integral_atoms_check.json": "7a8be1e68a75a3c98b5498484ca501204c2c2df80659c5bf490957f05cfbd6ef",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py": "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json": "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672",
    "lei_ren_part1_paper_compliant_current_core_first_interface.py": "3785967f555cb5669ae227cda4883dab7f084bd02ca87d6742462d78605627ff",
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.json": "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400",
    "lei_ren_part1_paper_compliant_current_core_interior_moments.py": "9c45b35c9fe1ef6c4eb54738d7c6323b1790b77d5ebda95fd052544ed1f05eaf",
    "lei_ren_part1_paper_compliant_current_core_interior_moments_check.json": "cf1297a115b12cee236c97a8b837eac4d57e1e880d87c104170056e61cad70b9",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py": "d09f1bd0e1a55b14fd8a4006fd38f0991ea7787b88bc55ebb136c576736abf9f",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json": "a1f9a3632a7765f6bc8d38ac1bd25c28850e128a8b25adaa0d42f769f3704d62",
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild.py": "d37bfa6fe89b8955e17b38a00e1c8ae53b6725bb5a73c4648067737ce3be990f",
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild_check.json": "02790017aeddc309624dcff46f2d4a2cd9d557075cc0e88096ff0c63cd9b65d0",
    "lei_ren_part1_paper_compliant_rooted_core_field.py": "d1830d0779bd0c83a7232596726a395d744168dc30ccb657a0edca6642bdc5e4",
    "lei_ren_part1_paper_compliant_rooted_core_field_check.json": "9e8b15a0cd4fedf4e7e9580a3caa8a39092f5c2cfab323dc6bc297ce37ccb31c",
    "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py": "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
    "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source_check.json": "9d8b746ef78502956cbb2b53f96d452be636eecbf645d021aaca3f9c2e21885d",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    "lei_ren_part1_paper_compliant_current_original_micro_functions.py": "d4ef2385ff4627c91f3fd99d1c6819e5598db3e789093f02fe1dc965ec8be779",
    "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": "54d14670660ea41e90ffc25f536e65747da2a3008d8c34e790780cd2d9df8fac",
}

DECODED_RECORD_HASHES = {}


def _digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True, allow_nan=False).encode("utf-8")).hexdigest()


def _get(value, path, default=None):
    current = value
    for part in path.split(".") if path else ():
        if not isinstance(current, Mapping) or part not in current:
            return default
        current = current[part]
    return current


def _record_text(records, name, gates):
    record = records.get(name)
    if not isinstance(record, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        gates.append("source_path_mismatch:" + name)
        return ""
    if record.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
        return ""
    try:
        if record.get("encoding") == "gzip":
            raw = base64.b64decode(record["raw_bytes_base64"], validate=True)
            if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
                gates.append("source_raw_hash_mismatch:" + name)
                return ""
            decoded = gzip.decompress(raw)
            declared = record.get("decoded_sha256")
            known = DECODED_RECORD_HASHES.get(name, declared)
            if known and hashlib.sha256(decoded).hexdigest() != known:
                gates.append("source_decoded_hash_mismatch:" + name)
            return decoded.decode("utf-8")
        text = record.get("source_text")
        if not isinstance(text, str):
            raise TypeError("plain source_text missing")
        if hashlib.sha256(text.encode("utf-8")).hexdigest() != RECORD_HASHES[name]:
            gates.append("source_hash_mismatch:" + name)
            return ""
        return text
    except Exception as error:  # keep a source-only negative result reviewable
        gates.append("source_decode_error:" + name + ":" + type(error).__name__)
        return ""


def _tree(text, name, gates):
    try:
        return ast.parse(text, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _json(texts, name, gates):
    try:
        return json.loads(texts[name])
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _function(tree, name, method=None):
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    if method is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == name:
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method:
                        return child
    return None


def _class_method(tree, cls, method):
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method:
                    return child
    return None


def _same_expr(node, expression):
    try:
        expected = ast.parse(expression, mode="eval").body
        return ast.dump(node, include_attributes=False) == ast.dump(expected, include_attributes=False)
    except (SyntaxError, TypeError):
        return False


def _assignments(fn, target):
    if fn is None:
        return []
    try:
        target_text = "".join(ast.unparse(ast.parse(target, mode="eval").body).split())
    except (SyntaxError, TypeError):
        target_text = target
    values = []
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign):
            if any("".join(ast.unparse(item).split()) == target_text
                   for item in node.targets):
                values.append(node.value)
        elif isinstance(node, ast.AnnAssign) and ast.unparse(node.target) == target:
            values.append(node.value)
    return values


def _assignment(fn, target, expression):
    return any(_same_expr(value, expression) for value in _assignments(fn, target))


def _keyword(fn, call_name, keyword, expression):
    if fn is None:
        return False
    return any(
        isinstance(node, ast.Call) and ast.unparse(node.func) == call_name and
        any(item.arg == keyword and _same_expr(item.value, expression)
            for item in node.keywords)
        for node in ast.walk(fn))


def _contains_call(fn, dotted):
    return fn is not None and any(
        isinstance(node, ast.Call) and ast.unparse(node.func) == dotted
        for node in ast.walk(fn))


def _contains_text(text, *needles):
    return all(needle in text for needle in needles)


def _receipt(receipts, name, gates, paths=()):
    value = receipts.get(name, {})
    ok = isinstance(value, Mapping) and value.get("all_passed") is True
    if not ok:
        gates.append("receipt_not_accepted:" + name)
    for path in paths:
        actual = _get(value, path)
        expected_values = {
            "source_bindings_and_incidence.actual_production_gauge_differences": ["0", "0"],
            "checks.source_rows.preserved_true_core_C5_sources": ["Phi", "Psi", "mean_Psi", "U0", "gradient", "chi", "G", "logF0", "F0"],
            "checks.source_rows.actual_seven_row_complete_moment_histories": 6,
        }
        expected = expected_values.get(path, True)
        if isinstance(expected, list):
            matches = isinstance(actual, list) and actual == expected
        elif type(expected) is int:
            matches = type(actual) is int and actual == expected
        else:
            matches = actual is True
        if not matches:
            gates.append("receipt_missing:" + name + ":" + path)
            ok = False
    return ok


def _source_record_set(source_records):
    return source_records if isinstance(source_records, Mapping) else {}


def _check_reports(texts, gates):
    parent = _json(texts, PARENT_REPORT, gates)
    positive = _json(texts, POSITIVE_REPORT, gates)
    if _get(positive, "parent_report_sha256") != PARENT_SHA256:
        gates.append("positive_parent_pin_changed")
    scope = _get(positive, "scope", {})
    required_scope = (
        "common_independent_normalized_P0_source_bound",
        "old_enclosures_same_unique_analytic_core_projection_bound",
        "C2_C5_native_m_transport_definitions_checked",
    )
    for path in required_scope:
        if scope.get(path) is not True:
            gates.append("positive_scope_missing:" + path)
    if scope.get("full_first_micro_cross_owner_source_identity_proved") is not False:
        gates.append("positive_scope_overclaims_full_micro")
    if scope.get("full_goal_complete") is not False:
        gates.append("positive_scope_overclaims_goal")
    native = _get(positive, "native_m_transport_definitions", {})
    for path in ("native_transport_definitions_checked",):
        if native.get(path) is not True:
            gates.append("native_definitions_missing:" + path)
    for path in ("nodes_digest",):
        if native.get(path) != EXPECTED_NODES_DIGEST:
            gates.append("native_nodes_digest_changed")
    first = _get(positive, "first_micro_source_identity", {})
    if first.get("nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("first_micro_nodes_digest_changed")
    for path in (
        "component_proofs.conditional_formula_witness",
        "component_proofs.common_core_projection_sublemma.passed",
        "component_proofs.common_pressure_P0_sublemma.passed",
        "component_proofs.c2_actual_source.passed",
        "component_proofs.prescribed_C6_actual_micro.passed",
    ):
        if _get(first, path) is not True:
            gates.append("positive_first_micro_missing:" + path)
    rows = first.get("rows") if isinstance(first, Mapping) else None
    if not (isinstance(rows, list) and len(rows) == 1 and isinstance(rows[0], Mapping)
            and rows[0].get("chart") == "first_micro"
            and rows[0].get("same_actual_source_identity") is False):
        gates.append("positive_first_micro_overclaims_identity")
    return parent, positive


def _check_receipts(receipts, gates):
    specs = {
        "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": (
            "full_finite_width_signed_bridge_integral_enclosures_available",
            "actual_own_six_moment_feedback_enclosures_available",
        ),
        "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": (
            "original_9_23_namespace_bound",
            "actual_atoms_and_formal_width_remainders_preserved",
        ),
        "lei_ren_part1_paper_compliant_core_integral_atoms_check.json": (),
        "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json": (
            "current_core_actual_operator_Banach_fixed_point_admitted",
            "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point",
            "current_core_actual_operator_contraction_and_tail_certified",
            "actual_Banach_and_common_coefficient_proof.actual_source_operator_defined_on_same_complete_Xh_pair",
            "actual_Banach_and_common_coefficient_proof.unique_actual_analytic_fixed_point_for_each_fixed_implicit_source",
            "actual_Banach_and_common_coefficient_proof.all_original_finite_rows_enclose_coefficients_of_this_unique_solution",
            "actual_Banach_and_common_coefficient_proof.no_interval_overlap_or_fitted_source_equality_argument",
            "actual_Banach_and_common_coefficient_proof.same_source_analytic_functions_not_just_shared_receipt_ids",
            "same_solution_all_order_tail_proof.same_actual_fixed_point_coefficient_norm_used",
            "same_solution_all_order_tail_proof.finite_and_infinite_parts_of_same_analytic_solution",
            "same_solution_all_order_tail_proof.ordinary_derivative_factorials_and_radial_units_preserved",
        ),
        "lei_ren_part1_paper_compliant_current_core_first_interface_check.json": (
            "current_core_bridge_functional_mixed4_join_certified",
            "current_common_core_bridge_graph.common_rebuild_is_actual_atom_rebuild",
            "current_common_core_bridge_graph.actual_comparison_uses_same_atoms",
            "original_boundary_source_bindings.actual_first_initial_values_and_moments_are_atoms_not_covers",
            "original_boundary_source_bindings.exact_first_zero_source_increment_branch_bound",
            "original_boundary_source_bindings.original_covering_integral_is_only_an_enclosure",
            "original_boundary_source_bindings.true_amplitude_product_derivatives_retained",
            "phase0_ODE_and_exact_positive_width_pullback.equality_by_functions_and_ODE_uniqueness_not_interval_overlap",
            "common_pressure_primitive_and_true_atom.raw_bridge_axial_V_is_Uz_not_pressure_primitive",
            "scaled_atoms_are_the_original_physical_primitives.actual_inlet_uses_scaled_Phi_Uz_atom_mapping",
            "scaled_atoms_are_the_original_physical_primitives.F0_axial_dependence_cancelled_exactly",
        ),
        "lei_ren_part1_paper_compliant_current_core_interior_moments_check.json": (
            "current_core_six_interior_moment_functions_and_tails_available",
            "current_core_interior_moment_rho4_axial6_and_Euler4_available",
            "current_core_interior_normalized_axis_and_inlet_source_identities_certified",
            "current_core_interior_source_theorem.source_functions_not_chosen_from_enclosures",
            "current_core_interior_source_theorem.passed",
        ),
        "lei_ren_part1_paper_compliant_rooted_core_field_check.json": (
            "independent_source_identities.same_primitive_divergence_identity",
            "independent_source_identities.actual_coefficientwise_radial_and_axial_primitive_identities",
            "independent_source_identities.pressure_squared_profile_primitive_and_epsilon_conversion",
            "exact_source_H_root_and_radial_model_tail_valuation_preserved",
            "original_positive_pressure_primitive_and_datum_preserved",
        ),
        "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json": (
            "source_bindings_and_incidence.actual_production_gauge_differences",
            "independent_production_extraction.unmodified_production_advance_scaled_one_executed",
            "independent_production_extraction.pressure_source_and_pressure_Z_coupling_checked",
            "independent_production_extraction.unavailable_axial_derivative_not_filled_with_zero",
            "independent_Green_and_models.Bessel_Phi_model_and_linear_Psi_model_preserved",
        ),
        "lei_ren_part1_paper_compliant_core_coefficient_rebuild_check.json": (),
        "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source_check.json": (
            "current_original_true_core_inlet_C6_generic_source_functions_installed",
            "current_original_true_core_inlet_axial6_numeric_source_installed",
            "checks.source_rows.preserved_true_core_C5_sources",
            "checks.source_rows.actual_seven_row_complete_moment_histories",
            "checks.physical_units.nonunit_radius_and_Pstar_retained",
            "checks.physical_units.raw_h_and_raw_k_primitive_normalizations_proved",
            "checks.actual_numeric_source.all_density_tails_retained",
        ),
        "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": (
            "current_original_true_prescribed_micro_macro_C6_functions_installed",
            "checks.actual_source.preserved_core_C6_before_actual_C7_comparison_source",
            "checks.actual_source.F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho",
            "checks.actual_source.same_original_both_micro_controls_and_comparison_direction_source",
            "checks.own_six_histories.all_nonzero_incoming_histories_and_sixth_derivatives_retained",
            "checks.own_six_histories.rates_H2_M1_K2_A1_B2_C1_and_B_forcing_not_halved",
        ),
        "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": (
            "current_original_Ra_micro_hydrated_whole_source_functions_installed",
            "independent_symbolic_controls_and_measure.physical_hb_measure_once_verified",
            "actual_live_source.comparison_owner_hydrated_without_ancestor_constructor",
            "actual_live_source.finite_N_boundary_and_global_moment_closure_still_unsupplied",
        ),
    }
    results = {}
    for name, paths in specs.items():
        results[name] = _receipt(receipts, name, gates, paths)
    return results


def _ast_witnesses(texts, gates):
    trees = {name: _tree(text, name, gates) for name, text in texts.items() if name.endswith(".py")}
    ok = {}

    bridge = trees.get("lei_ren_part1_paper_compliant_actual_bridge_integrals.py")
    packet = _function(bridge, "packet")
    prepare = _function(bridge, "prepare")
    feedback = _function(bridge, "own_moment_feedback")
    zero = ast.parse("if chart=='first' and hi==0:\n ell=ell*0;delta_phi=delta_phi*0;delta_V=delta_V*0").body[0]
    ok["actual_bridge.prepare_inlet"] = _assignment(prepare, "data", "self.comparison.inlet(Z,root)")
    ok["actual_bridge.packet_source"] = all((
        _assignment(packet, "p", "self.prepare(Z,root)"),
        _assignment(packet, "angular,axial", "self.contributions(p,chart,q)"),
        _assignment(packet, "ell,delta_phi,delta_V", "self.field_changes(p,angular,axial)"),
        _assignment(packet, "phi0", "p['data']['phi0'].truncate(5)"),
        _assignment(packet, "V0", "p['data']['V0'].truncate(5)"),
        _assignment(packet, "initial", "{n:j.truncate(5) for n,j in p['data']['moments'].items()}"),
        _assignment(packet, "history_delta_V", "p['uniform_delta_V'] if not (chart=='first' and hi==0) else delta_V*0"),
        _assignment(packet, "own", "own_moment_feedback(phi0,V0,initial,theta,delta_phi,history_delta_V)"),
        _assignment(packet, "point_delta_phi", "phi0*exp_minus_one"),
        _assignment(packet, "phi", "phi0+point_delta_phi"),
        _assignment(packet, "V", "V0+delta_V"),
        any(ast.dump(node, include_attributes=False) == ast.dump(zero, include_attributes=False)
            for node in ast.walk(packet)),
    ))
    ok["actual_bridge.feedback_zero_history"] = all((
        _assignment(feedback, "decay", "theta**rate"),
        _assignment(feedback, "baseline[name]", "initial[name]*decay+targets[name]*(1-decay)"),
        _assignment(feedback, "feedback[name]", "changes[name]*((1-decay)/rate)"),
        _assignment(feedback, "actual[name]", "baseline[name]+feedback[name]"),
    ))

    comparison = trees.get("lei_ren_part1_paper_compliant_comparison_point_integrals.py")
    inlet = _function(comparison, "inlet")
    profile = _function(comparison, "source_profile_jet")
    ok["comparison.inlet_atoms"] = all((
        _assignment(inlet, "packet", "self.atoms.field.build_root_rows(24,6) if root else self.atoms.rebuild.rebuild(Z,24,6)"),
        _assignment(inlet, "atoms", "self.atoms.atoms_from_packet(packet,6,shared_root=root)"),
        _assignment(inlet, "local", "[self.source_profile_jet(packet,4,i,root) for i in range(3)]"),
        _keyword(inlet, "dict", "moments", "{name:IntervalTaylor(c,values) for name,values in atoms['actual_core_atom_axial_coefficients'].items()}"),
    ))
    ok["comparison.profile_tail_consumption"] = _contains_text(
        texts["lei_ren_part1_paper_compliant_comparison_point_integrals.py"],
        "packet['rows'][name][n][k]", "self.core.correction*tail_factor(",
        "finite+symmetric(c,tail)")

    atoms = trees.get("lei_ren_part1_paper_compliant_core_integral_atoms.py")
    finite = _function(atoms, "finite_atom_coefficients")
    atoms_from = _function(atoms, "atoms_from_packet")
    ok["atoms.normalizations"] = _assignment(finite, "return", "dict(H=linear(phi,1,8),M=linear(uz,0,4),K=product(phi,uz,1,8),A=product(uz,uz,0,4),B=product(phi,phi,1,16),C=product(phi,phi,0,4))")
    # ``return`` is represented by Return, so check it explicitly as well.
    ok["atoms.normalizations"] = any(
        isinstance(node, ast.Return) and _same_expr(node.value,
            "dict(H=linear(phi,1,8),M=linear(uz,0,4),K=product(phi,uz,1,8),A=product(uz,uz,0,4),B=product(phi,phi,1,16),C=product(phi,phi,0,4))")
        for node in ast.walk(finite))
    ok["atoms.tail_addition"] = all((
        _assignment(atoms_from, "phi", "[row[:order+1] for row in rows['A'][:N+1]]"),
        _assignment(atoms_from, "uz", "[row[:order+1] for row in rows['Uz'][:N+1]]"),
        _assignment(atoms_from, "finite", "finite_atom_coefficients(c,phi,uz,order)"),
        _assignment(atoms_from, "correction", "self.core.correction*factor/math.factorial(k)"),
        _assignment(atoms_from, "atoms", "{name:[value+symmetric(c,errors[name][k]) for k,value in enumerate(values)] for name,values in finite.items()}"),
    ))

    interior = trees.get("lei_ren_part1_paper_compliant_current_core_interior_moments.py")
    source_theorem = _function(interior, "source_theorem")
    ok["interior.normalized_integral_change"] = _contains_text(
        texts["lei_ren_part1_paper_compliant_current_core_interior_moments.py"],
        "s.Integral(norm*t**w", "(t,0,1)", "phi(rho*t,z)", "uz(rho*t,z)",
        "d_rho^i d_Z^k J")
    ok["interior.theorem_function"] = source_theorem is not None and _contains_call(source_theorem, "s.Integral")
    ok["interior.spec"] = _contains_text(
        texts["lei_ren_part1_paper_compliant_current_core_interior_moments.py"],
        "'H':('Phi',None,1,2)", "'M':('Uz',None,0,1)", "'K':('Phi','Uz',1,2)",
        "'A':('Uz','Uz',0,1)", "'B':('Phi','Phi',1,1)", "'C':('Phi','Phi',0,1)")

    common_text = texts["lei_ren_part1_paper_compliant_current_core_common_fixed_point.py"]
    common = trees.get("lei_ren_part1_paper_compliant_current_core_common_fixed_point.py")
    consumer = _function(common, "consumer_source_bindings")
    ok["common.same_fixed_point_consumer"] = _contains_text(
        common_text,
        "def same_source_models_and_domains", "leading_Psi_is_B_rho_over_two=True",
        "leading_Phi_is_same_chi_Bessel_model=True",
        "old_new_envelope_is_only_an_upper_bound_not_a_source_definition=True",
        "advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)")
    ok["common.consumer_function"] = consumer is not None and _contains_text(
        ast.get_source_segment(common_text, consumer) or "",
        "advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)",
        "tail_factor(c,degree=N,radial_order=i,axial_order=k")

    first_text = texts["lei_ren_part1_paper_compliant_current_core_first_interface.py"]
    first = trees.get("lei_ren_part1_paper_compliant_current_core_first_interface.py")
    bindings = _function(first, "defining_boundary_bindings")
    ok["first_interface.operands"] = all((
        _contains_text(first_text, "def defining_boundary_bindings", "original_covering_integral_is_only_an_enclosure=True"),
        _assignment(bindings, "zero", "ast.parse(\"if chart=='first' and hi==0:\\n ell=ell*0;delta_phi=delta_phi*0;delta_V=delta_V*0\").body[0]"),
    ))
    # The helper stores several exact source expressions as arguments to
    # ``binding``.  Authentication plus these literals binds the operands.
    ok["first_interface.binding_literals"] = _contains_text(
        first_text,
        "self.comparison.inlet(Z,root)",
        "p['data']['phi0'].truncate(5)",
        "{n:j.truncate(5) for n,j in p['data']['moments'].items()}",
        "PI_core(4,Z)=4*dress(C)")

    inlet_c6_text = texts["lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"]
    inlet_c6 = trees.get("lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py")
    c6_init = _class_method(inlet_c6, "CurrentOriginalCoreInletAxial6GenericSource", "__init__")
    ok["c6_inlet.core_functions"] = all((
        _assignment(c6_init, "self.U", "a.add(self.core_functions['U0'],a.scale(self.core_functions['Psi'],self.core.parameters['epsilon']))"),
        _assignment(c6_init, "sample", "g.mul(rho,t)"),
        _assignment(c6_init, "phi", "a.at(self.Phi,rho,sample)"),
        _assignment(c6_init, "u", "a.at(self.U,rho,sample)"),
        _assignment(c6_init, "densities", "dict(H=a.scale(phi,g.mul(c(2),t)),M=u,K=a.scale(a.mul(phi,u),g.mul(c(2),t)),A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi))"),
        _assignment(c6_init, "self.inlet_moments", "{key:a.at(jet,rho,c(4)) for key,jet in self.moment_functions.items()}"),
    ))
    ok["c6_inlet.normalized_integrals"] = _contains_text(
        inlet_c6_text, "actual_core_normalized_moment=key",
        "same_unique_nonlinear_core_and_full_density_product=True",
        "self.radius=g.mul(c(4),self.core.parameters['epsilon'])",
        "self.Pstar=g.unary('exp',self.core.parameters['logP'])")

    bridge_c6_text = texts["lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"]
    bridge_c6 = trees.get("lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py")
    bridge_init = _class_method(bridge_c6, "CurrentOriginalPrescribedBridgeC6Functions", "__init__")
    actual_fields = _function(bridge_c6, "actual_fields")
    ok["c6_bridge.initial_projection"] = all((
        _assignment(bridge_init, "self.initial7", "{name:extend(d,jet,7) for name,jet in inlet_owner.inlet_moments.items()}"),
        _assignment(bridge_init, "phi0", "a.at(self.core7['Phi'],self.core.rho,c(4))"),
        _assignment(bridge_init, "U0", "a.add(self.core7['U0'],a.scale(a.at(self.core7['Psi'],self.core.rho,c(4)),self.core.parameters['epsilon']))"),
        _assignment(bridge_init, "actualJ", "histories(b,actual['phi'],actual['U'],{k:b.trim(v) for k,v in self.initial7.items()},s,hb,'actual_micro')"),
    ))
    ok["c6_bridge.distinct_width_and_signs"] = all((
        "self.hb=hb=" in bridge_c6_text and "self.Ra=Ra=" in bridge_c6_text,
        "hbB=hbS" not in bridge_c6_text,
        _assignment(actual_fields, "Vy", "a.scale(a.mul(ratio,full_drive),g.neg(chi))"),
        _assignment(actual_fields, "U", "a.sub(initial['U'],a.integral(a.scale(a.at(Vy,phase,v),g.neg(jacobian)),v,g.zero,phase,actual_prescribed_axial_equation=True,all_hydro_pressure_swirl_sectors_retained=True,source_measure_applied_once=True,prescribed_bridge_chart=label))"),
    ))
    ok["c6_bridge.geometry"] = _contains_text(
        bridge_c6_text,
        "rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))",
        "R=g.mul(Ra,g.unary('exp',g.mul(hb,s)))",
        "F0=b.trim(self.core7['F0'])",
        "Pstar=inlet_owner.Pstar")

    micro_text = texts["lei_ren_part1_paper_compliant_current_original_micro_functions.py"]
    micro = trees.get("lei_ren_part1_paper_compliant_current_original_micro_functions.py")
    micro_init = _class_method(micro, "MicroFunctions", "__init__")
    prefix = _function(micro, "prefix", None)
    moment = _function(micro, "moments", None)
    ok["c2.initial_consumer"] = all((
        _assignment(micro_init, "self.phi0", "flow.jet(jet(actual_core['actual_phi_axial5']))"),
        _assignment(micro_init, "self.V0", "flow.jet(jet(actual_core['actual_raw_V_axial5']))"),
        _assignment(micro_init, "self.initial", "{name:flow.jet(jet(row)) for name,row in actual_core['actual_own_six_moments_axial5'].items()}"),
        _assignment(micro_init, "self.Ra", "flow.factor((0,0,0,1,0))"),
        _contains_call(prefix, "f.multiply"),
        _contains_call(prefix, "f.add"),
    ))
    ok["c2.zero_initial_branch"] = _contains_text(
        micro_text,
        "if lo==hi==0:return self.initial,dict(exact_original_core_atoms=True)",
        "phi=f.multiply(self.phi0,exp)", "V=f.add(self.V0,*deltaV.values())")

    # The source text is authenticated, and each witness records only source
    # operands.  No witness evaluates a model or chooses a range endpoint.
    return ok


def _all_witnesses_true(witnesses, gates):
    for name, value in witnesses.items():
        if value is not True:
            gates.append("ast_witness_missing:" + name)
    return all(witnesses.values())


def _projection_record(receipts, witnesses):
    return {
        "point": {"micro_coordinate": "s=0", "core_radius": "rho=4", "physical_radius": "Ra"},
        "common_core_function": {
            "Phi4": "Phi(rho=4,Z)",
            "U0_plus_epsilon_Psi4": "U0(rho=4,Z)+epsilon*Psi(rho=4,Z)",
            "ownership": "same unique analytic fixed-point core functions",
        },
        "six_normalized_moments": {
            "labels": ["H", "M", "K", "A", "B", "C"],
            "change_of_variables": "rho=4*t, t in [0,1]",
            "recipes": {
                "H": "2*int_0^1 t*Phi(4*t,Z) dt",
                "M": "int_0^1 Uz(4*t,Z) dt",
                "K": "2*int_0^1 t*Phi(4*t,Z)*Uz(4*t,Z) dt",
                "A": "int_0^1 Uz(4*t,Z)^2 dt",
                "B": "int_0^1 t*Phi(4*t,Z)^2 dt",
                "C": "int_0^1 Phi(4*t,Z)^2 dt",
            },
            "ordinary_orders": [0, 1, 2],
            "factorials_and_radial_units_preserved": True,
        },
        "consumer_bindings": {
            "actual_bridge_packet_zero_branch": "chart='first' and hi==0 => source increments ell, delta_phi, delta_V vanish",
            "c2": "actual_phi_axial5, actual_raw_V_axial5, actual_own_six_moments_axial5",
            "c6_inlet": "same Phi/U/moment functions evaluated at rho=4",
            "c6_bridge": "initial7 extends the same six inlet moments; phi0 and U0 are the same core endpoint",
        },
        "projection_kind": "function/projection identity; no interval-value or graph-node equality",
        "tail_errors_retained": True,
        "ordinary_Z_factors_retained": True,
        "physical_Pstar_retained": True,
        "physical_radius_retained": True,
        "U_y_double_negative_retained": True,
    }

def derive(nodes, source_records, *, expected_nodes_digest):
    """Authenticate source records and return the bounded initial projection.

    The function intentionally accepts only serialized records.  It never
    imports any project module and never invokes a producer, checker, helper,
    quadrature, recurrence, norm, or source operator.
    """
    gates = []
    if not isinstance(expected_nodes_digest, str) or expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("expected_nodes_digest_mismatch")
    try:
        actual_nodes_digest = _digest(nodes)
    except Exception as error:
        actual_nodes_digest = ""
        gates.append("nodes_digest_error:" + type(error).__name__)
    if actual_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("nodes_digest_mismatch")
    if not isinstance(nodes, list) or len(nodes) != EXPECTED_NODE_COUNT:
        gates.append("node_count_mismatch")

    records = _source_record_set(source_records)
    texts = {name: _record_text(records, name, gates) for name in RECORD_HASHES}
    report_parent, report_positive = _check_reports(texts, gates)
    receipts = {name: _json(texts, name, gates) for name in RECORD_HASHES if name.endswith("_check.json")}
    receipt_status = _check_receipts(receipts, gates)
    witnesses = _ast_witnesses(texts, gates)
    ast_ok = _all_witnesses_true(witnesses, gates)
    positive_scope_ok = not any(item.startswith((
        "positive_", "native_", "first_micro_")) for item in gates)
    accepted = bool(actual_nodes_digest == EXPECTED_NODES_DIGEST and
                    isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT and
                    ast_ok and positive_scope_ok and all(receipt_status.values()) and not gates)
    # These are deliberate scope limits, not unresolved gates: the first
    # inlet projection is complete without an hbB=hbS premise or evolution.
    limitations = [
        "full_micro_evolution_identity_proved=False",
        "F0_functional_Z7_evolution_not_claimed",
        "bridge_width_hbB_and_switch_width_hbS_not_identified",
        "phase_cutoff_fullfield_density_and_global_closure_not_claimed",
    ]
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": accepted,
        "nodes_digest": actual_nodes_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "common_core_initial_projection_proved": accepted,
        "six_inlet_moment_projection_proved": accepted,
        "all_tail_enclosures_retained": accepted,
        "ordinary_orders": [0, 1, 2],
        "initial_projection": _projection_record(receipts, witnesses),
        "ast_operands": witnesses,
        "accepted_receipts": {name: value for name, value in receipt_status.items() if value},
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_SHA256,
        "accepted_source_stage_report": POSITIVE_REPORT,
        "accepted_source_stage_sha256": POSITIVE_SHA256,
        "full_first_micro_evolution_proved": False,
        "full_micro_evolution_identity_proved": False,
        "F0_identity_proved": False,
        "hb_identity_proved": False,
        "phase_cutoff_identity_proved": False,
        "full_field_identity_proved": False,
        "density_identity_proved": False,
        "global_closure_proved": False,
        "unresolved_gates": [] if accepted else sorted(set(gates)),
        "scope_limitations": limitations,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
    }
