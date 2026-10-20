"""Static source identity for the original two-chart micro evolution.

The consumer authenticates inherited source records and inspects their ASTs.  It
does not import a producer, materialise a graph value, replay an operator, or
turn an interval enclosure into a function value.  Its conclusion is a
function/ODE source correspondence; the density and finite-N seams stay open.
"""
from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction

SCHEMA = "original_micro_C2_C6_evolution_source_identity_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
POSITIVE_REPORT = "lei_ren_part1_paper_compliant_current_original_micro_affine_source_identity.json"
POSITIVE_SHA256 = "3ba7c5829ebc1243ef63d04f092c17613cfda7e76913980e7e16b71c810fe37c"

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    POSITIVE_REPORT: POSITIVE_SHA256,
    "lei_ren_part1_paper_compliant_current_original_micro_functions.py": "d4ef2385ff4627c91f3fd99d1c6819e5598db3e789093f02fe1dc965ec8be779",
    "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": "54d14670660ea41e90ffc25f536e65747da2a3008d8c34e790780cd2d9df8fac",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    "lei_ren_part1_paper_compliant_inner_bridge_profiles.py": "142945ac1dd8a5705c380e79da29ea6bd89374c18cfa50f0412f1ba70c520840",
    "lei_ren_part1_paper_compliant_inner_bridge_profiles_check.json": "ea360228f830af0dc78fdb38f43f8585cc68f6bf7b6e10b6bb3b4adce4af64ff",
    "lei_ren_part1_paper_compliant_actual_bridge_integrals.py": "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3",
    "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": "575fc6b5c1c78957bf9e0d178582bfcff26c167d4de037cd7aea41fc69f52145",
    "lei_ren_part1_paper_compliant_comparison_point_integrals.py": "a7103469587709ff80b087b99fcf5490c8445afcd25c2794a4097b73963d98b7",
    "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": "3b2cd05bd8732acaf3b53b4f97422c55721a4ab9db9c760758b5cbdd0e31171b",
    "lei_ren_part1_paper_compliant_frozen_comparison_field.py": "e38e99800be7ea5d6b142e9101e75e871c97e8dcf7a2b781927f729dcebd9bfb",
    "lei_ren_part1_paper_compliant_frozen_comparison_field_check.json": "68737b02622414f095ce8aecb6daf0dd5c222e0a4479a265cdb672249b00d163",
    "lei_ren_part1_paper_compliant_bridge_mixed_C4.py": "55569c8679beb43eaf4ba9681393a2604d34dfed1a83d66d53c6219cfa7e9187",
    "lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json": "9842df322ae1d84ef79760b6bfe4a01a2d94877639a2be7130571bac19b801e8",
    "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py": "f3ad62f07b7631ccbfe67c787b4989c0d70db2ff49c96325ca3ab2be284ad87c",
    "lei_ren_part1_paper_compliant_anchored_axis_amplitude_check.json": "299b52671d3ff05b91968b1b78e9ab57804abb60815b32d6d469c4d3e47a00a9",
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json": "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    "lei_ren_part1_paper_compliant_current_original_Rm_native_transport_pressure_sources.json": "7c88052348a78ab7e26728f492296fb57858cead1794b9d8cb6e57ad873a7b9b",
    "lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py": "74511ee00a655184021b0516478cffb7c8fdca830b5dc01ab44dfde2e9ad648b",
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json": "572c50c5210ec51c9998780765bdaf3ca9331ead8f4bc8e968245c6a556f8468",
    "lei_ren_part1_paper_compliant_outer_initial.py": "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py": "38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01",
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.json": "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400",
    "lei_ren_part1_paper_shared_analytic_tube.json": "1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6",
    "lei_ren_part1_paper_compliant_K1_ledger.json": "38cd8e76c0d5f5519b510d9cb0502b8915446ffa1c6338d3a3a98c4aca2cdf7d",
    "lei_ren_part1_paper_compliant_core_uniform_bounds.json": "7415cba5b9abc486b2078c450bd0f2e1b7798d1921f5863c00c04526620844e7",
    "lei_ren_part1_paper_compliant_core_uniform_bounds.py": "23d38d188d44e3169d9107bcf2eba201f6b866f4e408c4c83a84aa2f9d1420d7",
}
DECODED_RECORD_HASHES = {}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def _get(value, path, default=None):
    for part in path.split(".") if path else ():
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _bound_endpoint(packet, endpoint):
    """Read a saved directed binary endpoint, never a rounded display value."""
    if not isinstance(packet, Mapping):
        return None
    value = packet.get(endpoint + "_exact_mpf_tuple")
    if not isinstance(value, list) or len(value) != 4 or any(type(v) is not int for v in value):
        return None
    sign, mantissa, exponent, bits = value
    if sign not in (0, 1) or mantissa < 0 or bits != mantissa.bit_length() or abs(exponent) > 10000:
        return None
    result = Fraction(mantissa << exponent) if exponent >= 0 else Fraction(mantissa, 1 << -exponent)
    return -result if sign else result


def _text(records, name, gates):
    rec = records.get(name)
    if not isinstance(rec, Mapping):
        gates.append("missing_source_record:" + name); return ""
    if rec.get("path") != "experiments/root_st073/" + name:
        gates.append("source_path_mismatch:" + name); return ""
    if rec.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name); return ""
    if rec.get("encoding") is not None:
        gates.append("unexpected_encoded_record:" + name); return ""
    value = rec.get("source_text")
    if not isinstance(value, str):
        gates.append("plain_source_text_missing:" + name); return ""
    if hashlib.sha256(value.encode()).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_hash_mismatch:" + name); return ""
    if "\r" in value:
        gates.append("source_crlf:" + name)
    return value


def _tree(value, name, gates):
    try:
        return ast.parse(value, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _json(value, name, gates):
    try:
        return json.loads(value, parse_int=lambda raw: int(raw), parse_float=float)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _fn(tree, name, cls=None):
    if cls:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                        return child
        return None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def _norm(value):
    return "".join(value.split())


def _assign(fn, target, expression):
    if fn is None:
        return False
    target = _norm(target)
    expression_tree = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    for node in ast.walk(fn):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(_norm(ast.unparse(item)) == target for item in targets):
                value = node.value
                if ast.dump(value, include_attributes=False) == expression_tree:
                    return True
    return False


def _has(fn, *needles):
    if fn is None:
        return False
    source = ast.unparse(fn)
    return all(needle.replace(" ", "") in source.replace(" ", "") for needle in needles)


def _source_has(text, *needles):
    compact = text.replace(" ", "").replace("\n", "")
    return all(needle.replace(" ", "").replace("\n", "") in compact for needle in needles)


def _returns_mapping(fn, *keys):
    if fn is None:
        return False
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Dict):
            names = {key.value for key in node.value.keys if isinstance(key, ast.Constant)}
            if set(keys).issubset(names): return True
        if (isinstance(node, ast.Return) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == "dict"):
            names = {key.arg for key in node.value.keywords if key.arg is not None}
            if set(keys).issubset(names):
                return True
    return False


def _return_expression(fn, key, expression):
    """An exact returned operand, independently of source formatting."""
    if fn is None:
        return False
    expected = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Call):
            if isinstance(node.value.func, ast.Name) and node.value.func.id == "dict":
                for item in node.value.keywords:
                    if item.arg == key and ast.dump(item.value, include_attributes=False) == expected:
                        return True
    return False


def _imports(tree, module, symbol):
    return any(isinstance(node, ast.ImportFrom) and node.module == module
        and any(item.name == symbol and item.asname is None for item in node.names)
        for node in ast.walk(tree))


def _receipt(receipts, name, paths, gates):
    value = receipts.get(name, {})
    result = {path: _get(value, path) is True for path in paths}
    if value.get("all_passed") is not True:
        gates.append("receipt_not_accepted:" + name)
    for path, ok in result.items():
        if not ok:
            gates.append("receipt_missing:" + name + ":" + path)
    return bool(value.get("all_passed") is True and result and all(result.values()))


def _reports(texts, gates):
    parent = _json(texts[PARENT_REPORT], PARENT_REPORT, gates)
    positive = _json(texts[POSITIVE_REPORT], POSITIVE_REPORT, gates)
    if _get(positive, "parent_report_sha256") != PARENT_SHA256:
        gates.append("positive_parent_pin_changed")
    family = {
        "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
        "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
        "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
    }
    if _get(positive, "source_family") != family:
        gates.append("positive_source_family_changed")
    micro = _get(positive, "micro_initial_source_identity", {})
    affine = _get(positive, "native_affine_kernel_identity", {})
    for path in ("passed", "common_core_initial_projection_proved", "six_inlet_moment_projection_proved", "all_tail_enclosures_retained"):
        if micro.get(path) is not True: gates.append("accepted_initial_missing:" + path)
    for path in ("passed", "cross_archive_kernel_identity_proved", "global_offset_identity_proved"):
        if affine.get(path) is not True: gates.append("accepted_affine_missing:" + path)
    if micro.get("nodes_digest") != EXPECTED_NODES_DIGEST or affine.get("nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("accepted_stage_graph_changed")
    if (_get(affine, "ra_recipe.proved") is not True
            or _get(affine, "source_parameter_aliases.hbB") != 8
            or _get(affine, "source_parameter_aliases.logP") != 5):
        gates.append("accepted_exact_bridge_parameters_missing")
    rows = affine.get("rows", [])
    micro_rows = [row for row in rows if isinstance(row, Mapping) and row.get("chart") in ("first_micro", "second_micro")]
    if len(micro_rows) != 2 or any(row.get("actual_c2_J_equals_actual_c5_J") is not True
            or _get(row, "c2_handles.Jacobian") != 8 or _get(row, "c5_handles.physical_measure") != 8 for row in micro_rows):
        gates.append("accepted_actual_micro_hb_source_identity_missing")
    scope = positive.get("scope", {})
    for path in ("micro_background_s0_rho4_finite_and_tail_projection", "micro_six_nonzero_inlet_moments"):
        if scope.get(path) is not True: gates.append("accepted_scope_missing:" + path)
    for path in ("micro_finite_N_sc_over_2_correction", "offset_same_density_certificate",
                 "offset_finite_N_history_identity", "offset_physical_Q_y_Z1_seam", "full_micro_evolution_identity_proved",
                 "F0_identity_proved", "hb_identity_proved", "phase_cutoff_identity_proved", "density_identity_proved",
                 "full_goal_complete"):
        if scope.get(path) is not False: gates.append("accepted_scope_overclaims:" + path)
    if _get(scope, "offset_flags.same_density_certificate") is not False:
        gates.append("accepted_scope_overclaims:offset_flags.same_density_certificate")
    if scope.get("ordinary_Z_orders") != [0, 1, 2]: gates.append("ordinary_Z_orders_changed")
    if _get(parent, "source_family") != family: gates.append("parent_source_family_changed")
    pressure_name = "lei_ren_part1_paper_compliant_current_original_Rm_native_transport_pressure_sources.json"
    pressure = _json(texts[pressure_name], pressure_name, gates)
    if (positive.get("accepted_source_stage_sha256") != RECORD_HASHES[pressure_name]
            or _get(pressure, "scope.common_independent_normalized_P0_source_bound") is not True
            or _get(pressure, "first_micro_source_identity.component_proofs.common_pressure_P0_sublemma.passed") is not True):
        gates.append("accepted_common_normalized_P0_source_missing")
    return parent, positive


def _receipts(texts, gates):
    receipts = {name: _json(text, name, gates) for name, text in texts.items() if name.endswith("_check.json")}
    specs = {
        "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json": (
            "full_finite_width_signed_bridge_integral_enclosures_available", "actual_own_six_moment_feedback_enclosures_available"),
        "lei_ren_part1_paper_compliant_comparison_point_integrals_check.json": (
            "original_9_23_namespace_bound", "actual_atoms_and_formal_width_remainders_preserved"),
        "lei_ren_part1_paper_compliant_frozen_comparison_field_check.json": (
            "functional_identities.passed", "actual_core_exit_five_primitive_axial5_enclosures_available"),
        "lei_ren_part1_paper_compliant_inner_bridge_profiles_check.json": (
            "angular_log_sign_bound_to_accepted_same_source_comparison_theorem", "source_bound_smoothed_comparison_axial6_enclosures_available",
            "actual_prescribed_shear_bridge_axial5_enclosures_available"),
        "lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json": (
            "phase1_actual_two_sided_physical_grids_identical", "phase2_actual_histories_retained",
            "actual_core_exit_R100_histories_and_original_P0_retained", "R100_physical_zeroth_radial_source_rows_identical",
            "exact_comparison_integral_and_six_own_moment_histories_bound", "core_bridge_and_R100_local_functional_joins_certified"),
        "lei_ren_part1_paper_compliant_anchored_axis_amplitude_check.json": (
            "six_poles_and_original_unique_anchor_certified", "actual_anchored_G_resolved_with_directed_error",
            "positive_F0_retained_without_exp_logF0", "source_identities.exact_rational_pole_decomposition",
            "source_identities.exact_simple_pole_residue_identity"),
        "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json": (
            "current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed", "source_functions_and_directed_enclosures_are_distinct",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.unique_H_root_in_original_bracket",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.amplitude_uses_same_core",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.amplitude_uses_same_nonlinear_rebuild",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_source_scalar_and_unique_anchor_identity.anchor_leaf_uses_this_exact_source",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.exact_graph_logP_delta_logLambda_epsilon_recipes_verified",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_core_encloses_same_exact_recipes.logP",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_core_encloses_same_exact_recipes.logLambda",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.actual_core_encloses_same_exact_recipes.epsilon",
            "checks.actual_oracle.exact_parameter_and_native_P0_identity.interval_overlap_not_used_as_a_function_identity"),
        "lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json": (
            "analytic_original_beta_and_flat_envelope_checks.all_required_endpoint_flat_limits",
            "analytic_original_beta_and_flat_envelope_checks.sigma_termwise_constants",
            "original_flat_shapes_retained"),
        "lei_ren_part1_paper_compliant_current_core_first_interface_check.json": (
            "phase0_ODE_and_exact_positive_width_pullback.exact_global_positive_scalar_hb_definition_bound",
            "phase0_ODE_and_exact_positive_width_pullback.hb_axial_derivatives_are_zero_by_global_norm_choice",
            "phase0_ODE_and_exact_positive_width_pullback.source_width_cap_not_used_as_hb_value"),
        "lei_ren_part1_paper_compliant_current_original_micro_functions_check.json": (
            "current_original_Ra_micro_hydrated_whole_source_functions_installed", "independent_original_micro_ODE.passed",
            "independent_original_micro_ODE.exact_core_memory_and_three_micro_control_seams_checked", "actual_live_source.passed",
            "actual_live_source.comparison_owner_hydrated_without_ancestor_constructor", "actual_live_source.finite_N_boundary_and_global_moment_closure_still_unsupplied",
            "original_source_bindings.original_micro_controls.first_chi.verified", "original_source_bindings.original_micro_controls.second_and_macro_chi.verified",
            "original_source_bindings.original_micro_controls.comparison_direction.verified", "original_source_bindings.original_micro_controls.actual_F_source_integral.verified",
            "original_source_bindings.original_micro_controls.actual_V_source_integral.verified", "original_source_bindings.recovered_norms_are_upper_bounds_not_source_values"),
        "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check.json": (
            "current_original_true_prescribed_micro_macro_C6_functions_installed", "checks.actual_source.passed",
            "checks.actual_source.preserved_core_C6_before_actual_C7_comparison_source", "checks.actual_source.G_higher_rows_follow_original_anchored_FTC",
            "checks.actual_source.F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho", "checks.actual_source.original_P0_C6_preserved_and_executable_Z7_projection",
            "checks.actual_source.same_original_both_micro_controls_and_comparison_direction_source", "checks.original_direction.physical_swirl_F0_squared_dressing_and_external_scale_once",
            "checks.original_direction.original_linear_quadratic_pressure_and_meridional_terms_retained", "checks.signed_transport.nonunit_jacobian_chi_radius_Pstar_and_F0"),
    }
    status = {}
    for name, paths in specs.items(): status[name] = _receipt(receipts, name, paths, gates)
    mixed = receipts.get("lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json", {})
    for fixture in ("core_log_fixture", "varying_direction_fixture", "controls_fixture"):
        if _get(mixed, fixture + ".finite_fixture_only") is not True:
            gates.append("fixture_scope_changed:" + fixture + ".finite_fixture_only")
        if _get(mixed, fixture + ".actual_source_admission") is not False:
            gates.append("fixture_scope_changed:" + fixture + ".actual_source_admission")
    return receipts, status


def _witnesses(texts, trees, gates):
    micro = "lei_ren_part1_paper_compliant_current_original_micro_functions.py"
    c6 = "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
    inner = "lei_ren_part1_paper_compliant_inner_bridge_profiles.py"
    actual = "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
    comparison = "lei_ren_part1_paper_compliant_comparison_point_integrals.py"
    frozen = "lei_ren_part1_paper_compliant_frozen_comparison_field.py"
    mixed = "lei_ren_part1_paper_compliant_bridge_mixed_C4.py"
    anchored = "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py"
    n1 = "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
    macro = "lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py"
    flat = "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
    outer = "lei_ren_part1_paper_compliant_outer_initial.py"
    gamma = "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py"
    uniform = "lei_ren_part1_paper_compliant_core_uniform_bounds.py"
    mt, ct, it = trees[micro], trees[c6], trees[inner]
    at, pt, ft = trees[actual], trees[comparison], trees[frozen]
    bt, ant, nt, ut = trees[mixed], trees[anchored], trees[n1], trees[uniform]
    minit = _fn(mt, "__init__", "MicroFunctions"); mprefix = _fn(mt, "prefix"); mhist = _fn(mt, "histories"); meval = _fn(mt, "evaluate")
    cinit = _fn(ct, "__init__", "CurrentOriginalPrescribedBridgeC6Functions"); chist = _fn(ct, "histories"); cdirection = _fn(ct, "comparison_directions"); cfields = _fn(ct, "actual_fields")
    idirection = _fn(it, "direction"); ain = _fn(at, "__init__", "CompliantActualBridgeIntegrals"); aprepare = _fn(at, "prepare"); feedback = _fn(at, "own_moment_feedback"); iinputs = _fn(it, "inputs")
    pinlet = _fn(pt, "inlet"); pprofile = _fn(pt, "source_profile_jet"); dress = _fn(ft, "dress")
    formal = _fn(bt, "formal_comparison_source"); ane = _fn(ant, "evaluate"); attach = _fn(nt, "attach_oracle")
    n1init = _fn(nt, "__init__", "CurrentOriginalN1CoreFixedPointFunctions")
    macro_init = _fn(trees[macro], "__init__", "MacroFlow")
    macro_owner = _fn(trees[macro], "owner", "OriginalBridgeMacroFunctions")
    uniform_run = _fn(ut, "run")
    out = {}
    k1name = "lei_ren_part1_paper_compliant_K1_ledger.json"
    ubname = "lei_ren_part1_paper_compliant_core_uniform_bounds.json"
    k1 = _json(texts[k1name], k1name, gates)
    ub = _json(texts[ubname], ubname, gates)
    tube_name = "lei_ren_part1_paper_shared_analytic_tube.json"
    tube = _json(texts[tube_name], tube_name, gates)
    budget_upper = _get(k1, "decreasing_K_smallness_bounds.double_h.upper")
    upper_hb_budget = _bound_endpoint(_get(k1, "decreasing_K_smallness_bounds.double_h"), "upper")
    lower_phi = _bound_endpoint(_get(ub, "actual_Phi_lower"), "lower")
    upper_delta = _bound_endpoint(_get(tube, "delta"), "upper")
    lower_sigma = _bound_endpoint(_get(tube, "sigma"), "lower")
    out["accepted_hb_domain_budget"] = bool(
        _get(k1, "smallness_thresholds.double_h") == ".01"
        and _get(k1, "smallness_checks.double_h") is True
        and _get(k1, "parameter_definition.h_b") == "cstar*K^-100"
        and _get(k1, "parameter_definition.epsilon_b") == "h_b, the identical scalar"
        and _get(k1, "h_b_equals_epsilon_b_by_definition") is True
        and _get(k1, "positive_width_not_materialized") is True
        and budget_upper == "6.03985961041886660757579592153415107609035339111330482224815052001095236774274579994934865e-643"
        and _get(ub, "analytic_core_family_sha256") == "a4056322a3946c0ec0d43cc14d6f82d9a022410a5c077bab3577e87f8e9d7bdf"
        and _get(ub, "real_axial_domain") == ["-1", "1"]
        and _get(ub, "scaled_radial_domain") == ["0", "4.1"]
        and _get(ub, "unique_analytic_field_bounded_not_sampled") is True
        and _get(k1, "base_analytic_core_family_sha256") == _get(ub, "analytic_core_family_sha256")
        and _get(ub, "full_real_axis_normalized_swirl_positive") is True
        and _get(ub, "finite_coefficient_state_used") is False
        and _get(ub, "all_supplied_core_exit_inputs_certified") is True
        and upper_hb_budget is not None and 0 < upper_hb_budget < Fraction(1, 100)
        and lower_phi is not None and lower_phi > 0
        and upper_delta is not None and 0 < upper_delta < Fraction(1, 2)
        and lower_sigma is not None and lower_sigma > 0
        and _get(tube, "common_complex_axis_poles_excluded") is True
        and _get(tube, "analytic_core_family_sha256") == _get(ub, "analytic_core_family_sha256"))
    out["rho_continuation_theorem"] = bool(out["accepted_hb_domain_budget"])
    out["core_uniform_bound_ast"] = all((_assign(uniform_run, "model_lower", "1-qmax/2+qmax*qmax/12-qmax**3/144"),
        _assign(uniform_run, "correction_value", "correction_norm*embedding(c,h,radius,0,0)"),
        _assign(uniform_run, "phi_lower", "model_lower-correction_value"),
        _has(uniform_run, "ifendpoints(phi_lower)[0]<=0", "full_real_axis_normalized_swirl_positive=True", "unique_analytic_field_bounded_not_sampled=True"),
        _source_has(texts[uniform], "Alternating model bound", "Xh embedding/moment averaging", "finite_coefficient_state_used=False")))
    out["c2_initial_actual_core_values"] = all((_assign(minit, "self.phi0", "flow.jet(jet(actual_core['actual_phi_axial5']))"),
        _assign(minit, "self.V0", "flow.jet(jet(actual_core['actual_raw_V_axial5']))"),
        _assign(minit, "self.initial", "{name:flow.jet(jet(row)) for name,row in actual_core['actual_own_six_moments_axial5'].items()}")))
    out["c2_geometry_and_chi"] = all((_source_has(texts[micro], "actual_geometry='R=Ra*exp(hb*s), dy=hb*ds'", "actual_first_equations='phi_y=-((1-sigma)+hb*sigma)*Dbar*phi/2; V_y=-((1-sigma)+hb*sigma)*(phi/barphi)*drive'", "actual_second_equations='phi_y=-hb*Dbar*phi/2; V_y=-hb*(phi/barphi)*drive'"),
        _has(meval, "R=self.Ra*self.series.scalar_series(1,s)[0]", "chi=f.scalar(1-sig)+f.h*sig", "chi=f.h", "true_physical_measure='dy=hb*ds'", "source_interval_functions_not_endpoint_hulls=True")))
    out["c2_enclosure_prefix_not_values"] = all((_has(mprefix, "absolute_full_prefix_log_norm", "full_axial_nonlinear_errors", "source_control_terms"),
        _source_has(texts[micro], "recovered_norms_are_upper_bounds_not_source_values=True", "source_interval_functions_not_endpoint_hulls=True")))
    out["complete_two_chart_control_memory"] = all((
        _has(mprefix, "terms=[('first_micro',f.h*wc),('first_micro',f.h*f.h*ws)]",
            "terms.append(('second_micro',f.h*f.h*(s-1)))", "bridge.pulse_control_masses(c,s)",
            "wc,ws=(c.mpf('.5'),c.mpf('.5'))"),
        _has(mhist, "first_full=self.prefix('first_micro',c.mpf([0,1]))",
            "full={name:[hull(a,b)for a,b in zip(first_full[name],full[name])]for name in ('phi','V')}",
            "result[name]=f.add(f.scale(self.initial[name],decay),f.scale(rhs[name],mass))")))
    out["c2_own_six_history_ode"] = all((_has(mhist, "iflo==hi==0", "H=f.scale(phi,2)", "K=f.scale(f.multiply(phi,V),2)", "A=f.multiply(V,V)", "B=f.multiply(phi,phi)", "C=f.multiply(phi,phi)", "original_core_memory_retained=True", "source_rhs_covers_every_actual_prefix=True", "forname,rateinbridge.RATES.items()"),
        _source_has(texts[micro], "true_incoming_decay=decay", "true_positive_Volterra_mass=mass")))
    out["six_exact_rates_and_densities"] = all((
        _assign(at, "RATES", "dict(H=2,M=1,K=2,A=1,B=2,C=1)"),
        _assign(ct, "RATES", "bridge.RATES"),
        _assign(mhist, "rhs", "dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))"),
        all(_return_expression(_fn(ct, "densities"), key, expression) for key, expression in {
            "H": "a.scale(phi,2)", "M": "U", "K": "a.scale(a.mul(phi,U),2)",
            "A": "a.mul(U,U)", "B": "a.mul(phi,phi)", "C": "a.mul(phi,phi)"}.items())))
    out["c2_f0_ratio_bell_recurrence"] = all((_has(iinputs, "defratios(multiplier)", "out[n-1-j]", "math.factorial(k)", "F0_ratios=ratios(1)", "F0_squared_ratios=ratios(2)"),
        _source_has(texts[inner], "#Bell recurrence at order6; no materialization of F0 or Cstar", "out.append(sum((-multiplier*self.core.Lambda*gradient[j]*out[n-1-j]", "math.factorial(k)")))
    left = _fn(trees[flat], "_sigma_left")
    sigjets = _fn(trees[flat], "sigma_jets")
    point = _fn(trees[outer], "point")
    out["same_exact_flat_sigma"] = all((
        _imports(at, "lei_ren_part1_paper_compliant_flat_pulse_derivatives", "sigma_jets"),
        _assign(left, "odds", "1/(1-x)**2-1/x**2"), _assign(left, "value", "e/(1+e)"),
        _has(sigjets, "IntervalTaylor.constant(c,0,4)", "IntervalTaylor.constant(c,1,4)", "1-left[0]", "(-1)**(n+1)*left[n]"),
        _assign(point, "a", "c.exp(-1/x**2)"), _assign(point, "b", "c.exp(-1/(1-x)**2)"),
        _assign(point, "value", "a/(a+b)"),
        _has(point, "ifv<=0:", "return(c.mpf(0),c.mpf(0))", "ifv>=1:", "return(c.mpf(1),c.mpf(0))"),
        _has(_fn(trees[gamma], "sigma"), "source_binding=parent.active.current.ast_binding(flat_source.stable_sigma)",
            "exact_function_not_its_interval_cap=True", "'exact_original_source_flat_sigma'"),
        _has(meval, "bridge.sigma_jets(c,s)[0]", "sig=c.mpf(0)"),
        _assign(cinit, "sigma", "lambda q:inlet.previous.lower.own.jets.complete.Gamma.sigma(g,q)")))
    out["c6_alpha_rho_chi_radius_measure"] = all((_has(cinit, "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))", "rho=lambdaq:g.mul(c(4),g.unary('exp',g.mul(hb,q)))", "ratio=a.div(atcore(self.radial7['Phi'][1],u),atcore(self.core7['Phi'],u)", "barphi=a.exp(logbar)", "barU=a.add(U0", "chi=g.add(g.one,g.mul(g.sub(hb,g.one),sigma(s)))", "R=g.mul(Ra,g.unary('exp',g.mul(hb,s)))", "physical_measure=hb"),
        _source_has(texts[c6], "domain=[0,2]", "subcharts=dict(first=[0,1],second=[1,2])", "exact_micro_dy_hb_ds_and_macro_dy_length_dm_once=True", "true_prescribed_background_before_finite_N_modification=True")))
    out["c6_ode_history_rates_and_double_minus"] = all((_has(chist, "forname,rateinRATES.items()", "g.mul(jacobian,kernel)", "actual_six_history_rate=rate", "original_incoming_memory_retained=True", "source_measure_applied_once=True"),
        _has(cfields, "full_drive=a.add(a.scale(known['drive_hydro'],radius)", "a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar))", "a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl'])", "Vy=a.scale(a.mul(ratio,full_drive),g.neg(chi))", "a.sub(initial['U'],a.integral(a.scale(a.at(Vy,phase,v),g.neg(jacobian))")))
    out["c6_direction_physical_terms"] = all((_returns_mapping(cdirection, "D_over_R", "drive_hydro", "drive_pressure", "drive_swirl"),
        _has(cdirection, "f,h,k=(b.mul(F0,jet)forjetin(phi,J['H'],J['K']))", "F2=b.mul(F0,F0)", "B,C=(b.mul(F2,J[name])fornamein('B','C'))", "swirl=a.div(dressed_swirl,a.trim(F2)", "same strictly positive full F0 squared"),
        _has(cdirection, "pressure=divide", "hydro=divide", "dressed_swirl=divide")))
    out["c2_c6_direction_bridge"] = all((_returns_mapping(idirection, "D_over_R", "drive_hydro", "drive_pressure", "drive_swirl"),
        _has(idirection, "f=dress(phi,ratios)", "h=dress(moments[MTH],ratios)", "b=dress(moments[MZT]['swirl'],ratios2)", "p=dress(moments[MP],ratios2)", "drive_swirl"),
        _source_has(texts[inner], "including true amplitude derivatives", "drive_pressure=((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L)")))
    out["complete_direction_operands"] = all((
        _assign(idirection, "W", "1-(z*m)*(1-delta)-d*derivative(m)"),
        _assign(idirection, "angular", "h*(1-delta/2)-(z*derivative(h))*((1-delta)/2)-d*derivative(k)+(z*k)*(2*delta-1)"),
        _return_expression(idirection, "D_over_R", "(-W+angular/(2*f))/L"),
        _return_expression(idirection, "drive_hydro", "(-W*v+(m-z*derivative(m))*((1-delta)/2)+(z*a)*(2*delta)-d*derivative(a))/(2*L)"),
        _return_expression(idirection, "drive_pressure", "((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L)"),
        _return_expression(idirection, "drive_swirl", "(-(z*b)*(2*delta)+d*derivative(b)+(z*p)*(2*(1+delta))-d*derivative(p))/(2*L)"),
        _assign(cdirection, "W", "a.sub(one,a.add(a.mul(a.sub(one,de),a.mul(z,M)),a.mul(d,Mz)))"),
        _assign(cdirection, "angular", "a.add(a.mul(a.sub(one,a.scale(de,'1/2')),a.trim(h)),a.neg(a.scale(a.mul(a.sub(one,de),a.mul(z,a.shift(h))),'1/2')),a.neg(a.mul(d,a.shift(k))),a.mul(a.sub(a.scale(de,2),one),a.mul(z,a.trim(k))))"),
        _assign(cdirection, "hydro", "divide(a.add(a.neg(a.mul(W,a.trim(U))),a.scale(a.mul(a.sub(one,de),a.sub(M,a.mul(z,Mz))),'1/2'),a.scale(a.mul(de,a.mul(z,A)),2),a.neg(a.mul(d,Az))))"),
        _assign(cdirection, "pressure", "divide(a.add(a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(P0))),2),a.neg(a.mul(d,a.shift(P0)))))"),
        _assign(cdirection, "dressed_swirl", "divide(a.add(a.neg(a.scale(a.mul(de,a.mul(z,a.trim(B))),2)),a.mul(d,a.shift(B)),a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(C))),2),a.neg(a.mul(d,a.shift(C)))))")))
    out["swirl_physical_germ_statement"] = _source_has(texts[c6], "external physical drive multiplies F0^2 once", "dressed first, then normalized back by the SAME full F0^2") and _source_has(texts[inner], "b=dress(moments[MZT]['swirl'],ratios2)")
    out["c2_physical_amplitude_factor"] = all((
        _assign(macro_init, "self.logs", "tuple(c.mpf(v) for v in (logh,logP2,logF02,logRa,0))"),
        _assign(macro_owner, "logF02", "2*previous.read_interval(c,amplitude['logF0'])"),
        _assign(macro_owner, "logP2", "2*previous.read_interval(c,core['logPstar'])"),
        _assign(macro_owner, "flow", "MacroFlow(c,logh,logP2,logF02,logRa,c.ln(100)-logRa,weight)"),
        _assign(meval, "drives", "{part:f.scale(f.jet(dirs['drive_'+part].truncate(5)),R*R if part=='swirl' else R) for part in self.parts}"),
        _assign(meval, "Vy", "f.add(*[f.scale(f.multiply(quotient,row),-chi*f.factor((0,int(part=='pressure'),int(part=='swirl'),0,0))) for part,row in drives.items()])")))
    out["actual_parameter_consumer_chain"] = all((
        _assign(macro_init, "self.h", "self.factor((1,0,0,0,0))"),
        _assign(macro_owner, "logh", "previous.read_interval(c,original['source_log_hb_enclosure'])"),
        _assign(macro_owner, "logLambda", "previous.read_interval(c,core['logLambda'])"),
        _assign(macro_owner, "logRa", "c.ln(4)-logLambda"),
        _assign(minit, "self.Ra", "flow.factor((0,0,0,1,0))"),
        _assign(cinit, "self.hb", "inlet_owner.previous.lower.parameters['hb']"),
        _assign(cinit, "self.Ra", "inlet_owner.radius"), _assign(cinit, "Pstar", "inlet_owner.Pstar")))
    out["formal_comparison_exact_global_functions"] = all((_has(formal, "recipe=dict", "fields=fields", "moments=moments", "return dict(**recipe"),
        _has(formal, "alpha='1-sigma((t-hb)/hb)'", "core_radius='rho=4exp(t)'", "initial_weight='exp(-lambda*y)'", "kernel='exp(-lambda*(y-t))'", "numeric_covers_are_not_exact_point_values=True"),
        _source_has(texts[mixed], "fields=dict(logphi=dict(initial='log(Phi_core(4,Z))'", "H_y+2H=2phi_bar", "C_y+C=phi_bar^2")))
    point_init = _fn(pt, "__init__", "CompliantComparisonPointIntegrals")
    out["actual_comparison_consumers"] = all((
        _assign(point_init, "self.bridge", "CompliantBridgeMixedC4()"),
        _assign(point_init, "self.source", "self.bridge.formal_comparison_source()"),
        _assign(cinit, "phi0", "a.at(self.core7['Phi'],self.core.rho,c(4))"),
        _assign(cinit, "logbar", "a.add(a.log(phi0,'same positive core inlet Phi'),a.integral(a.scale(ratio,g.mul(alpha,hb,rho(u))),u,g.zero,s,exact_original_comparison_smoothing=True,comparison_not_actual_field=True))"),
        _assign(cinit, "barphi", "a.exp(logbar)"),
        _assign(cinit, "barU", "a.add(U0,a.integral(a.scale(atcore(self.radial7['Psi'][1],u),g.mul(alpha,hb,rho(u),self.core.parameters['epsilon'])),u,g.zero,s,exact_original_comparison_smoothing=True,comparison_not_actual_field=True))"),
        _assign(cinit, "self.comparison_micro", "dict(phase=s,phi=barphi,U=barU,logphi=logbar,moments=barJ,directions=known)")))
    out["actual_integral_source_definitions"] = all((_has(aprepare, "self.comparison.inlet(Z,root)"),
        _has(feedback, "decay=theta**rate", "baseline[name]=initial[name]*decay+targets[name]*(1-decay)", "feedback[name]=changes[name]*((1-decay)/rate)", "actual[name]=baseline[name]+feedback[name]"),
        _source_has(texts[inner], "F='F=f*exp(-.5*integral_0^y chi*Dbar dt)'", "V='V=v-integral_0^y chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dt'")))
    out["enclosure_tail_is_not_value"] = all((_source_has(texts[comparison], "finite+symmetric(c,tail)", "source_width_materialized=False", "cap_used_as_width_value=False"),
        _has(pinlet, "self.atoms.field.build_root_rows(24,6)", "self.atoms.atoms_from_packet(packet,6,shared_root=root)"), _has(pprofile, "symmetric(c,tail)")))
    out["dress_ast_and_pressure_scales"] = all((_has(dress, "c=jet.ctx", "jet*IntervalTaylor(c,[ratios[k]/math.factorial(k)for k in range(jet.order+1)])"),
        _source_has(texts[frozen], "P0='Pstar^2'", "PI='R*F0^2'", "R^2*F0^2*swirl")))
    out["f0_anchor_and_unique_core"] = all((_assign(ane, "logF0", "-self.core.logC-self.core.Lambda*actual_G"),
        _has(attach, "self.oracle_parameter_identity=self.parameter_identity(oracle)", "self._oracle=oracle"),
        _source_has(texts[anchored], "original_anchored_primitive='integral_a^Z L H/(H^2+sigma^2), H(a)=0", "exp_logF0_not_materialized=True", "positive_F0_retained_formally=True"),
        _source_has(texts[n1], "G_anchored_at_unique_H_root=True", "Lambda_epsilon_exactly_one=True", "same_implicit_native_P0=True")))
    out["f0_z_ode_recurrence"] = all((_assign(n1init, "logF", "a.add(a.fixed(g.neg(self.parameters['logC'])),a.neg(a.scale(Gjet,Lambda)))"),
        _assign(n1init, "Gjet", "c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)"),
        _assign(n1init, "F0value", "g.unary('exp',logF.value)"),
        _assign(n1init, "F0", "diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))"),
        _assign(n1init, "gradient", "a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))"),
        _has(ane, "G_gradient_taylor_coefficients=self.gradient_jets(z)"),
        _source_has(texts[micro], "gradient=[fields.previous.read_interval(c,v)forv in amplitude['G_gradient_taylor_coefficients']", "values[n-1-j]", "return[v*math.factorial(n)forn,vinenumerate(values)]")))
    out["positivity_existence_uniqueness"] = all((out["rho_continuation_theorem"], out["core_uniform_bound_ast"], _has(cinit, "same positive unique analytic Phi on rho4..4exp(2hb)", "same positive core inlet Phi"),
        _source_has(texts[anchored], "positive_F0_retained_formally=True"),
        _source_has(texts[n1], "unique_analytic_Banach_solution=True", "admitted_solution_not_a_model_or_enclosure_center=True"),
        out["same_exact_flat_sigma"], out["formal_comparison_exact_global_functions"], out["actual_comparison_consumers"]))
    return out


def _result(gates, witnesses, receipt_status, nodes, actual_digest, reports):
    graph_ok = actual_digest == EXPECTED_NODES_DIGEST and isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT
    initial_ok = all(witnesses.get(key) is True for key in ("c2_initial_actual_core_values", "c2_f0_ratio_bell_recurrence"))
    background_keys = ("c2_geometry_and_chi", "c2_enclosure_prefix_not_values", "c6_alpha_rho_chi_radius_measure",
        "c6_direction_physical_terms", "c2_c6_direction_bridge", "complete_direction_operands", "swirl_physical_germ_statement", "formal_comparison_exact_global_functions", "actual_comparison_consumers",
        "actual_integral_source_definitions", "dress_ast_and_pressure_scales", "c2_physical_amplitude_factor", "actual_parameter_consumer_chain", "complete_two_chart_control_memory", "f0_anchor_and_unique_core", "f0_z_ode_recurrence",
        "positivity_existence_uniqueness", "same_exact_flat_sigma")
    history_keys = ("c2_own_six_history_ode", "six_exact_rates_and_densities", "c6_ode_history_rates_and_double_minus", "actual_integral_source_definitions")
    receipts_ok = all(receipt_status.values()) and bool(receipt_status)
    background = bool(graph_ok and receipts_ok and initial_ok and not gates and all(witnesses.get(key) is True for key in background_keys))
    own = bool(background and all(witnesses.get(key) is True for key in history_keys))
    accepted = bool(graph_ok and receipts_ok and background and own and not gates)
    return {
        "schema": SCHEMA, "source_only": True, "passed": accepted,
        "nodes_digest": actual_digest, "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "full_micro_background_evolution_projection_proved": background,
        "full_micro_own_history_projection_proved": own,
        "full_micro_evolution_identity_proved": accepted,
        "full_micro_own_history_identity_proved": own,
        "full_micro_proved": accepted,
        "ordinary_Z_orders": [0, 1, 2], "ordinary_orders": [0, 1, 2], "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "same_density_certificate": False, "physical_Q_y_Z1": False,
        "numeric_values_computed": False, "old_operators_replayed": False,
        "mathematical_operator_replayed": False, "full_goal_complete": False,
        "ast_witnesses": witnesses,
        "accepted_receipts": {name: True for name, ok in receipt_status.items() if ok},
        "receipt_status": receipt_status,
        "fixture_limits": {fixture: {"finite_fixture_only": True, "actual_source_admission": False}
            for fixture in ("core_log_fixture", "varying_direction_fixture", "controls_fixture")},
        "parent_report": PARENT_REPORT, "parent_report_sha256": PARENT_SHA256,
        "accepted_source_stage_report": POSITIVE_REPORT, "accepted_source_stage_sha256": POSITIVE_SHA256,
        "source_family": _get(reports[1], "source_family", {}),
        "projection_kind": "exact source-function identity with enclosure remainders retained; no finite-prefix equality",
        "physical_swirl_identity": "F0_base^2*C2_drive_swirl_germ == [F0^2*C6_swirl]_germ; bare higher-Z swirl equality is not claimed",
        "ordinary_F0_germ_proof": {
            "definition": "r_m(zeta)=exp(-m*Lambda*(G(Zbase+zeta)-G(Zbase)))",
            "ODE_and_initial": "r_m'=-m*Lambda*G'(Zbase+zeta)*r_m; r_m(0)=1",
            "coefficient_rule": "r[n]=sum_j(-m*Lambda*Gprime[j]*r[n-1-j])/n; ordinary row=n!*r[n]",
            "uniqueness": "The same anchored G/logC/Lambda define one positive F0. The coefficient ODE uniquely determines both m=1 and m=2 normalized germs.",
            "angular": "Common F0/F0base cancels in angular/f; its differentiated terms remain in the numerator.",
            "swirl": "C2 uses raw F0^2 germ divided by constant F0base^2. C6 divides by variable F0^2 then multiplies back before the physical ODE.",
            "time_order_recursion": False},
        "domain_evidence": {"budget_upper_bound": "6.03985961041886660757579592153415107609035339111330482224815052001095236774274579994934865e-643",
            "accepted_threshold": "double_h <= .01", "parameter_definition": "h_b=cstar*K^-100=epsilon_b",
            "continuation": "0<=h_b*s<=2*h_b<=1/100; exp(x)<=1/(1-x)<=100/99; rho=4*exp(h_b*s)<=400/99<41/10 for 0<=s<=2",
            "source_domain": "real axial Z in [-1,1], scaled radial rho in [0,4.1]", "actual_hb_materialized": False,
            "comparison_arithmetic": "saved exact binary upper/lower endpoints; rounded decimal displays are diagnostic only",
            "positive_denominators": "Phi>=saved positive lower; L=1-delta*Z^2>=1/2; H^2+sigma^2>0; barPhi=exp(logbar)>0",
            "hb_positivity": "accepted exact global hb=cstar*K^-100 positive scalar theorem; the double_h budget is an upper bound only",
            "uniqueness": "continuous analytic coefficients on the full compact phase interval give a triangular linear field ODE and six nonreset linear history ODEs"},
        "sublemmas": {"initial": initial_ok, "background": background, "own_history": own,
            "F0_same_anchored_function": witnesses.get("f0_anchor_and_unique_core", False) and witnesses.get("f0_z_ode_recurrence", False),
            "same_exact_ODE_and_uniqueness": witnesses.get("positivity_existence_uniqueness", False),
            "same_six_history_rates_and_initials": own,
            "all_old_prefixes_are_enclosures": witnesses.get("c2_enclosure_prefix_not_values", False) and witnesses.get("enclosure_tail_is_not_value", False)},
        "scope_flags": {"same_density_certificate": False, "physical_Q_y_Z1": False, "finite_N_sc_over_2": False,
            "time_order_recursion": False, "full_goal_complete": False},
        "unresolved_gates": [] if accepted else sorted(set(gates)),
        "potential_gates": ["same_density_certificate", "physical_Q_y_Z1", "finite_N_sc_over_2_correction", "global_full_goal_closure"],
        "record_hashes": dict(RECORD_HASHES), "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "source_inspection_audit": {"source_records_only": True, "ast_and_saved_bound_comparisons_only": True, "nodes_hashed_only": True,
            "new_source_consumer_calls": 1, "accepted_ancestor_operator_calls": 0,
            "producer_provider_builder_checker_replayed": False, "quadrature_or_norm_calls": 0},
        "worker_role": "luna_max_worker", "runtime_model_metadata": None,
    }


def derive(nodes, source_records, *, expected_nodes_digest):
    """Return the fail-closed source-function projection for both micro charts."""
    gates = []
    if expected_nodes_digest != EXPECTED_NODES_DIGEST: gates.append("expected_nodes_digest_mismatch")
    try: actual_digest = _digest(nodes)
    except Exception as error:
        actual_digest = ""; gates.append("nodes_digest_error:" + type(error).__name__)
    if actual_digest != EXPECTED_NODES_DIGEST: gates.append("nodes_digest_mismatch")
    if not isinstance(nodes, list) or len(nodes) != EXPECTED_NODE_COUNT: gates.append("node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    texts = {name: _text(records, name, gates) for name in RECORD_HASHES}
    reports = _reports(texts, gates)
    _, receipt_status = _receipts(texts, gates)
    trees = {name: _tree(value, name, gates) for name, value in texts.items() if name.endswith(".py")}
    witnesses = _witnesses(texts, trees, gates)
    for name, ok in witnesses.items():
        if ok is not True: gates.append("ast_witness_missing:" + name)
    result = _result(gates, witnesses, receipt_status, nodes, actual_digest, reports)
    return result
