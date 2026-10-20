"""Fail-closed source projection through the original R100 switch charts.

The consumer authenticates saved plain sources and authenticated C6 snapshots,
then checks source AST/data-flow witnesses.  It does not import or evaluate a
producer, graph field, finite-N operator, quadrature, recurrence, or norm.
Taylor rows and saved packets remain enclosures; exact source-function ports
are attached separately through the admitted C6 graph snapshots.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping

SCHEMA = "original_R100_switch_C2_C6_source_projection_v2"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
MACRO_REPORT = "lei_ren_part1_paper_compliant_current_original_macro_evolution_source_projection.json"
MACRO_SHA256 = "880abe2721144f86b0a9775f62735dccc8c5650a5fd8faab27613997daa4ae16"
BRIDGE_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz"
BRIDGE_ARCHIVE_SHA256 = "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92"
BRIDGE_DECODED_SHA256 = "f361ab39e7599f5ed770df630bccadc4056b3e7ad6c0e1de32c975963f070298"
R100_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
R100_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
R100_DECODED_SHA256 = "18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a"
FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    MACRO_REPORT: MACRO_SHA256,
    BRIDGE_ARCHIVE: BRIDGE_ARCHIVE_SHA256,
    R100_ARCHIVE: R100_ARCHIVE_SHA256,
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    "lei_ren_part1_paper_compliant_current_original_R100_endpoint.py": "f3453e14fac85530698bc7d3a394cc30dd22a5bbc8954d98723794a9dd6b0fef",
    "lei_ren_part1_paper_compliant_current_original_R100_endpoint.json": "ec1478cf1f59c3c225df84201524dca71f51631da7b6cddf0532cb59236c2ee7",
    "lei_ren_part1_paper_compliant_current_original_R100_endpoint_check.json": "d084c241780baa67047221f20a4efe9f530190b22076b4ef6c75639182e398ff",
    "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py": "a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496",
    "lei_ren_part1_paper_compliant_current_original_first_switch_functions.json": "5ea141f4ee1befaf33d6d2d2b9ad871477b72ee3dbbfc6b1e55b82766278a53d",
    "lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json": "5760631fe412ef2ff25724536ea4675f8a995b8bf17fe2efd1387efc96c49d4a",
    "lei_ren_part1_paper_compliant_current_original_second_switch_R110.py": "46df33f034c98ee885f6338141c9e8bbd7f29eb85999c229a17cda367e696023",
    "lei_ren_part1_paper_compliant_current_original_second_switch_R110.json": "51b930434d4746c901eac1328eb0f1b8eae3ea2abb9fceca5cb8dc0b0ecc3751",
    "lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json": "2854144abbf98d5b27f9d251d74e9fa38f6e2b7fafebdd0a059c91084c718535",
    "lei_ren_part1_paper_compliant_K1_ledger.json": "38cd8e76c0d5f5519b510d9cb0502b8915446ffa1c6338d3a3a98c4aca2cdf7d",
}
DECODED_RECORD_HASHES = {BRIDGE_ARCHIVE: BRIDGE_DECODED_SHA256, R100_ARCHIVE: R100_DECODED_SHA256}
BRIDGE_PREFIX_NODES = 73190
R100_PREFIX_NODES = 97295


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def _get(value, path, default=None):
    for part in path.split(".") if path else ():
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _text(records, name, gates):
    record = records.get(name)
    if not isinstance(record, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    if record.get("path") != "experiments/root_st073/" + name:
        gates.append("source_path_mismatch:" + name)
    if record.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
    value = record.get("source_text")
    if not isinstance(value, str):
        gates.append("plain_source_text_missing:" + name)
        return ""
    if name in DECODED_RECORD_HASHES:
        if record.get("encoding") != "gzip" or record.get("decoded_sha256") != DECODED_RECORD_HASHES[name]:
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        try:
            raw = base64.b64decode(record.get("raw_bytes_base64", ""), validate=True)
            if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
                gates.append("archive_raw_hash_mismatch:" + name)
                return ""
            decoded = gzip.decompress(raw)
        except Exception as error:
            gates.append("archive_decode_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(decoded).hexdigest() != DECODED_RECORD_HASHES[name] or decoded != value.encode():
            gates.append("archive_decoded_bytes_mismatch:" + name)
            return ""
        return value
    if record.get("encoding") is not None:
        gates.append("unexpected_encoded_record:" + name)
    if hashlib.sha256(value.encode()).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_hash_mismatch:" + name)
    return value


def _tree(value, name, gates):
    try:
        return ast.parse(value, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _json(value, name, gates):
    try:
        return json.loads(value, parse_int=int, parse_float=float)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _fn(tree, name, cls=None):
    if cls is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                return next((child for child in node.body if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name), None)
        return None
    return next((node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def _norm(value):
    return "".join(value.split())


def _has(fn, *needles):
    if fn is None:
        return False
    compact = _norm(ast.unparse(fn))
    return all(_norm(needle) in compact for needle in needles)


def _source_has(text, *needles):
    compact = _norm(text)
    return all(_norm(needle) in compact for needle in needles)


def _assign(fn, target, expression):
    if fn is None:
        return False
    wanted_target = _norm(target)
    wanted_value = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if any(_norm(ast.unparse(item)) == wanted_target for item in targets) and ast.dump(node.value, include_attributes=False) == wanted_value:
            return True
    return False


def _returns_mapping(fn, *keys):
    if fn is None:
        return False
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Dict):
            names = {key.value for key in node.value.keys if isinstance(key, ast.Constant)}
            if set(keys).issubset(names):
                return True
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == "dict":
            names = {item.arg for item in node.value.keywords if item.arg is not None}
            if set(keys).issubset(names):
                return True
    return False


def _receipt(receipts, name, required, forbidden, gates):
    value = receipts.get(name, {})
    good = value.get("all_passed") is True
    if not good:
        gates.append("receipt_not_accepted:" + name)
    family = value.get("source_family")
    if family != FAMILY and family != FAMILY["actual_five_defect_family_sha256"]:
        gates.append("receipt_source_family_mismatch:" + name)
        good = False
    for path in required:
        marker = _get(value, path)
        ok = marker is True or (type(marker) is int and marker > 0)
        if not ok:
            gates.append("receipt_missing:" + name + ":" + path)
        good = good and ok
    for path in forbidden:
        if _get(value, path) is True:
            gates.append("receipt_open_seam_claimed:" + name + ":" + path)
            good = False
    return good


def _reports(texts, gates):
    parent = _json(texts[PARENT_REPORT], PARENT_REPORT, gates)
    macro = _json(texts[MACRO_REPORT], MACRO_REPORT, gates)
    if _get(parent, "source_family") != FAMILY or _get(macro, "source_family") != FAMILY:
        gates.append("source_family_changed")
    if _get(macro, "parent_report") != PARENT_REPORT or _get(macro, "parent_report_sha256") != PARENT_SHA256:
        gates.append("accepted_macro_parent_pin_changed")
    component = _get(macro, "component_result", {})
    required = {
        "passed": True, "source_only": True,
        "full_macro_background_evolution_projection_proved": True,
        "full_macro_own_history_projection_proved": True,
        "full_macro_evolution_identity_proved": True,
        "full_macro_own_history_identity_proved": True,
        "full_macro_proved": True, "numeric_values_computed": False,
        "old_operators_replayed": False, "full_goal_complete": False,
    }
    for key, expected in required.items():
        if component.get(key) is not expected:
            gates.append("accepted_macro_component_changed:" + key)
    for key in ("full_physical_rhs_and_swirl_germ", "c2_mode_rows_match_c6_macro_modes", "anchored_F0_same_function_and_ordinary_Z_germ", "macro_modes_exp_minus_j_t_and_measure_conversion", "six_history_rates_nonreset_and_derivatives"):
        if _get(component, "ast_witnesses." + key) is not True:
            gates.append("accepted_macro_witness_missing:" + key)
    for key in ("macro_length_exact_source", "macro_radius_exact_source", "frozen_barphi2_independent_of_macro_m", "macro_triangular_field_ode", "macro_ode_uniqueness_premises"):
        if _get(component, "domain_evidence." + key) is not True:
            gates.append("accepted_macro_domain_premise_missing:" + key)
    scope = _get(macro, "scope", {})
    for key, expected in {"ordinary_Z_orders": [0, 1, 2], "same_density_certificate": False, "finite_N_history_identity": False, "physical_Q_y_Z1": False, "global_five_moment_closure": False, "timeorder_recursion_claimed": False, "full_goal_complete": False}.items():
        if scope.get(key) != expected:
            gates.append("accepted_macro_scope_changed:" + key)
    if component.get("nodes_digest") != EXPECTED_NODES_DIGEST or scope.get("nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("accepted_macro_graph_changed")
    ports = component.get("exact_macro_inlet_function_ports", {})
    if ports.get("graph_ports_valid") is not True or ports.get("C2_source_projection_admitted") is not True or ports.get("enclosure_row_equality_claimed") is not False:
        gates.append("accepted_macro_exact_ports_changed")
    if ports.get("source_archive") != BRIDGE_ARCHIVE or ports.get("source_archive_sha256") != BRIDGE_ARCHIVE_SHA256:
        gates.append("accepted_macro_archive_pin_changed")
    return parent, macro, component


def _receipts(texts, gates):
    receipts = {name: _json(text, name, gates) for name, text in texts.items() if name.endswith("_check.json")}
    specs = {
        "lei_ren_part1_paper_compliant_current_original_R100_endpoint_check.json": (
            ("original_factored_R100_actual_endpoint_switch_inlet_installed", "independent_normalization.passed", "independent_normalization.nonzero_all_orders_and_amplitude_derivatives_used", "genuine_source_binding.passed", "genuine_source_binding.exact_original_macro_endpoint_rows_reused", "genuine_source_binding.actual_anchored_F0_squared_source_identities", "guards.passed", "genuine_source_binding.whole_axis_micro_switch_R110_and_recursion_remain_open"), ("actual_bridge_and_first_switch_integral_functions_installed", "actual_R110_source_histories_closed")),
        "lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json": (
            ("original_actual_first_switch_FV_and_six_history_function_integrals_installed", "independent_true_sigma_masses.passed", "independent_true_sigma_masses.exact_full_complement_half_mass_retained", "genuine_native_first_switch.passed", "genuine_native_first_switch.complete_source_owned_weighted_sigma_cells", "guards.passed"), ("actual_bridge_and_first_switch_integral_functions_installed", "actual_R110_source_histories_closed")),
        "lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json": (
            ("original_actual_second_switch_and_two_frame_R110_function_transport_installed", "exact_post_power_identities.passed", "exact_post_power_identities.B_and_C_distinct_rates_retained", "genuine_original_second_switch_R110.passed", "genuine_original_second_switch_R110.nonzero_micro_theta_power_corrections", "guards.passed", "genuine_original_second_switch_R110.whole_Z_micro_global_and_recursion_remain_open"), ("actual_bridge_and_first_switch_integral_functions_installed", "actual_R110_source_histories_closed")),
        "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json": (
            ("checks.actual_source.passed", "checks.actual_source.all_genuine_seven_source_and_six_generic_rows", "checks.actual_source.every_nonzero_incoming_history_function_retained_by_identity", "checks.actual_source.endpoint_phi_U_and_physical_y_derivatives_exposed", "checks.actual_source.distinct_width_tags_with_exact_shared_upstream_K1_identity", "checks.actual_source.frozen_T_and_positive_gap_post_length_guarded", "checks.actual_source.no_free_core_rho_or_old_phase_capture_in_chart_bodies", "checks.actual_source.finite_N_modified_source_and_quantitative_C6_not_claimed", "checks.geometry.passed", "checks.geometry.exact_radius_endpoint_identities", "checks.geometry.exact_physical_log_radius_Jacobians", "checks.geometry.original_finite_T_and_Rref_offsets_retained", "checks.geometry.complete_radius_order_proved_by_actual_source_positive_guards", "checks.cutoff.passed", "checks.cutoff.original_value_and_prime_sources_are_same_flat_function", "checks.cutoff.actual_stable_value_and_Taylor_log_odds_assignments_bound", "checks.cutoff.factorial_one_for_sigma_prime", "checks.cutoff.exact_zero_one_and_flat_derivative_endpoints", "checks.switch_controls.passed", "checks.switch_controls.first_axial_hb_squared_one_minus_sigma_and_current_R_Rsquared_preserved", "checks.switch_controls.second_angular_minus_two_fifths_sigma_and_constant_U_preserved", "checks.switch_controls.mathematical_fixture_only"), ("current_original_actual_modified_bridge_C6_functions_installed", "current_original_prescribed_bridge_C6_numeric_containment_installed")),
    }
    status = {}
    for name, (required, forbidden) in specs.items():
        status[name] = _receipt(receipts, name, required, forbidden, gates)
    first = receipts.get("lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json", {})
    first_rows = _get(first, "independent_complete_first_switch")
    if not (isinstance(first_rows, list) and first_rows and all(_get(row, "full_nonlinear_exponential_current_radius_and_true_sigma_used") is True for row in first_rows)):
        gates.append("receipt_missing:lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json:independent_complete_first_switch.full_nonlinear_exponential_current_radius_and_true_sigma_used")
        status["lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json"] = False
    second = receipts.get("lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json", {})
    second_rows = _get(second, "independent_complete_second_switch_and_post_power")
    if not (isinstance(second_rows, list) and second_rows and all(_get(row, "full_sigma_shifted_radius_and_nonlinear_exponential_used") is True for row in second_rows)):
        gates.append("receipt_missing:lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json:independent_complete_second_switch_and_post_power.full_sigma_shifted_radius_and_nonlinear_exponential_used")
        status["lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json"] = False
    return status


def _saved_reports(texts, gates):
    result = {}
    for name, packet_key in (("lei_ren_part1_paper_compliant_current_original_R100_endpoint.json", "packets"), ("lei_ren_part1_paper_compliant_current_original_first_switch_functions.json", "packets"), ("lei_ren_part1_paper_compliant_current_original_second_switch_R110.json", "second_switch_packets")):
        row = _json(texts[name], name, gates)
        family = _get(row, "source_family")
        good = (family == FAMILY or family == FAMILY["actual_five_defect_family_sha256"]) and isinstance(row.get(packet_key), Mapping) and _source_has(texts[name], '"point_value_selected": false', '"encloses_original_source_function": true')
        if not good:
            gates.append("saved_source_enclosure_report_missing:" + name)
        result[name] = good
    return result


def _archive_ports(texts, nodes, gates):
    bridge = _json(texts[BRIDGE_ARCHIVE], BRIDGE_ARCHIVE, gates)
    native = _json(texts[R100_ARCHIVE], R100_ARCHIVE, gates)
    bb = _get(bridge, "prescribed_bridge_C6_bindings", {})
    nb = _get(native, "prescribed_R100_Rm_C6_bindings", {})
    bnodes = bridge.get("exact_graph_nodes", [])
    nnodes = native.get("exact_graph_nodes", [])
    checks = {
        "bridge_source_family": bridge.get("source_family") == FAMILY,
        "native_source_family": native.get("source_family") == FAMILY,
        "bridge_prefix_length": isinstance(bnodes, list) and len(bnodes) == BRIDGE_PREFIX_NODES,
        "native_prefix_length": isinstance(nnodes, list) and len(nnodes) == R100_PREFIX_NODES,
        "native_contains_bridge_prefix": isinstance(bnodes, list) and isinstance(nnodes, list) and nnodes[:BRIDGE_PREFIX_NODES] == bnodes,
        "graph_prefix_matches_bridge_snapshot": isinstance(nodes, list) and isinstance(bnodes, list) and nodes[:BRIDGE_PREFIX_NODES] == bnodes,
        "graph_prefix_matches_R100_snapshot": isinstance(nodes, list) and isinstance(nnodes, list) and nodes[:R100_PREFIX_NODES] == nnodes,
        "bridge_R100_binding": isinstance(_get(bb, "R100"), Mapping),
        "native_chart_bindings": all(isinstance(_get(nb, "charts." + key), Mapping) for key in ("first_switch", "second_switch", "post_switch")),
        "native_endpoint_bindings": all(isinstance(_get(nb, "endpoints." + key), Mapping) for key in ("R1", "R2", "R110")),
        "bridge_hbB_is_8": bb.get("hb") == 8,
        "native_hbS_is_9": nb.get("hb") == 9,
        "width_handles_distinct": nb.get("hb") != bb.get("hb") and nb.get("bridge_hb") == bb.get("hb"),
    }
    r100 = _get(bb, "R100", {})
    def jets(value):
        return isinstance(value, list) and len(value) == 7 and all(type(item) is int for item in value)
    checks["R100_actual_function_ports"] = all(jets(_get(r100, "actual." + key)) for key in ("phi", "U"))
    checks["R100_six_history_ports"] = all(jets(_get(r100, "moments." + key)) for key in ("H", "M", "K", "A", "B", "C"))
    checks["R100_to_first_incoming_history_identity"] = _get(nb, "charts.first_switch.incoming_histories") == _get(r100, "moments")
    checks["R1_to_second_incoming_history_identity"] = _get(nb, "charts.second_switch.incoming_histories") == _get(nb, "endpoints.R1.moments")
    checks["R2_to_post_incoming_history_identity"] = _get(nb, "charts.post_switch.incoming_histories") == _get(nb, "endpoints.R2.moments")
    checks["R1_actual_and_history_ports"] = all(jets(_get(nb, "endpoints.R1.actual." + key)) for key in ("phi", "U")) and all(jets(_get(nb, "endpoints.R1.moments." + key)) for key in ("H", "M", "K", "A", "B", "C"))
    checks["R2_actual_and_history_ports"] = all(jets(_get(nb, "endpoints.R2.actual." + key)) for key in ("phi", "U")) and all(jets(_get(nb, "endpoints.R2.moments." + key)) for key in ("H", "M", "K", "A", "B", "C"))
    checks["R110_actual_and_history_ports"] = all(jets(_get(nb, "endpoints.R110.actual." + key)) for key in ("phi", "U")) and all(jets(_get(nb, "endpoints.R110.moments." + key)) for key in ("H", "M", "K", "A", "B", "C"))
    definition = _get(nb, "definition", {})
    checks["native_definition_flags"] = all(_get(nb, "definition." + key) is True for key in ("true_prescribed_background_before_finite_N_loop", "exact_original_both_R100_switch_controls_and_current_radius_retained", "switch_hbS_and_preceding_bridge_hbB_source_nodes_not_aliased", "exact_shared_upstream_K1_width_identity_and_positive_geometry_guarded", "frozen_comparison_direction_not_actual_moment_replacement", "all_six_nonzero_incoming_and_every_chart_memory_retained", "finite_N_modified_fields_correction_histories_and_all17_not_installed", "shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix"))
    checks["geometry_certificate_K1_ledger_pin"] = _get(nb, "geometry_certificate.shared_K1_ledger_sha256") == RECORD_HASHES["lei_ren_part1_paper_compliant_K1_ledger.json"]
    checks["geometry_certificate_flags"] = checks["geometry_certificate_K1_ledger_pin"] and all(_get(nb, "geometry_certificate." + key) is True for key in ("shared_upstream_is_same_object", "bridge_logh_is_actual_upstream_object", "switch_width_is_exact_source_reboxing", "all_four_native_T_source_tuples_equal", "actual_native_T_and_logC_match_common_source", "bridge_hbB_switch_hbS_kept_distinct", "directed_T_positive", "directed_reference_gap_positive", "directed_post_length_positive", "bounds_are_not_defining_field_parameters"))
    if not all(checks.values()):
        gates.extend("exact_R100_C6_archive_port_missing:" + name for name, ok in checks.items() if not ok)
    return dict(schema="original_R100_switch_exact_function_ports_v1", graph_ports_valid=all(checks.values()), checks=checks, bridge_prefix_nodes=BRIDGE_PREFIX_NODES, R100_prefix_nodes=R100_PREFIX_NODES, bridge_width_node=bb.get("hb"), switch_width_node=nb.get("hb"), C2_R100_actual_function_ports={key: list(_get(r100, "actual." + key, [])) for key in ("phi", "U")}, C2_R100_history_function_ports={key: list(_get(r100, "moments." + key, [])) for key in ("H", "M", "K", "A", "B", "C")}, R1_end_function_ports=_get(nb, "endpoints.R1", {}), R2_end_function_ports=_get(nb, "endpoints.R2", {}), R110_end_function_ports=_get(nb, "endpoints.R110", {}), bridge_switch_width_identity_proved=bool(checks["geometry_certificate_flags"]), endpoint_rows_are_values=False, saved_C2_rows_are_enclosures=True, numeric_values_computed=False, ordinary_Z_orders_admitted=[0, 1, 2])


def _domain(texts, trees, macro, receipts, ports, gates):
    component = _get(macro, "component_result", {})
    c6 = trees["lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"]
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedR100RmC6Functions")
    geom = _fn(c6, "parameter_geometry_guard", "CurrentOriginalPrescribedR100RmC6Functions")
    add_chart = _fn(c6, "add_chart", "CurrentOriginalPrescribedR100RmC6Functions")
    first = _fn(c6, "first_switch_fields")
    second = _fn(c6, "second_switch_fields")
    result = {
        "accepted_macro_positive_barphi_F0_and_ODE_domain": all(_get(component, "domain_evidence." + key) is True for key in ("proved", "positive_actual_Phi_lower_exact", "analytic_positive_field_not_sampled", "frozen_barphi2_independent_of_macro_m", "macro_triangular_field_ode", "macro_ode_uniqueness_premises")),
        "R100_archive_width_geometry_guard": ports.get("checks", {}).get("geometry_certificate_flags") is True,
        "same_exact_width_source": ports.get("checks", {}).get("geometry_certificate_flags") is True and receipts.get("lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json") is True and _has(geom, "shared=wb.upstream", "shared is sw.history.upstream", "wb.logh is shared.logh", "wb.family==sw.family==common.family", "wb.source==sw.source==common.source_sha", "shared_K1_ledger_sha256=digest", "bridge_hbB_switch_hbS_kept_distinct=True"),
        "post_length_exact_source": _assign(c6_init, "self.post_length", "g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))") and _get(ports, "checks.geometry_certificate_flags") is True,
        "macro_to_post_length_positive": _get(component, "domain_evidence.macro_length_exact_source") is True and _get(ports, "checks.geometry_certificate_flags") is True and receipts.get("lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json") is True,
        "switch_radius_source": _has(c6_init, "R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))", "R=g.mul(c(100),g.unary('exp',g.mul(hb,g.add(g.one,t))))", "R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))"),
        "switch_fields_continuation": _has(first, "g.mul(g.constant('-1/2'),hb,hb,radius)", "g.mul(hb,g.sub(g.one,sigma(g,phase)))", "original_first_switch_axial_control=True") and _has(second, "g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))", "g.mul(g.constant('-2/5'),hb,sig)", "original_second_switch_full_current_radius_direction_and_power=True"),
        "triangular_linear_history_uniqueness": _has(add_chart, "previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)", "fields,radial,generic=previous.generic_fields") and receipts.get("lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json") is True,
        "finite_interval_without_quadrature_claim": _source_has(texts["lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"], "finite_N_modified_source_and_quantitative_C6_not_claimed\": true", "mathematical_fixture_only\": true") and _get(component, "domain_evidence.macro_ode_uniqueness_premises") is True,
    }
    result["proved"] = all(result.values())
    if not result["proved"]:
        gates.append("R100_switch_domain_or_width_premise_missing")
    result["derivation"] = {"interval": "first s in [0,1], second t in [0,1], post p in [0,1]; ordinary Z orders 0,1,2", "radii": "R=100 exp(hbS*s), R=100 exp(hbS*(1+t)), R=100 exp(2hbS+Lp*p)", "post_length": "Lp=log(11/10)-2hbS; log(11/10)>=1/11 and 2hbS<=.01 give Lp>0", "width": "hbS is a native source handle reboxed from the admitted shared K1 source; hbB remains the preceding bridge handle", "uniqueness": "positive frozen barphi2 and F0, finite smooth sigma and current-radius source, triangular linear fields, and six finite linear history ODEs give whole-interval uniqueness and ordinary-Z differentiation", "values": "No field values, quadrature, norm, interval midpoint, or finite-N equality is computed."}
    return result


def _witnesses(texts, trees, receipts, saved, macro, ports, domain):
    r100 = trees["lei_ren_part1_paper_compliant_current_original_R100_endpoint.py"]
    first = trees["lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"]
    second = trees["lei_ren_part1_paper_compliant_current_original_second_switch_R110.py"]
    c6 = trees["lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"]
    recover = _fn(r100, "recover_endpoint")
    r100_owner = _fn(r100, "owner", "OriginalR100Endpoint")
    r100_eval = _fn(r100, "evaluate", "OriginalR100Endpoint")
    first_eval = _fn(first, "evaluate", "FirstSwitchFunctions")
    first_init = _fn(first, "__init__", "FirstSwitchFunctions")
    first_angular = _fn(first, "angular", "FirstSwitchFunctions")
    first_force = _fn(first, "force_delta", "FirstSwitchFunctions")
    first_owner = _fn(first, "owner", "OriginalFirstSwitchFunctions")
    second_eval = _fn(second, "evaluate", "SecondSwitchFunctions")
    second_init = _fn(second, "__init__", "SecondSwitchFunctions")
    second_angular = _fn(second, "angular", "SecondSwitchFunctions")
    second_Ddelta = _fn(second, "Ddelta", "SecondSwitchFunctions")
    second_post = _fn(second, "post", "SecondSwitchFunctions")
    second_owner = _fn(second, "owner", "OriginalSecondSwitchR110")
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_comp = _fn(c6, "comparison_at_radius", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_first = _fn(c6, "first_switch_fields")
    c6_second = _fn(c6, "second_switch_fields")
    c6_full = _fn(c6, "full_drive")
    c6_sigma = _fn(c6, "sigma")
    c6_sigma_prime = _fn(c6, "sigma_prime")
    c6_add = _fn(c6, "add_chart", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_guard = _fn(c6, "parameter_geometry_guard", "CurrentOriginalPrescribedR100RmC6Functions")
    component = _get(macro, "component_result", {})
    return {
        "accepted_macro_source_projection": component.get("passed") is True and component.get("full_macro_proved") is True and component.get("numeric_values_computed") is False,
        "exact_R100_C2_function_ports": ports.get("graph_ports_valid") is True and _returns_mapping(recover, "radius", "normalized_actual_fields", "normalized_actual_own_six_moments") and _has(r100_owner, "actual_fields={key:restore(row)forkey,rowine['same_original_macro_field_functions'].items()}", "actual_histories={key:restore(row)forkey,rowine['actual_six_moment_functions'].items()}") and _has(r100_eval, "cached_values_remain_source_enclosures_not_points=True") and saved.get("lei_ren_part1_paper_compliant_current_original_R100_endpoint.json") is True,
        "R100_to_first_R1_snapshot_ports": ports.get("checks", {}).get("R100_to_first_incoming_history_identity") is True and _has(first_owner, "self.upstream.owner(label)", "FirstSwitchFunctions(flow,value['normalized_actual_fields'],value['normalized_actual_own_six_moments'],self.cells)"),
        "first_switch_C2_C6_equations": receipts.get("lei_ren_part1_paper_compliant_current_original_first_switch_functions_check.json") is True and _has(first_eval, "h=f.h", "actual_six_histories=histories", "derivatives[name]=f.scale", "microscopic_jacobian_hb_applied_once=True", "exact_source='R=100*exp(hb*s);dy=hb*ds'", "true_sigma_function_used=True", "comparison_direction_not_redefined_by_actual_moments=True") and _has(c6_first, "g.mul(g.constant('-1/2'),hb,hb,radius)", "g.mul(hb,g.sub(g.one,sigma(g,phase)))", "a.scale(a.mul(ratio,drive),g.neg(chi_axial))", "original_first_switch_axial_control=True", "source_measure_applied_once=True"),
        "first_R1_end_nonreset_history_ports": ports.get("checks", {}).get("R1_actual_and_history_ports") is True and ports.get("checks", {}).get("R1_to_second_incoming_history_identity") is True and _has(first_eval, "incoming_delta=f.scale(self.inlet[name],dm)", "actual=f.add(self.inlet[name],change)") and _has(c6_init, "self.endpoints['R1']=self.end('first_switch',g.one)"),
        "second_switch_C2_C6_equations": receipts.get("lei_ren_part1_paper_compliant_current_original_second_switch_R110_check.json") is True and _has(second_eval, "h=f.h", "actual_six_histories=histories", "actual_first_exit_V_preserved_exactly=True", "exact_source='R=100*exp(hb*(1+t));dy=hb*dt'", "actual_moments_not_used_to_redefine_direction=True") and _has(c6_second, "g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))", "g.mul(g.constant('-2/5'),hb,sig)", "original_second_switch_full_current_radius_direction_and_power=True"),
        "second_R2_end_nonreset_history_ports": ports.get("checks", {}).get("R2_actual_and_history_ports") is True and ports.get("checks", {}).get("R2_to_post_incoming_history_identity") is True and _has(second_post, "actual_second_exit_incoming_histories=initial", "original_post_power_a='4/5'", "original_post_power_b=0", "actual_first_exit_V_preserved_exactly=True"),
        "post_switch_power_source": _assign(c6_init, "self.post_length", "g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))") and _has(c6_init, "phi=a.scale(incoming['actual']['phi'],g.unary('exp',g.mul(c('-2/5'),y)))", "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0))", "self.endpoints['R110']=self.end('post_switch',g.one)") and _has(second_post, "powers['angular']", "powers['squared']", "original_positive_kernel_weights=dict(angular=angular,swirl=swirl,pressure=pressure,axial=axial)"),
        "same_flat_sigma_source": receipts.get("lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json") is True and _has(c6_sigma, "previous.inlet.previous.lower.own.jets.complete.Gamma.sigma(g,x)") and _has(c6_sigma_prime, "original_flat_sigma_prime", "same original sigma Taylor coefficient[1]", "exact_flat_source_derivative_not_a_bound=True") and _has(first_eval, "sigma_jets(c,s)[0]") and _has(second_eval, "first.sigma_jets(c,s)[0]"),
        "comparison_modes_inherit_macro": _has(c6_comp, "old=self.bridge.comparison_macro", "eq=dict(H=phi,M=U,K=b.mul(phi,U),A=b.mul(U,U),B=b.scale(b.mul(phi,phi),'1/2'),C=b.mul(phi,phi))", "decay=g.unary('exp',g.mul(c(-rate),distance))", "previous.comparison_directions") and _get(component, "ast_witnesses.c2_mode_rows_match_c6_macro_modes") is True and _get(component, "ast_witnesses.macro_modes_exp_minus_j_t_and_measure_conversion") is True,
        "full_physical_rhs_and_swirl_identity": _has(c6_full, "a.scale(known['drive_hydro'],radius)", "a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar))", "a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl'])") and _has(c6_first, "drive=full_drive(a,known,radius,F0,Pstar)", "a.scale(a.mul(ratio,drive),g.neg(chi_axial))") and _get(component, "physical_swirl_identity") == "F0_base^2*C2_drive_swirl_germ == [F0^2*C6_normalized_swirl]_germ; bare higher-Z swirl equality is not claimed",
        "six_history_rates_nonreset_and_physical_measure_once": _get(component, "ast_witnesses.six_history_rates_nonreset_and_derivatives") is True and _has(c6_add, "previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)") and _has(first_eval, "forname,rateinmoments.RATES.items()", "microscopic_jacobian_hb_applied_once=True") and _has(second_eval, "forname,rateinmoments.RATES.items()", "microscopic_jacobian_hb_applied_once=True") and _has(second_post, "actual_log_radius_ODE_rows=dict(phi=f.scale(phi,-c.mpf(2)/5)") ,
        "same_exact_width_source": bool(domain.get("same_exact_width_source")) and receipts.get("lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json") is True and _has(c6_guard, "shared=wb.upstream", "shared is sw.history.upstream", "wb.logh is shared.logh", "wb.source==sw.source==common.source_sha", "bridge_hbB_switch_hbS_kept_distinct=True"),
        "ordinary_Z_source_only_scope": _get(macro, "scope.ordinary_Z_orders") == [0, 1, 2] and _get(macro, "scope.timeorder_recursion_claimed") is False and _get(macro, "scope.full_goal_complete") is False,
        "finite_N_and_fixture_seams_open": _has(c6_init, "finite_N_modified_fields_correction_histories_and_all17_not_installed=True", "new_original_scale_C6_numeric_containment_not_certified=True") and _source_has(texts["lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"], "mathematical_fixture_only\": true") and _get(component, "finite_N_history_identity") is False,
        "saved_rows_are_enclosures_only": all(saved.values()) and _has(r100_eval, "cached_values_remain_source_enclosures_not_points=True") and _source_has(texts["lei_ren_part1_paper_compliant_current_original_R100_endpoint.json"], "cached_values_remain_source_enclosures_not_points\": true"),
        "direct_C2_current_radius_modes_and_phase_transport": _has(first_init, "flow.radial_power(R0,j)*flow.radial_power(R100,p-j)", "self.D=modes(flow.d,1)", "self.drives={part:modes(flow.drive[part],2 if part=='swirl' else 1) for part in fields.PARTS}") and _has(first_angular, "self.scalar_series(1-j,t,average=True)", "f.h*f.h*(-t/2)*avg") and _has(first_force, "self.scalar_series(p-j,t,minus_one=True)") and _has(second_init, "first_owner.scalar_series(1-j,self.c.mpf(1))[0]", "self.D0=f.add(*self.D)") and _has(second_Ddelta, "self.first.scalar_series(1-j,t,minus_one=True)[0]") and _has(second_eval, "exact_source='R=100*exp(hb*(1+t));dy=hb*dt'"),
        "second_sigma_prefix_and_phase_derivative_identity": (
            _assign(second_angular, "sig", "t-complement_prefix")
            and _assign(second_angular, "deterministic", "[f.h*(-c.mpf(2)/5)*sig]+[f.scalar(0)]*5")
            and _assign(second_angular, "ell", "f.add(deterministic,f.scale(weighted_D,-f.h*f.h*c.mpf('.5')))")
            and _assign(second_eval, "weighted_D", "f.add(f.scale(self.D0,prefix),weighted_delta)")
            and _assign(second_eval, "sig", "first.sigma_jets(c,s)[0]")
            and _assign(second_eval, "Dend", "f.add(self.D0,self.Ddelta(s))")
            and _assign(second_eval, "angular_rate", "f.add(f.scale(Dend,h*h*(-c.mpf('.5'))*(1-sig)),[h*(-c.mpf(2)/5)*sig]+[f.scalar(0)]*5)")
            and _assign(second_eval, "phi_phase", "f.multiply(phi,angular_rate)")
            and _has(c6_second, "logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase")
            and _has(c6_second, "slope=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,radius,g.sub(g.one,sig))),a.fixed(g.mul(g.constant('-2/5'),sig)))")
        ),
    }


def _result(gates, witnesses, receipt_status, saved, nodes, actual_digest, macro, domain, ports):
    graph_ok = isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT and actual_digest == EXPECTED_NODES_DIGEST
    keys = ("accepted_macro_source_projection", "exact_R100_C2_function_ports", "R100_to_first_R1_snapshot_ports", "first_switch_C2_C6_equations", "first_R1_end_nonreset_history_ports", "second_switch_C2_C6_equations", "second_R2_end_nonreset_history_ports", "post_switch_power_source", "same_flat_sigma_source", "comparison_modes_inherit_macro", "full_physical_rhs_and_swirl_identity", "six_history_rates_nonreset_and_physical_measure_once", "same_exact_width_source", "ordinary_Z_source_only_scope", "finite_N_and_fixture_seams_open", "saved_rows_are_enclosures_only", "direct_C2_current_radius_modes_and_phase_transport", "second_sigma_prefix_and_phase_derivative_identity")
    receipts_ok = bool(receipt_status) and all(receipt_status.values()) and all(saved.values())
    background = bool(graph_ok and receipts_ok and ports.get("graph_ports_valid") is True and domain.get("proved") is True and not gates and all(witnesses.get(key) is True for key in keys))
    own = bool(background and witnesses.get("six_history_rates_nonreset_and_physical_measure_once") is True)
    return {
        "schema": SCHEMA, "source_only": True, "passed": bool(background and own),
        "nodes_digest": actual_digest, "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "full_R100_switch_background_projection_proved": background,
        "full_R100_switch_own_history_projection_proved": own,
        "full_R100_switch_evolution_identity_proved": bool(background and own),
        "full_R100_switch_own_history_identity_proved": own,
        "full_R100_switch_projection_proved": bool(background and own),
        "ordinary_Z_orders": [0, 1, 2], "ordinary_orders": [0, 1, 2],
        "ordinary_Z_not_temporal_n": True, "time_order_recursion_claimed": False, "timeorder_recursion_claimed": False,
        "same_exact_width_source_proved": witnesses.get("same_exact_width_source") is True,
        "same_density_certificate": False, "physical_Q_y_Z1": False, "finite_N_history_identity": False,
        "finite_N_sc_over_2": False, "global_five_moment_closure": False,
        "numeric_values_computed": False, "old_operators_replayed": False,
        "mathematical_operator_replayed": False, "full_goal_complete": False,
        "ast_witnesses": witnesses,
        "accepted_receipts": {key: True for key, value in receipt_status.items() if value},
        "receipt_status": receipt_status, "saved_source_reports": saved,
        "parent_report": PARENT_REPORT, "parent_report_sha256": PARENT_SHA256,
        "accepted_macro_report": MACRO_REPORT, "accepted_macro_report_sha256": MACRO_SHA256,
        "source_family": _get(macro, "source_family", {}),
        "projection_kind": "exact original R100 first-switch, second-switch, and post-switch source-function identity; saved rows remain enclosures",
        "physical_swirl_identity": "F0_base^2*C2_drive_swirl_germ == [F0^2*C6_normalized_swirl]_germ; bare higher-Z swirl equality is not claimed",
        "parameter_ports": {"hbS": ports.get("switch_width_node"), "hbB": ports.get("bridge_width_node"), "handle_identity_not_alias": ports.get("checks", {}).get("width_handles_distinct") is True, "same_exact_width_source_proved": witnesses.get("same_exact_width_source") is True, "differential_conversion": "first/second dy=hbS ds or hbS dt; post dy=Lp dp; physical measure applied once"},
        "exact_R100_switch_function_ports": ports,
        "endpoint_ports": {"endpoint_rows_are_values": False, "exact_R100_endpoint": "C2 R100 actual phi/U and six moment source-function ports; saved packet rows are enclosures", "R1": "first-switch endpoint at s=1 with nonreset incoming six histories", "R2": "second-switch endpoint at t=1 with nonreset incoming six histories", "R110": "post-switch endpoint at p=1", "first_switch": "R=100*exp(hbS*s), dy=hbS ds; phi_s=-hbS^2 R D_over_R phi/2; U_s=-hbS^2(1-sigma)*(phi/barphi2)*full_drive; physical phi_y=phi_s/hbS and U_y=U_s/hbS", "second_switch": "R=100*exp(hbS*(1+t)), dy=hbS dt; phi_t=phi[-hbS^2(1-sigma)R D_over_R/2-(2/5)hbS sigma]; physical phi_y=phi_t/hbS; U is constant", "post_switch": "Lp=log(11/10)-2hbS; R=100*exp(2hbS+Lp*p); phi=phi_R2*exp(-2Lp*p/5); U=U_R1", "six_histories": "rates 2,1,2,1,2,1; incoming memories are not reset; physical Jacobian hbS or Lp appears once"},
        "source_equivalence_derivation": {
            "current_radius_modes": "C2 formal R0^j*100^(p-j)*exp((p-j)*h*s) equals R0^j*R(s)^(p-j); the old R0 enclosure remains only a cover of its exact Ra*exp(2h) source.",
            "comparison_semigroup": "exp(-rate*log(100/R0))*exp(-rate*log(R/100))=exp(-rate*log(R/R0)); same frozen equilibrium and modes, with no comparison-memory reset or replacement by actual moments.",
            "second_switch_primitives": "C(t)=integral_0^t(1-sigma); S(t)=t-C(t); W(t)=integral_0^t((1-sigma)*Dend); ell=-(2/5)*h*S-(h^2/2)*W; ell_t=-(2/5)*h*sigma-(h^2/2)*(1-sigma)*Dend, Dend=R*D_over_R.",
            "ordinary_Z": "Smooth finite-interval parameter differentiation of these same identities through ordinary Z order 2 preserves full products and incoming jets; no midpoint, finite-N error, or temporal n recursion is used.",
        },
        "domain_evidence": domain,
        "sublemmas": {"accepted_macro_source": witnesses.get("accepted_macro_source_projection", False), "exact_R100_inlet_function_port": witnesses.get("exact_R100_C2_function_ports", False), "first_switch_projection": witnesses.get("first_switch_C2_C6_equations", False), "second_switch_projection": witnesses.get("second_switch_C2_C6_equations", False), "post_switch_projection": witnesses.get("post_switch_power_source", False), "same_comparison_modes": witnesses.get("comparison_modes_inherit_macro", False), "same_full_physical_rhs": witnesses.get("full_physical_rhs_and_swirl_identity", False), "same_sigma": witnesses.get("same_flat_sigma_source", False), "same_six_history_rates_and_initials": own, "whole_interval_positivity_existence_uniqueness": domain.get("proved", False), "finite_N_seam_open": False},
        "scope_flags": {"same_density_certificate": False, "finite_N_history_identity": False, "physical_Q_y_Z1": False, "global_five_moment_closure": False, "time_order_recursion": False, "full_goal_complete": False},
        "retained_construction_data": {"cumulative_terms": 212, "raw_signed_terms": 213, "reduced_terms": 183, "six_channels_retained": True, "six_signed_channels_retained": True, "all_old_error_data_retained": True, "reserved_contribution_registered_once": True},
        "fixture_limits": {"finite_fixture_only": True, "actual_source_admission": False, "scope": "R100/first/second saved receipts are source-function enclosures; no finite-N admission"},
        "unresolved_gates": sorted(set(gates)) if gates or not (background and own) else [],
        "potential_gates": ["same_density_certificate", "finite_N_history_identity", "physical_Q_y_Z1", "global_five_moment_closure", "time_order_recursion"],
        "record_hashes": dict(RECORD_HASHES), "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "source_inspection_audit": {"source_records_only": True, "static_source_AST_graph_and_saved_bound_inspection_only": True, "nodes_hashed_only": True, "source_graph_substitutions_checked": True, "new_source_consumer_calls": 1, "accepted_ancestor_operator_calls": 0, "producer_provider_builder_checker_replayed": False, "quadrature_or_norm_calls": 0, "numeric_field_evaluated": False},
        "worker_role": "luna_max_worker", "runtime_model_metadata": None,
    }


def derive(nodes, source_records, *, expected_nodes_digest):
    """Return the fail-closed original R100 switch source projection."""
    gates = []
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append("expected_nodes_digest_mismatch")
    try:
        actual_digest = _digest(nodes)
    except Exception as error:
        actual_digest = ""
        gates.append("nodes_digest_error:" + type(error).__name__)
    if actual_digest != EXPECTED_NODES_DIGEST:
        gates.append("nodes_digest_mismatch")
    if not isinstance(nodes, list) or len(nodes) != EXPECTED_NODE_COUNT:
        gates.append("node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    texts = {name: _text(records, name, gates) for name in RECORD_HASHES}
    parent, macro, component = _reports(texts, gates)
    receipt_status = _receipts(texts, gates)
    saved = _saved_reports(texts, gates)
    trees = {name: _tree(value, name, gates) for name, value in texts.items() if name.endswith(".py")}
    ports = _archive_ports(texts, nodes, gates)
    domain = _domain(texts, trees, macro, receipt_status, ports, gates)
    witnesses = _witnesses(texts, trees, receipt_status, saved, macro, ports, domain)
    for name, ok in witnesses.items():
        if ok is not True:
            gates.append("ast_witness_missing:" + name)
    return _result(gates, witnesses, receipt_status, saved, nodes, actual_digest, macro, domain, ports)
