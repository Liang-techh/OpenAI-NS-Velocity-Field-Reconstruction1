"""Static source projection from the accepted R100 switch stage to ``Rm``.

This module authenticates the saved source packets and checks source AST/data
flow witnesses.  It deliberately does not import or execute a producer,
provider, evaluator, builder, checker, recurrence, quadrature, norm, or
project target.  Saved graph nodes and authenticated source snapshots are
inspected as data only; enclosure rows remain enclosures and no field value is
selected.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping


SCHEMA = "original_R110_Rm_C2_C6_source_projection_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
SWITCH_REPORT = "lei_ren_part1_paper_compliant_current_original_R100_switch_source_projection.json"
SWITCH_SHA256 = "92581da55a9da2a7dec6dbca76a433ce988f15786092b0ba35897b39bdcfb9d4"

C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
C6_DECODED_SHA256 = "18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a"
C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
C6_CHECK = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"
C6_CHECK_SHA256 = "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"

# These two inherited owner sources carry the physical generic products used
# below.  They are pinned directly so a metadata label cannot stand in for the
# source-level W/FF and Pstar identities.
BRIDGE_C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
BRIDGE_C6_SOURCE_SHA256 = "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605"
CORE_INLET_SOURCE = "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
CORE_INLET_SOURCE_SHA256 = "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890"

LONG_SOURCE = "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py"
LONG_SOURCE_SHA256 = "9b41c9a587c52f3af5de35a39695c16e07de509155ce466680d984cabd54b231"
LONG_REPORT = "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.json"
LONG_REPORT_SHA256 = "fa243f3390c1d910eeb6231876e8ec37cec1404e9d5e05dadce882ebf0d24252"
LONG_CHECK = "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint_check.json"
LONG_CHECK_SHA256 = "161b972310323e9e0dc7cac4956127b666f9f5b4739b1aadf65cde83c50a2f39"

REF_SOURCE = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py"
REF_SOURCE_SHA256 = "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7"
REF_REPORT = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.json"
REF_REPORT_SHA256 = "47f21c94f6bb08ebef01288d7b506ef5135ecd2cf5e3fba9cfd15567f6cef440"
REF_CHECK = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions_check.json"
REF_CHECK_SHA256 = "94f5b78791bd571ae6759d2d53e6fdda47607b047138e9c4868e2611dff10b02"

LONG_PROFILE_SOURCE = "lei_ren_part1_paper_compliant_long_reshape_profiles.py"
LONG_PROFILE_SOURCE_SHA256 = "b12c7ef683033e3a0612ede04720d37b625a116ec8231c1d9a7a1806ae05c4c3"
LONG_PROFILE_REPORT = "lei_ren_part1_paper_compliant_long_reshape_profiles.json"
LONG_PROFILE_REPORT_SHA256 = "ae2612d307411e2eb8e2c4a93dc8f5b1ee834601d0360eea56f322e346f7f1bc"
LONG_PROFILE_CHECK = "lei_ren_part1_paper_compliant_long_reshape_profiles_check.json"
LONG_PROFILE_CHECK_SHA256 = "3f8ed1d3239da27b88279c90356aa4aa885690dee1a4d414fbc927f7f38a58dc"

REF_PROFILE_SOURCE = "lei_ren_part1_paper_compliant_reference_restore_profiles.py"
REF_PROFILE_SOURCE_SHA256 = "bee442990e9e52b04f37063d21e2da80ed5c35ec037941ac2af5b63860308765"
REF_PROFILE_REPORT = "lei_ren_part1_paper_compliant_reference_restore_profiles.json"
REF_PROFILE_REPORT_SHA256 = "8b21a26b79401a27ff174c68a670d08cbb10ff094ce34ac5772077449aaeb8f6"
REF_PROFILE_CHECK = "lei_ren_part1_paper_compliant_reference_restore_profiles_check.json"
REF_PROFILE_CHECK_SHA256 = "fecc91c0ba8368c681a2631db174b7f25da82ee16ceace2a744fb3f545f15087"

MOMENTS_SOURCE = "lei_ren_part1_paper_compliant_current_original_bridge_macro_moments.py"
MOMENTS_SOURCE_SHA256 = "ab0f25f7a9c9436f8fbaba9952d538a4c861695911d4e6596eebc219ac124b85"
K1_LEDGER = "lei_ren_part1_paper_compliant_K1_ledger.json"
K1_LEDGER_SHA256 = "38cd8e76c0d5f5519b510d9cb0502b8915446ffa1c6338d3a3a98c4aca2cdf7d"

FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
REFERENCE_JOIN_FAMILY_SHA256 = "aa754c733df332110edae66247ccb83a813c3c3a6f9ed90fac89cf07054f691f"
HISTORICAL_PROFILE_PROVENANCE = {
    "actual_five_defect_family_sha256": FAMILY["actual_five_defect_family_sha256"],
    "implicit_source_sha256": FAMILY["implicit_source_sha256"],
    "reference_join_family_sha256": REFERENCE_JOIN_FAMILY_SHA256,
}

CHARTS = ("first_switch", "second_switch", "post_switch", "long_reshape", "reference", "restoration", "post_restore")
ENDPOINTS = ("R1", "R2", "R110", "Rsh", "Rz", "Rrestore", "Rm")
RAW_HISTORY_RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
NORMALIZED_SHAPE_RATES = {
    "theta": "8/5",
    "theta_z": "8/5",
    "pressure": "1/5",
    "swirl": "6/5",
    "mean": "1",
    "axial": "1",
}

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    SWITCH_REPORT: SWITCH_SHA256,
    C6_ARCHIVE: C6_ARCHIVE_SHA256,
    C6_SOURCE: C6_SOURCE_SHA256,
    C6_CHECK: C6_CHECK_SHA256,
    BRIDGE_C6_SOURCE: BRIDGE_C6_SOURCE_SHA256,
    CORE_INLET_SOURCE: CORE_INLET_SOURCE_SHA256,
    LONG_SOURCE: LONG_SOURCE_SHA256,
    LONG_REPORT: LONG_REPORT_SHA256,
    LONG_CHECK: LONG_CHECK_SHA256,
    REF_SOURCE: REF_SOURCE_SHA256,
    REF_REPORT: REF_REPORT_SHA256,
    REF_CHECK: REF_CHECK_SHA256,
    LONG_PROFILE_SOURCE: LONG_PROFILE_SOURCE_SHA256,
    LONG_PROFILE_REPORT: LONG_PROFILE_REPORT_SHA256,
    LONG_PROFILE_CHECK: LONG_PROFILE_CHECK_SHA256,
    REF_PROFILE_SOURCE: REF_PROFILE_SOURCE_SHA256,
    REF_PROFILE_REPORT: REF_PROFILE_REPORT_SHA256,
    REF_PROFILE_CHECK: REF_PROFILE_CHECK_SHA256,
    MOMENTS_SOURCE: MOMENTS_SOURCE_SHA256,
    K1_LEDGER: K1_LEDGER_SHA256,
}
DECODED_RECORD_HASHES = {C6_ARCHIVE: C6_DECODED_SHA256}
PREFIX_NODES = 97295


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def _get(value, path, default=None):
    for part in path.split(".") if path else ():
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _norm(value):
    return "".join(value.split())


def _text(records, name, gates):
    record = records.get(name)
    if not isinstance(record, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        gates.append("source_path_mismatch:" + name)
    if record.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
    value = record.get("source_text")
    if not isinstance(value, str):
        gates.append("source_text_missing:" + name)
        return ""
    if name in DECODED_RECORD_HASHES:
        if record.get("encoding") != "gzip" or record.get("decoded_sha256") != DECODED_RECORD_HASHES[name]:
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        try:
            # Authenticate the raw gzip bytes before decompression.
            raw = base64.b64decode(record.get("raw_bytes_base64", ""), validate=True)
            raw_sha = hashlib.sha256(raw).hexdigest()
            if raw_sha != RECORD_HASHES[name]:
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


def _json(value, name, gates):
    try:
        return json.loads(value, parse_int=int, parse_float=float)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _tree(value, name, gates):
    try:
        return ast.parse(value, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _fn(tree, name, cls=None):
    if cls is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                return next((child for child in node.body if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name), None)
        return None
    return next((node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def _has(fn, *needles):
    if fn is None:
        return False
    compact = _norm(ast.unparse(fn))
    return all(_norm(needle) in compact for needle in needles)


def _source_has(text, *needles):
    compact = _norm(text)
    return all(_norm(needle) in compact for needle in needles)


def _assign(fn_or_tree, target, expression):
    """Compare assignment ASTs, preserving structural equality under reboxing."""
    if fn_or_tree is None:
        return False
    wanted_target = _norm(target)
    try:
        wanted_value = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    except Exception:
        return False
    for node in ast.walk(fn_or_tree):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if any(_norm(ast.unparse(item)) == wanted_target for item in targets) and ast.dump(node.value, include_attributes=False) == wanted_value:
            return True
    return False


def _family_ok(row):
    family = row.get("source_family") if isinstance(row, Mapping) else None
    return family == FAMILY or family == FAMILY["actual_five_defect_family_sha256"] or (
        isinstance(row, Mapping)
        and row.get("actual_five_defect_family_sha256") == FAMILY["actual_five_defect_family_sha256"]
        and row.get("implicit_source_sha256") == FAMILY["implicit_source_sha256"]
        and row.get("datum_enclosure_sha256") == FAMILY["datum_enclosure_sha256"]
    )


def _accepted_switch(texts, gates):
    parent = _json(texts[PARENT_REPORT], PARENT_REPORT, gates)
    switch = _json(texts[SWITCH_REPORT], SWITCH_REPORT, gates)
    if not _family_ok(parent) or not _family_ok(switch):
        gates.append("accepted_source_family_changed")
    if switch.get("parent_report") != PARENT_REPORT or switch.get("parent_report_sha256") != PARENT_SHA256:
        gates.append("accepted_switch_parent_pin_changed")
    component = switch.get("component_result", {})
    required = {
        "passed": True,
        "source_only": True,
        "full_R100_switch_background_projection_proved": True,
        "full_R100_switch_own_history_projection_proved": True,
        "same_exact_width_source_proved": True,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
    }
    for key, expected in required.items():
        if component.get(key) is not expected:
            gates.append("accepted_switch_component_changed:" + key)
    if component.get("nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("accepted_switch_graph_changed")
    scope = switch.get("scope", {})
    scope_required = {
        "ordinary_Z_orders": [0, 1, 2],
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "full_goal_complete": False,
    }
    for key, expected in scope_required.items():
        if scope.get(key) != expected:
            gates.append("accepted_switch_scope_changed:" + key)
    ports = component.get("exact_R100_switch_function_ports", {})
    if not isinstance(ports, Mapping):
        ports = {}
    for key, expected in {
        "graph_ports_valid": True,
        "endpoint_rows_are_values": False,
        "saved_C2_rows_are_enclosures": True,
        "numeric_values_computed": False,
        "ordinary_Z_orders_admitted": [0, 1, 2],
    }.items():
        if ports.get(key) != expected:
            gates.append("accepted_switch_ports_changed:" + key)
    for key in ("C2_R100_actual_function_ports", "C2_R100_history_function_ports", "R1_end_function_ports", "R2_end_function_ports", "R110_end_function_ports"):
        if not isinstance(ports.get(key), Mapping):
            gates.append("accepted_switch_port_missing:" + key)
    audit = switch.get("source_inspection_audit", {})
    for key, expected in {
        "source_records_only": True,
        "static_source_AST_graph_and_saved_bound_inspection_only": True,
        "nodes_hashed_only": True,
        "source_graph_substitutions_checked": True,
        "new_source_consumer_calls": 1,
        "accepted_ancestor_operator_calls": 0,
        "quadrature_or_norm_calls": 0,
        "numeric_field_evaluated": False,
    }.items():
        if audit.get(key) != expected:
            gates.append("accepted_switch_audit_changed:" + key)
    accepted = component.get("passed") is True and component.get("source_only") is True and not any(
        gate.startswith(("accepted_switch_", "accepted_source_family")) for gate in gates
    )
    return dict(parent=parent, switch=switch, component=component, scope=scope, ports=ports, accepted=accepted)


def _receipt(receipts, name, required, forbidden, gates, *, historical_profile=False):
    row = receipts.get(name, {})
    good = row.get("all_passed") is True
    if not good:
        gates.append("receipt_not_accepted:" + name)
    profile_family_ok = historical_profile and all(
        row.get(key) == value for key, value in HISTORICAL_PROFILE_PROVENANCE.items()
    )
    if not _family_ok(row) and not profile_family_ok:
        gates.append("receipt_source_family_mismatch:" + name)
        good = False
    if historical_profile and not profile_family_ok:
        gates.append("historical_profile_provenance_mismatch:" + name)
    for path in required:
        marker = _get(row, path)
        ok = marker is True or (type(marker) is int and marker > 0)
        if not ok:
            gates.append("receipt_missing:" + name + ":" + path)
        good = good and ok
    for path in forbidden:
        if _get(row, path) is True:
            gates.append("receipt_open_seam_claimed:" + name + ":" + path)
            good = False
    return good


def _receipts(texts, gates):
    receipts = {name: _json(text, name, gates) for name, text in texts.items() if name.endswith("_check.json")}
    specs = {
        C6_CHECK: (
            (
                "checks.actual_source.passed",
                "checks.actual_source.all_genuine_seven_source_and_six_generic_rows",
                "checks.actual_source.every_nonzero_incoming_history_function_retained_by_identity",
                "checks.actual_source.endpoint_phi_U_and_physical_y_derivatives_exposed",
                "checks.actual_source.distinct_width_tags_with_exact_shared_upstream_K1_identity",
                "checks.actual_source.frozen_T_and_positive_gap_post_length_guarded",
                "checks.actual_source.no_free_core_rho_or_old_phase_capture_in_chart_bodies",
                "checks.actual_source.finite_N_modified_source_and_quantitative_C6_not_claimed",
                "checks.geometry.passed",
                "checks.geometry.exact_radius_endpoint_identities",
                "checks.geometry.exact_physical_log_radius_Jacobians",
                "checks.geometry.original_finite_T_and_Rref_offsets_retained",
                "checks.geometry.complete_radius_order_proved_by_actual_source_positive_guards",
                "checks.long_restore_joins.passed",
                "checks.long_restore_joins.nonconstant_actual_R110_B_and_velocity_inputs",
                "checks.long_restore_joins.full_finite_T_and_offsets_minus8_minus7_minus6",
                "checks.long_restore_joins.restoration_end_actual_U_is_4Z_through_sixth_derivative",
                "checks.long_restore_joins.mathematical_fixture_only_not_original_scale_numeric_accuracy",
                "checks.switch_controls.passed",
                "checks.switch_controls.first_axial_hb_squared_one_minus_sigma_and_current_R_Rsquared_preserved",
                "checks.switch_controls.second_angular_minus_two_fifths_sigma_and_constant_U_preserved",
                "checks.switch_controls.mathematical_fixture_only",
            ),
            (
                "current_original_actual_modified_R100_Rm_C6_functions_installed",
                "current_original_prescribed_R100_Rm_C6_numeric_containment_installed",
                "temporal_recursion",
                "full_corrected_Navier_Stokes_solution",
            ),
        ),
        LONG_CHECK: (
            (
                "original_actual_two_frame_long_reshape_terminal_history_functions_installed",
                "exact_source_transport_identities.passed",
                "genuine_original_terminal_source_contracts.passed",
                "guards.passed",
            ),
            (
                "actual_active_patch_genuine_coefficients_installed",
                "actual_complete_patch_integral_installed",
                "numerical_original_source_point_or_integral_oracle_installed",
                "actual_coefficient_recursion_certified",
            ),
        ),
        REF_CHECK: (
            (
                "original_actual_two_frame_Rsh_reference_restoration_Rm_functions_installed",
                "exact_transport_and_source_identities.passed",
                "genuine_actual_source_joins_and_histories.passed",
                "guards.passed",
            ),
            (
                "actual_active_patch_genuine_coefficients_installed",
                "actual_complete_patch_integral_installed",
                "numerical_original_source_point_or_integral_oracle_installed",
                "actual_coefficient_recursion_certified",
            ),
        ),
        LONG_PROFILE_CHECK: (
            (
                "kernel_fixture.passed",
                "structural_identities.passed",
                "original_selected_A_T_retained",
                "exact_R110_histories_retained",
                "actual_B_C2_shear_sigma_and_kernel_rate_theorems_directly_bound",
                "actual_long_reshape_log_velocity_axial5_available",
                "actual_long_reshape_moment_axial5_enclosures_available",
                "actual_Rsh_exit_axial5_available",
            ),
            (
                "actual_source_admission",
                "reshape_radial_mixed4_certified",
                "reference_continuation_installed",
                "actual_five_moment_patch_connected",
                "temporal_recursion",
            ),
        ),
        REF_PROFILE_CHECK: (
            (
                "restoration_fixture.passed",
                "structural_identities.passed",
                "actual_source_binding.passed",
                "exact_Rsh_Rz_histories_and_P0_retained",
                "original_axial_restoration_reaches_exact_4Z",
                "negative_inherited_defect_terms_retained",
                "actual_reference_continuation_axial5_enclosures_available",
                "actual_axial_restoration_axial5_enclosures_available",
                "actual_unpatched_Rm_Rh_moment_axial5_enclosures_available",
                "actual_five_defect_profile_axial5_enclosures_available",
            ),
            (
                "actual_source_admission",
                "radial_mixed4_certified",
                "all_radial_endpoint_joins_certified",
                "actual_five_moment_patch_connected",
                "temporal_recursion",
            ),
        ),
    }
    status = {}
    historical_profiles = {LONG_PROFILE_CHECK, REF_PROFILE_CHECK}
    for name, (required, forbidden) in specs.items():
        status[name] = _receipt(
            receipts, name, required, forbidden, gates,
            historical_profile=name in historical_profiles,
        )
    # Current accepted endpoint owners must pin the historical profile source,
    # report and check bytes explicitly.  This keeps the profile packets in a
    # separate algorithm/fixture provenance scope and never promotes their
    # numerical rows into the current physical family.
    profile_links = {
        "long": (LONG_CHECK, (LONG_PROFILE_SOURCE, LONG_PROFILE_REPORT, LONG_PROFILE_CHECK)),
        "reference": (REF_CHECK, (REF_PROFILE_SOURCE, REF_PROFILE_REPORT, REF_PROFILE_CHECK)),
    }
    for label, (owner, names) in profile_links.items():
        owner_row = receipts.get(owner, {})
        owner_inputs = owner_row.get("input_hashes", {}) if isinstance(owner_row, Mapping) else {}
        linked = all(owner_inputs.get(name) == RECORD_HASHES[name] for name in names)
        status["historical_profile_" + label + "_linkage"] = linked
        if not linked:
            gates.append("historical_profile_input_linkage_missing:" + label)
    return status


def _archive_ports(texts, nodes, gates):
    archive = _json(texts[C6_ARCHIVE], C6_ARCHIVE, gates)
    bindings = archive.get("prescribed_R100_Rm_C6_bindings", {})
    charts = bindings.get("charts", {})
    endpoints = bindings.get("endpoints", {})
    saved_nodes = archive.get("exact_graph_nodes", [])
    chart_endpoint_pairs = {
        "R1": "first_switch",
        "R2": "second_switch",
        "R110": "post_switch",
        "Rsh": "long_reshape",
        "Rz": "reference",
        "Rrestore": "restoration",
        "Rm": "post_restore",
    }
    incoming_previous = {
        "second_switch": "R1",
        "post_switch": "R2",
        "long_reshape": "R110",
        "reference": "Rsh",
        "restoration": "Rz",
        "post_restore": "Rrestore",
    }

    def all_graph_ids(value):
        if type(value) is int:
            return 0 <= value < len(saved_nodes)
        if isinstance(value, list):
            return all(all_graph_ids(item) for item in value)
        if isinstance(value, Mapping):
            return all(all_graph_ids(item) for item in value.values())
        return True

    def graph_row(node_id):
        if type(node_id) is not int or not 0 <= node_id < len(saved_nodes):
            return None
        return saved_nodes[node_id]

    exact_chart_nodes = {
        name: {
            key: charts[name].get(key)
            for key in ("phase", "radius", "physical_measure", "actual", "incoming_histories", "moments", "fields", "radial", "generic")
        }
        for name in CHARTS
        if isinstance(charts.get(name), Mapping)
    }
    exact_endpoint_nodes = {
        name: {
            key: endpoints[name].get(key)
            for key in ("radius", "actual", "moments", "fields", "generic")
        }
        for name in ENDPOINTS
        if isinstance(endpoints.get(name), Mapping)
    }
    incoming_history_identity = {
        chart: isinstance(charts.get(chart), Mapping)
        and isinstance(endpoints.get(previous), Mapping)
        and charts[chart].get("incoming_histories") == endpoints[previous].get("moments")
        for chart, previous in incoming_previous.items()
    }
    endpoint_substitution_identity = {}
    for endpoint, chart in chart_endpoint_pairs.items():
        row = endpoints.get(endpoint, {})
        radius_id = row.get("radius") if isinstance(row, Mapping) else None
        graph = graph_row(radius_id)
        endpoint_substitution_identity[endpoint] = (
            isinstance(graph, Mapping)
            and graph.get("operation") == "function_substitution"
            and graph.get("expression") == charts.get(chart, {}).get("radius")
            and graph.get("variable") == charts.get(chart, {}).get("phase")
            and graph.get("value") == 1
            and graph.get("Z_independent_substitution") is True
        )
    expected_measures = {
        "first_switch": "hb",
        "second_switch": "hb",
        "post_switch": "post_length",
        "long_reshape": "T",
        "reference": "reference_length",
        "restoration": 1,
        "post_restore": 1,
    }
    fixed_measure_identity = {
        chart: charts.get(chart, {}).get("physical_measure") == (
            bindings.get(want) if isinstance(want, str) else want
        )
        for chart, want in expected_measures.items()
    }
    parameter_nodes = {
        key: bindings.get(key)
        for key in ("T", "F0", "P0", "Pstar", "hb", "post_length", "reference_length", "logref", "B", "E", "logq")
    }
    parameter_nodes_complete = all(
        isinstance(parameter_nodes.get(key), (int, list, tuple))
        for key in ("T", "F0", "P0", "Pstar", "hb", "post_length", "reference_length", "logref", "B", "E", "logq")
    )
    source_binding_identity = parameter_nodes_complete and all(all_graph_ids(value) for value in parameter_nodes.values())
    graph_id_integrity = all(all_graph_ids(value) for value in exact_chart_nodes.values()) and all(all_graph_ids(value) for value in exact_endpoint_nodes.values())
    endpoint_chain_shape = all(
        isinstance(exact_endpoint_nodes.get(endpoint), Mapping)
        and set(exact_endpoint_nodes[endpoint].get("actual", {})) == {"phi", "U", "phi_y", "U_y"}
        and set(exact_endpoint_nodes[endpoint].get("moments", {})) == {"H", "M", "K", "A", "B", "C"}
        for endpoint in ENDPOINTS
    )
    checks = {
        "source_family": archive.get("source_family") == FAMILY,
        "prefix_length": isinstance(saved_nodes, list) and len(saved_nodes) == PREFIX_NODES,
        "graph_prefix_matches_snapshot": isinstance(nodes, list) and isinstance(saved_nodes, list) and nodes[:PREFIX_NODES] == saved_nodes,
        "charts_present": all(isinstance(charts.get(name), Mapping) for name in CHARTS),
        "endpoints_present": all(isinstance(endpoints.get(name), Mapping) for name in ENDPOINTS),
        "chart_ports_complete": all(
            isinstance(charts.get(name), Mapping)
            and all(key in charts[name] for key in ("phase", "radius", "physical_measure", "actual", "incoming_histories", "moments", "fields", "radial", "generic"))
            for name in CHARTS
        ),
        "endpoint_ports_complete": all(
            isinstance(endpoints.get(name), Mapping)
            and all(key in endpoints[name] for key in ("radius", "actual", "moments", "fields", "generic"))
            for name in ENDPOINTS
        ),
        "geometry_certificate": isinstance(bindings.get("geometry_certificate"), Mapping),
        "definition": isinstance(bindings.get("definition"), Mapping),
        "graph_id_integrity": graph_id_integrity,
        "source_binding_identity": source_binding_identity,
        "parameter_nodes_complete": parameter_nodes_complete,
        "incoming_history_identity": all(incoming_history_identity.values()),
        "endpoint_substitution_identity": all(endpoint_substitution_identity.values()),
        "fixed_measure_identity": all(fixed_measure_identity.values()),
        "endpoint_chain_shape": endpoint_chain_shape,
    }
    geometry = bindings.get("geometry_certificate", {})
    definition = bindings.get("definition", {})
    for key in (
        "shared_upstream_is_same_object",
        "bridge_logh_is_actual_upstream_object",
        "switch_width_is_exact_source_reboxing",
        "all_four_native_T_source_tuples_equal",
        "actual_native_T_and_logC_match_common_source",
        "bridge_hbB_switch_hbS_kept_distinct",
        "directed_T_positive",
        "directed_reference_gap_positive",
        "directed_post_length_positive",
        "bounds_are_not_defining_field_parameters",
    ):
        checks["geometry_" + key] = geometry.get(key) is True
    for key in (
        "true_prescribed_background_before_finite_N_loop",
        "full_original_finite_T_400A_not_shortened_or_infinite",
        "original_B_from_actual_R110_and_anchored_G",
        "all_six_nonzero_incoming_and_every_chart_memory_retained",
        "original_axial_restoration_4Z_plus_actual_E_times_one_minus_sigma",
        "original_independent_P0_and_physical_Pstar_F0_factors_once",
        "all_seven_source_rows_and_six_generic_output_rows_genuine",
        "finite_N_modified_fields_correction_histories_and_all17_not_installed",
        "new_original_scale_C6_numeric_containment_not_certified",
    ):
        checks["definition_" + key] = definition.get(key) is True
    checks["offsets"] = definition.get("exact_Rref_Rz_Rrestore_Rm_offsets") == [-8, -7, -6]
    checks["radius_order"] = geometry.get("complete_radius_order") == "100<R1<R2<110<Rsh<Rz<Rrestore<Rm"
    for name, ok in checks.items():
        if ok is not True:
            gates.append("archive_port_missing:" + name)
    graph_ports_valid = all(value is True for value in checks.values())
    chart_shapes = {
        name: {
            "keys": sorted(charts[name]),
            "actual_keys": sorted(charts[name].get("actual", {})) if isinstance(charts[name].get("actual"), Mapping) else [],
            "incoming_history_keys": sorted(charts[name].get("incoming_histories", {})) if isinstance(charts[name].get("incoming_histories"), Mapping) else [],
            "moment_keys": sorted(charts[name].get("moments", {})) if isinstance(charts[name].get("moments"), Mapping) else [],
        }
        for name in CHARTS
        if isinstance(charts.get(name), Mapping)
    }
    endpoint_shapes = {
        name: {
            "keys": sorted(endpoints[name]),
            "actual_keys": sorted(endpoints[name].get("actual", {})) if isinstance(endpoints[name].get("actual"), Mapping) else [],
            "moment_keys": sorted(endpoints[name].get("moments", {})) if isinstance(endpoints[name].get("moments"), Mapping) else [],
        }
        for name in ENDPOINTS
        if isinstance(endpoints.get(name), Mapping)
    }
    return {
        "source_archive": C6_ARCHIVE,
        "source_archive_sha256": C6_ARCHIVE_SHA256,
        "decoded_source_archive_sha256": C6_DECODED_SHA256,
        "graph_ports_valid": graph_ports_valid,
        "endpoint_rows_are_values": False,
        "chart_ports": chart_shapes,
        "endpoint_ports": endpoint_shapes,
        # Preserve the authenticated node identities themselves.  The shape
        # summaries above are useful for inspection, but cannot establish
        # that the accepted R110 ports and the four incoming substitutions
        # are the exact saved graph objects.
        "exact_chart_nodes": exact_chart_nodes,
        "exact_endpoint_nodes": exact_endpoint_nodes,
        "parameter_nodes": parameter_nodes,
        "incoming_history_identity": incoming_history_identity,
        "endpoint_substitution_identity": endpoint_substitution_identity,
        "fixed_measure_identity": fixed_measure_identity,
        "charts": list(CHARTS),
        "endpoints": list(ENDPOINTS),
        "geometry_certificate": {key: geometry.get(key) for key in geometry if isinstance(geometry.get(key), (bool, str))},
        "definition": {key: definition.get(key) for key in definition if isinstance(definition.get(key), (bool, str, list))},
        "checks": checks,
    }


def _domain(texts, trees, accepted, archive, receipt_status, gates):
    c6 = trees.get(C6_SOURCE, ast.Module(body=[], type_ignores=[]))
    long_tree = trees.get(LONG_SOURCE, ast.Module(body=[], type_ignores=[]))
    ref_tree = trees.get(REF_SOURCE, ast.Module(body=[], type_ignores=[]))
    profile_tree = trees.get(LONG_PROFILE_SOURCE, ast.Module(body=[], type_ignores=[]))
    moments_tree = trees.get(MOMENTS_SOURCE, ast.Module(body=[], type_ignores=[]))
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_add = _fn(c6, "add_chart", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_phi = _fn(c6, "phi_from_Utheta", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_sigma = _fn(c6, "sigma")
    c6_sigma_prime = _fn(c6, "sigma_prime")
    long_eval = _fn(long_tree, "evaluate", "ActualLongReshapeEndpoint")
    ref_packet = _fn(ref_tree, "packet", "ActualReferenceRestoreFunctions")
    profile_inputs = _fn(profile_tree, "inputs", "CompliantLongReshapeProfiles")
    profile_eval = _fn(profile_tree, "evaluate", "CompliantLongReshapeProfiles")
    ref_profile_tree = trees.get(REF_PROFILE_SOURCE, ast.Module(body=[], type_ignores=[]))
    geometry = archive.get("geometry_certificate", {})
    definition = archive.get("definition", {})
    b_formula = _assign(
        c6_init,
        "self.B",
        "a.add(a.scale(a.trim(bridge_owner.core7['G']),g.neg(self.core.parameters['Lambda'])),a.log(incoming['actual']['phi'],'positive actual R110 phi'),a.fixed(g.mul(c('1/2'),g.unary('log',c(220)))),logq)",
    )
    e_formula = _assign(c6_init, "self.E", "a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))")
    normalization = _has(c6_phi, "returna.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R))", "F0isZ-only")
    long_source = texts.get(LONG_SOURCE, "")
    c6_source = texts.get(C6_SOURCE, "")
    ref_source = texts.get(REF_SOURCE, "")
    profile_source = texts.get(LONG_PROFILE_SOURCE, "")
    ref_profile_source = texts.get(REF_PROFILE_SOURCE, "")
    bridge_c6_source = texts.get(BRIDGE_C6_SOURCE, "")
    core_inlet_source = texts.get(CORE_INLET_SOURCE, "")
    moments_source = texts.get(MOMENTS_SOURCE, "")
    accepted_domain_evidence = accepted.get("component", {}).get("domain_evidence", {})
    source_shape_map = _source_has(
        long_source,
        "inlets=dict(theta=f.scale(f.multiply(histories['H'],invphi),c.mpf('.5'))",
        "theta_z=f.scale(f.multiply(histories['K'],invphi),c.mpf('.5'))",
        "pressure=f.multiply(histories['C'],invphi2)",
        "swirl=f.multiply(histories['B'],invphi2)",
        "mean=histories['M']",
        "axial=histories['A']",
    ) and _source_has(
        long_source,
        "rates=dict(theta=c.mpf('1.6'),theta_z=c.mpf('1.6'),pressure=c.mpf('.2'),swirl=c.mpf('1.2'),mean=1,axial=1)",
    )
    raw_rates = _assign(moments_tree, "RATES", "dict(H=2,M=1,K=2,A=1,B=2,C=1)") and _source_has(
        moments_source, "RATES=dict(H=2,M=1,K=2,A=1,B=2,C=1)"
    )
    rate_conjugation = raw_rates and source_shape_map and _source_has(
        long_source,
        "terminal_Utheta_squared_derivative_ratios",
        "terminal_Utheta_derivative_ratios",
        "original_nonzero_incoming_decay_rows",
    )
    full_phi_normalization = _has(
        profile_inputs,
        "normalized=dict(theta=moments[MTH]/(phi*2),theta_z=moments[MTHZ]/(phi*2),pressure=moments[MP]/square(phi),swirl=moments[MZT]['swirl']/square(phi),mean=moments[MZ],axial=moments[MZT]['axial'])",
    ) and _has(
        profile_eval,
        "logu=B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)",
        "amplitude2=relative_amplitude_jet(logu,2)",
        "true_shapes=dict(theta=moment_shapes['theta']*amplitude,theta_z=moment_shapes['theta_z']*amplitude,pressure=moment_shapes['pressure']*amplitude2,swirl=moment_shapes['swirl']*amplitude2,mean=mean,axial=axial)",
    )
    cutoff_conjugation = _has(
        profile_eval,
        "kernels={name:backward_kernel(c,B,self.T,y,m,low,high,self.proofs)forname,m,low,highin(('theta',1,'1.55','1.65'),('pressure',2,'.1','.3'),('swirl',2,'1.1','1.3'))}",
        "decays={name:normalized_amplitude_decay(c,B*(sig*m)-y*k,low,high,y,self.proofs)forname,m,k,low,highin(('theta',1,'1.6','1.55','1.65'),('pressure',2,'.2','.1','.3'),('swirl',2,'1.2','1.1','1.3'))}",
        "inherited={name:normalized_amplitude_decay(c,B*(sig*m)-y*k,low,high,y,self.proofs,initial=inp['moments'][name])forname,m,k,low,highin(('theta',1,'1.6','1.55','1.65'),('theta_z',1,'1.6','1.55','1.65'),('pressure',2,'.2','.1','.3'),('swirl',2,'1.2','1.1','1.3'))}",
    )
    physical_product_units = _source_has(
        profile_source,
        "moment_prefactors=dict(theta='sqrt(2)*R^(3/2)*Utheta',theta_z='sqrt(2)*R^(3/2)*Utheta',mean='R',axial='R',swirl='R*Utheta^2/2',pressure='Utheta^2/2')",
        "exact_physical_moment_definitions=dict(Mtheta='sqrt(2)*R^(3/2)*Utheta*theta_shape',Mtheta_z='sqrt(2)*R^(3/2)*Utheta*theta_z_shape',Mz='R*mean',Mztheta='R*axial-R*Utheta^2*swirl_shape/2',Mp='Utheta^2*pressure_shape/2',P='Pstar^2*original_P0_normalized+Mp',Ur='sqrt(R/2)*Q')",
    ) and _has(c6_init, "self.F0=F0=a.trim(bridge_owner.core7['F0'])", "phi=self.phi_from_Utheta(Utheta,R)")
    # The C6 generic source is the authoritative physical-unit map.  Keep the
    # actual W and FF constructions in the witness instead of accepting the
    # descriptive metadata emitted by the C2 profiles alone.
    c6_generic_products = _source_has(
        bridge_c6_source,
        "inverse=g.quotient(g.one,Pstar,'same strictly positive Pstar')",
        "W=a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))",
        "FF=a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))",
        "h=a.scale(a.mul(W,J['H']),'1/2')",
        "k=a.scale(a.mul(W,J['K']),g.mul(c('1/2'),inverse))",
        "e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B']))",
        "p=a.mul(FF,J['C'])",
        "fields['pressure']=a.add(fields['P0'],fields['p'])",
    )
    core_generic_products = _source_has(
        core_inlet_source,
        "self.Pstar=g.unary('exp',self.core.parameters['logP'])",
        "self.F0=self.core_functions['F0']",
        "FF=a.scale(a.mul(F,F),g.mul(self.radius,invS2))",
        "self.fields['pressure']=a.add(self.P0,self.fields['p'])",
    )
    c2_long_products = _source_has(
        long_source,
        "ratio=IntervalTaylor(c,[c.mpf(0)]+list(logU.coefficients[1:])).exp()",
        "ratio2=IntervalTaylor(c,[c.mpf(0)]+[2*v for v in logU.coefficients[1:]]).exp()",
        "Mtheta=f.scale(dressed['theta'],R*rroot*U*c.sqrt(2))",
        "Mtheta_z=f.scale(dressed['theta_z'],R*rroot*U*c.sqrt(2))",
        "Mztheta=f.add(f.scale(shapes['axial'],R),f.scale(dressed['swirl'],-R*U2*c.mpf('.5')))",
        "Mp=increment",
        "physical_cumulative_moment_axial5=moment",
    )
    c2_reference_products = _source_has(
        ref_source,
        "ratio=IntervalTaylor(c,[c.mpf(0)]+list(logu.coefficients[1:])).exp();ratio2=IntervalTaylor(c,[c.mpf(0)]+[2*v for v in logu.coefficients[1:]]).exp()",
        "Mtheta=[f.scale(dress(binrate(sy['theta'],'1.6',k),1),R*rootR*U*c.sqrt(2)) for k in range(5)]",
        "Mztheta=[f.add(f.scale(binrate(sy['axial'],1,k),R),f.scale(dress(binrate(sy['swirl'],'1.2',k),2),-R*U2*c.mpf('.5'))) for k in range(5)]",
        "Mp=increments",
        "physical_velocity_y_Z_mixed4",
        "physical_pressure_y_Z_mixed4",
        "physical_five_primitive_y_Z_mixed4",
    )
    c2_physical_product_sources = c6_generic_products and core_generic_products and c2_long_products and c2_reference_products
    c6_long_logu = _source_has(
        c6_source,
        "logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q))),a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
    )
    c6_reference_logu = _source_has(
        c6_source,
        "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
    )
    c6_restore_logu = _source_has(
        c6_source,
        "logUtheta=a.add(a.neg(logq),a.fixed(g.add(self.core.parameters['logP'],g.mul(c('1/10'),offset))))",
        "physical_Utheta=Utheta,log_Utheta=logUtheta",
    )
    c2_long_logu = _source_has(
        profile_source,
        "logu=B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)",
    )
    c2_reference_logu = _source_has(
        ref_profile_source,
        "logu=-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)",
        "logu=-logarithm(1+square(inp['z']))+(-8+t)/10",
        "logu=-logarithm(1+square(inp['z']))+offset/10",
    )
    pstar_recipe = _source_has(
        c6_source,
        "self.Pstar=Pstar=bridge_owner.inlet.Pstar",
        "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
    ) and _source_has(
        ref_source,
        "self.logP=f.logs[1]/2",
        "self.logref=10*(self.logC+self.logP)",
        "U=self.scalar_exp(logu[0]+self.logP)",
    )
    per_chart_physical_Utheta_source_equality = (
        c6_long_logu and c6_reference_logu and c6_restore_logu
        and c2_long_logu and c2_reference_logu and pstar_recipe
        and _source_has(c6_source, "def phi_from_Utheta", "physical_Utheta=Utheta", "phi=self.phi_from_Utheta(Utheta,R)")
    )
    phi_y_chart_qualification = (
        _assign(c6_init, "slope", "a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T,'same positive original T'))))")
        and _assign(c6_init, "actual", "dict(phi=phi,U=incoming['actual']['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0),physical_Utheta=Utheta,log_Utheta=logUtheta)")
        and _assign(c6_init, "actual", "dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),physical_Utheta=Utheta,log_Utheta=logUtheta)")
        and _assign(c6_init, "actual", "dict(phi=phi,U=U,phi_y=a.scale(phi,'-2/5'),U_y=a.scale(E,g.neg(sigma_prime(g,t))),physical_Utheta=Utheta,log_Utheta=logUtheta)")
        and _assign(c6_init, "actual", "dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),physical_Utheta=Utheta,log_Utheta=logUtheta)")
    )
    c6_u_alias_and_mean_source = _source_has(
        c6_source,
        "actual=dict(phi=phi,U=incoming['actual']['U']",
        "self.endpoints['R110']['actual']['U']",
    ) and _source_has(
        long_source,
        "incoming_mean=f.scale(f.add(self.inlets['mean'],f.scale(self.V,-1)),dmean)",
        "mean=f.add(self.V,incoming_mean)",
    )
    domain = {
        "finite_positive_T": geometry.get("directed_T_positive") is True and definition.get("full_original_finite_T_400A_not_shortened_or_infinite") is True and _has(c6_init, "self.T=T=self.native_parameter('T','same original T=400*A singleton; all-N long window_length exact tuple guards')"),
        "finite_positive_reference_length": geometry.get("directed_reference_gap_positive") is True and _has(c6_init, "self.reference_length=length=g.sub(g.sub(logref,T),c(8))") and _source_has(ref_source, "self.gap=self.logref-self.T-8"),
        "positive_R110_phi": geometry.get("positive_R110_phi") == "positive true core inlet Phi times finite real exponentials" and _has(c6_init, "phi=self.phi_from_Utheta(Utheta,R)", "Utheta=a.exp(logUtheta)"),
        "finite_positive_F0": accepted_domain_evidence.get("accepted_macro_positive_barphi_F0_and_ODE_domain") is True and _has(c6_init, "self.F0=F0=a.trim(bridge_owner.core7['F0'])", "phi_from_Utheta(Utheta,R)") and _source_has(core_inlet_source, "self.F0=self.core_functions['F0']"),
        "exponential_phi_source": _has(c6_init, "Utheta=a.exp(logUtheta)", "phi=self.phi_from_Utheta(Utheta,R)") and full_phi_normalization,
        "anchored_B_source": b_formula and definition.get("original_B_from_actual_R110_and_anchored_G") is True and _source_has(long_source, "G*(-Lambda)+logarithm(phi)+logarithm(1+z*z)+c.ln(2*R)/2"),
        "E_source": e_formula and definition.get("original_axial_restoration_4Z_plus_actual_E_times_one_minus_sigma") is True and _source_has(c6_source, "self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))"),
        "physical_Utheta_normalization": normalization and _source_has(c6_source, "def phi_from_Utheta", "sqrt(2R)*F0; F0 is Z-only"),
        "per_chart_physical_Utheta_source_equality": per_chart_physical_Utheta_source_equality,
        "pstar_exp_logP_recipe": pstar_recipe and _source_has(core_inlet_source, "self.Pstar=g.unary('exp',self.core.parameters['logP'])"),
        "c6_physical_generic_W_FF_source": c6_generic_products and core_generic_products,
        "c2_physical_product_source_assignments": c2_physical_product_sources,
        "smooth_sigma_and_derivative": _has(c6_sigma, "previous.inlet.previous.lower.own.jets.complete.Gamma.sigma(g,x)") and _has(c6_sigma_prime, "original_flat_sigma_prime", "same original sigma Taylor coefficient[1]", "exact_flat_source_derivative_not_a_bound=True"),
        "triangular_fields_and_histories": _has(c6_add, "previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)") and _has(long_eval, "actual_terminal_log_radius_ODE_rows=ODE", "original_nonzero_incoming_decay_rows=decays") and _has(ref_packet, "actual_centered_y_derivative_axial5=yr", "actual_normalized_shape_y_derivative_axial5=sy"),
        "six_history_rate_conjugation": rate_conjugation,
        "full_phi_history_normalization": full_phi_normalization,
        "m_half_cutoff_conjugation": cutoff_conjugation,
        "physical_product_units": physical_product_units and c2_physical_product_sources,
        "phi_y_chart_qualification": phi_y_chart_qualification,
        "c6_U_alias_and_mean_transport": c6_u_alias_and_mean_source,
        "common_Z_source_interval": _source_has(long_source, "if ep(self.Z)[0]<-1 or ep(self.Z)[1]>1") and _source_has(ref_source, "if ep(self.Z)[0]<-1 or ep(self.Z)[1]>1"),
        "ordinary_Z_product_dressing": _source_has(long_source, "terminal_Utheta_squared_derivative_ratios", "physical_cumulative_moment_axial5") and _source_has(ref_source, "physical_five_primitive_y_Z_mixed4", "radial_prefactors_differentiated_before_mixed_grid"),
        "finite_interval_chart_chain": all(isinstance(archive.get("chart_ports", {}).get(name), Mapping) for name in CHARTS),
        "accepted_source_receipts": bool(receipt_status) and all(receipt_status.values()),
        "historical_profile_receipts_algorithm_only": receipt_status.get("historical_profile_long_linkage") is True and receipt_status.get("historical_profile_reference_linkage") is True,
        "source_only_ancestor": accepted.get("accepted") is True and accepted.get("ports", {}).get("endpoint_rows_are_values") is False,
    }
    domain["proved"] = all(value is True for value in domain.values())
    domain["theorem"] = (
        "On the common finite Z source interval, inherited positive F0 and anchored B, finite-positive T=400A and "
        "Lref, the exponential phi source, smooth flat sigma and triangular six-history linear systems give the "
        "parameter-dependent source continuation and ordinary-Z differentiation through T, Lref, 1 and 1 Jacobians. "
        "The physical normalization is the source identity phi=Utheta/(F0*sqrt(2R)); actual C6 W/FF and C2 ratio/ratio2 "
        "products are differentiated before any enclosure. No finite-N, density, norm, high-radial join or temporal-recursion conclusion is made."
    )
    domain["rates"] = {
        "raw_C6": dict(RAW_HISTORY_RATES),
        "normalized_C2": dict(NORMALIZED_SHAPE_RATES),
        "conjugation": {
            "theta": {"raw": "H", "phi_power": "1", "identity": "2-2/5=8/5"},
            "theta_z": {"raw": "K", "phi_power": "1", "identity": "2-2/5=8/5"},
            "pressure": {"raw": "C", "phi_power": "2", "identity": "1-4/5=1/5"},
            "swirl": {"raw": "B", "phi_power": "2", "identity": "2-4/5=6/5"},
            "mean": {"raw": "M", "phi_power": "0", "identity": "1"},
            "axial": {"raw": "A", "phi_power": "0", "identity": "1"},
        },
        "source_map": "H_y=2phi-2H; K_y=2phi*U-2K; C_y=phi^2-C; B_y=phi^2-2B; M_y=U-M; A_y=U^2-A. Shapes are H/(2phi), K/(2phi), C/phi^2, B/phi^2, M, A respectively; C6 U aliases C2 V and incoming memories remain nonreset.",
    }
    return domain


def _witnesses(texts, trees, accepted, archive, domain, receipt_status):
    c6 = trees.get(C6_SOURCE, ast.Module(body=[], type_ignores=[]))
    bridge_c6 = trees.get(BRIDGE_C6_SOURCE, ast.Module(body=[], type_ignores=[]))
    core_inlet = trees.get(CORE_INLET_SOURCE, ast.Module(body=[], type_ignores=[]))
    long_tree = trees.get(LONG_SOURCE, ast.Module(body=[], type_ignores=[]))
    ref_tree = trees.get(REF_SOURCE, ast.Module(body=[], type_ignores=[]))
    long_profile_tree = trees.get(LONG_PROFILE_SOURCE, ast.Module(body=[], type_ignores=[]))
    ref_profile_tree = trees.get(REF_PROFILE_SOURCE, ast.Module(body=[], type_ignores=[]))
    c6_init = _fn(c6, "__init__", "CurrentOriginalPrescribedR100RmC6Functions")
    c6_add = _fn(c6, "add_chart", "CurrentOriginalPrescribedR100RmC6Functions")
    long_eval = _fn(long_tree, "evaluate", "ActualLongReshapeEndpoint")
    ref_reference = _fn(ref_tree, "reference", "ActualReferenceRestoreFunctions")
    ref_restoration = _fn(ref_tree, "restoration", "ActualReferenceRestoreFunctions")
    ref_post = _fn(ref_tree, "postrestore", "ActualReferenceRestoreFunctions")
    ref_packet = _fn(ref_tree, "packet", "ActualReferenceRestoreFunctions")
    c6_source = texts.get(C6_SOURCE, "")
    bridge_c6_source = texts.get(BRIDGE_C6_SOURCE, "")
    core_inlet_source = texts.get(CORE_INLET_SOURCE, "")
    long_source = texts.get(LONG_SOURCE, "")
    ref_source = texts.get(REF_SOURCE, "")
    long_profile_source = texts.get(LONG_PROFILE_SOURCE, "")
    ref_profile_source = texts.get(REF_PROFILE_SOURCE, "")
    accepted_ports = accepted.get("ports", {})
    exact_endpoint_nodes = archive.get("exact_endpoint_nodes", {})
    accepted_endpoint_port_identity = all(
        exact_endpoint_nodes.get(endpoint) == accepted_ports.get(port_name)
        for endpoint, port_name in (("R1", "R1_end_function_ports"), ("R2", "R2_end_function_ports"), ("R110", "R110_end_function_ports"))
    )
    return {
        "accepted_R100_switch_background_and_own_history": accepted.get("accepted") is True,
        "exact_endpoint_function_ports": archive.get("graph_ports_valid") is True and accepted.get("ports", {}).get("endpoint_rows_are_values") is False,
        "accepted_R110_R1_R2_endpoint_node_identity": accepted_endpoint_port_identity and archive.get("endpoint_rows_are_values") is False,
        "saved_chart_and_endpoint_node_lists": bool(archive.get("exact_chart_nodes")) and bool(archive.get("exact_endpoint_nodes")) and all(archive.get("incoming_history_identity", {}).values()) and all(archive.get("endpoint_substitution_identity", {}).values()),
        "accepted_endpoint_ports_are_saved_nodes": all(archive.get("endpoint_substitution_identity", {}).values()) and all(archive.get("fixed_measure_identity", {}).values()),
        "finite_T_400A": domain.get("finite_positive_T") is True,
        "anchored_B_is_nonconstant_source": domain.get("anchored_B_source") is True and _has(c6_init, "self.B=B=a.add", "a.log(incoming['actual']['phi']"),
        "E_is_U_R110_minus_4Z": domain.get("E_source") is True,
        "R110_to_long_reshape": _has(c6_init, "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])", "self.endpoints['Rsh']=self.end('long_reshape',g.one)") and domain.get("anchored_B_source") is True,
        "long_to_reference": _has(c6_init, "self.add_chart('reference',q,R,length,actual,incoming['moments'])", "self.endpoints['Rz']=self.end('reference',g.one)") and _has(ref_reference, "log_radius_Jacobian=1", "exact_source='Rsh*exp(phase*(10(logCstar+logPstar)-T-8))'") ,
        "reference_to_restoration": (
            _has(c6_init, "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])", "self.endpoints['Rrestore']=self.end('restoration',g.one)")
            and _has(ref_restoration, "self.E", "offset=-8+t", "logu=-logarithm(1+self.z*self.z)+offset/10", "logR=c.ln(110)+self.logref+offset", "exact_source='Rz*exp(t); Rz=Rref*exp(-8)'")
        ),
        "restoration_to_post_restore": _has(c6_init, "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])", "self.endpoints['Rm']=self.end('post_restore',g.one)") and _has(ref_post, "offsets[-7,-6] before active patch", "log_radius_Jacobian=1"),
        "exact_offsets_and_jacobians": domain.get("finite_interval_chart_chain") is True and _source_has(c6_source, "exact_Rref_Rz_Rrestore_Rm_offsets=[-8,-7,-6]") and _source_has(ref_source, "log_radius_Jacobian=1"),
        "physical_Utheta_F0_sqrt_2R_normalization": domain.get("physical_Utheta_normalization") is True and _has(_fn(c6, "phi_from_Utheta", "CurrentOriginalPrescribedR100RmC6Functions"), "self.F0", "sqrt", "2", "R") and _source_has(c6_source, "physical_Utheta=Utheta", "phi=self.phi_from_Utheta(Utheta,R)"),
        "physical_Utheta_source_normalization_AST": (
            _assign(_fn(c6, "__init__", "CurrentOriginalPrescribedR100RmC6Functions"), "self.Pstar", "bridge_owner.inlet.Pstar")
            and _assign(_fn(ref_tree, "__init__", "ActualReferenceRestoreFunctions"), "self.logP", "f.logs[1]/2")
            and _has(_fn(c6, "phi_from_Utheta", "CurrentOriginalPrescribedR100RmC6Functions"), "Utheta", "self.F0", "sqrt", "2", "R")
        ),
        "per_chart_physical_Utheta_Pstar_source_equality": domain.get("per_chart_physical_Utheta_source_equality") is True and _source_has(ref_source, "self.logP=f.logs[1]/2", "U=self.scalar_exp(logu[0]+self.logP)"),
        "physical_generic_W_FF_source_AST": domain.get("c6_physical_generic_W_FF_source") is True and _source_has(bridge_c6_source, "W=a.scale(F0", "FF=a.scale(a.mul(F0,F0)"),
        "core_Pstar_and_F0_source_AST": _source_has(core_inlet_source, "self.Pstar=g.unary('exp',self.core.parameters['logP'])", "self.F0=self.core_functions['F0']"),
        "C2_actual_ratio_ratio2_and_product_sources": domain.get("c2_physical_product_source_assignments") is True and _source_has(long_source, "terminal_Utheta_derivative_ratios", "terminal_Utheta_squared_derivative_ratios") and _source_has(ref_source, "physical_five_primitive_y_Z_mixed4", "physical_pressure_y_Z_mixed4"),
        "C2_profile_logu_sources": _source_has(long_profile_source, "logu=B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)") and _source_has(ref_profile_source, "logu=-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)", "logu=-logarithm(1+square(inp['z']))+(-8+t)/10", "logu=-logarithm(1+square(inp['z']))+offset/10"),
        "phi_y_chart_qualification": domain.get("phi_y_chart_qualification") is True,
        "c6_U_alias_and_mean_U_minus_mean_source": domain.get("c6_U_alias_and_mean_transport") is True,
        "full_phi_history_normalization": domain.get("full_phi_history_normalization") is True and domain.get("exponential_phi_source") is True,
        "m_half_cutoff_conjugation": domain.get("m_half_cutoff_conjugation") is True,
        "physical_product_units": domain.get("physical_product_units") is True,
        "six_history_rates_nonreset_and_conjugated": domain.get("six_history_rate_conjugation") is True and _has(c6_add, "incoming", "previous.histories") and _has(long_eval, "original_nonzero_incoming_decay_rows=decays", "actual_inherited_terminal_history_contributions=inherited"),
        "smooth_sigma_and_triangular_fields": domain.get("smooth_sigma_and_derivative") is True and domain.get("triangular_fields_and_histories") is True,
        "ordinary_Z_products_dressed_once": domain.get("ordinary_Z_product_dressing") is True and _source_has(ref_source, "radial_prefactors_differentiated_before_mixed_grid", "actual_source_coefficients_not_selected_from_caps=True"),
        "endpoint_rows_are_enclosures": accepted.get("ports", {}).get("endpoint_rows_are_values") is False and _source_has(long_source, "original_kernel_coordinate_change_is_directed_enclosure_only=True", "physical_positive_factors_not_materialized=True"),
        "historical_profile_receipts_algorithm_only": domain.get("historical_profile_receipts_algorithm_only") is True,
        "all_new_receipts_accepted": bool(receipt_status) and all(receipt_status.values()),
    }


def _result(gates, nodes, actual_digest, accepted, archive, domain, witnesses, receipt_status):
    graph_ok = isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT and actual_digest == EXPECTED_NODES_DIGEST
    witness_ok = all(value is True for value in witnesses.values())
    receipts_ok = bool(receipt_status) and all(receipt_status.values())
    background = bool(graph_ok and accepted.get("accepted") is True and archive.get("graph_ports_valid") is True and domain.get("proved") is True and receipts_ok and not gates and witness_ok)
    own = bool(background and witnesses.get("six_history_rates_nonreset_and_conjugated") is True)
    component = {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(background and own),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "full_R110_Rm_background_projection_proved": background,
        "full_R110_Rm_own_history_projection_proved": own,
        "full_R110_Rm_evolution_identity_proved": bool(background and own),
        "full_R110_Rm_own_history_identity_proved": own,
        "full_R110_Rm_projection_proved": bool(background and own),
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_orders": [0, 1, 2],
        "ordinary_Z_not_temporal_n": True,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "same_exact_width_source_proved": accepted.get("component", {}).get("same_exact_width_source_proved") is True,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "finite_N_sc_over_2": False,
        "global_five_moment_closure": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "mathematical_operator_replayed": False,
        "full_goal_complete": False,
        "ast_witnesses": witnesses,
        "accepted_receipts": {name: True for name, ok in receipt_status.items() if ok},
        "receipt_status": receipt_status,
        "historical_profile_receipts": {
            "algorithm_only": True,
            "actual_five_defect_family_sha256": HISTORICAL_PROFILE_PROVENANCE["actual_five_defect_family_sha256"],
            "implicit_source_sha256": HISTORICAL_PROFILE_PROVENANCE["implicit_source_sha256"],
            "reference_join_family_sha256": HISTORICAL_PROFILE_PROVENANCE["reference_join_family_sha256"],
            "pinned_hashes": {
                name: RECORD_HASHES[name]
                for name in (LONG_PROFILE_SOURCE, LONG_PROFILE_REPORT, LONG_PROFILE_CHECK,
                             REF_PROFILE_SOURCE, REF_PROFILE_REPORT, REF_PROFILE_CHECK)
            },
            "current_owner_links": {
                "long": receipt_status.get("historical_profile_long_linkage") is True,
                "reference": receipt_status.get("historical_profile_reference_linkage") is True,
            },
            "values_promoted": False,
        },
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_SHA256,
        "accepted_source_stage_report": SWITCH_REPORT,
        "accepted_source_stage_sha256": SWITCH_SHA256,
        "source_family": FAMILY,
        "projection_kind": "exact original R110 to finite long reshape, Rsh, Rz, restoration and Rm source-function identity; saved rows remain enclosures",
        # Public endpoint port object: these are saved FunctionJet node
        # identities, never enclosure values.  Keep it separate from the
        # legacy summary object so consumers cannot mistake shape metadata for
        # an exact port binding.
        "endpoint_ports": {
            "endpoint_rows_are_values": False,
            "accepted_R100_ports": accepted.get("ports", {}),
            "exact_chart_nodes": archive.get("exact_chart_nodes", {}),
            "exact_endpoint_nodes": archive.get("exact_endpoint_nodes", {}),
            "incoming_history_identity": archive.get("incoming_history_identity", {}),
            "endpoint_substitution_identity": archive.get("endpoint_substitution_identity", {}),
            "fixed_measure_identity": archive.get("fixed_measure_identity", {}),
            "parameter_nodes": archive.get("parameter_nodes", {}),
            "R110": archive.get("exact_endpoint_nodes", {}).get("R110", {}),
            "Rsh": archive.get("exact_endpoint_nodes", {}).get("Rsh", {}),
            "Rz": archive.get("exact_endpoint_nodes", {}).get("Rz", {}),
            "Rrestore": archive.get("exact_endpoint_nodes", {}).get("Rrestore", {}),
            "Rm": archive.get("exact_endpoint_nodes", {}).get("Rm", {}),
            "Jacobians": {"long_reshape": "T", "reference": "Lref", "restoration": 1, "post_restore": 1},
            "offsets": [-8, -7, -6],
        },
        "exact_endpoint_ports": {
            "endpoint_rows_are_values": False,
            "accepted_R100_ports": accepted.get("ports", {}),
            "C6_chart_ports": archive.get("chart_ports", {}),
            "C6_endpoint_ports": archive.get("endpoint_ports", {}),
            "C6_chart_nodes": archive.get("exact_chart_nodes", {}),
            "C6_endpoint_nodes": archive.get("exact_endpoint_nodes", {}),
            "incoming_history_identity": archive.get("incoming_history_identity", {}),
            "endpoint_substitution_identity": archive.get("endpoint_substitution_identity", {}),
            "fixed_measure_identity": archive.get("fixed_measure_identity", {}),
            "parameter_nodes": archive.get("parameter_nodes", {}),
            "R110": archive.get("endpoint_ports", {}).get("R110", {}),
            "Rsh": archive.get("endpoint_ports", {}).get("Rsh", {}),
            "Rz": archive.get("endpoint_ports", {}).get("Rz", {}),
            "Rrestore": archive.get("endpoint_ports", {}).get("Rrestore", {}),
            "Rm": archive.get("endpoint_ports", {}).get("Rm", {}),
            "Jacobians": {"long_reshape": "T", "reference": "Lref", "restoration": 1, "post_restore": 1},
            "offsets": [-8, -7, -6],
        },
        "domain_evidence": domain,
        "source_normalization_theorem": {
            "identity": "phi=Utheta/(F0*sqrt(2R))",
            "F0": "positive Z-only source factor",
            "full_phi": "log(phi)=log(Utheta)-log(F0)-1/2 log(2R); long_reshape phi_y/phi=-2/5-B*sigma_prime/T, while reference/restoration/post_restore phi_y/phi=-2/5",
            "per_chart_logUtheta_over_Pstar": {
                "long_reshape": "log(Utheta/Pstar)=B*(1-sigma(q))-log(1+Z^2)+y/10-logC-logP",
                "reference": "log(Utheta/Pstar)=-log(1+Z^2)+y/10-logC-logP",
                "restoration": "log(Utheta/Pstar)=-log(1+Z^2)+(-8+t)/10",
                "post_restore": "log(Utheta/Pstar)=-log(1+Z^2)+offset/10",
                "Pstar": "logP=f.logs[1]/2; Pstar=exp(logP)",
            },
            "physical_Utheta_source_identity": "C6 physical Utheta=Pstar*exp(log(Utheta/Pstar)); phi=Utheta/(F0*sqrt(2R))",
            "ordinary_Z": "all ordinary Z derivatives of Utheta, F0 and sqrt(2R) are retained before product formation",
            "normalized_history_shapes": {
                "theta": "H/(2phi)",
                "theta_z": "K/(2phi)",
                "pressure": "C/phi^2",
                "swirl": "B/phi^2",
                "mean": "M",
                "axial": "A",
            },
            "shape_transport": {
                "theta": "1+[-8/5+B*sigma_prime/T]*theta",
                "theta_z": "U+[ -8/5+B*sigma_prime/T]*theta_z",
                "pressure": "1+[-1/5+2B*sigma_prime/T]*pressure",
                "swirl": "1+[-6/5+2B*sigma_prime/T]*swirl",
                "mean": "U-mean",
                "axial": "U^2-axial",
            },
            "chart_qualification": {
                "long_reshape": "phi_y/phi=-2/5-B*sigma_prime/T; B*sigma_prime/T remains in the C2 kernel",
                "reference": "phi_y/phi=-2/5",
                "restoration": "phi_y/phi=-2/5",
                "post_restore": "phi_y/phi=-2/5",
            },
            "physical_products": {
                "Mtheta": "F0*R^2*H",
                "Mtheta_z": "F0*R^2*K",
                "Mp": "F0^2*R*C",
                "Mztheta": "R*A-F0^2*R^2*B",
                "Mz": "R*M",
                "pressure": "Pstar^2*P0+Mp",
            },
            "physical_product_source_witness": {
                "C6_generic": "W=F0*sqrt(2R)/Pstar; FF=F0^2*R/Pstar^2; h=W*H/2; k=W*K/(2*Pstar); e=A/Pstar^2-FF*B; p=FF*C; pressure=P0+p",
                "C2_long": "ratio and ratio2 dress the full F0/product derivatives before Mtheta, Mtheta_z, Mztheta and Mp",
                "C2_reference": "ratio and ratio2 dress the full F0/product derivatives before physical mixed y/Z rows",
            },
            "alias_and_transport": "C6 U is the C2 V source; mean transport is U-mean and axial transport is U^2-axial",
            "cutoff_factors": "m=1 for theta/theta_z and m=2 for pressure/swirl; incoming memories are convolved before enclosure",
            "finite_interval": "T=400A and Lref are finite positive inherited source lengths; restoration and post_restore Jacobians are 1",
            "analytic_domain": "On the common finite Z source interval, inherited analytic positive F0 and anchored B, exponential phi, smooth sigma, and triangular finite linear history systems provide the source-level parameter-dependent uniqueness theorem",
            "sigma": "smooth flat sigma and sigma_prime source identities",
            "scope": "source theorem only; no field values, norms, quadrature, density, finite-N or high-radial/mixed-join claim",
        },
        "history_rate_conjugation": domain.get("rates", {}),
        "retained_construction_data": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "six_channels_retained": True,
            "six_signed_channels_retained": True,
            "all_old_error_data_retained": True,
            "reserved_contribution_registered_once": True,
        },
        "scope_flags": {
            "same_density_certificate": False,
            "finite_N_history_identity": False,
            "physical_Q_y_Z1": False,
            "global_five_moment_closure": False,
            "time_order_recursion": False,
            "timeorder_recursion": False,
            "full_goal_complete": False,
        },
        "fixture_limits": {
            "finite_fixture_only": True,
            "actual_source_admission": False,
            "scope": "R110/Rsh/Rz/restoration/Rm source functions and saved enclosures; finite-N and global closure remain open",
        },
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "frozen_raw_sha256": {C6_ARCHIVE: C6_ARCHIVE_SHA256},
        "exact_scope": {
            "finite_T": "T=400*A",
            "B": "-Lambda*G_anchored+log(phi_R110)+1/2 log(220)+log(1+Z^2)",
            "E": "U_R110-4Z",
            "normalization": "phi=Utheta/(F0*sqrt(2R)); Pstar=exp(logP), logP=f.logs[1]/2",
            "physical_products": {
                "W": "F0*sqrt(2R)/Pstar",
                "FF": "F0^2*R/Pstar^2",
                "moments": ["F0*R^2*H", "F0*R^2*K", "R*A-F0^2*R^2*B", "R*M"],
                "pressure": "Pstar^2*P0+F0^2*R*C",
            },
            "offsets": [-8, -7, -6],
            "Jacobians": {"long_reshape": "T", "reference": "Lref", "restoration": 1, "post_restore": 1},
            "ordinary_Z_orders": [0, 1, 2],
            "source_only": True,
            "numeric_values_computed": False,
        },
        "source_inspection_audit": {
            "source_records_only": True,
            "static_source_AST_graph_and_saved_bound_inspection_only": True,
            "nodes_hashed_only": True,
            "source_graph_substitutions_checked": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "producer_provider_builder_checker_replayed": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "recurrence_calls": 0,
        },
        "worker_role": "luna_max_worker",
        "runtime_model_metadata": None,
        "unresolved_gates": sorted(set(gates)) if gates or not (background and own) else [],
        "potential_gates": ["same_density_certificate", "finite_N_history_identity", "physical_Q_y_Z1", "global_five_moment_closure", "time_order_recursion"],
    }
    return component


def derive(nodes, source_records, *, expected_nodes_digest):
    """Return a fail-closed, source-only R110-to-Rm projection."""
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
    accepted = _accepted_switch(texts, gates)
    receipt_status = _receipts(texts, gates)
    trees = {name: _tree(value, name, gates) for name, value in texts.items() if name.endswith(".py")}
    archive = _archive_ports(texts, nodes, gates)
    domain = _domain(texts, trees, accepted, archive, receipt_status, gates)
    witnesses = _witnesses(texts, trees, accepted, archive, domain, receipt_status)
    for name, ok in witnesses.items():
        if ok is not True:
            gates.append("ast_witness_missing:" + name)
    return _result(gates, nodes, actual_digest, accepted, archive, domain, witnesses, receipt_status)
