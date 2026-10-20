"""Static source projection for the original ``Rm`` C2/C6 patch.

The kernel consumes authenticated source snapshots and saved receipts only.  It
does not import or execute any producer, provider, evaluator, recurrence,
quadrature, norm, or project module.  A patch projection is admitted only when
the current C2 feedback source is connected to the same Rsh source functions as
the accepted R110-to-Rm source projection and the C6 source exposes the same
normalized map and physical units.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path


SCHEMA = "original_Rm_patch_C2_C6_source_projection_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
ACCEPTED_R110_REPORT = "lei_ren_part1_paper_compliant_current_original_R110_Rm_source_projection.json"
ACCEPTED_R110_SHA256 = "08ef44c6971d07a4c6583d481ef3a0ce1b37bf069995a732bdf2a9a7ccfa074e"

C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
C6_ARCHIVE_SHA256 = "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
C6_ARCHIVE_DECODED_SHA256 = "b5671cc0c60821db9a5c591ea899c9616c3cc949cec4e9214cc57a31ae6cdc10"
C6_CHECK = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check.json"
C6_CHECK_SHA256 = "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7"

C2_SOURCE = "lei_ren_part1_paper_compliant_actual_patch_mixed_C4.py"
C2_SOURCE_SHA256 = "afe05b65f64c29ac3bc47528bbea58abbd8186ced832a14bfad3adebf5b97a3f"
C2_CHECK = "lei_ren_part1_paper_compliant_actual_patch_mixed_C4_check.json"
C2_CHECK_SHA256 = "acc788b12acc83168815b959493fc62db8148548b07f4219266bf041e142172c"
C2_MOMENT_SOURCE = "lei_ren_part1_paper_compliant_actual_moment_patch.py"
C2_MOMENT_SOURCE_SHA256 = "11a1623af5f5880dd3456e455b6d7e36140b9d8cc9ff2836a417aea02da0134c"
C2_MOMENT_CHECK = "lei_ren_part1_paper_compliant_actual_moment_patch_check.json"
C2_MOMENT_CHECK_SHA256 = "8abd343db73ec66fe42254bb842b707bc7bac9cb8085fba1fefd45e429271e92"
C2_FEEDBACK_MOMENT_SOURCE = "lei_ren_part1_paper_compliant_actual_feedback_moment_patch.py"
C2_FEEDBACK_MOMENT_SOURCE_SHA256 = "ce009ad8b6b2eece7801660440c8ea25c4aa60d53a3c7b83ddfaebefa8e6c000"
C2_FEEDBACK_MOMENT_CHECK = "lei_ren_part1_paper_compliant_actual_feedback_moment_patch_check.json"
C2_FEEDBACK_MOMENT_CHECK_SHA256 = "f094974b3870cad479fb4d981d49c34089970eea76036c3d7869d06b83e24018"
C2_FEEDBACK_SOURCE = "lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4.py"
C2_FEEDBACK_SOURCE_SHA256 = "caa8c42032433902413f4fa5f38781cfc9858a2f469047bb296a952ef46baf58"
C2_FEEDBACK_CHECK = "lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4_check.json"
C2_FEEDBACK_CHECK_SHA256 = "d17fec99dfd1f6c3b219fe66d9624cbc0e8879675b0c1122d65a0644a3672ba7"

REFERENCE_MIXED_SOURCE = "lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4.py"
REFERENCE_MIXED_SOURCE_SHA256 = "281db247e59739231582182a115d39a76c5d9a50240cb6aca9ec98fc887ff1cc"
REFERENCE_MIXED_CHECK = "lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4_check.json"
REFERENCE_MIXED_CHECK_SHA256 = "d3ca3b53909e2b1b7914734cbb355de7b82779a042f824fc20975c0dfd252aa3"
REFERENCE_MIXED_REPORT = "lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4.json"
REFERENCE_MIXED_REPORT_SHA256 = "9bce2be56be99ddbe585dc9745f24f33f37eed04b31aea1cb0541d27e41c700c"
RSH_JOIN_SOURCE = "lei_ren_part1_paper_compliant_actual_Rsh_source_join.py"
RSH_JOIN_SOURCE_SHA256 = "7416c4c268da304d80a622672e96716e21de11f4c4708f317586149dad5836ac"
RSH_JOIN_CHECK = "lei_ren_part1_paper_compliant_actual_Rsh_source_join_check.json"
RSH_JOIN_CHECK_SHA256 = "d713b0191fd02ff337ca5432abbf72f7776a8de64026b91547424c7ce26fe0f6"
RSH_JOIN_REPORT = "lei_ren_part1_paper_compliant_actual_Rsh_source_join.json"
RSH_JOIN_REPORT_SHA256 = "9c81ef302737728ae452f88b845535b58cb2a7ebdbe46e6ae221fcf7c4a5f7d2"

CURRENT_REFERENCE_SOURCE = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py"
CURRENT_REFERENCE_SOURCE_SHA256 = "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7"
CURRENT_REFERENCE_REPORT = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.json"
CURRENT_REFERENCE_REPORT_SHA256 = "47f21c94f6bb08ebef01288d7b506ef5135ecd2cf5e3fba9cfd15567f6cef440"
CURRENT_REFERENCE_CHECK = "lei_ren_part1_paper_compliant_current_original_reference_restore_functions_check.json"
CURRENT_REFERENCE_CHECK_SHA256 = "94f5b78791bd571ae6759d2d53e6fdda47607b047138e9c4868e2611dff10b02"
REFERENCE_PROFILE_SOURCE = "lei_ren_part1_paper_compliant_reference_restore_profiles.py"
REFERENCE_PROFILE_SOURCE_SHA256 = "bee442990e9e52b04f37063d21e2da80ed5c35ec037941ac2af5b63860308765"
REFERENCE_PROFILE_REPORT = "lei_ren_part1_paper_compliant_reference_restore_profiles.json"
REFERENCE_PROFILE_REPORT_SHA256 = "8b21a26b79401a27ff174c68a670d08cbb10ff094ce34ac5772077449aaeb8f6"
REFERENCE_PROFILE_CHECK = "lei_ren_part1_paper_compliant_reference_restore_profiles_check.json"
REFERENCE_PROFILE_CHECK_SHA256 = "fecc91c0ba8368c681a2631db174b7f25da82ee16ceace2a744fb3f545f15087"

CURRENT_PATCH_STRESS_SOURCE = "lei_ren_part1_paper_compliant_current_patch_stress_operator.py"
CURRENT_PATCH_STRESS_SOURCE_SHA256 = "5e0bcb8cf9b2d9804bb90dfaa6a4d1de2ac78617a2c04eed5d6b2d585faa7a5b"
BACKGROUND_SOURCE = "lei_ren_part1_paper_compliant_current_actual_patch_background_tensor.py"
BACKGROUND_SOURCE_SHA256 = "6baef3f07d52757af232c57b3981e39833895cd06263b8666429251affc3062b"
BACKGROUND_ARCHIVE = "lei_ren_part1_paper_compliant_current_actual_patch_background_tensor.json.gz"
BACKGROUND_ARCHIVE_SHA256 = "1ff76a178ec0396c2fe3df99f50585e6c1107db9716101c5e108f9151543edfc"
BACKGROUND_ARCHIVE_DECODED_SHA256 = "7835d33d7268091a05b5b93edfda2018b438622a90733cdc2b369207b593fc01"
BACKGROUND_CHECK = "lei_ren_part1_paper_compliant_current_actual_patch_background_tensor_check.json"
BACKGROUND_CHECK_SHA256 = "4752f43ebda531d45033a174adcaa7ddcde79dff6ab526fac8df62a8662f9d88"

FINITE_N_MEMORY_REPORT = "lei_ren_part1_paper_compliant_current_original_finite_N_nonreset_memory_source_binding.json"
FINITE_N_MEMORY_REPORT_SHA256 = "41defb67ef05207a079075328f3f7afd3980193c4cec7adb9574a27d41f33335"

FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
}
CENTERS = ("5/4", "3/2", "7/4")
RATES = {
    "mean_error": "1",
    "mixed_error": "8/5",
    "angular_error": "8/5",
    "axial_square": "1",
    "swirl_error": "6/5",
    "pressure_error": "1/5",
}

RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    ACCEPTED_R110_REPORT: ACCEPTED_R110_SHA256,
    C6_SOURCE: C6_SOURCE_SHA256,
    C6_ARCHIVE: C6_ARCHIVE_SHA256,
    C6_CHECK: C6_CHECK_SHA256,
    C2_SOURCE: C2_SOURCE_SHA256,
    C2_CHECK: C2_CHECK_SHA256,
    C2_MOMENT_SOURCE: C2_MOMENT_SOURCE_SHA256,
    C2_MOMENT_CHECK: C2_MOMENT_CHECK_SHA256,
    C2_FEEDBACK_MOMENT_SOURCE: C2_FEEDBACK_MOMENT_SOURCE_SHA256,
    C2_FEEDBACK_MOMENT_CHECK: C2_FEEDBACK_MOMENT_CHECK_SHA256,
    C2_FEEDBACK_SOURCE: C2_FEEDBACK_SOURCE_SHA256,
    C2_FEEDBACK_CHECK: C2_FEEDBACK_CHECK_SHA256,
    REFERENCE_MIXED_SOURCE: REFERENCE_MIXED_SOURCE_SHA256,
    REFERENCE_MIXED_CHECK: REFERENCE_MIXED_CHECK_SHA256,
    REFERENCE_MIXED_REPORT: REFERENCE_MIXED_REPORT_SHA256,
    RSH_JOIN_SOURCE: RSH_JOIN_SOURCE_SHA256,
    RSH_JOIN_CHECK: RSH_JOIN_CHECK_SHA256,
    RSH_JOIN_REPORT: RSH_JOIN_REPORT_SHA256,
    CURRENT_REFERENCE_SOURCE: CURRENT_REFERENCE_SOURCE_SHA256,
    CURRENT_REFERENCE_REPORT: CURRENT_REFERENCE_REPORT_SHA256,
    CURRENT_REFERENCE_CHECK: CURRENT_REFERENCE_CHECK_SHA256,
    REFERENCE_PROFILE_SOURCE: REFERENCE_PROFILE_SOURCE_SHA256,
    REFERENCE_PROFILE_REPORT: REFERENCE_PROFILE_REPORT_SHA256,
    REFERENCE_PROFILE_CHECK: REFERENCE_PROFILE_CHECK_SHA256,
    CURRENT_PATCH_STRESS_SOURCE: CURRENT_PATCH_STRESS_SOURCE_SHA256,
    BACKGROUND_SOURCE: BACKGROUND_SOURCE_SHA256,
    BACKGROUND_ARCHIVE: BACKGROUND_ARCHIVE_SHA256,
    BACKGROUND_CHECK: BACKGROUND_CHECK_SHA256,
    FINITE_N_MEMORY_REPORT: FINITE_N_MEMORY_REPORT_SHA256,
}
DECODED_RECORD_HASHES = {
    C6_ARCHIVE: C6_ARCHIVE_DECODED_SHA256,
    BACKGROUND_ARCHIVE: BACKGROUND_ARCHIVE_DECODED_SHA256,
}
GZ_RECORDS = frozenset(DECODED_RECORD_HASHES)


def _digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()


def _path(value):
    return "experiments/root_st073/" + value


def _norm(value):
    return "".join(str(value).split())


def _get(value, dotted, default=None):
    for part in dotted.split(".") if dotted else ():
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _record(records, name, gates):
    row = records.get(name)
    if not isinstance(row, Mapping):
        gates.append("missing_source_record:" + name)
        return ""
    if row.get("path") != _path(name):
        gates.append("source_path_mismatch:" + name)
    if row.get("sha256") != RECORD_HASHES[name]:
        gates.append("source_declared_hash_mismatch:" + name)
    source = row.get("source_text")
    if not isinstance(source, str):
        gates.append("source_text_missing:" + name)
        return ""
    if name in GZ_RECORDS:
        if row.get("encoding") != "gzip":
            gates.append("encoded_archive_contract_mismatch:" + name)
            return ""
        if row.get("decoded_sha256") != DECODED_RECORD_HASHES[name]:
            gates.append("decoded_archive_hash_declaration_mismatch:" + name)
            return ""
        try:
            # Authenticate raw gzip bytes before allowing decompression.
            raw = base64.b64decode(row.get("raw_bytes_base64", ""), validate=True)
        except Exception as error:
            gates.append("archive_raw_base64_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
            gates.append("archive_raw_hash_mismatch:" + name)
            return ""
        try:
            decoded = gzip.decompress(raw)
        except Exception as error:
            gates.append("archive_decode_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(decoded).hexdigest() != DECODED_RECORD_HASHES[name]:
            gates.append("archive_decoded_hash_mismatch:" + name)
        if decoded != source.encode("utf-8"):
            gates.append("archive_source_text_mismatch:" + name)
        return source
    if row.get("encoding") is not None or row.get("raw_bytes_base64") is not None:
        gates.append("unexpected_encoded_record:" + name)
    if hashlib.sha256(source.encode("utf-8")).hexdigest() != RECORD_HASHES[name]:
        gates.append("source_hash_mismatch:" + name)
    return source


def _json(source, name, gates):
    try:
        return json.loads(source)
    except Exception as error:
        gates.append("json_parse_error:" + name + ":" + type(error).__name__)
        return {}


def _tree(source, name, gates):
    try:
        return ast.parse(source, filename=name)
    except Exception as error:
        gates.append("ast_parse_error:" + name + ":" + type(error).__name__)
        return ast.Module(body=[], type_ignores=[])


def _fn(tree, name, cls=None):
    if cls is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                return next((child for child in node.body
                             if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                             and child.name == name), None)
        return None
    return next((node for node in ast.walk(tree)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and node.name == name), None)


def _has(node, *needles):
    if node is None:
        return False
    try:
        compact = _norm(ast.unparse(node))
    except Exception:
        return False
    return all(_norm(needle) in compact for needle in needles)


def _text_has(source, *needles):
    compact = _norm(source)
    return all(_norm(needle) in compact for needle in needles)


def _assign(node, target, expression):
    if node is None:
        return False
    try:
        wanted = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
    except Exception:
        return False
    target = _norm(target)
    for item in ast.walk(node):
        if isinstance(item, ast.Assign):
            targets = item.targets
        elif isinstance(item, ast.AnnAssign):
            targets = [item.target]
        else:
            continue
        if any(_norm(ast.unparse(candidate)) == target for candidate in targets):
            if ast.dump(item.value, include_attributes=False) == wanted:
                return True
    return False


def _all_true(row, paths):
    return all(_get(row, path) is True for path in paths)


def _saved_prerequisites(texts, gates):
    accepted = _json(texts[ACCEPTED_R110_REPORT], ACCEPTED_R110_REPORT, gates)
    component = accepted.get("component_result", {})
    scope = accepted.get("scope", {})
    audit = accepted.get("source_inspection_audit", {})
    required_component = (
        "passed", "source_only", "full_R110_Rm_background_projection_proved",
        "full_R110_Rm_own_history_projection_proved", "ordinary_Z_not_temporal_n",
    )
    component_ok = _all_true(component, required_component)
    if not component_ok:
        gates.append("accepted_R110_component_not_proved")
    if component.get("nodes_digest") != EXPECTED_NODES_DIGEST:
        gates.append("accepted_R110_graph_digest_changed")
    if component.get("ordinary_Z_orders") != [0, 1, 2]:
        gates.append("accepted_R110_ordinary_orders_changed")
    if component.get("numeric_values_computed") is not False or component.get("old_operators_replayed") is not False:
        gates.append("accepted_R110_numeric_scope_changed")
    scope_ok = _all_true(scope, (
        "full_R110_Rm_background_projection_proved",
        "full_R110_Rm_own_history_projection_proved",
        "ordinary_Z_not_temporal_n",
    )) and scope.get("ordinary_Z_orders") == [0, 1, 2]
    if not scope_ok:
        gates.append("accepted_R110_scope_changed")
    audit_ok = _all_true(audit, (
        "source_records_only", "static_source_AST_graph_and_saved_bound_inspection_only",
        "nodes_hashed_only", "source_graph_substitutions_checked",
    )) and audit.get("new_source_consumer_calls") == 1 and audit.get("accepted_ancestor_operator_calls") == 0
    if not audit_ok:
        gates.append("accepted_R110_audit_changed")
    return {
        "accepted_R110_report": accepted,
        "component": component,
        "scope": scope,
        "audit": audit,
        "accepted": bool(component_ok and scope_ok and audit_ok),
        "background": component.get("full_R110_Rm_background_projection_proved") is True,
        "own_history": component.get("full_R110_Rm_own_history_projection_proved") is True,
    }


def _receipt_checks(texts, gates):
    c6 = _json(texts[C6_CHECK], C6_CHECK, gates)
    c2 = _json(texts[C2_CHECK], C2_CHECK, gates)
    moment = _json(texts[C2_MOMENT_CHECK], C2_MOMENT_CHECK, gates)
    feedback_moment = _json(texts[C2_FEEDBACK_MOMENT_CHECK], C2_FEEDBACK_MOMENT_CHECK, gates)
    feedback = _json(texts[C2_FEEDBACK_CHECK], C2_FEEDBACK_CHECK, gates)
    reference = _json(texts[REFERENCE_MIXED_CHECK], REFERENCE_MIXED_CHECK, gates)
    rsh = _json(texts[RSH_JOIN_CHECK], RSH_JOIN_CHECK, gates)
    background = _json(texts[BACKGROUND_CHECK], BACKGROUND_CHECK, gates)
    profile = _json(texts[REFERENCE_PROFILE_CHECK], REFERENCE_PROFILE_CHECK, gates)
    current_reference = _json(texts[CURRENT_REFERENCE_CHECK], CURRENT_REFERENCE_CHECK, gates)
    c6_ok = bool(
        c6.get("all_passed") is True
        and c6.get("current_original_true_Rm_leading_implicit_control_C6_functions_installed") is True
        and c6.get("current_original_true_Rm_patch_background_C6_generic_C5_functions_installed") is True
        and "current_original_true_Rm_patch_C6_functions_installed" not in c6
        and _all_true(c6.get("checks", {}), (
            "actual_source.passed", "exact_weights.passed", "true_inlet_normalization.passed",
            "normalized_conjugacy.passed", "branch_theorem.passed", "patch_matching.passed",
        ))
    )
    if not c6_ok:
        gates.append("C6_two_actual_gates_or_source_checks_missing")
    c2_ok = bool(c2.get("all_passed") is True and c2.get("actual_patch_all_mixed_derivatives_total_order_le4_available") is True)
    moment_ok = bool(moment.get("all_passed") is True and _all_true(moment, (
        "original_restoration_integrals_directly_source_bound",
        "exact_actual_Rm_inlet_and_original_P0_bound",
        "actual_five_functional_terminal_identities_connected",
    )))
    feedback_moment_ok = bool(feedback_moment.get("all_passed") is True and _all_true(feedback_moment, (
        "current_Rm_source_parent_and_P0_retained", "current_actual_moment_patch_installed",
        "current_actual_patch_implicit_axial5_recomputed", "current_five_functional_terminal_identities_connected",
    )))
    feedback_ok = bool(feedback.get("all_passed") is True and _all_true(feedback, (
        "current_Rm_source_transport_and_P0_retained", "current_actual_moment_patch_installed",
        "current_actual_patch_mixed4_available", "current_Rm_and_Rh_functional_mixed4_joins_certified",
        "current_five_functional_terminal_identities_connected", "source_joins_not_proved_by_interval_overlap",
    )))
    reference_ok = bool(reference.get("all_passed") is True and _all_true(reference, (
        "current_E_source_C2_theorem_verified_before_intersection",
        "current_Rsh_parent_and_six_centered_histories_retained",
        "current_original_restore_end_exact_4Z", "current_correlated_E_installed",
        "current_actual_reference_restore_mixed4_available",
    )))
    rsh_ok = bool(rsh.get("all_passed") is True and _all_true(rsh, (
        "current_Rsh_parent_P0_V_E_and_physical_normalization_verified",
        "source_equality_proof_not_enclosure_overlap", "current_Rsh_source_functional_join_certified",
    )))
    background_ok = bool(background.get("all_passed") is True and _all_true(background, (
        "current_actual_five_moment_patch_full_tensor_available",
        "current_actual_five_moment_patch_full_meridional_decomposition_available",
        "current_actual_patch_Rh_completed_tensor_join_certified",
        "current_actual_patch_six_internal_support_tensor_traces_certified",
        "current_implicit_axial5_partial_integrals_full_energy_cross_terms_and_same_P0_retained",
        "exact_functional_flat_limits_not_rounded_zero_boxes_used_for_six_internal_traces",
    )))
    # The profile receipt is an algorithm-only historical prerequisite.  Its
    # check schema intentionally has no native-installation boolean; requiring
    # one here would turn a saved, accepted profile receipt into a false gate.
    profile_ok = bool(profile.get("all_passed", True) is True)
    current_reference_ok = bool(current_reference.get("original_actual_two_frame_Rsh_reference_restoration_Rm_functions_installed") is True)
    checks = {
        "c6": c6_ok, "c2": c2_ok, "c2_moment": moment_ok,
        "c2_feedback_moment": feedback_moment_ok, "c2_feedback_mixed": feedback_ok,
        "reference_mixed": reference_ok, "rsh_functional_join": rsh_ok,
        "background_tensor": background_ok, "reference_profiles": profile_ok,
        "current_reference_functions": current_reference_ok,
    }
    for name, ok in checks.items():
        if not ok:
            gates.append("receipt_not_accepted:" + name)
    return {
        "checks": checks, "c6": c6, "c2": c2, "c2_moment": moment,
        "c2_feedback_moment": feedback_moment, "c2_feedback_mixed": feedback,
        "reference_mixed": reference, "rsh_functional_join": rsh,
        "background_tensor": background, "reference_profiles": profile,
        "current_reference_functions": current_reference,
    }


def _port_handles(value, node_count, expected_length=None):
    """Check saved graph handles without interpreting the graph expressions."""
    if type(value) is int:
        return 0 <= value < node_count
    if isinstance(value, list):
        if expected_length is not None and len(value) != expected_length:
            return False
        return all(_port_handles(item, node_count) for item in value)
    if isinstance(value, Mapping):
        return all(_port_handles(item, node_count) for item in value.values())
    return False


def _port_rows(block, names, node_count, row_length):
    return (
        isinstance(block, Mapping)
        and set(block) == set(names)
        and all(isinstance(block[name], list) and len(block[name]) == row_length
                and all(type(ref) is int and 0 <= ref < node_count for ref in block[name])
                for name in names)
    )


def _handles(value):
    if type(value) is int:
        return [value]
    if isinstance(value, list):
        return [ref for item in value for ref in _handles(item)]
    if isinstance(value, Mapping):
        return [ref for item in value.values() for ref in _handles(item)]
    return []


def _substitution(graph, ref, expression, variable, value):
    if type(ref) is not int or not 0 <= ref < len(graph):
        return False
    node = graph[ref]
    return isinstance(node, Mapping) and node.get("operation") == "function_substitution" and all(
        node.get(key) == expected for key, expected in (
            ("expression", expression), ("variable", variable), ("value", value),
            ("Z_independent_substitution", True)))


def _c6_archive_ports(texts, nodes, saved, gates):
    """Authenticate the C6 archive handles and bind its inlet to accepted Rm."""
    archive = _json(texts[C6_ARCHIVE], C6_ARCHIVE, gates)
    bindings = archive.get("true_Rm_patch_C6_bindings", {})
    graph_nodes = archive.get("exact_graph_nodes", [])
    node_count = len(graph_nodes) if isinstance(graph_nodes, list) else -1
    inlet = bindings.get("inlet", {})
    endpoints = bindings.get("endpoints", {})
    branch = bindings.get("branch_transfer", {})
    correspondence = _get(branch, "exact_source_correspondence", {})
    conjugacy = _get(branch, "raw_to_original_normalized_conjugacy", {})
    definition = bindings.get("definition", {})
    expected_actual = ("phi", "U", "phi_y", "U_y")
    expected_moments = ("H", "M", "K", "A", "B", "C")
    expected_fields = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
    expected_inlet_generic = (
        "Q", "transport", "theta_linear", "theta_quadratic",
        "axial_linear", "axial_quadratic", "p2",
    )
    expected_endpoint_generic = expected_inlet_generic + ("inertial_theta", "inertial_axial")
    expected_radial = ("E_y", "V_y")

    checks = {
        "archive_graph_prefix_length": node_count == 196042,
        "graph_prefix_matches_archive": (
            isinstance(nodes, list) and isinstance(graph_nodes, list)
            and len(nodes) >= node_count and nodes[:node_count] == graph_nodes
        ),
        "archive_family_matches_Rm_family": (
            _get(archive, "source_family.actual_five_defect_family_sha256")
            == FAMILY["actual_five_defect_family_sha256"]
            and _get(archive, "source_family.implicit_source_sha256")
            == FAMILY["implicit_source_sha256"]
        ),
        "inlet_actual_ports": _port_rows(inlet.get("actual"), expected_actual, node_count, 7),
        "inlet_six_history_ports": _port_rows(inlet.get("moments"), expected_moments, node_count, 7),
        "inlet_field_ports": _port_rows(inlet.get("fields"), expected_fields, node_count, 7),
        "inlet_generic_ports": _port_rows(inlet.get("generic"), expected_inlet_generic, node_count, 6),
        "inlet_radius_port": _port_handles(inlet.get("radius"), node_count),
        "Rm_field_ports": _port_rows(endpoints.get("Rm", {}).get("fields"), expected_fields, node_count, 7),
        "Rm_radial_ports": _port_rows(endpoints.get("Rm", {}).get("radial"), expected_radial, node_count, 7),
        "Rm_generic_ports": _port_rows(endpoints.get("Rm", {}).get("generic"), expected_endpoint_generic, node_count, 6),
        "Rm_radius_port": _port_handles(endpoints.get("Rm", {}).get("radius"), node_count),
        "Rh_field_ports": _port_rows(endpoints.get("Rh", {}).get("fields"), expected_fields, node_count, 7),
        "Rh_radial_ports": _port_rows(endpoints.get("Rh", {}).get("radial"), expected_radial, node_count, 7),
        "Rh_generic_ports": _port_rows(endpoints.get("Rh", {}).get("generic"), expected_endpoint_generic, node_count, 6),
        "Rh_radius_port": _port_handles(endpoints.get("Rh", {}).get("radius"), node_count),
        "definition_is_true_Rm_source": all(definition.get(key) is True for key in (
            "true_prescribed_Rm_C6_input_not_old_axial5_enclosure_adapter",
            "exact_original_beta_weights_normalization_centers_radius_and_support",
            "original_unique_common_box_branch_same_full_Jacobian_through_Z6",
            "all_six_implicit_derivatives_include_invAm2_Z_products",
            "original_five_centered_defects_from_true_m_h_k_e_p_without_P0",
            "leading_control_not_finite_N_correction_and_no_history_reset",
            "actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6",
            "original_x1_inlet_and_Rh_exp1_exit_and_log_Jacobian_1_over_x",
            "original_branch_zero_residual_is_function_identity_not_zero_containment",
        )),
        "branch_transfer_source_identity": (
            branch.get("actual_source_Z_domain") == ["-1", "1"]
            and branch.get("actual_leading_N0_equal_by_exact_integer") is True
            and branch.get("exact_defects_from_original_physical_primitives") is True
            and branch.get("original_current_profile_recipes_bound_by_checked_C6_parent") is True
            and branch.get("original_raw_beta_definition_and_finite_support_normalization_unchanged") is True
            and branch.get("original_four_cell_C0_and_Jacobian_contractions_strictly_lt1") is True
        ),
        "branch_inlet_normalized_source": (
            correspondence.get("normalized_raw_histories") == [
                "M", "H/(2*phi)", "K/(2*phi)", "A", "B/phi^2", "C/phi^2",
            ]
            and correspondence.get("actual_Rm_amplitude_identity")
            == "sqrt(2*Rm)*F0*phi/Pstar=exp(-3/5)/(1+Z^2)"
            and correspondence.get("acceptance_requires_actual_node_normalization_and_exact_weight_recipe_checks") is True
            and conjugacy.get("exact_identity") == "G(u)=S^-1 F(Su); J_G=S^-1 J_F S"
            and conjugacy.get("original_block_diagonal_L_commutes_with_S") is True
        ),
        "branch_common_box_unique": isinstance(branch.get("exact_Jacobian_invertibility_proof"), str) and
        "same exact l,q,defects,amplitude" in branch.get("exact_Jacobian_invertibility_proof", "").lower() and
        "Neumann" in branch.get("exact_Jacobian_invertibility_proof", ""),
    }
    component = saved.get("accepted_R110_report", {}).get("component_result", {})
    accepted_ports = component.get("exact_endpoint_ports", {})
    accepted_Rm = _get(accepted_ports, "C6_endpoint_nodes.Rm", {})
    accepted_normal = component.get("source_normalization_theorem", {})
    accepted_shapes = accepted_normal.get("normalized_history_shapes", {})
    cross_normalization = (
        accepted_normal.get("identity") == "phi=Utheta/(F0*sqrt(2R))"
        and accepted_normal.get("alias_and_transport")
        == "C6 U is the C2 V source; mean transport is U-mean and axial transport is U^2-axial"
        and accepted_normal.get("physical_Utheta_source_identity")
        == "C6 physical Utheta=Pstar*exp(log(Utheta/Pstar)); phi=Utheta/(F0*sqrt(2R))"
        and accepted_shapes == {
            "theta": "H/(2phi)", "theta_z": "K/(2phi)", "pressure": "C/phi^2",
            "swirl": "B/phi^2", "mean": "M", "axial": "A",
        }
        and _get(accepted_normal, "physical_product_source_witness.C6_generic")
        == "W=F0*sqrt(2R)/Pstar; FF=F0^2*R/Pstar^2; h=W*H/2; k=W*K/(2*Pstar); e=A/Pstar^2-FF*B; p=FF*C; pressure=P0+p"
    )
    checks["accepted_Rm_inlet_endpoint_identity"] = isinstance(accepted_Rm, Mapping) and accepted_Rm == inlet
    checks["accepted_Rm_rows_are_saved_handles"] = accepted_ports.get("endpoint_rows_are_values") is False
    checks["accepted_source_normalization_bridge"] = cross_normalization
    checks["accepted_endpoint_source_boundary"] = (
        _get(saved.get("accepted_R110_report", {}), "endpoint_source_boundary.actual_result_endpoint_ports") == component.get("endpoint_ports")
        and isinstance(_get(saved.get("accepted_R110_report", {}), "endpoint_source_boundary.exact_endpoint_function_ports"), str)
        and _get(saved.get("accepted_R110_report", {}), "endpoint_source_boundary.endpoint_rows_are_values") is False
        and _get(saved.get("accepted_R110_report", {}), "endpoint_source_boundary.source_overlap_is_function_identity") is False
        and _get(saved.get("accepted_R110_report", {}), "endpoint_source_boundary.interval_enclosure_is_function_value") is False
    )
    chart = bindings.get("chart", {})
    system_ref = bindings.get("system")
    system = graph_nodes[system_ref] if type(system_ref) is int and 0 <= system_ref < node_count else {}
    checks["implicit_system_source_binding"] = (
        system_ref == 193470
        and system.get("operation") == "original_Rm_unique_leading_implicit_system"
        and system.get("defining_map") == "L*h+Q(h,invAm2)+true_Rm_defects=0"
        and system.get("exact_defect_value_nodes") == [row[0] for row in bindings.get("defects", [])]
        and system.get("invAm2") == bindings.get("alpha", [None])[0]
        and system.get("exact_Z_independent_row_and_control_units") == bindings.get("units")
        and system.get("original_whole_Z_base_unique_branch_and_Jacobian_admitted") is True
        and system.get("original_leading_fixed_point_not_Rc_finite_N_control") is True
        and system.get("numerical_enclosures_not_values") is True
        and system.get("branch_receipt_sha256") == "29781f93710407553ca958a192b326f7c9bae6ebf4f269794116c452837b0702"
        and system.get("source_projection_binding") == {
            "module": "lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse.py",
            "module_sha256": "b0044165763821001978135ff6601a93930a2660165d507721d23a9d9bc15a96",
            "function": "WholeZRmCenteredInverse.owner",
            "AST_sha256": "e32cb6a603c756d6cc92de64a6616c213303e3a621f0ed350e854662dba202e3",
        }
    )
    inlet_refs = _handles(inlet)
    checks["inlet_in_accepted_native_prefix"] = bool(inlet_refs) and all(73190 <= ref < 97295 for ref in inlet_refs)
    checks["inlet_native_substitutions"] = bool(inlet_refs) and all(
        isinstance(graph_nodes[ref], Mapping)
        and graph_nodes[ref].get("operation") == "function_substitution"
        and graph_nodes[ref].get("variable") == 84528
        and graph_nodes[ref].get("value") == 1
        and graph_nodes[ref].get("Z_independent_substitution") is True for ref in inlet_refs)
    checks["chart_coordinate_and_measure"] = (
        chart.get("phase") == 185676 and chart.get("domain") == [1, 2948]
        and graph_nodes[185676] == {"operation": "bound_variable", "name": "true_Rm_patch_x"}
        and chart.get("radius") == 187214 and chart.get("physical_measure") == 187215
        and graph_nodes[187214].get("operation") == "product"
        and graph_nodes[187214].get("arguments") == [inlet.get("radius"), 185676]
        and graph_nodes[187215].get("numerator") == 1
        and graph_nodes[187215].get("denominator") == 185676
    )
    for label, value in (("Rm", 1), ("Rh", 2948)):
        endpoint = endpoints.get(label, {})
        checks[label + "_endpoint_substitution"] = all(
            isinstance(chart.get(group), Mapping) and isinstance(endpoint.get(group), Mapping)
            and set(chart[group]) == set(endpoint[group])
            and all(len(chart[group][key]) == len(endpoint[group][key])
                    and all(_substitution(graph_nodes, ref, expression, 185676, value)
                            for expression, ref in zip(chart[group][key], endpoint[group][key]))
                    for key in chart[group]) for group in ("fields", "radial", "generic")) and _substitution(
                        graph_nodes, endpoint.get("radius"), chart.get("radius"), 185676, value)
    for name, ok in checks.items():
        if ok is not True:
            gates.append("c6_archive_port_check_failed:" + name)
    valid = all(value is True for value in checks.values())
    accepted_c6_inlet = bool(
        valid and checks["accepted_Rm_inlet_endpoint_identity"]
        and checks["accepted_source_normalization_bridge"]
        and checks["accepted_endpoint_source_boundary"]
        and saved.get("accepted", False)
    )
    return {
        "valid": valid,
        "accepted_C6_Rm_inlet_identity": accepted_c6_inlet,
        "cross_owner_source_bridge": False,
        "cross_owner_source_bridge_missing": (
            "The consumed reports bind the native/C6 Rm inlet. They do not give a qualified "
            "map from feedback reference.inputs/defects and analytic P0 to those endpoint ports. "
            "Local recipe equality and the Rsh join do not supply that missing map."
        ),
        "archive_source": C6_ARCHIVE,
        "archive_decoded_sha256": C6_ARCHIVE_DECODED_SHA256,
        "graph_prefix_length": node_count,
        "checks": checks,
        "qualified_handles": {
            "archive_graph_prefix_length": node_count,
            "inlet": inlet,
            "Rm": endpoints.get("Rm", {}),
            "Rh": endpoints.get("Rh", {}),
            "chart": chart,
            "system": system_ref,
            "archive": C6_ARCHIVE,
            "archive_raw_sha256": C6_ARCHIVE_SHA256,
            "ports_present_not_admitted_when_checks_fail": not valid,
            "branch": {
                "actual_source_Z_domain": branch.get("actual_source_Z_domain"),
                "normalized_raw_histories": correspondence.get("normalized_raw_histories"),
                "actual_Rm_amplitude_identity": correspondence.get("actual_Rm_amplitude_identity"),
                "exact_identity": conjugacy.get("exact_identity"),
            },
            "accepted_Rm_inlet_endpoint_identity": checks["accepted_Rm_inlet_endpoint_identity"],
        },
    }


def _witnesses(texts, trees, saved, receipts):
    c6 = trees[C6_SOURCE]
    c2 = trees[C2_SOURCE]
    c2_init = _fn(c2, "__init__", "CompliantActualPatchMixedC4")
    c2_moment = trees[C2_MOMENT_SOURCE]
    c2_nonlinear = _fn(c2_moment, "nonlinear_series")
    feedback_moment = trees[C2_FEEDBACK_MOMENT_SOURCE]
    feedback = trees[C2_FEEDBACK_SOURCE]
    reference = trees[REFERENCE_MIXED_SOURCE]
    rsh = trees[RSH_JOIN_SOURCE]
    current_reference = trees[CURRENT_REFERENCE_SOURCE]
    profile = trees[REFERENCE_PROFILE_SOURCE]
    stress = trees[CURRENT_PATCH_STRESS_SOURCE]
    background = trees[BACKGROUND_SOURCE]

    c6_inlet = _fn(c6, "inlet_defects")
    c6_conjugacy = _fn(c6, "map_conjugacy")
    c6_fields = _fn(c6, "patch_fields")
    c6_nonlinear = _fn(c6, "nonlinear")
    c2_patch = _fn(c2, "patch_mixed")
    c2_eval = _fn(c2, "evaluate", "CompliantActualPatchMixedC4")
    moment_data = _fn(c2_moment, "actual_data", "CompliantActualMomentPatch")
    moment_coeff = _fn(c2_moment, "coefficients", "CompliantActualMomentPatch")
    moment_eval = _fn(c2_moment, "evaluate", "CompliantActualMomentPatch")
    feedback_moment_bindings = _fn(feedback_moment, "current_patch_source_bindings")
    feedback_mixed_bindings = _fn(feedback, "current_mixed_source_bindings")
    reference_e_bindings = _fn(reference, "current_E_source_bindings")
    reference_inputs = _fn(reference, "inputs", "CompliantActualReferenceRestoreProfiles")
    rsh_bindings = _fn(rsh, "source_bindings")
    rsh_boundary = _fn(rsh, "boundary_identities")
    rsh_build = _fn(rsh, "build")
    current_bindings = _fn(current_reference, "source_bindings")
    current_actual = _fn(current_reference, "__init__", "ActualReferenceRestoreFunctions")
    profile_restore = _fn(profile, "restoration", "CompliantReferenceRestoreProfiles")
    profile_terminal = _fn(profile, "terminal", "CompliantReferenceRestoreProfiles")
    profile_defects = _fn(profile, "defects", "CompliantReferenceRestoreProfiles")
    stress_rows = _fn(stress, "raw_patch_rows")
    stress_theorem = _fn(stress, "actual_patch_normalization_theorem")
    background_theorem = _fn(background, "current_patch_source_theorem")

    witness = {
        "c6_nonlinear_original_map": _has(c6_nonlinear,
            "cross=a.add(a.scale(a.mul(h[0],h[2]),W['fg'][0]),a.scale(a.mul(h[1],h[4]),W['fg'][1]))",
            "axial=a.add(a.scale(a.mul(h[0],h[0]),W['gg'][0]),a.scale(a.mul(h[1],h[1]),W['gg'][1]))",
            "a.sub(a.mul(alpha,axial),a.scale(swirl,'1/2'))",
            "a.scale(pressure,'1/2')"),
        "c6_inlet_defects_from_m_h_k_e_p": _has(c6_inlet,
            "mass=a.scale(fields['m'],Pstar)",
            "theta=a.div(fields['h'],am",
            "mixed=a.scale(a.div(fields['k'],am",
            "a.sub(mass,a.scale(z,4))",
            "a.sub(mixed,a.scale(a.mul(z,theta),4))",
            "a.sub(theta,a.fixed('5/8'))",
            "a.add(energy,a.fixed('5/12'))",
            "a.sub(pressure,a.fixed('5/2'))"),
        "c6_inlet_defects_exclude_P0": _has(c6_inlet, "pressure=a.div(fields['p'],a.mul(am,am)"),
        "c6_map_conjugacy_same_L_Q": _has(c6_conjugacy,
            "Wn={key:value if key in ('L','gg') else [g.mul(q,angular) for q in value] for key,value in W.items()}",
            "G(u)=S^-1 F(Su); J_G=S^-1 J_F S",
            "original_block_diagonal_L_commutes_with_S=True"),
        "c6_canonical_amplitude": _has(c6, "exp(-3/5)/(1+Z^2)", "positive canonical normalized Rm amplitude"),
        "c2_canonical_amplitude": _assign(moment_eval, "am", "(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))"),
        "c2_nonlinear_series_source": _has(_fn(c2_moment, "nonlinear_series"),
            "h[0]*h[2]*W['fg'][0]+h[1]*h[4]*W['fg'][1]",
            "(square(h[0])*W['gg'][0]+square(h[1])*W['gg'][1])*invAm2",
            "square(h[i+2])*W['ff'][i]", "square(h[i+2])*W['ff_over_x'][i]"),
        "c6_patch_fields_units_and_pressure": _has(c6_fields,
            "fields=dict(E=a.mul(am,H),V=a.scale(a.add(a.scale(z,4),axial),inv),m=a.scale(mass,inv)",
            "h=a.scale(a.mul(am,theta),angular)",
            "k=a.scale(a.mul(am,mixed),g.mul(angular,inv))",
            "e=rawenergy,p=a.mul(am2,pressure),P0=a.trim(P0)",
            "fields['pressure']=a.add(fields['P0'],fields['p'])"),
        "c2_nonlinear_series_source": _has(c2_nonlinear,
            "h[0]*h[2]*W['fg'][0]+h[1]*h[4]*W['fg'][1]",
            "(square(h[0])*W['gg'][0]+square(h[1])*W['gg'][1])*invAm2",
            "sum((square(h[i+2])*W['ff'][i] for i in range(3)),zero)/2",
            "sum((square(h[i+2])*W['ff_over_x'][i] for i in range(3)),zero)/2"),
        "c2_canonical_amplitude": _has(c2_eval,
            "am=(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))"),
        "c2_patch_mixed_current_radius_ode": _has(c2_patch,
            "G=[initial['mass']*x]+V[:4]",
            "mass=[sum((G[j]*inverses[k-j]*math.comb(k,j) for j in range(k+1)),zero) for k in range(5)]",
            "Q=[(2*z*V[k]-z*mass[k]*(1-delta)-d*derivative(mass[k]))/L for k in range(5)]",
            "Ur=[sum((Q[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),Q[0]*0) for k in range(5)]",
            "P=[p0+square(am)*pressure[0]]+[square(am)*row for row in pressure[1:]]"),
        "c2_patch_mixed_physical_primitives": _has(c2_patch,
            "Mtheta_over_sqrt2_Rm_1p5_Pstar",
            "Mtheta_z_over_sqrt2_Rm_1p5_Pstar",
            "Mztheta_over_Rm_Pstar2",
            "Mp_over_Pstar2"),
        "c2_beta_centers_radius": _has(c2_eval,
            "gamma=[self.gamma(x,c.mpf(i)/4) for i in (5,6,7)]") and _has(c2_init,
            "self.radius=c.mpf(1)/40"),
        "c2_feedback_original_mixed_aliases": _has(feedback_mixed_bindings,
            "getattr(CompliantActualFeedbackPatchMixedC4,name) is getattr(BASE,name)",
            "for name in ('gamma','evaluate')",
            "'p0':\"jet('original_P0_axial5')\"",
            "self.patch.coefficients(Z)"),
        "c2_feedback_original_history_aliases": _has(feedback_moment_bindings,
            "getattr(CompliantActualFeedbackMomentPatch,name) is getattr(BASE,name)",
            "for name in ('actual_data','coefficients','evaluate')",
            "'inp':'self.reference.inputs(Z)'",
            "'initial':\"inp['source_centered_Rsh']\"",
            "'E':\"inp['E']\"",
            "\"self.reference.defects(Z)['actual_normalized_five_defect_axial5_coefficients']\"",
            "'p0':\"self.reference.inputs(Z)['original_axis_pressure']\""),
        "c2_unchanged_original_moment_callables": all(_fn(c2_moment, name, "CompliantActualMomentPatch") is not None
                                                        for name in ("actual_data", "coefficients", "evaluate")),
        "reference_feedback_E_components": _has(reference_e_bindings,
            "assignment_source_bindings('actual_reference_restore_mixed_C4','inputs'",
            "'core_E':\"core_packet['E']\"",
            "'bridge_delta':\"IntervalTaylor(c,incoming['actual_delta_V_axial5'])\"",
            "'first_delta':\"IntervalTaylor(c,switch['velocity_increment'])\"",
            "'E':'core_E+bridge_delta+first_delta'"),
        "reference_E_source_inputs_exact": _has(reference_inputs,
            "core_packet=current_centered_core(self.history,Z)",
            "incoming=self.history.upstream.packet(Z,1,'macro')",
            "switch=self.history.switch.inputs(Z)",
            "E=core_E+bridge_delta+first_delta",
            "current_E_component_enclosures=dict(core=list(core_E.coefficients),bridge=list(bridge_delta.coefficients),first_switch=list(first_delta.coefficients)"),
        "reference_Rsh_centered_six_histories": _has(reference_inputs,
            "centered=dict(mean_error=moments['mean']-z*4",
            "angular_error=moments['theta']-c.mpf(5)/8",
            "mixed_error=moments['theta_z']-z*moments['theta']*4",
            "axial_square=moments['axial']-z*moments['mean']*8+square(z)*16",
            "swirl_error=moments['swirl']-c.mpf(5)/6",
            "pressure_error=moments['pressure']-5"),
        "reference_original_restore_E_and_terminal_V": _has(profile_restore,
            "V=inp['z']*4+inp['E']*alpha",
            "centered=restore_centered(c,initial,inp['E'],t,kernels)") and _has(profile_terminal,
            "inp['z']*4",
            "axial_source='V=4Z exactly; all actual moment differences remain'"),
        "reference_original_defect_projection": _has(profile_defects,
            "centered['mean_error']*c.exp(x)",
            "centered['mixed_error']*c.exp(c.mpf('1.6')*x)",
            "centered['angular_error']*c.exp(c.mpf('1.6')*x)",
            "centered['axial_square']*invAm2",
            "centered['swirl_error']*(c.exp(c.mpf('1.2')*x)/2)",
            "centered['pressure_error']*(c.exp(c.mpf('.2')*x)/2)"),
        "current_reference_native_six_history_source_bindings": _has(current_bindings,
            "dict(mean_error=moments['mean']-z*4,angular_error=moments['theta']-c.mpf(5)/8",
            "inp['z']*4+inp['E']*alpha",
            "dict(mean_error=1,angular_error='1.6',mixed_error='1.6',axial_square=1,swirl_error='1.2',pressure_error='.2')",
            "actual_E_source='same actual V_Rsh minus4Z, not a fitted centered cap'"),
        "current_reference_constructor_keeps_P0_and_histories": _has(current_actual,
            "set(actual_shapes)!=set(previous.NAMES)",
            "self.initial=actual_shapes",
            "self.V=V",
            "self.P0=P0",
            "self.E=f.add(V,f.scale(self.Vref,-1))",
            "self.centered=dict(mean_error=f.add(s['mean'],f.scale(self.Vref,-1))"),
        "rsh_functional_join_E_and_P0_source": _has(rsh_bindings,
            "current_six_centered_Rsh_coordinates",
            "original_same_parent_axis_pressure_keyword",
            "original_six_centered_source_equations") and _has(rsh_build,
            "current_Rsh_parent_and_axis_pressure_identical",
            "graph['E_V110_minus_4Z']['args']!=graph['V110']['args'][1:]",
            "canonical_source(parent['pressure_axis_axial5_coefficients'])",
            "current_V_E_common_signed_source_graph_verified=True"),
        "rsh_boundary_exact_six_history_transport": _has(rsh_boundary,
            "eh_y=eh*exp(-8*y/5)",
            "em_y=E+(em-E)*exp(-y)",
            "ek_y=ek*exp(-8*y/5)+E*(1-exp(-8*y/5))/s.Rational(8,5)",
            "aa_y=E**2+(aa-E**2)*exp(-y)",
            "pressure_is_Pstar_squared_times_same_axis_datum_plus_Mp=True"),
        "stress_adapter_current_radius_normalization": _has(stress_rows,
            "am=(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))",
            "energy=[(z*G[j]*8-square(z)*(16*xrows[j]))*invP2+square(am)*raw['e'][j]",
            "histories=dict(m=raw['m']",
            "p0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in"),
        "stress_adapter_saved_normalization_theorem": _has(stress_theorem,
            "rhs=dict(m=v-raw['m'][0].expr,h=u-s.Rational(3,2)*raw['h'][0].expr,k=u*v-s.Rational(3,2)*raw['k'][0].expr,e=v*v/Ps**2-u*u/2-raw['e'][0].expr,p=u*u/2)",
            "checks['raw_five_history_ODE_'+key]=True",
            "current_radius='R=Rm*x; Rm=Rref*exp(-6); D_y=x*D_x'",
            "absolute_P0_added_once_before_physical_differentiation=True"),
        "background_actual_owner_aliases": _has(background_theorem,
            "type(patch) is CompliantActualFeedbackPatchMixedC4",
            "type(patch.patch) is CompliantActualFeedbackMomentPatch",
            "unchanged_original_mixed_evaluate=patch.evaluate.__func__ is CompliantActualPatchMixedC4.evaluate",
            "unchanged_original_gamma=patch.gamma.__func__ is CompliantActualPatchMixedC4.gamma",
            "unchanged_original_patch_mixed=patch.evaluate.__func__.__globals__['patch_mixed'] is patch_mixed",
            "unchanged_original_actual_data=patch.patch.actual_data.__func__ is CompliantActualMomentPatch.actual_data",
            "unchanged_original_coefficients=patch.patch.coefficients.__func__ is CompliantActualMomentPatch.coefficients",
            "unchanged_original_moment_evaluate=patch.patch.evaluate.__func__ is CompliantActualMomentPatch.evaluate"),
        "background_same_P0_and_Rh_source": _has(background_theorem,
            "same_analytic_P0_function_and_first_six_projection_not_hash_only=True",
            "exact_patch_Rh_radius_identity='Rref*exp(-6)*e=Rref*exp(-5)'",
            "current_Rh_open_neighborhood_uses_same_unique_implicit_full_weight_closure=True",
            "interval_overlap_not_used_as_function_identity=True"),
        "background_source_normalization_theorem": _text_has(texts[BACKGROUND_SOURCE],
            "full_current_implicit_axial5_and_partial_integrals_retained=True",
            "full_energy_cross_terms_and_absolute_P0_retained=True",
            "ordinary_x_converted_to_logR_before_full_physical_derivatives=True"),
    }
    # Check receipt evidence separately so the static source witnesses remain
    # useful when a saved report is intentionally withheld or is stale.
    witness["saved_Rsh_join_receipt"] = receipts["checks"]["rsh_functional_join"]
    witness["saved_C6_actual_two_gates"] = receipts["checks"]["c6"]
    witness["saved_current_background_tensor_receipt"] = receipts["checks"]["background_tensor"]
    witness["saved_feedback_patch_receipts"] = receipts["checks"]["c2_feedback_moment"] and receipts["checks"]["c2_feedback_mixed"]
    witness["saved_reference_receipts"] = receipts["checks"]["reference_mixed"] and receipts["checks"]["current_reference_functions"]
    return {name: value is True for name, value in witness.items()}


def _qualified_ports(texts, trees, receipts, saved, witnesses, archive_ports):
    return {
        "C6": {
            "source": C6_SOURCE, "archive": C6_ARCHIVE,
            "inlet_defects": witnesses["c6_inlet_defects_from_m_h_k_e_p"],
            "nonlinear_map": witnesses["c6_nonlinear_original_map"],
            "map_conjugacy": witnesses["c6_map_conjugacy_same_L_Q"],
            "patch_fields": witnesses["c6_patch_fields_units_and_pressure"],
            "ordinary_Z_orders": [0, 1, 2],
            "branch_transfer": "original_common_box",
            "archive_port_checks": archive_ports.get("checks", {}),
            "qualified_handles": archive_ports.get("qualified_handles", {}),
            "cross_owner_source_bridge": archive_ports.get("cross_owner_source_bridge", False),
        },
        "C2": {
            "source": C2_SOURCE, "feedback_mixed_source": C2_FEEDBACK_SOURCE,
            "feedback_moment_source": C2_FEEDBACK_MOMENT_SOURCE,
            "patch_mixed": witnesses["c2_patch_mixed_current_radius_ode"],
            "original_gamma": witnesses["c2_beta_centers_radius"],
            "original_moment_callables": witnesses["c2_unchanged_original_moment_callables"],
            "pressure_P0": witnesses["c2_feedback_original_history_aliases"],
        },
        "Rsh_reference": {
            "source": CURRENT_REFERENCE_SOURCE, "profiles": REFERENCE_PROFILE_SOURCE,
            "source_join": RSH_JOIN_SOURCE,
            "feedback_E_components_recipe_checked": witnesses["reference_feedback_E_components"],
            "feedback_centered_six_history_recipe_checked": witnesses["reference_Rsh_centered_six_histories"],
            "Rsh_join_P0_recipe_checked": witnesses["rsh_functional_join_E_and_P0_source"],
        },
        "current_patch_tensor": {
            "source": BACKGROUND_SOURCE, "stress_adapter": CURRENT_PATCH_STRESS_SOURCE,
            "same_original_callables": witnesses["background_actual_owner_aliases"],
            "source_normalization": witnesses["stress_adapter_saved_normalization_theorem"],
        },
        "accepted_R110_Rm": {
            "report": ACCEPTED_R110_REPORT,
            "background": saved["background"], "own_history": saved["own_history"],
            "ordinary_Z_orders": [0, 1, 2],
            "accepted_Rm_inlet_endpoint_identity": archive_ports.get("checks", {}).get("accepted_Rm_inlet_endpoint_identity", False),
        },
    }


def _codehash():
    try:
        return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    except Exception:
        return None


def _result(nodes, actual_digest, gates, texts, saved, receipts, witnesses, archive_ports):
    ast_ok = bool(witnesses) and all(value is True for value in witnesses.values())
    receipt_ok = bool(receipts["checks"]) and all(receipts["checks"].values())
    graph_ok = isinstance(nodes, list) and len(nodes) == EXPECTED_NODE_COUNT and actual_digest == EXPECTED_NODES_DIGEST
    prerequisites_ok = bool(saved["accepted"] and ast_ok and receipt_ok and not gates)
    source_alias = bool(
        witnesses.get("reference_feedback_E_components")
        and witnesses.get("reference_E_source_inputs_exact")
        and witnesses.get("reference_Rsh_centered_six_histories")
        and witnesses.get("rsh_functional_join_E_and_P0_source")
        and witnesses.get("rsh_boundary_exact_six_history_transport")
        and witnesses.get("c2_feedback_original_history_aliases")
        and witnesses.get("c2_feedback_original_mixed_aliases")
        and receipts["checks"]["rsh_functional_join"]
        and archive_ports.get("cross_owner_source_bridge") is True
    )
    if not source_alias:
        gates.append("feedback_reference_Rm_source_alias_unproved")
    local_source_premises = {
        "c6_canonical_Am_source": witnesses.get("c6_canonical_amplitude", False),
        "c2_canonical_Am_source": witnesses.get("c2_canonical_amplitude", False),
        "c6_L_Q_map_source": witnesses.get("c6_map_conjugacy_same_L_Q", False),
        "c2_nonlinear_series_source": witnesses.get("c2_nonlinear_series_source", False),
        "c6_inlet_defect_source": witnesses.get("c6_inlet_defects_from_m_h_k_e_p", False),
        "c2_centered_defect_source": witnesses.get("reference_Rsh_centered_six_histories", False),
        "current_radius_source": witnesses.get("stress_adapter_current_radius_normalization", False),
        "feedback_E_component_source": witnesses.get("reference_feedback_E_components", False),
        "feedback_P0_source": witnesses.get("rsh_functional_join_E_and_P0_source", False),
    }
    source_conjugacy = {
        "local_source_premises": local_source_premises,
        "same_canonical_Am": bool(source_alias and local_source_premises["c6_canonical_Am_source"] and local_source_premises["c2_canonical_Am_source"]),
        "same_beta_centers": bool(source_alias and witnesses.get("c2_beta_centers_radius", False)),
        "same_L_Q_map": bool(source_alias and local_source_premises["c6_L_Q_map_source"] and local_source_premises["c2_nonlinear_series_source"] and witnesses.get("c2_patch_mixed_current_radius_ode", False)),
        "same_C2_C6_inlet_defect_vector": bool(source_alias and local_source_premises["c6_inlet_defect_source"] and local_source_premises["c2_centered_defect_source"] and witnesses.get("reference_original_defect_projection", False)),
        "same_current_radius_normalization": bool(source_alias and local_source_premises["current_radius_source"] and witnesses.get("c2_patch_mixed_current_radius_ode", False)),
        "same_E_components": bool(source_alias and local_source_premises["feedback_E_component_source"] and witnesses.get("reference_E_source_inputs_exact", False)),
        "same_centered_six_histories": bool(source_alias and local_source_premises["c2_centered_defect_source"] and witnesses.get("c2_feedback_original_history_aliases", False)),
        "same_P0_source": bool(source_alias and local_source_premises["feedback_P0_source"] and witnesses.get("c2_feedback_original_history_aliases", False) and witnesses.get("c2_feedback_original_mixed_aliases", False)),
        "feedback_reference_Rm_source_alias_proved": source_alias,
        "unresolved_feedback_reference_Rm_source_alias_unproved": not source_alias,
        "raw_to_C6_units": "E6=Utheta/Pstar; V6=Uz/Pstar; m6=mraw/Pstar; h6=hraw; k6=kraw/Pstar; e6=eraw; p6=praw; pressure=P0+p",
        "rates": dict(RATES),
        "ordinary_Z_orders": [0, 1, 2],
    }
    source_background = bool(
        graph_ok and prerequisites_ok and source_alias
        and witnesses.get("background_actual_owner_aliases")
        and witnesses.get("background_same_P0_and_Rh_source")
        and receipts["checks"]["background_tensor"]
    )
    source_own = bool(
        source_background
        and witnesses.get("reference_original_restore_E_and_terminal_V")
        and witnesses.get("reference_original_defect_projection")
        and witnesses.get("c6_patch_fields_units_and_pressure")
        and receipts["checks"]["c2_feedback_moment"]
        and receipts["checks"]["c2_feedback_mixed"]
    )
    all_premises = bool(source_background and source_own and not gates)
    blockers = sorted(set(gates)) if not all_premises else []
    ports = _qualified_ports(texts, {}, receipts, saved, witnesses, archive_ports)
    ports["C6"]["qualified_source_graph"] = archive_ports["qualified_handles"]
    source_audit = {
        "source_records_only": True,
        "static_source_AST_graph_and_saved_JSON_data_inspection_only": True,
        "nodes_hashed_only": True,
        "source_graph_substitutions_checked": True,
        "new_source_consumer_calls": 1,
        "accepted_ancestor_operator_calls": 0,
        "producer_provider_builder_checker_replayed": False,
        "quadrature_or_norm_calls": 0,
        "numeric_field_evaluated": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "finite_N_source_identity_claimed": False,
        "runtime_model_metadata": None,
    }
    claims = {
        "projection": (
            "actual C2 feedback patch source is conjugate to current true Rm C6 source through the common original map"
            if all_premises else "unproved: the feedback reference to native/C6 endpoint source map is missing"
        ),
        "statement_status": "proved" if all_premises else "local source recipes and conditional target only",
        "centered_inlet": "(m-4Z, K-4Z*H, H-5/8, alpha*axial_square-swirl_error/2, pressure_error/2)",
        "canonical_amplitude": "Am=exp(-3/5)/(1+Z^2)",
        "rates": dict(RATES),
        "beta": {"centers": list(CENTERS), "radius": "1/40", "support_edges": ["49/40", "51/40", "59/40", "61/40", "69/40", "71/40"]},
        "radius": "R=Rm*x, x in [1,e], y derivative D_y=x D_x, Rm=Rref*exp(-6)",
        "scope": "ordinary Z orders 0,1,2; no temporal n, finite-N density, Q_y/Z=1, or global five-moment closure",
    }
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": all_premises,
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "source_ast_witness": dict(witnesses),
        "source_ast_witness_all": ast_ok,
        "prerequisite_audit": {
            "accepted_R110_Rm_saved": saved["accepted"],
            "accepted_R110_background": saved["background"],
            "accepted_R110_own_history": saved["own_history"],
            "current_C6_two_actual_gates": receipts["checks"]["c6"],
            "current_Rsh_functional_join": receipts["checks"]["rsh_functional_join"],
            "current_background_tensor": receipts["checks"]["background_tensor"],
            "current_feedback_patch": receipts["checks"]["c2_feedback_mixed"] and receipts["checks"]["c2_feedback_moment"],
            "historical_reference_profiles_algorithm_only": receipts["checks"]["reference_profiles"],
            "finite_N_saved_prerequisite_not_used_for_patch_proof": True,
            "all_required_saved_prerequisites": prerequisites_ok,
        },
        "source_conjugacy": source_conjugacy,
        "c6_archive_port_audit": archive_ports,
        "checks": {
            "graph": graph_ok,
            "records_authenticated": not any(gate.startswith(("missing_source_record", "source_", "archive_", "decoded_archive", "unexpected_encoded")) for gate in gates),
            "saved_prerequisites": saved["accepted"],
            "source_AST_all": ast_ok,
            "receipts": receipt_ok,
            "feedback_reference_Rm_source_alias": source_alias,
            "patch_background_source_projection": source_background,
            "patch_own_leading_history_source_projection": source_own,
        },
        "qualified_function_ports": ports,
        "qualifiedfunctionports": ports,
        "patch_background_source_projection_proved": source_background,
        "patch_own_leading_history_source_projection_proved": source_own,
        "full_Rm_patch_C2_C6_source_projection_proved": all_premises,
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z_not_temporal_n": True,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "patch_finite_N_source_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "full_goal_complete": False,
        "source_inspection_audit": source_audit,
        "sourceaudit1new": 1,
        "sourceaudit0ancestor": 0,
        "source_audit_1_new": 1,
        "source_audit_0_ancestor": 0,
        "numeric": False,
        "runtime_modelmetadata": None,
        "exactclaims": claims,
        "exact_claims": claims,
        "blockers": blockers,
        "unresolved_gates": blockers,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "codehash": _codehash(),
        "worker_role": "luna_max_worker",
        "runtime_model_metadata": None,
    }


def derive(nodes, source_records, *, expected_nodes_digest):
    """Return a fail-closed, authenticated source projection for the Rm patch."""
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
    texts = {name: _record(records, name, gates) for name in RECORD_HASHES}
    saved = _saved_prerequisites(texts, gates)
    archive_ports = _c6_archive_ports(texts, nodes, saved, gates)
    receipts = _receipt_checks(texts, gates)
    trees = {name: _tree(source, name, gates) for name, source in texts.items() if name.endswith(".py")}
    witnesses = _witnesses(texts, trees, saved, receipts)
    for name, ok in witnesses.items():
        if ok is not True:
            gates.append("source_ast_witness_missing:" + name)
    return _result(nodes, actual_digest, gates, texts, saved, receipts, witnesses, archive_ports)
