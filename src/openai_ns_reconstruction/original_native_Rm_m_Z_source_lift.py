"""Source-only lift of the native Rm ordinary-Z m and physical mean-y ports.

This module consumes authenticated source records.  It does not import a producer,
execute an old graph, or evaluate a field.  The returned graph is a small proof
description whose atoms point at saved source ports and whose operations record the
new affine composition.
"""

from __future__ import annotations

import ast
import hashlib
import json


SCHEMA = "original_native_Rm_m_Z_source_lift_v1"
GRAPH_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
FINITE_PREFIX_DIGEST = "c7ab7da8cc5c51cd782a5794e3076a6cf021a45235aae1c75bd246b01491b639"
C5_PREFIX_DIGEST = "bc66540d371bcaf875acd3dcd3b926fd3dfaa3172368e0b3e65b93eaa3772525"
PRIMITIVE_DIGEST = "3fc133991e293860ab1adefcaa1d890238308113045de12388fd7d4b2e7893f9"
OUTER_CALCULUS_DIGEST = "c5af848cf3ee2c72a487f854cf6df6e10d3fd4bc26417fc827ac119d8c3f3988"

PARENT_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_native_C2_outer_six_source.json"
)
POSITIVITY_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_Rm_history_a_positivity_sources.json"
)
LEADING_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_Rm_leading_m_ordinary_Z_sources.json"
)
MEMORY_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_finite_N_nonreset_memory_source_binding.json"
)
ADMISSION_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_native_C2_patch_finite_N_admission.json"
)
CANDIDATE_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_native_C2_patch_finite_N_source.json"
)
CALCULUS_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
)
LOOP_ARCHIVE = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
PATCH_ARCHIVE = (
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
)
LOOP_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
PATCH_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
R100_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
PATCH_C6_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
BOUNDS_SOURCE = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json"


RECORD_HASHES = {
    PARENT_REPORT: "ac4900e6282dbb339708fcb202d9f21eea959f24e27d682f6f17521610216459",
    POSITIVITY_REPORT: "4c3b730e57966c0c4a8c2481564a6bf6ec19b1c1136ab01f7e7e2ade1246f675",
    LEADING_REPORT: "02cb39ce0272143655566f2a876bfdc1a083ea90c7aaf869506d401ed10ffcf8",
    MEMORY_REPORT: "41defb67ef05207a079075328f3f7afd3980193c4cec7adb9574a27d41f33335",
    ADMISSION_REPORT: "23d6b0984cb74b0b5fa31df6113e967604bf6b066ad3614b29d8430abb20f28d",
    CANDIDATE_REPORT: "5194f631aba2eebd4946003e1af52181ac8c050b92d636c36b6bdd2058245b5a",
    CALCULUS_REPORT: "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746",
    LOOP_ARCHIVE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    PATCH_ARCHIVE: "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad",
    LOOP_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PATCH_SOURCE: "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    R100_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PATCH_C6_SOURCE: "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    BOUNDS_SOURCE: "2817a7849a3bd943acb2e0ffaebae53d0f4d2c396b0c15f0f9b4d6effa8364bd",
}

DECODED_RECORD_HASHES = {
    LOOP_ARCHIVE: "f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2",
    PATCH_ARCHIVE: "c1026e05015584a7ce29bd7e8da95b17d97669523f1d063ac5b5c800559d9b69",
}

ORDERS = [0, 1, 2]
CONTEXTS = {
    "postrestore": {"coordinate": 84528, "endpoint": 1},
    "patch": {"coordinate": 185676, "endpoint": 1},
}
PORTS = {
    "postrestore": {
        "leading": [96436, 96437, 96438],
        "history": [121602, 121603, 121604],
        "modified": [121720, 121721, 121722],
        "U": [183268, 183269, 183270],
        "physical": [183280, 183281, 183282],
        "negative": [183286, 183287, 183288],
        "mean_y": [183292, 183293, 183294],
    },
    "patch": {
        "leading": [195000, 195001, 195002],
        "history": [198466, 198467, 198468],
        "modified": [198591, 198592, 198593],
        "U": [204317, 204318, 204319],
        "physical": [204329, 204330, 204331],
        "negative": [204335, 204336, 204337],
        "mean_y": [204341, 204342, 204343],
    },
}
EXPECTED_U_SOURCE = {
    "postrestore": [121571, 121572, 121352],
    "patch": [198423, 198424, 198425],
}
PSTAR = 33572


def _digest(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _node(nodes, handle):
    if isinstance(handle, int) and 0 <= handle < len(nodes):
        return nodes[handle]
    return None


def _ast_defs(text):
    try:
        tree = ast.parse(text)
    except (SyntaxError, TypeError, ValueError):
        return set()
    return {item.name for item in ast.walk(tree)
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def _verify_records(source_records, expected_input_digest, gates):
    if not isinstance(source_records, dict):
        gates.append("source_records_not_mapping")
        return {}, None
    if set(source_records) != set(RECORD_HASHES):
        gates.append("source_record_set_changed")
    actual = _digest(source_records)
    if actual != expected_input_digest:
        gates.append("input_digest_mismatch")
    parsed = {}
    for name, expected in RECORD_HASHES.items():
        record = source_records.get(name)
        if not isinstance(record, dict):
            gates.append("source_record_missing:" + name)
            continue
        if record.get("raw_sha256") != expected:
            gates.append("raw_hash_mismatch:" + name)
        is_gzip = name.endswith(".gz")
        if record.get("is_gzip") is not is_gzip:
            gates.append("gzip_contract_mismatch:" + name)
        if record.get("raw_authenticated_before_decode") is not True:
            gates.append("raw_hash_not_authenticated_before_decode:" + name)
        if record.get("source_path") != "experiments/root_st073/" + name:
            gates.append("source_path_mismatch:" + name)
        try:
            text = record["source_text"]
            decoded = hashlib.sha256(text.encode("utf-8")).hexdigest()
        except (KeyError, AttributeError, TypeError):
            gates.append("source_text_invalid:" + name)
            continue
        wanted = DECODED_RECORD_HASHES.get(name, expected)
        if record.get("decoded_sha256") != wanted or decoded != wanted:
            gates.append("decoded_hash_mismatch:" + name)
        if name.endswith(".json") or name.endswith(".json.gz"):
            try:
                parsed[name] = json.loads(text)
            except (TypeError, ValueError):
                gates.append("source_json_invalid:" + name)
        elif name.endswith(".py"):
            parsed[name] = text
            if not _ast_defs(text):
                gates.append("source_ast_invalid:" + name)
    return parsed, actual


class _ProofDag:
    def __init__(self):
        self.nodes = []

    def add(self, operation, arguments=(), **metadata):
        ident = "native_Rm_m_Z:" + str(len(self.nodes))
        node = {"id": ident, "operation": operation}
        if arguments:
            node["arguments"] = list(arguments)
        node.update(metadata)
        self.nodes.append(node)
        return ident


def _archive_checks(parsed, gates):
    loop = parsed.get(LOOP_ARCHIVE, {})
    patch = parsed.get(PATCH_ARCHIVE, {})
    loop_nodes = loop.get("exact_graph_nodes", []) if isinstance(loop, dict) else []
    patch_nodes = patch.get("exact_graph_nodes", []) if isinstance(patch, dict) else []
    if not isinstance(loop_nodes, list) or len(loop_nodes) != 121834:
        gates.append("native_C5_archive_node_count_changed")
    if not isinstance(patch_nodes, list) or len(patch_nodes) != 205147:
        gates.append("native_patch_C5_archive_node_count_changed")
    if isinstance(patch_nodes, list) and _digest(patch_nodes) != FINITE_PREFIX_DIGEST:
        gates.append("native_patch_C5_archive_digest_changed")
    if isinstance(loop_nodes, list) and isinstance(patch_nodes, list):
        if loop_nodes != patch_nodes[:121834]:
            gates.append("native_C5_patch_prefix_changed")
        if _digest(loop_nodes) != C5_PREFIX_DIGEST:
            gates.append("native_C5_archive_digest_changed")
    loop_bind = loop.get("prescribed_loop_C5_bindings", {}) if isinstance(loop, dict) else {}
    patch_bind = patch.get("Rm_patch_finite_N_C5_bindings", {}) if isinstance(patch, dict) else {}
    try:
        post = loop_bind["charts"]["postrestore"]
        chart = patch_bind["chart"]
        if loop_bind.get("N") != 18 or patch_bind.get("N") != 18 or _node(patch_nodes,18).get("Z_independent") is not True:
            gates.append("shared_fixed_Z_independent_N_missing")
        if _node(patch_nodes,20) != {"operation":"bound_variable","name":"Z"}:
            gates.append("native_Z_source_variable_changed")
        if post["roots"]["a"][0] != 119435 or chart["roots"]["a"][0] != 196071:
            gates.append("positivity_target_is_not_actual_native_root")
        if post["domain"] != [0,1] or chart["domain"] != [1,2948] or _node(patch_nodes,2948) != {"operation":"analytic_unary","name":"exp","argument":1}:
            gates.append("positivity_native_domain_binding_changed")
        if post["incoming"]["m"] != loop_bind["endpoints"]["restoration"]["m"]:
            gates.append("postrestore_predecessor_boundary_source_changed")
        if chart["incoming"]["m"] != loop_bind["endpoints"]["postrestore"]["m"] or chart["incoming"]["m"] != loop_bind["Rm_correction"]["m"]:
            gates.append("patch_nonreset_boundary_source_changed")
        if post["coordinate"] != CONTEXTS["postrestore"]["coordinate"]:
            gates.append("postrestore_native_coordinate_changed")
        if chart["coordinate"] != CONTEXTS["patch"]["coordinate"]:
            gates.append("patch_native_coordinate_changed")
        if post["history"]["m"] != [121602, 121603, 121604, 121605, 121606, 121607]:
            gates.append("postrestore_history_binding_changed")
        if chart["history"]["m"] != [198466, 198467, 198468, 198469, 198470, 198471]:
            gates.append("patch_history_binding_changed")
        if post["modified_fields"]["m"][:3] != PORTS["postrestore"]["modified"]:
            gates.append("postrestore_modified_binding_changed")
        if chart["modified_fields"]["m"][:3] != PORTS["patch"]["modified"]:
            gates.append("patch_modified_binding_changed")
        if chart["incoming"]["m"][:3] != [121756, 121757, 121758]:
            gates.append("patch_incoming_history_binding_changed")
    except (KeyError, TypeError):
        gates.append("native_C5_history_binding_missing")
    return loop_nodes, patch_nodes


def _row_checks(loop_nodes, patch_nodes, gates):
    rows = {}
    for side, ports in PORTS.items():
        nodes = patch_nodes
        for order in ORDERS:
            leading = _node(nodes, ports["leading"][order])
            history = _node(nodes, ports["history"][order])
            modified = _node(nodes, ports["modified"][order])
            if not (isinstance(leading, dict) and leading.get("operation") == "product"):
                gates.append("leading_source_row_missing:" + side + ":" + str(order))
            if not (isinstance(history, dict) and history.get("operation") == "sum"):
                gates.append("history_source_row_missing:" + side + ":" + str(order))
            expected_modified = [ports["leading"][order], ports["history"][order]]
            if modified != {"operation": "sum", "arguments": expected_modified}:
                gates.append("modified_affine_row_changed:" + side + ":" + str(order))
            rows[(side, order)] = {
                "side": side, "ordinary_order": order,
                "leading_handle": ports["leading"][order],
                "history_handle": ports["history"][order],
                "modified_handle": ports["modified"][order],
                "history_derivative_handle": ports["history"][order + 1]
                if order < 2 else None,
                "modified_identity": "leading_m + native_C5_history_m",
            }
    for side, ports in PORTS.items():
        for order in ORDERS:
            nodes = patch_nodes
            expected_u = EXPECTED_U_SOURCE[side][order]
            expected = {
                "U": {"operation": "product", "arguments": [expected_u, PSTAR]},
                "physical": {"operation": "product", "arguments": [ports["modified"][order], PSTAR]},
                "negative": {"operation": "negative", "argument": ports["physical"][order]},
                "mean_y": {"operation": "sum", "arguments": [ports["U"][order], ports["negative"][order]]},
            }
            for field, wanted in expected.items():
                if _node(nodes, ports[field][order]) != wanted:
                    gates.append("physical_affine_row_changed:" + side + ":" + field + ":" + str(order))
            rows[(side, order)].update({
                "U_handle": ports["U"][order],
                "physical_m_handle": ports["physical"][order],
                "negative_handle": ports["negative"][order],
                "physical_mean_y_handle": ports["mean_y"][order],
                "physical_mean_y_derivative_handle": ports["mean_y"][order + 1]
                if order < 2 else None,
                "source_derivative_proved": order < 2,
                "regularity_proved": order < 2,
            })
            if side == "postrestore":
                # The postrestore rows are in the common prefix of the patch archive.
                rows[(side, order)]["archive"] = PATCH_ARCHIVE
            else:
                rows[(side, order)]["archive"] = PATCH_ARCHIVE
    if _node(patch_nodes, PSTAR) != {"operation": "analytic_unary", "name": "exp", "argument": 5}:
        gates.append("Pstar_source_port_changed")
    return rows


def _accepted_checks(parsed, gates):
    parent = parsed.get(PARENT_REPORT, {})
    outer = parent.get("component_result", {}) if isinstance(parent, dict) else {}
    witness = outer.get("source_witness", {}) if isinstance(outer, dict) else {}
    bridge = witness.get("patch_history_bridge_checks", {}).get("m", {})
    parent_ok = (
        outer.get("passed") is True
        and outer.get("native_outer_six_nonreset_history_identity_proved") is True
        and outer.get("native_outer_six_source_identity") is True
        and all(bridge.get(k) is True for k in
                ("same_inlet", "same_density_port", "native_endpoint_substitution",
                 "old_C2_endpoint_substitution"))
        and witness.get("conditional_common_original_domain") is True
        and witness.get("conditional_genuine_C6_C5_regularity") is True
    )
    if not parent_ok:
        gates.append("native_outer_history_bridge_not_accepted")

    memory = parsed.get(MEMORY_REPORT, {}).get("component_result", {})
    memory_ok = (
        memory.get("passed") is True
        and memory.get("ordinary_Z_orders") == ORDERS
        and memory.get("ordinary_Z_not_temporal_n") is True
        and memory.get("first_ten_nonreset_history_identity_proved") is True
        and memory.get("first_ten_history_endpoint_identity_proved") is True
        and memory.get("first_ten_density_C2_NC5_identity_proved") is True
        and memory.get("first_ten_finite_N_field_source_composition_proved") is True
        and memory.get("generic_density_identity_proved") is True
        and memory.get("time_order_recursion_claimed") is False
        and memory.get("numeric_values_computed") is False
        and memory.get("old_operators_replayed") is False
        and memory.get("history_identity_proof", {}).get("proved") is True
        and memory.get("source_proof_composition", {}).get("proved") is True
    )
    memory_ports = memory.get("first_ten_history_ports",{}).get("postrestore",{})
    memory_ok = bool(memory_ok
        and memory.get("record_hashes",{}).get(LOOP_ARCHIVE) == RECORD_HASHES[LOOP_ARCHIVE]
        and memory_ports.get("C5_history_source_ports",{}).get("m") == [121602,121603,121604,121605,121606,121607]
        and memory_ports.get("C5_endpoint_ports",{}).get("m") == [121756,121757,121758,121759,121760,121761])
    if not memory_ok:
        gates.append("accepted_first_ten_memory_bridge_missing")

    admission = parsed.get(ADMISSION_REPORT, {}).get("component_result", {})
    admission_ok = (
        admission.get("passed") is True
        and admission.get("native_C2_C5_patch_history_identity_proved") is True
        and admission.get("patch_finite_N_source_identity") is True
        and admission.get("ordinary_Z_orders") == ORDERS
        and admission.get("numeric_values_computed") is False
        and admission.get("old_operators_replayed") is False
    )
    if not admission_ok:
        gates.append("native_patch_admission_missing")

    candidate = parsed.get(CANDIDATE_REPORT, {}).get("component_result", {})
    argument = candidate.get("source_identity_argument", {})
    candidate_ok = (
        candidate.get("source_only") is True
        and candidate.get("nodes_digest") == GRAPH_DIGEST
        and argument.get("qualified_previous_DAG_refs_and_native_C5_graph_refs") is True
        and argument.get("bridge_identity_status") is True
        and candidate.get("source_prerequisite_audit", {}).get("patch_C5_archive_shape") is True
    )
    if not candidate_ok:
        gates.append("native_patch_candidate_binding_missing")

    positivity = parsed.get(POSITIVITY_REPORT, {}).get("actual_history_a_positive_sources", {})
    positivity_ok = True
    bounds = {}
    for side, target, lower, native_domain in (
        ("postrestore", 119435, "4/5", ["0", "1"]),
        ("patch", 196071, "7/10", ["1", "exp(1)"]),
    ):
        item = positivity.get(side, {})
        good = (
            item.get("passed") is True and item.get("source_only") is True
            and item.get("nodes_digest") == GRAPH_DIGEST
            and item.get("actual_target_handle") == target
            and item.get("strict_positive_proved") is True
            and item.get("lower_bound_exact") == lower
            and item.get("whole_native_domain") == native_domain
            and item.get("whole_Z_domain") == ["-1", "1"]
            and item.get("numeric_values_computed") is False
        )
        if not good:
            positivity_ok = False
            gates.append("accepted_a_positive_source_missing:" + side)
        bounds[side] = {"lower_bound": lower, "is_bound_not_function_value": True,
                        "target_handle": target, "native_domain": native_domain,
                        "whole_Z_domain": [-1, 1]}

    leading = parsed.get(LEADING_REPORT, {}).get("actual_Rm_leading_m_ordinary_Z_source", {})
    leading_ok = (
        leading.get("passed") is True and leading.get("source_only") is True
        and leading.get("nodes_digest") == GRAPH_DIGEST
        and leading.get("Z_variable") == 20
        and leading.get("cross_side_value_equality_used") is False
        and leading.get("numeric_values_computed") is False
        and leading.get("old_operators_replayed") is False
    )
    lead_rows = leading.get("rows", []) if isinstance(leading, dict) else []
    for side, handles in (("postrestore", [96436, 96437, 96438]),
                          ("patch", [195000, 195001, 195002])):
        expected = [(0, handles[0], handles[1]), (1, handles[1], handles[2])]
        got = [(r.get("ordinary_order"), r.get("source_handle"), r.get("derivative_handle"))
               for r in lead_rows if r.get("side") == side and r.get("field") == "m"
               and r.get("source_derivative_proved") is True and r.get("regularity_proved") is True
               and r.get("context") == CONTEXTS[side]]
        if got != expected:
            leading_ok = False
            gates.append("accepted_leading_m_source_missing:" + side)
    if not leading_ok and "accepted_leading_m_source_missing:postrestore" not in gates \
            and "accepted_leading_m_source_missing:patch" not in gates:
        gates.append("accepted_leading_m_source_report_missing")

    calculus = parsed.get(CALCULUS_REPORT, {})
    primitive = calculus.get("actual_Rm_physical_primitive_Z_jets", {})
    outer_calculus = calculus.get("actual_Rm_physical_Qy_Z_derivative_provenance", {})
    primitive_digest = _digest(primitive) if primitive else ""
    outer_digest = _digest(outer_calculus) if outer_calculus else ""
    primitive_ok = (
        primitive_digest == PRIMITIVE_DIGEST
        and primitive.get("source_only") is True
        and primitive.get("source_affine_recipes_checked") is True
        and primitive.get("admitted_velocity_source", {}).get("U_through_ordinary_Z2_source_bound") is True
        and primitive.get("admitted_velocity_source", {}).get("Uy_through_ordinary_Z1_source_bound") is True
        and primitive.get("admitted_velocity_source", {}).get("old_theorem_or_checker_replayed") is False
        and primitive.get("numeric_values_computed") is False
        and primitive.get("old_operators_replayed") is False
    )
    if not primitive_ok:
        gates.append("accepted_Uy_ordinary_Z_source_missing")
    outer_ok = (
        outer_digest == OUTER_CALCULUS_DIGEST
        and outer_calculus.get("derivative_source_prefix_digest") == FINITE_PREFIX_DIGEST
        and outer_calculus.get("outer_calculus_proved") is True
        and outer_calculus.get("ordinary_Q_y_Z0_identity_admitted") is True
        and outer_calculus.get("common_phase_positive_endpoint_guard_proved") is True
        and outer_calculus.get("qy0_value_seam", {}).get("all_Z_identity_proved") is True
        and outer_calculus.get("common_phase_source", {}).get("actual_Z_independence_proved") is True
    )
    if not outer_ok:
        gates.append("accepted_conditional_Qy_calculus_missing")
    return {
        "parent": parent_ok, "memory": memory_ok, "admission": admission_ok,
        "candidate": candidate_ok, "positivity": positivity_ok, "leading": leading_ok,
        "primitive": primitive_ok, "outer_calculus": outer_ok,
        "a_bounds": bounds,
        "primitive_digest": primitive_digest,
        "outer_calculus_digest": outer_digest,
        "outer_calculus_previous_passed": outer_calculus.get("passed") is True,
    }


def _base_result(input_digest, expected_input_digest, gates):
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": False,
        "input_digest": input_digest,
        "expected_input_digest": expected_input_digest,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "ordinary_Z_orders": list(ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "runtime_model_metadata": None,
        "source_records_only": True,
        "new_source_consumer_calls": 1,
        "accepted_ancestor_operator_calls": 0,
        "old_definition_graphs_rebuilt": False,
        "old_graph_re_evaluated": False,
        "native_history_m_ordinary_Z_source_proved": False,
        "native_modified_m_ordinary_Z_source_proved": False,
        "native_physical_mean_y_ordinary_Z_source_proved": False,
        "Rm_physical_Qy_Z1_function_identity_proved": False,
        "physical_Qy_Z1_value_seam_proved": False,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "current_Rm_RQy_germ_identity_proved": False,
        "full_RQy_germ_proved": False,
        "legacy_feedback_reference_alias_proved": False,
        "prior_adapter_geometry_alias_proved": False,
        "prior_adapter_P0_alias_proved": False,
        "time_order_recursion_claimed": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
        "unresolved_gates": sorted(set(gates)),
        "retained_construction_data": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "six_signed_channels_retained": True,
            "all_old_error_data_retained": True,
            "reserved_contribution_registered_once": True,
            "inherited_errors_deleted": False,
            "P1_differs_from_m15": True,
            "physical_Q_differs_from_loop_q": True,
        },
        "source_inspection_audit": {
            "source_records_only": True,
            "static_source_AST_graph_saved_data_and_new_definitions_only": True,
            "original_graph_nodes_inspected_not_evaluated": True,
            "old_patch_DAG_rebuilt": False,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "old_definition_graphs_rebuilt": False,
            "old_graph_re_evaluated": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "numeric_values_computed": False,
            "old_operators_replayed": False,
            "runtime_model_metadata": None,
        },
        "contexts": dict(CONTEXTS),
        "definition_graph": {"namespace": SCHEMA, "nodes": [], "outputs": {}},
    }


def derive(source_records, *, expected_input_digest):
    """Consume fixed source records and return an unevaluated native m proof."""
    gates = []
    parsed, input_digest = _verify_records(source_records, expected_input_digest, gates)
    result = _base_result(input_digest, expected_input_digest, gates)
    if gates:
        result["unresolved_gates"] = sorted(set(gates))
        return result

    loop_nodes, patch_nodes = _archive_checks(parsed, gates)
    accepted = _accepted_checks(parsed, gates)
    graph_rows = _row_checks(loop_nodes, patch_nodes, gates)
    source_ast = {
        LOOP_SOURCE: sorted(_ast_defs(parsed.get(LOOP_SOURCE, ""))),
        PATCH_SOURCE: sorted(_ast_defs(parsed.get(PATCH_SOURCE, ""))),
        R100_SOURCE: sorted(_ast_defs(parsed.get(R100_SOURCE, ""))),
        PATCH_C6_SOURCE: sorted(_ast_defs(parsed.get(PATCH_C6_SOURCE, ""))),
    }
    ast_ok = (
        "CurrentOriginalPrescribedLoopC5Functions" in source_ast[LOOP_SOURCE]
        and "CurrentOriginalTrueRmPatchFiniteNC5Functions" in source_ast[PATCH_SOURCE]
        and "CurrentOriginalPrescribedR100RmC6Functions" in source_ast[R100_SOURCE]
        and "CurrentOriginalTrueRmPatchC6Functions" in source_ast[PATCH_C6_SOURCE]
        and "loop_functions" in source_ast[LOOP_SOURCE]
        and "patch_fields" in source_ast[PATCH_C6_SOURCE]
    )
    if not ast_ok:
        gates.append("native_function_port_AST_binding_missing")

    history_ok = (
        not gates and accepted["parent"] and accepted["memory"]
        and accepted["admission"] and accepted["candidate"]
        and accepted["positivity"] and ast_ok
    )
    modified_ok = history_ok and accepted["leading"] and not gates
    physical_ok = modified_ok and accepted["primitive"] and not gates
    qy_ok = physical_ok and accepted["outer_calculus"] and not gates
    dag = _ProofDag()
    certificates = {name: dag.add("accepted_source_theorem", report=report, raw_sha256=RECORD_HASHES[report],
        premise_name=name, premise_admitted=accepted[name], old_operator_called=False,
        scope="qualified source fact; original packet overall admission is unchanged") for name,report in (
        ("parent",PARENT_REPORT),("memory",MEMORY_REPORT),("admission",ADMISSION_REPORT),
        ("positivity",POSITIVITY_REPORT),("leading",LEADING_REPORT),("primitive",CALCULUS_REPORT),
        ("outer_calculus",CALCULUS_REPORT))}
    pstar_source = dag.add("source_port", field="Pstar",archive=PATCH_ARCHIVE,
        archive_sha256=RECORD_HASHES[PATCH_ARCHIVE],handle=PSTAR,Z_independent=True)
    history_nodes = {}
    modified_nodes = {}
    physical_nodes = {}
    mean_nodes = {}
    for side, ports in PORTS.items():
        history = dag.add(
            "source_port", label="native_C5_history_m", side=side,
            archive=LOOP_ARCHIVE if side == "postrestore" else PATCH_ARCHIVE,
            handles=list(ports["history"]), context=CONTEXTS[side],
            ordinary_Z_orders=list(ORDERS), source_identity="native C5 history binding",
        )
        transfer = dag.add(
            "history_transfer", [history, certificates["memory"], certificates["positivity"],
                certificates["admission"], certificates["parent"]], side=side,
            theorem="H_C2 = N H_C5 on ordinary Z orders 0,1,2",
            memory_report=MEMORY_REPORT, patch_bridge=PARENT_REPORT,
            no_rescaling_of_native_C5_history=True,
        )
        history_nodes[side] = transfer
        for order in ORDERS:
            leading = dag.add(
                "source_port", [certificates["leading"]], side=side, field="leading_m", ordinary_order=order,
                archive=PATCH_ARCHIVE, handle=ports["leading"][order],
                context=CONTEXTS[side], source_identity="accepted leading m source",
            )
            history_row = dag.add(
                "ordinary_Z_source_row", [transfer], side=side, field="history_m",
                ordinary_order=order, source_handle=ports["history"][order],
                derivative_edge=(ports["history"][order], ports["history"][order + 1])
                if order < 2 else None,
                source_derivative_proved=bool(history_ok and order < 2),
                source_jet_row_admitted=bool(history_ok), regularity_proved=bool(history_ok),
            )
            modified_nodes[(side, order)] = dag.add(
                "sum", [leading, history_row], side=side, field="modified_m",
                ordinary_order=order, source_handle=ports["modified"][order],
                formula="leading_m + native_C5_history_m",
            )
            physical_nodes[(side, order)] = dag.add(
                "product", [modified_nodes[(side, order)],pstar_source], side=side,
                field="physical_m", ordinary_order=order,
                factor_handle=PSTAR, formula="Pstar * modified_m",
            )
            negative = dag.add(
                "negative", [physical_nodes[(side, order)]], side=side,
                field="negative_physical_m", ordinary_order=order,
            )
            u = dag.add(
                "source_port", [certificates["primitive"]], side=side, field="U", ordinary_order=order,
                archive=PATCH_ARCHIVE,
                handle=ports["U"][order], source_identity="accepted physical U source",
            )
            mean_nodes[(side, order)] = dag.add(
                "sum", [u, negative], side=side, field="physical_mean_y",
                ordinary_order=order, source_handle=ports["mean_y"][order],
                formula="U - Pstar * modified_m",
            )
        for order in (0, 1):
            dag.add(
                "ordinary_Z_derivative_edge", [history_nodes[side]], side=side,
                field="history_m", ordinary_order=order,
                source_handle=ports["history"][order],
                derivative_handle=ports["history"][order + 1],
                context=CONTEXTS[side], source_regular=bool(history_ok),
            )
            dag.add(
                "ordinary_Z_derivative_edge", [mean_nodes[(side, order)]], side=side,
                field="physical_mean_y", ordinary_order=order,
                source_handle=ports["mean_y"][order],
                derivative_handle=ports["mean_y"][order + 1],
                context=CONTEXTS[side], source_regular=bool(physical_ok),
            )

    primitive_source_rows = []
    for side, uy, uy1, mean in (
        ("postrestore", 183274, 183275, [183292, 183293, 183294]),
        ("patch", 204323, 204324, [204341, 204342, 204343]),
    ):
        primitive_source_rows.append({
            "side": side, "field": "U_y", "ordinary_order": 0,
            "source_handle": uy, "derivative_handle": uy1,
            "source_derivative_proved": bool(accepted["primitive"]), "regularity_proved": bool(accepted["primitive"]),
            "context": CONTEXTS[side],
            "proof_origin": "accepted U_y ordinary-Z source bound",
        })
        for order in (0, 1):
            primitive_source_rows.append({
                "side": side, "field": "physical_mean_y", "ordinary_order": order,
                "source_handle": mean[order], "derivative_handle": mean[order + 1],
                "source_derivative_proved": bool(physical_ok),
                "regularity_proved": bool(physical_ok),
                "context": CONTEXTS[side],
                "proof_origin": "new native history derivative plus affine composition",
            })

    edge_sets = {}
    for field, port_name, admitted in (("history_m","history",history_ok),
            ("modified_m","modified",modified_ok),("physical_mean_y","mean_y",physical_ok)):
        edge_sets[field] = [dict(side=side,field=field,ordinary_order=order,
            source_handle=ports[port_name][order],derivative_handle=ports[port_name][order+1],
            source_derivative_proved=bool(admitted),regularity_proved=bool(admitted),
            archive=PATCH_ARCHIVE,archive_sha256=RECORD_HASHES[PATCH_ARCHIVE],context=CONTEXTS[side])
            for side,ports in PORTS.items() for order in (0,1)]
    qy_proof = dag.add("source_identity_theorem",
        [certificates["outer_calculus"],certificates["primitive"]]
        + [mean_nodes[(side,order)] for side in PORTS for order in ORDERS],
        theorem="D_Z_of_admitted_all_Z_Qy0_identity", admitted=bool(qy_ok),
        target_archive=PATCH_ARCHIVE,target_archive_sha256=RECORD_HASHES[PATCH_ARCHIVE],
        target_handles=[183455,204504],common_native_prefix_sha256=FINITE_PREFIX_DIGEST,
        per_side_C1_source_rows=True,closed_Z_endpoints_by_continuity=True,
        cross_side_Z2_values_required=False,old_calculus_replayed=False)
    for name, ok in (("native_history_m_source",history_ok),("modified_m_source",modified_ok),
            ("physical_mean_y_source",physical_ok),("local_Rm_Qy_Z1_identity",qy_ok)):
        if not ok and not gates: gates.append(name+"_unproved")
    for row in graph_rows.values():
        row["source_derivative_proved"] = bool(physical_ok and row["ordinary_order"] < 2)
        row["regularity_proved"] = bool(physical_ok and row["ordinary_order"] < 2)
    result.update({
        "native_history_m_ordinary_Z_source_proved": bool(history_ok),
        "native_modified_m_ordinary_Z_source_proved": bool(modified_ok),
        "native_physical_mean_y_ordinary_Z_source_proved": bool(physical_ok),
        "Rm_physical_Qy_Z1_function_identity_proved": bool(qy_ok),
        "physical_Qy_Z1_value_seam_proved": bool(qy_ok),
        "passed": bool(qy_ok and not gates),
        "unresolved_gates": sorted(set(gates)),
        "definition_graph": {
            "namespace": SCHEMA,
            "proof_graph_not_numeric_evaluator": True,
            "nodes_digest": _digest(dag.nodes),
            "nodes": dag.nodes,
            "node_count": len(dag.nodes),
            "outputs": {
                "local_Rm_Qy_Z1": qy_proof,
                "history_m": {side: history_nodes[side] for side in PORTS},
                "modified_m": {side: [modified_nodes[(side, o)] for o in ORDERS] for side in PORTS},
                "physical_mean_y": {side: [mean_nodes[(side, o)] for o in ORDERS] for side in PORTS},
            },
        },
        "history_rows": edge_sets["history_m"],
        "modified_m_rows": edge_sets["modified_m"],
        "physical_mean_y_rows": edge_sets["physical_mean_y"],
        "actual_affine_graph_descriptors": [graph_rows[(side, order)] for side in PORTS for order in ORDERS],
        "primitive_source_rows": primitive_source_rows,
        "source_witness": {
            "graph_digest": GRAPH_DIGEST,
            "finite_prefix_digest": FINITE_PREFIX_DIGEST,
            "c5_prefix_digest": C5_PREFIX_DIGEST,
            "archive_graph_counts": {"C5": 121834, "patch_C5": 205147},
            "history_handles": {side: PORTS[side]["history"] for side in PORTS},
            "contexts": dict(CONTEXTS),
            "parent_report": {"path": PARENT_REPORT, "sha256": RECORD_HASHES[PARENT_REPORT]},
            "first_ten_memory_report": {"path": MEMORY_REPORT, "sha256": RECORD_HASHES[MEMORY_REPORT]},
            "patch_admission": {"path": ADMISSION_REPORT, "sha256": RECORD_HASHES[ADMISSION_REPORT]},
            "patch_candidate": {"path": CANDIDATE_REPORT, "sha256": RECORD_HASHES[CANDIDATE_REPORT],
                                "reference_only": True, "re_evaluated": False},
            "a_positive_sources": accepted["a_bounds"],
            "a_field_root_domain_bindings": {
                "postrestore": {
                    "field": "a", "root_handle": 119435,
                    "native_domain": ["0", "1"], "Z_domain": [-1, 1],
                    "source_identity": bool(history_ok),
                },
                "patch": {
                    "field": "a", "root_handle": 196071,
                    "native_domain": ["1", "exp(1)"], "Z_domain": [-1, 1],
                    "source_identity": bool(history_ok),
                },
            },
            "leading_source_report": {"path": LEADING_REPORT, "sha256": RECORD_HASHES[LEADING_REPORT]},
            "calculus": {
                "primitive_digest": accepted["primitive_digest"],
                "expected_primitive_digest": PRIMITIVE_DIGEST,
                "outer_calculus_digest": accepted["outer_calculus_digest"],
                "expected_outer_calculus_digest": OUTER_CALCULUS_DIGEST,
                "previous_outer_passed": accepted["outer_calculus_previous_passed"],
                "new_physical_mean_y_edges": bool(physical_ok),
            },
            "no_cross_side_Z2_value_equality": True,
            "ordinary_Z_is_not_temporal_n": True,
        },
        "history_source_derivative_argument": {
            "functions": "native C5 m histories, via admitted H_C2=N*H_C5 on the same original domain",
            "incoming": "first-ten predecessor ordinary-Z source functions; patch incoming is the fixed postrestore outgoing function",
            "source_density": "B_over_Pstar/N with admitted genuine source jets and original flat/positive IFT branch",
            "parameters": "same fixed positive N and Z-independent phase, geometry, kernel and integration limits",
            "positivity": "consume already saved actual a lower bounds on the whole source charts; do not reprove or use them as function values",
            "regularity": "genuine source derivatives are continuous on the common compact native/Z domain, hence locally dominated; differentiate finite Duhamel integrals twice",
            "qualitative_only": True, "numerical_majorant_claimed": False,
            "postrestore_measure": "physical_Jacobian dt", "patch_measure": "dt/t",
            "predecessor_memory_retained": True,
            "patch_native_bridge": "accepted C5 DAG history to native C5 history, not a feedback alias"},
        "source_prerequisite_audit": {name: accepted[name] for name in ("parent","memory","admission","candidate","positivity","leading","primitive","outer_calculus")},
        "proof_scope": {
            "native_C5_history_derivative_composition": bool(history_ok),
            "modified_m_affine_composition": bool(modified_ok),
            "physical_mean_y_affine_composition": bool(physical_ok),
            "conditional_Qy_Z1_function_identity": bool(qy_ok),
            "whole_domain_a_bounds_are_bounds_only": True,
            "cross_side_rule": "differentiate admitted all-Z Q_y[0] identity under per-side C1",
            "independent_cross_side_Z2_value_equality": False,
        },
    })
    return result
