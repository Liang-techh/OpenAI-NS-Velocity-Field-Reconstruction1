"""Build an unevaluated local DAG for the twenty unresolved frame-gap roles.

This draft deliberately has no dependency on the historical graph builders.
It reads the authenticated saved lookup packet and the archived JSON graph,
then copies only the missing gap and endpoint roots into a new JSON namespace.
Saved ``rewritten_nodes`` entries are opaque boundary references. They are
never followed, evaluated, or compared with P/F/G expressions.
"""
from __future__ import annotations

import base64
import copy
import gzip
import hashlib
import json
from pathlib import Path


ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
DEFERRED_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_historical_frame_gap_deferred_roles.json"
)
DEFERRED_REPORT_SHA256 = (
    "caca47450da722fa290d5cb8d29d845e2c747b2d81b311317966ab38a2e33e25"
)
DEFERRED_RESULT_SHA256 = (
    "3ba2e27cb66459ca75592de04bf6359fa3295733cf87aa1b02f4efcf62d30ec6"
)
LOOKUP_REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_historical_frame_gap_saved_rewrite_lookup.json"
LOOKUP_REPORT_SHA256 = "f9e6eff7457afbe591323268655be122455de98b14bf320a235c80137be75b00"
LOOKUP_RESULT_SHA256 = (
    "44da94d04c52e020799692049e8a6c58f1585a1ea4f794a0146ba933c71cfc4d"
)

SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"


# This is the explicit operand schema used by the frozen source importer,
# copied locally so this draft does not import or call that importer.
REFERENCE_KEYS = {
    "mathematical_positive_infinity": (),
    "exact_original_source_heat_phi": ("argument",),
    "exact_rational": (),
    "bound_variable": (),
    "mathematical_pi": (),
    "current_original_source_parameter": (),
    "original_source_parameter": (),
    "exact_positive_integer_expression": (),
    "shared_positive_integer": (),
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "exact_real_comparison": ("left", "right"),
    "logical_or": ("arguments",),
    "compact_raw_beta": ("argument",),
    "exact_original_source_flat_sigma": ("argument",),
    "current_original_leading_function_recipe": ("coordinate", "Z_variable"),
    "current_original_physical_source_scale": (),
    "current_original_physical_source_quantity": ("coordinate", "Z_variable"),
    "current_original_physical_radius": ("log_radius", "free_query"),
    "substitute_original_native_coordinate": ("body", "query_coordinate"),
    "substitute_free_query_at_native_coordinate": ("body", "free_query"),
    "substitute_original_inverse_angle": ("body", "inverse_angle"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "lower", "upper"),
    "original_lazy_flat_branch": (
        "active_body",
        "Delta",
        "eta",
        "flat_predicate",
        "original_collar_coordinate",
        "original_collar_upper",
    ),
    "original_monotone_phase_inverse": (
        "phase_function",
        "target_phase",
        "angle_lower",
        "angle_upper",
        "Delta",
        "eta",
    ),
}

ROLE_ORDER = (
    ("gap_m15", "gap_m15"),
    ("lower_R_q_y", "lower_R_q_y"),
    ("upper_R_q_y", "upper_R_q_y"),
    ("lower_q_plus_q_y", "lower_q_plus_q_y"),
    ("upper_q_plus_q_y", "upper_q_plus_q_y"),
)


def _read_hashed_bytes(root: Path, relative: str, expected: str) -> bytes:
    raw = (root / relative).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise ValueError("authenticated input hash mismatch: " + relative)
    return raw


def _read_lookup(root: Path) -> dict:
    report = json.loads(_read_hashed_bytes(root, LOOKUP_REPORT, LOOKUP_REPORT_SHA256))
    raw = base64.b64decode(report["artifacts"]["tmp/historical_frame_gap_saved_rewrite_lookup_v1/result.json"]["base64"], validate=True)
    if hashlib.sha256(raw).hexdigest() != LOOKUP_RESULT_SHA256:
        raise ValueError("embedded lookup result hash mismatch")
    packet = json.loads(raw)
    if packet.get("status") != "returned" or packet.get("consumer_calls") != 1:
        raise ValueError("saved lookup packet is not the single returned result")
    result = packet.get("result")
    if not isinstance(result, dict):
        raise ValueError("saved lookup packet has no result object")
    if result.get("schema") != "historical_frame_gap_saved_rewrite_lookup_v1":
        raise ValueError("saved lookup schema is not the expected leaf result")
    if result.get("status") != "saved_rewrite_map_lookups_not_numerical_values":
        raise ValueError("saved lookup status is not the expected direct-map status")
    return result


def _read_deferred_result(root: Path) -> dict:
    raw = _read_hashed_bytes(root, DEFERRED_REPORT, DEFERRED_REPORT_SHA256)
    report = json.loads(raw)
    artifact = (
        report.get("artifacts", {})
        .get("tmp/historical_frame_gap_deferred_roles_v1/result.json")
    )
    if not isinstance(artifact, dict) or not isinstance(artifact.get("base64"), str):
        raise ValueError("deferred roles report has no embedded result artifact")
    encoded = base64.b64decode(artifact["base64"], validate=True)
    if hashlib.sha256(encoded).hexdigest() != DEFERRED_RESULT_SHA256:
        raise ValueError("embedded deferred roles result hash mismatch")
    packet = json.loads(encoded)
    if packet.get("status") != "returned" or not isinstance(packet.get("result"), dict):
        raise ValueError("embedded deferred roles result is not returned")
    result = packet["result"]
    if result.get("schema") != "historical_frame_gap_deferred_roles_v1":
        raise ValueError("embedded deferred roles schema is not expected")
    return result


def _opaque_int(value, label: str) -> int:
    if type(value) is not int:
        raise ValueError(label + " must be an opaque integer reference")
    return value


def _assert_lookup_matches_deferred(deferred: dict, lookup: dict) -> None:
    if lookup.get("prior_result_sha256") != DEFERRED_RESULT_SHA256:
        raise ValueError("lookup packet is not pinned to the deferred roles result")
    roles = lookup.get("roles")
    if not isinstance(roles, dict):
        raise ValueError("lookup packet has no roles object")

    def lookup_refs(name: str) -> list[int]:
        rows = roles.get(name)
        if not isinstance(rows, list):
            raise ValueError("lookup role is not a row list: " + name)
        return [_opaque_int(item.get("original_ref"), name + ".original_ref") for item in rows]

    endpoint = deferred.get("original_gap_endpoint_roles")
    if not isinstance(endpoint, dict):
        raise ValueError("deferred result has no endpoint roles")
    expected = {
        "gap_m15": deferred.get("original_gap_m15_rows"),
        "lower_R_q_y": endpoint.get("lower", {}).get("R_q_y"),
        "upper_R_q_y": endpoint.get("upper", {}).get("R_q_y"),
        "lower_q_plus_q_y": endpoint.get("lower", {}).get("q_plus_q_y"),
        "upper_q_plus_q_y": endpoint.get("upper", {}).get("q_plus_q_y"),
    }
    for name, expected_rows in expected.items():
        if not isinstance(expected_rows, list):
            raise ValueError("deferred role is not a row list: " + name)
        expected_refs = [_opaque_int(item, name + ".deferred_ref") for item in expected_rows]
        if lookup_refs(name) != expected_refs:
            raise ValueError("lookup and deferred role references disagree: " + name)
    template = deferred.get("original_own_gap_field_template")
    radius = roles.get("field_radius")
    if not isinstance(template, dict) or not isinstance(radius, dict):
        raise ValueError("missing deferred radius role")
    if _opaque_int(template.get("radius"), "deferred.radius") != _opaque_int(
        radius.get("original_ref"), "lookup.radius.original_ref"
    ):
        raise ValueError("lookup and deferred radius references disagree")


def _read_archive(root: Path, lookup: dict) -> dict:
    if lookup.get("frame_archive_raw_sha256") != ARCHIVE_SHA256:
        raise ValueError("lookup packet is pinned to a different frame archive")
    raw = _read_hashed_bytes(root, ARCHIVE, ARCHIVE_SHA256)
    archive = json.loads(gzip.decompress(raw))
    if not isinstance(archive, dict):
        raise ValueError("frame archive is not a JSON object")
    return archive


def _ref(namespace: str, node: int | str) -> dict:
    return {"namespace": namespace, "node": node}


def _missing_root_rows(result: dict) -> dict:
    roles = result.get("roles")
    if not isinstance(roles, dict):
        raise ValueError("saved lookup has no roles object")
    radius = roles.get("field_radius")
    if (
        not isinstance(radius, dict)
        or radius.get("status") != "present_in_saved_rewrite_map"
        or _opaque_int(radius.get("original_ref"), "field_radius.original_ref")
        != _opaque_int(radius.get("rewritten_ref"), "field_radius.rewritten_ref")
    ):
        raise ValueError("field radius is not the saved identity boundary")
    if result.get("roles_present_count") != 1 or result.get("roles_missing_count") != 20:
        raise ValueError("saved lookup does not contain exactly one present and twenty missing roles")

    def rows_for(name: str) -> list[dict]:
        rows = roles.get(name)
        if not isinstance(rows, list) or len(rows) != 4:
            raise ValueError(name + " must contain four source rows")
        out = []
        for index, item in enumerate(rows):
            if not isinstance(item, dict):
                raise ValueError(name + " row is not an object")
            if item.get("status") != "not_in_saved_rewrite_map":
                raise ValueError(name + " contains a non-missing row")
            if item.get("rewritten_ref") is not None:
                raise ValueError(name + " missing row has a rewritten target")
            out.append(
                {
                    "row_order": index,
                    "original_ref": _opaque_int(item.get("original_ref"), name + ".original_ref"),
                }
            )
        return out

    missing = {name: rows_for(name) for name, _ in ROLE_ORDER}
    if sum(len(rows) for rows in missing.values()) != 20:
        raise ValueError("missing role scope is not exactly twenty rows")
    return {
        "gap_m15": missing["gap_m15"],
        "endpoints": {
            "lower": {
                "R_q_y": missing["lower_R_q_y"],
                "q_plus_q_y": missing["lower_q_plus_q_y"],
            },
            "upper": {
                "R_q_y": missing["upper_R_q_y"],
                "q_plus_q_y": missing["upper_q_plus_q_y"],
            },
        },
        "radius_boundary": {
            "original_ref": radius["original_ref"],
            "frame_ref": radius["rewritten_ref"],
            "same_opaque_id": True,
        },
    }


def _make_substitution(exact_graph_nodes: list, rewritten_nodes: dict, roots: dict) -> dict:
    if not isinstance(exact_graph_nodes, list) or not exact_graph_nodes:
        raise ValueError("archive exact_graph_nodes must be a nonempty list")
    if not isinstance(rewritten_nodes, dict) or not rewritten_nodes:
        raise ValueError("archive rewritten_nodes must be a nonempty mapping")

    local_nodes: list[dict] = []
    local_by_original: dict[int, str] = {}
    visiting: set[int] = set()
    root_context: dict[str, object] = {}

    def check_id(value, label):
        value = _opaque_int(value, label)
        if value < 0 or value >= len(exact_graph_nodes):
            raise ValueError(label + " outside authenticated archive graph: " + str(value))
        return value

    def saved_boundary(original: int) -> dict | None:
        check_id(original, "original ref")
        key = str(original)
        if key not in rewritten_nodes:
            return None
        target = check_id(rewritten_nodes[key], "saved target")
        return _ref(SAVED_NAMESPACE, target)

    def resolve(original: int) -> dict:
        original = check_id(original, "root")
        stack = [(original, False)]
        while stack:
            node, ready = stack.pop()
            boundary = saved_boundary(node)
            if boundary is not None or node in local_by_original:
                continue
            if ready:
                row = exact_graph_nodes[node]
                operation = row["operation"]
                rewritten_row = copy.deepcopy(row)
                def resolved_ref(child):
                    child = check_id(child, "operand")
                    saved = saved_boundary(child)
                    if saved is not None:
                        return saved
                    if child not in local_by_original:
                        raise ValueError("unresolved child in postorder: " + str(child))
                    return _ref(FRESH_NAMESPACE, local_by_original[child])
                for key in REFERENCE_KEYS[operation]:
                    if key in row and row[key] is not None:
                        value = row[key]
                        rewritten_row[key] = [resolved_ref(v) for v in value] if isinstance(value, list) else resolved_ref(value)
                visiting.remove(node)
                node_id = "n%08d" % len(local_nodes)
                local_by_original[node] = node_id
                local_nodes.append({"namespace": FRESH_NAMESPACE, "node_id": node_id,
                    "original_ref": node, "root_context": dict(root_context), "operand_schema": list(REFERENCE_KEYS[operation]),
                    "operator_row": rewritten_row})
                if len(local_nodes) > 100000:
                    raise ValueError("fresh DAG exceeds 100000 nodes")
                continue
            if node in visiting:
                raise ValueError("cycle in fresh source DAG at original ref " + str(node))
            row = exact_graph_nodes[node]
            if not isinstance(row, dict):
                raise ValueError("source node is not an object: " + str(node))
            operation = row.get("operation")
            if not isinstance(operation, str) or operation not in REFERENCE_KEYS:
                raise ValueError("unspecified graph operand schema at " + str(node) + ": " + repr(operation))
            visiting.add(node)
            stack.append((node, True))
            children = []
            for key in REFERENCE_KEYS[operation]:
                if key in row and row[key] is not None:
                    value = row[key]
                    children.extend(value if isinstance(value, list) else [value])
            for child in reversed(children):
                child = check_id(child, "operand")
                if saved_boundary(child) is None and child not in local_by_original:
                    stack.append((child, False))
        saved = saved_boundary(original)
        return saved if saved is not None else _ref(FRESH_NAMESPACE, local_by_original[original])

    def resolve_rows(rows: list[dict], role: str) -> list[dict]:
        out = []
        for item in rows:
            # Independent local memo for each role/order: no sharing or
            # comparison across source roots or P/F graph families.
            local_by_original.clear()
            visiting.clear()
            root_context.clear()
            root_context.update({"role": role, "row_order": item["row_order"]})
            if str(item["original_ref"]) in rewritten_nodes:
                raise ValueError("previously missing root is now present in saved map")
            out.append({"row_order": item["row_order"],
                "original_ref": item["original_ref"],
                "node_ref": resolve(item["original_ref"])})
        return out

    radius = roots["radius_boundary"]
    radius_ref = saved_boundary(radius["original_ref"])
    if radius_ref is None or radius_ref["node"] != radius["frame_ref"]:
        raise ValueError("radius saved boundary differs from prior lookup")
    resolved = {
        "field_radius": {"original_ref": radius["original_ref"], "node_ref": radius_ref},
        "gap_m15": resolve_rows(roots["gap_m15"], "gap_m15"),
        "endpoints": {
            direction: {
                field: resolve_rows(rows, direction + "." + field)
                for field, rows in fields.items()
            }
            for direction, fields in roots["endpoints"].items()
        },
    }
    return {
        "namespace": FRESH_NAMESPACE,
        "existing_frame_namespace": SAVED_NAMESPACE,
        "roots": resolved,
        "nodes": local_nodes,
        "fresh_node_count": len(local_nodes),
        "missing_root_count": 20,
        "graph_domain_contract": "original and saved target refs are indices in the same authenticated Rp archive exact_graph_nodes",
        "fresh_memo_scope": "one independent root at a time",
        "boundary_cache": {
            "namespace": SAVED_NAMESPACE,
            "owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
            "semantics": "opaque existing frame node references; not followed",
        },
    }


def recover(root: Path) -> dict:
    """Return a pure JSON, unevaluated local substitution packet."""
    deferred = _read_deferred_result(root)
    lookup = _read_lookup(root)
    _assert_lookup_matches_deferred(deferred, lookup)
    archive = _read_archive(root, lookup)
    frame = archive.get("Rp_native_C5_pulse_source_bindings")
    if not isinstance(frame, dict):
        raise ValueError("frame archive is missing Rp_native_C5_pulse_source_bindings")
    if frame.get("input_substitution") != deferred.get("archive_input_substitution"):
        raise ValueError("saved input substitution differs from prior source packet")
    nodes = archive.get("exact_graph_nodes")
    mapping = frame.get("rewritten_nodes")
    bindings = frame.get("input_substitution")
    if not isinstance(nodes, list) or not isinstance(mapping, dict) or not isinstance(bindings, dict) or not bindings:
        raise ValueError("missing authenticated graph, saved mapping or inlet bindings")
    for key, target in bindings.items():
        if not isinstance(key, str) or not key.isdecimal() or str(int(key)) != key or int(key) in (0, 1):
            raise ValueError("invalid saved inlet binding key")
        if type(target) is not int or not (0 <= target < len(nodes)) or not (0 <= int(key) < len(nodes)):
            raise ValueError("saved inlet binding outside archive graph")
        if mapping.get(key) != target or type(mapping.get(key)) is not int:
            raise ValueError("saved inlet binding does not agree with saved rewrite mapping")
    if deferred.get("archived_receipts_passed") is not True or deferred.get("same_source_family_metadata") is not True:
        raise ValueError("prior archive provenance gates are not recorded")
    roots = _missing_root_rows(lookup)
    dag = _make_substitution(
        archive.get("exact_graph_nodes"),
        frame.get("rewritten_nodes"),
        roots,
    )
    dag["schema"] = "frame_gap_fresh_substitution_v1"
    dag["status"] = "unevaluated_local_source_dag"
    dag["deferred_result_sha256"] = DEFERRED_RESULT_SHA256
    dag["lookup_result_sha256"] = LOOKUP_RESULT_SHA256
    dag["frame_archive_raw_sha256"] = ARCHIVE_SHA256
    dag["radius_boundary"] = roots["radius_boundary"]
    dag["source_values_computed"] = False
    dag["integrals_evaluated"] = False
    dag["old_builders_called"] = False
    dag["old_error_terms_removed"] = False
    dag["extra_ledger_terms_added"] = False
    dag["P_equals_F_proved"] = False
    dag["G_equals_F_proved"] = False
    dag["admissible_for_numeric_use"] = False
    dag["full_goal_complete"] = False
    return dag
