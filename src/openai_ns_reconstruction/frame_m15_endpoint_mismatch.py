"""Unexecuted source-bound m15 endpoint seam expressions.

This leaf reads the authenticated saved Rp archive and the already-published
fresh frame-gap receipt only when :func:`recover` is explicitly called.  It
does not import a project module, construct a graph, follow graph edges,
evaluate a source value, evaluate an integral, or assert an endpoint identity.

The packet keeps the eight endpoint records (main/gap/end and global, lower
and upper) and the four finite m15 boundary fields at ordinary Z orders zero
and one.  Missing saved-map entries remain an opaque original reference plus
an unexecuted saved-inlet-substitution descriptor.  The two gap ``R_q_y``
rows may reuse their previously accepted fresh references, with their original
source references and unexecuted substitution provenance retained.

The local-minus-global expression is conditional on additive integral
blocks, algebra-preserving simultaneous substitution, and commutation of the
ordinary first-Z row with substitution.  No equality, error reduction, or
numeric use is claimed by this draft.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


ARCHIVE_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
ARCHIVE_RAW_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

C3_ARCHIVE_NAME = "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
C3_ARCHIVE_SHA256 = "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"

FRESH_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution.json"
)
FRESH_REPORT_SHA256 = "780a667e601422a1dfb3b13864f8ffa2e8ea12f8c2e875183f025e145d95ff95"
FRESH_RESULT_SHA256 = "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
FRESH_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_v1/result.json"
AUDIT_REPORT_RELATIVE = Path("experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution_audit.json")
AUDIT_REPORT_SHA256 = "de6de458582dab9b77eee5ed81ef6d1669bff245416aa118536d3fadb68f5432"
AUDIT_RESULT_SHA256 = "a3de519abfdc42fe18a0fa6830105c9c2c39a7bb4188d27693123117123087a3"
AUDIT_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_audit_v1/result.json"

C3_INPUT_KEY = C3_ARCHIVE_NAME
SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
DEFERRED_NAMESPACE = "unexecuted_saved_inlet_substitution_v1"

M15_ENDPOINT_FIELDS = (
    "R_squared_q",
    "R_squared_q_squared",
    "R_squared_U_squared",
    "R_q_y",
)
ORDINARY_Z_ROWS = (0, 1)

EXPECTED_GAP_R_Q_Y = {
    "lower": (34762, 34763),
    "upper": (34791, 34792),
}

ENDPOINT_LAYOUT = (
    ("main.lower", "chart", "main", 0),
    ("main.upper", "chart", "main", 1),
    ("gap.lower", "chart", "gap", 0),
    ("gap.upper", "chart", "gap", 1),
    ("end.lower", "chart", "end", 0),
    ("end.upper", "chart", "end", 1),
    ("global.lower", "global", "global", 0),
    ("global.upper", "global", "global", 1),
)

SEAM_TERMS = (
    (1, "global.lower"),
    (-1, "main.lower"),
    (1, "main.upper"),
    (-1, "gap.lower"),
    (1, "gap.upper"),
    (-1, "end.lower"),
    (1, "end.upper"),
    (-1, "global.upper"),
)


def _read_hashed_bytes(root: Path, relative: Path, expected: str) -> bytes:
    raw = (root / relative).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise ValueError("authenticated input hash mismatch: " + str(relative))
    return raw


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(label + " must be an object")
    return value


def _rows(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(label + " must contain at least Z0 and Z1 rows")
    return value


def _source_ref(
    original_ref: Any,
    rewritten_nodes: Mapping[str, Any],
    source_count: int,
    substitution_digest: str,
    *,
    fresh_ref: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Represent an opaque source row without constructing its expression."""

    if not _is_int(original_ref) or not (0 <= original_ref < source_count):
        raise ValueError("source row is not an in-range opaque archive index")
    key = str(original_ref)
    if key in rewritten_nodes:
        target = rewritten_nodes[key]
        if not _is_int(target) or not (0 <= target < source_count):
            raise ValueError("saved rewritten target is not an in-range archive index")
        return {
            "namespace": SAVED_NAMESPACE,
            "node": target,
            "original_ref": original_ref,
            "status": "saved_frame_reference",
        }

    substitution = {
        "namespace": DEFERRED_NAMESPACE,
        "operation": "simultaneous_saved_inlet_substitution",
        "mapping_pointer": "$.Rp_native_C5_pulse_source_bindings.input_substitution",
        "mapping_sha256": substitution_digest,
        "executed": False,
    }
    if fresh_ref is not None:
        if fresh_ref.get("namespace") != FRESH_NAMESPACE or not isinstance(
            fresh_ref.get("node"), str
        ):
            raise ValueError("reused fresh reference has the wrong shape")
        return {
            "namespace": FRESH_NAMESPACE,
            "node": fresh_ref["node"],
            "original_ref": original_ref,
            "status": "reused_committed_gap_R_q_y_reference",
            "substitution_provenance": {
                "operation": "simultaneous_saved_inlet_substitution",
                "mapping_pointer": substitution["mapping_pointer"],
                "mapping_sha256": substitution_digest,
                "fresh_source_DAG_constructed": True,
                "source_values_evaluated": False,
            },
            "fresh_report_sha256": FRESH_REPORT_SHA256,
            "fresh_result_sha256": FRESH_RESULT_SHA256,
        }
    return {
        "namespace": DEFERRED_NAMESPACE,
        "original_ref": original_ref,
        "status": "missing_saved_reference",
        "saved_inlet_substitution": substitution,
    }


def _load_fresh_gap_R_q_y(root: Path) -> dict[tuple[str, int], dict[str, Any]]:
    """Read only the two direct gap root lists from the frozen fresh report."""

    audit_report = json.loads(_read_hashed_bytes(root, AUDIT_REPORT_RELATIVE, AUDIT_REPORT_SHA256))
    audit_raw = base64.b64decode(audit_report["artifacts"][AUDIT_RESULT_ARTIFACT]["base64"], validate=True)
    if hashlib.sha256(audit_raw).hexdigest() != AUDIT_RESULT_SHA256:
        raise ValueError("fresh structural audit result hash mismatch")
    audit = json.loads(audit_raw)
    if audit.get("status") != "returned" or audit.get("result", {}).get("status") != "passed":
        raise ValueError("fresh structural audit is not passed")
    report_raw = _read_hashed_bytes(
        root, FRESH_REPORT_RELATIVE, FRESH_REPORT_SHA256
    )
    report = json.loads(report_raw.decode("utf-8"))
    artifact = _mapping(report.get("artifacts"), "fresh report artifacts").get(
        FRESH_RESULT_ARTIFACT
    )
    artifact = _mapping(artifact, "fresh result artifact")
    encoded = artifact.get("base64")
    if not isinstance(encoded, str):
        raise ValueError("fresh result artifact has no base64 payload")
    raw = base64.b64decode(encoded, validate=True)
    if hashlib.sha256(raw).hexdigest() != FRESH_RESULT_SHA256:
        raise ValueError("fresh result raw hash mismatch")
    packet = json.loads(raw.decode("utf-8"))
    if packet.get("status") != "returned":
        raise ValueError("fresh result is not a returned packet")
    body = _mapping(packet.get("result"), "fresh result body")
    if body.get("schema") != "frame_gap_fresh_substitution_v1":
        raise ValueError("fresh result schema mismatch")
    roots = _mapping(body.get("roots"), "fresh result roots")
    endpoints = _mapping(roots.get("endpoints"), "fresh result endpoints")

    out: dict[tuple[str, int], dict[str, Any]] = {}
    for side in ("lower", "upper"):
        expected = EXPECTED_GAP_R_Q_Y[side]
        side_rows = _mapping(endpoints.get(side), "fresh gap " + side)
        entries = side_rows.get("R_q_y")
        if not isinstance(entries, list) or len(entries) != 4:
            raise ValueError("fresh gap " + side + " R_q_y has wrong row count")
        for order in ORDINARY_Z_ROWS:
            entry = _mapping(entries[order], "fresh gap R_q_y row")
            if entry.get("row_order") != order or entry.get("original_ref") != expected[order]:
                raise ValueError("fresh gap R_q_y source row mismatch")
            node_ref = _mapping(entry.get("node_ref"), "fresh gap node_ref")
            if node_ref.get("namespace") != FRESH_NAMESPACE or not isinstance(
                node_ref.get("node"), str
            ):
                raise ValueError("fresh gap R_q_y node_ref is not fresh")
            out[(side, order)] = {
                "namespace": FRESH_NAMESPACE,
                "node": node_ref["node"],
                "original_ref": expected[order],
            }
    return out


def _load_archive(root: Path) -> tuple[dict[str, Any], bytes]:
    raw = _read_hashed_bytes(root, ARCHIVE_RELATIVE, ARCHIVE_RAW_SHA256)
    decoded = gzip.decompress(raw)
    archive = json.loads(decoded.decode("utf-8"))
    return _mapping(archive, "Rp archive"), raw


def _endpoint_packet(
    label: str,
    endpoint: Mapping[str, Any],
    rewritten_nodes: Mapping[str, Any],
    source_count: int,
    substitution_digest: str,
    fresh_gap_R_q_y: Mapping[tuple[str, int], Mapping[str, Any]],
) -> dict[str, Any]:
    fields: dict[str, Any] = {}
    side = label.split(".", 1)[1] if label.startswith("gap.") else None
    for field in M15_ENDPOINT_FIELDS:
        raw_rows = _rows(endpoint.get(field), label + "." + field)
        row_bindings: dict[str, Any] = {}
        original_refs: list[int] = []
        for order in ORDINARY_Z_ROWS:
            original = raw_rows[order]
            if not _is_int(original):
                raise ValueError(label + "." + field + " has a non-integer source row")
            original_refs.append(original)
            fresh = None
            if label == "gap." + str(side) and field == "R_q_y":
                fresh = fresh_gap_R_q_y.get((str(side), order))
                if str(original) != str(EXPECTED_GAP_R_Q_Y[str(side)][order]):
                    raise ValueError("archive and fresh gap R_q_y rows disagree")
            row_bindings["Z" + str(order)] = _source_ref(
                original,
                rewritten_nodes,
                source_count,
                substitution_digest,
                fresh_ref=fresh,
            )
        fields[field] = {
            "original_refs_Z0_Z1": original_refs,
            "rows": row_bindings,
        }
    return {"endpoint": label, "fields": fields}


def _formula_packet(endpoint: Mapping[str, Any], delta_ref: Mapping[str, Any]) -> dict[str, Any]:
    """Return exact rational forms; A/B/C/D are bound below to source rows."""

    fields = _mapping(endpoint["fields"], endpoint["endpoint"] + ".fields")
    aliases = {
        "A0": fields["R_squared_q"]["rows"]["Z0"],
        "A1": fields["R_squared_q"]["rows"]["Z1"],
        "B0": fields["R_squared_q_squared"]["rows"]["Z0"],
        "B1": fields["R_squared_q_squared"]["rows"]["Z1"],
        "C0": fields["R_squared_U_squared"]["rows"]["Z0"],
        "C1": fields["R_squared_U_squared"]["rows"]["Z1"],
        "D0": fields["R_q_y"]["rows"]["Z0"],
        "D1": fields["R_q_y"]["rows"]["Z1"],
    }
    return {
        "endpoint": endpoint["endpoint"],
        "delta_source_ref": dict(delta_ref),
        "formal_symbol": "Z",
        "weights": {
            "L": "1-delta*Z^2",
            "L_Z": "-2*delta*Z",
        },
        "aliases": aliases,
        "Z0": {
            "K": "A0/(2*L) + B0/4 - Z^2*C0/L^2 - D0",
            "note": "B0/4 is intentionally not divided by L",
        },
        "Z1": {
            "K_Z": (
                "A1/(2*L) - A0*L_Z/(2*L^2) + B1/4 "
                "- (2*Z*C0 + Z^2*C1)/L^2 "
                "+ 2*Z^2*C0*L_Z/L^3 - D1"
            ),
            "expanded_L_Z": "L_Z=-2*delta*Z",
            "note": "Z1 includes both L_Z and the derivative of Z^2*C0",
        },
        "unvalued": True,
        "integrals_evaluated": False,
    }


def _seam_expression(component: str) -> dict[str, Any]:
    labels = [
        ("+" if coefficient > 0 else "-") + "K" + component + "(" + label + ")"
        for coefficient, label in SEAM_TERMS
    ]
    return {
        "component": component,
        "terms": [
            {"coefficient": coefficient, "endpoint": label, "component": "K" + component}
            for coefficient, label in SEAM_TERMS
        ],
        "rational_form": " ".join(labels).replace("+ -", "- "),
        "zero_claimed": False,
    }


def recover(root: Path) -> dict[str, Any]:
    """Return a source-reference packet; this function is intentionally opt-in."""

    archive, archive_raw = _load_archive(root)
    fresh_gap_R_q_y = _load_fresh_gap_R_q_y(root)
    frame = _mapping(
        archive.get("Rp_native_C5_pulse_source_bindings"),
        "Rp frame bindings",
    )
    input_hashes = _mapping(archive.get("input_hashes"), "Rp archive input hashes")
    if input_hashes.get(C3_ARCHIVE_NAME) != C3_ARCHIVE_SHA256:
        raise ValueError("Rp archive does not pin the C3 input archive")
    inputs = _mapping(frame.get("inputs"), "Rp frame inputs")
    c3 = _mapping(inputs.get(C3_INPUT_KEY), "saved C3 pulse input")
    endpoints = _mapping(c3.get("endpoints"), "saved C3 endpoints")
    global_endpoints = c3.get("global_endpoints")
    if not isinstance(global_endpoints, list) or len(global_endpoints) != 2:
        raise ValueError("saved C3 global_endpoints must contain lower and upper")

    rewritten_nodes = _mapping(frame.get("rewritten_nodes"), "saved rewritten_nodes")
    substitution = _mapping(frame.get("input_substitution"), "saved input_substitution")
    source_nodes = archive.get("exact_graph_nodes")
    if not isinstance(source_nodes, list) or not source_nodes:
        raise ValueError("saved exact_graph_nodes must be a nonempty list")
    source_count = len(source_nodes)
    substitution_digest = hashlib.sha256(
        json.dumps(substitution, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    geometry = _mapping(frame.get("geometry"), "saved frame geometry")
    delta = geometry.get("delta")
    if not _is_int(delta):
        raise ValueError("saved frame geometry delta must be an opaque source ref")
    if not (0 <= delta < source_count):
        raise ValueError("saved frame geometry delta is outside archive graph")
    delta_ref = {"namespace": SAVED_NAMESPACE, "node": delta,
        "status": "direct_saved_frame_geometry_reference"}

    endpoint_packets: dict[str, Any] = {}
    for label, kind, container, index in ENDPOINT_LAYOUT:
        if kind == "global":
            endpoint = _mapping(global_endpoints[index], label)
        else:
            chart = endpoints.get(container)
            if isinstance(chart, list) and len(chart) == 2:
                endpoint = _mapping(chart[index], label)
            elif isinstance(chart, Mapping):
                chart_rows = chart.get("lower" if index == 0 else "upper")
                endpoint = _mapping(chart_rows, label)
            else:
                raise ValueError(container + " endpoint list is missing")
        endpoint_packets[label] = _endpoint_packet(
            label,
            endpoint,
            rewritten_nodes,
            source_count,
            substitution_digest,
            fresh_gap_R_q_y,
        )

    formulas = {
        label: _formula_packet(packet, delta_ref)
        for label, packet in endpoint_packets.items()
    }
    row_count = sum(
        len(packet["fields"][field]["rows"])
        for packet in endpoint_packets.values()
        for field in M15_ENDPOINT_FIELDS
    )
    if len(endpoint_packets) != 8 or row_count != 64:
        raise ValueError("endpoint packet is not exactly eight endpoints and 64 rows")
    status_counts: dict[str, int] = {}
    for packet in endpoint_packets.values():
        for field in M15_ENDPOINT_FIELDS:
            for binding in packet["fields"][field]["rows"].values():
                status = binding["status"]
                status_counts[status] = status_counts.get(status, 0) + 1

    return {
        "schema": "frame_m15_endpoint_mismatch_v1",
        "status": "unexecuted_source_bound_endpoint_mismatch",
        "endpoint_layout": [label for label, _, _, _ in ENDPOINT_LAYOUT],
        "endpoint_fields": list(M15_ENDPOINT_FIELDS),
        "ordinary_Z_rows": list(ORDINARY_Z_ROWS),
        "endpoint_count": len(endpoint_packets),
        "source_row_count": row_count,
        "binding_status_counts": status_counts,
        "source_references_are_opaque": True,
        "saved_frame_namespace": SAVED_NAMESPACE,
        "fresh_gap_namespace": FRESH_NAMESPACE,
        "unexecuted_substitution_namespace": DEFERRED_NAMESPACE,
        "provenance": {
            "archive_relative": str(ARCHIVE_RELATIVE),
            "archive_raw_sha256": hashlib.sha256(archive_raw).hexdigest(),
            "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
            "C3_archive_name": C3_ARCHIVE_NAME,
            "C3_archive_sha256": C3_ARCHIVE_SHA256,
            "fresh_report_relative": str(FRESH_REPORT_RELATIVE),
            "fresh_report_sha256": FRESH_REPORT_SHA256,
            "fresh_result_sha256": FRESH_RESULT_SHA256,
            "fresh_structural_audit_report_sha256": AUDIT_REPORT_SHA256,
            "fresh_structural_audit_result_sha256": AUDIT_RESULT_SHA256,
            "saved_rewritten_nodes_count": len(rewritten_nodes),
            "saved_input_substitution_sha256": substitution_digest,
        },
        "delta_source_ref": delta_ref,
        "endpoints": endpoint_packets,
        "unvalued_K_forms": formulas,
        "conditional_identity": {
            "name": "sum_local_m15_minus_global_m15",
            "rational_form": (
                "K(global.lower)-K(main.lower)+K(main.upper)-K(gap.lower)+"
                "K(gap.upper)-K(end.lower)+K(end.upper)-K(global.upper)"
            ),
            "Z1_rational_form": (
                "K_Z(global.lower)-K_Z(main.lower)+K_Z(main.upper)-K_Z(gap.lower)+"
                "K_Z(gap.upper)-K_Z(end.lower)+K_Z(end.upper)-K_Z(global.upper)"
            ),
            "Z0_terms": _seam_expression(""),
            "Z1_terms": _seam_expression("_Z"),
            "requires": [
                "main/gap/end integral blocks add under the same source substitution",
                "simultaneous saved-inlet substitution preserves the source algebra",
                "ordinary Z1 differentiation commutes with that substitution",
                "delta is independent of Z and L=1-delta*Z^2 is nonzero",
            ],
            "zero_requires": ["the four displayed seam contributions vanish by separately established identities"],
            "endpoint_equalities_proved": False,
            "residual_zero_proved": False,
        },
        "source_values_computed": False,
        "integrals_evaluated": False,
        "numeric_values_computed": False,
        "old_builders_called": False,
        "old_operators_replayed": False,
        "old_error_terms_removed": False,
        "extra_ledger_terms_added": False,
        "admissible_for_numeric_use": False,
        "full_goal_complete": False,
    }

