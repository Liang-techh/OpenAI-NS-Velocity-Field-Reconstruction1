"""Independent structural audit for the deferred m15 source DAG receipt.

This module is intentionally independent of the project constructors and
checkers.  ``recover`` is opt-in: it authenticates the published receipt and
the source archive, then checks the new 168-node packet by direct indexed
references.  It does not evaluate values, integrals, P/F expressions, or
physical/semantic identities.

The prior 1131-node frame-gap DAG is an authenticated opaque boundary.  This
audit checks only its signed report/result status, hashes, root references,
and an ``original_ref`` index.  It deliberately does not re-audit prior-node
operator rows or operands.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


AUDIT_SCHEMA = "frame_m15_deferred_endpoint_substitution_audit_v1"

RESULT_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_m15_deferred_endpoint_substitution.json"
)
RESULT_REPORT_SHA256 = (
    "65be948ec715056f310e17c827b3132f37f06b6439d932743cfb3b57e2c50d39"
)
RESULT_ARTIFACT = "tmp/frame_m15_deferred_endpoint_substitution_v1/result.json"
RESULT_SHA256 = (
    "bf07cdce2ba58cd9beae53fc576a3ff4d023f56e4a027de495f194cd0047d5f9"
)

ARCHIVE_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
ARCHIVE_RAW_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)

MISMATCH_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_m15_endpoint_mismatch.json"
)
MISMATCH_REPORT_SHA256 = (
    "cb8d8aa0cf5a877b03bcda35f7cae32c6761171d9b739e4035ca2302f938d884"
)
MISMATCH_RESULT_ARTIFACT = "tmp/frame_m15_endpoint_mismatch_v1/result.json"
MISMATCH_RESULT_SHA256 = (
    "6adc3c1a272b2fcc6f6259c8ad844b8560a9893ac4efdcd8dc689e7bcd177d93"
)

FRESH_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution.json"
)
FRESH_REPORT_SHA256 = (
    "780a667e601422a1dfb3b13864f8ffa2e8ea12f8c2e875183f025e145d95ff95"
)
FRESH_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_v1/result.json"
FRESH_RESULT_SHA256 = (
    "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
)

AUDIT_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution_audit.json"
)
AUDIT_REPORT_SHA256 = (
    "de6de458582dab9b77eee5ed81ef6d1669bff245416aa118536d3fadb68f5432"
)
AUDIT_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_audit_v1/result.json"
AUDIT_RESULT_SHA256 = (
    "a3de519abfdc42fe18a0fa6830105c9c2c39a7bb4188d27693123117123087a3"
)

C3_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
)
C3_ARCHIVE_SHA256 = (
    "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"
)

SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
PRIOR_FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
FRESH_NAMESPACE = "fresh_frame_m15_deferred_endpoint_substitution_v1"
DEFERRED_NAMESPACE = "unexecuted_saved_inlet_substitution_v1"

M15_ENDPOINT_FIELDS = (
    "R_squared_q",
    "R_squared_q_squared",
    "R_squared_U_squared",
    "R_q_y",
)
ORDINARY_Z_ROWS = (0, 1)
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
EXPECTED_GAP_R_Q_Y = {
    "lower": (34762, 34763),
    "upper": (34791, 34792),
}

EXPECTED_DEFERRED_OCCURRENCES = 44
EXPECTED_RETAINED_OCCURRENCES = 20
EXPECTED_NEW_NODE_COUNT = 168
EXPECTED_NEW_ROOT_COUNT = 32
EXPECTED_PRIOR_ROOT_COUNT = 12
EXPECTED_SAVED_RETAINED_COUNT = 16
EXPECTED_PRIOR_RETAINED_COUNT = 4

# This is a frozen source-row operand schema.  It is deliberately local and
# is not imported from the constructor or from the prior audit.
FROZEN_OPERAND_SCHEMA: dict[str, tuple[str, ...]] = {
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
LIST_OPERAND_FIELDS = frozenset({"arguments"})


class _AuditAbort(Exception):
    """Internal fail-closed stop after an authenticated input is malformed."""


class _Audit:
    def __init__(self) -> None:
        self.failures: list[dict[str, Any]] = []
        self.failure_count = 0
        self.counts: dict[str, int] = {
            "deferred_occurrences_checked": 0,
            "retained_occurrences_checked": 0,
            "new_nodes_checked": 0,
            "new_operand_refs_checked": 0,
            "saved_boundary_refs_checked": 0,
            "prior_boundary_refs_checked": 0,
            "new_boundary_refs_checked": 0,
            "source_rows_checked": 0,
            "root_context_checks": 0,
            "postorder_checks": 0,
            "prior_index_entries": 0,
        }

    def fail(self, code: str, detail: str, **context: Any) -> None:
        self.failure_count += 1
        if len(self.failures) >= 100:
            return
        item: dict[str, Any] = {"code": code, "detail": detail}
        if context:
            item["context"] = context
        self.failures.append(item)

    def abort(self, code: str, detail: str, **context: Any) -> None:
        self.fail(code, detail, **context)
        raise _AuditAbort(detail)

    def report(
        self,
        *,
        result_digest: str | None = None,
        archive_digest: str | None = None,
    ) -> dict[str, Any]:
        flags = {
            "semantic_identity_proved": False,
            "ordinary_Z1_commutation_proved": False,
            "Z1_commutation_proved": False,
            "seam_identity_proved": False,
            "endpoint_equalities_proved": False,
            "P_equals_F_proved": False,
            "G_equals_F_proved": False,
            "physical_equivalence_proved": False,
            "physical_equalities_authenticated": False,
            "source_values_computed": False,
            "numeric_values_computed": False,
            "integrals_evaluated": False,
            "admissible_for_numeric_use": False,
            "old_builders_called": False,
            "old_constructors_called": False,
            "old_operators_replayed": False,
            "old_error_terms_removed": False,
            "extra_ledger_terms_added": False,
            "old_twenty_roots_replayed": False,
            "full_goal_complete": False,
        }
        return {
            "schema": AUDIT_SCHEMA,
            "status": "passed" if self.failure_count == 0 else "failed",
            "result_sha256": result_digest,
            "expected_result_sha256": RESULT_SHA256,
            "archive_raw_sha256": archive_digest,
            "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
            "counts": dict(self.counts),
            "failure_count": self.failure_count,
            "failure_details": list(self.failures),
            "failure_details_truncated": self.failure_count > len(self.failures),
            "semantic_flags": dict(flags),
            **flags,
        }


def _is_int(value: Any) -> bool:
    return type(value) is int


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise _AuditAbort(label + " must be an object")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise _AuditAbort(label + " must be a list")
    return value


def _index(value: Any, label: str, source_count: int) -> int:
    if not _is_int(value) or value < 0 or value >= source_count:
        raise _AuditAbort(label + " is outside the authenticated source graph")
    return value


def _read_hashed_bytes(root: Path, relative: Path, expected: str) -> tuple[bytes, str]:
    raw = (root / relative).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected:
        raise _AuditAbort(
            "hash mismatch for %s: %s != %s" % (relative, digest, expected)
        )
    return raw, digest


def _load_json_result(
    root: Path,
    report_relative: Path,
    report_sha256: str,
    artifact_name: str,
    result_sha256: str,
    label: str,
) -> tuple[dict[str, Any], str]:
    report_raw, _ = _read_hashed_bytes(root, report_relative, report_sha256)
    try:
        report = json.loads(report_raw.decode("utf-8"))
    except Exception as exc:
        raise _AuditAbort(label + " report JSON is invalid") from exc
    report_obj = _mapping(report, label + " report")
    artifacts = _mapping(report_obj.get("artifacts"), label + " artifacts")
    artifact = _mapping(artifacts.get(artifact_name), label + " artifact")
    encoded = artifact.get("base64")
    if not isinstance(encoded, str):
        raise _AuditAbort(label + " artifact has no base64 payload")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except Exception as exc:
        raise _AuditAbort(label + " artifact base64 is invalid") from exc
    digest = hashlib.sha256(raw).hexdigest()
    if digest != result_sha256:
        raise _AuditAbort(label + " result hash mismatch")
    metadata_digest = artifact.get("sha256")
    if metadata_digest not in (None, result_sha256):
        raise _AuditAbort(label + " artifact metadata hash mismatch")
    try:
        packet = json.loads(raw.decode("utf-8"))
    except Exception as exc:
        raise _AuditAbort(label + " result JSON is invalid") from exc
    packet_obj = _mapping(packet, label + " packet")
    if packet_obj.get("status") != "returned":
        raise _AuditAbort(label + " packet is not returned")
    body = _mapping(packet_obj.get("result"), label + " result")
    return dict(body), digest


def _load_archive(root: Path) -> tuple[dict[str, Any], str]:
    raw, digest = _read_hashed_bytes(root, ARCHIVE_RELATIVE, ARCHIVE_RAW_SHA256)
    try:
        archive = json.loads(gzip.decompress(raw).decode("utf-8"))
    except Exception as exc:
        raise _AuditAbort("authenticated Rp archive cannot be decoded") from exc
    return dict(_mapping(archive, "Rp archive")), digest


def _canonical_mapping(
    value: Any, label: str, source_count: int
) -> dict[str, int]:
    mapping = _mapping(value, label)
    if not mapping:
        raise _AuditAbort(label + " must be nonempty")
    output: dict[str, int] = {}
    for key, target in mapping.items():
        if not isinstance(key, str) or not key.isdecimal() or str(int(key)) != key:
            raise _AuditAbort(label + " contains a noncanonical source index")
        if key in output:
            raise _AuditAbort(label + " contains a duplicate source index")
        output[key] = _index(target, label + "[" + key + "]", source_count)
    return output


def _substitution_digest(value: Mapping[str, Any]) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _endpoint_from_c3(c3: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    for candidate, kind, container, index in ENDPOINT_LAYOUT:
        if candidate != label:
            continue
        if kind == "global":
            rows = c3.get("global_endpoints")
            if not isinstance(rows, list) or len(rows) != 2:
                raise _AuditAbort("saved C3 global_endpoints has wrong shape")
            return _mapping(rows[index], label + " C3 endpoint")
        charts = c3.get("endpoints")
        if isinstance(charts, list):
            if len(charts) != 3:
                raise _AuditAbort("saved C3 endpoints list has wrong chart count")
            chart = charts[("main", "gap", "end").index(container)]
        elif isinstance(charts, Mapping):
            chart = charts.get(container)
        else:
            raise _AuditAbort("saved C3 endpoints are missing")
        if isinstance(chart, list):
            if len(chart) != 2:
                raise _AuditAbort(label + " C3 endpoint list has wrong length")
            return _mapping(chart[index], label + " C3 endpoint")
        if isinstance(chart, Mapping):
            return _mapping(
                chart.get("lower" if index == 0 else "upper"),
                label + " C3 endpoint",
            )
        raise _AuditAbort(label + " C3 endpoint is missing")
    raise _AuditAbort("unknown endpoint label: " + label)


def _load_source_context(
    archive: Mapping[str, Any], audit: _Audit
) -> tuple[list[Any], dict[str, int], dict[str, int], str, Mapping[str, Any]]:
    source_nodes = _list(archive.get("exact_graph_nodes"), "Rp exact_graph_nodes")
    if not source_nodes:
        raise _AuditAbort("Rp exact_graph_nodes is empty")
    source_count = len(source_nodes)
    frame = _mapping(
        archive.get("Rp_native_C5_pulse_source_bindings"),
        "Rp frame bindings",
    )
    substitution = _canonical_mapping(
        frame.get("input_substitution"),
        "saved input_substitution",
        source_count,
    )
    rewritten = _canonical_mapping(
        frame.get("rewritten_nodes"),
        "saved rewritten_nodes",
        source_count,
    )
    if len(rewritten) < len(substitution):
        raise _AuditAbort("saved rewritten_nodes is smaller than input_substitution")
    for key, target in substitution.items():
        if key in ("0", "1"):
            raise _AuditAbort("saved input_substitution rewrites reserved root")
        if rewritten.get(key) != target:
            raise _AuditAbort("saved input_substitution disagrees with rewritten_nodes")
    input_hashes = _mapping(archive.get("input_hashes"), "Rp input hashes")
    if input_hashes.get(C3_ARCHIVE_NAME) != C3_ARCHIVE_SHA256:
        raise _AuditAbort("Rp archive does not pin the C3 input archive")
    c3 = _mapping(
        _mapping(frame.get("inputs"), "Rp frame inputs").get(C3_ARCHIVE_NAME),
        "saved C3 pulse input",
    )
    digest = _substitution_digest(_mapping(frame.get("input_substitution"), "input_substitution"))
    audit.counts["source_rows_checked"] = 0
    return source_nodes, rewritten, substitution, digest, c3


def _expected_scope(
    mismatch: Mapping[str, Any],
    c3: Mapping[str, Any],
    source_count: int,
    rewritten: Mapping[str, int],
    substitution_digest: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if mismatch.get("schema") != "frame_m15_endpoint_mismatch_v1":
        raise _AuditAbort("mismatch schema mismatch")
    if mismatch.get("status") != "unexecuted_source_bound_endpoint_mismatch":
        raise _AuditAbort("mismatch status mismatch")
    if mismatch.get("source_references_are_opaque") is not True:
        raise _AuditAbort("mismatch result is not source-opaque")
    provenance = _mapping(mismatch.get("provenance"), "mismatch provenance")
    expected_provenance = {
        "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
        "archive_raw_sha256": ARCHIVE_RAW_SHA256,
        "C3_archive_name": C3_ARCHIVE_NAME,
        "C3_archive_sha256": C3_ARCHIVE_SHA256,
        "saved_input_substitution_sha256": substitution_digest,
    }
    for key, expected in expected_provenance.items():
        if provenance.get(key) != expected:
            raise _AuditAbort("mismatch provenance disagrees: " + key)
    counts = _mapping(mismatch.get("binding_status_counts"), "mismatch counts")
    if dict(counts) != {
        "missing_saved_reference": 44,
        "saved_frame_reference": 16,
        "reused_committed_gap_R_q_y_reference": 4,
    }:
        raise _AuditAbort("mismatch status counts are not 44/16/4")
    endpoint_map = _mapping(mismatch.get("endpoints"), "mismatch endpoints")
    expected_labels = {entry[0] for entry in ENDPOINT_LAYOUT}
    if set(endpoint_map) != expected_labels:
        raise _AuditAbort("mismatch endpoint labels disagree")
    missing: list[dict[str, Any]] = []
    retained: list[dict[str, Any]] = []
    for label, _, _, _ in ENDPOINT_LAYOUT:
        endpoint = _mapping(endpoint_map.get(label), label + " mismatch endpoint")
        if endpoint.get("endpoint") != label:
            raise _AuditAbort(label + " mismatch endpoint label disagrees")
        fields = _mapping(endpoint.get("fields"), label + ".fields")
        archive_endpoint = _endpoint_from_c3(c3, label)
        for field in M15_ENDPOINT_FIELDS:
            selected = _mapping(fields.get(field), label + "." + field)
            originals = _list(
                selected.get("original_refs_Z0_Z1"),
                label + "." + field + ".refs",
            )
            if len(originals) != 2:
                raise _AuditAbort(label + "." + field + " does not have Z0/Z1")
            rows = _mapping(selected.get("rows"), label + "." + field + ".rows")
            archive_rows = _list(
                archive_endpoint.get(field),
                label + " C3 " + field,
            )
            if len(archive_rows) < 2:
                raise _AuditAbort(label + " C3 field lacks Z0/Z1")
            for order in ORDINARY_Z_ROWS:
                z_label = "Z" + str(order)
                original = _index(
                    originals[order],
                    label + "." + field + "." + z_label,
                    source_count,
                )
                row = _mapping(
                    rows.get(z_label),
                    label + "." + field + "." + z_label,
                )
                if _index(row.get("original_ref"), "mismatch original_ref", source_count) != original:
                    raise _AuditAbort(label + "." + field + "." + z_label + " ref disagrees")
                if _index(archive_rows[order], "C3 source row", source_count) != original:
                    raise _AuditAbort("C3 source row disagrees with mismatch row")
                item = {
                    "endpoint": label,
                    "field": field,
                    "z_order": order,
                    "z_label": z_label,
                    "original_ref": original,
                    "report_status": row.get("status"),
                }
                status = row.get("status")
                if status == "missing_saved_reference":
                    if str(original) in rewritten:
                        raise _AuditAbort("missing row is present in saved rewritten_nodes")
                    deferred = _mapping(
                        row.get("saved_inlet_substitution"),
                        "missing substitution",
                    )
                    if deferred.get("namespace") != DEFERRED_NAMESPACE:
                        raise _AuditAbort("missing substitution namespace disagrees")
                    if deferred.get("operation") != "simultaneous_saved_inlet_substitution":
                        raise _AuditAbort("missing substitution operation disagrees")
                    if deferred.get("mapping_pointer") != "$.Rp_native_C5_pulse_source_bindings.input_substitution":
                        raise _AuditAbort("missing substitution mapping pointer disagrees")
                    if deferred.get("mapping_sha256") != substitution_digest:
                        raise _AuditAbort("missing substitution digest disagrees")
                    if deferred.get("executed") is not False:
                        raise _AuditAbort("missing substitution is marked executed")
                    missing.append(item)
                elif status == "saved_frame_reference":
                    target = _index(row.get("node"), "saved report target", source_count)
                    if rewritten.get(str(original)) != target:
                        raise _AuditAbort("saved report target disagrees with rewritten_nodes")
                    retained.append(item | {"node": target})
                elif status == "reused_committed_gap_R_q_y_reference":
                    if label not in ("gap.lower", "gap.upper") or field != "R_q_y":
                        raise _AuditAbort("prior fresh status appears on unexpected field")
                    retained.append(item | {"binding": dict(row)})
                else:
                    raise _AuditAbort("unknown mismatch endpoint binding status")
    if len(missing) != EXPECTED_DEFERRED_OCCURRENCES:
        raise _AuditAbort("published mismatch scope is not exactly 44 missing rows")
    if len(retained) != EXPECTED_RETAINED_OCCURRENCES:
        raise _AuditAbort("published mismatch retained scope is not exactly 20 rows")
    return missing, retained


def _load_prior_index(
    root: Path,
    source_count: int,
    rewritten: Mapping[str, int],
) -> tuple[dict[int, list[str]], set[str], set[str], set[int], dict[str, Any]]:
    """Load only signed prior metadata and an original_ref index.

    No prior ``operator_row`` or prior operand is inspected here.  The signed
    audit result is treated as the prior DAG's structural attestation.
    """

    fresh, _ = _load_json_result(
        root,
        FRESH_REPORT_RELATIVE,
        FRESH_REPORT_SHA256,
        FRESH_RESULT_ARTIFACT,
        FRESH_RESULT_SHA256,
        "prior fresh substitution",
    )
    if fresh.get("schema") != "frame_gap_fresh_substitution_v1":
        raise _AuditAbort("prior fresh schema mismatch")
    if fresh.get("status") != "unevaluated_local_source_dag":
        raise _AuditAbort("prior fresh status mismatch")
    if fresh.get("fresh_node_count") != 1131 or fresh.get("missing_root_count") != 20:
        raise _AuditAbort("prior fresh count metadata mismatch")
    if fresh.get("frame_archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        raise _AuditAbort("prior fresh result is pinned to another archive")
    nodes = _list(fresh.get("nodes"), "prior fresh nodes")
    if len(nodes) != 1131:
        raise _AuditAbort("prior fresh node list count mismatch")
    prior_by_original: dict[int, list[str]] = {}
    prior_ids: set[str] = set()
    for raw_node in nodes:
        node = _mapping(raw_node, "prior fresh node")
        node_id = node.get("node_id")
        if not isinstance(node_id, str) or not node_id or node_id in prior_ids:
            raise _AuditAbort("prior fresh node id is not unique")
        original = _index(node.get("original_ref"), "prior fresh original_ref", source_count)
        if node.get("namespace") != PRIOR_FRESH_NAMESPACE:
            raise _AuditAbort("prior fresh node namespace mismatch")
        prior_ids.add(node_id)
        prior_by_original.setdefault(original, []).append(node_id)

    audit, _ = _load_json_result(
        root,
        AUDIT_REPORT_RELATIVE,
        AUDIT_REPORT_SHA256,
        AUDIT_RESULT_ARTIFACT,
        AUDIT_RESULT_SHA256,
        "prior fresh structural audit",
    )
    if audit.get("schema") != "frame_gap_fresh_substitution_audit_v1":
        raise _AuditAbort("prior structural audit schema mismatch")
    if audit.get("status") != "passed":
        raise _AuditAbort("prior structural audit is not passed")
    if audit.get("archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        raise _AuditAbort("prior audit archive hash mismatch")
    if audit.get("expected_archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        raise _AuditAbort("prior audit expected archive hash mismatch")
    if audit.get("expected_result_sha256") != FRESH_RESULT_SHA256:
        raise _AuditAbort("prior audit result hash mismatch")
    counts = _mapping(audit.get("counts"), "prior audit counts")
    if counts.get("fresh_node_count") != 1131 or counts.get("fresh_nodes_checked") != 1131:
        raise _AuditAbort("prior audit node count mismatch")
    if counts.get("expected_root_count") != 20 or counts.get("roots_checked") != 20:
        raise _AuditAbort("prior audit root count mismatch")
    if audit.get("failure_count") != 0 or audit.get("failure_details") != []:
        raise _AuditAbort("prior audit contains failures")
    semantic_flags = _mapping(audit.get("semantic_flags"), "prior audit semantic flags")
    required_prior = {"G_equals_F_proved", "P_equals_F_proved", "semantic_equivalence_proved", "physical_equivalence_proved", "source_values_computed", "numeric_values_computed", "integrals_evaluated", "admissible_for_numeric_use", "full_goal_complete"}
    if not required_prior.issubset(semantic_flags) or any(value is not False for value in semantic_flags.values()):
        raise _AuditAbort("prior audit semantic flags missing or not false")
    for key in (
        "source_values_computed",
        "integrals_evaluated",
        "numeric_values_computed",
        "old_builders_called",
        "old_constructors_called",
        "old_operators_replayed",
        "old_error_terms_removed",
        "extra_ledger_terms_added",
        "full_goal_complete",
    ):
        if audit.get(key) is not False:
            raise _AuditAbort("prior audit flag is not false: " + key)

    prior_root_ids: set[str] = set()
    prior_root_originals: set[int] = set()
    roots = _mapping(fresh.get("roots"), "prior fresh roots")

    def check_root_rows(value: Any, label: str) -> None:
        rows = _list(value, label)
        if len(rows) != 4:
            raise _AuditAbort(label + " must contain four rows")
        for position, raw_row in enumerate(rows):
            row = _mapping(raw_row, label + " row")
            original = _index(row.get("original_ref"), label + ".original_ref", source_count)
            if row.get("row_order") != position:
                raise _AuditAbort(label + " row order mismatch")
            node_ref = _mapping(row.get("node_ref"), label + ".node_ref")
            if node_ref.get("namespace") != PRIOR_FRESH_NAMESPACE:
                raise _AuditAbort(label + " root namespace mismatch")
            node_id = node_ref.get("node")
            if not isinstance(node_id, str) or node_id not in prior_ids:
                raise _AuditAbort(label + " root node is unknown")
            if node_id not in prior_by_original.get(original, []):
                raise _AuditAbort(label + " root original/node pair disagrees")
            prior_root_ids.add(node_id)
            prior_root_originals.add(original)

    check_root_rows(roots.get("gap_m15"), "prior gap_m15 roots")
    endpoint_roots = _mapping(roots.get("endpoints"), "prior endpoint roots")
    for side in ("lower", "upper"):
        side_roots = _mapping(endpoint_roots.get(side), "prior " + side + " roots")
        check_root_rows(side_roots.get("R_q_y"), "prior " + side + " R_q_y roots")
        check_root_rows(
            side_roots.get("q_plus_q_y"),
            "prior " + side + " q_plus_q_y roots",
        )
    radius = _mapping(roots.get("field_radius"), "prior radius root")
    radius_original = _index(radius.get("original_ref"), "prior radius original_ref", source_count)
    radius_ref = _mapping(radius.get("node_ref"), "prior radius node_ref")
    if radius_ref.get("namespace") != SAVED_NAMESPACE:
        raise _AuditAbort("prior radius boundary namespace mismatch")
    if rewritten.get(str(radius_original)) != radius_ref.get("node"):
        raise _AuditAbort("prior radius boundary is not the saved opaque node")
    if len(prior_root_ids) != 20:
        raise _AuditAbort("prior root node count is not 20")

    receipt = {
        "report_sha256": AUDIT_REPORT_SHA256,
        "result_sha256": AUDIT_RESULT_SHA256,
        "schema": audit.get("schema"),
        "status": audit.get("status"),
        "fresh_node_count": len(nodes),
        "fresh_root_count": len(prior_root_originals),
        "semantic_flags": dict(semantic_flags),
    }
    return prior_by_original, prior_ids, prior_root_ids, prior_root_originals, receipt


def _source_operand_values(value: Any) -> tuple[list[Any], bool]:
    if isinstance(value, list):
        return list(value), True
    return [value], False


def _check_boundary_ref(
    audit: _Audit,
    ref: Any,
    source_operand: Any,
    *,
    parent_label: str,
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
    new_ids: set[str],
    new_original_by_id: Mapping[str, int],
    new_context_by_id: Mapping[str, Mapping[str, Any]],
    parent_context: Mapping[str, Any],
    parent_position: int,
    position_by_id: Mapping[str, int],
) -> None:
    audit.counts["new_operand_refs_checked"] += 1
    if not _is_int(source_operand):
        audit.fail("source_operand_not_index", "source operand is not an integer", node=parent_label)
        return
    if not isinstance(ref, Mapping):
        audit.fail("operand_ref_shape", "operand replacement is not an object", node=parent_label)
        return
    namespace = ref.get("namespace")
    target = ref.get("node")
    if namespace == SAVED_NAMESPACE:
        audit.counts["saved_boundary_refs_checked"] += 1
        if not _is_int(target):
            audit.fail("saved_target_shape", "saved operand target is not an integer", node=parent_label)
            return
        if rewritten.get(str(source_operand)) != target:
            audit.fail(
                "saved_mapping_mismatch",
                "saved operand target differs from rewritten_nodes",
                node=parent_label,
                original_ref=source_operand,
                target=target,
            )
        return
    if namespace == PRIOR_FRESH_NAMESPACE:
        audit.counts["prior_boundary_refs_checked"] += 1
        if not isinstance(target, str) or target not in prior_ids:
            audit.fail("prior_target_missing", "prior operand target is unknown", node=parent_label)
            return
        if target not in prior_by_original.get(source_operand, []):
            audit.fail(
                "prior_original_mismatch",
                "prior target is not indexed under the source operand original_ref",
                node=parent_label,
                target=target,
                original_ref=source_operand,
            )
        return
    if namespace == FRESH_NAMESPACE:
        audit.counts["new_boundary_refs_checked"] += 1
        if not isinstance(target, str) or target not in new_ids:
            audit.fail("new_target_missing", "new operand target is unknown", node=parent_label)
            return
        if new_original_by_id.get(target) != source_operand:
            audit.fail(
                "new_original_mismatch",
                "new operand target original_ref differs from source operand",
                node=parent_label,
                target=target,
                expected=source_operand,
                actual=new_original_by_id.get(target),
            )
        if dict(new_context_by_id.get(target, {})) != dict(parent_context):
            audit.fail("new_root_context_mismatch", "new child crosses occurrence root_context", node=parent_label, target=target)
        child_position = position_by_id.get(target)
        audit.counts["postorder_checks"] += 1
        if child_position is None or child_position >= parent_position:
            audit.fail("new_postorder_violation", "new child is not before its parent", node=parent_label, target=target)
        return
    audit.fail("operand_namespace", "operand reference namespace is not saved, prior, or new", node=parent_label, namespace=namespace)


def _check_operator_row(
    audit: _Audit,
    node: Mapping[str, Any],
    position: int,
    source_nodes: list[Any],
    source_count: int,
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
    occurrence_node_ids: set[str],
    new_ids: set[str],
    new_original_by_id: Mapping[str, int],
    new_context_by_id: Mapping[str, Mapping[str, Any]],
    position_by_id: Mapping[str, int],
    expected_context: Mapping[str, Any],
) -> None:
    node_id = node.get("node_id")
    label = node_id if isinstance(node_id, str) else "<malformed>"
    audit.counts["new_nodes_checked"] += 1
    original = _index(node.get("original_ref"), label + ".original_ref", source_count)
    source_row = _mapping(source_nodes[original], label + " source row")
    audit.counts["source_rows_checked"] += 1
    operation = source_row.get("operation")
    if operation not in FROZEN_OPERAND_SCHEMA:
        audit.fail("operation_not_frozen", "source operation is not in frozen schema", node=label, operation=operation)
        return
    if node.get("namespace") != FRESH_NAMESPACE:
        audit.fail("new_node_namespace", "new node namespace mismatch", node=label)
    if node.get("root_context") != dict(expected_context):
        audit.fail("root_context_mismatch", "new node root_context differs from occurrence", node=label)
    operator_row = node.get("operator_row")
    if not isinstance(operator_row, Mapping):
        audit.fail("operator_row_shape", "new operator_row is not an object", node=label)
        return
    if operator_row.get("operation") != operation:
        audit.fail("operation_mismatch", "new operation differs from source", node=label)
    schema = FROZEN_OPERAND_SCHEMA[operation]
    if node.get("operand_schema") != list(schema):
        audit.fail("operand_schema_mismatch", "new operand_schema differs from frozen source schema", node=label)
    if set(operator_row) != set(source_row):
        audit.fail("operator_key_shape", "new operator row does not preserve source key shape", node=label)
    for key in set(source_row) - set(schema) - {"operation"}:
        if operator_row.get(key) != source_row.get(key):
            audit.fail("nonoperand_metadata_mismatch", "non-operand source attribute changed", node=label, field=key)

    parent_context = _mapping(node.get("root_context"), label + ".root_context")
    for field in schema:
        source_has = field in source_row
        new_has = field in operator_row
        if source_has != new_has:
            audit.fail("operand_field_presence", "operand field presence changed", node=label, field=field)
            continue
        if not source_has:
            continue
        source_value = source_row.get(field)
        new_value = operator_row.get(field)
        if source_value is None:
            if new_value is not None:
                audit.fail("operand_null_shape", "null operand changed", node=label, field=field)
            continue
        source_values, source_is_list = _source_operand_values(source_value)
        new_values, new_is_list = _source_operand_values(new_value)
        if source_is_list and field not in LIST_OPERAND_FIELDS:
            audit.fail("unexpected_list_operand", "non-list operand field is list-valued", node=label, field=field)
        if source_is_list != new_is_list:
            audit.fail("operand_list_shape", "scalar/list operand shape changed", node=label, field=field)
            continue
        if len(source_values) != len(new_values):
            audit.fail("operand_list_length", "operand list length changed", node=label, field=field)
            continue
        for source_operand, new_operand in zip(source_values, new_values):
            _check_boundary_ref(
                audit,
                new_operand,
                source_operand,
                parent_label=label,
                rewritten=rewritten,
                prior_by_original=prior_by_original,
                prior_ids=prior_ids,
                new_ids=occurrence_node_ids,
                new_original_by_id=new_original_by_id,
                new_context_by_id=new_context_by_id,
                parent_context=parent_context,
                parent_position=position,
                position_by_id=position_by_id,
            )


def _check_boundary_events(
    audit: _Audit,
    occurrence: Mapping[str, Any],
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
) -> None:
    events = occurrence.get("boundary_events")
    if not isinstance(events, list):
        audit.fail("boundary_events_shape", "occurrence boundary_events is not a list", occurrence=occurrence.get("occurrence_id"))
        return
    for event in events:
        if not isinstance(event, Mapping):
            audit.fail("boundary_event_shape", "boundary event is not an object")
            continue
        original = event.get("original_ref")
        if not _is_int(original):
            audit.fail("boundary_event_original", "boundary event original_ref is not an integer")
            continue
        kind = event.get("kind")
        namespace = event.get("namespace")
        node = event.get("node")
        if kind == "saved_rewritten_nodes":
            if namespace != SAVED_NAMESPACE or rewritten.get(str(original)) != node:
                audit.fail("saved_boundary_event", "saved boundary event disagrees with saved map", original_ref=original)
        elif kind == "prior_fresh_boundary":
            if namespace != PRIOR_FRESH_NAMESPACE or not isinstance(node, str) or node not in prior_ids:
                audit.fail("prior_boundary_event", "prior boundary event target is unknown", original_ref=original)
            elif node not in prior_by_original.get(original, []):
                audit.fail("prior_boundary_event_original", "prior boundary event original_ref disagrees", original_ref=original, node=node)
        else:
            audit.fail("boundary_event_kind", "unknown boundary event kind", kind=kind)


def _check_retained(
    audit: _Audit,
    retained_output: Any,
    expected_retained: list[dict[str, Any]],
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
) -> None:
    retained = _list(retained_output, "retained_occurrences")
    if len(retained) != EXPECTED_RETAINED_OCCURRENCES:
        audit.fail("retained_count", "retained occurrence count is not 20")
        return
    expected_roles = [
        (item["endpoint"], item["field"], item["z_order"], item["original_ref"], item["report_status"])
        for item in expected_retained
    ]
    actual_roles = []
    for position, raw in enumerate(retained):
        audit.counts["retained_occurrences_checked"] += 1
        item = _mapping(raw, "retained occurrence")
        role = (item.get("endpoint"), item.get("field"), item.get("z_order"), item.get("original_ref"), item.get("status"))
        actual_roles.append(role)
        if role != expected_roles[position]:
            audit.fail("retained_order", "retained occurrence role/order differs from mismatch receipt", position=position)
        original = item.get("original_ref")
        if not _is_int(original):
            continue
        source_ref = _mapping(item.get("source_ref"), "retained source_ref")
        status = item.get("status")
        if status == "saved_frame_reference":
            if source_ref.get("namespace") != SAVED_NAMESPACE or rewritten.get(str(original)) != source_ref.get("node"):
                audit.fail("retained_saved_ref", "retained saved ref disagrees with saved map", position=position)
        elif status == "reused_committed_gap_R_q_y_reference":
            endpoint = item.get("endpoint")
            side = endpoint.split(".", 1)[1] if isinstance(endpoint, str) and "." in endpoint else None
            expected = EXPECTED_GAP_R_Q_Y.get(side, ())
            if endpoint not in ("gap.lower", "gap.upper") or item.get("field") != "R_q_y" or item.get("z_order") not in ORDINARY_Z_ROWS or original != expected[item.get("z_order")]:
                audit.fail("retained_prior_scope", "retained prior fresh ref is outside the four gap R_q_y rows", position=position)
            if source_ref.get("namespace") != PRIOR_FRESH_NAMESPACE or not isinstance(source_ref.get("node"), str) or source_ref.get("node") not in prior_by_original.get(original, []) or source_ref.get("node") not in prior_ids:
                audit.fail("retained_prior_ref", "retained prior ref disagrees with prior index", position=position)
        else:
            audit.fail("retained_status", "unknown retained status", position=position)
    if actual_roles != expected_roles:
        audit.fail("retained_scope", "retained occurrence scope/order differs from published mismatch")


def _check_occurrences(
    audit: _Audit,
    result: Mapping[str, Any],
    expected_missing: list[dict[str, Any]],
    source_nodes: list[Any],
    source_count: int,
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
    prior_root_ids: set[str],
    prior_root_originals: set[int],
) -> None:
    occurrences = _list(result.get("occurrences"), "deferred occurrences")
    if len(occurrences) != EXPECTED_DEFERRED_OCCURRENCES:
        audit.fail("occurrence_count", "deferred occurrence count is not 44")
        return
    actual_new_roots = 0
    actual_prior_roots = 0
    seen_root_nodes: set[str] = set()
    retained_prior_ids: set[str] = set()
    for raw_retained in _list(result.get("retained_occurrences"), "retained occurrences"):
        retained_item = _mapping(raw_retained, "retained occurrence")
        if retained_item.get("status") != "reused_committed_gap_R_q_y_reference":
            continue
        retained_ref = retained_item.get("source_ref")
        if (
            isinstance(retained_ref, Mapping)
            and retained_ref.get("namespace") == PRIOR_FRESH_NAMESPACE
            and isinstance(retained_ref.get("node"), str)
        ):
            retained_prior_ids.add(retained_ref["node"])
    for position, raw in enumerate(occurrences):
        audit.counts["deferred_occurrences_checked"] += 1
        occurrence = _mapping(raw, "deferred occurrence")
        expected = expected_missing[position]
        expected_role = (expected["endpoint"], expected["field"], expected["z_order"], expected["z_label"], expected["original_ref"], expected["report_status"])
        actual_role = (occurrence.get("endpoint"), occurrence.get("field"), occurrence.get("z_order"), occurrence.get("z_label"), occurrence.get("original_ref"), occurrence.get("report_status"))
        if actual_role != expected_role:
            audit.fail("occurrence_role_order", "occurrence role/order differs from mismatch receipt", position=position)
        occurrence_id = "occ%04d" % position
        if occurrence.get("occurrence_id") != occurrence_id:
            audit.fail("occurrence_id", "occurrence id is not the frozen sequence id", position=position)
        expected_context = {
            "occurrence_id": occurrence_id,
            "endpoint": expected["endpoint"],
            "field": expected["field"],
            "z_order": expected["z_order"],
            "z_label": expected["z_label"],
        }
        if occurrence.get("role") != expected["endpoint"] + "." + expected["field"] + "." + expected["z_label"]:
            audit.fail("occurrence_role", "occurrence role string mismatch", position=position)
        if occurrence.get("dag_namespace") != FRESH_NAMESPACE:
            audit.fail("occurrence_namespace", "occurrence dag namespace mismatch", position=position)
        if occurrence.get("original_ref") != expected["original_ref"]:
            audit.fail("occurrence_original_ref", "occurrence original_ref mismatch", position=position)
        root_ref = _mapping(occurrence.get("root_ref"), "occurrence root_ref")
        source_ref = _mapping(occurrence.get("source_ref"), "occurrence source_ref")
        if dict(root_ref) != dict(source_ref):
            audit.fail("root_source_ref", "root_ref and source_ref disagree", position=position)
        root_status = occurrence.get("root_status")
        root_node = root_ref.get("node")
        node_refs = _list(occurrence.get("dag_node_refs"), "occurrence dag_node_refs")
        if occurrence.get("dag_node_count") != len(node_refs):
            audit.fail("occurrence_node_count", "dag_node_count differs from dag_node_refs", position=position)
        if root_status == "new_fresh_dag":
            actual_new_roots += 1
            if root_ref.get("namespace") != FRESH_NAMESPACE or not isinstance(root_node, str):
                audit.fail("new_root_ref", "new root_ref is not a fresh node", position=position)
            last_ref = node_refs[-1] if node_refs else None
            if (
                not isinstance(last_ref, Mapping)
                or not node_refs
                or not isinstance(root_node, str)
                or last_ref.get("node") != root_node
            ):
                audit.fail("new_root_postorder", "new root is not the final occurrence node", position=position)
        elif root_status == "prior_fresh_boundary":
            actual_prior_roots += 1
            if root_ref.get("namespace") != PRIOR_FRESH_NAMESPACE or not isinstance(root_node, str):
                audit.fail("prior_root_ref", "prior root_ref is not a prior fresh node", position=position)
            # A previously audited internal node can be a newly requested
            # root; require the exact original-ref pairing and prior membership.
            if not isinstance(root_node, str) or root_node not in prior_ids:
                audit.fail("prior_root_scope", "deferred root is not an audited prior node", position=position)
            elif root_node not in prior_by_original.get(expected["original_ref"], []):
                audit.fail("prior_root_original_pair", "deferred prior node does not match its original_ref", position=position)
            if node_refs:
                audit.fail("prior_root_nodes", "prior boundary root unexpectedly owns new nodes", position=position)
        else:
            audit.fail("root_status", "unknown deferred root status", position=position)
        if isinstance(root_node, str):
            if root_node in seen_root_nodes:
                audit.fail("duplicate_root_node", "same root node is used by multiple occurrences", position=position, node=root_node)
            seen_root_nodes.add(root_node)
        _check_boundary_events(audit, occurrence, rewritten, prior_by_original, prior_ids)
    if actual_new_roots != EXPECTED_NEW_ROOT_COUNT:
        audit.fail("new_root_count", "deferred new root count is not 32", actual=actual_new_roots)
    if actual_prior_roots != EXPECTED_PRIOR_ROOT_COUNT:
        audit.fail("prior_root_count", "deferred prior root count is not 12", actual=actual_prior_roots)
    if result.get("prior_fresh_root_reuses") != EXPECTED_PRIOR_ROOT_COUNT:
        audit.fail("prior_root_reuse_count", "declared prior_fresh_root_reuses is not 12")


def _check_new_nodes(
    audit: _Audit,
    result: Mapping[str, Any],
    expected_missing: list[dict[str, Any]],
    source_nodes: list[Any],
    source_count: int,
    rewritten: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    prior_ids: set[str],
) -> None:
    nodes = _list(result.get("nodes"), "new deferred nodes")
    if len(nodes) != EXPECTED_NEW_NODE_COUNT or result.get("fresh_node_count") != EXPECTED_NEW_NODE_COUNT:
        audit.fail("new_node_count", "new deferred node count is not 168", actual=len(nodes))
    nodes_by_id: dict[str, Mapping[str, Any]] = {}
    position_by_id: dict[str, int] = {}
    for position, raw in enumerate(nodes):
        node = _mapping(raw, "new deferred node")
        node_id = node.get("node_id")
        if not isinstance(node_id, str) or not node_id:
            audit.fail("new_node_id", "new node id is not a nonempty string", position=position)
            continue
        if node_id in nodes_by_id:
            audit.fail("new_node_id_duplicate", "new node id is duplicated", node=node_id)
            continue
        nodes_by_id[node_id] = node
        position_by_id[node_id] = position
    occurrence_ids = {"occ%04d" % index for index in range(EXPECTED_DEFERRED_OCCURRENCES)}
    node_ids_by_occurrence: dict[str, set[str]] = {key: set() for key in occurrence_ids}
    node_context_by_id: dict[str, Mapping[str, Any]] = {}
    node_original_by_id: dict[str, int] = {}
    for node_id, node in nodes_by_id.items():
        context = node.get("root_context")
        if isinstance(context, Mapping) and isinstance(context.get("occurrence_id"), str):
            node_ids_by_occurrence.setdefault(context["occurrence_id"], set()).add(node_id)
            node_context_by_id[node_id] = context
        original = node.get("original_ref")
        if _is_int(original):
            node_original_by_id[node_id] = original
    for occurrence_index, expected in enumerate(expected_missing):
        occurrence_id = "occ%04d" % occurrence_index
        occurrence_context = {
            "occurrence_id": occurrence_id,
            "endpoint": expected["endpoint"],
            "field": expected["field"],
            "z_order": expected["z_order"],
            "z_label": expected["z_label"],
        }
        output_occurrences = _list(result.get("occurrences"), "deferred occurrences")
        occurrence = _mapping(output_occurrences[occurrence_index], "deferred occurrence")
        refs = _list(occurrence.get("dag_node_refs"), "occurrence dag_node_refs")
        occurrence_node_ids: set[str] = set()
        for raw_ref in refs:
            ref = _mapping(raw_ref, "occurrence dag node ref")
            if ref.get("namespace") != FRESH_NAMESPACE or not isinstance(ref.get("node"), str):
                audit.fail("occurrence_node_ref", "occurrence node ref namespace/shape mismatch", occurrence=occurrence_id)
                continue
            node_id = ref["node"]
            occurrence_node_ids.add(node_id)
            if node_id not in nodes_by_id:
                audit.fail("occurrence_node_missing", "occurrence references unknown new node", occurrence=occurrence_id, node=node_id)
        if occurrence_node_ids != node_ids_by_occurrence.get(occurrence_id, set()):
            audit.fail("occurrence_node_membership", "occurrence dag refs do not match node root_context membership", occurrence=occurrence_id)
        ordered_ids: list[str] = []
        for raw_ref in refs:
            if not isinstance(raw_ref, Mapping):
                continue
            if raw_ref.get("namespace") == FRESH_NAMESPACE and isinstance(raw_ref.get("node"), str):
                ordered_ids.append(raw_ref["node"])
        ordered_positions = [
            position_by_id[node_id]
            for node_id in ordered_ids
            if node_id in position_by_id
        ]
        if len(ordered_ids) != len(set(ordered_ids)):
            audit.fail("occurrence_duplicate_node_ref", "occurrence dag refs contain a duplicate", occurrence=occurrence_id)
        if ordered_positions and ordered_positions != sorted(ordered_positions):
            audit.fail("occurrence_postorder_sequence", "occurrence dag refs are not in saved node order", occurrence=occurrence_id)
        if ordered_ids and ordered_ids[-1] in nodes_by_id:
            root_ref = _mapping(occurrence.get("root_ref"), "occurrence root_ref")
            if root_ref.get("namespace") == FRESH_NAMESPACE and root_ref.get("node") != ordered_ids[-1]:
                audit.fail("root_ref_postorder", "new root_ref is not the final occurrence node", occurrence=occurrence_id)
        root_ref = _mapping(occurrence.get("root_ref"), "occurrence root_ref")
        if root_ref.get("namespace") == FRESH_NAMESPACE and isinstance(root_ref.get("node"), str):
            if node_original_by_id.get(root_ref["node"]) != expected["original_ref"]:
                audit.fail("new_root_original_pair", "new root node original_ref differs from occurrence", occurrence=occurrence_id)
        for node_id in occurrence_node_ids:
            node = nodes_by_id.get(node_id)
            if node is None:
                continue
            position = position_by_id[node_id]
            _check_operator_row(
                audit,
                node,
                position,
                source_nodes,
                source_count,
                rewritten,
                prior_by_original,
                prior_ids,
                occurrence_node_ids,
                set(nodes_by_id),
                node_original_by_id,
                node_context_by_id,
                position_by_id,
                occurrence_context,
            )
    referenced = {
        ref["node"]
        for raw in _list(result.get("occurrences"), "deferred occurrences")
        for ref in _list(_mapping(raw, "deferred occurrence").get("dag_node_refs"), "dag refs")
        if isinstance(ref, Mapping) and ref.get("namespace") == FRESH_NAMESPACE and isinstance(ref.get("node"), str)
    }
    if referenced != set(nodes_by_id):
        audit.fail("unowned_new_nodes", "new node list and occurrence references are not identical")


def _check_result_contract(
    audit: _Audit,
    result: Mapping[str, Any],
    archive_digest: str,
    substitution_digest: str,
    prior_receipt: Mapping[str, Any],
) -> None:
    if result.get("schema") != "frame_m15_deferred_endpoint_substitution_v1":
        audit.fail("result_schema", "deferred result schema mismatch")
    if result.get("status") != "unevaluated_deferred_source_dag":
        audit.fail("result_status", "deferred result status mismatch")
    if result.get("namespace") != FRESH_NAMESPACE:
        audit.fail("result_namespace", "deferred result namespace mismatch")
    if result.get("endpoint_fields") != list(M15_ENDPOINT_FIELDS):
        audit.fail("endpoint_fields", "deferred endpoint fields mismatch")
    if result.get("ordinary_Z_rows") != list(ORDINARY_Z_ROWS):
        audit.fail("ordinary_Z_rows", "deferred ordinary Z rows mismatch")
    if result.get("deferred_occurrence_count") != EXPECTED_DEFERRED_OCCURRENCES:
        audit.fail("deferred_count", "deferred occurrence count mismatch")
    if result.get("retained_occurrence_counts") != {
        "missing_saved_reference_deferred": EXPECTED_DEFERRED_OCCURRENCES,
        "saved_frame_reference_retained": EXPECTED_SAVED_RETAINED_COUNT,
        "prior_fresh_reference_retained": EXPECTED_PRIOR_RETAINED_COUNT,
    }:
        audit.fail("retained_counts", "retained occurrence counts mismatch")
    provenance = _mapping(result.get("provenance"), "deferred provenance")
    for key, expected in {
        "archive_raw_sha256": archive_digest,
        "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
        "C3_archive_name": C3_ARCHIVE_NAME,
        "C3_archive_sha256": C3_ARCHIVE_SHA256,
        "saved_input_substitution_sha256": substitution_digest,
        "fresh_structural_audit_report_sha256": AUDIT_REPORT_SHA256,
        "fresh_structural_audit_result_sha256": AUDIT_RESULT_SHA256,
    }.items():
        if provenance.get(key) != expected:
            audit.fail("provenance", "deferred provenance mismatch", field=key)
    if provenance.get("prior_audit_receipt") != dict(prior_receipt):
        audit.fail("prior_receipt", "deferred prior audit receipt differs")
    input_selection = _mapping(result.get("input_selection"), "input_selection")
    if input_selection.get("source_frame_inputs_checked") is not True or input_selection.get("source_frame_endpoint_rows_checked") != 64 or input_selection.get("missing_status_selected_from_published_result") is not True:
        audit.fail("input_selection", "deferred input selection flags mismatch")
    flags = _mapping(result.get("semantic_flags"), "deferred semantic flags")
    required_flags = {"semantic_identity_proved", "ordinary_Z1_commutation_proved", "Z1_commutation_proved", "seam_identity_proved", "endpoint_equalities_proved", "P_equals_F_proved", "G_equals_F_proved", "source_values_computed", "numeric_values_computed", "integrals_evaluated", "physical_equivalence_proved", "admissible_for_numeric_use", "old_builders_called", "old_constructors_called", "old_operators_replayed", "old_error_terms_removed", "extra_ledger_terms_added", "full_goal_complete"}
    if not required_flags.issubset(flags) or any(value is not False for value in flags.values()):
        audit.fail("semantic_flags", "deferred semantic flags missing or not false")
    for key in (
        "old_twenty_roots_replayed",
        "old_builders_called",
        "old_operators_replayed",
        "old_error_terms_removed",
        "extra_ledger_terms_added",
        "source_values_computed",
        "integrals_evaluated",
        "numeric_values_computed",
        "admissible_for_numeric_use",
        "full_goal_complete",
    ):
        if result.get(key) is not False:
            audit.fail("safety_flag", "deferred safety flag is not false", flag=key)


def recover(root: Path) -> dict[str, Any]:
    """Return a failure-collecting structural receipt; this is opt-in."""

    audit = _Audit()
    result_digest: str | None = None
    archive_digest: str | None = None
    try:
        result, result_digest = _load_json_result(
            root,
            RESULT_REPORT_RELATIVE,
            RESULT_REPORT_SHA256,
            RESULT_ARTIFACT,
            RESULT_SHA256,
            "deferred m15 substitution",
        )
        archive, archive_digest = _load_archive(root)
        source_nodes, rewritten, _substitution, substitution_digest, c3 = _load_source_context(archive, audit)
        source_count = len(source_nodes)
        mismatch, _ = _load_json_result(
            root,
            MISMATCH_REPORT_RELATIVE,
            MISMATCH_REPORT_SHA256,
            MISMATCH_RESULT_ARTIFACT,
            MISMATCH_RESULT_SHA256,
            "frame m15 mismatch",
        )
        expected_missing, expected_retained = _expected_scope(
            mismatch,
            c3,
            source_count,
            rewritten,
            substitution_digest,
        )
        prior_by_original, prior_ids, prior_root_ids, prior_root_originals, prior_receipt = _load_prior_index(root, source_count, rewritten)
        _check_result_contract(audit, result, archive_digest, substitution_digest, prior_receipt)
        _check_retained(
            audit,
            result.get("retained_occurrences"),
            expected_retained,
            rewritten,
            prior_by_original,
            prior_ids,
        )
        _check_occurrences(
            audit,
            result,
            expected_missing,
            source_nodes,
            source_count,
            rewritten,
            prior_by_original,
            prior_ids,
            prior_root_ids,
            prior_root_originals,
        )
        _check_new_nodes(
            audit,
            result,
            expected_missing,
            source_nodes,
            source_count,
            rewritten,
            prior_by_original,
            prior_ids,
        )
    except _AuditAbort as exc:
        if not audit.failures or audit.failures[-1].get("detail") != str(exc):
            audit.fail("audit_abort", str(exc))
    except Exception as exc:  # pragma: no cover - defensive failure receipt
        audit.fail("audit_exception", "unexpected audit exception: " + str(exc))
    return audit.report(result_digest=result_digest, archive_digest=archive_digest)

