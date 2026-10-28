"""Independent, unentered source-only bridge for the physical ``log Rp`` port.

This draft is deliberately narrower than the earlier geometry bridge.  It
consumes only frozen metadata and four literal graph rows.  It does not import
the prior bridge, traverse a graph, compare chart roots, inspect measures or
endpoints, evaluate a source value, or consume a Gamma expectation.

The only conditional statement admitted here is:

``log(current Rp) = old log Rp``

under the saved ``absolute-radius log recipe only`` source identity, the
literal old row ``old Rp = exp(old log Rp)``, and the recorded real-positivity
premise.  The real rule ``log(exp(x)) = x`` is recorded as a rule premise; no
numeric evaluation is performed.  The Gamma report contributes only its
authenticated ``imported_logRp`` / ``schedule.logRp`` source-port witness.

Raw gzip bytes are authenticated before decompression.  Every consumed value
is read through an exact JSON pointer and missing or malformed pointers block
the result.  Hashes pin provenance; they never establish function identity.
There is intentionally no command entry point, test call, or import of an
ancestor/source consumer in this file.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from typing import Any, Mapping


SCHEMA = "original_native_logRp_source_bridge_v1_unentered"

BASE_GRAPH = "5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e"

GEOMETRY_RESULT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_pulse_geometry_literal_bridge.json"
)
SCALE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json"
)
TERMINAL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_terminal_log_repair.json"
)
MASTER_RAW = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
PHYSICAL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json"
)
GAMMA_RAW = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.json.gz"
)
GAMMA_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions_check.json"
)

# These hashes are source pins only.  The two gzip records are decompressed
# only after their raw-byte hash passes.
EXPECTED_RAW_HASHES = {
    GEOMETRY_RESULT: "0bd4c50246150cf8595ad2c3bbf4e7a1866e607c4e64de6624f72a940ba939f9",
    SCALE_REPORT: "b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760",
    TERMINAL_REPORT: "9566d63baaf92856b8210e1ce07824a5b12c97b6af7a641bfa4c5cab36f5f040",
    MASTER_RAW: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    PHYSICAL_REPORT: "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    GAMMA_RAW: "dbc3dd702930b85599e7967b5e14d934f3c326a731cd1b9528990158ff74c612",
    GAMMA_CHECK: "ee6fc38317768a68b5775f6224eeb250f6383ee3f732d83c68b30316c720f927",
}

EXPECTED_DECODED_HASHES = {
    MASTER_RAW: "56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1",
    GAMMA_RAW: "957473c9479a4e0945af4966a99d5c2b95e7a186f63c521304d8c141e15aa07e",
}

OLD_LOG_RP = 31161
OLD_RP = 31162
CURRENT_RP = 263349
CURRENT_LOG_RP = 280182
GAMMA_IMPORTED_LOG_RP = 31161

EXPECTED_RP_ALIAS = {
    "source_graph_digest": BASE_GRAPH,
    "native_node": CURRENT_RP,
    "old_node": OLD_RP,
    "identity_scope": "absolute-radius log recipe only",
    "adapter_alias": False,
}

EXPECTED_OLD_LOG_RP_ROW = {
    "operation": "sum",
    "arguments": [31157, 31160],
}
EXPECTED_OLD_RP_ROW = {
    "operation": "analytic_unary",
    "name": "exp",
    "argument": OLD_LOG_RP,
}
EXPECTED_CURRENT_RP_ROW = {
    "operation": "function_substitution",
    "expression": 263348,
    "variable": 263344,
    "value": 263346,
    "Z_independent_substitution": True,
}
EXPECTED_CURRENT_LOG_RP_ROW = {
    "operation": "analytic_unary",
    "name": "log",
    "argument": CURRENT_RP,
    "source_positive_certificate": "actual native Rp>0",
    "ordinary_Z1_exact_zero": True,
}

LOG_EXP_RULE = {
    "name": "real_principal_log_exp",
    "statement": "log(exp(x)) = x",
    "domain": "x is real and exp(x) is in the positive real domain of principal log",
    "numeric_evaluation": False,
}


class SourceAuthenticationError(ValueError):
    """Raised internally when a frozen source record cannot be authenticated."""


_MISSING = object()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _workspace_path(workspace_root: os.PathLike[str] | str, relative_path: str) -> str:
    root = os.fspath(workspace_root)
    if not isinstance(root, str) or not isinstance(relative_path, str):
        raise TypeError("workspace_root and relative_path must be strings")
    return os.path.join(root, *relative_path.split("/"))


def _authenticated_document(
    workspace_root: os.PathLike[str] | str,
    relative_path: str,
    *,
    gzip_input: bool,
) -> dict[str, Any]:
    """Read one frozen JSON record, hashing raw bytes before any decode."""

    path = _workspace_path(workspace_root, relative_path)
    with open(path, "rb") as stream:
        raw = stream.read()
    raw_sha256 = _sha256(raw)
    expected_raw = EXPECTED_RAW_HASHES.get(relative_path)
    if expected_raw is None or raw_sha256 != expected_raw:
        raise SourceAuthenticationError("raw SHA-256 mismatch: " + relative_path)

    decoded = gzip.decompress(raw) if gzip_input else raw
    decoded_sha256 = _sha256(decoded)
    expected_decoded = EXPECTED_DECODED_HASHES.get(relative_path)
    if expected_decoded is not None and decoded_sha256 != expected_decoded:
        raise SourceAuthenticationError("decoded SHA-256 mismatch: " + relative_path)

    try:
        document = json.loads(decoded.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SourceAuthenticationError("JSON decode failed: " + relative_path) from error
    if not isinstance(document, dict):
        raise SourceAuthenticationError("source document is not an object: " + relative_path)
    return {
        "path": relative_path,
        "raw_sha256": raw_sha256,
        "decoded_sha256": decoded_sha256,
        "raw_authenticated_before_decode": True,
        "gzip": gzip_input,
        "document": document,
    }


def _pointer_value(document: Any, pointer: str) -> Any:
    """Read one exact RFC 6901 pointer; never fall back to a nearby key."""

    if not isinstance(pointer, str) or (pointer and not pointer.startswith("/")):
        raise KeyError("invalid JSON pointer")
    value = document
    if not pointer:
        return value
    for escaped in pointer[1:].split("/"):
        segment = escaped.replace("~1", "/").replace("~0", "~")
        if isinstance(value, dict):
            if segment not in value:
                raise KeyError(pointer)
            value = value[segment]
        elif isinstance(value, list):
            if segment == "-" or not segment.isdigit():
                raise KeyError(pointer)
            index = int(segment)
            if index >= len(value):
                raise KeyError(pointer)
            value = value[index]
        else:
            raise KeyError(pointer)
    return value


def _record_ok(record: Any, path: str, *, gzip_input: bool) -> bool:
    if not isinstance(record, Mapping):
        return False
    if record.get("path") != path:
        return False
    if record.get("raw_sha256") != EXPECTED_RAW_HASHES.get(path):
        return False
    if record.get("raw_authenticated_before_decode") is not True:
        return False
    if record.get("gzip") is not gzip_input:
        return False
    expected_decoded = EXPECTED_DECODED_HASHES.get(path)
    if expected_decoded is not None and record.get("decoded_sha256") != expected_decoded:
        return False
    return isinstance(record.get("document"), dict)


def _read(record: Any, pointer: str) -> Any:
    if not isinstance(record, Mapping) or not isinstance(record.get("document"), dict):
        return _MISSING
    try:
        return _pointer_value(record["document"], pointer)
    except (KeyError, TypeError, ValueError):
        return _MISSING


def _optional_read(record: Any, pointer: str) -> Any:
    value = _read(record, pointer)
    return None if value is _MISSING else value


def _require(
    checks: dict[str, dict[str, Any]],
    failures: list[str],
    name: str,
    passed: bool,
    detail: str,
) -> bool:
    ok = passed is True
    checks[name] = {"passed": ok, "detail": detail}
    if not ok:
        failures.append(name)
    return ok


def _preserved_scope(geometry: Any) -> dict[str, Any]:
    """Keep prior limits explicit; this bridge cannot remove any old error."""

    review_limitations = _optional_read(geometry, "/review_limitations")
    return {
        "old_errors_removed": False,
        "chart_measures_proved": False,
        "chart_endpoints_proved": False,
        "gamma_expectation_equality_proved": False,
        "m15_identity_proved": False,
        "physical_pulse_correction_zero_proved": False,
        "finite_N": {
            "materialized": False,
            "numeric_value_used": False,
            "node_handle_semantics": "preserved; N is not consumed by this bridge",
        },
        "prior_geometry_result": {
            "all_old_error_data_retained": _optional_read(
                geometry, "/all_old_error_data_retained"
            ),
            "reserved_contribution_registered_once": _optional_read(
                geometry, "/reserved_contribution_registered_once"
            ),
            "review_limitations": review_limitations,
        },
    }


def build_logrp_bridge(source_records: Mapping[str, Any]) -> dict[str, Any]:
    """Build a fail-closed source witness from already authenticated records.

    This function does not load files, import source modules, evaluate rows, or
    traverse any graph.  ``load_logrp_bridge_once`` is the only function that
    authenticates and parses the six frozen records.
    """

    checks: dict[str, dict[str, Any]] = {}
    failures: list[str] = []
    records = source_records if isinstance(source_records, Mapping) else {}
    geometry = records.get("geometry_result")

    record_specs = {
        "geometry_result": (GEOMETRY_RESULT, False),
        "scale": (SCALE_REPORT, False),
        "terminal": (TERMINAL_REPORT, False),
        "master": (MASTER_RAW, True),
        "physical": (PHYSICAL_REPORT, False),
        "gamma": (GAMMA_RAW, True),
        "gamma_check": (GAMMA_CHECK, False),
    }
    for role, (path, gzip_input) in record_specs.items():
        _require(
            checks,
            failures,
            "authenticated_record_" + role,
            _record_ok(records.get(role), path, gzip_input=gzip_input),
            "raw SHA and source document attestation are required; gzip is decoded only after raw SHA",
        )

    geometry_status = _read(geometry, "/result_status")
    _require(
        checks,
        failures,
        "accepted_geometry_only_result",
        geometry_status == "geometry-only"
        and _read(geometry, "/checks_passed") == 51
        and _read(geometry, "/source_only") is True
        and _read(geometry, "/prior_Rp_identity_is_a_scoped_premise") is True
        and _read(geometry, "/chart_measures_proved") is False
        and _read(geometry, "/chart_endpoints_proved") is False
        and _read(geometry, "/full_pulse_pressure_identity_proved") is False
        and _read(geometry, "/physical_pulse_correction_zero_proved") is False,
        "consume the frozen geometry-only result as a scoped premise only",
    )

    alias = _read(records.get("scale"), "/component_result/source_interface/source_aliases/Rp")
    _require(
        checks,
        failures,
        "scoped_Rp_alias_source_premise",
        alias == EXPECTED_RP_ALIAS
        and _read(records.get("scale"), "/component_result/source_only") is True
        and _read(records.get("scale"), "/component_result/native_original_absolute_Rc_Rp_identity_proved") is True,
        "only the saved absolute-radius log recipe alias is admitted; no global node alias",
    )

    radius_audit = _read(
        records.get("terminal"), "/component_result/absolute_radius_source_audit"
    )
    _require(
        checks,
        failures,
        "terminal_radius_scope_and_positive_laws",
        isinstance(radius_audit, dict)
        and radius_audit.get("source_graph_digest") == BASE_GRAPH
        and radius_audit.get("native_Rp") == CURRENT_RP
        and radius_audit.get("old_Rp") == OLD_RP
        and radius_audit.get("exact_positive_radius_log_laws_only") is True
        and radius_audit.get("source_values_or_bounds_evaluated") is False
        and radius_audit.get("old_math_consumers_called") == 0
        and radius_audit.get("control_alias_proved") is False,
        "terminal audit supplies positivity and scope; its unproved control alias is not inferred",
    )

    old_log_row = _read(records.get("master"), "/exact_graph_nodes/31161")
    old_rp_row = _read(records.get("master"), "/exact_graph_nodes/31162")
    current_rp_row = _read(records.get("master"), "/exact_graph_nodes/263349")
    current_log_handle = _read(records.get("physical"), "/component_result/physical_pulse_bindings/current_logRp")
    current_log_row = _read(
        records.get("physical"), "/component_result/definition_graph/appended_nodes/0"
    )
    _require(
        checks,
        failures,
        "literal_logRp_and_Rp_rows",
        old_log_row == EXPECTED_OLD_LOG_RP_ROW
        and old_rp_row == EXPECTED_OLD_RP_ROW
        and current_rp_row == EXPECTED_CURRENT_RP_ROW
        and current_log_handle == CURRENT_LOG_RP
        and current_log_row == EXPECTED_CURRENT_LOG_RP_ROW,
        "31161, 31162, 263349 and physical current_logRp 280182 match exact literal rows",
    )

    _require(
        checks,
        failures,
        "real_log_exp_rule_domain",
        isinstance(current_log_row, dict)
        and current_log_row.get("source_positive_certificate") == "actual native Rp>0"
        and radius_audit is not _MISSING
        and isinstance(radius_audit, dict)
        and radius_audit.get("exact_positive_radius_log_laws_only") is True,
        "principal real log is used only on the recorded positive-radius source domain",
    )

    gamma_imported = _read(
        records.get("gamma"), "/function_bindings/native_binding/imported_logRp"
    )
    gamma_schedule = _read(records.get("gamma"), "/function_bindings/schedule/logRp")
    _require(
        checks,
        failures,
        "Gamma_imported_logRp_source_port",
        _read(records.get("gamma_check"), "/all_passed") is True
        and _read(records.get("gamma_check"), "/current_original_n1_exact_Gamma_three_future_C1_functions_installed") is True
        and _read(records.get("gamma_check"), "/input_hashes/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.json.gz") == EXPECTED_RAW_HASHES[GAMMA_RAW]
        and gamma_imported == GAMMA_IMPORTED_LOG_RP
        and gamma_schedule == GAMMA_IMPORTED_LOG_RP,
        "Gamma report contributes imported/schedule logRp only; downstream expectation is not consumed",
    )

    passed = not failures
    preserved = _preserved_scope(geometry)
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "passed" if passed else "blocked",
        "passed": passed,
        "source_only": True,
        "checks": checks,
        "failures": failures,
        "rule": dict(LOG_EXP_RULE),
        "source_identity_scope": "absolute-radius log recipe only",
        "source_identity_is_global_alias": False,
        "proof": {
            "statement": "Under the saved scoped Rp source identity, log(263349)=31161",
            "conditional": True,
            "numeric_evaluation": False,
            "steps": [
                "current Rp source function = old Rp source function (saved scoped alias)",
                "old Rp source row = exp(old log Rp source row)",
                "current log Rp source row = log(current Rp source function)",
                "real principal log(exp(x)) = x on the recorded positive domain",
            ],
            "current_logRp_node": CURRENT_LOG_RP,
            "old_logRp_node": OLD_LOG_RP,
            "current_Rp_node": CURRENT_RP,
            "old_Rp_node": OLD_RP,
        },
        "Gamma_source_witness": {
            "imported_logRp": gamma_imported,
            "schedule_logRp": gamma_schedule,
            "expectation_equality_proved": False,
            "m15_or_correction_claimed": False,
        },
        "source_pointers": {
            "geometry_result": {"path": GEOMETRY_RESULT, "pointer": "/result_status"},
            "Rp_alias": {
                "path": SCALE_REPORT,
                "pointer": "/component_result/source_interface/source_aliases/Rp",
            },
            "terminal_radius_audit": {
                "path": TERMINAL_REPORT,
                "pointer": "/component_result/absolute_radius_source_audit",
            },
            "old_logRp": {"path": MASTER_RAW, "pointer": "/exact_graph_nodes/31161"},
            "old_Rp": {"path": MASTER_RAW, "pointer": "/exact_graph_nodes/31162"},
            "current_Rp": {"path": MASTER_RAW, "pointer": "/exact_graph_nodes/263349"},
            "current_logRp": {
                "path": PHYSICAL_REPORT,
                "pointer": "/component_result/physical_pulse_bindings/current_logRp",
            },
            "current_logRp_row": {
                "path": PHYSICAL_REPORT,
                "pointer": "/component_result/definition_graph/appended_nodes/0",
            },
            "Gamma_imported_logRp": {
                "path": GAMMA_RAW,
                "pointer": "/function_bindings/native_binding/imported_logRp",
            },
            "Gamma_schedule_logRp": {
                "path": GAMMA_RAW,
                "pointer": "/function_bindings/schedule/logRp",
            },
        },
        "preserved_scope": preserved,
        "source_inspection_audit": {
            "new_source_consumer_calls": 0,
            "accepted_ancestor_operator_calls": 0,
            "old_math_consumers_called": 0,
            "numeric_values_computed": False,
            "graph_traversal": False,
            "chart_root_measure_endpoint_comparison": False,
            "Gamma_expectation_consumed": False,
        },
    }
    return result


def load_logrp_bridge_once(workspace_root: os.PathLike[str] | str) -> dict[str, Any]:
    """Authenticate the six frozen records, then hand them to the pure kernel."""

    try:
        records = {
            "geometry_result": _authenticated_document(
                workspace_root, GEOMETRY_RESULT, gzip_input=False
            ),
            "scale": _authenticated_document(workspace_root, SCALE_REPORT, gzip_input=False),
            "terminal": _authenticated_document(
                workspace_root, TERMINAL_REPORT, gzip_input=False
            ),
            "master": _authenticated_document(workspace_root, MASTER_RAW, gzip_input=True),
            "physical": _authenticated_document(
                workspace_root, PHYSICAL_REPORT, gzip_input=False
            ),
            "gamma": _authenticated_document(workspace_root, GAMMA_RAW, gzip_input=True),
            "gamma_check": _authenticated_document(workspace_root, GAMMA_CHECK, gzip_input=False),
        }
    except (OSError, SourceAuthenticationError, TypeError, ValueError) as error:
        return {
            "schema": SCHEMA,
            "status": "blocked",
            "passed": False,
            "source_only": True,
            "checks": {},
            "failures": ["authenticated_source_loading"],
            "reason": str(error),
            "rule": dict(LOG_EXP_RULE),
            "proof": {
                "conditional": True,
                "numeric_evaluation": False,
                "statement": "source authentication failed; no logRp identity is admitted",
            },
            "preserved_scope": _preserved_scope(None),
            "source_inspection_audit": {
                "new_source_consumer_calls": 0,
                "numeric_values_computed": False,
                "Gamma_expectation_consumed": False,
            },
        }
    return build_logrp_bridge(records)

