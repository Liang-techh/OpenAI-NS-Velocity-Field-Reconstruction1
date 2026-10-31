"""Source-reference bindings for the old log-radius seam.

This leaf is deliberately smaller than a source evaluator.  A runner supplies
two already authenticated/decoded mappings: the ``own`` archive and the
previous signed-seam result.  Only the named endpoint records are read.  An
old ``log_radius`` or ``query`` integer is an opaque source handle, never a
numeric argument.  Likewise, the previous ``F.R`` and ``F.R_Z`` records are
kept as opaque candidate expressions; no graph row is looked up and no source
function is called.

The returned objects are JSON-safe expression records for a later evaluator::

    D_native = source_at(old.log_radius, old.query, 11) - log(F.R)
    D        = D_native + rho_ref
    D_Z      = partial_Z_at_fixed_own_query(source_at(old.log_radius,
                    old.query, 11)) - F.R_Z/F.R + rho_ref_Z

``rho_ref`` and ``rho_ref_Z`` are caller-owned unit/reference-convention
functions.  They are intentionally unresolved rather than silently set to
zero.  The conditions needed for physical use are recorded as open
obligations, even when an audit packet carries a boolean evidence flag.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_logradius_source_binding_v2"
COMPONENT = SCHEMA
STATUS = "source-reference candidate"
CHARTS = ("main", "gap")
SIDES = ("main.upper", "gap.lower")
_SIDE_INFO = {
    "main.upper": ("main", "upper"),
    "gap.lower": ("gap", "lower"),
}


class LogRadiusSourceBindingInputError(ValueError):
    """Raised when a required authenticated record is absent or unsafe."""


def _expr(operator: str, *arguments: Any, **metadata: Any) -> dict[str, Any]:
    """Build one unevaluated JSON expression node."""

    result: dict[str, Any] = {
        "operator": operator,
        "arguments": list(arguments),
        "source_values_selected": False,
    }
    result.update(metadata)
    return result


def _literal(text: str) -> dict[str, Any]:
    """Return syntax, not a source value."""

    return {
        "operator": "exact_literal",
        "text": text,
        "source_values_selected": False,
        "node_id_is_not_a_numeric_value": True,
    }


def _json_safe(value: Any, *, depth: int = 0) -> Any:
    """Copy JSON-shaped metadata without evaluating or following source data.

    A candidate is expected to be a small JSON DAG record.  The depth and item
    limits are defensive boundaries: they also prevent an accidentally supplied
    graph/map from becoming a traversal performed by this leaf.
    """

    if depth > 32:
        raise LogRadiusSourceBindingInputError("metadata depth limit exceeded")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for index, (key, item) in enumerate(value.items()):
            if index >= 128:
                raise LogRadiusSourceBindingInputError("metadata mapping limit exceeded")
            out[str(key)] = _json_safe(item, depth=depth + 1)
        return out
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) > 128:
            raise LogRadiusSourceBindingInputError("metadata sequence limit exceeded")
        return [_json_safe(item, depth=depth + 1) for item in value]
    # A non-JSON opaque token is represented, not called or coerced to a value.
    return {"kind": "opaque_python_token", "python_type": type(value).__name__}


def _required(mapping: Mapping[str, Any], key: str, path: str) -> Any:
    """Read exactly one named path component and fail closed if it is absent."""

    try:
        value = mapping[key]
    except (KeyError, IndexError, TypeError) as exc:
        raise LogRadiusSourceBindingInputError(f"missing required path: {path}") from exc
    if value is None:
        raise LogRadiusSourceBindingInputError(f"null required path: {path}")
    return value


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise LogRadiusSourceBindingInputError(f"path is not a mapping: {path}")
    return value


def _handle_id(value: Any) -> Any:
    """Extract only a conventional id for diagnostics; it remains opaque."""

    if type(value) is int:
        return value
    if isinstance(value, Mapping):
        for key in ("source_node_id", "node_id", "node", "handle", "root"):
            candidate = value.get(key)
            if type(candidate) is int or isinstance(candidate, str):
                return candidate
    return None


def _namespace(value: Any, fallback: str) -> str:
    if isinstance(value, Mapping):
        for key in ("namespace", "source_namespace", "owner_namespace"):
            candidate = value.get(key)
            if isinstance(candidate, str) and candidate:
                return candidate
    return fallback


def _opaque_ref(value: Any, *, owner: str, namespace: str, role: str, path: str) -> dict[str, Any]:
    """Retain an opaque handle without allowing it to become arithmetic."""

    return {
        "operator": "opaque_source_ref",
        "owner": owner,
        "namespace": _namespace(value, namespace),
        "source_role": role,
        "source_path": path,
        "source_handle": _json_safe(value),
        "source_node_id": _handle_id(value),
        "node_id_is_not_a_numeric_value": True,
        "handle_is_opaque": True,
        "source_values_selected": False,
    }


def _candidate_ref(value: Any, *, side: str, field: str, path: str) -> dict[str, Any]:
    """Wrap a saved F candidate as an expression reference, never as a value."""

    candidate = _mapping(value, path)
    if candidate.get("source_values_selected") is True or candidate.get("source_value_selected") is True:
        raise LogRadiusSourceBindingInputError(
            f"saved candidate has selected source values: {path}"
        )
    copied = _json_safe(candidate)
    return _expr(
        "saved_source_candidate_ref",
        copied,
        owner="F",
        namespace=_namespace(candidate, f"prior_seam_result.endpoints.{side}"),
        source_role=f"F.{side}.{field}",
        source_path=path,
        candidate=copied,
        source_candidate=copied,
        candidate_is_expression=True,
        candidate_value_selected=False,
        node_id_is_not_a_numeric_value=True,
    )


def _exact_xi_11() -> dict[str, Any]:
    return {
        "operator": "exact_rational",
        "numerator": "11",
        "denominator": "1",
        "text": "11",
        "typed_exact_rational_11": True,
        "source_values_selected": False,
        "node_id_is_not_a_numeric_value": True,
    }


def _source_at(log_radius: Mapping[str, Any], query: Mapping[str, Any], *, side: str) -> dict[str, Any]:
    xi = _exact_xi_11()
    return _expr(
        "source_at",
        log_radius,
        query,
        xi,
        source_role=f"old.{side}.log_radius_at_query",
        expression_ref=log_radius,
        query_ref=query,
        endpoint_ref=xi,
        old_source_path=log_radius.get("source_path"),
        old_query_path=query.get("source_path"),
        source_executed=False,
        capture_safe_endpoint_proved=False,
        source_values_selected=False,
    )


def _rho_ref(name: str, formula: str) -> dict[str, Any]:
    """An explicit caller-owned unit/reference function, never zero."""

    return _expr(
        name,
        name=name,
        symbol=name,
        formula_text=formula,
        caller_supplied=True,
        unresolved=True,
        defaulted=False,
        zero_assumption=False,
        source_values_selected=False,
    )


def _audit_copy(audit: Mapping[str, Any]) -> dict[str, Any]:
    # Keep the prior audit as an external source reference.  It can contain
    # deeply nested expression records, and is not needed to build D or D_Z.
    # Re-copying it would traverse unrelated source expressions and could
    # silently truncate evidence; the prior authenticated receipt retains it.
    return {"source_path": "prior_seam_result['audit']", "retained_in_prior_receipt": True,
            "entire_audit_copied": False}


def _condition(audit: Mapping[str, Any], names: Sequence[str], statement: str) -> dict[str, Any]:
    """Record caller evidence without turning it into a proof."""

    evidence: Any = None
    evidence_key: str | None = None
    for name in names:
        value = audit.get(name)
        if value is not None:
            evidence = _json_safe(value)
            evidence_key = name
            break
    return {
        "condition": statement,
        "required_for_physical_use": True,
        "caller_evidence_key": evidence_key,
        "caller_evidence": evidence,
        "evidence_supplied_true": evidence is True,
        "proved": False,
        "unresolved": True,
        "source_values_selected": False,
    }


def _endpoint(
    own_root: Mapping[str, Any],
    prior_endpoints: Mapping[str, Any],
    audit: Mapping[str, Any],
    side: str,
) -> dict[str, Any]:
    chart, endpoint_name = _SIDE_INFO[side]
    own_fields = _mapping(
        _required(
            own_root,
            chart,
            f"own_archive['own_main_gap_end_fields']['{chart}']",
        ),
        f"own_archive['own_main_gap_end_fields']['{chart}']",
    )
    old_log = _opaque_ref(
        _required(
            own_fields,
            "log_radius",
            f"own_archive['own_main_gap_end_fields']['{chart}']['log_radius']",
        ),
        owner="old",
        namespace=f"own_main_gap_end_fields.{chart}",
        role=f"old.{side}.log_radius",
        path=f"own_archive['own_main_gap_end_fields']['{chart}']['log_radius']",
    )
    old_query = _opaque_ref(
        _required(
            own_fields,
            "query",
            f"own_archive['own_main_gap_end_fields']['{chart}']['query']",
        ),
        owner="old",
        namespace=f"own_main_gap_end_fields.{chart}",
        role=f"old.{side}.query",
        path=f"own_archive['own_main_gap_end_fields']['{chart}']['query']",
    )
    prior_side = _mapping(
        _required(prior_endpoints, side, f"prior_seam_result['endpoints']['{side}']"),
        f"prior_seam_result['endpoints']['{side}']",
    )
    f_r = _candidate_ref(
        _required(prior_side, "R", f"prior_seam_result['endpoints']['{side}']['R']"),
        side=side,
        field="R",
        path=f"prior_seam_result['endpoints']['{side}']['R']",
    )
    f_r_z = _candidate_ref(
        _required(prior_side, "R_Z", f"prior_seam_result['endpoints']['{side}']['R_Z']"),
        side=side,
        field="R_Z",
        path=f"prior_seam_result['endpoints']['{side}']['R_Z']",
    )

    old_source_at = _source_at(old_log, old_query, side=side)
    d_native = _expr(
        "subtract",
        old_source_at,
        _expr(
            "log",
            f_r,
            source_role=f"F.{side}.R.log",
            F_source_path=f"prior_seam_result['endpoints']['{side}']['R']",
            source_values_selected=False,
        ),
        formula_text="D_native=source_at(old.log_radius,old.query,11)-log(F.R)",
        old_source_path=old_log["source_path"],
        F_source_path=f_r["source_path"],
        source_values_selected=False,
    )
    rho = _rho_ref(
        "rho_ref",
        "rho_ref=log(old_reference_radius/F_reference_radius)",
    )
    d_physical = _expr(
        "add",
        d_native,
        rho,
        formula_text="D=D_native+rho_ref",
        rho_ref=rho,
        caller_unit_function_unresolved=True,
        source_values_selected=False,
    )
    partial = _expr(
        "partial_Z_at_fixed_own_query",
        old_source_at,
        fixed_query=old_query,
        source_role=f"old.{side}.log_radius.partial_Z_at_fixed_own_query",
        derivative_defaulted_zero=False,
        source_values_selected=False,
    )
    quotient = _expr(
        "quotient",
        f_r_z,
        f_r,
        formula_text="F.R_Z/F.R",
        denominator_positive_required=True,
        source_values_selected=False,
    )
    rho_z = _rho_ref(
        "rho_ref_Z",
        "rho_ref_Z=partial_Z log(old_reference_radius/F_reference_radius)",
    )
    d_z = _expr(
        "add",
        _expr(
            "subtract",
            partial,
            quotient,
            formula_text="partial_Z_at_fixed_own_query(source_at(old.log_radius,11))-F.R_Z/F.R",
            source_values_selected=False,
        ),
        rho_z,
        formula_text="D_Z=partial_Z_at_fixed_own_query(source_at(old.log_radius,11))-F.R_Z/F.R+rho_ref_Z",
        rho_ref_Z=rho_z,
        old_source_path=old_log["source_path"],
        F_source_path=f_r["source_path"],
        F_R_Z_source_path=f_r_z["source_path"],
        derivative_defaulted_zero=False,
        source_values_selected=False,
    )
    return {
        "side": side,
        "chart": chart,
        "endpoint": endpoint_name,
        "old": {
            "log_radius": old_log,
            "query": old_query,
            "source_values_selected": False,
        },
        "source_at_old_log_radius": old_source_at,
        "old_log_radius_source_at": old_source_at,
        "F": {
            "R": f_r,
            "R_Z": f_r_z,
            "source_values_selected": False,
        },
        "F_R": f_r,
        "F_R_Z": f_r_z,
        "rho_ref": rho,
        "rho_ref_Z": rho_z,
        "D_native": d_native,
        "D": d_physical,
        "D_Z": d_z,
        "conditions": {
            "old_log_radius_dimensionless_actual_source": _condition(
                audit,
                ("old_log_radius_is_dimensionless_actual_source", "old_logR_dimensionless_actual_source"),
                "old logR is an actual source log-radius and is dimensionless",
            ),
            "F_R_positive_same_reference_convention": _condition(
                audit,
                ("F_R_positive_same_reference_convention", "F_R_positive_and_same_reference_convention"),
                "F.R is positive and uses the same numerical reference convention",
            ),
            "edge_substitution_capture_safe": _condition(
                audit,
                ("edge_substitution_capture_safe", "capture_safe_edge_substitution"),
                "xi=11 edge substitution is capture-safe",
            ),
            "ordinary_Z1_source_row_semantics": _condition(
                audit,
                ("ordinary_Z1_source_row_semantics", "ordinary_Z1_semantics"),
                "ordinary Z1 is the intended source-row derivative",
            ),
            "F_R_Z_is_derivative": _condition(
                audit,
                ("F_R_Z_is_derivative", "F_R_Z_true_derivative"),
                "F.R_Z is truly the ordinary derivative of F.R",
            ),
        },
        "source_values_selected": False,
        "source_executed": False,
    }


def build(own_archive: Mapping[str, Any], prior_seam_result: Mapping[str, Any]) -> dict[str, Any]:
    """Build independent main/gap source-reference candidates.

    The two public inputs must already be authenticated/parsed by the caller.
    This function reads only the explicitly named paths documented in the
    module docstring and fails closed on any missing path or selected value.
    """

    own = _mapping(own_archive, "own_archive")
    prior = _mapping(prior_seam_result, "prior_seam_result")
    own_root = _mapping(
        _required(own, "own_main_gap_end_fields", "own_archive['own_main_gap_end_fields']"),
        "own_archive['own_main_gap_end_fields']",
    )
    prior_endpoints = _mapping(
        _required(prior, "endpoints", "prior_seam_result['endpoints']"),
        "prior_seam_result['endpoints']",
    )
    source_family = _json_safe(
        _required(prior, "source_family", "prior_seam_result['source_family']")
    )
    audit = _mapping(
        _required(prior, "audit", "prior_seam_result['audit']"),
        "prior_seam_result['audit']",
    )
    endpoints = {
        side: _endpoint(own_root, prior_endpoints, audit, side)
        for side in SIDES
    }
    d_native = {side: endpoints[side]["D_native"] for side in SIDES}
    d_physical = {side: endpoints[side]["D"] for side in SIDES}
    d_z = {side: endpoints[side]["D_Z"] for side in SIDES}
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": STATUS,
        "status_code": "source_reference_candidate",
        "source_reference_candidate": True,
        "constructed_means_expression_records_only": True,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "source_executed": False,
        "graph_accessed": False,
        "graph_traversed": False,
        "exact_graph_nodes_read": False,
        "gzip_read": False,
        "source_family": source_family,
        "prior_source_family": source_family,
        "prior_audit_evidence": _audit_copy(audit),
        "endpoints": endpoints,
        "D_native": d_native,
        "D": d_physical,
        "D_Z": d_z,
        "D_native_by_side": d_native,
        "D_by_side": d_physical,
        "D_Z_by_side": d_z,
        "independent_endpoint_sides": True,
        "main_gap_combined_or_cancelled": False,
        "old_errors_removed": False,
        "historical_error_terms_deleted": False,
        "rho_ref_defaulted_zero": False,
        "rho_ref_Z_defaulted_zero": False,
        "derivatives_defaulted_zero": False,
        "conditions_are_unresolved_obligations": True,
        "audit": {
            "prior_audit_preserved": True,
            "old_log_radius_and_F_candidates_retained": True,
            "old_errors_removed": False,
            "rho_ref_resolved": False,
            "rho_ref_Z_resolved": False,
            "capture_safe_edge_substitution_proved": False,
            "ordinary_Z1_source_row_semantics_proved": False,
            "F_R_positive_and_reference_convention_proved": False,
            "F_R_Z_derivative_proved": False,
            "source_values_computed": False,
        },
    }


build_source_reference_candidate = build
construct = build


__all__ = [
    "CHARTS",
    "COMPONENT",
    "LogRadiusSourceBindingInputError",
    "SCHEMA",
    "SIDES",
    "STATUS",
    "build",
    "build_source_reference_candidate",
    "construct",
]
