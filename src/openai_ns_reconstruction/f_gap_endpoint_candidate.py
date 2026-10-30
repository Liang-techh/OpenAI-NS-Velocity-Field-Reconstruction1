"""Candidate F gap endpoint source packet.

This source-only leaf authenticates both frozen prior manifests and every raw
file named by them before it reads either prior ``result.json``.  It consumes
only F ``leading_fields_gap`` rows and the saved historical-radius contract.
Saved ids remain opaque metadata.  The radius owner, capture-free endpoint
substitution, fixed-coordinate xi partials, and ordinary-Z transport remain
explicitly unresolved; no value, graph edge, quadrature, old helper, or old
error is consumed.  The emitted lower/upper eight-field, K, K_Z and
``-(upper-lower)`` objects are candidate expression DAGs, not proofs.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "f_gap_endpoint_candidate_v1"
COMPONENT = SCHEMA
ACTUAL_MANIFEST = Path("tmp/pulse_gap_actual_source_program_v1/source.json")
ACTUAL_MANIFEST_SHA256 = "fa460ab8417828dda6b49be7dec6d3fa824af2fe4ed170f23a8e74ff5233390f"
ACTUAL_RESULT = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
ACTUAL_RESULT_SHA256 = "e611cdc7e04940e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"
PF_MANIFEST = Path("tmp/pf_noncommon_moment_v1/source.json")
PF_MANIFEST_SHA256 = "7f72faa2af671e3e14ee1ed135b9cdd85153131e2144079beaed7733dc5f8060"
PF_RESULT = Path("tmp/pf_noncommon_moment_v1/result.json")
PF_RESULT_SHA256 = "fa5bfea33cfb082964772de2f113b54171c5942a55cdcdb75c9cd0eecde8ee20"
ORDERS = (0, 1, 2)
ENDPOINT_FIELDS = (
    "R_q", "R_squared_q", "R_q_squared", "R_squared_q_squared",
    "R_U_squared", "R_squared_U_squared", "q_plus_q_y", "R_q_y",
)
F_PRODUCER_POINTER = "$.Rp_native_C5_pulse_source_bindings"
F_LEADING_POINTER = F_PRODUCER_POINTER + ".leading_fields.gap"


class FGapEndpointCandidateInputError(ValueError):
    """A frozen manifest/result boundary is absent or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapEndpointCandidateInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _rel(path: Path) -> str:
    return path.as_posix()


def _expr(operator: str, *args: Any, **meta: Any) -> dict[str, Any]:
    result = {"operator": operator, "arguments": list(args), "source_values_selected": False}
    result.update(meta)
    return result


def _literal(value: str) -> dict[str, Any]:
    return {"kind": "exact_scalar_literal", "value": value, "source_values_selected": False}


def _formal(name: str, **meta: Any) -> dict[str, Any]:
    result = {"kind": "formal_source_symbol", "name": name, "source_values_selected": False}
    result.update(meta)
    return result


def _unresolved(name: str, role: str, reason: str) -> dict[str, Any]:
    return {
        "status": "unresolved", "kind": "unresolved_source_role", "name": name,
        "source_role": role, "reason": reason, "source_values_selected": False,
    }


def _read(root: Path, relative: Path) -> bytes:
    path = root / relative
    _require(path.is_file(), f"missing frozen input: {_rel(relative)}")
    return path.read_bytes()


def _manifest_path(root: Path, raw_path: Any) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "manifest path missing")
    path = Path(raw_path)
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    _require(path.is_relative_to(root.resolve()), f"manifest path escapes workspace: {raw_path}")
    _require(path.is_file(), f"manifest file missing: {raw_path}")
    return path


def _authenticate_manifest(root: Path, relative: Path, expected: str, label: str) -> dict[str, Any]:
    """Hash the manifest and every listed raw file before reading prior results."""

    raw = _read(root, relative)
    actual = _sha(raw)
    _require(actual == expected, f"{label}: manifest SHA-256 changed")
    try:
        manifest = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapEndpointCandidateInputError(f"{label}: invalid manifest JSON") from exc
    _require(isinstance(manifest, Mapping), f"{label}: manifest is not an object")
    files = manifest.get("files")
    _require(isinstance(files, list) and files, f"{label}: manifest files missing")
    seen: set[str] = set()
    authenticated: list[dict[str, Any]] = []
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"{label}: malformed files[{index}]")
        raw_path, expected_file = item.get("path"), item.get("sha256")
        _require(isinstance(expected_file, str) and len(expected_file) == 64,
                 f"{label}: malformed files[{index}] SHA")
        key = str(raw_path)
        _require(key not in seen, f"{label}: duplicate manifest path {key}")
        seen.add(key)
        path = _manifest_path(root, raw_path)
        actual_file = _sha(path.read_bytes())
        _require(actual_file == expected_file, f"{label}: raw SHA changed for {key}")
        authenticated.append({
            "manifest_path": key,
            "workspace_path": path.relative_to(root.resolve()).as_posix(),
            "expected_raw_sha256": expected_file,
            "raw_sha256": actual_file,
            "raw_sha256_authenticated": True,
        })
    return {
        "path": _rel(relative), "raw_sha256": actual, "schema": manifest.get("schema"),
        "files": authenticated, "all_raw_sha256_authenticated": True,
    }


def _prior(root: Path, relative: Path, expected: str, label: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    raw = _read(root, relative)
    actual = _sha(raw)
    _require(actual == expected, f"{label}: result SHA-256 changed")
    try:
        packet = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapEndpointCandidateInputError(f"{label}: invalid JSON") from exc
    _require(isinstance(packet, Mapping) and packet.get("status") == "returned",
             f"{label}: result is not returned")
    result = packet.get("result")
    _require(isinstance(result, Mapping), f"{label}: result payload is not an object")
    return result, {"path": _rel(relative), "raw_sha256": actual,
                    "raw_hash_authenticated_before_parse": True,
                    "source_values_selected": False}


def _actual_ref(value: Any, role: str, gates: list[str]) -> dict[str, Any]:
    """Admit only the exact direct-ref schema emitted by the F producer."""

    valid = (
        isinstance(value, Mapping)
        and value.get("status") == "direct_source_ref"
        and value.get("source_role") == role
        and type(value.get("source_node_id")) is int
        and value.get("node_id_is_not_a_numeric_value") is True
        and value.get("source_value_selected") is False
    )
    if valid:
        return {
            "status": "candidate_source_ref",
            "owner": "F",
            "namespace": "frame_Rp_native_C5_archive",
            "source_role": role,
            "source_node_id": value["source_node_id"],
            "node_id_is_not_a_numeric_value": True,
            "source_values_selected": False,
        }
    gates.append(f"{role}:direct_source_ref_contract_mismatch")
    return _unresolved(role, role, "F producer direct-source-ref contract did not match")


def _template_ref(value: Any, role: str, gates: list[str]) -> dict[str, Any]:
    """Admit only the historical contract's saved-opaque-ref schema."""

    valid = (
        isinstance(value, Mapping)
        and value.get("status") == "saved_opaque_ref"
        and value.get("source_role") == role
        and type(value.get("source_node_id")) is int
        and value.get("node_id_is_not_a_numeric_value") is True
        and value.get("source_value_selected") is False
    )
    if valid:
        return {
            "status": "saved_opaque_ref",
            "owner": "historical_template",
            "namespace": "historical_frame_gap_template",
            "source_role": role,
            "source_node_id": value["source_node_id"],
            "node_id_is_not_a_numeric_value": True,
            "source_values_selected": False,
        }
    gates.append(f"{role}:saved_opaque_ref_contract_mismatch")
    return _unresolved(role, role, "historical template ref contract did not match")


def _rows(record: Any, field: str, gates: list[str]) -> dict[str, Any]:
    role = f"{F_LEADING_POINTER}.{field}"
    group_valid = (
        isinstance(record, Mapping)
        and record.get("status") == "direct_source_rows"
        and record.get("source_role") == role
        and record.get("ordinary_Z_orders") == list(ORDERS)
        and record.get("source_values_selected") is False
    )
    if not group_valid:
        gates.append(f"{role}:direct_source_rows_contract_mismatch")
    # The producer's ordinary_Z0_to_Z2 field is the authoritative saved row
    # group.  all_saved_rows is deliberately not a fallback.
    raw_rows = record.get("ordinary_Z0_to_Z2") if isinstance(record, Mapping) else None
    if not isinstance(raw_rows, list):
        gates.append(f"{role}:ordinary_Z0_to_Z2_missing")
        raw_rows = []
    out = []
    for order in ORDERS:
        item = raw_rows[order] if order < len(raw_rows) else None
        if item is None:
            gates.append(f"{role}[{order}]:missing")
        row_role = f"{role}[{order}]"
        row = _actual_ref(item, row_role, gates) if group_valid else _unresolved(
            f"F.{field}.Z{order}", row_role,
            "group source-role contract did not match; row kept unresolved",
        )
        row["ordinary_Z_order"] = order
        out.append(row)
    return {"status": "direct_source_rows" if group_valid and all(x.get("status") == "candidate_source_ref" for x in out) else "unresolved",
            "source_role": role, "ordinary_Z_orders": list(ORDERS), "rows": out,
            "saved_row_count": len(raw_rows), "source_values_selected": False}


def _row_expr(rows: Mapping[str, Any], order: int, label: str) -> dict[str, Any]:
    values = rows.get("rows")
    if not isinstance(values, list) or order >= len(values):
        return _expr("unresolved_source_row", source_role=f"{label}[{order}]", ordinary_Z_order=order)
    return _expr("opaque_source_ref", dict(values[order]),
                 source_role=f"{label}[{order}]", node_id_is_not_a_numeric_value=True)


def _source_at(function_ref: Any, xi: Any, role: str) -> dict[str, Any]:
    return _expr("source_function_at", function_ref, xi, source_role=role,
                 capture_free_substitution_proved=False)


def _parameter(name: str, role: str, saved: Any | None = None) -> dict[str, Any]:
    result = _formal(name, source_role=role, parameter_owner_unproved=True)
    if saved is not None:
        result["saved_candidate_ref"] = saved
    return result


def _endpoint_parameters(endpoint: str, f_side: Mapping[str, Any], contract: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    geometry = f_side.get("geometry") if isinstance(f_side.get("geometry"), Mapping) else {}
    delta = _actual_ref(geometry.get("delta"), F_PRODUCER_POINTER + ".geometry.delta", gates)
    mu_saved = _actual_ref(geometry.get("mu"), F_PRODUCER_POINTER + ".geometry.mu", gates)
    radius_info = contract.get("historical_F_radius") if isinstance(contract.get("historical_F_radius"), Mapping) else {}
    template = radius_info.get("source_template") if isinstance(radius_info.get("source_template"), Mapping) else {}
    if not radius_info:
        gates.append("historical_F_radius:contract_missing")
    if not template:
        gates.append("historical_F_radius.source_template:missing")
    mu = _parameter("mu_F", "historical_F.mu", mu_saved)
    log_rp = _parameter("logRp_F", "historical_F.logRp",
                        _template_ref(template.get("log_radius"),
                                      "original_own_gap_field_template.log_radius", gates))
    mu_z = _unresolved("mu_F_Z", "historical_F.mu_Z", "radius contract leaves mu_Z unbound")
    log_rp_z = _unresolved("logRp_F_Z", "historical_F.logRp_Z", "radius contract leaves logRp_Z unbound")
    delta_z = _unresolved("delta_F_Z", "historical_F.delta_Z", "retain delta_Z in L_Z")
    if endpoint == "lower":
        xi = _literal("11")
        # The candidate coordinate is the exact literal xi_lower=11, so its
        # derivative is exactly zero.  This does not bind the actual F
        # endpoint owner or prove capture-free substitution; those gates stay
        # false below.
        xi_z = _literal("0")
    else:
        xi = _expr("subtract", _literal("13"), _expr("multiply", _literal("4"), mu),
                    formula_text="xi_upper=13-4*mu_F")
        xi_z = _expr("negate", _expr("multiply", _literal("4"), mu_z),
                     formula_text="xi_upper_Z=-4*mu_F_Z")
    ell = _expr("add", log_rp, _expr("divide", xi, mu),
                formula_text="ell=logRp_F+xi/mu_F")
    radius = _expr("exp", ell, formula_text="R=exp(logRp_F+xi/mu_F)", F_radius_owner_bound=False)
    radius_z = _expr("multiply", radius, _expr(
        "add", log_rp_z, _expr("divide", xi_z, mu),
        _expr("negate", _expr("divide", _expr("multiply", xi, mu_z), _expr("power", mu, _literal("2")))),
        formula_text="R_Z=R*(logRp_Z+xi_Z/mu-xi*mu_Z/mu^2)"),
        formula_text="R_Z=R*(logRp_Z+xi_Z/mu-xi*mu_Z/mu^2)", F_radius_owner_bound=False)
    Z = _formal("Z")
    L = _expr("subtract", _literal("1"), _expr("multiply", delta, _expr("power", Z, _literal("2"))),
               formula_text="L=1-delta*Z^2")
    L_z = _expr("negate", _expr("add", _expr("multiply", delta_z, _expr("power", Z, _literal("2"))),
                                      _expr("multiply", _literal("2"), delta, Z)),
                 formula_text="L_Z=-(delta_Z*Z^2+2*delta*Z)")
    return {"endpoint": endpoint, "xi": xi, "xi_Z": xi_z, "mu": mu, "mu_Z": mu_z,
            "logRp": log_rp, "logRp_Z": log_rp_z, "delta": delta, "delta_Z": delta_z,
            "ell": ell, "R": radius, "R_Z": radius_z, "L": L, "L_Z": L_z,
            "F_radius_owner_bound": False, "capture_free_endpoint_substitution_proved": False,
            "historical_radius_formula": radius_info.get("formula"),
            "source_values_selected": False}


def _field(field: str, rows: Mapping[str, Any], params: Mapping[str, Any], endpoint: str) -> dict[str, Any]:
    xi = params["xi"]
    value = _source_at(_row_expr(rows, 0, rows["source_role"]), xi, f"F.{field}.{endpoint}.Z0")
    fixed_z = _source_at(_row_expr(rows, 1, rows["source_role"]), xi, f"F.{field}.{endpoint}.fixed_Z1")
    partial = _unresolved(f"{field}_xi", f"F.{field}.{endpoint}.partial_xi",
                          "fixed-coordinate xi partial is absent from prior result")
    total_z = _expr("add", fixed_z, _expr("multiply", params["xi_Z"], partial),
                    formula_text=f"{field}_Z_total={field}_Z_fixed+xi_Z*partial_xi({field})")
    return {"value": value, "fixed_coordinate_Z1": fixed_z, "fixed_coordinate_xi_partial": partial,
            "total_Z1": total_z, "ordinary_Z_rows": dict(rows),
            "ordinary_Z_total_transport_proved": False, "source_values_selected": False}


def _prod_z1(values: Sequence[Any], derivs: Sequence[Any], text: str) -> dict[str, Any]:
    return _expr("add", *[_expr("multiply", *(derivs[i] if i == j else values[j] for j in range(len(values))))
                           for i in range(len(values))], formula_text=text)


def _power(value: Any, exponent: str) -> dict[str, Any]:
    return _expr("power", value, _literal(exponent))


def _endpoint(endpoint: str, f_side: Mapping[str, Any], rows: Mapping[str, Any],
              contract: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    params = _endpoint_parameters(endpoint, f_side, contract, gates)
    fields = {field: _field(field, rows[field], params, endpoint) for field in ("U", "Q", "Q_y")}
    R, R_z = params["R"], params["R_Z"]
    Q, Q_z = fields["Q"]["value"], fields["Q"]["total_Z1"]
    U, U_z = fields["U"]["value"], fields["U"]["total_Z1"]
    Q_y, Q_y_z = fields["Q_y"]["value"], fields["Q_y"]["total_Z1"]
    R2, R2_z = _power(R, "2"), _prod_z1((R, R), (R_z, R_z), "(R^2)_Z")
    Q2, Q2_z = _power(Q, "2"), _prod_z1((Q, Q), (Q_z, Q_z), "(Q^2)_Z")
    U2, U2_z = _power(U, "2"), _prod_z1((U, U), (U_z, U_z), "(U^2)_Z")
    values = {
        "R_q": _expr("multiply", R, Q),
        "R_squared_q": _expr("multiply", R2, Q),
        "R_q_squared": _expr("multiply", R, Q2),
        "R_squared_q_squared": _expr("multiply", R2, Q2),
        "R_U_squared": _expr("multiply", R, U2),
        "R_squared_U_squared": _expr("multiply", R2, U2),
        "q_plus_q_y": _expr("add", Q, Q_y),
        "R_q_y": _expr("multiply", R, Q_y),
    }
    derivatives = {
        "R_q": _prod_z1((R, Q), (R_z, Q_z), "(R*Q)_Z"),
        "R_squared_q": _prod_z1((R2, Q), (R2_z, Q_z), "(R^2*Q)_Z"),
        "R_q_squared": _prod_z1((R, Q2), (R_z, Q2_z), "(R*Q^2)_Z"),
        "R_squared_q_squared": _prod_z1((R2, Q2), (R2_z, Q2_z), "(R^2*Q^2)_Z"),
        "R_U_squared": _prod_z1((R, U2), (R_z, U2_z), "(R*U^2)_Z"),
        "R_squared_U_squared": _prod_z1((R2, U2), (R2_z, U2_z), "(R^2*U^2)_Z"),
        "q_plus_q_y": _expr("add", Q_z, Q_y_z, formula_text="(Q+Q_y)_Z"),
        "R_q_y": _prod_z1((R, Q_y), (R_z, Q_y_z), "(R*Q_y)_Z"),
    }
    eight = {name: {"candidate": values[name], "ordinary_Z1_candidate": derivatives[name],
                    "source_identity_proved": False, "source_values_selected": False}
             for name in ENDPOINT_FIELDS}
    Z = _formal("Z")
    A, A_z = _expr("multiply", R2, Q, formula_text="A=R^2*Q"), _prod_z1((R2, Q), (R2_z, Q_z), "A_Z")
    B, B_z = _expr("multiply", R2, Q2, formula_text="B=R^2*Q^2"), _prod_z1((R2, Q2), (R2_z, Q2_z), "B_Z")
    C, C_z = _expr("multiply", R2, U2, formula_text="C=R^2*U^2"), _prod_z1((R2, U2), (R2_z, U2_z), "C_Z")
    D, D_z = _expr("multiply", R, Q_y, formula_text="D=R*Q_y"), _prod_z1((R, Q_y), (R_z, Q_y_z), "D_Z")
    L, L_z = params["L"], params["L_Z"]
    K = _expr("add",
        _expr("divide", A, _expr("multiply", _literal("2"), L)),
        _expr("divide", B, _literal("4")),
        _expr("negate", _expr("divide", _expr("multiply", _power(Z, "2"), C), _power(L, "2"))),
        _expr("negate", D),
        formula_text="K=R^2*Q/(2*L)+R^2*Q^2/4-Z^2*R^2*U^2/L^2-R*Q_y",
        candidate=True, source_identity_proved=False)
    K_z = _expr("add",
        _expr("divide", A_z, _expr("multiply", _literal("2"), L)),
        _expr("negate", _expr("divide", _expr("multiply", A, L_z),
                                _expr("multiply", _literal("2"), _power(L, "2")))),
        _expr("divide", B_z, _literal("4")),
        _expr("negate", _expr("divide", _expr("add", _expr("multiply", _literal("2"), Z, C),
                                               _expr("multiply", _power(Z, "2"), C_z)), _power(L, "2"))),
        _expr("divide", _expr("multiply", _literal("2"), _power(Z, "2"), C, L_z), _power(L, "3")),
        _expr("negate", D_z),
        formula_text="K_Z=A_Z/(2L)-A*L_Z/(2L^2)+B_Z/4-(2Z*C+Z^2*C_Z)/L^2+2Z^2*C*L_Z/L^3-D_Z",
        candidate=True, source_identity_proved=False, delta_Z_retained=True)
    return {"endpoint": endpoint, "parameters": params, "source_rows": dict(rows), "fields": fields,
            "eight_fields": eight, "intermediate": {"A": A, "A_Z": A_z, "B": B, "B_Z": B_z,
            "C": C, "C_Z": C_z, "D": D, "D_Z": D_z}, "K": K, "K_Z": K_z,
            "candidate_not_proved": True, "source_values_selected": False}


def _contract_summary(contract: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    """Retain only sanitized radius-role metadata, never raw id-bearing trees."""

    info = contract.get("historical_F_radius")
    info = info if isinstance(info, Mapping) else {}
    template = info.get("source_template")
    template = template if isinstance(template, Mapping) else {}
    return {
        "schema": contract.get("schema"),
        "status": contract.get("status"),
        "formula": info.get("formula"),
        "owner_bound": False,
        "candidate_template_only": True,
        "source_template": {
            name: _template_ref(
                template.get(name),
                f"original_own_gap_field_template.{name}",
                gates,
            )
            for name in ("radius", "query", "log_radius")
        },
        "moving_endpoint_formulas": {
            "lower": "xi_lower=11; xi_lower_Z remains source-owned",
            "upper": "xi_upper=13-4*mu_F; xi_upper_Z=-4*mu_F_Z",
        },
        "mu_Z_source_role": "historical_F.mu_Z",
        "logRp_Z_source_role": "historical_F.logRp_Z",
        "source_values_selected": False,
    }


def build_candidate(repo_root: str | Path = ".", source: Mapping[str, Any] | None = None,
                    hashes: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Build the unevaluated candidate after complete prior-manifest auth."""

    root = Path(repo_root).resolve()
    actual_manifest = _authenticate_manifest(root, ACTUAL_MANIFEST, ACTUAL_MANIFEST_SHA256, "actual-source")
    pf_manifest = _authenticate_manifest(root, PF_MANIFEST, PF_MANIFEST_SHA256, "pf-noncommon")
    actual, actual_input = _prior(root, ACTUAL_RESULT, ACTUAL_RESULT_SHA256, "actual-source prior")
    pf, pf_input = _prior(root, PF_RESULT, PF_RESULT_SHA256, "pf-noncommon prior")
    _require(actual.get("schema") == "pulse_gap_actual_source_program_v1", "actual prior schema changed")
    _require(actual.get("source_only") is True and actual.get("source_values_selected") is not True,
             "actual prior source scope changed")
    actual_flags = actual.get("flags")
    _require(isinstance(actual_flags, Mapping), "actual prior source flags missing")
    for flag in ("numeric_values_computed", "source_graph_traversed",
                 "source_function_values_selected", "old_source_constructor_imported",
                 "old_source_checker_imported", "incident_pair_DAG_reentered",
                 "root_mapping_performed", "rewritten_node_mapping_consumed"):
        _require(actual_flags.get(flag) is False, f"actual prior source gate is not false: {flag}")
    sides = actual.get("sides")
    f_side = sides.get("F") if isinstance(sides, Mapping) and isinstance(sides.get("F"), Mapping) else {}
    contract = pf.get("historical_radius_contract") if isinstance(pf.get("historical_radius_contract"), Mapping) else {}
    packet = pf.get("source_packet")
    _require(isinstance(packet, Mapping), "PF prior source_packet missing")
    _require(packet.get("source_only") is True, "PF prior source-only gate changed")
    for flag in ("source_values_computed", "source_graph_traversed",
                 "numeric_use_admissible", "certified", "old_errors_removed"):
        _require(packet.get(flag) is False, f"PF prior source gate is not false: {flag}")
    _require(isinstance(contract, Mapping), "PF historical radius contract missing")
    for flag in ("source_values_computed", "graph_traversed", "admissible_for_numeric_use"):
        _require(contract.get(flag) is False, f"PF radius source gate is not false: {flag}")
    gates: list[str] = []
    if not contract:
        gates.append("historical_radius_contract:missing")
    elif contract.get("schema") != "pf_historical_radius_contract_v1":
        gates.append("historical_radius_contract:schema_changed")
    contract_summary = _contract_summary(contract, gates)
    leading = f_side.get("leading_fields_gap") if isinstance(f_side.get("leading_fields_gap"), Mapping) else {}
    rows = {field: _rows(leading.get(field), field, gates) for field in ("U", "Q", "Q_y")}
    lower = _endpoint("lower", f_side, rows, contract, gates)
    upper = _endpoint("upper", f_side, rows, contract, gates)
    direction = {"orientation": "lower_to_upper", "formula": "-(upper-lower)",
                 "K": _expr("negate", _expr("subtract", upper["K"], lower["K"]), formula_text="-(K_upper-K_lower)"),
                 "K_Z": _expr("negate", _expr("subtract", upper["K_Z"], lower["K_Z"]), formula_text="-(K_Z_upper-K_Z_lower)"),
                 "direction_identity_proved": False, "source_values_selected": False}
    retained = actual.get("retained_212_scope")
    if isinstance(retained, Mapping):
        retained = dict(retained)
        retained["old_errors_removed"] = False
        retained["source_values_selected"] = False
    else:
        retained = {"status": "unresolved", "source_values_selected": False}
        gates.append("retained_212_scope:missing")
    source_gates = {
        "manifest_all_raw_sha256_authenticated_before_prior_json_read": True,
        "F_leading_U_Q_Qy_source_rows_present": all(item["status"] == "direct_source_rows" for item in rows.values()),
        "historical_F_radius_owner_bound": False, "historical_F_endpoint_owner_bound": False,
        "fixed_coordinate_xi_partials_bound": False, "ordinary_Z_total_transport_proved": False,
        "capture_free_endpoint_substitution": False, "delta_Z_source_bound": False,
        "eight_field_source_identity_proved": False, "K_source_identity_proved": False,
        "K_Z_source_identity_proved": False, "endpoint_direction_identity_proved": False,
        "numeric_use_admissible": False, "old_errors_removed": False,
        "ledger_212_213_183_changed": False,
    }
    for gate in gates:
        source_gates[gate] = False
    return {
        "schema": SCHEMA, "component": COMPONENT,
        "status": "candidate_source_owned_F_gap_endpoint_unproved",
        "source_only": True, "source_values_computed": False, "source_values_selected": False,
        "graph_traversed": False, "certified": False, "old_leaf_called": False,
        "old_helper_imported": False, "quadrature_called": False,
        "old_errors_removed": False, "error_transport_certified": False,
        "provenance": {"runner_source_supplied": isinstance(source, Mapping),
                       "runner_hashes_supplied": isinstance(hashes, Mapping),
                       "actual_manifest": actual_manifest, "pf_noncommon_manifest": pf_manifest,
                       "prior_results": {"actual_source": actual_input, "pf_noncommon": pf_input},
                       "all_manifest_raw_sha256_authenticated_before_prior_json_read": True},
        "F_source_roles": {"owner": "F", "namespace": "frame_Rp_native_C5_archive",
                           "leading_fields_gap": rows,
                           "saved_geometry": {name: _actual_ref(
                               f_side.get("geometry", {}).get(name) if isinstance(f_side.get("geometry"), Mapping) else None,
                               F_PRODUCER_POINTER + f".geometry.{name}", []) for name in ("mu", "delta")},
                           "source_values_selected": False},
        "historical_radius_contract": contract_summary,
        "endpoints": {"lower": lower, "upper": upper, "ordinary_Z_orders": list(ORDERS),
                      "source_values_selected": False},
        "direction": direction, "retained_212_scope": retained,
        "retained_legacy": {"physical_pulse_corrections_preserved": True,
                            "independent_reference_offset_preserved": True,
                            "old_errors_removed": False, "source_values_selected": False},
        "source_gates": source_gates,
        "unresolved": sorted(set(gates + ["F historical radius owner", "F finite endpoint owner",
                                           "fixed-coordinate xi partials", "capture-free endpoint substitution",
                                           "ordinary-Z total derivative transport", "certified signed error transport"])),
        "admissible_for_numeric_use": False, "full_goal_complete": False,
    }


def build_f_gap_endpoint_source_candidate(root: str | Path = ".", source: Mapping[str, Any] | None = None,
                                          hashes: Mapping[str, str] | None = None) -> dict[str, Any]:
    return build_candidate(root, source=source, hashes=hashes)


def _fixture_fields(z: Any, xi: Any) -> dict[str, Any]:
    return {"Q": 3 + 2*z + 5*xi + z*xi,
            "U": 1 + z + 2*xi + z*xi**2,
            "Q_y": 4 + 3*z + xi**2 + z*xi}


def _fixture_eight(R: Any, Q: Any, U: Any, Q_y: Any) -> dict[str, Any]:
    return {"R_q": R*Q, "R_squared_q": R**2*Q, "R_q_squared": R*Q**2,
            "R_squared_q_squared": R**2*Q**2, "R_U_squared": R*U**2,
            "R_squared_U_squared": R**2*U**2, "q_plus_q_y": Q+Q_y,
            "R_q_y": R*Q_y}


def _fixture_eight_z(R: Any, R_z: Any, Q: Any, Q_z: Any, U: Any, U_z: Any,
                     Q_y: Any, Q_y_z: Any) -> dict[str, Any]:
    return {"R_q": R_z*Q + R*Q_z,
            "R_squared_q": 2*R*R_z*Q + R**2*Q_z,
            "R_q_squared": R_z*Q**2 + 2*R*Q*Q_z,
            "R_squared_q_squared": 2*R*R_z*Q**2 + 2*R**2*Q*Q_z,
            "R_U_squared": R_z*U**2 + 2*R*U*U_z,
            "R_squared_U_squared": 2*R*R_z*U**2 + 2*R**2*U*U_z,
            "q_plus_q_y": Q_z + Q_y_z, "R_q_y": R_z*Q_y + R*Q_y_z}


def _fixture_json_safe(value: Any) -> Any:
    """Serialize SymPy expressions without making them source values."""

    if isinstance(value, Mapping):
        return {str(key): _fixture_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_fixture_json_safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def synthetic_fixture() -> dict[str, Any]:
    """Independent exact nonzero-moving-endpoint fixture; never run on import."""

    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover
        return {"schema": "f_gap_endpoint_candidate_synthetic_v1", "checks_passed": False,
                "synthetic_only": True, "source_values_computed": False,
                "error": f"sympy unavailable: {exc}"}
    z, x = sp.symbols("z x")
    mu, log_rp = 2 + z/sp.Integer(5), 1 + z/sp.Integer(7)
    delta = sp.Rational(1, 3) + z/sp.Integer(11)
    mu_z, log_rp_z, delta_z = sp.diff(mu, z), sp.diff(log_rp, z), sp.diff(delta, z)
    checks = {"nonzero_mu_Z": mu_z != 0, "nonzero_logRp_Z": log_rp_z != 0,
              "nonzero_delta_Z": delta_z != 0}
    data = {}
    for name, xi in (("lower", 11 + 2*z), ("upper", 13 - 4*mu)):
        xi_z = sp.diff(xi, z)
        checks[f"moving_{name}_endpoint_nonzero"] = xi_z != 0
        raw = _fixture_fields(z, x)
        values = {key: value.subs(x, xi) for key, value in raw.items()}
        fixed_z = {key: sp.diff(value, z).subs(x, xi) for key, value in raw.items()}
        partial_xi = {key: sp.diff(value, x).subs(x, xi) for key, value in raw.items()}
        total_z = {key: sp.simplify(fixed_z[key] + xi_z*partial_xi[key]) for key in raw}
        R = sp.exp(log_rp + xi/mu)
        R_z = R*(log_rp_z + xi_z/mu - xi*mu_z/mu**2)
        L, L_z = 1-delta*z**2, -(delta_z*z**2 + 2*delta*z)
        Q, U, Q_y = values["Q"], values["U"], values["Q_y"]
        Q_z, U_z, Q_y_z = total_z["Q"], total_z["U"], total_z["Q_y"]
        A, B, C, D = R**2*Q, R**2*Q**2, R**2*U**2, R*Q_y
        A_z = 2*R*R_z*Q + R**2*Q_z
        B_z = 2*R*R_z*Q**2 + 2*R**2*Q*Q_z
        C_z = 2*R*R_z*U**2 + 2*R**2*U*U_z
        D_z = R_z*Q_y + R*Q_y_z
        K = A/(2*L) + B/4 - z**2*C/L**2 - D
        K_z = (A_z/(2*L) - A*L_z/(2*L**2) + B_z/4
               - (2*z*C + z**2*C_z)/L**2 + 2*z**2*C*L_z/L**3 - D_z)
        eight, eight_z = _fixture_eight(R, Q, U, Q_y), _fixture_eight_z(
            R, R_z, Q, Q_z, U, U_z, Q_y, Q_y_z)
        data[name] = {"xi_Z": xi_z, "mu_Z": mu_z, "logRp_Z": log_rp_z,
                      "delta_Z": delta_z, "K": K, "K_Z": K_z}
        checks[f"{name}_radius_chain_rule"] = sp.simplify(
            R_z - sp.diff(R, z)) == 0
        checks[f"{name}_total_Z_chain_rule"] = all(
            sp.simplify(total_z[key] - sp.diff(values[key], z)) == 0 for key in raw)
        checks[f"{name}_K_Z_exact"] = sp.simplify(K_z - sp.diff(K, z)) == 0
        checks[f"{name}_eight_fields_exact"] = all(
            sp.simplify(eight_z[key] - sp.diff(eight[key], z)) == 0
            for key in ENDPOINT_FIELDS)
    direction = data["lower"]["K"] - data["upper"]["K"]
    direction_z = data["lower"]["K_Z"] - data["upper"]["K_Z"]
    checks["direction_minus_upper_minus_lower"] = sp.simplify(
        -(data["upper"]["K"] - data["lower"]["K"]) - direction) == 0
    checks["direction_Z_matches"] = sp.simplify(
        -(data["upper"]["K_Z"] - data["lower"]["K_Z"]) - direction_z) == 0
    checks["direction_nonzero_at_zero"] = sp.simplify(direction.subs(z, 0)) != 0
    checks = {name: bool(value) for name, value in checks.items()}
    return {"schema": "f_gap_endpoint_candidate_synthetic_v1", "checks": checks,
            "checks_passed": bool(all(checks.values())), "synthetic_only": True,
            "exact_arithmetic": True, "source_values_computed": False,
            "fixture_scope": "independent_formula_reimplementation_only",
            "emitted_candidate_dag_evaluated": False,
            "endpoint_data": _fixture_json_safe(data)}


def verify_synthetic_fixture() -> dict[str, Any]:
    return synthetic_fixture()


def recover(repo_root: str | Path = ".", source: Mapping[str, Any] | None = None,
            hashes: Mapping[str, str] | None = None) -> dict[str, Any]:
    return build_candidate(repo_root, source=source, hashes=hashes)


__all__ = ["ACTUAL_MANIFEST", "ACTUAL_MANIFEST_SHA256", "ACTUAL_RESULT", "ACTUAL_RESULT_SHA256",
           "COMPONENT", "ENDPOINT_FIELDS", "FGapEndpointCandidateInputError", "ORDERS",
           "PF_MANIFEST", "PF_MANIFEST_SHA256", "PF_RESULT", "PF_RESULT_SHA256", "SCHEMA",
           "build_candidate", "build_f_gap_endpoint_source_candidate", "recover",
           "synthetic_fixture", "verify_synthetic_fixture"]
