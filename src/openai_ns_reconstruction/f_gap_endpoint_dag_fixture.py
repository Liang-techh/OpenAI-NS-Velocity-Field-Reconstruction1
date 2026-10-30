"""Independent synthetic evaluator for the emitted F-gap endpoint DAG.

This module is deliberately a consumer-side fixture.  It does not open the
published result, import the candidate producer, walk an ancestor graph, or
turn an opaque source node id into a value.  A one-shot authenticated runner
passes the already-frozen ``prior_result`` to :func:`verify_emitted_candidate`.

The evaluator has two independent paths:

* ``_evaluate_emitted_dag`` is a small, allow-listed interpreter for the
  expression objects emitted by ``f_gap_endpoint_candidate_v2``.  It obtains
  all scalar values, ordinary-Z rows, and xi partials through a caller-owned
  synthetic value table and source-role map.
* ``_independent_sympy_expectation`` reconstructs the toy functions and takes
  fresh SymPy derivatives.  It never reads or rewrites the emitted expression
  tree.

The fixture is therefore an algebra/transport check only.  Its successful
synthetic comparison must not be interpreted as an F-owner, capture-free,
ordinary-Z transport, source identity, numerical, or error certificate.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping
from typing import Any


SCHEMA = "f_gap_endpoint_dag_fixture_v1"
COMPONENT = "f_gap_endpoint_dag_fixture"
F_PRODUCER_POINTER = "$.Rp_native_C5_pulse_source_bindings"
F_LEADING_POINTER = F_PRODUCER_POINTER + ".leading_fields.gap"
ENDPOINT_FIELDS = (
    "R_q",
    "R_squared_q",
    "R_q_squared",
    "R_squared_q_squared",
    "R_U_squared",
    "R_squared_U_squared",
    "q_plus_q_y",
    "R_q_y",
)
FIELD_NAMES = ("U", "Q", "Q_y")
ENDPOINT_NAMES = ("lower", "upper")
ORDERS = (0, 1, 2)


class FGapEndpointFixtureInputError(ValueError):
    """Raised when a candidate DAG, role map, or synthetic input is unsafe."""


class _EvaluationContext:
    """Private state for one allow-listed DAG evaluation."""

    def __init__(self, sp: Any, values: Mapping[str, Any], refs: Mapping[str, Any]):
        self.sp = sp
        self.values = values
        self.refs = refs
        self.xi_symbol = values["xi_symbol"]
        self.used_roles: set[str] = set()
        self.used_operator_counts: dict[str, int] = {}
        self.source_function_calls = 0
        self.opaque_source_refs_seen = 0
        self.opaque_ids_used_as_values = False
        self.current_endpoint: str | None = None
        self.current_xi: Any | None = None
        self.endpoint_xis: dict[str, Any] = {}

    def count_operator(self, operator: str) -> None:
        self.used_operator_counts[operator] = self.used_operator_counts.get(operator, 0) + 1


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapEndpointFixtureInputError(message)


def _sympy() -> Any:
    """Load SymPy lazily so importing this source-only fixture has no effect."""

    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise FGapEndpointFixtureInputError(f"SymPy is required for the fixture: {exc}") from exc
    return sp


def _is_int(value: Any) -> bool:
    return type(value) is int


def _to_scalar(value: Any, sp: Any, label: str) -> Any:
    """Admit scalar synthetic data, never a mapping carrying an opaque id."""

    _require(not isinstance(value, (Mapping, list, tuple, set)), f"{label}: scalar expected")
    try:
        converted = sp.sympify(value)
    except Exception as exc:  # pragma: no cover - defensive malformed input path
        raise FGapEndpointFixtureInputError(f"{label}: not a SymPy scalar") from exc
    # Do not rely on optional/top-level SymPy class exports (notably
    # ``sympy.Relational`` is not present in every supported 1.x build).
    _require(getattr(converted, "is_Matrix", False) is not True,
             f"{label}: scalar expected")
    _require(getattr(converted, "is_Relational", False) is not True,
             f"{label}: scalar expected")
    return converted


def _default_synthetic_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the independent nonzero moving-endpoint synthetic data.

    The source rows are data supplied to the evaluator, not values recovered
    from a source graph.  Row 0 is a value function, row 1 is a fixed-
    coordinate ordinary-Z1 derivative, and row 2 is an unused ordinary-Z2
    derivative retained to expose the distinction.
    """

    sp = _sympy()
    z, xi = sp.symbols("z xi")
    mu = sp.Integer(2) + z / sp.Integer(5)
    log_rp = sp.Integer(1) + z / sp.Integer(7)
    delta = sp.Rational(1, 3) + z / sp.Integer(11)
    functions = {
        "Q": 3 + 2 * z + 5 * xi + z * xi,
        "U": 1 + z + 2 * xi + z * xi**2,
        "Q_y": 4 + 3 * z + xi**2 + z * xi,
    }
    ordinary_z1 = {name: sp.diff(expr, z) for name, expr in functions.items()}
    ordinary_z2 = {name: sp.diff(expr, z, 2) for name, expr in functions.items()}
    xi_partials = {name: sp.diff(expr, xi) for name, expr in functions.items()}

    values: dict[str, Any] = {
        "schema": "f_gap_endpoint_synthetic_inputs_v1",
        "z_symbol": z,
        "xi_symbol": xi,
        # The emitted candidate's formal Z is the ordinary derivative
        # coordinate z in this exact fixture.
        "Z": z,
        "scalars": {
            "historical_F.mu": mu,
            "historical_F.mu_Z": sp.diff(mu, z),
            "historical_F.logRp": log_rp,
            "historical_F.logRp_Z": sp.diff(log_rp, z),
            F_PRODUCER_POINTER + ".geometry.delta": delta,
            "historical_F.delta_Z": sp.diff(delta, z),
            F_PRODUCER_POINTER + ".geometry.mu": mu,
        },
        "field_functions": functions,
        "field_ordinary_Z1": ordinary_z1,
        "field_ordinary_Z2": ordinary_z2,
        "field_xi_partials": xi_partials,
    }
    return values, _default_source_role_ref_map()


def _default_source_role_ref_map() -> dict[str, dict[str, Any]]:
    """Return role-to-semantic-key metadata; no values or ids are stored here."""

    refs: dict[str, dict[str, Any]] = {}
    for field in FIELD_NAMES:
        for order in ORDERS:
            semantics = {
                0: "field_value_function",
                1: "fixed_coordinate_ordinary_Z1",
                2: "fixed_coordinate_ordinary_Z2",
            }
            refs[f"{F_LEADING_POINTER}.{field}[{order}]"] = {
                "kind": "field_row",
                "field": field,
                "ordinary_Z_order": order,
                "semantic": semantics[order],
            }
        for endpoint in ENDPOINT_NAMES:
            refs[f"F.{field}.{endpoint}.partial_xi"] = {
                "kind": "field_xi_partial",
                "field": field,
                "semantic": "fixed_coordinate_xi_partial",
            }

    scalar_roles = {
        "historical_F.mu": "historical_F.mu",
        "historical_F.mu_Z": "historical_F.mu_Z",
        "historical_F.logRp": "historical_F.logRp",
        "historical_F.logRp_Z": "historical_F.logRp_Z",
        F_PRODUCER_POINTER + ".geometry.delta": F_PRODUCER_POINTER + ".geometry.delta",
        "historical_F.delta_Z": "historical_F.delta_Z",
        F_PRODUCER_POINTER + ".geometry.mu": F_PRODUCER_POINTER + ".geometry.mu",
    }
    for role, scalar_key in scalar_roles.items():
        refs[role] = {
            "kind": "scalar",
            "scalar_key": scalar_key,
            "semantic": "synthetic_scalar",
        }
    return refs


def build_synthetic_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    """Return default synthetic values and role metadata for a one-shot runner.

    The tuple contains SymPy expressions and is intentionally not a result
    packet.  :func:`verify_emitted_candidate` converts all report material to
    JSON-safe strings before returning it.
    """

    return _default_synthetic_inputs()


def _normalise_inputs(
    synthetic_values: Mapping[str, Any] | None,
    source_role_ref_map: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any], Any]:
    sp = _sympy()
    if synthetic_values is None and source_role_ref_map is None:
        values, refs = _default_synthetic_inputs()
    else:
        _require(isinstance(synthetic_values, Mapping), "synthetic_values must be supplied with source_role_ref_map")
        values = dict(synthetic_values)
        if source_role_ref_map is None:
            _, refs = _default_synthetic_inputs()
        else:
            _require(isinstance(source_role_ref_map, Mapping), "source_role_ref_map must be an object")
            refs = dict(source_role_ref_map)
    _validate_synthetic_inputs(values, sp)
    _validate_role_map(refs)
    return values, refs, sp


def _validate_synthetic_inputs(values: Mapping[str, Any], sp: Any) -> None:
    required = (
        "Z",
        "z_symbol",
        "xi_symbol",
        "scalars",
        "field_functions",
        "field_ordinary_Z1",
        "field_ordinary_Z2",
        "field_xi_partials",
    )
    for name in required:
        _require(name in values, f"synthetic_values.{name}: missing")
    _require(isinstance(values["scalars"], Mapping), "synthetic_values.scalars: object expected")
    _require(isinstance(values["field_functions"], Mapping), "synthetic_values.field_functions: object expected")
    _require(isinstance(values["field_ordinary_Z1"], Mapping), "synthetic_values.field_ordinary_Z1: object expected")
    _require(isinstance(values["field_ordinary_Z2"], Mapping), "synthetic_values.field_ordinary_Z2: object expected")
    _require(isinstance(values["field_xi_partials"], Mapping), "synthetic_values.field_xi_partials: object expected")
    _to_scalar(values["Z"], sp, "synthetic_values.Z")
    z_symbol = values["z_symbol"]
    _require(getattr(z_symbol, "is_Symbol", False) is True, "synthetic_values.z_symbol must be a SymPy symbol")
    _require(bool(sp.simplify(values["Z"] - z_symbol) == 0),
             "synthetic_values.Z must be the ordinary z_symbol used for derivatives")
    xi_symbol = values["xi_symbol"]
    _require(getattr(xi_symbol, "is_Symbol", False) is True, "synthetic_values.xi_symbol must be a SymPy symbol")
    for role in (
        "historical_F.mu",
        "historical_F.mu_Z",
        "historical_F.logRp",
        "historical_F.logRp_Z",
        F_PRODUCER_POINTER + ".geometry.delta",
        "historical_F.delta_Z",
    ):
        _require(role in values["scalars"], f"synthetic_values.scalars.{role}: missing")
        _to_scalar(values["scalars"][role], sp, f"synthetic_values.scalars.{role}")
    for field in FIELD_NAMES:
        for key in ("field_functions", "field_ordinary_Z1", "field_ordinary_Z2", "field_xi_partials"):
            table = values[key]
            _require(field in table, f"synthetic_values.{key}.{field}: missing")
            _to_scalar(table[field], sp, f"synthetic_values.{key}.{field}")


def _validate_role_map(refs: Mapping[str, Any]) -> None:
    _require(isinstance(refs, Mapping), "source_role_ref_map must be an object")
    for role, entry in refs.items():
        _require(isinstance(role, str) and role, "source_role_ref_map keys must be nonempty roles")
        _require(isinstance(entry, Mapping), f"source_role_ref_map[{role}]: object expected")
        # Values belong in synthetic_values.  An opaque id may be copied as
        # metadata, but it can never select the value table.
        if "source_node_id" in entry:
            _require(_is_int(entry["source_node_id"]), f"{role}: malformed opaque source_node_id")
            _require(entry.get("node_id_is_not_a_numeric_value") is True,
                     f"{role}: opaque id must be explicitly nonnumeric")
        for forbidden in ("value", "numeric_value", "expression", "source_value"):
            _require(forbidden not in entry, f"{role}: values must be outside source_role_ref_map")
        kind = entry.get("kind")
        _require(kind in {"scalar", "field_row", "field_xi_partial"}, f"{role}: unsupported role-map kind")
        semantic = entry.get("semantic")
        _require(isinstance(semantic, str) and semantic, f"{role}: semantic row meaning missing")
        if kind == "scalar":
            _require(semantic == "synthetic_scalar", f"{role}: scalar semantic mismatch")
            _require(isinstance(entry.get("scalar_key"), str) and entry["scalar_key"],
                     f"{role}: scalar_key missing")
            _require("node_id" not in entry["scalar_key"].lower(),
                     f"{role}: scalar_key cannot disguise an opaque node id")
        elif kind == "field_row":
            _require(entry.get("field") in FIELD_NAMES, f"{role}: unknown field row")
            _require(_is_int(entry.get("ordinary_Z_order")), f"{role}: ordinary_Z_order must be int")
            order = entry["ordinary_Z_order"]
            _require(order in ORDERS, f"{role}: ordinary_Z_order out of range")
            _require(role == f"{F_LEADING_POINTER}.{entry['field']}[{order}]",
                     f"{role}: source row role does not match field/order")
            expected = {
                0: "field_value_function",
                1: "fixed_coordinate_ordinary_Z1",
                2: "fixed_coordinate_ordinary_Z2",
            }[order]
            _require(semantic == expected, f"{role}: row semantic does not match ordinary_Z_order")
        else:
            _require(entry.get("field") in FIELD_NAMES, f"{role}: unknown xi-partial field")
            _require(role in {
                f"F.{entry['field']}.lower.partial_xi",
                f"F.{entry['field']}.upper.partial_xi",
            }, f"{role}: xi-partial role does not match field/endpoint")
            _require(semantic == "fixed_coordinate_xi_partial", f"{role}: xi-partial semantic mismatch")


def _opaque_ref_from_node(node: Any, ctx: _EvaluationContext) -> str:
    _require(isinstance(node, Mapping), "source_function_at: opaque source ref object expected")
    _require(node.get("operator") == "opaque_source_ref", "source_function_at: malformed ref operator")
    arguments = node.get("arguments")
    _require(isinstance(arguments, list) and len(arguments) == 1,
             "opaque_source_ref: exactly one metadata argument required")
    ref = arguments[0]
    _require(isinstance(ref, Mapping), "opaque_source_ref: metadata object expected")
    _require(ref.get("status") in {"candidate_source_ref", "saved_opaque_ref"},
             "opaque_source_ref: unsupported ref status")
    _require(_is_int(ref.get("source_node_id")), "opaque_source_ref: source_node_id must be opaque int")
    _require(ref.get("node_id_is_not_a_numeric_value") is True,
             "opaque_source_ref: source_node_id cannot be used as a numeric value")
    _require(ref.get("source_values_selected") is False,
             "opaque_source_ref: source values must remain unselected")
    role = ref.get("source_role")
    _require(isinstance(role, str) and role, "opaque_source_ref: source_role missing")
    ctx.opaque_source_refs_seen += 1
    # Only the semantic role leaves this function.  In particular, the id is
    # not returned and cannot enter arithmetic.
    return role


def _raw_opaque_role(node: Mapping[str, Any], ctx: _EvaluationContext) -> str:
    """Validate a direct source-role ref used by the candidate geometry leaf."""

    _require(node.get("status") in {"candidate_source_ref", "saved_opaque_ref"},
             "direct source ref: unsupported ref status")
    _require(_is_int(node.get("source_node_id")),
             "direct source ref: source_node_id must be opaque int")
    _require(node.get("node_id_is_not_a_numeric_value") is True,
             "direct source ref: source_node_id cannot be numeric")
    _require(node.get("source_values_selected") is False,
             "direct source ref: source values must remain unselected")
    role = node.get("source_role")
    _require(isinstance(role, str) and role, "direct source ref: source_role missing")
    ctx.opaque_source_refs_seen += 1
    return role


def _role_value(role: str, ctx: _EvaluationContext) -> Any:
    _require(role in ctx.refs, f"source role has no caller-provided synthetic mapping: {role}")
    entry = ctx.refs[role]
    _require(isinstance(entry, Mapping), f"source role mapping malformed: {role}")
    ctx.used_roles.add(role)
    if "source_node_id" in entry:
        _require(_is_int(entry["source_node_id"]), f"{role}: malformed opaque source_node_id")
        _require(entry.get("node_id_is_not_a_numeric_value") is True,
                 f"{role}: opaque id must not be numeric")
    for forbidden in ("value", "numeric_value", "expression", "source_value"):
        _require(forbidden not in entry, f"{role}: source-role map cannot carry a value")
    kind = entry.get("kind")
    if kind == "scalar":
        key = entry.get("scalar_key")
        _require(isinstance(key, str), f"{role}: scalar_key malformed")
        _require("node_id" not in key.lower(), f"{role}: scalar_key cannot use an opaque node id")
        _require(key in ctx.values["scalars"], f"{role}: scalar synthetic value missing: {key}")
        return _to_scalar(ctx.values["scalars"][key], ctx.sp, f"synthetic scalar {key}")
    if kind == "field_xi_partial":
        field = entry.get("field")
        _require(field in FIELD_NAMES, f"{role}: xi-partial field malformed")
        partial = _to_scalar(ctx.values["field_xi_partials"][field], ctx.sp,
                             f"synthetic xi partial {field}")
        # The emitted unresolved partial is endpoint-specific by source role;
        # consume its caller-supplied function at the current endpoint, while
        # keeping the row-1 derivative fixed-coordinate by construction.
        # Direction roots inline both endpoint trees after the endpoint loop;
        # use the endpoint encoded in the source role rather than whichever
        # endpoint happened to be evaluated last.
        role_parts = role.split(".")
        role_endpoint = role_parts[2] if len(role_parts) >= 4 else None
        substitution_xi = ctx.endpoint_xis.get(role_endpoint, ctx.current_xi)
        if substitution_xi is not None:
            partial = ctx.sp.sympify(partial).subs(ctx.xi_symbol, substitution_xi)
        return partial
    raise FGapEndpointFixtureInputError(
        f"{role}: a field row is only valid inside source_function_at"
    )


def _row_value(role: str, ctx: _EvaluationContext, xi: Any) -> Any:
    _require(role in ctx.refs, f"source row has no caller-provided mapping: {role}")
    entry = ctx.refs[role]
    _require(isinstance(entry, Mapping) and entry.get("kind") == "field_row",
             f"{role}: source_function_at requires a field-row mapping")
    field = entry.get("field")
    order = entry.get("ordinary_Z_order")
    _require(field in FIELD_NAMES and _is_int(order) and order in ORDERS,
             f"{role}: malformed field row mapping")
    expected_semantic = {
        0: "field_value_function",
        1: "fixed_coordinate_ordinary_Z1",
        2: "fixed_coordinate_ordinary_Z2",
    }[order]
    _require(entry.get("semantic") == expected_semantic,
             f"{role}: field source-row semantic mismatch")
    ctx.used_roles.add(role)
    table_name = {0: "field_functions", 1: "field_ordinary_Z1", 2: "field_ordinary_Z2"}[order]
    expression = _to_scalar(ctx.values[table_name][field], ctx.sp,
                            f"synthetic {field} row {order}")
    # x/xi substitution is the only operation allowed for a synthetic source
    # function.  No source graph or id lookup occurs here.
    return ctx.sp.sympify(expression).subs(ctx.xi_symbol, xi)


def _literal(value: Any, ctx: _EvaluationContext) -> Any:
    _require(isinstance(value, str), "exact_scalar_literal.value must be a string")
    try:
        return ctx.sp.Rational(value)
    except Exception as exc:
        raise FGapEndpointFixtureInputError(f"malformed exact scalar literal: {value!r}") from exc


def _eval_node(node: Any, ctx: _EvaluationContext) -> Any:
    _require(isinstance(node, Mapping), "expression node must be an object")
    operator = node.get("operator")
    if operator is not None:
        _require(isinstance(operator, str), "expression operator must be a string")
        ctx.count_operator(operator)
        arguments = node.get("arguments")
        _require(isinstance(arguments, list), f"{operator}: arguments must be a list")
        if operator == "opaque_source_ref":
            # An opaque handle is legal only as the first operand of
            # source_function_at and is never a scalar expression.
            ctx.opaque_ids_used_as_values = True
            raise FGapEndpointFixtureInputError("opaque_source_ref cannot be evaluated as a numeric value")
        if operator == "source_function_at":
            _require(len(arguments) == 2, "source_function_at: expected ref and xi")
            role = _opaque_ref_from_node(arguments[0], ctx)
            xi = _eval_node(arguments[1], ctx)
            ctx.source_function_calls += 1
            return _row_value(role, ctx, xi)
        if operator == "add":
            _require(len(arguments) >= 1, "add: at least one argument required")
            return ctx.sp.Add(*[_eval_node(item, ctx) for item in arguments])
        if operator == "subtract":
            _require(len(arguments) == 2, "subtract: exactly two arguments required")
            return _eval_node(arguments[0], ctx) - _eval_node(arguments[1], ctx)
        if operator == "multiply":
            _require(len(arguments) >= 1, "multiply: at least one argument required")
            return ctx.sp.Mul(*[_eval_node(item, ctx) for item in arguments])
        if operator == "divide":
            _require(len(arguments) == 2, "divide: exactly two arguments required")
            return _eval_node(arguments[0], ctx) / _eval_node(arguments[1], ctx)
        if operator == "negate":
            _require(len(arguments) == 1, "negate: exactly one argument required")
            return -_eval_node(arguments[0], ctx)
        if operator == "power":
            _require(len(arguments) == 2, "power: exactly two arguments required")
            return ctx.sp.Pow(_eval_node(arguments[0], ctx), _eval_node(arguments[1], ctx))
        if operator == "exp":
            _require(len(arguments) == 1, "exp: exactly one argument required")
            return ctx.sp.exp(_eval_node(arguments[0], ctx))
        raise FGapEndpointFixtureInputError(f"unsupported/malformed operator: {operator}")

    # ``_actual_ref`` in the published candidate is intentionally a direct
    # source-role metadata object for geometry.delta (not an operator wrapper).
    # Admit that one safe shape only as a role lookup; the id is discarded.
    if node.get("status") in {"candidate_source_ref", "saved_opaque_ref"}:
        return _role_value(_raw_opaque_role(node, ctx), ctx)
    kind = node.get("kind")
    _require(isinstance(kind, str), "expression node must have an operator or kind")
    if kind == "exact_scalar_literal":
        return _literal(node.get("value"), ctx)
    if kind == "formal_source_symbol":
        name = node.get("name")
        _require(isinstance(name, str) and name, "formal_source_symbol.name missing")
        if name == "Z":
            return _to_scalar(ctx.values["Z"], ctx.sp, "synthetic Z")
        role = node.get("source_role")
        _require(isinstance(role, str) and role, f"formal source symbol {name}: source_role missing")
        return _role_value(role, ctx)
    if kind == "unresolved_source_role":
        role = node.get("source_role")
        _require(isinstance(role, str) and role, "unresolved_source_role.source_role missing")
        return _role_value(role, ctx)
    raise FGapEndpointFixtureInputError(f"unsupported/malformed expression kind: {kind}")


def _candidate_root(prior_result: Mapping[str, Any]) -> Mapping[str, Any]:
    """Locate the frozen result payload without opening or importing anything."""

    _require(isinstance(prior_result, Mapping), "prior_result must be a mapping")
    current: Any = prior_result
    if isinstance(current.get("result"), Mapping):
        current = current["result"]
    if isinstance(current.get("source_role_candidate"), Mapping):
        current = current["source_role_candidate"]
    _require(isinstance(current, Mapping), "source_role_candidate payload missing")
    # The fixture may consume only the source-only candidate packet.  A packet
    # that has selected values, traversed a graph, or claimed certification is
    # outside this evaluator's authority and is rejected before any DAG node.
    _require(current.get("source_only") is True, "candidate source_only gate is not true")
    for flag in ("source_values_selected", "source_values_computed", "source_graph_traversed",
                 "graph_traversed", "certified", "numeric_use_admissible"):
        if flag in current:
            _require(current.get(flag) is False, f"candidate scope gate is not false: {flag}")
    _require(current.get("source_values_selected") is False,
             "candidate source_values_selected gate is not false")
    _require(current.get("certified") is False, "candidate certified gate is not false")
    _require(isinstance(current.get("endpoints"), Mapping), "candidate endpoints missing")
    _require(isinstance(current.get("direction"), Mapping), "candidate direction missing")
    return current


def _evaluate_emitted_dag(
    candidate: Mapping[str, Any],
    values: Mapping[str, Any],
    refs: Mapping[str, Any],
    sp: Any,
) -> tuple[dict[str, Any], _EvaluationContext]:
    endpoints = candidate["endpoints"]
    ctx = _EvaluationContext(sp, values, refs)
    actual: dict[str, Any] = {"endpoints": {}, "direction": {}}
    for endpoint in ENDPOINT_NAMES:
        packet = endpoints.get(endpoint)
        _require(isinstance(packet, Mapping), f"candidate endpoints.{endpoint}: object missing")
        params = packet.get("parameters")
        fields = packet.get("fields")
        eight = packet.get("eight_fields")
        _require(isinstance(params, Mapping), f"{endpoint}.parameters: object missing")
        _require(isinstance(fields, Mapping), f"{endpoint}.fields: object missing")
        _require(isinstance(eight, Mapping), f"{endpoint}.eight_fields: object missing")
        evaluated_params = {}
        for key in ("xi", "xi_Z", "R", "R_Z", "L", "L_Z"):
            _require(key in params, f"{endpoint}.parameters.{key}: missing")
            evaluated_params[key] = _eval_node(params[key], ctx)
        ctx.current_endpoint = endpoint
        ctx.current_xi = evaluated_params["xi"]
        ctx.endpoint_xis[endpoint] = evaluated_params["xi"]
        evaluated_fields: dict[str, Any] = {}
        for field in FIELD_NAMES:
            field_packet = fields.get(field)
            _require(isinstance(field_packet, Mapping), f"{endpoint}.fields.{field}: object missing")
            evaluated_fields[field] = {}
            for key in ("value", "fixed_coordinate_Z1", "fixed_coordinate_xi_partial", "total_Z1"):
                _require(key in field_packet, f"{endpoint}.fields.{field}.{key}: missing")
                evaluated_fields[field][key] = _eval_node(field_packet[key], ctx)
        evaluated_eight: dict[str, Any] = {}
        for name in ENDPOINT_FIELDS:
            field_packet = eight.get(name)
            _require(isinstance(field_packet, Mapping), f"{endpoint}.eight_fields.{name}: object missing")
            for key in ("candidate", "ordinary_Z1_candidate"):
                _require(key in field_packet, f"{endpoint}.eight_fields.{name}.{key}: missing")
            evaluated_eight[name] = {
                "value": _eval_node(field_packet["candidate"], ctx),
                "ordinary_Z1": _eval_node(field_packet["ordinary_Z1_candidate"], ctx),
            }
        for key in ("K", "K_Z"):
            _require(key in packet, f"{endpoint}.{key}: missing")
        actual["endpoints"][endpoint] = {
            "parameters": evaluated_params,
            "fields": evaluated_fields,
            "eight_fields": evaluated_eight,
            "K": _eval_node(packet["K"], ctx),
            "K_Z": _eval_node(packet["K_Z"], ctx),
        }
    direction = candidate["direction"]
    for key in ("K", "K_Z"):
        _require(key in direction, f"direction.{key}: missing")
    actual["direction"] = {
        "K": _eval_node(direction["K"], ctx),
        "K_Z": _eval_node(direction["K_Z"], ctx),
    }
    return actual, ctx


def _independent_sympy_expectation(values: Mapping[str, Any], sp: Any) -> dict[str, Any]:
    """Derive the expected endpoint algebra independently of the emitted DAG."""

    z = values["z_symbol"] if "z_symbol" in values else values["Z"]
    xi_symbol = values["xi_symbol"]
    scalars = values["scalars"]
    mu = _to_scalar(scalars["historical_F.mu"], sp, "expected mu")
    log_rp = _to_scalar(scalars["historical_F.logRp"], sp, "expected logRp")
    delta = _to_scalar(scalars[F_PRODUCER_POINTER + ".geometry.delta"], sp, "expected delta")
    expected: dict[str, Any] = {"endpoints": {}, "direction": {}}
    for endpoint in ENDPOINT_NAMES:
        xi = sp.Integer(11) if endpoint == "lower" else sp.Integer(13) - 4 * mu
        xi_z = sp.diff(xi, z)
        radius = sp.exp(log_rp + xi / mu)
        radius_z = sp.diff(radius, z)
        L = 1 - delta * z**2
        L_z = sp.diff(L, z)
        fields: dict[str, Any] = {}
        for field in FIELD_NAMES:
            function = _to_scalar(values["field_functions"][field], sp, f"expected {field} function")
            value = function.subs(xi_symbol, xi)
            fixed_z = sp.diff(function, z).subs(xi_symbol, xi)
            partial = sp.diff(function, xi_symbol).subs(xi_symbol, xi)
            total_z = sp.diff(function.subs(xi_symbol, xi), z)
            fields[field] = {
                "value": value,
                "fixed_coordinate_Z1": fixed_z,
                "fixed_coordinate_xi_partial": partial,
                "total_Z1": total_z,
            }
        R = radius
        Q = fields["Q"]["value"]
        U = fields["U"]["value"]
        Q_y = fields["Q_y"]["value"]
        Q_z = fields["Q"]["total_Z1"]
        U_z = fields["U"]["total_Z1"]
        Q_y_z = fields["Q_y"]["total_Z1"]
        values8 = {
            "R_q": R * Q,
            "R_squared_q": R**2 * Q,
            "R_q_squared": R * Q**2,
            "R_squared_q_squared": R**2 * Q**2,
            "R_U_squared": R * U**2,
            "R_squared_U_squared": R**2 * U**2,
            "q_plus_q_y": Q + Q_y,
            "R_q_y": R * Q_y,
        }
        values8_z = {name: sp.diff(value, z) for name, value in values8.items()}
        A = R**2 * Q
        B = R**2 * Q**2
        C = R**2 * U**2
        D = R * Q_y
        K = A / (2 * L) + B / 4 - z**2 * C / L**2 - D
        expected["endpoints"][endpoint] = {
            "parameters": {"xi": xi, "xi_Z": xi_z, "R": R, "R_Z": radius_z, "L": L, "L_Z": L_z},
            "fields": fields,
            "eight_fields": {
                name: {"value": values8[name], "ordinary_Z1": values8_z[name]}
                for name in ENDPOINT_FIELDS
            },
            "K": K,
            "K_Z": sp.diff(K, z),
            "independent_checks": {
                "mu_Z": sp.diff(mu, z),
                "logRp_Z": sp.diff(log_rp, z),
                "delta_Z": sp.diff(delta, z),
                "xi_partial": {field: fields[field]["fixed_coordinate_xi_partial"] for field in FIELD_NAMES},
                "source_function_total_Z": {field: fields[field]["total_Z1"] for field in FIELD_NAMES},
            },
        }
    expected["direction"] = {
        "K": expected["endpoints"]["lower"]["K"] - expected["endpoints"]["upper"]["K"],
        "K_Z": sp.diff(
            expected["endpoints"]["lower"]["K"] - expected["endpoints"]["upper"]["K"], z
        ),
    }
    return expected


def _simplifies_equal(actual: Any, expected: Any, sp: Any) -> bool:
    try:
        return bool(sp.simplify(actual - expected) == 0)
    except Exception:
        return False


def _is_nonzero(value: Any, sp: Any) -> bool:
    try:
        return bool(sp.simplify(value) != 0)
    except Exception:
        return False


def _comparison(actual: Any, expected: Any, sp: Any) -> dict[str, Any]:
    return {
        "matches": _simplifies_equal(actual, expected, sp),
        "actual": str(sp.simplify(actual)),
        "expected": str(sp.simplify(expected)),
    }


def _compare_results(actual: Mapping[str, Any], expected: Mapping[str, Any], values: Mapping[str, Any], sp: Any) -> dict[str, Any]:
    checks: dict[str, Any] = {"endpoints": {}, "direction": {}}
    for endpoint in ENDPOINT_NAMES:
        a_ep, e_ep = actual["endpoints"][endpoint], expected["endpoints"][endpoint]
        ep_checks: dict[str, Any] = {
            "parameters": {
                key: _comparison(a_ep["parameters"][key], e_ep["parameters"][key], sp)
                for key in ("xi", "xi_Z", "R", "R_Z", "L", "L_Z")
            },
            "fields": {},
            "eight_fields": {},
            "K": _comparison(a_ep["K"], e_ep["K"], sp),
            "K_Z": _comparison(a_ep["K_Z"], e_ep["K_Z"], sp),
        }
        for field in FIELD_NAMES:
            ep_checks["fields"][field] = {
                key: _comparison(a_ep["fields"][field][key], e_ep["fields"][field][key], sp)
                for key in ("value", "fixed_coordinate_Z1", "fixed_coordinate_xi_partial", "total_Z1")
            }
        for name in ENDPOINT_FIELDS:
            ep_checks["eight_fields"][name] = {
                "value": _comparison(a_ep["eight_fields"][name]["value"], e_ep["eight_fields"][name]["value"], sp),
                "ordinary_Z1": _comparison(
                    a_ep["eight_fields"][name]["ordinary_Z1"], e_ep["eight_fields"][name]["ordinary_Z1"], sp
                ),
            }
        ep_checks["all_eight_fields_match"] = all(
            ep_checks["eight_fields"][name][part]["matches"]
            for name in ENDPOINT_FIELDS
            for part in ("value", "ordinary_Z1")
        )
        ep_checks["all_fields_match"] = all(
            check["matches"]
            for field in FIELD_NAMES
            for check in ep_checks["fields"][field].values()
        )
        ep_checks["K_and_K_Z_match"] = ep_checks["K"]["matches"] and ep_checks["K_Z"]["matches"]
        checks["endpoints"][endpoint] = ep_checks
    checks["direction"] = {
        "K": _comparison(actual["direction"]["K"], expected["direction"]["K"], sp),
        "K_Z": _comparison(actual["direction"]["K_Z"], expected["direction"]["K_Z"], sp),
    }
    checks["direction"]["matches"] = checks["direction"]["K"]["matches"] and checks["direction"]["K_Z"]["matches"]

    supplied_scalars = values["scalars"]
    scalar_checks = {
        "mu_Z": _comparison(supplied_scalars["historical_F.mu_Z"], sp.diff(supplied_scalars["historical_F.mu"], values["z_symbol"]), sp),
        "logRp_Z": _comparison(supplied_scalars["historical_F.logRp_Z"], sp.diff(supplied_scalars["historical_F.logRp"], values["z_symbol"]), sp),
        "delta_Z": _comparison(
            supplied_scalars["historical_F.delta_Z"],
            sp.diff(supplied_scalars[F_PRODUCER_POINTER + ".geometry.delta"], values["z_symbol"]),
            sp,
        ),
    }
    field_jet_checks: dict[str, Any] = {}
    for field in FIELD_NAMES:
        function = values["field_functions"][field]
        field_jet_checks[field] = {
            "ordinary_Z1_supplied_matches_SymPy": _comparison(
                values["field_ordinary_Z1"][field], sp.diff(function, values["z_symbol"]), sp
            ),
            "ordinary_Z2_supplied_matches_SymPy": _comparison(
                values["field_ordinary_Z2"][field], sp.diff(function, values["z_symbol"], 2), sp
            ),
            "xi_partial_supplied_matches_SymPy": _comparison(
                values["field_xi_partials"][field], sp.diff(function, values["xi_symbol"]), sp
            ),
            "xi_partial_nonzero": _is_nonzero(values["field_xi_partials"][field], sp),
        }
    checks["independent_supplied_derivative_checks"] = {
        "geometry": scalar_checks,
        "fields": field_jet_checks,
    }
    lower_xi = actual["endpoints"]["lower"]["parameters"]["xi"]
    lower_xi_z = actual["endpoints"]["lower"]["parameters"]["xi_Z"]
    upper_xi_z = actual["endpoints"]["upper"]["parameters"]["xi_Z"]
    checks["endpoint_assumptions"] = {
        "lower_xi_is_fixed_11": _simplifies_equal(lower_xi, sp.Integer(11), sp),
        "lower_xi_Z_is_zero": _simplifies_equal(lower_xi_z, sp.Integer(0), sp),
        "upper_endpoint_moves": _is_nonzero(upper_xi_z, sp),
        "mu_Z_nonzero": _is_nonzero(supplied_scalars["historical_F.mu_Z"], sp),
        "logRp_Z_nonzero": _is_nonzero(supplied_scalars["historical_F.logRp_Z"], sp),
        "delta_Z_nonzero": _is_nonzero(supplied_scalars["historical_F.delta_Z"], sp),
    }
    return checks


def _negative_tests(values: Mapping[str, Any], refs: Mapping[str, Any], sp: Any) -> dict[str, Any]:
    """Exercise safety rejects without touching a producer or source graph."""

    def run(label: str, fn: Any) -> dict[str, Any]:
        try:
            fn()
        except FGapEndpointFixtureInputError as exc:
            return {"rejected": True, "detail": str(exc)}
        except Exception as exc:  # pragma: no cover - defensive negative path
            return {"rejected": True, "detail": f"unexpected rejection: {exc}"}
        return {"rejected": False, "detail": f"negative case accepted: {label}"}

    def context() -> _EvaluationContext:
        return _EvaluationContext(sp, values, refs)

    opaque = {
        "operator": "opaque_source_ref",
        "arguments": [{
            "status": "candidate_source_ref",
            "source_role": f"{F_LEADING_POINTER}.Q[0]",
            "source_node_id": 991,
            "node_id_is_not_a_numeric_value": True,
            "source_values_selected": False,
        }],
        "source_values_selected": False,
    }
    malformed_opaque = copy.deepcopy(opaque)
    malformed_opaque["arguments"][0]["node_id_is_not_a_numeric_value"] = False

    def unknown_operator() -> None:
        _eval_node({"operator": "sqrt", "arguments": [{"kind": "exact_scalar_literal", "value": "1"}]}, context())

    def opaque_as_number() -> None:
        _eval_node({"operator": "add", "arguments": [opaque, {"kind": "exact_scalar_literal", "value": "1"}]}, context())

    def malformed_ref() -> None:
        _eval_node({"operator": "source_function_at", "arguments": [malformed_opaque, {"kind": "exact_scalar_literal", "value": "11"}]}, context())

    def missing_role() -> None:
        bad = copy.deepcopy(opaque)
        bad["arguments"][0]["source_role"] = "unmapped.opaque.role"
        _eval_node({"operator": "source_function_at", "arguments": [bad, {"kind": "exact_scalar_literal", "value": "11"}]}, context())

    def row_semantic_swap() -> None:
        bad_refs = copy.deepcopy(dict(refs))
        role = f"{F_LEADING_POINTER}.Q[0]"
        bad_refs[role]["ordinary_Z_order"] = 1
        _eval_node({"operator": "source_function_at", "arguments": [opaque, {"kind": "exact_scalar_literal", "value": "11"}]}, _EvaluationContext(sp, values, bad_refs))

    def source_id_value_alias() -> None:
        bad_refs = copy.deepcopy(dict(refs))
        bad_refs["bad.scalar"] = {
            "kind": "scalar",
            "scalar_key": "source_node_id",
            "source_node_id": 991,
            "node_id_is_not_a_numeric_value": True,
            "semantic": "synthetic_scalar",
        }
        _eval_node({"kind": "unresolved_source_role", "source_role": "bad.scalar"}, _EvaluationContext(sp, values, bad_refs))

    cases = {
        "unknown_operator_rejected": unknown_operator,
        "opaque_ref_cannot_be_numeric": opaque_as_number,
        "malformed_opaque_ref_rejected": malformed_ref,
        "unmapped_source_role_rejected": missing_role,
        "row_semantic_swap_rejected": row_semantic_swap,
        "source_node_id_value_alias_rejected": source_id_value_alias,
    }
    result = {name: run(name, fn) for name, fn in cases.items()}
    result["all_negative_cases_rejected"] = all(item["rejected"] for item in result.values())
    return result


def _gate_summary(candidate: Mapping[str, Any]) -> dict[str, Any]:
    retained = candidate.get("retained_legacy")
    retained_legacy_observed = isinstance(retained, Mapping)
    return {
        "F_radius_owner_bound": False,
        "F_endpoint_owner_bound": False,
        "capture_free_endpoint_substitution": False,
        "fixed_coordinate_xi_partials_bound": False,
        "ordinary_Z_total_transport_proved": False,
        "delta_Z_source_bound": False,
        "eight_field_source_identity_proved": False,
        "K_source_identity_proved": False,
        "K_Z_source_identity_proved": False,
        "endpoint_direction_identity_proved": False,
        "numeric_use_admissible": False,
        "source_values_selected": False,
        "old_errors_removed": False,
        "error_transport_certified": False,
        "historical_error_state_mutated_by_fixture": False,
        "retained_legacy_metadata_observed": retained_legacy_observed,
        "retained_legacy_metadata_forwarded_without_reinterpretation": retained_legacy_observed,
    }


def _json_safe(value: Any) -> Any:
    """Convert expressions and nested metadata to JSON-safe values."""

    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def verify_emitted_candidate(
    prior_result: Mapping[str, Any] | None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate a caller-supplied frozen candidate and compare it to SymPy.

    Parameters are intentionally data-only.  ``prior_result`` is the parsed
    authenticated packet from the one-shot runner; this function never reads a
    path.  If custom synthetic inputs are supplied, they must contain the
    fields documented by :func:`build_synthetic_inputs`, and the role map must
    contain semantic entries only, never values selected by opaque ids.
    """

    base: dict[str, Any] = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "synthetic_only": True,
        "source_only": True,
        # Scope flags are deliberately top-level so a runner cannot confuse
        # this algebra check with actual source evaluation or certification.
        "source_values_computed": False,
        "source_graph_traversed": False,
        "graph_traversed": False,
        "source_values_selected": False,
        "source_function_values_selected": False,
        "numeric_use_admissible": False,
        "certified": False,
        "actual_source_function_called": False,
        "opaque_source_node_ids_used_as_values": False,
        "admissible_for_numeric_use": False,
        "numeric_use_admissible": False,
        "old_errors_removed": False,
        "error_transport_certified": False,
        "errors": [],
    }
    if not isinstance(prior_result, Mapping):
        base.update({"status": "rejected_missing_candidate", "candidate_dag_evaluated": False})
        return _json_safe(base)
    try:
        candidate = _candidate_root(prior_result)
        values, refs, sp = _normalise_inputs(synthetic_values, source_role_ref_map)
        actual, context = _evaluate_emitted_dag(candidate, values, refs, sp)
        expected = _independent_sympy_expectation(values, sp)
        checks = _compare_results(actual, expected, values, sp)
        negatives = _negative_tests(values, refs, sp)
        endpoint_checks = checks["endpoints"]
        all_formula_checks = all(
            endpoint_checks[endpoint]["all_eight_fields_match"]
            and endpoint_checks[endpoint]["all_fields_match"]
            and endpoint_checks[endpoint]["K_and_K_Z_match"]
            for endpoint in ENDPOINT_NAMES
        )
        all_assumptions = all(checks["endpoint_assumptions"].values()) and all(
            item["ordinary_Z1_supplied_matches_SymPy"]["matches"]
            and item["ordinary_Z2_supplied_matches_SymPy"]["matches"]
            and item["xi_partial_supplied_matches_SymPy"]["matches"]
            and item["xi_partial_nonzero"]
            for item in checks["independent_supplied_derivative_checks"]["fields"].values()
        ) and all(item["matches"] for item in checks["independent_supplied_derivative_checks"]["geometry"].values())
        all_checks = all_formula_checks and checks["direction"]["matches"] and all_assumptions
        gate_summary = _gate_summary(candidate)
        base.update({
            "status": "synthetic_comparison_passed" if all_checks and negatives["all_negative_cases_rejected"] else "synthetic_comparison_failed",
            "candidate_dag_evaluated": True,
            "candidate_source_values_selected": candidate.get("source_values_selected") is True,
            "actual_source_function_called": False,
            "opaque_source_node_ids_used_as_values": context.opaque_ids_used_as_values,
            "source_role_ref_map_consumed": True,
            "used_source_roles": sorted(context.used_roles),
            "source_function_at_calls": context.source_function_calls,
            "opaque_source_refs_seen": context.opaque_source_refs_seen,
            "operator_counts": dict(sorted(context.used_operator_counts.items())),
            "field_source_row_semantics": {
                "row_0": "caller synthetic field value function, substituted at xi; not a numeric node id",
                "row_1": "caller synthetic fixed-coordinate ordinary Z1 derivative; not total endpoint derivative",
                "row_2": "caller synthetic fixed-coordinate ordinary Z2 derivative; retained but not consumed by v2 emitted endpoint roots",
                "xi_partial": "caller synthetic fixed-coordinate xi partial; consumed only through total_Z1 chain rule",
                "lower_endpoint": "candidate xi=11 exact literal and xi_Z=0",
                "upper_endpoint": "candidate xi=13-4*mu_F and xi_Z=-4*mu_F_Z",
            },
            "checks": checks,
            "negative_tests": negatives,
            "gate_preservation": gate_summary,
            "actual": actual,
            "independent_sympy_expectation": expected,
            "all_formula_checks_passed": all_formula_checks,
            "all_assumption_checks_passed": all_assumptions,
            "all_checks_passed": all_checks,
            "checks_passed": all_checks and negatives["all_negative_cases_rejected"],
            "limitations": [
                "Synthetic data are caller supplied and are not actual F source values.",
                "Opaque source_node_id metadata is validated and discarded, never evaluated.",
                "A successful comparison proves only the toy algebra and supplied jet semantics.",
                "F owner identity, capture-free substitution, ordinary-Z transport, same-source seams, numeric admissibility, and certified error transport remain false.",
                "Historical error metadata is not rewritten or removed by this fixture.",
            ],
        })
    except Exception as exc:
        base.update({
            "status": "synthetic_comparison_rejected",
            "candidate_dag_evaluated": False,
            "errors": [str(exc)],
            "checks_passed": False,
            "gate_preservation": _gate_summary(prior_result if isinstance(prior_result, Mapping) else {}),
            "negative_tests": {"all_negative_cases_rejected": False, "not_run": True},
        })
    return _json_safe(base)


def build_and_test(
    prior_result: Mapping[str, Any] | None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """One-shot runner alias retained for explicit build/test call sites."""

    return verify_emitted_candidate(prior_result, synthetic_values, source_role_ref_map)


def test_emitted_candidate(
    prior_result: Mapping[str, Any] | None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Test-oriented alias for :func:`verify_emitted_candidate`."""

    return verify_emitted_candidate(prior_result, synthetic_values, source_role_ref_map)


def verify_emitted_dag_synthetic_fixture(
    emitted_dag: Mapping[str, Any] | None = None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Runner-discovery entry point; ``emitted_dag`` is already authenticated."""

    return verify_emitted_candidate(emitted_dag, synthetic_values, source_role_ref_map)


def verify_emitted_dag_fixture(
    emitted_dag: Mapping[str, Any] | None = None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Compatibility alias for one-shot emitted-DAG runners."""

    return verify_emitted_candidate(emitted_dag, synthetic_values, source_role_ref_map)


def verify_synthetic_fixture(
    emitted_dag: Mapping[str, Any] | None = None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Compatibility alias used by fixture-oriented runners."""

    return verify_emitted_candidate(emitted_dag, synthetic_values, source_role_ref_map)


def synthetic_fixture(
    emitted_dag: Mapping[str, Any] | None = None,
    synthetic_values: Mapping[str, Any] | None = None,
    source_role_ref_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Short compatibility alias; no producer or result file is read."""

    return verify_emitted_candidate(emitted_dag, synthetic_values, source_role_ref_map)


__all__ = [
    "COMPONENT",
    "ENDPOINT_FIELDS",
    "FGapEndpointFixtureInputError",
    "F_PRODUCER_POINTER",
    "F_LEADING_POINTER",
    "SCHEMA",
    "build_and_test",
    "build_synthetic_inputs",
    "synthetic_fixture",
    "test_emitted_candidate",
    "verify_emitted_dag_fixture",
    "verify_emitted_dag_synthetic_fixture",
    "verify_synthetic_fixture",
    "verify_emitted_candidate",
]
