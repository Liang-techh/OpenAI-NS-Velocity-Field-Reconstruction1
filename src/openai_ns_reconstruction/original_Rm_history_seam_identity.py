"""Bounded source consumer for the post-restore ``Rm`` history seam.

This module consumes saved graph ports.  It does not import a transport
operator, replay an accepted graph, decode an archive, evaluate a function,
or select a numerical value.  The graph checks below are deliberately local:
they authenticate the exact post-restore Duhamel source and the exact patch
incoming substitutions before appending the first native-y source rows.

The returned rows are source identities only.  In particular, a returned
``density(1)`` substitution is a guarded graph definition, not a value of a
density.  The only endpoint source used by the new derivative row is the
already authenticated patch incoming handle, whose expression is the current
post-restore ``Rm`` history at ``p84528 = 1``.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


FIELDS = ("m", "h", "k", "e", "p")
ORDERS = (0, 1)
RATES = {
    "m": Fraction(1, 1),
    "h": Fraction(3, 2),
    "k": Fraction(3, 2),
    "e": Fraction(1, 1),
    "p": Fraction(0, 1),
}

POST_COORDINATE = 84528
PATCH_COORDINATE = 185676
ENDPOINT_ONE = 1
LOOP_DENSITY_ARCHIVE_SHA256 = (
    "29c4fec78d3c622b36bb10867792b3d4df1b62f332cf3bc93cd5ed684706fe29"
)
PATCH_INCOMING_ROOTS = {
    "m": (121756, 121757),
    "h": (121762, 121763),
    "k": (121768, 121769),
    "e": (121774, 121775),
    "p": (121780, 121781),
}

_GUARD_KEYS = (
    "Z_independent_substitution",
    "native_substitution_guard",
    "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard",
    "guarded",
)
_FIELD_ALIASES = {field: (field, field.upper()) for field in FIELDS}


def canonical_digest(value: Any) -> str:
    """Hash JSON-shaped source data without normalizing graph algebra."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def static_source_proof() -> dict[str, Any]:
    """Return static source facts for the runner's AST/LF/EOF audit.

    This reads only this source file.  It does not inspect or execute a saved
    graph.  The runner may repeat the checks independently.
    """

    raw = Path(__file__).read_bytes()
    text = raw.decode("utf-8")
    ast.parse(text, filename=str(Path(__file__)))
    return {
        "ast_parse": True,
        "utf8": True,
        "lf_only": b"\r" not in raw,
        "single_eof": raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
        "source_bytes_sha256": hashlib.sha256(raw).hexdigest(),
        "source_only": True,
        "no_provider_execution": True,
        "no_archive_decode": True,
    }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sequence(value: Any, label: str) -> list[Any]:
    _require(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
        f"{label} must be a sequence",
    )
    return list(value)


def _handle(value: Any, label: str, count: int) -> int:
    if isinstance(value, Mapping):
        value = value.get("node", value.get("handle", value.get("root")))
    _require(
        isinstance(value, int) and not isinstance(value, bool),
        f"{label} must be an integer graph handle",
    )
    _require(0 <= value < count, f"{label} is outside the source graph")
    return value


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        text = value.strip()
        try:
            if "/" in text:
                parts = text.split("/")
                _require(len(parts) == 2, f"{label} is not an exact rational")
                return Fraction(int(parts[0]), int(parts[1]))
            return Fraction(int(text), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{label} is not an exact rational") from exc
    raise ValueError(f"{label} is not an exact rational")


def _rate_text(rate: Fraction) -> str:
    return str(rate.numerator) if rate.denominator == 1 else (
        f"{rate.numerator}/{rate.denominator}"
    )


def _op(nodes: Sequence[Mapping[str, Any]], index: int) -> str:
    row = nodes[index]
    _require(isinstance(row, Mapping), f"node {index} is not a mapping")
    operation = row.get("operation", row.get("op"))
    _require(isinstance(operation, str) and operation, f"node {index} has no operation")
    return operation


def _node(nodes: Sequence[Mapping[str, Any]], index: int, label: str = "node") -> Mapping[str, Any]:
    _handle(index, label, len(nodes))
    row = nodes[index]
    _require(isinstance(row, Mapping), f"{label} {index} is not a mapping")
    return row


def _args(nodes: Sequence[Mapping[str, Any]], index: int) -> list[int]:
    row = _node(nodes, index)
    raw = row.get("arguments", row.get("operands"))
    values = _sequence(raw, f"node {index}.arguments")
    return [
        _handle(value, f"node {index}.arguments[{offset}]", len(nodes))
        for offset, value in enumerate(values)
    ]


def _exact(nodes: Sequence[Mapping[str, Any]], index: int) -> Fraction | None:
    if _op(nodes, index) != "exact_rational":
        return None
    row = _node(nodes, index)
    try:
        numerator = _fraction(row.get("numerator", row.get("value")), "rational numerator")
        denominator = _fraction(row.get("denominator", 1), "rational denominator")
        _require(denominator != 0, "rational denominator is zero")
        return numerator / denominator
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _guarded(row: Mapping[str, Any]) -> bool:
    return any(row.get(key) is True for key in _GUARD_KEYS)


def _first(owner: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in owner:
            return owner[key]
    return None


def _field_raw(container: Mapping[str, Any], field: str) -> Any:
    for key in _FIELD_ALIASES[field]:
        if key in container:
            return container[key]
    roots = container.get("roots")
    if isinstance(roots, Mapping):
        for key in _FIELD_ALIASES[field]:
            if key in roots:
                return roots[key]
    return None


def _rows(
    raw: Any,
    label: str,
    count: int,
    orders: Sequence[int] = ORDERS,
) -> dict[int, int]:
    result: dict[int, int] = {}
    if isinstance(raw, Mapping):
        for order in orders:
            value = raw.get(order, raw.get(str(order)))
            _require(value is not None, f"{label} lacks ordinary order {order}")
            result[order] = _handle(value, f"{label}[{order}]", count)
        return result
    values = _sequence(raw, label)
    for order in orders:
        _require(order < len(values), f"{label} lacks ordinary order {order}")
        result[order] = _handle(values[order], f"{label}[{order}]", count)
    return result


def _field_rows(
    container: Mapping[str, Any],
    field: str,
    label: str,
    count: int,
) -> dict[int, int]:
    raw = _field_raw(container, field)
    _require(raw is not None, f"{label} is missing field {field}")
    return _rows(raw, f"{label}.{field}", count)


def _flatten_product(nodes: Sequence[Mapping[str, Any]], root: int) -> list[int]:
    if _op(nodes, root) not in {"product", "mul"}:
        return [root]
    arguments = _args(nodes, root)
    _require(arguments, f"product node {root} has no factors")
    result: list[int] = []
    for argument in arguments:
        result.extend(_flatten_product(nodes, argument))
    return result


def _unary_argument(nodes: Sequence[Mapping[str, Any]], index: int, label: str) -> int:
    row = _node(nodes, index, label)
    value = row.get("argument", row.get("operand", row.get("value")))
    return _handle(value, f"{label}.argument", len(nodes))


def _affine(
    nodes: Sequence[Mapping[str, Any]],
    root: int,
    leaves: set[int],
    stack: tuple[int, ...] = (),
) -> tuple[dict[int, Fraction], Fraction]:
    """Parse only the exact affine source algebra used by the kernel."""

    _require(root not in stack, f"cyclic kernel source at node {root}")
    operation = _op(nodes, root)
    if operation == "exact_rational":
        value = _exact(nodes, root)
        _require(value is not None, f"node {root} is not an exact rational")
        return {}, value
    if operation == "bound_variable":
        _require(root in leaves, f"kernel uses an out-of-scope variable node {root}")
        return {root: Fraction(1, 1)}, Fraction(0, 1)
    if operation in {"negative", "neg"}:
        return _scale(
            Fraction(-1, 1),
            _affine(nodes, _unary_argument(nodes, root, "negative source"), leaves, stack + (root,)),
        )
    if operation in {"sum", "add"}:
        values = _args(nodes, root)
        _require(values, f"sum node {root} has no terms")
        result: tuple[dict[int, Fraction], Fraction] = ({}, Fraction(0, 1))
        for value in values:
            result = _add(result, _affine(nodes, value, leaves, stack + (root,)))
        return result
    if operation in {"product", "mul"}:
        values = _args(nodes, root)
        _require(values, f"product node {root} has no factors")
        scalar = Fraction(1, 1)
        affine: tuple[dict[int, Fraction], Fraction] | None = None
        for value in values:
            current = _affine(nodes, value, leaves, stack + (root,))
            if current[0]:
                _require(affine is None, f"kernel node {root} is nonlinear")
                affine = current
            else:
                scalar *= current[1]
        if affine is None:
            return {}, scalar
        return _scale(scalar, affine)
    raise ValueError(f"node {root} is outside the exact affine kernel source")


def _add(
    left: tuple[dict[int, Fraction], Fraction],
    right: tuple[dict[int, Fraction], Fraction],
) -> tuple[dict[int, Fraction], Fraction]:
    coefficients = dict(left[0])
    for key, value in right[0].items():
        coefficients[key] = coefficients.get(key, Fraction(0, 1)) + value
        if coefficients[key] == 0:
            del coefficients[key]
    return coefficients, left[1] + right[1]


def _scale(
    scalar: Fraction,
    value: tuple[dict[int, Fraction], Fraction],
) -> tuple[dict[int, Fraction], Fraction]:
    return (
        {key: scalar * item for key, item in value[0].items() if scalar * item},
        scalar * value[1],
    )


def _unwrap_left(value: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any] | None]:
    current: Mapping[str, Any] = value
    if isinstance(current.get("prescribed_loop_C5_bindings"), Mapping):
        current = current["prescribed_loop_C5_bindings"]
    if isinstance(current.get("charts"), Mapping) and isinstance(current["charts"].get("postrestore"), Mapping):
        chart = current["charts"]["postrestore"]
    elif isinstance(current.get("postrestore"), Mapping):
        chart = current["postrestore"]
    else:
        chart = current
    outlet = chart.get("outlet") if isinstance(chart.get("outlet"), Mapping) else None
    return chart, outlet


def _unwrap_right(value: Mapping[str, Any]) -> Mapping[str, Any]:
    current: Mapping[str, Any] = value
    if isinstance(current.get("Rm_patch_finite_N_C5_bindings"), Mapping):
        current = current["Rm_patch_finite_N_C5_bindings"]
    if isinstance(current.get("chart"), Mapping):
        current = current["chart"]
    return current


def _chart_container(chart: Mapping[str, Any], *keys: str) -> Mapping[str, Any]:
    for key in keys:
        candidate = chart.get(key)
        if isinstance(candidate, Mapping):
            return candidate
    raise ValueError(f"chart is missing {'/'.join(keys)}")


def _density_container(chart: Mapping[str, Any]) -> Mapping[str, Any]:
    finite = chart.get("finite_N")
    owners = (finite, chart) if isinstance(finite, Mapping) else (chart,)
    for owner in owners:
        for key in (
            "five_signed_density_increments",
            "signed_density_increments",
            "density",
            "source_density",
        ):
            candidate = owner.get(key)
            if isinstance(candidate, Mapping):
                return candidate
    raise ValueError("postrestore chart is missing five signed density increments")


def _coordinate(nodes: Sequence[Mapping[str, Any]], chart: Mapping[str, Any], label: str) -> int:
    value = _first(chart, "coordinate", "native_coordinate", "x")
    return _handle(value, f"{label} coordinate", len(nodes))


def _domain(chart: Mapping[str, Any], expected: Sequence[int], label: str) -> list[int]:
    raw = _sequence(chart.get("domain"), f"{label}.domain")
    _require(raw == list(expected), f"{label}.domain is not {list(expected)!r}")
    return [int(item) for item in raw]


def _physical_info(
    nodes: Sequence[Mapping[str, Any]],
    chart: Mapping[str, Any],
    coordinate: int,
    label: str,
    *,
    postrestore: bool,
) -> dict[str, Any]:
    value = _first(chart, "physical_measure", "physical_Jacobian", "jacobian")
    physical = _handle(value, f"{label} physical measure", len(nodes))
    if postrestore:
        _require(_exact(nodes, physical) == 1, f"{label} physical measure is not exact one")
        return {
            "node": physical,
            "endpoint_value": "1",
            "endpoint_J_exact_one": True,
            "factor_count": "metadata_only_once",
        }
    row = _node(nodes, physical, f"{label} physical measure")
    _require(_op(nodes, physical) in {"positive_quotient", "quotient"},
             f"{label} physical measure is not a positive quotient")
    numerator = _handle(row.get("numerator"), f"{label} physical numerator", len(nodes))
    denominator = _handle(row.get("denominator"), f"{label} physical denominator", len(nodes))
    _require(_exact(nodes, numerator) == 1, f"{label} physical numerator is not exact one")
    _require(denominator == coordinate, f"{label} physical denominator changed")
    return {
        "node": physical,
        "formula": "1/x",
        "numerator": numerator,
        "denominator": denominator,
        "lower_endpoint_node": ENDPOINT_ONE,
        "endpoint_J_exact_one": True,
        "factor_count": "metadata_only_once",
    }


def _validate_postrestore(
    nodes: Sequence[Mapping[str, Any]],
    chart: Mapping[str, Any],
) -> tuple[dict[tuple[str, int], dict[str, int]], dict[str, Any]]:
    coordinate = _coordinate(nodes, chart, "postrestore")
    _require(coordinate == POST_COORDINATE, "postrestore coordinate is not p84528")
    _domain(chart, (0, 1), "postrestore")
    physical = _physical_info(nodes, chart, coordinate, "postrestore", postrestore=True)
    density = _density_container(chart)
    incoming = _chart_container(chart, "incoming")
    history = _chart_container(chart, "history")
    contributions = _chart_container(chart, "contributions", "contribution")
    memory = _chart_container(chart, "memory")
    finite = chart.get("finite_N")
    _require(isinstance(finite, Mapping), "postrestore chart is missing finite_N")
    densities = _density_container(chart)

    ports: dict[tuple[str, int], dict[str, int]] = {}
    evidence: list[dict[str, Any]] = []
    for field in FIELDS:
        drows = _field_rows(densities, field, "postrestore finite_N density", len(nodes))
        irows = _field_rows(incoming, field, "postrestore incoming", len(nodes))
        hrows = _field_rows(history, field, "postrestore history", len(nodes))
        crows = _field_rows(contributions, field, "postrestore contribution", len(nodes))
        raw_memory = _field_raw(memory, field)
        _require(raw_memory is not None, f"postrestore memory is missing field {field}")
        mroot = _handle(raw_memory, f"postrestore memory.{field}", len(nodes))
        rate = RATES[field]
        memory_info = _validate_memory(nodes, mroot, coordinate, rate, field, "postrestore")
        for order in ORDERS:
            density_root = drows[order]
            incoming_root = irows[order]
            history_root = hrows[order]
            contribution_root = crows[order]
            _require(density_root != incoming_root, f"{field} Z{order} density aliases incoming")
            _require(_exact(nodes, incoming_root) != Fraction(0, 1),
                     f"{field} Z{order} incoming is an exact zero")
            contribution_info = _validate_contribution(
                nodes,
                contribution_root,
                density_root,
                coordinate,
                rate,
                field,
                order,
            )
            history_info = _validate_history(
                nodes,
                history_root,
                contribution_root,
                incoming_root,
                mroot,
                rate,
                field,
                order,
            )
            key = (field, order)
            ports[key] = {
                "density": density_root,
                "history": history_root,
                "incoming": incoming_root,
                "memory": mroot,
                "contribution": contribution_root,
            }
            evidence.append({
                "field": field,
                "ordinary_Z_order": order,
                "density_root": density_root,
                "history_root": history_root,
                "incoming_root": incoming_root,
                "memory_root": mroot,
                "contribution_root": contribution_root,
                "rate": _rate_text(rate),
                "history_body": history_info,
                "contribution_source": contribution_info,
                "memory_source": memory_info,
            })
    return ports, {
        "coordinate": coordinate,
        "domain": [0, 1],
        "physical": physical,
        "ports": {field: {str(order): dict(ports[(field, order)]) for order in ORDERS} for field in FIELDS},
        "rows": evidence,
        "full_topology": True,
    }


def _validate_memory(
    nodes: Sequence[Mapping[str, Any]],
    root: int,
    coordinate: int,
    rate: Fraction,
    field: str,
    chart_name: str,
) -> dict[str, Any]:
    if rate == 0:
        _require(_exact(nodes, root) == 1, f"{chart_name} {field} memory is not exact one")
        return {"node": root, "kind": "exact_one", "endpoint_unit": True}
    _require(_op(nodes, root) in {"analytic_unary", "exp"},
             f"{chart_name} {field} memory is not exponential")
    row = _node(nodes, root, "memory")
    name = row.get("name", row.get("function", "exp"))
    _require(name == "exp", f"{chart_name} {field} memory is not exp")
    argument = _unary_argument(nodes, root, "memory")
    affine = _affine(nodes, argument, {coordinate})
    expected = {coordinate: -rate}, Fraction(0, 1)
    _require(affine == expected, f"{chart_name} {field} memory exponential changed")
    return {
        "node": root,
        "kind": "exp_affine_width",
        "argument": argument,
        "rate": _rate_text(rate),
        "affine": {str(coordinate): _rate_text(-rate)},
        "endpoint_unit": False,
        "endpoint_value_formula": "exp(-" + _rate_text(rate) + ") at p=1",
        "outlet_memory_never_reset": True,
    }


def _validate_contribution(
    nodes: Sequence[Mapping[str, Any]],
    root: int,
    density_root: int,
    coordinate: int,
    rate: Fraction,
    field: str,
    order: int,
) -> dict[str, Any]:
    row = _node(nodes, root, "postrestore contribution")
    _require(_op(nodes, root) == "definite_integral",
             f"postrestore {field} Z{order} contribution is not definite_integral")
    _require(row.get("lower") == 0, f"postrestore {field} Z{order} lower endpoint changed")
    _require(row.get("upper") == coordinate,
             f"postrestore {field} Z{order} upper coordinate changed")
    _require(_fraction(row.get("original_rate"), f"postrestore {field} rate") == rate,
             f"postrestore {field} Z{order} rate changed")
    required = {
        "actual_finite_N_signed_correction": True,
        "endpoint_independent_of_Z": True,
        "exact_function_integral": True,
        "numerical_cap_not_defining_integrand": True,
        "numerical_value_not_installed": True,
        "original_chart": "postrestore",
        "physical_Jacobian_applied_once": True,
        "predecessor_memory_never_reset": True,
        "true_core_inlet_complete_history": True,
        "zero_correction_inlet_only_at_original_rminus": True,
    }
    for key, expected in required.items():
        _require(row.get(key) == expected,
                 f"postrestore {field} Z{order} contribution metadata {key} changed")
    _require(row.get("ordinary_slow_Z_derivative_order") == order,
             f"postrestore {field} Z{order} derivative order changed")
    _require(row.get("variable") == "true_loop_C5_transport_postrestore",
             f"postrestore {field} Z{order} transport variable changed")
    integrand = _handle(row.get("integrand"), "postrestore integrand", len(nodes))
    factors = _flatten_product(nodes, integrand)
    density_substitutions: list[int] = []
    for factor in factors:
        if _op(nodes, factor) != "function_substitution":
            continue
        candidate = _node(nodes, factor, "postrestore density substitution")
        if (
            candidate.get("expression") == density_root
            and candidate.get("variable") == coordinate
            and isinstance(candidate.get("value"), int)
        ):
            variable = _handle(candidate.get("value"), "postrestore integration variable", len(nodes))
            variable_row = _node(nodes, variable, "postrestore integration variable")
            _require(_op(nodes, variable) == "bound_variable",
                     "postrestore integration variable is not bound_variable")
            _require(variable_row.get("name") == "true_loop_C5_transport_postrestore",
                     "postrestore integration variable name changed")
            _require(_guarded(candidate),
                     f"postrestore {field} Z{order} density substitution is unguarded")
            density_substitutions.append(factor)
    _require(len(density_substitutions) == 1,
             f"postrestore {field} Z{order} density source count changed")
    non_density = [factor for factor in factors if factor not in density_substitutions]
    if rate == 0:
        _require(len(non_density) == 0,
                 f"postrestore {field} Z{order} rate-zero source has extra factors")
        kernel_info = {"kind": "exact_one", "node": None}
    else:
        _require(len(non_density) == 1,
                 f"postrestore {field} Z{order} source has extra factors")
        kernel = non_density[0]
        _require(_op(nodes, kernel) in {"analytic_unary", "exp"},
                 f"postrestore {field} Z{order} kernel is not exp")
        kernel_row = _node(nodes, kernel, "postrestore kernel")
        _require(kernel_row.get("name", kernel_row.get("function", "exp")) == "exp",
                 f"postrestore {field} Z{order} kernel name changed")
        argument = _unary_argument(nodes, kernel, "postrestore kernel")
        variable = _handle(_node(nodes, density_substitutions[0]).get("value"),
                           "postrestore integration variable", len(nodes))
        affine = _affine(nodes, argument, {coordinate, variable})
        expected = {coordinate: -rate, variable: rate}, Fraction(0, 1)
        _require(affine == expected,
                 f"postrestore {field} Z{order} kernel is not exact exp[-rate*(p-t)]")
        kernel_info = {
            "kind": "exact_exponential",
            "node": kernel,
            "argument": argument,
            "integration_variable": variable,
            "affine": {
                str(coordinate): _rate_text(-rate),
                str(variable): _rate_text(rate),
            },
        }
    return {
        "node": root,
        "integrand": integrand,
        "density_substitution": density_substitutions[0],
        "density_root": density_root,
        "lower": 0,
        "upper": coordinate,
        "rate": _rate_text(rate),
        "kernel": kernel_info,
        "physical_Jacobian_applied_once": True,
        "factor_topology_exact": True,
    }


def _validate_history(
    nodes: Sequence[Mapping[str, Any]],
    root: int,
    contribution: int,
    incoming: int,
    memory: int,
    rate: Fraction,
    field: str,
    order: int,
) -> dict[str, Any]:
    _require(_op(nodes, root) in {"sum", "add"},
             f"postrestore {field} Z{order} history is not a sum")
    arguments = _args(nodes, root)
    _require(arguments.count(contribution) == 1,
             f"postrestore {field} Z{order} history contribution branch changed")
    others = [value for value in arguments if value != contribution]
    _require(len(others) == 1,
             f"postrestore {field} Z{order} history has extra branches")
    incoming_branch = others[0]
    if rate == 0:
        _require(incoming_branch == incoming,
                 f"postrestore {field} Z{order} rate-zero incoming branch changed")
        return {
            "node": root,
            "arguments": arguments,
            "incoming_branch": incoming_branch,
            "incoming_branch_factors": [incoming],
            "incoming_branch_kind": "incoming",
            "incoming_retained": True,
            "memory_retained": False,
        }
    factors = _flatten_product(nodes, incoming_branch)
    _require(len(factors) == 2 and sorted(factors) == sorted((incoming, memory)),
             f"postrestore {field} Z{order} memory*incoming branch changed")
    return {
        "node": root,
        "arguments": arguments,
        "incoming_branch": incoming_branch,
        "incoming_branch_factors": factors,
        "incoming_branch_kind": "memory_times_incoming",
        "incoming_retained": True,
        "memory_retained": True,
    }


def _unwrap_patch_receipt(value: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in (
        "actual_Rm_history_native_source",
        "accepted_patch_history",
        "patch_history",
        "source",
    ):
        candidate = value.get(key)
        if isinstance(candidate, Mapping) and candidate.get("schema") == "original_Rm_history_native_source_v1":
            return candidate
    if value.get("schema") == "original_Rm_history_native_source_v1":
        return value
    raise ValueError("accepted_patch_history lacks original_Rm_history_native_source_v1")


def _packet_rows(packet: Mapping[str, Any]) -> dict[tuple[str, int], Mapping[str, Any]]:
    raw = packet.get("rows")
    values = _sequence(raw, "accepted patch history rows")
    result: dict[tuple[str, int], Mapping[str, Any]] = {}
    for value in values:
        _require(isinstance(value, Mapping), "accepted patch history row is not a mapping")
        field = value.get("field")
        order = value.get("ordinary_Z_order", value.get("order"))
        _require(field in FIELDS and order in ORDERS,
                 "accepted patch history row has unsupported field/order")
        key = (field, order)
        _require(key not in result, "accepted patch history has duplicate field/order")
        result[key] = value
    _require(set(result) == {(field, order) for field in FIELDS for order in ORDERS},
             "accepted patch history does not contain exactly ten ordinary rows")
    return result


def _packet_port(
    packet: Mapping[str, Any],
    field: str,
    order: int,
    label: str,
    count: int,
) -> Mapping[str, Any]:
    ports = packet.get("source_ports")
    _require(isinstance(ports, Mapping), f"{label} lacks source_ports")
    owner = ports.get(field, ports.get(field.upper()))
    _require(isinstance(owner, Mapping), f"{label} source_ports lacks field {field}")
    value = owner.get(order, owner.get(str(order)))
    _require(isinstance(value, Mapping), f"{label} source_ports lacks {field} Z{order}")
    return value


def _validate_patch_history(
    nodes: Sequence[Mapping[str, Any]],
    chart: Mapping[str, Any],
    post_ports: Mapping[tuple[str, int], Mapping[str, int]],
    packet_value: Mapping[str, Any],
    *,
    physical_measure_node: int,
) -> tuple[dict[tuple[str, int], dict[str, int]], dict[str, Any]]:
    packet = _unwrap_patch_receipt(packet_value)
    _require(packet.get("passed") is True, "accepted patch history is not passed")
    _require(packet.get("source_only") is True, "accepted patch history is not source-only")
    _require(packet.get("quadratures_replayed") is False,
             "accepted patch history reports quadrature replay")
    _require(packet.get("accepted_builders_or_owners_replayed") is False,
             "accepted patch history reports builder replay")
    _require(packet.get("accepted_checkers_replayed") is False,
             "accepted patch history reports checker replay")
    _require(packet.get("accepted_reducers_or_providers_replayed") is False,
             "accepted patch history reports provider replay")
    incoming = _chart_container(chart, "incoming")
    history = _chart_container(chart, "history")
    contributions = _chart_container(chart, "contributions", "contribution")
    memory = _chart_container(chart, "memory")
    density = _density_container(chart)
    endpoint_rows = _patch_endpoint_rows(chart)
    packet_rows = _packet_rows(packet)
    ports: dict[tuple[str, int], dict[str, int]] = {}
    evidence: list[dict[str, Any]] = []
    for field in FIELDS:
        drows = _field_rows(density, field, "patch density", len(nodes))
        irows = _field_rows(incoming, field, "patch incoming", len(nodes))
        hrows = _field_rows(history, field, "patch history", len(nodes))
        crows = _field_rows(contributions, field, "patch contribution", len(nodes))
        mroot = _handle(_field_raw(memory, field), f"patch memory.{field}", len(nodes))
        for order in ORDERS:
            key = (field, order)
            density_root = drows[order]
            incoming_root = irows[order]
            history_root = hrows[order]
            contribution_root = crows[order]
            endpoint_root = endpoint_rows.get(key)
            _require(incoming_root == PATCH_INCOMING_ROOTS[field][order],
                     f"patch {field} Z{order} incoming handle changed")
            incoming_row = _node(nodes, incoming_root, "patch incoming")
            _require(_op(nodes, incoming_root) == "function_substitution",
                     f"patch {field} Z{order} incoming is not a substitution")
            _require(_handle(incoming_row.get("expression"), "patch incoming expression", len(nodes))
                     == post_ports[key]["history"],
                     f"patch {field} Z{order} incoming is not postrestore history")
            _require(_handle(incoming_row.get("variable"), "patch incoming variable", len(nodes))
                     == POST_COORDINATE,
                     f"patch {field} Z{order} incoming variable changed")
            value = _handle(incoming_row.get("value"), "patch incoming endpoint", len(nodes))
            _require(_exact(nodes, value) == 1,
                     f"patch {field} Z{order} incoming endpoint is not exact one")
            _require(_guarded(incoming_row), f"patch {field} Z{order} incoming is unguarded")
            if endpoint_root is not None:
                _require(endpoint_root == _packet_field(packet_rows[key], "history_endpoint_root", len(nodes)),
                         f"patch {field} Z{order} endpoint root changed")
            packet_row = packet_rows[key]
            expected = {
                "density_root": density_root,
                "history_root": history_root,
                "incoming_root": incoming_root,
                "memory_root": mroot,
                "contribution_root": contribution_root,
            }
            for name, expected_value in expected.items():
                actual = _packet_field(packet_row, name, len(nodes))
                _require(actual == expected_value,
                         f"accepted patch {field} Z{order} {name} changed")
            _require(packet_row.get("status", "proved") in {"proved", "passed"},
                     f"accepted patch {field} Z{order} is not proved")
            _require(_fraction(packet_row.get("rate"), f"accepted patch {field} Z{order} rate")
                     == RATES[field],
                     f"accepted patch {field} Z{order} rate changed")
            _require(packet_row.get("first_native_y_formula")
                     == "D_y H(1)=density(1)-rate*incoming",
                     f"accepted patch {field} Z{order} derivative formula changed")
            identity = packet_row.get("history_endpoint_identity")
            _require(isinstance(identity, Mapping),
                     f"accepted patch {field} Z{order} lacks endpoint identity")
            _require(identity.get("formula") == "H(1)=incoming",
                     f"accepted patch {field} Z{order} endpoint formula changed")
            wrapper = identity.get("endpoint_wrapper")
            if isinstance(wrapper, Mapping):
                _require(_packet_field(wrapper, "expression", len(nodes)) == history_root,
                         f"accepted patch {field} Z{order} endpoint expression changed")
                _require(_packet_field(wrapper, "variable", len(nodes)) == PATCH_COORDINATE,
                         f"accepted patch {field} Z{order} endpoint variable changed")
                _require(_packet_field(wrapper, "endpoint_node", len(nodes)) == ENDPOINT_ONE,
                         f"accepted patch {field} Z{order} endpoint value changed")
                _require(wrapper.get("guarded") is True,
                         f"accepted patch {field} Z{order} endpoint is unguarded")
            # The accepted first-native source is consumed as a saved graph
            # receipt.  Validate each referenced node exactly; do not rebuild
            # or execute the old operator.
            density_endpoint = _packet_field(
                packet_row, "density_at_endpoint_root", len(nodes)
            )
            density_endpoint_row = _node(
                nodes, density_endpoint, "accepted patch density endpoint"
            )
            _require(_op(nodes, density_endpoint) == "function_substitution",
                     f"accepted patch {field} Z{order} density endpoint is not a substitution")
            _require(
                _handle(density_endpoint_row.get("expression"),
                        "accepted patch density endpoint expression", len(nodes))
                == density_root,
                f"accepted patch {field} Z{order} density endpoint expression changed",
            )
            _require(
                _handle(density_endpoint_row.get("variable"),
                        "accepted patch density endpoint variable", len(nodes))
                == PATCH_COORDINATE,
                f"accepted patch {field} Z{order} density endpoint variable changed",
            )
            endpoint_value = _handle(
                density_endpoint_row.get("value"),
                "accepted patch density endpoint value",
                len(nodes),
            )
            _require(_exact(nodes, endpoint_value) == 1,
                     f"accepted patch {field} Z{order} density endpoint is not one")
            _require(_guarded(density_endpoint_row),
                     f"accepted patch {field} Z{order} density endpoint is unguarded")
            _require(density_endpoint_row.get("source_native_substitution") is True,
                     f"accepted patch {field} Z{order} density endpoint lacks source guard")
            coefficient = _packet_field(packet_row, "rate_coefficient_root", len(nodes))
            _require(_exact(nodes, coefficient) == -RATES[field],
                     f"accepted patch {field} Z{order} rate coefficient changed")
            coefficient_row = _node(nodes, coefficient, "accepted patch rate coefficient")
            _require(coefficient_row.get("source_formula") == "-rate",
                     f"accepted patch {field} Z{order} coefficient source changed")
            incoming_term = _packet_field(packet_row, "incoming_term_root", len(nodes))
            incoming_term_row = _node(nodes, incoming_term, "accepted patch incoming term")
            _require(_op(nodes, incoming_term) in {"product", "mul"},
                     f"accepted patch {field} Z{order} incoming term is not a product")
            _require(_args(nodes, incoming_term) == [coefficient, incoming_root],
                     f"accepted patch {field} Z{order} incoming term factors changed")
            _require(incoming_term_row.get("incoming_retained") is True,
                     f"accepted patch {field} Z{order} incoming was not retained")
            _require(incoming_term_row.get("rate_zero_product_retained") is (RATES[field] == 0),
                     f"accepted patch {field} Z{order} rate-zero retention changed")
            _require(incoming_term_row.get("source_formula") == "-rate*incoming",
                     f"accepted patch {field} Z{order} incoming source changed")
            first_native = _packet_field(packet_row, "first_native_y_root", len(nodes))
            first_native_row = _node(nodes, first_native, "accepted patch first native-y")
            _require(_op(nodes, first_native) in {"sum", "add"},
                     f"accepted patch {field} Z{order} first native-y is not a sum")
            _require(_args(nodes, first_native) == [density_endpoint, incoming_term],
                     f"accepted patch {field} Z{order} first native-y factors changed")
            _require(first_native_row.get("physical_Jacobian_at_endpoint") == 1,
                     f"accepted patch {field} Z{order} endpoint Jacobian changed")
            _require(first_native_row.get("physical_Jacobian_node") == physical_measure_node,
                     f"accepted patch {field} Z{order} Jacobian source changed")
            _require(first_native_row.get("source_formula")
                     == "D_y H(1)=density(1)-rate*incoming",
                     f"accepted patch {field} Z{order} native-y source changed")
            ports[key] = {
                "density": density_root,
                "history": history_root,
                "incoming": incoming_root,
                "memory": mroot,
                "contribution": contribution_root,
                "history_endpoint": _packet_field(packet_row, "history_endpoint_root", len(nodes)),
                "density_at_endpoint": density_endpoint,
                "rate_coefficient": coefficient,
                "incoming_term": incoming_term,
                "first_native_y": first_native,
            }
            evidence.append({
                "field": field,
                "ordinary_Z_order": order,
                "rate": _rate_text(RATES[field]),
                "density_root": density_root,
                "history_root": history_root,
                "incoming_root": incoming_root,
                "memory_root": mroot,
                "contribution_root": contribution_root,
                "history_endpoint_root": ports[key]["history_endpoint"],
                "accepted_patch_native": {
                    "first_native_y_function_handle": ports[key]["first_native_y"],
                    "density_at_endpoint_handle": ports[key]["density_at_endpoint"],
                    "rate_coefficient_handle": ports[key]["rate_coefficient"],
                    "incoming_term_handle": ports[key]["incoming_term"],
                    "history_endpoint_function_handle": ports[key]["history_endpoint"],
                    "formula": "D_y H(1)=density(1)-rate*incoming",
                    "source_only_receipt": True,
                },
            })
    return ports, {
        "schema": packet.get("schema"),
        "rows": evidence,
        "source_digest": canonical_digest(packet),
        "accepted_packet_source_only": True,
        "full_topology": True,
    }


def _packet_field(value: Mapping[str, Any], name: str, count: int) -> int:
    candidate = value.get(name)
    return _handle(candidate, name, count)


def _patch_endpoint_rows(chart: Mapping[str, Any]) -> dict[tuple[str, int], int]:
    endpoints = chart.get("endpoints")
    if not isinstance(endpoints, Mapping):
        return {}
    rm = endpoints.get("Rm", endpoints.get("rm"))
    if not isinstance(rm, Mapping):
        return {}
    correction = rm.get("correction", rm.get("history"))
    if not isinstance(correction, Mapping):
        return {}
    result: dict[tuple[str, int], int] = {}
    for field in FIELDS:
        raw = _field_raw(correction, field)
        if raw is None:
            return {}
        rows = _rows(raw, f"patch endpoint {field}", 10**9)
        for order in ORDERS:
            result[(field, order)] = rows[order]
    return result


def _validate_density_admission(
    density_value: Mapping[str, Any],
    nodes: Sequence[Mapping[str, Any]],
    expected_nodes_digest: str,
    post_ports: Mapping[tuple[str, int], Mapping[str, int]],
    patch_ports: Mapping[tuple[str, int], Mapping[str, int]],
) -> dict[str, Any]:
    _require(density_value.get("passed") is True, "density admission is not passed")
    _require(density_value.get("source_only") is True, "density admission is not source-only")
    for key in ("mathematical_operator_executed", "numerical_values_computed", "full_goal_complete"):
        if key in density_value:
            _require(density_value.get(key) is False,
                     f"density admission reports forbidden {key}")
    archive_sha = density_value.get("archive_sha256", density_value.get("source_archive_sha256"))
    _require(str(archive_sha).lower() == LOOP_DENSITY_ARCHIVE_SHA256,
             "density admission archive SHA-256 changed")
    _require(str(density_value.get("nodes_digest", "")).lower() == expected_nodes_digest.lower(),
             "density admission nodes digest changed")
    if density_value.get("graph_node_count") is not None:
        _require(density_value.get("graph_node_count") == len(nodes),
                 "density admission graph node count changed")
    _require(density_value.get("ordinary_orders") == list(ORDERS),
             "density admission ordinary orders changed")
    values = _sequence(density_value.get("rows"), "density admission rows")
    selected: dict[tuple[str, int], Mapping[str, Any]] = {}
    for value in values:
        _require(isinstance(value, Mapping), "density admission row is not a mapping")
        field = value.get("field")
        order = value.get("order", value.get("ordinary_Z_order"))
        if field not in FIELDS:
            continue
        _require(order in ORDERS, "density admission row has unsupported ordinary order")
        key = (field, order)
        _require(key not in selected, "density admission has duplicate selected row")
        selected[key] = value
    _require(set(selected) == {(field, order) for field in FIELDS for order in ORDERS},
             "density admission does not contain exactly ten selected rows")
    rows: list[dict[str, Any]] = []
    row_source_equality: dict[str, bool | None] = {}
    source_equality_checked = 0
    for field in FIELDS:
        for order in ORDERS:
            key = (field, order)
            value = selected[key]
            _require(value.get("status") == "proved",
                     f"density admission {field} Z{order} is not proved")
            post = value.get("postrestore")
            patch = value.get("patch")
            _require(isinstance(post, Mapping) and isinstance(patch, Mapping),
                     f"density admission {field} Z{order} lacks postrestore/patch roots")
            post_root = _handle(post.get("node"), "density admission postrestore node", len(nodes))
            patch_root = _handle(patch.get("node"), "density admission patch node", len(nodes))
            _require(post_root == post_ports[key]["density"],
                     f"density admission {field} Z{order} postrestore root changed")
            _require(patch_root == patch_ports[key]["density"],
                     f"density admission {field} Z{order} patch root changed")
            post_value = post.get("value")
            patch_value = patch.get("value")
            _require(isinstance(post_value, Mapping) and isinstance(patch_value, Mapping),
                     f"density admission {field} Z{order} lacks admitted canonical source forms")
            value_equal = canonical_digest(post_value) == canonical_digest(patch_value)
            _require(value_equal, f"density admission {field} Z{order} source forms differ")
            source_equality_checked += 1
            row_source_equality[f"{field}:{order}"] = value_equal
            rows.append({
                "field": field,
                "ordinary_Z_order": order,
                "postrestore_density_root": post_root,
                "patch_density_root": patch_root,
                "status": "proved",
                "value_consumed": False,
                "source_value_equality_checked": value_equal is not None,
            })
    return {
        "source_archive_sha256": str(archive_sha).lower(),
        "nodes_digest": expected_nodes_digest,
        "graph_node_count": len(nodes),
        "rows": rows,
        "bounded_rows_consumed": len(rows),
        "values_consumed": False,
        "row_source_equality": row_source_equality,
        "source_value_equality_checked": source_equality_checked == len(rows),
        "source_value_equality_rows": source_equality_checked,
        "source_only": True,
        "digest": canonical_digest(density_value),
    }


def _validate_expected_digest(nodes: Sequence[Mapping[str, Any]], expected: Any) -> str:
    _require(isinstance(expected, str) and len(expected) == 64,
             "expected_nodes_digest must be SHA-256")
    try:
        int(expected, 16)
    except ValueError as exc:
        raise ValueError("expected_nodes_digest must be SHA-256") from exc
    actual = canonical_digest(list(nodes))
    _require(actual == expected.lower(), "nodes digest does not match expected_nodes_digest")
    return actual


def _append(delta: list[dict[str, Any]], row: Mapping[str, Any]) -> int:
    index = len(delta)
    delta.append(dict(row))
    return index


def derive(
    nodes: Sequence[Mapping[str, Any]],
    left_chart: Mapping[str, Any],
    right_chart: Mapping[str, Any],
    accepted_patch_history: Mapping[str, Any],
    density_admission: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Consume bounded saved ports and append first-native postrestore rows.

    ``nodes`` is copied only as a shallow list.  No graph node is evaluated.
    The returned ``graph_delta`` contains forty source rows: one guarded
    ``density(1)`` substitution, one exact ``-rate`` coefficient, one
    retained incoming product, and one native-y sum for each of ten rows.
    """

    base = _sequence(nodes, "nodes")
    _require(all(isinstance(row, Mapping) for row in base), "nodes must contain mappings")
    for value, label in (
        (left_chart, "left_chart"),
        (right_chart, "right_chart"),
        (accepted_patch_history, "accepted_patch_history"),
        (density_admission, "density_admission"),
    ):
        _require(isinstance(value, Mapping), f"{label} must be a mapping")
    actual_digest = _validate_expected_digest(base, expected_nodes_digest)
    left, outlet = _unwrap_left(left_chart)
    right = _unwrap_right(right_chart)
    post_ports, post_proof = _validate_postrestore(base, left)
    patch_coordinate = _coordinate(base, right, "patch")
    _require(patch_coordinate == PATCH_COORDINATE, "patch coordinate changed")
    _domain(right, (1, 2948), "patch")
    patch_physical = _physical_info(base, right, patch_coordinate, "patch", postrestore=False)
    patch_density = _density_container(right)
    patch_incoming = _chart_container(right, "incoming")
    patch_history = _chart_container(right, "history")
    patch_contributions = _chart_container(right, "contributions", "contribution")
    patch_memory = _chart_container(right, "memory")
    patch_ports_layout: dict[tuple[str, int], dict[str, int]] = {}
    for field in FIELDS:
        drows = _field_rows(patch_density, field, "patch density", len(base))
        irows = _field_rows(patch_incoming, field, "patch incoming", len(base))
        hrows = _field_rows(patch_history, field, "patch history", len(base))
        crows = _field_rows(patch_contributions, field, "patch contribution", len(base))
        mroot = _handle(_field_raw(patch_memory, field), f"patch memory.{field}", len(base))
        for order in ORDERS:
            patch_ports_layout[(field, order)] = {
                "density": drows[order],
                "history": hrows[order],
                "incoming": irows[order],
                "memory": mroot,
                "contribution": crows[order],
            }
    patch_ports, patch_proof = _validate_patch_history(
        base,
        right,
        post_ports,
        accepted_patch_history,
        physical_measure_node=patch_physical["node"],
    )
    # The packet and chart must describe the same patch layout.  This check is
    # intentionally exact and does not accept ancestry labels or containment.
    for key, port in patch_ports_layout.items():
        for name in ("density", "history", "incoming", "memory", "contribution"):
            _require(patch_ports[key][name] == port[name],
                     f"patch port {key!r} {name} differs from chart")
    density_proof = _validate_density_admission(
        density_admission,
        base,
        expected_nodes_digest.lower(),
        post_ports,
        patch_ports,
    )

    # If a caller supplied an explicit outlet/Rm-correction port map, consume
    # it as an equality receipt.  The patch substitutions remain authoritative
    # and are still checked above when no optional wrapper is present.
    outlet_proof: dict[str, Any] = {"provided": False, "checked": 0}
    if isinstance(outlet, Mapping):
        candidate = outlet.get("Rm_correction", outlet.get("correction", outlet.get("history")))
        if isinstance(candidate, Mapping):
            checked = 0
            for field in FIELDS:
                raw = _field_raw(candidate, field)
                if raw is None:
                    continue
                rows = _rows(raw, f"outlet {field}", len(base))
                for order in ORDERS:
                    _require(rows[order] == patch_ports[(field, order)]["incoming"],
                             f"outlet {field} Z{order} differs from patch incoming")
                    checked += 1
            _require(checked == 10, "outlet correction map is incomplete")
            outlet_proof = {"provided": True, "checked": checked, "exact_match": True}

    row_by_key = {
        (row["field"], row["ordinary_Z_order"]): row
        for row in post_proof["rows"]
    }
    delta: list[dict[str, Any]] = []
    first_native: dict[str, dict[str, int]] = {field: {} for field in FIELDS}
    density_endpoint: dict[str, dict[str, int]] = {field: {} for field in FIELDS}
    incoming_term: dict[str, dict[str, int]] = {field: {} for field in FIELDS}

    for field in FIELDS:
        for order in ORDERS:
            key = (field, order)
            post = post_ports[key]
            patch = patch_ports[key]
            rate = RATES[field]
            density_equal = density_proof["row_source_equality"].get(f"{field}:{order}")
            endpoint_density = _append(delta, {
                "operation": "function_substitution",
                "expression": post["density"],
                "variable": POST_COORDINATE,
                "value": ENDPOINT_ONE,
                "Z_independent_substitution": True,
                "source_native_substitution": True,
                "guarded": True,
                "source_formula": "density_post(1)",
                "field": field,
                "ordinary_Z_order": order,
                "original_chart": "postrestore",
            })
            density_endpoint[field][str(order)] = len(base) + endpoint_density
            coefficient = _append(delta, {
                "operation": "exact_rational",
                "numerator": -rate.numerator,
                "denominator": rate.denominator,
                "source_formula": "-rate",
                "rate": _rate_text(rate),
                "field": field,
                "ordinary_Z_order": order,
            })
            incoming_product = _append(delta, {
                "operation": "product",
                "arguments": [len(base) + coefficient, patch["incoming"]],
                "source_formula": "-rate*post_history_endpoint",
                "incoming_retained": True,
                "rate_zero_product_retained": rate == 0,
                "post_history_endpoint_source": patch["incoming"],
                "field": field,
                "ordinary_Z_order": order,
            })
            incoming_term[field][str(order)] = len(base) + incoming_product
            derivative = _append(delta, {
                "operation": "sum",
                "arguments": [
                    len(base) + endpoint_density,
                    len(base) + incoming_product,
                ],
                "source_formula": "D_y H(1)=density(1)-rate*incoming",
                "normal_formula": "D_y H(1)=density_post(1)-rate*post_history_endpoint",
                "physical_Jacobian_at_endpoint": 1,
                "physical_Jacobian_applied_once": True,
                "post_history_endpoint_source": patch["incoming"],
                "patch_incoming_root": patch["incoming"],
                "field": field,
                "ordinary_Z_order": order,
            })
            first_native[field][str(order)] = len(base) + derivative
            row_by_key[key]["post_density_endpoint_root"] = len(base) + endpoint_density
            row_by_key[key]["post_native_y_root"] = len(base) + derivative
            row_by_key[key]["post_endpoint_source_root"] = patch["incoming"]
            row_by_key[key]["patch_history_endpoint_root"] = patch["history_endpoint"]
            row_by_key[key]["accepted_patch_native_function_handle"] = patch["first_native_y"]
            row_by_key[key]["accepted_patch_density_endpoint_root"] = patch["density_at_endpoint"]
            row_by_key[key]["accepted_patch_rate_coefficient_root"] = patch["rate_coefficient"]
            row_by_key[key]["accepted_patch_incoming_term_root"] = patch["incoming_term"]
            row_by_key[key]["post_history_body"] = post["history"]
            row_by_key[key]["patch_history_body"] = patch["history"]

    full_nodes = base + delta
    result_rows: list[dict[str, Any]] = []
    for field in FIELDS:
        for order in ORDERS:
            key = (field, order)
            post = post_ports[key]
            patch = patch_ports[key]
            density_equal = density_proof["row_source_equality"][f"{field}:{order}"]
            _require(density_equal is True, "paired density source identity missing")
            result_rows.append({
                "field": field,
                "ordinary_Z_order": order,
                "rate": _rate_text(RATES[field]),
                "post_history_body": post["history"],
                "post_history_endpoint_source": patch["incoming"],
                "post_endpoint_source_is_patch_incoming": True,
                "patch_history_body": patch["history"],
                "patch_history_endpoint_function_handle": patch["history_endpoint"],
                "accepted_patch_native_function_handle": patch["first_native_y"],
                "accepted_patch_density_endpoint_handle": patch["density_at_endpoint"],
                "accepted_patch_rate_coefficient_handle": patch["rate_coefficient"],
                "accepted_patch_incoming_term_handle": patch["incoming_term"],
                "post_density_endpoint_definition": density_endpoint[field][str(order)],
                "post_native_y_definition": first_native[field][str(order)],
                "post_new_native_y_handle": first_native[field][str(order)],
                "accepted_patch_native_y_handle": patch["first_native_y"],
                "incoming_term_definition": incoming_term[field][str(order)],
                "normal_formula": "D_y H(1)=density_post(1)-rate*post_history_endpoint",
                "source_formula": "D_y H(1)=density(1)-rate*incoming",
                "physical_J_at_endpoint": 1,
                "rate_zero_incoming_retained": RATES[field] == 0,
                "native_pairing": {
                    "same_density_source": density_equal is True,
                    "density_source_value_equality_checked": density_equal is not None,
                    "post_density_handle_binding_exact": True,
                    "accepted_patch_density_handle_binding_exact": True,
                    "same_rate": True,
                    "rate": _rate_text(RATES[field]),
                    "same_incoming_source": post["history"] == _node(
                        nodes, patch["incoming"], "paired patch incoming"
                    ).get("expression"),
                    "reason": (
                        "paired because the bounded density admission records equal "
                        "postrestore/patch source data (when present), while the "
                        "accepted patch density endpoint, exact -rate coefficient, "
                        "and retained incoming term are each graph-equal to the "
                        "saved patch native source bindings"
                    ),
                },
                "status": "source_identity_proved",
            })

    static = static_source_proof()
    source_proof = {
        "postrestore": post_proof,
        "patch": patch_proof,
        "density_admission": density_proof,
        "outlet": outlet_proof,
        "paired_rows": len(result_rows),
        "ordinary_orders": list(ORDERS),
        "fields": list(FIELDS),
        "exact_endpoint": "p84528=1",
        "physical_J_endpoint": 1,
        "full_topology": True,
        "no_ancestry_or_contains_shortcut": True,
        "all_factors_and_Jacobian_once": True,
    }
    return {
        "schema": "original_Rm_history_seam_identity_v1",
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "actual_saved_graph_run": False,
        "actual_source_graph_inspected": True,
        "accepted_graph_replayed": False,
        "accepted_operator_replayed": False,
        "provider_execution": False,
        "history_archive_decoded": False,
        "numerical_values_computed": False,
        "full_goal_complete": False,
        "modified_moment_value_native_rows_claimed": False,
        "expected_nodes_digest": expected_nodes_digest.lower(),
        "nodes_digest": actual_digest,
        "module_sha256": static["source_bytes_sha256"],
        "ordinary_orders": list(ORDERS),
        "fields": list(FIELDS),
        "source_proof": source_proof,
        "sourceproof": source_proof,
        "rows": result_rows,
        "postrestore_endpoint": {
            "coordinate_node": POST_COORDINATE,
            "endpoint_node": ENDPOINT_ONE,
            "formula": "H_post(1)=post_history_endpoint",
            "physical_J_at_endpoint": 1,
        },
        "patch_inlet": {
            "coordinate_node": PATCH_COORDINATE,
            "incoming_roots": {
                field: {str(order): patch_ports[(field, order)]["incoming"] for order in ORDERS}
                for field in FIELDS
            },
            "formula": "patch_incoming=postrestore_Rm_correction(1)",
            "endpoint_source_is_current_loop_history": True,
        },
        "density_admission": density_proof,
        "graph_delta": {
            "base_graph_node_count": len(base),
            "base_nodes_sha256": actual_digest,
            "appended_nodes": delta,
            "appended_node_count": len(delta),
            "full_node_count": len(full_nodes),
            "full_nodes_sha256": canonical_digest(full_nodes),
            "post_density_endpoint_roots": density_endpoint,
            "post_native_y_history_derivative_roots": first_native,
        },
        "static_proof": static,
        "limitations": [
            "Source identities only; no density or derivative value is selected.",
            "No modified-moment VALUE/native rows are claimed.",
            "No physical Q_y, pressure, or global closure identity is claimed.",
            "The caller must run independent AST/LF/single-EOF checks on this module.",
        ],
        "scope": {
            "source_graph_only": True,
            "actual_source_graph_checked": True,
            "accepted_graph_replayed": False,
            "no_archive_decode": True,
            "no_numeric_integration": True,
            "no_numeric_values_selected": True,
            "exact_postrestore_topology_verified": True,
            "exact_patch_incoming_topology_verified": True,
            "first_native_y_history_source_installed": True,
            "rate_zero_pressure_incoming_retained": True,
            "modified_moment_value_native_rows_claimed": False,
            "full_goal_complete": False,
        },
    }


__all__ = [
    "FIELDS",
    "ORDERS",
    "RATES",
    "POST_COORDINATE",
    "PATCH_COORDINATE",
    "PATCH_INCOMING_ROOTS",
    "canonical_digest",
    "static_source_proof",
    "derive",
]
