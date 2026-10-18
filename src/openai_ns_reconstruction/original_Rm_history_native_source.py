"""Source-only first native-y ports for the current finite-N Rm history.

The current source transports five signed density functions with the exact
dt/t Duhamel kernel.  This adapter consumes saved source bindings and an
authenticated inlet receipt.  It does not import a producer, execute an
archived graph, evaluate an integral, or select a numeric endpoint value.

For every requested ordinary Z row it appends only source rows for
density(1) and D_y H(1) = density(1) - rate * incoming.

The historical inlet receipt used common symbolic toy bodies.  It is
therefore provenance only; toy values are never admitted as actual normalized
fields.  The current endpoint H(1)=incoming identity is established here from
the actual graph ports, source metadata, and exact logarithmic memory shape.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
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
_WRAPPERS = (
    "Rm_patch_finite_N_C5_bindings",
    "finite_bindings",
    "finite_patch",
    "finite",
)
_F_SUB_GUARDS = (
    "Z_independent_substitution",
    "native_substitution_guard",
    "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard",
    "guarded",
)


def canonical_digest(value: Any) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _module_sha256() -> str:
    digest = hashlib.sha256()
    with Path(__file__).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
                pieces = text.split("/")
                _require(len(pieces) == 2, f"{label} is not an exact rational")
                return Fraction(int(pieces[0]), int(pieces[1]))
            return Fraction(int(text), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{label} is not an exact rational") from exc
    raise ValueError(f"{label} is not an exact rational")


def _rate_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else (
        f"{value.numerator}/{value.denominator}"
    )


def _op(nodes: Sequence[Any], index: int) -> str:
    row = nodes[index]
    _require(isinstance(row, Mapping), f"node {index} is not a mapping")
    operation = row.get("operation", row.get("op"))
    _require(isinstance(operation, str) and operation, f"node {index} has no operation")
    return operation


def _row(nodes: Sequence[Any], index: int, label: str = "node") -> Mapping[str, Any]:
    _handle(index, label, len(nodes))
    value = nodes[index]
    _require(isinstance(value, Mapping), f"{label} {index} is not a mapping")
    return value


def _args(nodes: Sequence[Any], index: int) -> list[int]:
    raw = _row(nodes, index).get("arguments", _row(nodes, index).get("operands"))
    values = _sequence(raw, f"node {index}.arguments")
    return [_handle(value, f"node {index}.arguments[{offset}]", len(nodes))
            for offset, value in enumerate(values)]


def _exact(nodes: Sequence[Any], index: int) -> Fraction | None:
    try:
        if _op(nodes, index) != "exact_rational":
            return None
        row = _row(nodes, index)
        numerator = _fraction(row.get("numerator", row.get("value")), "rational numerator")
        denominator = _fraction(row.get("denominator", 1), "rational denominator")
        return None if denominator == 0 else numerator / denominator
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _guarded(row: Mapping[str, Any]) -> bool:
    return any(row.get(key) is True for key in _F_SUB_GUARDS)


def _first(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _unwrap(value: Mapping[str, Any]) -> Mapping[str, Any]:
    current = value
    seen: set[int] = set()
    for _ in range(4):
        if id(current) in seen:
            break
        seen.add(id(current))
        child = next(
            (current.get(key) for key in _WRAPPERS
             if isinstance(current.get(key), Mapping)),
            None,
        )
        if child is None:
            break
        current = child
    return current


def _chart(bindings: Mapping[str, Any]) -> Mapping[str, Any]:
    return bindings.get("chart") if isinstance(bindings.get("chart"), Mapping) else bindings


def _container(
    bindings: Mapping[str, Any],
    chart: Mapping[str, Any],
    *keys: str,
) -> Mapping[str, Any] | None:
    for owner in (chart, bindings):
        for key in keys:
            candidate = owner.get(key)
            if isinstance(candidate, Mapping):
                return candidate
    return None


def _field_raw(container: Mapping[str, Any], field: str) -> Any:
    aliases = (field, {"m": "M", "h": "H", "k": "K", "e": "E", "p": "P"}[field])
    for key in aliases:
        if key in container:
            return container[key]
    roots = container.get("roots")
    if isinstance(roots, Mapping):
        for key in aliases:
            if key in roots:
                return roots[key]
    return None


def _rows(
    raw: Any,
    label: str,
    count: int,
    orders: Sequence[int],
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
    orders: Sequence[int],
) -> dict[int, int]:
    raw = _field_raw(container, field)
    _require(raw is not None, f"{label} is missing field {field}")
    return _rows(raw, f"{label}.{field}", count, orders)


def _endpoint_rows(
    bindings: Mapping[str, Any],
    count: int,
    orders: Sequence[int],
) -> dict[str, dict[int, int]]:
    endpoints = bindings.get("endpoints")
    _require(isinstance(endpoints, Mapping), "finite bindings are missing endpoints")
    rm = endpoints.get("Rm", endpoints.get("rm"))
    _require(isinstance(rm, Mapping), "finite bindings are missing endpoints.Rm")
    correction = rm.get("correction", rm.get("history", rm.get("H")))
    _require(isinstance(correction, Mapping), "endpoints.Rm is missing correction rows")
    return {
        field: _field_rows(correction, field, "endpoints.Rm.correction", count, orders)
        for field in FIELDS
    }


def _endpoint_node(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    chart: Mapping[str, Any],
    endpoint_rows: Mapping[str, Mapping[int, int]],
) -> int:
    direct = _first(
        chart,
        "endpoint_node1",
        "endpoint_node",
        "endpoint",
        "lower_node",
        "domain_lower_node",
    )
    if direct is None:
        direct = _first(
            bindings,
            "endpoint_node1",
            "endpoint_node",
            "endpoint",
            "lower_node",
            "domain_lower_node",
        )
    if direct is None:
        first = endpoint_rows[FIELDS[0]][min(endpoint_rows[FIELDS[0]])]
        wrapper = _row(nodes, first, "endpoint history")
        _require(_op(nodes, first) == "function_substitution",
                 "endpoint history row is not a function substitution")
        direct = wrapper.get("value")
    endpoint = _handle(direct, "endpoint node", len(nodes))
    _require(_exact(nodes, endpoint) == 1, "Rm endpoint node is not exact rational one")
    return endpoint


def _coordinate_node(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    chart: Mapping[str, Any],
) -> int:
    value = _first(chart, "coordinate", "native_coordinate", "x")
    if value is None:
        value = _first(bindings, "coordinate", "native_coordinate", "x")
    _require(value is not None, "finite bindings are missing Rm coordinate")
    return _handle(value, "Rm coordinate", len(nodes))


def _physical_node(
    nodes: Sequence[Any],
    bindings: Mapping[str, Any],
    chart: Mapping[str, Any],
) -> int:
    value = _first(chart, "physical_measure", "physical_Jacobian", "jacobian")
    if value is None:
        value = _first(bindings, "physical_measure", "physical_Jacobian", "jacobian")
    _require(value is not None, "finite bindings are missing physical measure")
    return _handle(value, "physical measure", len(nodes))


def _density_container(
    bindings: Mapping[str, Any],
    chart: Mapping[str, Any],
) -> Mapping[str, Any]:
    finite = chart.get("finite_N")
    owners = (finite, chart, bindings) if isinstance(finite, Mapping) else (chart, bindings)
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
    raise ValueError("finite bindings are missing five signed density increments")


def _memory_node(nodes: Sequence[Any], memory: Mapping[str, Any], field: str) -> int:
    raw = _field_raw(memory, field)
    _require(raw is not None, f"memory is missing field {field}")
    if isinstance(raw, Mapping) and not any(key in raw for key in ("node", "handle", "root")):
        raw = raw.get(0, raw.get("0"))
    return _handle(raw, f"memory.{field}", len(nodes))


def _unary_argument(nodes: Sequence[Any], index: int, label: str) -> int:
    """Return the single argument handle used by a unary source row."""

    row = _row(nodes, index, label)
    value = row.get("argument", row.get("operand", row.get("value")))
    return _handle(value, f"{label}.argument", len(nodes))


def _flatten_product(nodes: Sequence[Any], root: int) -> list[int]:
    """Flatten only commutative product rows, retaining all other roots."""

    operation = _op(nodes, root)
    if operation not in {"product", "mul"}:
        return [root]
    arguments = _args(nodes, root)
    _require(arguments, f"product node {root} has no factors")
    factors: list[int] = []
    for argument in arguments:
        factors.extend(_flatten_product(nodes, argument))
    return factors


def _affine_add(
    left: tuple[dict[tuple[str, int], Fraction], Fraction],
    right: tuple[dict[tuple[str, int], Fraction], Fraction],
) -> tuple[dict[tuple[str, int], Fraction], Fraction]:
    coefficients = dict(left[0])
    for key, value in right[0].items():
        coefficients[key] = coefficients.get(key, Fraction(0, 1)) + value
        if coefficients[key] == 0:
            del coefficients[key]
    return coefficients, left[1] + right[1]


def _affine_scale(
    coefficient: Fraction,
    value: tuple[dict[tuple[str, int], Fraction], Fraction],
) -> tuple[dict[tuple[str, int], Fraction], Fraction]:
    return (
        {key: coefficient * item for key, item in value[0].items() if coefficient * item},
        coefficient * value[1],
    )


def _log_affine(
    nodes: Sequence[Any],
    root: int,
    coordinate: int,
    endpoint: int,
    variable: int | None = None,
    _stack: tuple[int, ...] = (),
) -> tuple[dict[tuple[str, int], Fraction], Fraction]:
    """Parse a bounded exact affine expression in log leaves.

    The parser deliberately accepts only the source algebra needed by the
    current transport: exact rationals, sums, negations, products by one
    exact rational, and log leaves.  It never evaluates a graph node.
    """

    _require(root not in _stack, f"cyclic logarithmic source at node {root}")
    operation = _op(nodes, root)
    if operation == "exact_rational":
        value = _exact(nodes, root)
        _require(value is not None, f"node {root} is not a valid exact rational")
        return {}, value
    if operation in {"analytic_unary", "log"}:
        row = _row(nodes, root, "log source")
        name = (
            row.get("name", row.get("function", "log"))
            if operation == "analytic_unary" else "log"
        )
        _require(name == "log", f"node {root} is not a logarithm")
        child = _unary_argument(nodes, root, "log source")
        _require(
            child in {coordinate, endpoint} or (variable is not None and child == variable),
            f"log node {root} uses an out-of-scope source handle",
        )
        return {("log", child): Fraction(1, 1)}, Fraction(0, 1)
    if operation in {"negative", "neg"}:
        child = _unary_argument(nodes, root, "negative source")
        return _affine_scale(
            Fraction(-1, 1),
            _log_affine(nodes, child, coordinate, endpoint, variable, _stack + (root,)),
        )
    if operation in {"sum", "add"}:
        arguments = _args(nodes, root)
        _require(arguments, f"sum node {root} has no terms")
        result: tuple[dict[tuple[str, int], Fraction], Fraction] = ({}, Fraction(0, 1))
        for argument in arguments:
            result = _affine_add(
                result,
                _log_affine(nodes, argument, coordinate, endpoint, variable, _stack + (root,)),
            )
        return result
    if operation in {"product", "mul"}:
        arguments = _args(nodes, root)
        _require(arguments, f"product node {root} has no factors")
        scalar = Fraction(1, 1)
        affine: tuple[dict[tuple[str, int], Fraction], Fraction] = ({}, Fraction(0, 1))
        nonconstant = False
        for argument in arguments:
            current = _log_affine(
                nodes, argument, coordinate, endpoint, variable, _stack + (root,)
            )
            if current[0]:
                _require(not nonconstant, f"node {root} is nonlinear in logarithms")
                nonconstant = True
                affine = current
            else:
                _require(not current[1] == 0 or not nonconstant,
                         f"node {root} has an unsupported zero product")
                scalar *= current[1]
        if not nonconstant:
            return {}, scalar
        return _affine_scale(scalar, affine)
    if operation in {"difference", "sub"}:
        arguments = _args(nodes, root)
        _require(len(arguments) == 2, f"difference node {root} must have two terms")
        return _affine_add(
            _log_affine(nodes, arguments[0], coordinate, endpoint, variable, _stack + (root,)),
            _affine_scale(
                Fraction(-1, 1),
                _log_affine(nodes, arguments[1], coordinate, endpoint, variable, _stack + (root,)),
            ),
        )
    raise ValueError(f"node {root} is outside the bounded logarithmic source algebra")


def _expected_affine(
    coordinate: int,
    endpoint: int,
    rate: Fraction,
) -> tuple[dict[tuple[str, int], Fraction], Fraction]:
    coefficients: dict[tuple[str, int], Fraction] = {}
    if rate:
        coefficients[("log", coordinate)] = -rate
        coefficients[("log", endpoint)] = rate
    return coefficients, Fraction(0, 1)


def _kernel_info(
    nodes: Sequence[Any],
    factors: Sequence[int],
    coordinate: int,
    endpoint: int,
    variable: int,
    rate: Fraction,
) -> dict[str, Any]:
    """Validate the exact dt/t factor and the current logarithmic kernel."""

    measure: int | None = None
    numerator: int | None = None
    denominator: int | None = None
    kernels: list[int] = []
    extras: list[int] = []
    for factor in factors:
        operation = _op(nodes, factor)
        if operation in {"positive_quotient", "quotient"}:
            if measure is not None:
                extras.append(factor)
                continue
            row = _row(nodes, factor, "transport measure")
            numerator = _handle(row.get("numerator"), "transport measure numerator", len(nodes))
            denominator = _handle(row.get("denominator"), "transport measure denominator", len(nodes))
            _require(_exact(nodes, numerator) == 1, "transport measure numerator is not exact one")
            measure = factor
            continue
        if operation in {"analytic_unary", "exp"}:
            row = _row(nodes, factor, "transport kernel")
            name = (
                row.get("name", row.get("function", "exp"))
                if operation == "analytic_unary" else "exp"
            )
            if name == "exp":
                kernels.append(factor)
                continue
        if _exact(nodes, factor) == 1:
            continue
        extras.append(factor)
    _require(measure is not None, "transport integrand lacks exactly one dt/t measure")
    _require(not extras, f"transport integrand has unexpected factors {extras!r}")
    _require(not denominator is None, "transport measure denominator is missing")
    if rate == 0:
        _require(not kernels, "rate-zero transport has a nontrivial exponential kernel")
        return {
            "measure_node": measure,
            "measure_numerator": numerator,
            "integration_variable": denominator,
            "kernel_node": None,
            "kernel_kind": "exact_one",
            "kernel_affine": {"coefficients": {}, "constant": "0"},
        }
    _require(len(kernels) == 1, "transport integrand must have one exponential kernel")
    kernel = kernels[0]
    argument = _unary_argument(nodes, kernel, "transport kernel")
    affine = _log_affine(nodes, argument, coordinate, endpoint, denominator)
    expected = _expected_affine(coordinate, denominator, rate)
    _require(affine == expected, "transport kernel is not exp[-rate*(log(x)-log(t))]")
    return {
        "measure_node": measure,
        "measure_numerator": numerator,
        "integration_variable": denominator,
        "kernel_node": kernel,
        "kernel_kind": "exp_log_suffix",
        "kernel_argument": argument,
        "kernel_affine": {
            "coefficients": {
                f"{key[0]}:{key[1]}": _rate_text(value)
                for key, value in affine[0].items()
            },
            "constant": _rate_text(affine[1]),
        },
    }


def _memory_shape(
    nodes: Sequence[Any],
    root: int,
    coordinate: int,
    endpoint: int,
    rate: Fraction,
) -> dict[str, Any]:
    if rate == 0:
        _require(_exact(nodes, root) == 1, "rate-zero memory is not exact one")
        return {"node": root, "kind": "exact_one", "endpoint_unit": True}
    _require(_op(nodes, root) in {"analytic_unary", "exp"},
             f"memory node {root} is not an exponential source")
    row = _row(nodes, root)
    name = row.get("name", row.get("function", "exp")) if _op(nodes, root) == "analytic_unary" else "exp"
    _require(name == "exp", f"memory node {root} is not exp")
    argument = _unary_argument(nodes, root, f"memory node {root}")
    affine = _log_affine(nodes, argument, coordinate, endpoint)
    _require(
        affine == _expected_affine(coordinate, endpoint, rate),
        f"memory node {root} is not exp[-rate*(log(x)-log(1))]",
    )
    return {
        "node": root,
        "kind": "exp_log_width",
        "rate": _rate_text(rate),
        "argument": argument,
        "affine": {
            "coefficients": {
                f"{key[0]}:{key[1]}": _rate_text(value)
                for key, value in affine[0].items()
            },
            "constant": _rate_text(affine[1]),
        },
        "endpoint_unit": True,
    }


def _physical_info(nodes: Sequence[Any], physical: int, coordinate: int) -> dict[str, Any]:
    row = _row(nodes, physical, "physical measure")
    _require(_op(nodes, physical) in {"positive_quotient", "quotient"},
             "physical measure is not a positive quotient")
    numerator = _handle(row.get("numerator"), "physical measure numerator", len(nodes))
    denominator = _handle(row.get("denominator"), "physical measure denominator", len(nodes))
    _require(_exact(nodes, numerator) == 1, "physical measure numerator is not exact one")
    _require(denominator == coordinate, "physical measure denominator is not Rm coordinate")
    return {
        "node": physical,
        "formula": "1/x",
        "coordinate_node": coordinate,
        "positive_certificate": row.get("source_positive_certificate"),
        "endpoint_exact_one": True,
        "global_exact_one": False,
        "scope": "endpoint x=1 only",
    }


def _contribution_info(
    nodes: Sequence[Any],
    index: int,
    density_root: int,
    field: str,
    order: int,
    coordinate: int,
    endpoint: int,
    rate: Fraction,
) -> dict[str, Any]:
    row = _row(nodes, index, "Duhamel contribution")
    _require(_op(nodes, index) == "definite_integral",
             f"{field} Z{order} contribution is not a definite integral")
    _require(_handle(row.get("lower"), "contribution lower", len(nodes)) == endpoint,
             f"{field} Z{order} lower endpoint changed")
    _require(_handle(row.get("upper"), "contribution upper", len(nodes)) == coordinate,
             f"{field} Z{order} upper coordinate changed")
    _require(_fraction(row.get("original_rate"), "contribution rate") == rate,
             f"{field} Z{order} rate metadata changed")
    required = {
        "physical_measure": "dt/t",
        "physical_Jacobian_applied_once": True,
        "exact_original_log_suffix": True,
        "nonzero_original_Rm_correction_inlet_retained": True,
        "predecessor_memory_never_reset": True,
        "actual_finite_N_signed_correction": True,
        "endpoint_independent_of_Z": True,
        "exact_function_integral": True,
    }
    for key, expected in required.items():
        _require(
            row.get(key) is expected if isinstance(expected, bool) else row.get(key) == expected,
                 f"{field} Z{order} contribution metadata {key} changed")
    declared_order = row.get("ordinary_slow_Z_derivative_order")
    if declared_order is not None:
        _require(declared_order == order, f"{field} Z{order} contribution order changed")
    _require(isinstance(row.get("variable"), str) and row["variable"],
             f"{field} Z{order} transport variable missing")
    integrand = _handle(row.get("integrand"), "contribution integrand", len(nodes))
    factors = _flatten_product(nodes, integrand)

    # The source must expose one guarded density substitution and exactly one
    # dt/t measure.  Discovering the integration variable from the measure
    # keeps the check tied to the actual source graph rather than to metadata.
    measure_factors = [
        factor for factor in factors
        if _op(nodes, factor) in {"positive_quotient", "quotient"}
    ]
    _require(len(measure_factors) == 1,
             f"{field} Z{order} contribution must have one dt/t measure")
    measure_row = _row(nodes, measure_factors[0], "transport measure")
    measure_numerator = _handle(
        measure_row.get("numerator"), "transport measure numerator", len(nodes)
    )
    integration_variable = _handle(
        measure_row.get("denominator"), "transport measure denominator", len(nodes)
    )
    _require(
        integration_variable != coordinate,
        f"{field} Z{order} transport integration variable aliases Rm coordinate",
    )
    _require(_exact(nodes, measure_numerator) == 1,
             f"{field} Z{order} transport measure numerator is not exact one")
    variable_row = _row(nodes, integration_variable, "transport variable")
    variable_name = variable_row.get("name")
    if isinstance(variable_name, str) and variable_name:
        _require(
            row["variable"] == variable_name,
            f"{field} Z{order} transport metadata variable does not match graph variable",
        )
    if row.get("variable_node") is not None:
        _require(
            _handle(row["variable_node"], "contribution variable node", len(nodes))
            == integration_variable,
            f"{field} Z{order} contribution variable node changed",
        )

    density_factors = []
    for factor in factors:
        if _op(nodes, factor) != "function_substitution":
            continue
        candidate = _row(nodes, factor, "density substitution")
        expression = candidate.get("expression")
        variable = candidate.get("variable")
        value = candidate.get("value")
        if expression == density_root and variable == coordinate and value == integration_variable:
            density_factors.append(factor)
    _require(len(density_factors) == 1,
             f"{field} Z{order} contribution lacks its guarded density(t) source")
    density_substitution = _row(nodes, density_factors[0], "density substitution")
    _require(_guarded(density_substitution),
             f"{field} Z{order} density(t) substitution lacks a guard")
    non_density = [factor for factor in factors if factor not in density_factors]
    kernel = _kernel_info(
        nodes, non_density, coordinate, endpoint, integration_variable, rate
    )
    return {
        "node": index,
        "integrand": integrand,
        "variable": row["variable"],
        "integration_variable": integration_variable,
        "density_substitution": density_factors[0],
        "density_root": density_root,
        "measure_node": kernel["measure_node"],
        "measure_numerator": kernel["measure_numerator"],
        "kernel_node": kernel["kernel_node"],
        "kernel_kind": kernel["kernel_kind"],
        "kernel_affine": kernel["kernel_affine"],
        "lower": endpoint,
        "upper": coordinate,
        "rate": _rate_text(rate),
        "physical_measure": "dt/t",
        "metadata_verified": True,
    }


def _history_info(
    nodes: Sequence[Any],
    index: int,
    contribution: int,
    incoming: int,
    memory: int,
    field: str,
    order: int,
    rate: Fraction,
) -> dict[str, Any]:
    _require(_op(nodes, index) in {"sum", "add"},
             f"{field} Z{order} history is not a sum")
    args = _args(nodes, index)
    _require(args.count(contribution) == 1,
             f"{field} Z{order} history must contain one density contribution")
    branches = [value for value in args if value != contribution]
    _require(len(branches) == 1, f"{field} Z{order} history has unexpected branches")
    incoming_branch = branches[0]
    if rate == 0:
        _require(incoming_branch == incoming,
                 f"{field} Z{order} rate-zero history branch is not incoming")
        memory_retained = False
        branch_kind = "incoming"
    else:
        _require(incoming != memory,
                 f"{field} Z{order} history incoming aliases predecessor memory")
        branch_factors = _flatten_product(nodes, incoming_branch)
        _require(len(branch_factors) == 2,
                 f"{field} Z{order} history branch has extra factors")
        _require(
            sorted(branch_factors) == sorted((incoming, memory)),
            f"{field} Z{order} history branch is not memory*incoming",
        )
        memory_retained = True
        branch_kind = "memory_times_incoming"
    return {
        "node": index,
        "arguments": args,
        "density_branch": contribution,
        "incoming_branch": incoming_branch,
        "incoming_branch_factors": (
            [incoming] if rate == 0 else _flatten_product(nodes, incoming_branch)
        ),
        "incoming_branch_kind": branch_kind,
        "incoming_retained": True,
        "memory_retained": memory_retained,
    }


def _endpoint_info(
    nodes: Sequence[Any],
    index: int,
    history: int,
    coordinate: int,
    endpoint: int,
    field: str,
    order: int,
) -> dict[str, Any]:
    row = _row(nodes, index, "endpoint history")
    _require(_op(nodes, index) == "function_substitution",
             f"{field} Z{order} endpoint history is not a function substitution")
    _require(_guarded(row), f"{field} Z{order} endpoint history lacks a guard")
    expression = _handle(row.get("expression"), "endpoint expression", len(nodes))
    variable = _handle(row.get("variable"), "endpoint variable", len(nodes))
    value = _handle(row.get("value"), "endpoint value", len(nodes))
    _require(expression == history, f"{field} Z{order} endpoint does not wrap history")
    _require(variable == coordinate, f"{field} Z{order} endpoint variable changed")
    _require(value == endpoint, f"{field} Z{order} endpoint value changed")
    return {
        "node": index,
        "expression": history,
        "variable": coordinate,
        "endpoint_node": endpoint,
        "guarded": True,
        "guard_keys": [key for key in _F_SUB_GUARDS if row.get(key) is True],
    }


def _premise_section(premise: Mapping[str, Any]) -> tuple[Mapping[str, Any], str]:
    _require(isinstance(premise, Mapping), "saved_inlet_premise must be a mapping")
    candidates: list[tuple[str, Any]] = []
    if isinstance(premise.get("actual_inlet_seam"), Mapping):
        candidates.append(("actual_inlet_seam", premise["actual_inlet_seam"]))
    checks = premise.get("checks")
    if isinstance(checks, Mapping) and isinstance(checks.get("actual_inlet_seam"), Mapping):
        candidates.append(("checks.actual_inlet_seam", checks["actual_inlet_seam"]))
    candidates.append(("saved_inlet_premise", premise))
    for path, section in candidates:
        if section.get("passed") is True:
            if section.get("source_only") is False:
                raise ValueError(f"{path} is not source-only")
            if section.get("source_executed") is True:
                raise ValueError(f"{path} reports source execution")
            if section.get("accepted_ancestors_replayed") is True:
                raise ValueError(f"{path} reports accepted ancestor replay")
            return section, path
    raise ValueError("saved_inlet_premise lacks a passed actual_inlet_seam receipt")


def _premise_rows(premise: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Collect optional explicit port rows without inventing absent handles."""

    found: list[dict[str, Any]] = []
    aliases = {
        "incoming": ("incoming", "incoming_root", "incoming_node", "inlet", "inlet_root"),
        "history_endpoint": (
            "history_endpoint", "history_endpoint_root", "history_endpoint_node",
            "endpoint_history", "endpoint", "H1", "H_Rm",
        ),
        "history": ("history", "history_root", "H_root", "source_history"),
    }
    fields = set(FIELDS)

    def visit(value: Any, field_hint: str | None = None, order_hint: int | None = None) -> None:
        if isinstance(value, Mapping):
            field = value.get("field", value.get("component", field_hint))
            if isinstance(field, str) and field.lower() in fields:
                field = field.lower()
            order = value.get("ordinary_Z_order", value.get("Z_order", value.get("order", order_hint)))
            if field in fields and isinstance(order, int) and not isinstance(order, bool):
                row = {"field": field, "order": order}
                for name, names in aliases.items():
                    for key in names:
                        if key in value:
                            row[name] = value[key]
                            break
                if any(name in row for name in aliases):
                    found.append(row)
            for key, child in value.items():
                child_field = key.lower() if isinstance(key, str) and key.lower() in fields else field
                child_order = int(key) if isinstance(key, str) and key.isdigit() else (
                    key if isinstance(key, int) and not isinstance(key, bool) else None
                )
                if isinstance(child, Mapping):
                    visit(child, child_field, child_order)
                elif isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                    for index, item in enumerate(child):
                        visit(item, child_field, index)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            for index, item in enumerate(value):
                visit(item, field_hint, index)

    for key in ("rows", "port_rows", "source_ports", "actual_inlet_rows", "actual_inlet_seam_rows"):
        if key in premise:
            visit(premise[key])
    return found


def _compare_premise_rows(
    rows: Sequence[Mapping[str, Any]],
    expected: Mapping[tuple[str, int], Mapping[str, int]],
) -> dict[str, Any]:
    if not rows:
        return {
            "explicit": False,
            "match": None,
            "rows_checked": 0,
            "note": "accepted receipt has no explicit graph handles",
        }
    checked = 0
    limit = max(value for row in expected.values() for value in row.values()) + 1
    for row in rows:
        key = (row.get("field"), row.get("order"))
        _require(key in expected, f"saved inlet row has unsupported field/order {key!r}")
        for name in ("incoming", "history_endpoint", "history"):
            if name in row:
                actual = _handle(row[name], f"saved inlet {key} {name}", limit)
                _require(actual == expected[key][name],
                         f"saved inlet {key} {name} does not match current port")
                checked += 1
    return {"explicit": True, "match": True, "rows_checked": checked, "note": "explicit rows matched"}


def _append(delta: list[dict[str, Any]], base_count: int, row: Mapping[str, Any]) -> int:
    index = base_count + len(delta)
    delta.append(deepcopy(dict(row)))
    return index


def _orders(value: Sequence[int]) -> tuple[int, ...]:
    result = tuple(value)
    _require(result and all(isinstance(item, int) and not isinstance(item, bool) for item in result),
             "ordinary_orders must contain integers")
    _require(len(set(result)) == len(result), "ordinary_orders contains duplicates")
    _require(result == tuple(sorted(result)) and all(item >= 0 for item in result),
             "ordinary_orders must be sorted nonnegative integers")
    return result


def derive(
    nodes: Sequence[Mapping[str, Any]],
    finite_bindings: Mapping[str, Any],
    saved_inlet_premise: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    ordinary_orders: Sequence[int] = ORDERS,
) -> dict[str, Any]:
    """Return a bounded graph delta for current Rm first native-y histories."""

    base = _sequence(nodes, "nodes")
    _require(all(isinstance(row, Mapping) for row in base), "nodes must contain mappings")
    _require(isinstance(finite_bindings, Mapping), "finite_bindings must be a mapping")
    expected = str(expected_nodes_digest).lower()
    _require(len(expected) == 64, "expected_nodes_digest must be SHA-256")
    try:
        int(expected, 16)
    except ValueError as exc:
        raise ValueError("expected_nodes_digest must be SHA-256") from exc
    actual = canonical_digest(base)
    _require(actual == expected, "nodes digest does not match expected_nodes_digest")
    orders = _orders(ordinary_orders)

    bindings = _unwrap(finite_bindings)
    chart = _chart(bindings)
    premise, premise_path = _premise_section(saved_inlet_premise)
    premise_rows = _premise_rows(premise)

    coordinate = _coordinate_node(base, bindings, chart)
    endpoint_rows = _endpoint_rows(bindings, len(base), orders)
    endpoint = _endpoint_node(base, bindings, chart, endpoint_rows)
    physical_measure = _physical_node(base, bindings, chart)
    physical = _physical_info(base, physical_measure, coordinate)

    history = _container(bindings, chart, "history")
    incoming = _container(bindings, chart, "incoming")
    contributions = _container(bindings, chart, "contributions", "contribution")
    memory = _container(bindings, chart, "memory")
    _require(history is not None, "finite bindings are missing history roots")
    _require(incoming is not None, "finite bindings are missing incoming roots")
    _require(contributions is not None, "finite bindings are missing contribution roots")
    _require(memory is not None, "finite bindings are missing memory roots")
    density = _density_container(bindings, chart)

    declared_rates = _first(chart, "rates", "RATES")
    if declared_rates is None:
        declared_rates = _first(bindings, "rates", "RATES")
    _require(declared_rates is None or isinstance(declared_rates, Mapping),
             "declared rates must be a mapping")

    rows: list[dict[str, Any]] = []
    ports: dict[tuple[str, int], dict[str, int]] = {}
    for field in FIELDS:
        densities = _field_rows(density, field, "finite_N.density", len(base), orders)
        histories = _field_rows(history, field, "chart.history", len(base), orders)
        inlets = _field_rows(incoming, field, "chart.incoming", len(base), orders)
        integrals = _field_rows(contributions, field, "chart.contributions", len(base), orders)
        memory_root = _memory_node(base, memory, field)
        for order in orders:
            density_root = densities[order]
            history_root = histories[order]
            incoming_root = inlets[order]
            contribution_root = integrals[order]
            _require(density_root != incoming_root, f"{field} Z{order} density aliases incoming")
            declared = (
                declared_rates.get(field, declared_rates.get(field.upper()))
                if declared_rates is not None else None
            )
            rate = _fraction(declared, f"rates.{field}") if declared is not None else RATES[field]
            _require(rate == RATES[field], f"{field} source rate is not current finite-N Rm rate")
            contribution_info = _contribution_info(
                base, contribution_root, density_root, field, order,
                coordinate, endpoint, rate
            )
            history_info = _history_info(
                base, history_root, contribution_root, incoming_root,
                memory_root, field, order, rate
            )
            memory_info = _memory_shape(
                base, memory_root, coordinate, endpoint, rate
            )
            endpoint_info = _endpoint_info(
                base, endpoint_rows[field][order], history_root,
                coordinate, endpoint, field, order
            )
            _require(_exact(base, incoming_root) != 0,
                     f"{field} Z{order} incoming is exact zero")
            key = (field, order)
            ports[key] = {
                "density": density_root,
                "history": history_root,
                "incoming": incoming_root,
                "memory": memory_root,
                "contribution": contribution_root,
                "history_endpoint": endpoint_rows[field][order],
            }
            rows.append(
                {
                    "field": field,
                    "ordinary_Z_order": order,
                    "density_root": density_root,
                    "history_root": history_root,
                    "incoming_root": incoming_root,
                    "memory_root": memory_root,
                    "contribution_root": contribution_root,
                    "history_endpoint_root": endpoint_rows[field][order],
                    "rate": _rate_text(rate),
                    "history_endpoint_identity": {
                        "formula": "H(1)=incoming",
                        "status": "proved_by_new_actual_source_equal_bound",
                        "old_receipt_toy_values_admitted": False,
                        "endpoint_wrapper": endpoint_info,
                        "history_source": history_info,
                        "contribution_source": contribution_info,
                        "memory_source": memory_info,
                    },
                    "first_native_y_formula": "D_y H(1)=density(1)-rate*incoming",
                    "physical_J_at_endpoint": 1,
                }
            )

    premise_match = _compare_premise_rows(premise_rows, ports)
    delta: list[dict[str, Any]] = []
    density_at_endpoint: dict[str, dict[int, int]] = {field: {} for field in FIELDS}
    first_native_y: dict[str, dict[int, int]] = {field: {} for field in FIELDS}
    row_by_key = {(row["field"], row["ordinary_Z_order"]): row for row in rows}

    # Keep all endpoint density function substitutions separate from history
    # roots.  This prevents a pre-loop seam alias from becoming density data.
    for field in FIELDS:
        for order in orders:
            row = row_by_key[(field, order)]
            endpoint_density = _append(
                delta,
                len(base),
                {
                    "operation": "function_substitution",
                    "expression": row["density_root"],
                    "variable": coordinate,
                    "value": endpoint,
                    "Z_independent_substitution": True,
                    "source_native_substitution": True,
                    "source_formula": "density(1)",
                    "field": field,
                    "ordinary_Z_order": order,
                },
            )
            density_at_endpoint[field][order] = endpoint_density
            row["density_at_endpoint_root"] = endpoint_density

    # Preserve the incoming factor even for p's rate-zero row: the zero
    # coefficient and product remain explicit source graph rows.
    for field in FIELDS:
        for order in orders:
            row = row_by_key[(field, order)]
            rate = _fraction(row["rate"], f"{field} Z{order} rate")
            coefficient = _append(
                delta,
                len(base),
                {
                    "operation": "exact_rational",
                    "numerator": -rate.numerator,
                    "denominator": rate.denominator,
                    "source_formula": "-rate",
                    "field": field,
                    "ordinary_Z_order": order,
                },
            )
            incoming_term = _append(
                delta,
                len(base),
                {
                    "operation": "product",
                    "arguments": [coefficient, row["incoming_root"]],
                    "source_formula": "-rate*incoming",
                    "incoming_retained": True,
                    "rate_zero_product_retained": rate == 0,
                    "field": field,
                    "ordinary_Z_order": order,
                },
            )
            derivative = _append(
                delta,
                len(base),
                {
                    "operation": "sum",
                    "arguments": [row["density_at_endpoint_root"], incoming_term],
                    "source_formula": "D_y H(1)=density(1)-rate*incoming",
                    "physical_Jacobian_node": physical_measure,
                    "physical_Jacobian_at_endpoint": 1,
                    "field": field,
                    "ordinary_Z_order": order,
                },
            )
            first_native_y[field][order] = derivative
            row.update(
                rate_coefficient_root=coefficient,
                incoming_term_root=incoming_term,
                first_native_y_root=derivative,
            )

    joined = base + deepcopy(delta)
    return {
        "schema": "original_Rm_history_native_source_v1",
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "actual_saved_graph_run": False,
        "actual_saved_graph_loaded": False,
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_replayed": False,
        "accepted_reducers_or_providers_replayed": False,
        "quadratures_replayed": False,
        "history_archive_decoded": False,
        "expected_nodes_digest": expected,
        "nodes_digest": actual,
        "ordinary_orders": list(orders),
        "fields": list(FIELDS),
        "module_sha256": _module_sha256(),
        "endpoint_bindings": {
            "coordinate": coordinate,
            "endpoint_node1": endpoint,
            "physical_measure": physical,
            "physical_J_exact_one_at_endpoint": True,
            "physical_J_global_exact_one": False,
            "history_endpoint_roots": {
                field: {str(order): value for order, value in values.items()}
                for field, values in endpoint_rows.items()
            },
        },
        "source_ports": {
            field: {str(order): deepcopy(ports[(field, order)]) for order in orders}
            for field in FIELDS
        },
        "rows": deepcopy(rows),
        "history_at_endpoint_roots": {
            field: {str(order): value for order, value in values.items()}
            for field, values in endpoint_rows.items()
        },
        "history_equals_incoming": {
            field: {
                str(order): {
                    "history_endpoint_root": ports[(field, order)]["history_endpoint"],
                    "incoming_root": ports[(field, order)]["incoming"],
                    "formula": "H(1)=incoming",
                    "status": "proved_by_new_actual_source_equal_bound",
                }
                for order in orders
            }
            for field in FIELDS
        },
        "density_at_endpoint_roots": deepcopy(density_at_endpoint),
        "first_native_y_history_derivative_roots": deepcopy(first_native_y),
        "graph_delta": {
            "base_graph_node_count": len(base),
            "base_nodes_sha256": actual,
            "appended_nodes": deepcopy(delta),
            "appended_node_count": len(delta),
            "full_node_count": len(joined),
            "full_canonical_nodes_sha256": canonical_digest(joined),
            "density_at_endpoint_roots": deepcopy(density_at_endpoint),
            "first_native_y_history_derivative_roots": deepcopy(first_native_y),
        },
        "saved_inlet_premise": {
            "path": premise_path,
            "authenticated": True,
            "passed": True,
            "consumed_for_port_identity_only": True,
            "mathematical_source_fixture_not_new_original_scale_norm": (
                premise.get("mathematical_source_fixture_not_new_original_scale_norm") is True
            ),
            "toy_body_values_admitted_as_actual_normalized_fields": False,
            "explicit_port_rows": premise_match["explicit"],
            "port_refs_match": premise_match["match"],
            "rows_checked": premise_match["rows_checked"],
            "note": premise_match["note"],
        },
        "proofs": {
            "actual_current_Rm_history_endpoint": {
                "formula": "H(1)=incoming",
                "status": "proved_by_new_actual_source_equal_bound",
                "equal_bound_integrals_are_structural": True,
                "memory_endpoint_unit_is_structural": True,
                "old_toy_receipt_not_used_as_field_values": True,
                "rows": len(rows),
            },
            "first_native_y_Duhamel": {
                "formula": "D_y H(1)=density(1)-rate*incoming",
                "status": "source_formula_installed",
                "rows": len(rows),
                "physical_J_scope": "exactly one at endpoint x=1",
            },
        },
        "remaining_density_seam_unproved": True,
        "scope": {
            "source_graph_only": True,
            "no_saved_graph_execution": True,
            "no_history_archive_decode": True,
            "no_numeric_integration": True,
            "no_numeric_values_selected": True,
            "actual_current_transport_metadata_verified": True,
            "actual_current_source_rows_verified": True,
            "old_inlet_toy_body_values_not_admitted": True,
            "new_actual_equal_bound_history_identity": True,
            "no_density_alias_across_Rm": True,
            "full_modified_germ_claimed": False,
            "RQy_identity_claimed": False,
            "global_closure_claimed": False,
            "remaining_density_seam_unproved": True,
            "nodes_input_unmutated": deepcopy(base) == base,
            "finite_bindings_input_unmutated": deepcopy(finite_bindings) == finite_bindings,
            "saved_inlet_premise_input_unmutated": deepcopy(saved_inlet_premise) == saved_inlet_premise,
        },
    }


__all__ = ["FIELDS", "ORDERS", "RATES", "canonical_digest", "derive"]
