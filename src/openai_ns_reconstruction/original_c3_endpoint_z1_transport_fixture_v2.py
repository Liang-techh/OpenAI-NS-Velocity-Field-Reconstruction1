"""Independent synthetic dual fixture for the C3 endpoint Z1 transport.

The fixture receives a ``build(refs)`` callable and never imports the
transport producer.  It supplies JSON-shaped opaque refs, evaluates only the
returned expression nodes with a local exact ``Fraction`` dual number, and
derives the comparison target from fresh analytic functions
``R(Z), delta(Z), Q(Z), U(Z), Q_y(Z)``.

The synthetic field rows deliberately model the old fixed-geometry records:
``A1=R^2 Q1``, ``B1=R^2 (Q^2)_1``, ``C1=R^2 (U^2)_1``, and
``D1=R Q_y1``.  Radius and ``delta`` derivatives are then checked through
the two explicit correction terms.  This is an algebra check only: no source
graph, archive, gzip member, source value, cancellation, or certification is
used.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "original_c3_endpoint_z1_transport_fixture_v2"


class _Dual:
    """Exact value/ordinary-Z1 pair over rational numbers."""

    __slots__ = ("value", "derivative")

    def __init__(self, value: Fraction | int, derivative: Fraction | int = 0) -> None:
        self.value = Fraction(value)
        self.derivative = Fraction(derivative)

    def __add__(self, other: _Dual | Fraction | int) -> _Dual:
        right = _as_dual(other)
        return _Dual(self.value + right.value, self.derivative + right.derivative)

    def __radd__(self, other: _Dual | Fraction | int) -> _Dual:
        return self + other

    def __neg__(self) -> _Dual:
        return _Dual(-self.value, -self.derivative)

    def __sub__(self, other: _Dual | Fraction | int) -> _Dual:
        return self + (-_as_dual(other))

    def __rsub__(self, other: _Dual | Fraction | int) -> _Dual:
        return _as_dual(other) - self

    def __mul__(self, other: _Dual | Fraction | int) -> _Dual:
        right = _as_dual(other)
        return _Dual(
            self.value * right.value,
            self.derivative * right.value + self.value * right.derivative,
        )

    def __rmul__(self, other: _Dual | Fraction | int) -> _Dual:
        return self * other

    def __truediv__(self, other: _Dual | Fraction | int) -> _Dual:
        right = _as_dual(other)
        if right.value == 0:
            raise ZeroDivisionError("synthetic dual denominator is zero")
        return _Dual(
            self.value / right.value,
            (self.derivative * right.value - self.value * right.derivative)
            / (right.value * right.value),
        )

    def __rtruediv__(self, other: _Dual | Fraction | int) -> _Dual:
        return _as_dual(other) / self

    def power(self, exponent: int) -> _Dual:
        if exponent == 0:
            return _Dual(1)
        return _Dual(
            self.value**exponent,
            exponent * self.value ** (exponent - 1) * self.derivative,
        )


def _as_dual(value: _Dual | Fraction | int) -> _Dual:
    return value if isinstance(value, _Dual) else _Dual(value)


def _analytic_family(z: _Dual, moving_geometry: bool = True) -> dict[str, _Dual]:
    """Return fresh nondegenerate analytic synthetic functions of ``Z``."""

    if moving_geometry:
        radius = 3 + z + z.power(2) / 5
        delta = Fraction(1, 5) + z / 7 + z.power(2) / 11
    else:
        # Keep the same values while setting only the geometric derivatives to
        # zero.  The field functions still have ordinary Z1 rows.
        moving = _analytic_family(z, moving_geometry=True)
        radius = _Dual(moving["R"].value)
        delta = _Dual(moving["delta"].value)

    return {
        "R": radius,
        "delta": delta,
        "Q": 4 + 2 * z + z.power(2) / 3,
        "U": 1 + 3 * z - z.power(2) / 4,
        "Q_y": -2 + z / 2 + z.power(2) / 6,
    }


def _products(family: Mapping[str, _Dual]) -> dict[str, _Dual]:
    """Build A, B, C, D and the old fixed-geometry first rows."""

    R, Q, U, Q_y = family["R"], family["Q"], family["U"], family["Q_y"]
    R0 = R.value
    return {
        "A0": R * R * Q,
        "A1": _Dual(R0 * R0 * Q.derivative),
        "B0": R * R * Q * Q,
        "B1": _Dual(R0 * R0 * 2 * Q.value * Q.derivative),
        "C0": R * R * U * U,
        "C1": _Dual(R0 * R0 * 2 * U.value * U.derivative),
        "D0": R * Q_y,
        "D1": _Dual(R0 * Q_y.derivative),
    }


def _direct_K(products: Mapping[str, _Dual], z: _Dual, delta: _Dual) -> _Dual:
    """Independent full K expression; intentionally not copied from the DAG."""

    L = 1 - delta * z.power(2)
    return products["A0"] / (2 * L) + products["B0"] / 4 - z.power(2) * products["C0"] / L.power(2) - products["D0"]


def _opaque_refs(tag: str) -> dict[str, dict[str, str]]:
    names = (
        "A0", "A1", "B0", "B1", "C0", "C1", "D0", "D1",
        "Z", "delta", "rho_Z", "delta_Z",
    )
    return {name: {"fixture": tag, "role": name} for name in names}


def _literal(node: Mapping[str, Any]) -> Fraction:
    text = node.get("text")
    if not isinstance(text, str):
        raise ValueError("exact_literal text is missing")
    return Fraction(text)


def _evaluate(node: Any, values: Mapping[str, _Dual]) -> _Dual:
    """Evaluate the small allow-listed DAG language used by the producer."""

    if not isinstance(node, Mapping):
        raise ValueError("DAG node must be a mapping")
    operator = node.get("operator")
    arguments = node.get("arguments", ())
    if not isinstance(arguments, Sequence) or isinstance(arguments, (str, bytes, bytearray)):
        raise ValueError("DAG arguments must be a sequence")

    if operator == "exact_literal":
        return _Dual(_literal(node))
    if operator == "opaque_source_ref":
        role = node.get("source_role")
        if not isinstance(role, str) or role not in values:
            raise ValueError("opaque source role is not admitted")
        if len(arguments) != 1:
            raise ValueError("opaque source ref must have one metadata argument")
        # Deliberately do not inspect arguments[0].  It is an opaque JSON ref,
        # never a numeric value for the evaluator.
        return values[role]
    if operator == "formal_source_symbol":
        role = node.get("name")
        if not isinstance(role, str) or role not in values:
            raise ValueError("formal source symbol is not admitted")
        return values[role]

    if operator in ("sum", "add"):
        total = _Dual(0)
        for argument in arguments:
            total = total + _evaluate(argument, values)
        return total
    if operator in ("product", "multiply"):
        total = _Dual(1)
        for argument in arguments:
            total = total * _evaluate(argument, values)
        return total
    if operator in ("negative", "negate"):
        if len(arguments) != 1:
            raise ValueError("negative expects one argument")
        return -_evaluate(arguments[0], values)
    if operator in ("quotient", "divide"):
        if len(arguments) != 2:
            raise ValueError("quotient expects two arguments")
        return _evaluate(arguments[0], values) / _evaluate(arguments[1], values)
    if operator == "subtract":
        if len(arguments) < 2:
            raise ValueError("subtract expects at least two arguments")
        result = _evaluate(arguments[0], values)
        for argument in arguments[1:]:
            result = result - _evaluate(argument, values)
        return result
    if operator == "power":
        if len(arguments) != 2:
            raise ValueError("power expects base and exponent")
        base = _evaluate(arguments[0], values)
        exponent = _evaluate(arguments[1], values)
        if exponent.derivative != 0 or exponent.value.denominator != 1:
            raise ValueError("power exponent must be a constant integer")
        return base.power(exponent.value.numerator)
    raise ValueError(f"unsupported DAG operator: {operator!r}")


def _safe_evaluate(node: Any, values: Mapping[str, _Dual]) -> _Dual | None:
    try:
        return _evaluate(node, values)
    except (ArithmeticError, KeyError, TypeError, ValueError):
        return None


def _value_table(family: Mapping[str, _Dual], products: Mapping[str, _Dual], z: _Dual) -> dict[str, _Dual]:
    radius = family["R"]
    delta = family["delta"]
    return {
        **products,
        "Z": z,
        "delta": delta,
        "rho_Z": _Dual(radius.derivative / radius.value),
        "delta_Z": _Dual(delta.derivative),
    }


def _check_zero_geometry(build: Callable[[Mapping[str, Any]], Mapping[str, Any]]) -> bool:
    """With R_Z=delta_Z=0, the full row must regress to K1_fixed."""

    z = _Dual(Fraction(1, 3), 1)
    family = _analytic_family(z, moving_geometry=False)
    products = _products(family)
    refs = _opaque_refs("zero-geometry")
    built = build(refs)
    values = _value_table(family, products, z)
    # Explicitly replace the geometry derivatives with zero while retaining
    # all nonzero field rows.
    values["rho_Z"] = _Dual(0)
    values["delta_Z"] = _Dual(0)
    fixed = _safe_evaluate(built.get("K1_fixed"), values)
    full = _safe_evaluate(built.get("K1_full"), values)
    direct = _direct_K(products, z, family["delta"])
    return fixed is not None and full is not None and fixed.value == full.value == direct.derivative


def run_fixture(build: Callable[[Mapping[str, Any]], Mapping[str, Any]]) -> dict[str, Any]:
    """Run the independent exact synthetic transport checks."""

    checks: dict[str, bool] = {}
    try:
        z = _Dual(Fraction(1, 3), 1)
        family = _analytic_family(z, moving_geometry=True)
        products = _products(family)
        refs = _opaque_refs("moving-geometry")
        built = build(refs)
        checks["build_result_mapping"] = isinstance(built, Mapping)
        if not checks["build_result_mapping"]:
            return {
                "schema": SCHEMA,
                "synthetic_only": True,
                "passed": False,
                "checks": checks,
                "source_values_computed": False,
                "numeric_values_computed": False,
            }

        checks["required_DAG_nodes_present"] = all(
            key in built for key in ("inputs", "K", "K1_fixed", "K1_full", "intermediates", "scope")
        )
        checks["opaque_refs_retained"] = (
            isinstance(built.get("inputs"), Mapping)
            and all(built["inputs"].get(name, {}).get("arguments", [None])[0] == refs[name]
                    for name in refs)
        )
        try:
            json.dumps(built, allow_nan=False)
            checks["DAG_JSON_safe"] = True
        except (TypeError, ValueError):
            checks["DAG_JSON_safe"] = False

        values = _value_table(family, products, z)
        direct = _direct_K(products, z, family["delta"])
        fixed_radius = _Dual(family["R"].value)
        fixed_delta = _Dual(family["delta"].value)
        fixed_products = dict(products)
        fixed_products["A0"] = fixed_radius * fixed_radius * family["Q"]
        fixed_products["B0"] = fixed_radius * fixed_radius * family["Q"] * family["Q"]
        fixed_products["C0"] = fixed_radius * fixed_radius * family["U"] * family["U"]
        fixed_products["D0"] = fixed_radius * family["Q_y"]
        fixed_direct = _direct_K(fixed_products, z, fixed_delta)

        K = _safe_evaluate(built.get("K"), values)
        K1_fixed = _safe_evaluate(built.get("K1_fixed"), values)
        K1_full = _safe_evaluate(built.get("K1_full"), values)
        rho_correction = _safe_evaluate(
            built.get("intermediates", {}).get("rho_correction"), values
        )
        delta_correction = _safe_evaluate(
            built.get("intermediates", {}).get("delta_correction"), values
        )
        checks["all_synthetic_DAG_nodes_evaluated"] = all(
            item is not None for item in (K, K1_fixed, K1_full, rho_correction, delta_correction)
        )
        checks["synthetic_R_Z_nonzero"] = family["R"].derivative != 0
        checks["synthetic_delta_Z_nonzero"] = family["delta"].derivative != 0
        checks["synthetic_field_rows_nonzero"] = all(
            family[name].value != 0 and family[name].derivative != 0
            for name in ("Q", "U", "Q_y")
        )
        checks["positive_R_and_nonzero_L"] = (
            family["R"].value > 0
            and (1 - family["delta"].value * z.value * z.value) != 0
        )
        if all(checks[name] for name in (
            "all_synthetic_DAG_nodes_evaluated",
            "synthetic_R_Z_nonzero",
            "synthetic_delta_Z_nonzero",
            "synthetic_field_rows_nonzero",
        )):
            assert K is not None
            assert K1_fixed is not None
            assert K1_full is not None
            assert rho_correction is not None
            assert delta_correction is not None
            checks["K_value_matches_full_dual"] = K.value == direct.value
            checks["fixed_candidate_matches_fixed_geometry_dual"] = (
                K1_fixed.value == fixed_direct.derivative
            )
            checks["corrected_K1_full_matches_moving_dual"] = (
                K1_full.value == direct.derivative
            )
            checks["rho_delta_correction_matches_difference"] = (
                rho_correction.value + delta_correction.value
                == direct.derivative - fixed_direct.derivative
            )
        else:
            checks["K_value_matches_full_dual"] = False
            checks["fixed_candidate_matches_fixed_geometry_dual"] = False
            checks["corrected_K1_full_matches_moving_dual"] = False
            checks["rho_delta_correction_matches_difference"] = False

        checks["zero_geometry_derivatives_recover_old_candidate"] = _check_zero_geometry(build)
    except (ArithmeticError, AssertionError, KeyError, TypeError, ValueError):
        checks["fixture_execution_without_source_access"] = False
        return {
            "schema": SCHEMA,
            "synthetic_only": True,
            "passed": False,
            "checks": checks,
            "source_values_computed": False,
            "numeric_values_computed": False,
            "source_graph_accessed": False,
            "old_leaf_called": False,
        }

    checks["fixture_execution_without_source_access"] = True
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "schema": SCHEMA,
        "synthetic_only": True,
        "passed": bool(all(checks.values())),
        "checks": checks,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "source_graph_accessed": False,
        "source_graph_traversed": False,
        "old_leaf_called": False,
    }


run = run_fixture


__all__ = ["SCHEMA", "run", "run_fixture"]
