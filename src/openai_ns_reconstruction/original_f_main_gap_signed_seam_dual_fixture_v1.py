"""Synthetic Fraction dual-number fixture for the signed main/gap seam.

This module is deliberately independent of the seam leaf.  The caller passes
the already selected ``build`` callable and an ``archives_factory``; this
fixture calls them once, then evaluates only the returned symbolic IR with
synthetic values.  No source handle is dereferenced and no source value is
treated as an archive value.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any, Callable


_Z = Fraction(1, 3)
_DELTA = Fraction(1, 5)
_DELTA_Z = Fraction(2, 7)


class _Dual:
    """Exact value/ordinary-Z derivative pair over ``Fraction``."""

    __slots__ = ("value", "derivative")

    def __init__(self, value: Fraction, derivative: Fraction = Fraction(0)) -> None:
        self.value = Fraction(value)
        self.derivative = Fraction(derivative)

    def __add__(self, other: "_Dual | int | Fraction") -> "_Dual":
        right = _as_dual(other)
        return _Dual(self.value + right.value, self.derivative + right.derivative)

    def __radd__(self, other: "_Dual | int | Fraction") -> "_Dual":
        return self + other

    def __neg__(self) -> "_Dual":
        return _Dual(-self.value, -self.derivative)

    def __sub__(self, other: "_Dual | int | Fraction") -> "_Dual":
        return self + (-_as_dual(other))

    def __rsub__(self, other: "_Dual | int | Fraction") -> "_Dual":
        return _as_dual(other) - self

    def __mul__(self, other: "_Dual | int | Fraction") -> "_Dual":
        right = _as_dual(other)
        return _Dual(
            self.value * right.value,
            self.derivative * right.value + self.value * right.derivative,
        )

    def __rmul__(self, other: "_Dual | int | Fraction") -> "_Dual":
        return self * other

    def __truediv__(self, other: "_Dual | int | Fraction") -> "_Dual":
        right = _as_dual(other)
        if right.value == 0:
            raise ZeroDivisionError("synthetic dual denominator is zero")
        denominator = right.value * right.value
        return _Dual(
            self.value / right.value,
            (self.derivative * right.value - self.value * right.derivative)
            / denominator,
        )

    def __rtruediv__(self, other: "_Dual | int | Fraction") -> "_Dual":
        return _as_dual(other) / self

    def power(self, exponent: int) -> "_Dual":
        if exponent == 0:
            return _Dual(Fraction(1))
        value = self.value**exponent
        return _Dual(value, exponent * (self.value ** (exponent - 1)) * self.derivative)


def _as_dual(value: _Dual | int | Fraction) -> _Dual:
    return value if isinstance(value, _Dual) else _Dual(Fraction(value))


# The tuples are (Z0 value, Z1 value).  Their values are synthetic only; they
# are intentionally not copied from any archive or source graph.
_SYNTHETIC: dict[str, dict[str, tuple[Fraction, Fraction]]] = {
    "main.upper": {
        "R": (Fraction(5, 2), Fraction(3, 7)),
        "U": (Fraction(4, 3), Fraction(-2, 5)),
        "Q": (Fraction(7, 6), Fraction(5, 8)),
        "Q_y": (Fraction(-9, 10), Fraction(11, 13)),
    },
    "gap.lower": {
        "R": (Fraction(8, 3), Fraction(-4, 9)),
        "U": (Fraction(5, 4), Fraction(7, 11)),
        "Q": (Fraction(-6, 7), Fraction(9, 10)),
        "Q_y": (Fraction(13, 12), Fraction(-5, 14)),
    },
}


def _fraction_literal(node: Mapping[str, Any]) -> Fraction:
    text = node.get("text")
    if not isinstance(text, str):
        raise ValueError("exact_literal text is missing")
    return Fraction(text)


def _source_role_parts(role: Any) -> tuple[str, str, str]:
    if not isinstance(role, str):
        raise ValueError("source role is missing")
    parts = role.split(".")
    if len(parts) < 5 or parts[0] != "F":
        raise ValueError("unexpected source_at role")
    field = parts[1]
    side = parts[2] + "." + parts[3]
    suffix = ".".join(parts[4:])
    if side not in _SYNTHETIC or field not in _SYNTHETIC[side]:
        raise ValueError("synthetic source role is outside the fixture")
    return field, side, suffix


def _source_at(node: Mapping[str, Any]) -> _Dual:
    field, side, suffix = _source_role_parts(node.get("source_role"))
    if suffix in ("Z0", "ordinary_Z0", "0"):
        value, derivative = _SYNTHETIC[side][field]
        return _Dual(value, derivative)
    if suffix in ("Z1", "ordinary_Z1", "1"):
        _, value = _SYNTHETIC[side][field]
        # The fixture only needs the value of a supplied first jet when it is
        # consumed by K_Z; no second jet is part of this synthetic contract.
        return _Dual(value)
    raise ValueError("unexpected source_at ordinary-Z suffix")


def _partial_z(node: Mapping[str, Any]) -> _Dual:
    role = node.get("source_role")
    if not isinstance(role, str):
        raise ValueError("partial-Z role is missing")
    if "delta_Z" in role:
        return _Dual(_DELTA_Z)
    # The leaf uses this fallback role when the direct radius shape is not
    # structurally admitted: F.R.<chart>.<side>.partial_Z_at_fixed_query.
    if role.startswith("F.R.") and (
        "partial_Z_at_fixed_query" in role or role.endswith(".R_Z")
    ):
        parts = role.split(".")
        side = parts[2] + "." + parts[3] if len(parts) >= 4 else ""
        if side in _SYNTHETIC:
            return _Dual(_SYNTHETIC[side]["R"][1])
    raise ValueError("unsupported synthetic partial-Z role")


def _eval(node: Any) -> _Dual:
    if not isinstance(node, Mapping):
        raise ValueError("symbolic IR node is not a mapping")
    operator = node.get("operator")
    arguments = node.get("arguments", ())
    if not isinstance(arguments, Sequence) or isinstance(arguments, (str, bytes, bytearray)):
        raise ValueError("symbolic IR arguments are not a sequence")

    if operator == "exact_literal":
        return _Dual(_fraction_literal(node))
    if operator == "formal_source_symbol":
        if node.get("name") != "Z":
            raise ValueError("fixture only admits formal source symbol Z")
        return _Dual(_Z, Fraction(1))
    if operator == "formal_source_atom":
        if node.get("name") != "delta":
            raise ValueError("fixture only admits formal source atom delta")
        return _Dual(_DELTA, _DELTA_Z)
    if operator == "source_at":
        return _source_at(node)
    if operator == "partial_Z_at_fixed_query":
        return _partial_z(node)

    if operator == "sum":
        total = _Dual(Fraction(0))
        for argument in arguments:
            total = total + _eval(argument)
        return total
    if operator == "product":
        total = _Dual(Fraction(1))
        for argument in arguments:
            total = total * _eval(argument)
        return total
    if operator == "negative":
        if len(arguments) != 1:
            raise ValueError("negative expects one argument")
        return -_eval(arguments[0])
    if operator == "quotient":
        if len(arguments) != 2:
            raise ValueError("quotient expects two arguments")
        return _eval(arguments[0]) / _eval(arguments[1])
    if operator == "power":
        if len(arguments) != 2:
            raise ValueError("power expects a base and exponent")
        base = _eval(arguments[0])
        exponent = _eval(arguments[1])
        if exponent.derivative != 0 or exponent.value.denominator != 1:
            raise ValueError("fixture power exponent is not a constant integer")
        return base.power(exponent.value.numerator)
    if operator == "subtract":
        if len(arguments) < 2:
            raise ValueError("subtract expects at least two arguments")
        value = _eval(arguments[0])
        for argument in arguments[1:]:
            value = value - _eval(argument)
        return value
    raise ValueError("unsupported synthetic IR operator")


def _safe_eval(node: Any) -> _Dual | None:
    try:
        return _eval(node)
    except (ArithmeticError, KeyError, TypeError, ValueError):
        return None


def _factory_pair(value: Any) -> tuple[Any, Any] | None:
    if isinstance(value, Mapping) and "F" in value and "own" in value:
        return value["F"], value["own"]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return value[0], value[1]
    return None


def run_fixture(
    build: Callable[[Mapping[str, Any], Mapping[str, Any]], Mapping[str, Any]],
    archives_factory: Callable[[], Any],
) -> dict[str, Any]:
    """Run one synthetic Fraction dual check against one symbolic build."""

    checks: dict[str, bool] = {
        "synthetic_R_Z_nonzero": all(
            values["R"][1] != 0 for values in _SYNTHETIC.values()
        ),
        "synthetic_delta_Z_nonzero": _DELTA_Z != 0,
        "synthetic_jets_nonzero_and_side_distinct": all(
            values[field][order] != 0
            for values in _SYNTHETIC.values()
            for field in ("R", "U", "Q", "Q_y")
            for order in (0, 1)
        )
        and all(
            _SYNTHETIC["main.upper"][field] != _SYNTHETIC["gap.lower"][field]
            for field in ("R", "U", "Q", "Q_y")
        ),
    }
    try:
        pair = _factory_pair(archives_factory())
        checks["archives_factory_pair"] = pair is not None
        if pair is None:
            return {"passed": False, "checks": checks}
        F_archive, own_archive = pair
        built = build(F_archive, own_archive)
        checks["build_result_mapping"] = isinstance(built, Mapping)
        if not isinstance(built, Mapping):
            return {"passed": False, "checks": checks}

        endpoints = built.get("endpoints")
        signed = built.get("signed_seam")
        checks["endpoint_and_signed_seam_present"] = isinstance(endpoints, Mapping) and isinstance(signed, Mapping)
        if not isinstance(endpoints, Mapping) or not isinstance(signed, Mapping):
            return {"passed": False, "checks": checks}

        main = endpoints.get("main.upper")
        gap = endpoints.get("gap.lower")
        checks["main_gap_endpoint_records_present"] = isinstance(main, Mapping) and isinstance(gap, Mapping)
        if not isinstance(main, Mapping) or not isinstance(gap, Mapping):
            return {"passed": False, "checks": checks}

        main_k = _safe_eval(main.get("K"))
        main_kz = _safe_eval(main.get("K_Z"))
        gap_k = _safe_eval(gap.get("K"))
        gap_kz = _safe_eval(gap.get("K_Z"))
        signed_k = _safe_eval(signed.get("K"))
        signed_kz = _safe_eval(signed.get("K_Z"))
        checks["all_synthetic_IR_evaluated"] = all(
            item is not None for item in (main_k, main_kz, gap_k, gap_kz, signed_k, signed_kz)
        )
        if not checks["all_synthetic_IR_evaluated"]:
            return {"passed": False, "checks": checks}

        assert main_k is not None
        assert main_kz is not None
        assert gap_k is not None
        assert gap_kz is not None
        assert signed_k is not None
        assert signed_kz is not None

        checks.update(
            {
                "main_K_Z_matches_dual_derivative": main_k.derivative == main_kz.value,
                "gap_K_Z_matches_dual_derivative": gap_k.derivative == gap_kz.value,
                "signed_K_is_main_minus_gap": signed_k.value == main_k.value - gap_k.value,
                "signed_K_Z_is_main_minus_gap": signed_kz.value == main_kz.value - gap_kz.value,
            }
        )
    except (ArithmeticError, KeyError, TypeError, ValueError):
        checks["fixture_execution_without_source_access"] = False
        return {"passed": False, "checks": checks}

    checks["fixture_execution_without_source_access"] = True
    return {"passed": all(checks.values()), "checks": checks}


__all__ = ["run_fixture"]
