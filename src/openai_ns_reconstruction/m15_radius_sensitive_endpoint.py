"""Caller-owned ordinary-Z1 endpoint algebra with a nonconstant radius jet."""

from __future__ import annotations

from fractions import Fraction
from decimal import Decimal
from numbers import Number, Rational
import math
from typing import Any, Mapping

COMPONENT = "m15_radius_sensitive_endpoint_v1"
SCHEMA = COMPONENT
_NAMES = ("R", "R_Z", "U", "U_Z", "Q", "Q_Z", "Q_y", "Q_yZ", "Z", "delta", "delta_Z")


class EndpointInputError(ValueError):
    """Caller data cannot be evaluated safely."""


def _check_scalar(value: Any, name: str) -> None:
    # This interface is point arithmetic. Reject intervals and symbolic or
    # undecidable inputs; interval enclosures require a separate implementation.
    if getattr(value, "_mpi_", None) is not None:
        raise EndpointInputError(f"{name} interval inputs are unsupported")
    if isinstance(value, bool) or isinstance(value, complex) or not isinstance(value, Number):
        raise EndpointInputError(f"{name} must be a finite real numeric scalar")
    if isinstance(value, (int, Rational)):
        return
    if isinstance(value, Decimal):
        if not value.is_finite():
            raise EndpointInputError(f"{name} must be finite")
        return
    if getattr(value, "_mpf_", None) is not None:
        import mpmath
        if not mpmath.isfinite(value):
            raise EndpointInputError(f"{name} must be finite")
        return
    try:
        if not math.isfinite(value):
            raise EndpointInputError(f"{name} must be finite")
    except (TypeError, ValueError, OverflowError) as exc:
        raise EndpointInputError(f"{name} must be a finite real numeric scalar") from exc


def _zero_status(value: Any) -> tuple[bool, bool]:
    # Comparison is exact for Fraction; for rounded inputs it is only a point
    # check at caller precision and is not a certified separation from zero.
    try:
        return bool(value == 0), isinstance(value, (int, Rational))
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise EndpointInputError("L zero status is undecidable") from exc


def _terms(
    *, R: Any, R_Z: Any, U: Any, U_Z: Any, Q: Any, Q_Z: Any,
    Q_y: Any, Q_yZ: Any, Z: Any, delta: Any, delta_Z: Any,
) -> dict[str, Any]:
    try:
        R2, Q2, U2, Z2 = R * R, Q * Q, U * U, Z * Z
        A, A_Z = R2 * Q, 2 * R * R_Z * Q + R2 * Q_Z
        B, B_Z = R2 * Q2, 2 * R * R_Z * Q2 + 2 * R2 * Q * Q_Z
        C, C_Z = R2 * U2, 2 * R * R_Z * U2 + 2 * R2 * U * U_Z
        D, D_Z = R * Q_y, R_Z * Q_y + R * Q_yZ
        L, L_Z = 1 - delta * Z2, -delta_Z * Z2 - 2 * delta * Z
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise EndpointInputError("caller endpoint values are not one scalar domain") from exc
    contains_zero, exact = _zero_status(L)
    if contains_zero:
        raise EndpointInputError("L=0 at the supplied point")
    try:
        L2, L3 = L * L, L * L * L
        K = A / (2 * L) + B / 4 - Z2 * C / L2 - D
        zterms = {
            "A_Z_over_2L": A_Z / (2 * L),
            "minus_A_L_Z_over_2L_squared": -A * L_Z / (2 * L2),
            "B_Z_over_4": B_Z / 4,
            "minus_2ZC_plus_Z_squared_C_Z_over_L_squared": -(2 * Z * C + Z2 * C_Z) / L2,
            "plus_2Z_squared_C_L_Z_over_L_cubed": 2 * Z2 * C * L_Z / L3,
            "minus_D_Z": -D_Z,
        }
        K_Z = sum(zterms.values(), 0)
        fA_Z, fB_Z, fC_Z, fD_Z = R2 * Q_Z, 2 * R2 * Q * Q_Z, 2 * R2 * U * U_Z, R * Q_yZ
        fixed_terms = {
            "A_Z_over_2L": fA_Z / (2 * L),
            "minus_A_L_Z_over_2L_squared": -A * L_Z / (2 * L2),
            "B_Z_over_4": fB_Z / 4,
            "minus_2ZC_plus_Z_squared_C_Z_over_L_squared": -(2 * Z * C + Z2 * fC_Z) / L2,
            "plus_2Z_squared_C_L_Z_over_L_cubed": 2 * Z2 * C * L_Z / L3,
            "minus_D_Z": -fD_Z,
        }
        fixed_K_Z = sum(fixed_terms.values(), 0)
        omitted = R_Z * (2 * R * (Q / (2 * L) + Q2 / 4 - Z2 * U2 / L2) - Q_y)
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise EndpointInputError("caller endpoint quotient arithmetic failed") from exc
    fixed = {"A": A, "A_Z": fA_Z, "B": B, "B_Z": fB_Z, "C": C, "C_Z": fC_Z,
             "D": D, "D_Z": fD_Z, "K": K, "K_Z": fixed_K_Z, "K_Z_terms": fixed_terms,
             "radius_jet_used": "R_Z=0", "delta_Z_retained": True}
    return {"A": A, "A_Z": A_Z, "B": B, "B_Z": B_Z, "C": C, "C_Z": C_Z,
            "D": D, "D_Z": D_Z, "L": L, "L_Z": L_Z, "K": K, "K_Z": K_Z,
            "K_Z_terms": zterms, "fixed_radius": fixed, "fixed_radius_K": K,
            "fixed_radius_K_Z": fixed_K_Z, "radius_omitted_term": omitted,
            "radius_omission_formula": "R_Z*(2*R*(Q/(2*L)+Q^2/4-Z^2*U^2/L^2)-Q_y)",
            "delta_Z_retained": True, "L_zero_decision_exact": exact,
            "certified": False, "source_equality_proved": False, "source_values_bound": False,
            "graph_traversed": False, "old_functions_called": False,
            "fixed_common_endpoint_assumed": True, "moving_endpoint_chain_rule_supplied_by_caller": "not_verified_caller_responsibility",
            "source_identity_scope": "caller_input_math_only"}


def evaluate_m15_endpoint(
    R: Any, R_Z: Any, U: Any, U_Z: Any, Q: Any, Q_Z: Any,
    Q_y: Any, Q_yZ: Any, Z: Any, delta: Any, delta_Z: Any,
) -> dict[str, Any]:
    """Evaluate K and ordinary fixed-common-endpoint K_Z from all 11 inputs."""
    values = dict(zip(_NAMES, (R, R_Z, U, U_Z, Q, Q_Z, Q_y, Q_yZ, Z, delta, delta_Z)))
    for name in _NAMES:
        _check_scalar(values[name], name)
    if R <= 0:
        raise EndpointInputError("physical radius R must be positive")
    result = _terms(**values)
    result["schema"], result["status"] = SCHEMA, "numerical_returned"
    result["interval_inputs_present"] = False
    result["inputs"] = values
    return result


ordinary_z1_endpoint = evaluate_m15_endpoint
evaluate_endpoint = evaluate_m15_endpoint
compute_m15_endpoint = evaluate_m15_endpoint


class _Polynomial:
    """Exact Fraction polynomial used only by the source-free fixture."""

    def __init__(self, coefficients: tuple[Fraction, ...] | list[Fraction]) -> None:
        values = tuple(Fraction(value) for value in coefficients)
        while len(values) > 1 and values[-1] == 0:
            values = values[:-1]
        self.coefficients = values or (Fraction(0),)

    def __add__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        right = _poly(other)
        width = max(len(self.coefficients), len(right.coefficients))
        return type(self)([
            (self.coefficients[i] if i < len(self.coefficients) else 0)
            + (right.coefficients[i] if i < len(right.coefficients) else 0)
            for i in range(width)
        ])

    def __radd__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        return self + other

    def __neg__(self) -> "_Polynomial":
        return type(self)([-value for value in self.coefficients])

    def __sub__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        return self + (-_poly(other))

    def __rsub__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        return _poly(other) - self

    def __mul__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        right = _poly(other)
        values = [Fraction(0)] * (len(self.coefficients) + len(right.coefficients) - 1)
        for i, left in enumerate(self.coefficients):
            for j, value in enumerate(right.coefficients):
                values[i + j] += left * value
        return type(self)(values)

    def __rmul__(self, other: "_Polynomial | int | Fraction") -> "_Polynomial":
        return self * other

    def derivative(self) -> "_Polynomial":
        return type(self)([
            i * value for i, value in enumerate(self.coefficients[1:], 1)
        ] or [Fraction(0)])

    def evaluate(self, value: Fraction) -> Fraction:
        result = Fraction(0)
        for coefficient in reversed(self.coefficients):
            result = result * value + coefficient
        return result


def _poly(value: _Polynomial | int | Fraction) -> _Polynomial:
    return value if isinstance(value, _Polynomial) else _Polynomial((Fraction(value),))


def _fixture_polynomials() -> dict[str, _Polynomial]:
    return {
        "R": _Polynomial((2, Fraction(1, 3), Fraction(-1, 5))),
        "U": _Polynomial((1, Fraction(-2, 7), Fraction(1, 6))),
        "Q": _Polynomial((Fraction(3, 5), Fraction(-1, 4), Fraction(1, 9))),
        "Q_y": _Polynomial((Fraction(-1, 6), Fraction(1, 8), Fraction(1, 10))),
    }


def _reference(
    p: Mapping[str, _Polynomial], *, z: Fraction, delta: Fraction, delta_z: Fraction
) -> dict[str, Fraction]:
    """Independent hand expansion of every displayed component."""
    r, rz = p["R"].evaluate(z), p["R"].derivative().evaluate(z)
    u, uz = p["U"].evaluate(z), p["U"].derivative().evaluate(z)
    q, qz = p["Q"].evaluate(z), p["Q"].derivative().evaluate(z)
    qy, qyz = p["Q_y"].evaluate(z), p["Q_y"].derivative().evaluate(z)
    r2, q2, u2, z2 = r * r, q * q, u * u, z * z
    A, Az = r2 * q, 2 * r * rz * q + r2 * qz
    B, Bz = r2 * q2, 2 * r * rz * q2 + 2 * r2 * q * qz
    C, Cz = r2 * u2, 2 * r * rz * u2 + 2 * r2 * u * uz
    D, Dz = r * qy, rz * qy + r * qyz
    L, Lz = 1 - delta * z2, -delta_z * z2 - 2 * delta * z
    K = A / (2 * L) + B / 4 - z2 * C / (L * L) - D
    terms = {
        "A_Z_over_2L": Az / (2 * L),
        "minus_A_L_Z_over_2L_squared": -A * Lz / (2 * L * L),
        "B_Z_over_4": Bz / 4,
        "minus_2ZC_plus_Z_squared_C_Z_over_L_squared": -(2 * z * C + z2 * Cz) / (L * L),
        "plus_2Z_squared_C_L_Z_over_L_cubed": 2 * z2 * C * Lz / (L * L * L),
        "minus_D_Z": -Dz,
    }
    return {"R": r, "R_Z": rz, "U": u, "U_Z": uz, "Q": q, "Q_Z": qz,
            "Q_y": qy, "Q_yZ": qyz, "Z": z, "delta": delta, "delta_Z": delta_z,
            "A": A, "A_Z": Az, "B": B, "B_Z": Bz, "C": C, "C_Z": Cz,
            "D": D, "D_Z": Dz, "L": L, "L_Z": Lz, "K": K,
            "K_Z": sum(terms.values(), Fraction(0)), **{
                "term_" + key: value for key, value in terms.items()
            }}


def _k_from_polynomials(p: Mapping[str, _Polynomial], z: Fraction, delta: Fraction) -> Fraction:
    r, u, q, qy = (p[name].evaluate(z) for name in ("R", "U", "Q", "Q_y"))
    z2, L = z * z, 1 - delta * z * z
    if L == 0:
        raise AssertionError("finite-difference sample is singular")
    return (r * r * q) / (2 * L) + (r * r * q * q) / 4 - z2 * (r * r * u * u) / (L * L) - r * qy


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run one exact Fraction fixture; it binds no physical/source values."""
    p = _fixture_polynomials()
    z, delta, delta_z = Fraction(1, 5), Fraction(1, 7), Fraction(2, 11)
    ref = _reference(p, z=z, delta=delta, delta_z=delta_z)
    result = evaluate_m15_endpoint(
        ref["R"], ref["R_Z"], ref["U"], ref["U_Z"], ref["Q"], ref["Q_Z"],
        ref["Q_y"], ref["Q_yZ"], z, delta, delta_z,
    )
    expected_omitted = ref["R_Z"] * (
        2 * ref["R"] * (ref["Q"] / (2 * ref["L"]) + ref["Q"] ** 2 / 4
        - z * z * ref["U"] ** 2 / (ref["L"] ** 2)) - ref["Q_y"]
    )
    h = Fraction(1, 10_000)
    delta_at = lambda x: delta + delta_z * (x - z)
    finite_difference = (
        _k_from_polynomials(p, z - 2 * h, delta_at(z - 2 * h))
        - 8 * _k_from_polynomials(p, z - h, delta_at(z - h))
        + 8 * _k_from_polynomials(p, z + h, delta_at(z + h))
        - _k_from_polynomials(p, z + 2 * h, delta_at(z + 2 * h))
    ) / (12 * h)
    rz0 = dict(ref, R_Z=Fraction(0))
    rz0_result = evaluate_m15_endpoint(
        rz0["R"], rz0["R_Z"], rz0["U"], rz0["U_Z"], rz0["Q"], rz0["Q_Z"],
        rz0["Q_y"], rz0["Q_yZ"], z, delta, delta_z,
    )
    dz0_result = evaluate_m15_endpoint(
        ref["R"], ref["R_Z"], ref["U"], ref["U_Z"], ref["Q"], ref["Q_Z"],
        ref["Q_y"], ref["Q_yZ"], z, delta, Fraction(0),
    )
    try:
        evaluate_m15_endpoint(ref["R"], ref["R_Z"], ref["U"], ref["U_Z"], ref["Q"], ref["Q_Z"],
                              ref["Q_y"], ref["Q_yZ"], Fraction(1), Fraction(1), delta_z)
    except EndpointInputError:
        singular_rejected = True
    else:
        singular_rejected = False
    def rejects(**overrides):
        params = {name: ref[name] for name in _NAMES}
        params.update(overrides)
        try:
            evaluate_m15_endpoint(**params)
        except EndpointInputError:
            return True
        return False
    terms = tuple(result["K_Z_terms"].values())
    controls = {
        "nonfinite_input_rejected": rejects(Q_Z=float("nan")),
        "nonpositive_radius_rejected": rejects(R=Fraction(0)),
        "fraction_exact_hand_reference": result["K"] == ref["K"] and result["K_Z"] == ref["K_Z"],
        "all_radius_leibniz_components": all(result[key] == ref[key] for key in
            ("A_Z", "B_Z", "C_Z", "D_Z")),
        "finite_difference_polynomial_expression": abs(finite_difference - ref["K_Z"]) < Fraction(1, 10**10),
        "nonzero_R_Z_omission": ref["R_Z"] != 0 and expected_omitted != 0,
        "fixed_radius_omission_identity": result["K_Z"] - result["fixed_radius_K_Z"] == expected_omitted,
        "R_Z_zero_special_case": rz0_result["radius_omitted_term"] == 0 and rz0_result["K_Z"] == rz0_result["fixed_radius_K_Z"],
        "delta_Z_nonzero_retained": result["L_Z"] != dz0_result["L_Z"] and result["K_Z"] != dz0_result["K_Z"],
        "positive_and_negative_terms": any(value > 0 for value in terms) and any(value < 0 for value in terms),
        "L_zero_rejected": singular_rejected,
        "caller_math_not_certified": result["certified"] is False and result["source_equality_proved"] is False,
    }
    return {
        "schema": "m15_radius_sensitive_endpoint_synthetic_fixture_v1", "synthetic_only": True,
        "source_values_selected": False, "source_equality_proved": False,
        "accepted_as_source_binding": False, "fixed_common_endpoint_assumed": True,
        "moving_endpoint_chain_rule_supplied_by_caller": "not_verified_caller_responsibility",
        "Z": str(z), "delta": str(delta), "delta_Z": str(delta_z),
        "K": str(result["K"]), "K_Z": str(result["K_Z"]),
        "independent_K": str(ref["K"]), "independent_K_Z": str(ref["K_Z"]),
        "finite_difference_K_Z": str(finite_difference),
        "finite_difference_error": str(abs(finite_difference - ref["K_Z"])),
        "radius_omitted_term": str(expected_omitted), "fixed_radius_K_Z": str(result["fixed_radius_K_Z"]),
        "K_Z_terms": {key: str(value) for key, value in result["K_Z_terms"].items()},
        "controls": controls, "passed": all(controls.values()), "runtime_model_metadata": None,
    }


run_fixture = verify_synthetic_fixture


__all__ = ["COMPONENT", "SCHEMA", "EndpointInputError", "evaluate_m15_endpoint",
           "ordinary_z1_endpoint", "evaluate_endpoint", "compute_m15_endpoint",
           "verify_synthetic_fixture", "run_fixture"]


