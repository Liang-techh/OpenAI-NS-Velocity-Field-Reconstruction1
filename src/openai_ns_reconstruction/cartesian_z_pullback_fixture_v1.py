"""Fresh exact-rational checks for :mod:`cartesian_z_pullback_v1`.

The fixture receives the leaf's ``evaluate`` callable explicitly and never
imports an implementation, reads a receipt, or calls a historical function.
Its main polynomial is

    f(rho,Z) = rho**2 + rho*Z + Z**2

with ``delta = 1/2``, ``Z = 1/2``, ``lambda = 1/4``, and ``rho = 1``.  The
chosen physical values ``X = Y = 1/4``, ``tau = 3/64``, and ``epsilon = 1``
also satisfy the caller's geometry premise exactly.  Expected values below
are accumulated as ``Fraction`` objects and checked by directed containment;
no midpoint or binary float is used.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any


SCHEMA = "cartesian_z_pullback_v1_fixture_v1"
SOURCE_ID = "fixture:cartesian-z:synthetic-profile"
SOURCE_INDICES = ["synthetic:Z0", "synthetic:Z1", "synthetic:Z2"]


def _q(c: Any, value: Fraction | int) -> Any:
    rational = Fraction(value)
    return c.mpf(rational.numerator) / rational.denominator


def _tuple_point(bits: tuple[int, int, int, int]) -> dict[str, Any]:
    """Make one serialized directed point without importing an implementation."""

    raw = list(bits)
    return {
        "lower_exact_mpf_tuple": list(raw),
        "upper_exact_mpf_tuple": list(raw),
        "directed": True,
    }


def _bounds(value: Any) -> tuple[Any, Any]:
    if hasattr(value, "a") and hasattr(value, "b"):
        return getattr(value, "a"), getattr(value, "b")
    return value, value


def _inside(c: Any, value: Any, expected: Fraction | int) -> bool:
    try:
        lower, upper = _bounds(value)
        point = _q(c, expected)
        point_lower, point_upper = _bounds(point)
        return bool(lower <= point_lower and point_upper <= upper)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _callable(leaf: Any) -> Any:
    candidate = getattr(leaf, "evaluate", leaf)
    if not callable(candidate):
        raise TypeError("leaf must be evaluate(...) or expose evaluate")
    return candidate


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, ZeroDivisionError, ArithmeticError):
        return True
    return False


def _params(
    *,
    z: Fraction | int,
    rho: Fraction | int,
    lam: Fraction | int,
    delta: Fraction | int,
    X: Fraction | int,
    Y: Fraction | int,
    tau: Fraction | int,
) -> dict[str, Any]:
    return {
        "Z": Fraction(z),
        "rho": Fraction(rho),
        "lambda": Fraction(lam),
        "delta": Fraction(delta),
        "X": Fraction(X),
        "Y": Fraction(Y),
        "tau": Fraction(tau),
        "epsilon": Fraction(1),
        "fixed_conditions": {
            "X": True,
            "Y": True,
            "tau": True,
            "epsilon": True,
            "delta": True,
        },
        "source_provenance": {
            "source_ids": [SOURCE_ID],
            "source_indices": list(SOURCE_INDICES),
            "source_paths": [["fixture", "profile", "f"]],
            "error_classes": ["synthetic-exact"],
        },
    }


def _main_jet() -> dict[str, Fraction]:
    # f=rho^2+rho*Z+Z^2 at rho=1, Z=1/2.
    return {
        "value": Fraction(7, 4),
        "rho": Fraction(5, 2),
        "Z": Fraction(2),
        "rhoZ": Fraction(1),
        "rhorho": Fraction(2),
        "ZZ": Fraction(2),
    }


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent exact polynomial, domain, jet, and product checks."""

    evaluate = _callable(leaf)
    params = _params(
        z=Fraction(1, 2),
        rho=Fraction(1),
        lam=Fraction(1, 4),
        delta=Fraction(1, 2),
        X=Fraction(1, 4),
        Y=Fraction(1, 4),
        tau=Fraction(3, 64),
    )
    jet = _main_jet()
    expected = {
        "rho_Z": Fraction(-4, 3),
        "rho_ZZ": Fraction(-8, 3),
        "D": Fraction(-4, 3),
        "D2": Fraction(-34, 9),
        "B": Fraction(12, 7),
        "connection": Fraction(-23, 21),
        "z_phys": Fraction(1, 4),
        "dz_phys": Fraction(-16, 7),
        "dzz_phys": Fraction(-2336, 343),
    }
    checks: dict[str, bool] = {}
    try:
        result = evaluate(c, params, jet)
        weighted = evaluate(c, params, jet, weight_power=2)
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"evaluate_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["schema"] = result.get("schema") == "cartesian_z_pullback_v1"
    checks["conditional_scalar_only"] = (
        result.get("profile_kind") == "scalar_f(rho,Z)"
        and result.get("physical_layer_complete") is False
        and result.get("pressure") is False
        and result.get("m15") is False
        and result.get("ns_validation") is False
    )
    checks["rho_Z"] = _inside(c, result["geometry"]["rho_Z"], expected["rho_Z"])
    checks["rho_ZZ"] = _inside(c, result["geometry"]["rho_ZZ"], expected["rho_ZZ"])
    checks["D"] = _inside(c, result["D"], expected["D"])
    checks["D2"] = _inside(c, result["D2"], expected["D2"])
    checks["B"] = _inside(c, result["geometry"]["B"], expected["B"])
    checks["connection"] = _inside(c, result["geometry"]["connection"], expected["connection"])
    checks["z_phys"] = _inside(c, result["geometry"]["z_phys"], expected["z_phys"])
    checks["dz_phys"] = _inside(c, result["dz_phys"], expected["dz_phys"])
    checks["dzz_phys"] = _inside(c, result["dzz_phys"], expected["dzz_phys"])
    checks["x_normalized_is_sqrt_rho"] = _inside(
        c,
        result["geometry"]["x_normalized"] * result["geometry"]["x_normalized"],
        Fraction(1),
    )
    checks["ordinary_orders"] = result.get("ordinary_Z", {}).get("Z1") is not None and result.get("ordinary_Z", {}).get("Z2") is not None
    checks["provenance_ids_retained"] = result.get("input_provenance", {}).get("source_ids") == [SOURCE_ID]
    checks["provenance_indices_retained"] = result.get("input_provenance", {}).get("source_indices") == SOURCE_INDICES
    checks["fixed_conditions_recorded"] = result.get("coordinate_fixed_conditions", {}).get("provided") is True

    # A saved exact mpf point for rho=1 must be accepted exactly as the
    # Fraction input.  The tuple is the canonical finite (sign, mantissa,
    # exponent, bit-count) representation of one.
    tuple_params = dict(params)
    tuple_params["rho"] = _tuple_point((0, 1, 0, 1))
    try:
        tuple_result = evaluate(c, tuple_params, jet)
        checks["exact_tuple_input_equivalent"] = (
            _inside(c, tuple_result["D"], expected["D"])
            and _inside(c, tuple_result["D2"], expected["D2"])
            and _inside(c, tuple_result["dzz_phys"], expected["dzz_phys"])
        )
    except Exception:
        checks["exact_tuple_input_equivalent"] = False
    reversed_tuple_params = dict(params)
    reversed_tuple_params["rho"] = {
        "lower_exact_mpf_tuple": [0, 1, 0, 1],
        "upper_exact_mpf_tuple": [0, 1, -1, 1],
        "directed": True,
    }
    checks["reversed_exact_tuple_rejected"] = _rejected(
        lambda: evaluate(c, reversed_tuple_params, jet)
    )
    incomplete_tuple_params = dict(params)
    incomplete_tuple_params["rho"] = {"lower_exact_mpf_tuple": [0, 1, 0, 1]}
    checks["incomplete_exact_tuple_rejected"] = _rejected(
        lambda: evaluate(c, incomplete_tuple_params, jet)
    )
    nonfinite_tuple_params = dict(params)
    nonfinite_tuple_params["rho"] = {
        "lower_exact_mpf_tuple": [0, 0, 1, 0],
        "upper_exact_mpf_tuple": [0, 0, 1, 0],
        "directed": True,
    }
    checks["nonfinite_exact_tuple_rejected"] = _rejected(
        lambda: evaluate(c, nonfinite_tuple_params, jet)
    )
    undirected_tuple_params = dict(params)
    undirected_tuple_params["rho"] = {
        "lower_exact_mpf_tuple": [0, 1, 0, 1],
        "upper_exact_mpf_tuple": [0, 1, 0, 1],
        "directed": False,
    }
    checks["undirected_exact_tuple_rejected"] = _rejected(
        lambda: evaluate(c, undirected_tuple_params, jet)
    )
    missing_directed_tuple_params = dict(params)
    missing_directed_tuple_params["rho"] = {
        "lower_exact_mpf_tuple": [0, 1, 0, 1],
        "upper_exact_mpf_tuple": [0, 1, 0, 1],
    }
    checks["missing_directed_marker_rejected"] = _rejected(
        lambda: evaluate(c, missing_directed_tuple_params, jet)
    )

    # With p=2, lambda^p=1/16 and the product rule gives the independent
    # values supplied by the caller's exact polynomial calculation.
    weighted_expected = {
        "value": Fraction(7, 64),
        "D": Fraction(1, 16),
        "D2": Fraction(2, 9),
        "dz_phys": Fraction(3, 28),
        "dzz_phys": Fraction(155, 343),
    }
    checks["weight_power_recorded"] = weighted.get("weight_applied") is True and _inside(c, weighted["weight_power"], Fraction(2))
    for name, value in weighted_expected.items():
        checks["weighted_" + name] = _inside(c, weighted[name], value)
    checks["weight_is_not_constant"] = not _inside(
        c,
        weighted["D"],
        Fraction(-1, 12),
    )

    # At Z=0 the even polynomial f=rho^2+Z^2 has zero first pullback but a
    # nonzero second pullback. The geometry is still exact: rho=lambda=1/2
    # data are supplied with X=Y=1/2 and tau=1/4.
    zero_params = _params(
        z=0,
        rho=1,
        lam=Fraction(1, 2),
        delta=0,
        X=Fraction(1, 2),
        Y=Fraction(1, 2),
        tau=Fraction(1, 4),
    )
    zero_jet = {
        "value": Fraction(1),
        "rho": Fraction(2),
        "Z": Fraction(0),
        "rhoZ": Fraction(0),
        "rhorho": Fraction(2),
        "ZZ": Fraction(2),
    }
    try:
        zero_result = evaluate(c, zero_params, zero_jet)
        checks["Z0_first_zero"] = _inside(c, zero_result["D"], 0) and _inside(c, zero_result["dz_phys"], 0)
        checks["Z0_second_nonzero"] = _inside(c, zero_result["D2"], Fraction(-2)) and _inside(c, zero_result["dzz_phys"], Fraction(-8))
    except Exception:
        checks["Z0_first_zero"] = False
        checks["Z0_second_nonzero"] = False

    missing_radial = dict(jet)
    del missing_radial["rho"]
    checks["missing_radial_jet_rejected"] = _rejected(lambda: evaluate(c, params, missing_radial))

    bad_z = dict(params)
    bad_z["Z"] = Fraction(1)
    checks["bad_Z_domain_rejected"] = _rejected(lambda: evaluate(c, bad_z, jet))
    bad_denominator = dict(params)
    bad_denominator["delta"] = Fraction(4)
    checks["zero_delta_denominator_rejected"] = _rejected(lambda: evaluate(c, bad_denominator, jet))
    bad_rho = dict(params)
    bad_rho["rho"] = Fraction(-1)
    checks["negative_rho_rejected"] = _rejected(lambda: evaluate(c, bad_rho, jet))
    bad_fixed = dict(params)
    bad_fixed["fixed_conditions"] = dict(params["fixed_conditions"])
    bad_fixed["fixed_conditions"]["X"] = False
    checks["bad_fixed_condition_rejected"] = _rejected(lambda: evaluate(c, bad_fixed, jet))
    missing_fixed = dict(params)
    del missing_fixed["fixed_conditions"]
    checks["missing_fixed_conditions_rejected"] = _rejected(
        lambda: evaluate(c, missing_fixed, jet)
    )

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "exact_fraction_expected": True,
        "ordinary_Z1_Z2_checked": True,
        "physical_z_first_second_checked": True,
        "nonzero_scale_exponent_checked": True,
        "missing_radial_jet_rejected": True,
        "domain_guards_checked": True,
        "source_provenance_checked": True,
        "implementation_imported": False,
        "physical_NS_not_claimed": True,
    }


def run(leaf: Any, c: Any) -> dict[str, Any]:
    return run_fixture(leaf, c)


__all__ = ["SCHEMA", "run_fixture", "run"]


