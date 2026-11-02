"""Fresh synthetic interval checks for ``reserved_partial_wholez_v3``.

Inputs are rational synthetic packets only.  The fixture covers a clipped
cell crossing ``Z=0`` and a positive nonzero cell, ``A``/``A_Z`` and the
ordinary-Z product rule, independent three-term containment, source
path/error retention, and exact synthetic ``r2=1`` pulse cancellation.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

SCHEMA = "reserved_partial_wholez_fixture_v3"


def _q(c: Any, value: Fraction | int) -> Any:
    return c.mpf(value.numerator) / c.mpf(value.denominator) if isinstance(value, Fraction) else c.mpf(value)


def _box(c: Any, lower: Fraction | int | str, upper: Fraction | int | str) -> Any:
    lo = c.mpf(lower) if isinstance(lower, str) else _q(c, lower)
    hi = c.mpf(upper) if isinstance(upper, str) else _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _point(c: Any, value: Fraction | int | str) -> Any:
    return c.mpf(value) if isinstance(value, str) else _q(c, value)


def _contains(interval: Any, value: Any) -> bool:
    return interval.a <= value <= interval.b


def _zero(interval: Any) -> bool:
    return interval.a == 0 and interval.b == 0


def synthetic_inputs(c: Any) -> dict[str, Any]:
    """Build fresh synthetic source/J/P/response packets in context ``c``."""

    cells = [
        {
            "cell_id": "synthetic:whole-Z:zero", "source_id": "synthetic-v3-source-zero",
            "source_path": ["fixture", "whole_Z", 0, "saved"],
            "clipped_Z": _box(c, Fraction(-1, 5), Fraction(1, 5)),
            "normalized_theta_inputs": {
                "Z0": (_box(c, Fraction(1, 2), Fraction(3, 5)), _box(c, Fraction(-2, 5), Fraction(-1, 3))),
                "Z1": (_box(c, Fraction(1, 3), Fraction(2, 5)), _box(c, Fraction(1, 5), Fraction(3, 10))),
            },
            "source_packet": {
                "id": "synthetic-v3-source-zero", "path": ["fixture", "whole_Z", 0, "saved"],
                "error_classes": ["fixture_zero_cell_error"], "errors": {"radius": "retained"},
            },
            "error_classes": ["fixture_zero_cell_error"], "error_packet": {"radius": "retained"},
        },
        {
            "cell_id": "synthetic:whole-Z:nonzero", "source_id": "synthetic-v3-source-nonzero",
            "source_path": ["fixture", "whole_Z", 1, "saved"],
            "clipped_Z": _box(c, Fraction(1, 3), Fraction(1, 2)),
            "normalized_theta_inputs": {
                "Z0": (_box(c, Fraction(-3, 5), Fraction(-1, 2)), _box(c, Fraction(2, 5), Fraction(1, 2))),
                "Z1": (_box(c, Fraction(1, 4), Fraction(1, 3)), _box(c, Fraction(-1, 2), Fraction(-2, 5))),
            },
            "source_packet": {
                "id": "synthetic-v3-source-nonzero", "path": ["fixture", "whole_Z", 1, "saved"],
                "error_classes": ["fixture_nonzero_cell_error"], "errors": {"theta": "retained"},
            },
            "error_classes": ["fixture_nonzero_cell_error"], "error_packet": {"theta": "retained"},
        },
    ]
    return {
        "response": {
            "r0": _point(c, Fraction(2, 3)), "r1": _point(c, Fraction(-4, 7)), "r2": _point(c, Fraction(3, 5)),
            "z_independent": True, "x": _point(c, "1.27"), "source_path": ["fixture", "saved_response", "x=1.27"],
        },
        "cells": cells, "mu": _box(c, Fraction(-3, 5), Fraction(-3, 5)),
        "joint": {"Z0": _box(c, Fraction(-3, 10), Fraction(-1, 5)), "Z1": _box(c, Fraction(1, 6), Fraction(1, 5)), "definition": "J=-nonpulse", "source_path": ["fixture", "saved_joint", "whole_Z"]},
        "pulse": {"Z0": _box(c, Fraction(2, 5), Fraction(1, 2)), "Z1": _box(c, Fraction(1, 3), Fraction(2, 5)), "positive": True, "is_P1": False, "source_path": ["fixture", "independent_positive_pulse", "m15"]},
    }


def _point_checks(result: dict[str, Any], mp: Any) -> dict[str, bool]:
    """Check containment at rational points, never interval endpoint equality."""

    r0, r1, r2 = mp.mpf(2) / 3, -mp.mpf(4) / 7, mp.mpf(3) / 5
    e = mp.exp(-mp.mpf(1) / 5)
    z_values = (mp.mpf(1) / 20, mp.mpf(5) / 12)
    theta_centres = (
        ((mp.mpf(11) / 20, mp.mpf(11) / 30), (-mp.mpf(7) / 20, mp.mpf(1) / 4)),
        ((-mp.mpf(11) / 20, mp.mpf(7) / 24), (mp.mpf(9) / 20, -mp.mpf(9) / 20)),
    )
    J = (-mp.mpf(1) / 4, mp.mpf(11) / 60)
    P = (mp.mpf(9) / 20, mp.mpf(11) / 30)
    checks = {"A_AZ_enclose": True, "Z1_product_rule_enclose": True, "direct_difference_enclose": True}
    for index, cell in enumerate(result["cell_results"]):
        z = z_values[index]
        den = 1 + z * z
        A, AZ = e / den, -2 * z * e / den**2
        checks["A_AZ_enclose"] &= _contains(cell["A"], A) and _contains(cell["A_Z"], AZ)
        theta0, theta1 = theta_centres[index]
        p0 = AZ * theta0[0] + A * theta0[1]
        p1 = AZ * theta1[0] + A * theta1[1]
        checks["Z1_product_rule_enclose"] &= _contains(cell["theta_products"]["theta0"]["Z1"], p0) and _contains(cell["theta_products"]["theta1"]["Z1"], p1)
        own0 = -(r0 * A * theta0[0] + r1 * A * theta1[0])
        own1 = -(r0 * p0 + r1 * p1)
        checks["direct_difference_enclose"] &= _contains(cell["result"]["Z0"], own0 + r2 * J[0] + (1 - r2) * P[0])
        checks["direct_difference_enclose"] &= _contains(cell["result"]["Z1"], own1 + r2 * J[1] + (1 - r2) * P[1])
    return checks


def run_fixture(build: Any = None, c: Any = None) -> dict[str, Any]:
    """Run synthetic checks through ``build`` (or the sibling leaf)."""

    if build is None:
        from . import reserved_partial_wholez_v3 as leaf
        build = leaf.build
    elif hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")
    if c is None:
        from mpmath.ctx_iv import MPIntervalContext
        c = MPIntervalContext()
        c.dps = 100
    inputs = synthetic_inputs(c)
    try:
        result = build(c, **inputs)
    except Exception as exc:  # pragma: no cover - diagnostic mutation path
        import traceback
        return {"schema": SCHEMA, "synthetic_only": True, "passed": False, "checks": {"build": False}, "check_count": 1, "error": repr(exc), "traceback": traceback.format_exc()}

    checks = {
        "schema": result.get("schema") == "reserved_partial_wholez_v3",
        "response_three": result.get("response_count") == 3 and len(result.get("response_intervals", ())) == 3,
        "fixed_x": result.get("x_decimal") == "1.27" and result.get("x_fixed") is True,
        "open_Z_conditional": result.get("open_Z_domain") is True and result.get("conditional") is True,
        "not_Rp_or_finite_N": result.get("not_physical_Rp") is True and result.get("not_finite_N_closure") is True,
        "no_212_append": result.get("appended_to_212_ledger") is False,
        "P_not_P1_channel": result.get("P1_distinct") is True and result.get("pulse_channel") == "m15 pressure pulse channel",
        "no_cross_cell_hull": result.get("cross_cell_hull") is False and result.get("cell_count") == 2,
        "uncertainty_not_separated": result.get("signed_terms_separate_from_errors") is False and result.get("source_errors_retained") is True,
        "source_path_id_retained": all(cell.get("source_path") == cell.get("source_packet", {}).get("path") and cell.get("source_id") == cell.get("source_packet", {}).get("id") for cell in result.get("cell_results", [])),
        "response_fixed_premise": result.get("response_Z_independent") is True,
        "old_helpers_not_called": result.get("old_helpers_called") is not True,
    }
    import mpmath as mpmath_module
    mp = mpmath_module.mp.clone()
    mp.dps = 90
    checks.update(_point_checks(result, mp))

    signed_inputs = dict(inputs)
    signed_inputs["pulse"] = {"Z0": _box(c, Fraction(-1, 2), Fraction(1, 3)), "Z1": _box(c, Fraction(-2, 5), Fraction(-1, 3)), "is_P1": False}
    signed_result = build(c, **signed_inputs)
    checks["signed_pulse_and_derivative_allowed"] = signed_result.get("cell_count") == 2

    response_one = dict(inputs["response"])
    response_one["r2"] = _point(c, 1)
    one_inputs = dict(inputs)
    one_inputs["response"] = response_one
    one_result = build(c, **one_inputs)
    checks["r2_one_exact_pulse_cancel"] = all(_zero(cell["pulse_term"]["Z0"]) and _zero(cell["pulse_term"]["Z1"]) for cell in one_result["cell_results"])
    checks["range_zero_and_nonzero"] = result["cell_results"][0]["clipped_Z"].a < 0 < result["cell_results"][0]["clipped_Z"].b and result["cell_results"][1]["clipped_Z"].a > 0 and result["cell_results"][1]["A_Z"].b < 0
    checks["derivative_metadata"] = result.get("A_Z_formula") == "-2*Z*exp(-17/10-5*mu/2)/(1+Z^2)^2"
    return {
        "schema": SCHEMA, "synthetic_only": True, "independent_point_probes": True,
        "passed": all(checks.values()), "check_count": len(checks), "checks": checks,
        "result": result, "r2_one_result": one_result, "source_receipts_read": False,
        "prior_helpers_imported": False, "old_212_ledger_appended": False,
    }


def run(build: Any = None, c: Any = None) -> dict[str, Any]:
    return run_fixture(build, c)


__all__ = ["SCHEMA", "synthetic_inputs", "run_fixture", "run"]

