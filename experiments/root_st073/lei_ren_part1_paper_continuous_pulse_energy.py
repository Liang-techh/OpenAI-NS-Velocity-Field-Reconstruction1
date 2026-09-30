"""Continuous MP energy atom for the Section 7.5 axial pulse.

The source definition is delegated to ``ContinuousAxialPulse.value_jet``.
Only the plateau integral is reduced analytically; startup and cutoff use MP
Gauss-Legendre nodes.  Startup values retain the source class's nested MP
primitive evaluation.  This module exposes ``continuous_pulse_energy`` for a
future continuous solver but does not install it into that solver.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_axial_pulse import pulse_K_p
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse


OUT_JSON = Path(__file__).with_suffix(".json")
# Keep decimal boundaries as text until an active MP work context exists.
# Constructing these at module import would round .02 at the default precision.
STARTUP_END = "0.02"
PLATEAU_END = "10"
SUPPORT_END = "11"


def _n(value: mp.mpf, digits: int = 80) -> str:
    return mp.nstr(mp.mpf(value), digits)


def _signed_log(value: mp.mpf, digits: int = 80) -> dict[str, object]:
    value = mp.mpf(value)
    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return {
        "sign": 1 if value > 0 else -1,
        "log_abs": _n(mp.log(abs(value)), digits),
        "arbitrary_exponent_value": _n(value, digits),
    }


def _relative(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    left = mp.mpf(left)
    right = mp.mpf(right)
    scale = max(abs(left), abs(right))
    return abs(left - right) / scale if scale else abs(left - right)


def _gauss_rule(order: int) -> tuple[list[mp.mpf], list[mp.mpf]]:
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    return (
        [mp.mpf(nodes[index]) for index in range(order)],
        [mp.mpf(weights[index]) for index in range(order)],
    )


def _gauss_integral(function, left: mp.mpf, right: mp.mpf, order: int) -> mp.mpf:
    """MP Gauss-Legendre integral with MP nodes and weights only."""

    nodes, weights = _gauss_rule(order)
    midpoint = (left + right) / 2
    half_width = (right - left) / 2
    return mp.fsum(
        weight * function(midpoint + half_width * node)
        for node, weight in zip(nodes, weights)
    ) * half_width


def _plateau_energy() -> mp.mpf:
    """Analytic integral on [.02, 10], where gp(xi) = xi - .01."""

    shift = mp.mpf(".01")
    startup_end = mp.mpf(STARTUP_END)
    plateau_end = mp.mpf(PLATEAU_END)

    def antiderivative(xi: mp.mpf) -> mp.mpf:
        y = xi - shift
        polynomial = y * y / 2 + y / 2 + mp.mpf(".25")
        return -mp.exp(-2 * xi) * polynomial

    return antiderivative(plateau_end) - antiderivative(startup_end)


def _continuous_parts(*, precision: int, quadrature_order: int) -> dict[str, mp.mpf]:
    with mp.workdps(precision):
        pulse = ContinuousAxialPulse(precision=precision)
        startup_end = mp.mpf(STARTUP_END)
        plateau_end = mp.mpf(PLATEAU_END)
        support_end = mp.mpf(SUPPORT_END)

        # Calling value_jet here is intentional: startup uses the exact same
        # MP primitive and cutoff convention as the shared source definition.
        def weighted_value(xi: mp.mpf) -> mp.mpf:
            value = pulse.value_jet(xi)["value"]
            return mp.exp(-2 * xi) * value * value

        startup = _gauss_integral(
            weighted_value, mp.mpf(0), startup_end, quadrature_order
        )
        plateau = _plateau_energy()
        cutoff = _gauss_integral(
            weighted_value, plateau_end, support_end, quadrature_order
        )
        return {
            "startup": startup,
            "plateau": plateau,
            "cutoff": cutoff,
            "total": startup + plateau + cutoff,
        }


def continuous_pulse_energy(
    precision: int = 100, *, quadrature_order: int | None = None
) -> mp.mpf:
    """Return MP ``int_0^11 exp(-2 xi) gp(xi)^2 dxi``.

    ``gp`` is exactly the source pulse represented by
    ``ContinuousAxialPulse.value_jet``.  The optional order is useful for
    bounded refinement; callers normally need only ``precision``.
    """

    try:
        precision = int(precision)
    except (TypeError, ValueError) as exc:
        raise ValueError("precision must be an integer at least 40") from exc
    if precision < 40:
        raise ValueError("precision must be at least 40")
    if quadrature_order is None:
        # MP order grows with requested precision while remaining bounded for
        # practical source diagnostics.  The startup primitive is nested, so
        # this is intentionally conservative rather than a global proof order.
        quadrature_order = max(48, min(128, precision * 4 // 5))
    try:
        quadrature_order = int(quadrature_order)
    except (TypeError, ValueError) as exc:
        raise ValueError("quadrature_order must be an integer at least 16") from exc
    if quadrature_order < 16:
        raise ValueError("quadrature_order must be at least 16")
    with mp.workdps(precision):
        return _continuous_parts(
            precision=precision, quadrature_order=quadrature_order
        )["total"]


def _receipt_value(parts: dict[str, mp.mpf], digits: int) -> dict[str, object]:
    return {
        "startup": _signed_log(parts["startup"], digits),
        "plateau": _signed_log(parts["plateau"], digits),
        "cutoff": _signed_log(parts["cutoff"], digits),
        "total": _signed_log(parts["total"], digits),
    }


def run() -> dict[str, object]:
    print("continuous pulse energy: precision 70", flush=True)
    with mp.workdps(120):
        # Keep order fixed here so this comparison isolates MP precision.
        parts70 = _continuous_parts(precision=70, quadrature_order=80)
        print("continuous pulse energy: precision 100", flush=True)
        parts100 = _continuous_parts(precision=100, quadrature_order=80)
        print("continuous pulse energy: precision 100 refined order", flush=True)
        parts100_refined = _continuous_parts(precision=100, quadrature_order=112)

        old_float = mp.mpf(str(pulse_K_p(quadrature_order=192)))
        value70 = parts70["total"]
        value100 = parts100["total"]
        value100_refined = parts100_refined["total"]
        digits = 90
        receipt: dict[str, object] = {
            "kind": "continuous_axial_pulse_energy",
            "definition": "integral_0^11 exp(-2*xi) * gp(xi)^2 dxi",
            "source_definition": "ContinuousAxialPulse.value_jet",
            "support": ["0", "11"],
            "splits": {
                "startup": "[0, .02] via MP Gauss-Legendre of value_jet (nested MP primitive)",
                "plateau": "[.02, 10] analytic antiderivative of exp(-2*xi)*(xi-.01)^2",
                "cutoff": "[10, 11] via MP Gauss-Legendre of value_jet",
            },
            "nodes": "All quadrature nodes and weights are mpmath values; no float nodes.",
            "precision_70_order_80": _receipt_value(parts70, digits),
            "precision_100_order_80": _receipt_value(parts100, digits),
            "precision_100_order_112": _receipt_value(parts100_refined, digits),
            "precision_refinement": {
                "log_total_precision_70_minus_100_same_order_80": _n(
                    mp.log(value70) - mp.log(value100), digits
                ),
                "relative_total_precision_70_vs_100_same_order_80": _signed_log(
                    _relative(value70, value100), digits
                ),
                "log_total_order_80_minus_112": _n(
                    mp.log(value100) - mp.log(value100_refined), digits
                ),
                "relative_total_order_80_vs_112": _signed_log(
                    _relative(value100, value100_refined), digits
                ),
            },
            "legacy_float_comparison": {
                "pulse_K_p_order_192": _signed_log(old_float, digits),
                "continuous_order_112_minus_legacy_float": _n(
                    value100_refined - old_float, digits
                ),
                "relative_continuous_vs_legacy_float": _signed_log(
                    _relative(value100_refined, old_float), digits
                ),
            },
            "public_function": {
                "name": "continuous_pulse_energy",
                "signature": "continuous_pulse_energy(precision=100, quadrature_order=None)",
                "returned_value_type": "mpmath.mpf",
                "installed_in_global_profile": False,
            },
            "scope": {
                "measured": [
                    "bounded MP quadrature/refinement for the source pulse energy atom",
                    "startup, analytic plateau, and cutoff contributions",
                    "comparison with legacy float pulse_K_p(order192)",
                ],
                "not_certified": [
                    "rigorous quadrature enclosure",
                    "global finite energy",
                    "continuous incoming moment or axial solve",
                    "global profile installation or recursive closure",
                ],
            },
        }

    OUT_JSON.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_JSON}", flush=True)
    print("global installation: false; quadrature enclosure: uncertified", flush=True)
    return receipt


if __name__ == "__main__":
    run()
