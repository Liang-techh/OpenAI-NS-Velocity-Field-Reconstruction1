"""Synthetic checks for the signed physical pulse correction kernel.

This fixture is intentionally source-free.  It uses two independent, smooth,
nonzero field germs and unequal physical radii; it is not a saved source,
source binding, Pstar value, unit assignment, or pulse-equality proof.  The
script is meant to be run once by the parent after the real source entry is
frozen.  It does not import any previous builder, provider, or checker.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np
from scipy.integrate import quad_vec

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from openai_ns_reconstruction.pulse_signed_m15_correction import (  # noqa: E402
    ChartSpec,
    FieldJets,
    FormulaParameters,
    PulseSide,
    QuadratureSettings,
    evaluate_signed_pulse_correction,
)


_CHART_SCALE = {"main": 1.0, "gap": 0.83, "end": 1.21}


def _chart_y(chart: str, x: float, mu: float) -> float:
    """Synthetic physical y coordinate for each chart's source query."""

    if chart in ("main", "gap"):
        return x / mu
    return 13.0 / mu + x


def _make_field(side_scale: float, mu: float) -> tuple[Callable[..., FieldJets], Callable[..., float]]:
    def field(chart: str, x: float, Z: float) -> FieldJets:
        c = _CHART_SCALE[chart]
        y = _chart_y(chart, x, mu)
        # U is an explicit smooth polynomial in physical y and slow Z.
        u0 = side_scale * c * (0.17 + 0.041 * y + 0.013 * y * y)
        u1 = side_scale * c * (0.11 - 0.027 * y + 0.009 * y * y)
        u2 = side_scale * c * 0.019
        # Q is a physical-y polynomial.  The Q_y, Q_Z, Q_ZZ and Q_yZ below
        # are its actual analytic partial derivatives, not independent jets.
        q0 = side_scale * c * (0.39 + 0.063 * y + 0.017 * y * y)
        q1 = side_scale * c * (0.16 + 0.031 * y - 0.008 * y * y)
        q2 = side_scale * c * 0.027
        qy0 = side_scale * c * (0.063 + 0.034 * y)
        qy1 = side_scale * c * (0.031 - 0.016 * y)
        qy2 = 0.0
        return FieldJets(
            U=u0 + Z * u1 + Z * Z * u2,
            U_Z=u1 + 2.0 * Z * u2,
            Q=q0 + Z * q1 + Z * Z * q2,
            Q_Z=q1 + 2.0 * Z * q2,
            Q_ZZ=2.0 * q2,
            Q_y=qy0 + Z * qy1 + Z * Z * qy2,
            Q_yZ=qy1 + 2.0 * Z * qy2,
        )

    def q_yy(chart: str, x: float, Z: float) -> float:
        # This is supplied separately to the independent N1p checker; it is
        # not read back from the kernel's FieldJets interface.
        return side_scale * _CHART_SCALE[chart] * (0.034 - 0.016 * Z)

    return field, q_yy


def _make_radius(base: float, mu: float):
    def radius(chart: str, x: float) -> float:
        # dR/dy=R exactly: y=x/mu on main/gap and y=13/mu+x on end.
        return base * math.exp(_chart_y(chart, x, mu))

    return radius


def _zero_field(chart: str, x: float, Z: float) -> FieldJets:
    return FieldJets(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _zero_q_yy(chart: str, x: float, Z: float) -> float:
    return 0.0


def _charts(mu: float) -> tuple[ChartSpec, ...]:
    return (
        ChartSpec(
            "main",
            "xi",
            mu=mu,
            mu_z_independent=True,
            lower=0.0,
            upper=11.0,
            endpoints_z_independent=True,
        ),
        ChartSpec(
            "gap",
            "xi",
            mu=mu,
            mu_z_independent=True,
            lower=11.0,
            upper=13.0 - 4.0 * mu,
            endpoints_z_independent=True,
        ),
        ChartSpec(
            "end",
            "s",
            mu=None,
            mu_z_independent=True,
            lower=-4.0,
            upper=0.0,
            endpoints_z_independent=True,
        ),
    )


def _evaluate(
    physical: PulseSide,
    frame: PulseSide,
    Z: float,
    *,
    mu: float = 2.5,
):
    return evaluate_signed_pulse_correction(
        physical_side=physical,
        frame_side=frame,
        charts=_charts(mu),
        parameters=FormulaParameters(delta=0.23),
        Z=Z,
        settings=QuadratureSettings(epsabs=1.0e-12, epsrel=1.0e-12, limit=300),
    )


def _close(actual: float, expected: float, *, tolerance: float = 2.0e-8) -> bool:
    return math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance)


def _assert_close(actual: float, expected: float, label: str) -> None:
    if not _close(actual, expected):
        raise AssertionError(f"{label}: {actual!r} != {expected!r}")


def _direct_source_density(jets: FieldJets, Z: float, delta: float = 0.23) -> float:
    """Independent C0 source formula; deliberately not the kernel helper."""

    h = 1.0 - delta
    L = 1.0 - delta * Z * Z
    A = jets.Q + 0.5 * h * Z * jets.Q_Z + jets.Q_y
    B = jets.Q * (0.5 * jets.Q + jets.Q_y)
    C = -2.0 * Z * jets.Q + (1.0 - Z * Z) * jets.Q_Z - 2.0 * Z * jets.Q_y
    return A / L + B + jets.U * C / L


def _direct_n1p_side(
    side: PulseSide,
    q_yy: Callable[[str, float, float], float],
    chart: str,
    Z: float,
    mu: float,
    lower: float,
    upper: float,
) -> float:
    """Directly integrate one side's ``-R^2 N1p`` expression."""

    weight = 1.0 / mu if chart in ("main", "gap") else 1.0

    def integrand(x: float) -> np.ndarray:
        x = float(x)
        jets = side.field(chart, x, Z)
        radius = side.radius(chart, x)
        density = _direct_source_density(jets, Z)
        n1p = -0.5 * density + (q_yy(chart, x, Z) + jets.Q_y) / radius
        return np.asarray((-weight * radius * radius * n1p,), dtype=float)

    value, _ = quad_vec(
        integrand,
        lower,
        upper,
        epsabs=1.0e-11,
        epsrel=1.0e-11,
        norm="max",
        limit=300,
    )
    return float(value[0])


def _direct_n1p_pair(
    physical: PulseSide,
    physical_q_yy: Callable[[str, float, float], float],
    frame: PulseSide,
    frame_q_yy: Callable[[str, float, float], float],
    chart: str,
    Z: float,
    mu: float,
    lower: float,
    upper: float,
) -> float:
    """Pair direct integral, formed from independently integrated sides."""

    return _direct_n1p_side(
        physical, physical_q_yy, chart, Z, mu, lower, upper
    ) - _direct_n1p_side(frame, frame_q_yy, chart, Z, mu, lower, upper)


def _direct_stable_volume(
    physical: PulseSide,
    frame: PulseSide,
    chart: str,
    Z: float,
    mu: float,
    lower: float,
    upper: float,
) -> float:
    """Independent unit-measure volume integral for the weight control."""

    def integrand(x: float) -> np.ndarray:
        x = float(x)
        p = physical.field(chart, x, Z)
        f = frame.field(chart, x, Z)
        sp = _direct_source_density(p, Z)
        sf = _direct_source_density(f, Z)
        rp = physical.radius(chart, x)
        rf = frame.radius(chart, x)
        return np.asarray(
            (0.5 * (rp * rp * (sp - sf) + (rp - rf) * (rp + rf) * sf),),
            dtype=float,
        )

    value, _ = quad_vec(
        integrand,
        lower,
        upper,
        epsabs=1.0e-11,
        epsrel=1.0e-11,
        norm="max",
        limit=300,
    )
    return float(value[0])


def run_fixture() -> dict[str, Any]:
    # This is a resolved synthetic scalar.  No source-node identifier is
    # interpreted as a numerical mu handle.
    mu = 2.5
    physical_field, physical_q_yy = _make_field(1.19, mu)
    frame_field, frame_q_yy = _make_field(0.73, mu)
    physical = PulseSide(
        "P_synthetic",
        physical_field,
        _make_radius(1.71, mu),
        radius_z_independent=True,
        source_radius_exponential_law=True,
    )
    frame = PulseSide(
        "F_synthetic",
        frame_field,
        _make_radius(0.86, mu),
        radius_z_independent=True,
        source_radius_exponential_law=True,
    )
    Z = 0.29
    result = _evaluate(physical, frame, Z, mu=mu)
    zero = PulseSide(
        "zero_synthetic",
        _zero_field,
        _make_radius(1.0, mu),
        radius_z_independent=True,
        source_radius_exponential_law=True,
    )
    physical_only = _evaluate(physical, zero, Z, mu=mu)
    frame_only_signed = _evaluate(zero, frame, Z, mu=mu)
    chart_map = result.chart_map
    if tuple(sorted(chart_map)) != ("end", "gap", "main"):
        raise AssertionError("all three charts must be returned")

    controls: dict[str, bool] = {}
    controls["synthetic_jets_are_distinct"] = (
        physical.field("main", 0.13, Z) != frame.field("main", 0.13, Z)
    )
    controls["radii_are_unequal"] = physical.radius("main", 0.13) != frame.radius("main", 0.13)
    controls["radius_exponential_law_declared"] = (
        physical.source_radius_exponential_law
        and frame.source_radius_exponential_law
    )
    if not controls["synthetic_jets_are_distinct"] or not controls["radii_are_unequal"]:
        raise AssertionError("fixture must exercise two unequal nonzero sides")
    controls["directed_endpoints_supported"] = (
        chart_map["main"].source_lower == 0.0
        and chart_map["main"].source_upper == 11.0
        and chart_map["gap"].source_lower == 11.0
        and chart_map["gap"].source_upper == 13.0 - 4.0 * mu
        and chart_map["end"].source_lower == -4.0
        and chart_map["end"].source_upper == 0.0
    )

    # The three chart contributions remain separate and sum to the signed
    # aggregate.  Each chart keeps its own physical endpoint term.
    _assert_close(
        sum(chart.m15_X for chart in result.charts),
        result.signed_m15,
        "three-chart signed split",
    )
    _assert_close(
        sum(chart.m15_X_Z1 for chart in result.charts),
        result.signed_m15_Z1,
        "three-chart Z1 split",
    )
    controls["three_chart_split"] = len(result.charts) == 3
    controls["main_gap_weight_is_inverse_mu"] = (
        chart_map["main"].physical_weight == 1.0 / mu
        and chart_map["gap"].physical_weight == 1.0 / mu
    )
    controls["end_weight_is_ds"] = chart_map["end"].physical_weight == 1.0

    # Independently integrate the original N1p expression.  Since dR/dy=R,
    # -integral R^2*N1p equals the endpoint-subtracted m15 formula.
    for chart in result.charts:
        direct_physical = _direct_n1p_side(
            physical,
            physical_q_yy,
            chart.chart,
            Z,
            mu,
            chart.source_lower,
            chart.source_upper,
        )
        direct_frame = _direct_n1p_side(
            frame,
            frame_q_yy,
            chart.chart,
            Z,
            mu,
            chart.source_lower,
            chart.source_upper,
        )
        direct_n1p = _direct_n1p_pair(
            physical,
            physical_q_yy,
            frame,
            frame_q_yy,
            chart.chart,
            Z,
            mu,
            chart.source_lower,
            chart.source_upper,
        )
        _assert_close(
            physical_only.chart_map[chart.chart].m15_X,
            direct_physical,
            f"{chart.chart} physical-side direct original N1p",
        )
        _assert_close(
            -frame_only_signed.chart_map[chart.chart].m15_X,
            direct_frame,
            f"{chart.chart} frame-side direct original N1p",
        )
        _assert_close(
            direct_physical - direct_frame,
            direct_n1p,
            f"{chart.chart} independently subtracted N1p sides",
        )
        _assert_close(
            chart.m15_X,
            direct_n1p,
            f"{chart.chart} direct original N1p integration",
        )
        direct_unit_measure = _direct_stable_volume(
            physical,
            frame,
            chart.chart,
            Z,
            mu,
            chart.source_lower,
            chart.source_upper,
        )
        _assert_close(
            chart.integral,
            chart.physical_weight * direct_unit_measure,
            f"{chart.chart} physical measure",
        )
    controls["physical_weights_checked"] = True
    controls["integration_by_parts_checked"] = True

    for chart in result.charts:
        _assert_close(
            chart.integral,
            chart.delta_integral + chart.radius_integral,
            f"{chart.chart} integral decomposition",
        )
        _assert_close(
            chart.boundary,
            chart.boundary_delta_upper
            + chart.boundary_radius_upper
            - chart.boundary_delta_lower
            - chart.boundary_radius_lower,
            f"{chart.chart} signed boundary",
        )
        _assert_close(
            chart.m15_X,
            chart.integral - chart.boundary,
            f"{chart.chart} m15 endpoint subtraction",
        )
    controls["boundary_contribution_checked"] = all(
        abs(chart.boundary) > 1.0e-8 for chart in result.charts
    )
    if not controls["boundary_contribution_checked"]:
        raise AssertionError("synthetic endpoint contribution unexpectedly vanished")

    # Independent central differences verify the analytic ordinary Z1 result
    # for each chart and for the signed aggregate.
    step = 1.0e-6
    plus = _evaluate(physical, frame, Z + step)
    minus = _evaluate(physical, frame, Z - step)
    derivative_errors = []
    for chart in result.charts:
        finite_difference = (
            plus.chart_map[chart.chart].m15_X - minus.chart_map[chart.chart].m15_X
        ) / (2.0 * step)
        derivative_errors.append(abs(finite_difference - chart.m15_X_Z1))
        _assert_close(
            finite_difference,
            chart.m15_X_Z1,
            f"{chart.chart} independent central difference",
        )
    aggregate_finite_difference = (plus.signed_m15 - minus.signed_m15) / (2.0 * step)
    _assert_close(
        aggregate_finite_difference,
        result.signed_m15_Z1,
        "aggregate independent central difference",
    )
    controls["ordinary_Z1_center_difference"] = True

    # The stable pair must vanish exactly when both ordered sides are the same
    # callable bundle.  This is an algebraic control, not a seam assumption.
    same = _evaluate(physical, physical, Z)
    controls["same_side_exact_zero"] = (
        same.signed_m15 == 0.0
        and same.signed_m15_Z1 == 0.0
        and all(chart.m15_X == 0.0 and chart.m15_X_Z1 == 0.0 for chart in same.charts)
    )
    if not controls["same_side_exact_zero"]:
        raise AssertionError("same-side stable signed pair must be exact zero")

    # A missing Z-independence declaration is rejected before any quadrature;
    # this prevents silently dropping radius Leibniz terms.
    rejected = PulseSide(
        "unbound_radius",
        physical.field,
        physical.radius,
        radius_z_independent=False,
        source_radius_exponential_law=False,
    )
    try:
        _evaluate(rejected, frame, Z)
    except ValueError:
        controls["radius_Z_independence_required"] = True
    else:
        raise AssertionError("undeclared radius Z independence was accepted")

    if abs(result.signed_m15) < 1.0e-8 or abs(result.signed_m15_Z1) < 1.0e-8:
        raise AssertionError("fixture must retain nonzero value and ordinary Z1")
    if all(
        result.charts[index].m15_X == result.charts[0].m15_X
        for index in range(1, len(result.charts))
    ):
        raise AssertionError("chart numbers should differ in the synthetic fixture")

    report = {
        "passed": all(controls.values()),
        "synthetic_only": True,
        "source_values_selected": False,
        "pstar_repeated": False,
        "pulse_equality_proved": result.pulse_equality_proved,
        "charts": [chart.chart for chart in result.charts],
        "signed_m15": result.signed_m15,
        "signed_m15_Z1": result.signed_m15_Z1,
        "maximum_Z1_center_difference_error": max(derivative_errors),
        "quadrature_is_certification": result.quadrature_is_certification,
        "controls": controls,
    }
    if not report["passed"]:
        raise AssertionError(f"synthetic pulse fixture failed: {report}")
    return report


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))
