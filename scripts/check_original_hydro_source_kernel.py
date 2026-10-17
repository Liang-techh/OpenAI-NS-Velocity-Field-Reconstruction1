"""Focused independent checks for the pure hydro ordinary-derivative kernel."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_hydro_source_kernel import (  # noqa: E402
    evaluate_ordinary_sectors,
    hydro_drive_ordinary,
    hydro_drive_ordinary_epsilon,
    hydro_drive_ordinary_epsilon_sectors,
    hydro_moment_difference_numerator_ordinary,
    hydro_moment_difference_ordinary_epsilon_sectors,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (TypeError, ValueError, ArithmeticError):
        return True
    return False


def _ordinary_rows(ctx: Any, function: Callable[[Any], Any], z: Any, order: int) -> list[Any]:
    return [ctx.diff(function, z, k) for k in range(order + 1)]


def _full_hydro(
    ctx: Any,
    z: Any,
    delta: Any,
    v: Callable[[Any], Any],
    m: Callable[[Any], Any],
    a: Callable[[Any], Any],
) -> Any:
    """Independent point formula used only as a differentiation oracle."""
    d = 1 - z * z
    L = 1 - delta * z * z
    m_z = ctx.diff(m, z)
    a_z = ctx.diff(a, z)
    W = 1 - (1 - delta) * z * m(z) - d * m_z
    numerator = (
        -W * v(z)
        + (1 - delta) * (m(z) - z * m_z) / 2
        + 2 * delta * z * a(z)
        - d * a_z
    )
    return numerator / (2 * L)


def _close(ctx: Any, left: Any, right: Any, tolerance: str = "1e-90") -> bool:
    return abs(left - right) < ctx.mpf(tolerance)


def _inside(point: Any, enclosure: Any) -> bool:
    return point in enclosure


def _box(ctx: Any, interval_ctx: Any, value: Any, radius: str = "1e-105") -> Any:
    width = ctx.mpf(radius)
    return interval_ctx.mpf(
        [ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)]
    )


def _nonzero(ctx: Any, row: list[Any], tolerance: str = "1e-70") -> bool:
    return any(abs(value) > ctx.mpf(tolerance) for value in row)


def run_fixture() -> dict[str, Any]:
    ctx = mp.mp.clone()
    ctx.dps = 140
    interval_ctx = MPIntervalContext()
    interval_ctx.dps = 180
    controls: list[dict[str, Any]] = []
    differences: list[Any] = []

    def passed(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    delta = ctx.mpf(".23")

    def v(z: Any) -> Any:
        return ctx.exp(z / 5) + ctx.sin(z / 3) + z**3 / 17

    def m(z: Any) -> Any:
        return 1 + ctx.exp(-z / 4) + z**2 / 11

    def a(z: Any) -> Any:
        return ctx.sin(z / 2) + ctx.exp(z / 7) / 9 + z**3 / 13

    for z_text in ("-.31", ".07", ".29"):
        z = ctx.mpf(z_text)
        v_row = _ordinary_rows(ctx, v, z, 2)
        m_row = _ordinary_rows(ctx, m, z, 3)
        a_row = _ordinary_rows(ctx, a, z, 3)
        got = hydro_drive_ordinary(z, delta, v_row, m_row, a_row)
        oracle = [ctx.diff(lambda q: _full_hydro(ctx, q, delta, v, m, a), z, k) for k in range(3)]
        for order, (value, expected) in enumerate(zip(got, oracle)):
            difference = abs(value - expected)
            differences.append(difference)
            passed(f"ordinary_full_formula_{z_text}_Z{order}", difference < ctx.mpf("1e-80"))

        iv_z = _box(ctx, interval_ctx, z)
        iv_delta = _box(ctx, interval_ctx, delta)
        iv_v = [_box(ctx, interval_ctx, value) for value in v_row]
        iv_m = [_box(ctx, interval_ctx, value) for value in m_row]
        iv_a = [_box(ctx, interval_ctx, value) for value in a_row]
        directed = hydro_drive_ordinary(iv_z, iv_delta, iv_v, iv_m, iv_a)
        passed(
            f"directed_containment_{z_text}",
            all(_inside(value, enclosure) for value, enclosure in zip(got, directed)),
        )

    z = ctx.mpf(".17")
    v_row = _ordinary_rows(ctx, v, z, 2)
    m_row = _ordinary_rows(ctx, m, z, 3)
    a_row = _ordinary_rows(ctx, a, z, 3)
    passed("missing_v_Z2_rejected", _raises(lambda: hydro_drive_ordinary(z, delta, v_row[:2], m_row, a_row)))
    passed("missing_m_Z3_rejected", _raises(lambda: hydro_drive_ordinary(z, delta, v_row, m_row[:3], a_row)))
    passed("missing_a_Z3_rejected", _raises(lambda: hydro_drive_ordinary(z, delta, v_row, m_row, a_row[:3])))
    passed("missing_m_rejected", _raises(lambda: hydro_drive_ordinary(z, delta, v_row, None, a_row)))
    passed("missing_a_rejected", _raises(lambda: hydro_drive_ordinary(z, delta, v_row, m_row, None)))

    # The baseline is nonzero, but its hydro numerator cancels identically by
    # construction.  The perturbation sectors then retain a genuine degree-2
    # m*v cross term.
    def m0(q: Any) -> Any:
        return 1 + q / 9 + q**2 / 31

    def v0(q: Any) -> Any:
        d = 1 - q * q
        m_z = ctx.diff(m0, q)
        w = 1 - (1 - delta) * q * m0(q) - d * m_z
        return (1 - delta) * (m0(q) - q * m_z) / (2 * w)

    def a0(q: Any) -> Any:
        return ctx.mpf(0)

    def v1(q: Any) -> Any:
        return ctx.exp(q / 6) + q / 8

    def m1(q: Any) -> Any:
        return ctx.sin(q / 4) + q**2 / 19

    def a1(q: Any) -> Any:
        return ctx.exp(-q / 5) + ctx.sin(q / 7)

    def make_sectors(z_value: Any) -> tuple[dict[int, list[Any]], dict[int, list[Any]], dict[int, list[Any]]]:
        return (
            {0: _ordinary_rows(ctx, v0, z_value, 2), 1: _ordinary_rows(ctx, v1, z_value, 2)},
            {0: _ordinary_rows(ctx, m0, z_value, 3), 1: _ordinary_rows(ctx, m1, z_value, 3)},
            {0: _ordinary_rows(ctx, a0, z_value, 3), 1: _ordinary_rows(ctx, a1, z_value, 3)},
        )

    z = ctx.mpf(".21")
    v_sectors, m_sectors, a_sectors = make_sectors(z)
    decomposition = hydro_drive_ordinary_epsilon_sectors(
        z,
        delta,
        v_sectors,
        m_sectors,
        a_sectors,
    )
    numerator_sectors = decomposition["numerator_sectors"]
    denominator_sectors = decomposition["denominator_sectors"]
    passed("sector_result_has_numerator_and_denominator", set(decomposition) == {"numerator_sectors", "denominator_sectors"})
    passed("quadratic_cross_sector_is_retained", 2 in numerator_sectors and _nonzero(ctx, numerator_sectors[2]))
    passed("vanishing_baseline_cancellation", 0 in numerator_sectors and all(abs(value) < ctx.mpf("1e-105") for value in numerator_sectors[0]))
    difference_decomposition = hydro_moment_difference_ordinary_epsilon_sectors(
        z,
        delta,
        v_sectors,
        m_sectors,
        a_sectors,
    )
    difference_epsilon = ctx.mpf(".35")
    direct_difference = hydro_moment_difference_numerator_ordinary(
        z,
        delta,
        evaluate_ordinary_sectors(v_sectors, difference_epsilon),
        evaluate_ordinary_sectors(m_sectors, difference_epsilon),
        evaluate_ordinary_sectors(a_sectors, difference_epsilon),
    )
    sector_difference = evaluate_ordinary_sectors(
        difference_decomposition["numerator_sectors"], difference_epsilon
    )
    passed(
        "moment_difference_sector_matches_direct_numerator",
        all(_close(ctx, left, right, "1e-100") for left, right in zip(sector_difference, direct_difference)),
    )
    passed(
        "moment_difference_keeps_common_denominator",
        all(
            _close(
                ctx,
                left,
                right,
                "1e-100",
            )
            for left, right in zip(
                evaluate_ordinary_sectors(difference_decomposition["denominator_sectors"], difference_epsilon),
                evaluate_ordinary_sectors(denominator_sectors, difference_epsilon),
            )
        ),
    )

    def v_eps(epsilon: Any, q: Any) -> Any:
        return v0(q) + epsilon * v1(q)

    def m_eps(epsilon: Any, q: Any) -> Any:
        return m0(q) + epsilon * m1(q)

    def a_eps(epsilon: Any, q: Any) -> Any:
        return a0(q) + epsilon * a1(q)

    for epsilon in (ctx.mpf("-.4"), ctx.mpf("0"), ctx.mpf(".35"), ctx.mpf(".9")):
        evaluated_v = evaluate_ordinary_sectors(v_sectors, epsilon)
        evaluated_m = evaluate_ordinary_sectors(m_sectors, epsilon)
        evaluated_a = evaluate_ordinary_sectors(a_sectors, epsilon)
        direct = hydro_drive_ordinary(z, delta, evaluated_v, evaluated_m, evaluated_a)
        sector = hydro_drive_ordinary_epsilon(
            z,
            delta,
            v_sectors,
            m_sectors,
            a_sectors,
            epsilon,
        )
        full = [
            ctx.diff(
                lambda q: _full_hydro(
                    ctx,
                    q,
                    delta,
                    lambda r: v_eps(epsilon, r),
                    lambda r: m_eps(epsilon, r),
                    lambda r: a_eps(epsilon, r),
                ),
                z,
                k,
            )
            for k in range(3)
        ]
        passed(
            f"sector_matches_direct_kernel_{epsilon}",
            all(_close(ctx, left, right, "1e-100") for left, right in zip(sector, direct)),
        )
        passed(
            f"sector_matches_independent_full_formula_{epsilon}",
            all(_close(ctx, left, right, "1e-100") for left, right in zip(sector, full)),
        )

    # Exercise the source-correlated path delta=epsilon*exp(970) at several
    # finite negative epsilons, which keep the common denominator positive.
    correlated_scale = ctx.exp(970)
    correlated_delta_sectors = {1: correlated_scale}
    for epsilon in (ctx.mpf("-.001"), ctx.mpf("-.0007"), ctx.mpf("-.0015")):
        correlated_delta = epsilon * correlated_scale
        sector = hydro_drive_ordinary_epsilon(
            z,
            None,
            v_sectors,
            m_sectors,
            a_sectors,
            epsilon,
            delta_sectors=correlated_delta_sectors,
        )
        direct = hydro_drive_ordinary(
            z,
            correlated_delta,
            evaluate_ordinary_sectors(v_sectors, epsilon),
            evaluate_ordinary_sectors(m_sectors, epsilon),
            evaluate_ordinary_sectors(a_sectors, epsilon),
        )
        passed(
            f"correlated_delta_sector_matches_direct_{epsilon}",
            all(_close(ctx, left, right, "1e-100") for left, right in zip(sector, direct)),
        )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "maximum_supporting_difference": ctx.nstr(max(differences), 25),
        "ordinary_formula": "Dhydro=[-W*v+(1-delta)*(m-z*m_Z)/2+2*delta*z*a-d*a_Z]/(2*L)",
        "sector_contract": {
            "fields": "degree -> ordinary jet; v through Z2, m/a through Z3",
            "delta_default": "{0: delta}",
            "correlated_delta_example": "{1: exp(970)}",
            "returned": ["numerator_sectors", "denominator_sectors"],
            "division": "common 2L is divided only after sector evaluation; no series truncation",
        },
        "ancestor_constructor_producer_checker_fixture_run": False,
        "source_loading_or_adapter": False,
        "quadrature_or_integration": False,
        "temporal_recursion": False,
        "NS_residual": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
