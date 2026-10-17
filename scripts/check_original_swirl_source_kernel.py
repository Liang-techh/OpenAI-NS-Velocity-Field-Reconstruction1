"""Focused independent checks for the pure ordinary swirl kernel."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_swirl_source_kernel import (  # noqa: E402
    dress_ordinary_jet,
    swirl_moment_mode_numerator_delta_sectors_ordinary,
    swirl_numerator_delta_sectors_ordinary,
    swirl_drive_ordinary,
    swirl_moment_modes_ordinary,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (TypeError, ValueError, ArithmeticError):
        return True
    return False


def _ordinary_rows(ctx: Any, function: Callable[[Any], Any], z: Any, order: int) -> list[Any]:
    return [ctx.diff(function, z, k) for k in range(order + 1)]


def _full_swirl(ctx: Any, z: Any, delta: Any, b: Callable[[Any], Any], c: Callable[[Any], Any]) -> Any:
    d = 1 - z * z
    L = 1 - delta * z * z
    numerator = (
        -2 * delta * z * b(z)
        + d * ctx.diff(b, z)
        + 2 * (1 + delta) * z * c(z)
        - d * ctx.diff(c, z)
    )
    return numerator / (2 * L)


def _full_swirl_numerator(
    ctx: Any,
    z: Any,
    delta: Any,
    b: Callable[[Any], Any],
    c: Callable[[Any], Any],
) -> Any:
    d = 1 - z * z
    return (
        -2 * delta * z * b(z)
        + d * ctx.diff(b, z)
        + 2 * (1 + delta) * z * c(z)
        - d * ctx.diff(c, z)
    )


def _box(ctx: Any, interval_ctx: Any, value: Any, radius: str = "1e-105") -> Any:
    width = ctx.mpf(radius)
    return interval_ctx.mpf(
        [ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)]
    )


def _inside(point: Any, enclosure: Any) -> bool:
    return point in enclosure


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

    # Independent source rows: exp, sin and polynomial terms ensure all three
    # returned derivatives exercise product-rule cross terms.
    def b(z: Any) -> Any:
        return ctx.exp(z / 5) + ctx.sin(z / 3) + z**3 / 17

    def c(z: Any) -> Any:
        return 1 + ctx.exp(-z / 4) + z**2 / 11

    for z_text in ("-.31", ".07", ".29"):
        z = ctx.mpf(z_text)
        b_row = _ordinary_rows(ctx, b, z, 3)
        c_row = _ordinary_rows(ctx, c, z, 3)
        got = swirl_drive_ordinary(z, delta, b_row, c_row)
        oracle = [ctx.diff(lambda q: _full_swirl(ctx, q, delta, b, c), z, k) for k in range(3)]
        for order, (value, expected) in enumerate(zip(got, oracle)):
            difference = abs(value - expected)
            differences.append(difference)
            passed(f"ordinary_full_formula_{z_text}_Z{order}", difference < ctx.mpf("1e-80"))

        iv_z = _box(ctx, interval_ctx, z)
        iv_delta = _box(ctx, interval_ctx, delta)
        iv_b = [_box(ctx, interval_ctx, value) for value in b_row]
        iv_c = [_box(ctx, interval_ctx, value) for value in c_row]
        directed = swirl_drive_ordinary(iv_z, iv_delta, iv_b, iv_c)
        passed(
            f"directed_containment_{z_text}",
            all(_inside(value, enclosure) for value, enclosure in zip(got, directed)),
        )

        sectors = swirl_numerator_delta_sectors_ordinary(z, b_row, c_row)
        finite_delta_numerator = [
            sectors[0][order] + delta * sectors[1][order] for order in range(3)
        ]
        numerator_oracle = [
            ctx.diff(lambda q: _full_swirl_numerator(ctx, q, delta, b, c), z, order)
            for order in range(3)
        ]
        passed(
            f"delta_sector_reconstruction_{z_text}",
            all(
                abs(value - expected) < ctx.mpf("1e-80")
                for value, expected in zip(finite_delta_numerator, numerator_oracle)
            ),
        )

    z = ctx.mpf(".17")
    b_row = _ordinary_rows(ctx, b, z, 3)
    c_row = _ordinary_rows(ctx, c, z, 3)
    passed("missing_b_Z3_rejected", _raises(lambda: swirl_drive_ordinary(z, delta, b_row[:3], c_row)))
    passed("missing_c_Z3_rejected", _raises(lambda: swirl_drive_ordinary(z, delta, b_row, c_row[:3])))
    passed("missing_b_rejected", _raises(lambda: swirl_drive_ordinary(z, delta, None, c_row)))
    passed("missing_c_rejected", _raises(lambda: swirl_drive_ordinary(z, delta, b_row, None)))

    # A nonconstant relative amplitude is dressed exactly once.  The oracle
    # differentiates the already-dressed functions, so a drive that adds an
    # external F0 derivative (or doubles the amplitude) fails this control.
    def relative_amplitude(q: Any) -> Any:
        return ctx.exp(q / 7) + q**2 / 13 + ctx.sin(q / 5) / 9

    def raw_b(q: Any) -> Any:
        return 2 + ctx.sin(q / 4) + q**3 / 19

    def raw_c(q: Any) -> Any:
        return 1 + ctx.exp(-q / 6) + q**2 / 23

    def dressed_b(q: Any) -> Any:
        return relative_amplitude(q) * raw_b(q)

    def dressed_c(q: Any) -> Any:
        return relative_amplitude(q) * raw_c(q)

    z = ctx.mpf(".13")
    amplitude_row = _ordinary_rows(ctx, relative_amplitude, z, 3)
    raw_b_row = _ordinary_rows(ctx, raw_b, z, 3)
    raw_c_row = _ordinary_rows(ctx, raw_c, z, 3)
    dressed_b_row = dress_ordinary_jet(amplitude_row, raw_b_row)
    dressed_c_row = dress_ordinary_jet(amplitude_row, raw_c_row)
    dressed = swirl_drive_ordinary(z, delta, dressed_b_row, dressed_c_row)
    dressed_oracle = [
        ctx.diff(lambda q: _full_swirl(ctx, q, delta, dressed_b, dressed_c), z, k)
        for k in range(3)
    ]
    passed(
        "nonconstant_amplitude_dressed_once",
        all(abs(value - expected) < ctx.mpf("1e-80") for value, expected in zip(dressed, dressed_oracle)),
    )
    wrong_external_amplitude = swirl_drive_ordinary(
        z,
        delta,
        dress_ordinary_jet(amplitude_row, dressed_b_row),
        dress_ordinary_jet(amplitude_row, dressed_c_row),
    )
    passed(
        "external_amplitude_derivative_not_doubled",
        max(abs(value - expected) for value, expected in zip(wrong_external_amplitude, dressed))
        > ctx.mpf("1e-6"),
    )

    # Same-source modal directions.  Keep one target row for the base and use
    # only its correlated delta-C / delta-B directions, so there is no
    # independent interval subtraction of the base target.
    def target_square(q: Any) -> Any:
        return 1 + q / 6 + ctx.exp(-q / 8) + q**2 / 29

    def d_b(q: Any) -> Any:
        return ctx.sin(q / 5) + q**2 / 31

    def d_c(q: Any) -> Any:
        return ctx.exp(q / 9) / 7 + q**3 / 37

    z = ctx.mpf(".24")
    target_row = _ordinary_rows(ctx, target_square, z, 3)
    d_b_row = _ordinary_rows(ctx, d_b, z, 3)
    d_c_row = _ordinary_rows(ctx, d_c, z, 3)
    modes = swirl_moment_modes_ordinary(z, delta, target_row, d_b_row, d_c_row)
    passed("three_modal_direction_numerators", len(modes) == 3 and all(len(row) == 3 for row in modes))

    mode_functions = (
        (lambda q: target_square(q) / 2, target_square),
        (lambda q: ctx.mpf(0), d_c),
        (d_b, lambda q: ctx.mpf(0)),
    )
    for mode, ((b_mode, c_mode), row) in enumerate(zip(mode_functions, modes)):
        oracle = [
            ctx.diff(lambda q: _full_swirl_numerator(ctx, q, delta, b_mode, c_mode), z, k)
            for k in range(3)
        ]
        passed(
            f"same_source_mode_{mode}",
            all(abs(value - expected) < ctx.mpf("1e-80") for value, expected in zip(row, oracle)),
        )

    mode_sectors = swirl_moment_mode_numerator_delta_sectors_ordinary(
        z, target_row, d_b_row, d_c_row
    )
    passed("three_modal_delta_sector_numerators", len(mode_sectors) == 3)
    for mode, ((b_mode, c_mode), sectors) in enumerate(zip(mode_functions, mode_sectors)):
        finite_delta_numerator = [
            sectors[0][order] + delta * sectors[1][order] for order in range(3)
        ]
        oracle = [
            ctx.diff(lambda q: _full_swirl_numerator(ctx, q, delta, b_mode, c_mode), z, order)
            for order in range(3)
        ]
        passed(
            f"same_source_mode_delta_sector_{mode}",
            all(
                abs(value - expected) < ctx.mpf("1e-80")
                for value, expected in zip(finite_delta_numerator, oracle)
            ),
        )

    passed(
        "modal_target_rows_require_Z3",
        _raises(lambda: swirl_moment_modes_ordinary(z, delta, target_row[:3], d_b_row, d_c_row)),
    )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "maximum_supporting_difference": ctx.nstr(max(differences), 25),
        "ordinary_formula": (
            "Dsw=[-2*delta*z*b+(1-z^2)*b_Z+2*(1+delta)*z*c-"
            "(1-z^2)*c_Z]/[2*(1-delta*z^2)]"
        ),
        "api": {
            "drive": "swirl_drive_ordinary(z, delta, b, c)",
            "modes": "swirl_moment_modes_ordinary(z, delta, target_square, dB, dC)",
            "numerator_delta_sectors": "swirl_numerator_delta_sectors_ordinary(z, b, c)",
            "mode_numerator_delta_sectors": (
                "swirl_moment_mode_numerator_delta_sectors_ordinary("
                "z, target_square, dB, dC)"
            ),
            "dress": "dress_ordinary_jet(relative_amplitude, raw_moment)",
            "rows": "b,c,target_square,dB,dC are complete ordinary Z0..Z3 rows",
            "modes_return": "[base_numerator, j1_deltaC_numerator, j2_deltaB_numerator], each Z0..Z2",
        },
        "ancestor_constructor_producer_checker_fixture_run": False,
        "source_loading_or_adapter": False,
        "quadrature_or_integration": False,
        "temporal_recursion": False,
        "NS_residual": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
