"""Bounded independent checks for the frozen macro radial kernel.

The production kernel is compared with fresh scalar quadrature of
``exp(-(L-u))*B(u,Z)`` for nonlinear ordinary ``Z`` source functions.  A
separate finite-length toy check integrates an exact own-mean ODE and verifies
the endpoint identity through ordinary ``Z1`` differentiation.  No repository
source report, owner, producer, width materializer, or accepted numerical
fixture is loaded.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any, Callable, Sequence

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_macro_weighted_kernel import (  # noqa: E402
    apply_c_delta,
    evaluate_macro_modes,
    macro_log_length,
    radial_kernel,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (ArithmeticError, TypeError, ValueError):
        return True
    return False


def _row(ctx: Any, function: Callable[[Any], Any], z: Any) -> list[Any]:
    return [function(z), ctx.diff(function, z), ctx.diff(function, z, 2)]


def _box(ctx: Any, iv: Any, value: Any, radius: str = "1e-48") -> Any:
    width = ctx.mpf(radius)
    return iv.mpf([ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)])


def _bounds(value: Any) -> tuple[Any, Any]:
    if hasattr(value, "_mpi_"):
        return mp.make_mpf(value._mpi_[0]), mp.make_mpf(value._mpi_[1])
    return value, value


def _contains(value: Any, enclosure: Any) -> bool:
    lower, upper = _bounds(enclosure)
    return lower <= value <= upper


def _subset(inner: Any, outer: Any) -> bool:
    ilo, ihi = _bounds(inner)
    olo, ohi = _bounds(outer)
    return olo <= ilo <= ihi <= ohi


def _source_functions(ctx: Any) -> tuple[Callable[[Any], Any], ...]:
    # Each source function is nonlinear in Z.  Their signs are deliberately
    # mixed, including the negative-rate mode in the radial kernel.
    return (
        lambda z: ctx.mpf(".71") + ctx.mpf(".12") * ctx.sin(z) + ctx.mpf(".025") * z**2,
        lambda z: -ctx.mpf(".23") + ctx.mpf(".09") * ctx.cos(2 * z) - ctx.mpf(".04") * z,
        lambda z: ctx.mpf(".16") - ctx.mpf(".05") * ctx.sin(3 * z) + ctx.mpf(".015") * z**3,
    )


def _oracle(ctx: Any, functions: Sequence[Callable[[Any], Any]], rates: Sequence[int], L: Any) -> Callable[[Any], Any]:
    def value(z: Any) -> Any:
        result = ctx.mpf("0")
        for function, rate in zip(functions, rates):
            result += function(z) * ctx.quad(
                lambda u, rate=rate: ctx.exp(rate * (u - L)),
                [0, L],
            )
        return result

    return value


def _own_mean_functions(ctx: Any, delta: Any):
    """Return V, exact m, Qy, H and E for a finite-length toy source."""

    def V(y: Any, z: Any) -> Any:
        v0 = 1 + z / 5 + z**2 / 17
        v1 = ctx.mpf(".31") + z / 7 + ctx.sin(z) / 9
        return v0 + v1 * ctx.exp(-y / 3)

    def Vy(y: Any, z: Any) -> Any:
        return -(
            ctx.mpf(".31") + z / 7 + ctx.sin(z) / 9
        ) * ctx.exp(-y / 3) / 3

    def m(y: Any, z: Any) -> Any:
        # m_y = V - m, m(0)=1+z/11.  The formula is exact for the two
        # exponential modes in V and avoids a numerical ODE solve.
        v0 = 1 + z / 5 + z**2 / 17
        v1 = ctx.mpf(".31") + z / 7 + ctx.sin(z) / 9
        m0 = 1 + z / 11
        return (
            v0
            + v1 * ctx.mpf("1.5") * ctx.exp(-y / 3)
            + (m0 - v0 - v1 * ctx.mpf("1.5")) * ctx.exp(-y)
        )

    def Qy(y: Any, z: Any) -> Any:
        vy = Vy(y, z)
        vz = ctx.diff(lambda zz: V(y, zz), z)
        mz = ctx.diff(lambda zz: m(y, zz), z)
        denominator = 1 - delta * z * z
        return (
            2 * z * vy
            - (1 - delta) * z * (V(y, z) - m(y, z))
            - (1 - z * z) * (vz - mz)
        ) / denominator

    def H(y: Any, z: Any) -> Any:
        return ctx.exp(y) * Qy(y, z)

    def E(y: Any, z: Any) -> Any:
        return ctx.exp(y) * 2 * z * Vy(y, z) / (1 - delta * z * z)

    return V, m, Qy, H, E


def run_fixture() -> dict[str, Any]:
    ctx = mp.mp.clone()
    ctx.dps = 90
    iv = MPIntervalContext()
    iv.dps = 130
    controls: list[dict[str, Any]] = []

    def passed(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    functions = _source_functions(ctx)
    rates = (2, 0, -1)
    L = ctx.mpf(".83")
    Li_coarse = iv.mpf([ctx.nstr(L - ctx.mpf("1e-18"), ctx.dps), ctx.nstr(L + ctx.mpf("1e-18"), ctx.dps)])
    Li_fine = iv.mpf([ctx.nstr(L - ctx.mpf("1e-30"), ctx.dps), ctx.nstr(L + ctx.mpf("1e-30"), ctx.dps)])
    z = ctx.mpf(".17")
    rows = [_row(ctx, function, z) for function in functions]
    modes = [
        {"ordinary_rows": row, "rate": rate}
        for row, rate in zip(rows, rates)
    ]

    point_packet = evaluate_macro_modes(
        ctx,
        modes,
        L,
        source_proved=True,
        L_Z_independent=True,
    )
    oracle = _oracle(ctx, functions, rates, L)
    oracle_rows = [ctx.diff(oracle, z, order) for order in range(3)]
    for order, (actual, expected) in enumerate(zip(point_packet["response"], oracle_rows)):
        passed(f"point_oracle_Z{order}", abs(actual - expected) < ctx.mpf("1e-70"))

    # The C_delta row uses the same independent oracle and its first two
    # ordinary derivatives.  The source functions are differentiated by mp.diff
    # through the fresh quadrature, rather than through production rows.
    delta = ctx.mpf(".23")
    C_point = apply_c_delta(z, delta, oracle_rows)
    C_expected = [
        ((1 - delta) * z * oracle(z) + (1 - z * z) * ctx.diff(oracle, z))
        / (1 - delta * z * z),
        ctx.diff(
            lambda zz: (
                (1 - delta) * zz * oracle(zz)
                + (1 - zz * zz) * ctx.diff(oracle, zz)
            )
            / (1 - delta * zz * zz),
            z,
        ),
    ]
    for order, (actual, expected) in enumerate(zip(C_point, C_expected)):
        passed(f"C_delta_oracle_Z{order}", abs(actual - expected) < ctx.mpf("1e-66"))

    # Directed rows provide a uniform source enclosure for the radial angular
    # factor.  They are intentionally wider than the point rows and are never
    # reduced to a midpoint.
    truth_modes = [
        {
            "ordinary_rows": [_box(ctx, iv, value) for value in row],
            "rate": rate,
        }
        for row, rate in zip(rows, rates)
    ]
    coarse = evaluate_macro_modes(
        iv, truth_modes, Li_coarse, source_proved=True, L_Z_independent=True
    )
    fine = evaluate_macro_modes(
        iv, truth_modes, Li_fine, source_proved=True, L_Z_independent=True
    )
    for order, value in enumerate(point_packet["response"]):
        passed(f"directed_contains_Z{order}", _contains(value, coarse["response"][order]))
        passed(f"directed_fine_contains_Z{order}", _contains(value, fine["response"][order]))
        passed(f"directed_refinement_Z{order}", _subset(fine["response"][order], coarse["response"][order]))
    for index, kernel in enumerate(coarse["kernels"]):
        lower, upper = _bounds(kernel)
        passed(f"positive_kernel_{index}", lower > 0 and upper >= lower)

    # Exact zero length and the exp-1 fallback are both part of the scalar
    # kernel contract.  The fallback context deliberately has no expm1 method.
    class ExpOnly:
        @staticmethod
        def mpf(value: Any) -> Any:
            return ctx.mpf(value)

        @staticmethod
        def exp(value: Any) -> Any:
            return ctx.exp(value)

        @staticmethod
        def log(value: Any) -> Any:
            return ctx.log(value)

    passed("zero_length_rate0_exact", radial_kernel(ctx, 0, ctx.mpf("0")) == 0)
    passed("zero_length_negative_rate_exact", radial_kernel(ctx, -2, ctx.mpf("0")) == 0)
    passed(
        "native2_to_R100_actual_length",
        abs(macro_log_length(ctx, ctx.mpf("4")) - ctx.log(ctx.mpf("25")))
        < ctx.mpf("1e-80"),
    )
    passed(
        "expm1_fallback",
        abs(radial_kernel(ExpOnly, -1, L) - radial_kernel(ctx, -1, L)) < ctx.mpf("1e-80"),
    )

    # A finite-length own-mean toy supplies an independent endpoint identity.
    # With H=R*Qy and E=R*2Z*Vy/L_delta, the source is independently defined
    # as R*C_delta(Vy).  The own-mean equation then implies E_y-H_y=source;
    # both that local identity and its endpoint FTC are checked through
    # ordinary Z0/Z1 derivatives at one nonlinear Z point.
    delta_ode = ctx.mpf(".19")
    V, mean, Qy, H, E = _own_mean_functions(ctx, delta_ode)
    lo = ctx.mpf(".11")
    hi = ctx.mpf(".83")
    passed(
        "own_mean_ode_exact",
        abs(
            ctx.diff(lambda yy: mean(yy, z), lo)
            - ((1 + z / 5 + z**2 / 17) + (ctx.mpf(".31") + z / 7 + ctx.sin(z) / 9) * ctx.exp(-lo / 3) - mean(lo, z))
        )
        < ctx.mpf("1e-70"),
    )

    def Vy(yy: Any, zz: Any) -> Any:
        return ctx.diff(lambda tt: V(tt, zz), yy)

    def c_source(yy: Any, zz: Any) -> Any:
        vy = Vy(yy, zz)
        vy_z = ctx.diff(lambda zz2: Vy(yy, zz2), zz)
        return ctx.exp(yy) * (
            (1 - delta_ode) * zz * vy + (1 - zz * zz) * vy_z
        ) / (1 - delta_ode * zz * zz)

    def endpoint_source(yy: Any, zz: Any) -> Any:
        return ctx.diff(lambda tt: E(tt, zz), yy) - ctx.diff(
            lambda tt: H(tt, zz), yy
        )

    interior_y = ctx.mpf(".37")
    for order in (0, 1):
        local_left = ctx.diff(lambda zz: endpoint_source(interior_y, zz), z, order)
        local_right = ctx.diff(lambda zz: c_source(interior_y, zz), z, order)
        passed(
            f"own_mean_local_Cdelta_source_Z{order}",
            abs(local_left - local_right) < ctx.mpf("1e-62"),
        )

    def endpoint_residual(zz: Any) -> Any:
        return H(lo, zz) - H(hi, zz) + E(hi, zz) - E(lo, zz) - ctx.quad(
            lambda yy: c_source(yy, zz), [lo, hi]
        )

    for order in (0, 1):
        passed(
            f"own_mean_endpoint_identity_Z{order}",
            abs(ctx.diff(endpoint_residual, z, order)) < ctx.mpf("1e-62"),
        )
    source_row = [c_source(lo, z), ctx.diff(lambda zz: c_source(lo, zz), z)]
    passed("own_mean_source_Z1_finite", all(ctx.isfinite(value) for value in source_row))

    # Fail closed on missing ordinary Z2, missing proof, and a Z-dependent
    # length declaration.
    passed(
        "missing_Z2_rejected",
        _raises(
            lambda: evaluate_macro_modes(
                ctx,
                [{"ordinary_rows": rows[0][:2], "rate": rates[0]}],
                L,
                source_proved=True,
                L_Z_independent=True,
            )
        ),
    )
    passed(
        "missing_source_proof_rejected",
        _raises(lambda: evaluate_macro_modes(ctx, modes, L, source_proved=False)),
    )
    passed(
        "missing_length_independence_rejected",
        _raises(
            lambda: evaluate_macro_modes(
                ctx, modes, L, source_proved=True, L_Z_independent=False
            )
        ),
    )
    passed(
        "noninteger_rate_rejected",
        _raises(
            lambda: evaluate_macro_modes(
                ctx,
                [{"ordinary_rows": rows[0], "rate": ctx.mpf(".5")}],
                L,
                source_proved=True,
            )
        ),
    )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "independent_oracle": "fresh ctx.quad of exp(-(L-u))*B(u,Z), with mp.diff through ordinary Z2",
        "own_mean_oracle": "exact m_y=V-m toy and direct endpoint source integral through ordinary Z1",
        "source_loaders_or_owners_called": False,
        "quadrature_in_production": False,
        "width_materialized": False,
        "physical_R100_epsilon_hbB_factors_attached": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
