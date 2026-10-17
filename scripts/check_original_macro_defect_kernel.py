"""Independent checks for the native/macro endpoint-defect kernel.

The production function is pure ordinary-row arithmetic.  This checker uses
fresh nonlinear ``Z`` functions and differentiates the full ``Q_y`` formula
independently.  A separate finite-``y`` toy solves ``m_y=V-m`` exactly and
uses fresh quadrature only for its endpoint identity.  No repository report,
owner, producer, constructor, or saved fixture is loaded.
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

from openai_ns_reconstruction.original_macro_defect_kernel import (  # noqa: E402
    evaluate_seam_defect,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (ArithmeticError, TypeError, ValueError):
        return True
    return False


def _ordinary_rows(
    ctx: Any, function: Callable[[Any], Any], z: Any, order: int = 1
) -> list[Any]:
    return [ctx.diff(function, z, derivative_order) for derivative_order in range(order + 1)]


def _box(ctx: Any, interval_ctx: Any, value: Any, radius: str = "1e-45") -> Any:
    width = ctx.mpf(radius)
    return interval_ctx.mpf(
        [ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)]
    )


def _bounds(value: Any) -> tuple[Any, Any]:
    if hasattr(value, "_mpi_"):
        lower, upper = value._mpi_
        return mp.make_mpf(lower), mp.make_mpf(upper)
    return value, value


def _contains(point: Any, enclosure: Any) -> bool:
    lower, upper = _bounds(enclosure)
    return lower <= point <= upper


def _full_qy(
    ctx: Any,
    z: Any,
    delta: Any,
    U: Callable[[Any], Any],
    mean: Callable[[Any], Any],
    Uy: Callable[[Any], Any],
) -> Any:
    """Independent full ``Q_y`` expression, including value terms."""

    denominator = 1 - delta * z * z
    return (
        2 * z * Uy(z)
        - (1 - delta) * z * (U(z) - mean(z))
        - (1 - z * z) * (ctx.diff(U, z) - ctx.diff(mean, z))
    ) / denominator


def _own_mean_toy(ctx: Any):
    """Return nonlinear ``V,m`` with an exact solution of ``m_y=V-m``."""

    def a(z: Any) -> Any:
        return 1 + z / 5 + z * z / 17 + ctx.sin(z) / 13

    def b(z: Any) -> Any:
        return ctx.mpf(".31") + z / 7 + ctx.sin(z) / 9

    def m0(z: Any) -> Any:
        return 1 + z / 11 + ctx.sin(z) / 19

    def V(y: Any, z: Any) -> Any:
        return a(z) + b(z) * ctx.exp(-y / 3)

    def Vy(y: Any, z: Any) -> Any:
        return -b(z) * ctx.exp(-y / 3) / 3

    def mean(y: Any, z: Any) -> Any:
        return (
            a(z)
            + ctx.mpf("1.5") * b(z) * ctx.exp(-y / 3)
            + (m0(z) - a(z) - ctx.mpf("1.5") * b(z)) * ctx.exp(-y)
        )

    return V, Vy, mean


def run_fixture() -> dict[str, Any]:
    """Run only new bounded arithmetic and independent toy controls."""

    ctx = mp.mp.clone()
    ctx.dps = 100
    interval_ctx = MPIntervalContext()
    interval_ctx.dps = 140
    controls: list[dict[str, Any]] = []

    def passed(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    delta = ctx.mpf(".23")
    epsilon = ctx.mpf(".37")
    hbB = ctx.mpf(".29")
    epsilon_hbB = epsilon * hbB

    # Nonlinear source rows exercise both the signed value and derivative
    # before the endpoint operator is applied.
    def native_B_function(z: Any) -> Any:
        return ctx.mpf(".71") + ctx.sin(z / 3) + ctx.exp(z / 7) + z * z / 11

    def macro_B_function(z: Any) -> Any:
        return -ctx.mpf(".23") + ctx.cos(z / 5) + z**3 / 17

    # Independent full-Qy values: native and macro share U,m, while their Uy
    # rows differ.  The normalized rows are Uy/(epsilon*hbB).
    def common_U(z: Any) -> Any:
        return 1 + ctx.sin(z / 4) + z * z / 19 + ctx.exp(-z / 6)

    def common_m(z: Any) -> Any:
        return ctx.mpf(".4") + ctx.cos(z / 5) + z / 13 + z**3 / 23

    def native_Uy(z: Any) -> Any:
        return epsilon_hbB * native_B_function(z)

    def macro_Uy(z: Any) -> Any:
        return epsilon_hbB * macro_B_function(z)

    for z_text in ("-.31", ".07", ".29"):
        z = ctx.mpf(z_text)
        native_rows = _ordinary_rows(ctx, lambda q: native_Uy(q) / epsilon_hbB, z)
        macro_rows = _ordinary_rows(ctx, lambda q: macro_Uy(q) / epsilon_hbB, z)
        result = evaluate_seam_defect(
            ctx,
            z,
            delta,
            native_rows,
            macro_rows,
            values_source_proved=True,
            same_radius_proved=True,
        )
        full_difference = lambda q: _full_qy(
            ctx, q, delta, common_U, common_m, native_Uy
        ) - _full_qy(ctx, q, delta, common_U, common_m, macro_Uy)
        independent_rows = [ctx.diff(full_difference, z, order) for order in range(2)]
        for order, (actual, expected) in enumerate(
            zip(independent_rows, [epsilon_hbB * value for value in result["defect_C0_C1"]])
        ):
            passed(
                f"full_Qy_difference_matches_Edefect_{z_text}_Z{order}",
                abs(actual - expected) < ctx.mpf("1e-75"),
            )

        passed(
            f"signed_source_rows_{z_text}",
            all(
                abs(actual - expected) < ctx.mpf("1e-80")
                for actual, expected in zip(
                    result["source_difference_rows"],
                    [
                        native_B_function(z) - macro_B_function(z),
                        ctx.diff(native_B_function, z) - ctx.diff(macro_B_function, z),
                    ],
                )
            ),
        )
        certificate = result["certificate"]
        passed(
            f"certificate_{z_text}",
            certificate["passed"] is True
            and certificate["source_Z_orders"] == [0, 1]
            and certificate["values_source_proved"] is True
            and certificate["same_radius_proved"] is True
            and certificate["widths_materialized"] is False,
        )

    # A nonzero physical scale stays outside the dimensionless operator.  The
    # external rho/4 role is also kept outside, as required by the same-radius
    # coupling owned by the caller.
    z = ctx.mpf(".17")
    native_rows = _ordinary_rows(ctx, lambda q: native_Uy(q) / epsilon_hbB, z)
    macro_rows = _ordinary_rows(ctx, lambda q: macro_Uy(q) / epsilon_hbB, z)
    result = evaluate_seam_defect(
        ctx,
        z,
        delta,
        native_rows,
        macro_rows,
        values_source_proved=True,
        same_radius_proved=True,
    )
    actual_R0 = ctx.mpf("2.75")
    external_rho = ctx.mpf("3.6")
    physical = actual_R0 * (external_rho / 4) * epsilon_hbB
    expected_physical = [physical * value for value in result["defect_C0_C1"]]
    independent_physical = [
        actual_R0
        * (external_rho / 4)
        * value
        for value in (
            ctx.diff(
                lambda q: _full_qy(ctx, q, delta, common_U, common_m, native_Uy)
                - _full_qy(ctx, q, delta, common_U, common_m, macro_Uy),
                z,
                order,
            )
            for order in (0, 1)
        )
    ]
    passed(
        "nonzero_R_epsilon_hbB_rho_over_4_scale",
        all(
            abs(actual - expected) < ctx.mpf("1e-75")
            for actual, expected in zip(independent_physical, expected_physical)
        )
        and actual_R0 != 0
        and epsilon_hbB != 0
        and external_rho / 4 != 0,
    )

    # Perturbing common U/m values adds the mean-value jump to full Qy.  The
    # resulting mismatch is evidence that the seam proof is missing; it is not
    # treated as an operator failure.
    def perturbed_macro_m(z_value: Any) -> Any:
        return common_m(z_value) + ctx.mpf(".17") + z_value / 9

    perturbed_difference = lambda q: _full_qy(
        ctx, q, delta, common_U, common_m, native_Uy
    ) - _full_qy(ctx, q, delta, common_U, perturbed_macro_m, macro_Uy)
    perturbed_rows = [ctx.diff(perturbed_difference, z, order) for order in range(2)]
    passed(
        "common_value_perturbation_adds_mean_jump",
        max(
            abs(actual - expected)
            for actual, expected in zip(
                perturbed_rows,
                [epsilon_hbB * value for value in result["defect_C0_C1"]],
            )
        )
        > ctx.mpf("1e-8"),
    )

    # Directed input rows must enclose the point result.  No interval midpoint
    # is selected by the production function.
    z_iv = _box(ctx, interval_ctx, z, "1e-35")
    delta_iv = _box(ctx, interval_ctx, delta, "1e-35")
    native_iv = [_box(ctx, interval_ctx, value) for value in native_rows]
    macro_iv = [_box(ctx, interval_ctx, value) for value in macro_rows]
    directed = evaluate_seam_defect(
        interval_ctx,
        z_iv,
        delta_iv,
        native_iv,
        macro_iv,
        values_source_proved=True,
        same_radius_proved=True,
    )
    passed(
        "directed_interval_contains_point",
        all(
            _contains(value, enclosure)
            for value, enclosure in zip(result["defect_C0_C1"], directed["defect_C0_C1"])
        ),
    )

    # A finite y toy independently checks the own-mean endpoint identity.  The
    # native and macro U,m values agree, but their Uy endpoint rows have a
    # deliberate offset.  The source is written directly as R*C_delta(Vy),
    # never defined as E_y-H_y.
    delta_ode = ctx.mpf(".19")
    V, Vy, mean = _own_mean_toy(ctx)
    lower_y = ctx.mpf(".11")
    upper_y = ctx.mpf(".83")
    seam_z = ctx.mpf(".17")

    def native_uy_y(y: Any, q: Any) -> Any:
        return Vy(y, q)

    def macro_uy_y(y: Any, q: Any) -> Any:
        return Vy(y, q) + ctx.mpf(".07") + ctx.sin(q) / 13

    def qy_y(y: Any, q: Any, uy: Callable[[Any, Any], Any]) -> Any:
        return _full_qy(
            ctx,
            q,
            delta_ode,
            lambda zz: V(y, zz),
            lambda zz: mean(y, zz),
            lambda zz: uy(y, zz),
        )

    def E_y(y: Any, q: Any, uy: Callable[[Any, Any], Any]) -> Any:
        return ctx.exp(y) * 2 * q * uy(y, q) / (1 - delta_ode * q * q)

    def H_y(y: Any, q: Any, uy: Callable[[Any, Any], Any]) -> Any:
        return ctx.exp(y) * qy_y(y, q, uy)

    def source_Cdelta_Vy(y: Any, q: Any) -> Any:
        denominator = 1 - delta_ode * q * q
        value = Vy(y, q)
        value_z = ctx.diff(lambda qq: Vy(y, qq), q)
        C_delta = ((1 - delta_ode) * q * value + (1 - q * q) * value_z) / denominator
        return ctx.exp(y) * C_delta

    def endpoint_residual(q: Any) -> Any:
        return (
            H_y(lower_y, q, native_uy_y)
            - H_y(upper_y, q, macro_uy_y)
            + E_y(upper_y, q, macro_uy_y)
            - E_y(lower_y, q, native_uy_y)
            - ctx.quad(lambda yy: source_Cdelta_Vy(yy, q), [lower_y, upper_y])
        )

    ode_error = ctx.diff(lambda yy: mean(yy, seam_z), lower_y) - (
        V(lower_y, seam_z) - mean(lower_y, seam_z)
    )
    passed("own_mean_ode_exact", abs(ode_error) < ctx.mpf("1e-75"))
    local_y = ctx.mpf(".37")
    independently_defined_source = source_Cdelta_Vy(local_y, seam_z)
    derivative_source = ctx.diff(
        lambda yy: E_y(yy, seam_z, macro_uy_y) - H_y(yy, seam_z, macro_uy_y),
        local_y,
    )
    passed(
        "source_is_independent_Cdelta_not_definition",
        abs(independently_defined_source - derivative_source) < ctx.mpf("1e-65"),
    )
    for order in (0, 1):
        passed(
            f"own_mean_endpoint_identity_Z{order}",
            abs(ctx.diff(endpoint_residual, seam_z, order)) < ctx.mpf("1e-55"),
        )

    # Fail closed on all seam prerequisites and dangerous/insufficient rows.
    passed(
        "missing_native_Bz_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows[:1],
                macro_rows,
                values_source_proved=True,
                same_radius_proved=True,
            )
        ),
    )
    passed(
        "missing_macro_Bz_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows,
                macro_rows[:1],
                values_source_proved=True,
                same_radius_proved=True,
            )
        ),
    )
    passed(
        "missing_values_source_proof_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows,
                macro_rows,
                values_source_proved=False,
                same_radius_proved=True,
            )
        ),
    )
    passed(
        "nonboolean_values_source_proof_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows,
                macro_rows,
                values_source_proved=1,
                same_radius_proved=True,
            )
        ),
    )
    passed(
        "missing_same_radius_proof_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows,
                macro_rows,
                values_source_proved=True,
                same_radius_proved=False,
            )
        ),
    )
    passed(
        "nonboolean_same_radius_proof_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                ctx,
                z,
                delta,
                native_rows,
                macro_rows,
                values_source_proved=True,
                same_radius_proved=ctx.mpf(1),
            )
        ),
    )
    passed(
        "denominator_crossing_rejected",
        _raises(
            lambda: evaluate_seam_defect(
                interval_ctx,
                interval_ctx.mpf(["-2", "2"]),
                interval_ctx.mpf(".5"),
                native_iv,
                macro_iv,
                values_source_proved=True,
                same_radius_proved=True,
            )
        ),
    )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "ordinary_formula": "E(D)=2*Z*D/(1-delta*Z^2), D=native_B-macro_B",
        "full_Qy_formula": "[2*Z*Uy-(1-delta)*Z*(U-m)-(1-Z^2)*(Uz-mz)]/(1-delta*Z^2)",
        "own_mean_source": "R*C_delta(Vy), independently defined and quadrature checked",
        "production_quadrature": False,
        "checker_quadrature": True,
        "repository_json_loaded": False,
        "ancestor_constructor_producer_fixture_replayed": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
