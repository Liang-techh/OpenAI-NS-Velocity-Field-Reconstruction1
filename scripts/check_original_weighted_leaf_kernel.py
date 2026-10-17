"""Independent checks for the pure factored weighted-leaf row algebra.

The oracle below integrates the finite, nonzero-width source directly.  Its
inner angular integral uses the exact ``expm1`` primitive, and each outer K
integral uses fresh Gauss-Legendre quadrature.  No accepted source producer,
saved moment packet, width logarithm, or ancestor checker is loaded.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_weighted_leaf_kernel import (  # noqa: E402
    build_factored_leaf_records,
    c_delta,
    combine_coefficients,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (ArithmeticError, TypeError, ValueError):
        return True
    return False


def _ordinary_product(left: list[Any], right: list[Any]) -> list[Any]:
    """Independent ordinary product through Z2, including the factor two."""

    return [
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
    ]


def _angular_inner(ctx: Any, functions: tuple[Callable[[Any], Any], ...], z: Any, t: Any, h: Any) -> Any:
    value = ctx.mpf("0")
    for index, function in enumerate(functions):
        rate = (1 - index) * h
        primitive = t if rate == 0 else ctx.expm1(rate * t) / rate
        value += function(z) * primitive
    return value


def _finite_source(
    ctx: Any,
    phase: Any,
    h: Any,
    q: int,
    angular: tuple[Callable[[Any], Any], ...],
    axial: tuple[Callable[[Any], Any], ...],
) -> Callable[[Any], Any]:
    """Return the direct finite-h weighted primitive J(z)."""

    def kappa(t: Any) -> Any:
        return 1 - t / 4

    def mass() -> Any:
        return phase - phase**2 / 8

    def kernel(z: Any, index: int) -> Any:
        rate = (q - index + 1) * h

        def integrand(t: Any) -> Any:
            ell = -h**2 / 2 * _angular_inner(ctx, angular, z, t, h)
            return kappa(t) * ctx.exp(rate * t + ell)

        return ctx.quadgl(integrand, [0, phase])

    def primitive(z: Any) -> Any:
        return sum(
            (function(z) * kernel(z, index) for index, function in enumerate(axial)),
            ctx.mpf("0"),
        )

    primitive.mass = mass  # type: ignore[attr-defined]
    return primitive


def _source_rows(ctx: Any, function: Callable[[Any], Any], z: Any) -> list[Any]:
    return [function(z), ctx.diff(function, z), ctx.diff(function, z, 2)]


def _evaluate_terms(ctx: Any, result: dict[str, Any], h: Any, order: int) -> Any:
    combined = result["combined_factored_jet"]
    scale = result["raw46"]["physical_scale"]
    terms = combined["terms"]
    return scale * h**2 * sum(
        h**term["width_power_inside_factor"] * term["coefficient_rows"][order]
        for term in terms
    )


def _box(ctx: Any, interval_ctx: Any, value: Any, radius: str = "1e-55") -> Any:
    width = ctx.mpf(radius)
    return interval_ctx.mpf(
        [ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)]
    )


def _inside(value: Any, enclosure: Any) -> bool:
    return value in enclosure


def run_fixture() -> dict[str, Any]:
    ctx = mp.mp.clone()
    ctx.dps = 70
    interval_ctx = MPIntervalContext()
    interval_ctx.dps = 110
    controls: list[dict[str, Any]] = []
    differences: list[Any] = []

    def passed(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    # Nonconstant signed source rows make both ordinary C0 and C1 visible.
    angular = (
        lambda z: ctx.mpf(".71") + ctx.mpf(".09") * z + ctx.sin(z / 3) / 8,
        lambda z: ctx.mpf("-.18") + ctx.mpf(".07") * z**2 + ctx.cos(z / 4) / 9,
        lambda z: ctx.mpf(".27") - ctx.mpf(".05") * z + ctx.sin(z / 5) / 11,
    )
    axial = (
        lambda z: ctx.mpf(".83") + ctx.mpf(".12") * z + ctx.sin(z / 6) / 10,
        lambda z: ctx.mpf("-.41") + ctx.mpf(".08") * z**2 + ctx.cos(z / 5) / 12,
        lambda z: ctx.mpf(".29") - ctx.mpf(".06") * z + ctx.sin(z / 7) / 13,
    )
    z = ctx.mpf(".23")
    delta = ctx.mpf(".031")
    q = 3
    R100 = ctx.mpf("100")
    epsilon = ctx.mpf(".17")

    for phase_text in (".37", ".71"):
        phase = ctx.mpf(phase_text)
        mass = phase - phase**2 / 8
        for h_text in (".002", ".001"):
            h = ctx.mpf(h_text)
            direct_j = _finite_source(ctx, phase, h, q, angular, axial)
            bsum = lambda x: sum((function(x) for function in axial), ctx.mpf("0"))
            remainder = lambda x: (direct_j(x) - bsum(x) * mass) / h
            b_rows = _source_rows(ctx, bsum, z)
            r_rows = _source_rows(ctx, remainder, z)
            result = combine_coefficients(
                z,
                delta,
                b_rows,
                r_rows,
                mass,
                R100=R100,
                epsilon=epsilon,
                width_symbol="h",
                source_proved=True,
                shared_source_proved=True,
                mass_Z_independent=True,
            )
            passed(
                f"{phase_text}_{h_text}_shape",
                result["raw46"]["source_ordinary_orders"] == [0, 1, 2]
                and result["raw47"]["source_ordinary_orders"] == [0, 1, 2]
                and result["raw46"]["ordinary_coefficient_orders"] == [0, 1]
                and result["raw47"]["ordinary_coefficient_orders"] == [0, 1]
                and result["width_materialized"] is False,
            )

            def c_of(function: Callable[[Any], Any], x: Any) -> Any:
                return (
                    (1 - delta) * x * function(x)
                    + (1 - x * x) * ctx.diff(function, x)
                ) / (1 - delta * x * x)

            # This is the unsplit physical primitive evaluated directly from
            # J.  It is independent from the leaf algebra and the R rows.
            physical = lambda x: R100 * epsilon * h**2 * c_of(direct_j, x)
            for order in (0, 1):
                factored = _evaluate_terms(ctx, result, h, order)
                direct = ctx.diff(physical, z, order)
                difference = abs(factored - direct)
                differences.append(difference)
                passed(
                    f"{phase_text}_{h_text}_physical_ordinary_Z{order}",
                    difference < ctx.mpf("1e-43"),
                )

            expected_main = [mass * value for value in c_delta(z, delta, b_rows)]
            passed(
                f"{phase_text}_{h_text}_raw46_main_factor",
                all(
                    abs(got - R100 * epsilon * expected) < ctx.mpf("1e-55")
                    for got, expected in zip(
                        result["raw46"]["coefficient_rows"], expected_main
                    )
                ),
            )
            passed(
                f"{phase_text}_{h_text}_raw47_signed_factor",
                all(
                    abs(got - R100 * epsilon * expected) < ctx.mpf("1e-55")
                    for got, expected in zip(
                        result["raw47"]["coefficient_rows"],
                        c_delta(z, delta, r_rows),
                    )
                ),
            )

    # The endpoint is an exact empty integration interval.  The source rows
    # remain nonzero and are still required by the API; only the independently
    # computed M0 and R vanish.
    phase = ctx.mpf("0")
    h = ctx.mpf(".001")
    endpoint_mass = ctx.mpf("0")
    endpoint_j = _finite_source(ctx, phase, h, q, angular, axial)
    endpoint_bsum = lambda x: sum((function(x) for function in axial), ctx.mpf("0"))
    endpoint_remainder = lambda x: (endpoint_j(x) - endpoint_bsum(x) * endpoint_mass) / h
    endpoint = combine_coefficients(
        z,
        delta,
        _source_rows(ctx, endpoint_bsum, z),
        _source_rows(ctx, endpoint_remainder, z),
        endpoint_mass,
        R100=R100,
        epsilon=epsilon,
        source_proved=True,
        shared_source_proved=True,
        mass_Z_independent=True,
    )
    passed(
        "phase_zero_exact",
        _evaluate_terms(ctx, endpoint, h, order=0) == 0,
    )
    passed(
        "phase_zero_source_not_replaced",
        any(abs(value) > ctx.mpf("1e-8") for value in endpoint["raw46"]["ordinary_source_rows"]),
    )

    # Independent product identity through Z2 checks the ordinary derivative
    # factor two retained by any caller assembling modal products.
    left = [ctx.mpf("1.2"), ctx.mpf("-.4"), ctx.mpf(".7")]
    right = [ctx.mpf("-.3"), ctx.mpf(".8"), ctx.mpf("-.2")]
    product = lambda x: (ctx.mpf("1.2") + ctx.mpf("-.4") * (x - z) + ctx.mpf(".7") * (x - z) ** 2 / 2) * (
        ctx.mpf("-.3") + ctx.mpf(".8") * (x - z) + ctx.mpf("-.2") * (x - z) ** 2 / 2
    )
    direct_product = [ctx.diff(product, z, order) for order in range(3)]
    passed(
        "ordinary_product_Z2_factor_two",
        all(abs(a - b) < ctx.mpf("1e-60") for a, b in zip(_ordinary_product(left, right), direct_product)),
    )

    # Re-run one nonzero finite-width row with narrow directed input boxes.
    # Every point coefficient must remain inside the independently propagated
    # interval result, including the two factored terms.
    iv_b = [_box(ctx, interval_ctx, value) for value in b_rows]
    iv_r = [_box(ctx, interval_ctx, value) for value in r_rows]
    iv_result = combine_coefficients(
        _box(ctx, interval_ctx, z),
        _box(ctx, interval_ctx, delta),
        iv_b,
        iv_r,
        _box(ctx, interval_ctx, mass),
        R100=_box(ctx, interval_ctx, R100),
        epsilon=_box(ctx, interval_ctx, epsilon),
        source_proved=True,
        shared_source_proved=True,
        mass_Z_independent=True,
    )
    for leaf_name in ("raw46", "raw47"):
        for field in ("ordinary_C0_C1", "coefficient_rows", "unscaled_coefficient_rows"):
            for order, point in enumerate(result[leaf_name][field]):
                passed(
                    f"directed_{leaf_name}_{field}_Z{order}",
                    _inside(point, iv_result[leaf_name][field][order]),
                )
    for term_index, term in enumerate(result["combined_factored_jet"]["terms"]):
        interval_term = iv_result["combined_factored_jet"]["terms"][term_index]
        for order, point in enumerate(term["coefficient_rows"]):
            passed(
                f"directed_combined_term{term_index}_Z{order}",
                _inside(point, interval_term["coefficient_rows"][order]),
            )
    passed(
        "directed_denominator_crossing_rejected",
        _raises(
            lambda: combine_coefficients(
                interval_ctx.mpf([".099", ".101"]),
                interval_ctx.mpf(["99", "101"]),
                iv_b,
                iv_r,
                _box(ctx, interval_ctx, mass),
                source_proved=True,
                shared_source_proved=True,
                mass_Z_independent=True,
            )
        ),
    )
    passed(
        "directed_negative_mass_rejected",
        _raises(
            lambda: combine_coefficients(
                _box(ctx, interval_ctx, z),
                _box(ctx, interval_ctx, delta),
                iv_b,
                iv_r,
                interval_ctx.mpf(["-.1", ".1"]),
                source_proved=True,
                shared_source_proved=True,
                mass_Z_independent=True,
            )
        ),
    )

    good_rows = [ctx.mpf("1"), ctx.mpf(".2"), ctx.mpf("-.1")]
    bad_common = dict(source_proved=True, mass_Z_independent=True)
    passed(
        "missing_source_proof_rejected",
        _raises(lambda: combine_coefficients(z, delta, good_rows, good_rows, ctx.mpf(".1"), mass_Z_independent=True)),
    )
    passed(
        "missing_mass_independence_proof_rejected",
        _raises(lambda: combine_coefficients(z, delta, good_rows, good_rows, ctx.mpf(".1"), source_proved=True)),
    )
    passed(
        "shared_source_proof_rejected",
        _raises(
            lambda: combine_coefficients(
                z,
                delta,
                good_rows,
                good_rows,
                ctx.mpf(".1"),
                source_proved=True,
                shared_source_proved=False,
                mass_Z_independent=True,
            )
        ),
    )
    passed(
        "missing_Z2_rejected",
        _raises(
            lambda: combine_coefficients(
                z,
                delta,
                good_rows[:2],
                good_rows,
                ctx.mpf(".1"),
                **bad_common,
            )
        ),
    )
    passed(
        "malformed_source_rejected",
        _raises(
            lambda: combine_coefficients(
                z,
                delta,
                [good_rows, [ctx.mpf("1"), None, ctx.mpf("0")]],
                good_rows,
                ctx.mpf(".1"),
                **bad_common,
            )
        ),
    )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "maximum_physical_difference": mp.nstr(max(differences), 25),
        "independent_oracle": "finite nonzero h direct Bj*K_a outer Gauss-Legendre quadrature with exact expm1 angular primitive",
        "source_producer_or_ancestor_checker_called": False,
        "width_selected_or_materialized": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
