"""Independent finite-width checks for the new weighted remainder kernel.

The actual primitive in this checker uses an independent exact inner
``expm1`` expression for ``ell`` and fresh outer Gauss-Legendre quadrature for
each K integral; it is not reconstructed from the kernel's coefficient
formulas.  The two positive test weights have closed moments, while the
analytic angular and axial rows are nonconstant in Z.  ``mp.diff``
differentiates the independently integrated finite-h primitive through
ordinary Z order two.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_weighted_remainder_kernel import (  # noqa: E402
    enclose_weighted_remainder,
    signed_remainder_coefficient,
    signed_remainder_coefficients,
)


def _rows(ctx, function):
    return [function(ctx.mpf(".17")), ctx.diff(function, ctx.mpf(".17")), ctx.diff(function, ctx.mpf(".17"), 2)]


def _interval_row(iv, row):
    return [iv.mpf(value) for value in row]


def _bounds(ctx, value):
    lower, upper = value._mpi_
    return ctx.make_mpf(lower), ctx.make_mpf(upper)


def _fixtures(ctx):
    angular_functions = (
        lambda z: ctx.mpf(".71") + ctx.mpf(".12") * ctx.sin(z) + ctx.mpf(".025") * z**2,
        lambda z: ctx.mpf(".23") - ctx.mpf(".09") * z + ctx.mpf(".04") * ctx.cos(z),
        lambda z: -ctx.mpf(".16") + ctx.mpf(".05") * ctx.sin(2 * z),
    )
    axial_functions = (
        lambda z: ctx.mpf("1.1") + ctx.mpf(".18") * ctx.sin(z),
        lambda z: ctx.mpf(".37") + ctx.mpf(".08") * z**2,
        lambda z: -ctx.mpf(".24") + ctx.mpf(".07") * ctx.cos(z),
    )
    rates = (ctx.mpf(".62"), -ctx.mpf(".31"), ctx.mpf(".44"))
    return angular_functions, axial_functions, rates


def _moment_formula(ctx, phase, kind, power):
    if kind == "one":
        return phase ** (power + 1) / (power + 1)
    # kappa(t)=1-t
    return phase ** (power + 1) / (power + 1) - phase ** (power + 2) / (power + 2)


def _kappa(kind):
    return (lambda t: 1) if kind == "one" else (lambda t: 1 - t)


def _actual_remainder(ctx, phase, h, kind, angular_functions, axial_functions, rates):
    kappa = _kappa(kind)
    # Fixed Gauss-Legendre quadrature keeps the independent nested oracle
    # quick enough for all ordinary Z derivatives.  The integrands are smooth
    # on the compact toy phases and this is deliberately separate from the
    # interval and coefficient arithmetic in the production kernel.
    # M0 is supplied by the explicit fixture identity, never recomputed by a
    # legacy or adaptive mass producer.
    mass = phase if kind == "one" else phase - phase**2 / 2

    def angular_sum(t, z):
        return sum(
            (
                angular_functions[index](z) * ctx.exp((1 - index) * h * t)
                for index in range(len(angular_functions))
            ),
            ctx.mpf(0),
        )

    def local_ell(t, z):
        # The inner integral is evaluated independently by its exact
        # expm1/rate primitive.  The outer K integral remains fresh
        # Gauss-Legendre quadrature, and mp.diff still differentiates the full
        # finite-h primitive through Z2.
        integral = ctx.mpf(0)
        for index, function in enumerate(angular_functions):
            rate = (1 - index) * h
            integral += function(z) * (t if rate == 0 else ctx.expm1(rate * t) / rate)
        return -h**2 / 2 * integral

    def primitive(z):
        kernels = [
            ctx.quadgl(
                lambda t, rate=rate: kappa(t)
                * ctx.exp(rate * h * t + local_ell(t, z)),
                [0, phase],
            )
            for rate in rates
        ]
        return (
            sum((function(z) * kernel for function, kernel in zip(axial_functions, kernels)), ctx.mpf(0))
            - sum((function(z) * mass for function in axial_functions), ctx.mpf(0))
        ) / h

    return primitive


def run_fixture() -> dict[str, object]:
    ctx = mp.mp.clone()
    ctx.dps = 65
    iv = MPIntervalContext()
    iv.dps = 110
    controls: list[dict[str, object]] = []

    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    angular_functions, axial_functions, rates = _fixtures(ctx)
    z = ctx.mpf(".17")
    angular_rows = [_rows(ctx, function) for function in angular_functions]
    axial_rows = [_rows(ctx, function) for function in axial_functions]
    axial_modes = [
        {"ordinary_rows": row, "rate": rate}
        for row, rate in zip(axial_rows, rates)
    ]
    check("source_rows_are_nonconstant", any(abs(row[1]) > ctx.mpf("1e-8") for row in angular_rows + axial_rows))

    for kind in ("one", "one_minus_t"):
        for phase_text in ("0", ".35", ".7"):
            phase = ctx.mpf(phase_text)
            mass = _moment_formula(ctx, phase, kind, 0)
            W1 = _moment_formula(ctx, phase, kind, 1)
            W2 = _moment_formula(ctx, phase, kind, 2)
            leading = signed_remainder_coefficient(
                axial_modes,
                angular_rows,
                mass,
                W1,
                source_proved=True,
            )
            packet = signed_remainder_coefficients(
                axial_modes,
                angular_rows,
                mass,
                W1,
                W2,
                True,
            )
            check(f"{kind}_{phase_text}_leading_matches_packet", leading == packet["leading"])

            # The finite-width toy uses eta=.01 and a modest cap.  The h values
            # are all positive and satisfy the source smallness proofs.
            for h_text in (".002", ".001", ".0005"):
                h = ctx.mpf(h_text)
                actual = _actual_remainder(ctx, phase, h, kind, angular_functions, axial_functions, rates)
                actual_rows = [
                    ctx.diff(actual, z, order) for order in range(3)
                ]
                axial_truth = [
                    {"ordinary_rows": _interval_row(iv, row), "rate": iv.mpf(rate)}
                    for row, rate in zip(axial_rows, rates)
                ]
                angular_truth = [_interval_row(iv, row) for row in angular_rows]
                packet_truth = {
                    "leading": _interval_row(iv, packet["leading"]),
                    "next": _interval_row(iv, packet["next"]),
                }
                truth, certificate = enclose_weighted_remainder(
                    ctx,
                    iv,
                    coeffpackettruth=packet_truth,
                    axial_modestruth=axial_truth,
                    angular_modestruth=angular_truth,
                    W2positiveIV=iv.mpf(W2),
                    weight=iv.mpf(".1"),
                    logh=iv.log(iv.mpf(h)),
                    source_proved=True,
                    error_cap=".1",
                    eta=".01",
                )
                for order, value in enumerate(actual_rows):
                    lower, upper = _bounds(ctx, truth[order])
                    check(
                        f"{kind}_{phase_text}_h{h_text}_Z{order}_enclosed",
                        lower <= value <= upper,
                    )
                check(
                    f"{kind}_{phase_text}_h{h_text}_proofs",
                    certificate["passed"] is True
                    and certificate["h_selected"] is False
                    and certificate["hb_materialized"] is False
                    and len(certificate["width_log_proofs"]) == 3
                    and len(certificate["E3_log_proof_per_order"]) == 3,
                )
                if phase != 0:
                    leading_rows = [ctx.mpf(value) for value in packet["leading"]]
                    next_rows = [ctx.mpf(value) for value in packet["next"]]
                    e3_bounds = certificate["ordinary_error_bound_per_order_exact_mpf_tuples"]
                    for order, value in enumerate(actual_rows):
                        e3_upper = ctx.make_mpf(tuple(e3_bounds[order][1]))
                        center = leading_rows[order] + h * next_rows[order]
                        check(
                            f"{kind}_{phase_text}_h{h_text}_Z{order}_third_order_residual",
                            abs(value - center) <= h**2 * e3_upper + ctx.mpf("1e-45"),
                        )
                    # Independent finite-h evaluation has the predicted signed
                    # leading limit; the packet's R1 is retained in the center.
                    check(
                        f"{kind}_{phase_text}_h{h_text}_signed_leading",
                        abs(actual_rows[0] - leading_rows[0]) < ctx.mpf(".1"),
                    )

            if phase == 0:
                check(f"{kind}_phase_zero_W2_exact", W2 == 0)
                _, zero_certificate = enclose_weighted_remainder(
                    ctx,
                    iv,
                    coeffpackettruth={
                        "leading": _interval_row(iv, packet["leading"]),
                        "next": _interval_row(iv, packet["next"]),
                    },
                    axial_modestruth=[
                        {"ordinary_rows": _interval_row(iv, row), "rate": iv.mpf(rate)}
                        for row, rate in zip(axial_rows, rates)
                    ],
                    angular_modestruth=[_interval_row(iv, row) for row in angular_rows],
                    W2positiveIV=iv.mpf(0),
                    weight=iv.mpf(".1"),
                    logh=iv.log(iv.mpf(".001")),
                    source_proved=True,
                    error_cap=".1",
                    eta=".01",
                )
                check(
                    f"{kind}_phase_zero_errors_exact",
                    zero_certificate["phase_zero_exact"] is True
                    and all(proof["exact_zero"] for proof in zero_certificate["E3_log_proof_per_order"]),
                )

    for label, call in (
        (
            "missing_source_proof",
            lambda: signed_remainder_coefficient(
                axial_modes, angular_rows, ctx.mpf(".1"), ctx.mpf(".01"), source_proved=False
            ),
        ),
        (
            "missing_axial_rows",
            lambda: signed_remainder_coefficient(
                [{"rate": rates[0]}], angular_rows, ctx.mpf(".1"), ctx.mpf(".01"), source_proved=True
            ),
        ),
        (
            "bad_moment",
            lambda: signed_remainder_coefficient(
                axial_modes, angular_rows, ctx.mpf("-.1"), ctx.mpf(".01"), source_proved=True
            ),
        ),
        (
            "inconsistent_moment",
            lambda: signed_remainder_coefficients(
                axial_modes,
                angular_rows,
                ctx.mpf(".01"),
                ctx.mpf(".02"),
                ctx.mpf(".001"),
                True,
            ),
        ),
    ):
        try:
            call()
        except ValueError:
            check(label, True)
        else:
            check(label, False)

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "independent_oracle": "exact expm1 inner ell primitive, fresh outer Gauss-Legendre K quadrature and ctx.diff through ordinary Z2",
        "kappa_fixtures": ["1", "1-t"],
        "source_proved": True,
        "actual_field_zeroed": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
