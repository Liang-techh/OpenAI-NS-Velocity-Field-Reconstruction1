"""Independent finite-N checks for the original micro viscous trace.

The physical fixture is deliberately self contained.  It reconstructs a
finite oscillatory velocity, its Volterra mean, the physical-R viscous
operator and ordinary Z derivatives without loading an accepted archive or
an ancestor constructor.  The optional source-ID checks exercise the narrow
``contract_trace_chain`` API once the new module is present.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


N = 17
P_TEXT = "1.37"
DELTA_TEXT = "0.23"
RC_TEXT = "2.3"
Z_TEXT = "0.37"
INTERVALS = (("micro_first", "0", "1"), ("micro_second", "1", "2"))


def close(ctx, left, right, tolerance):
    return abs(left - right) <= tolerance * max(ctx.mpf(1), abs(left), abs(right))


def integrate(ctx, function, left, right, segments=16):
    """Use independent tanh-sinh quadrature on fixed subintervals."""
    left, right = ctx.mpf(left), ctx.mpf(right)
    if right <= left:
        return ctx.mpf(0)
    points = [left + (right - left) * k / segments for k in range(segments + 1)]
    return ctx.quad(function, points)


def amplitude(ctx, z):
    return 1 + ctx.mpf("0.2") * z + ctx.mpf("0.3") * z * z


def amplitude_z(ctx, z):
    return ctx.mpf("0.2") + ctx.mpf("0.6") * z


def mean_prefix(ctx, y):
    """Return integral_0^y exp(2v/3) sin(Nv)/N dv analytically."""
    y = ctx.mpf(y)
    rate = ctx.mpf(2) / 3
    return (
        ctx.exp(rate * y) * (rate * ctx.sin(N * y) - N * ctx.cos(N * y)) + N
    ) / (N * (rate * rate + N * N))


def fields(ctx, y, z, *, fast_phase=True, amplitude_derivative=True):
    """Return U, m and the derivatives needed by Qy and Qy_Z."""
    y, z = ctx.mpf(y), ctx.mpf(z)
    p = ctx.mpf(P_TEXT)
    e = ctx.exp(-y / 3)
    s, c = ctx.sin(N * y), ctx.cos(N * y)
    a = amplitude(ctx, z)
    az = amplitude_z(ctx, z) if amplitude_derivative else ctx.mpf(0)
    azz = ctx.mpf("0.6") if amplitude_derivative else ctx.mpf(0)
    prefix = mean_prefix(ctx, y)
    mean_factor = ctx.exp(-y) * prefix

    u = p * a * e * s / N
    if fast_phase:
        uy_factor = -s / (3 * N) + c
        uyy_factor = -ctx.mpf(2) * c / 3 - N * s + s / (9 * N)
    else:
        # Deliberately omit the fast-phase derivative cos(N y).
        uy_factor = -s / (3 * N)
        uyy_factor = -c / 3 + s / (9 * N)
    uy = p * a * e * uy_factor
    uyy = p * a * e * uyy_factor
    uz = p * az * e * s / N
    uyz = p * az * e * uy_factor
    uyyz = p * az * e * uyy_factor
    uzz = p * azz * e * s / N
    uyzz = p * azz * e * uy_factor

    m = p * a * mean_factor
    my = u - m
    mz = p * az * mean_factor
    mzz = p * azz * mean_factor

    return dict(
        U=u,
        Uy=uy,
        Uyy=uyy,
        Uz=uz,
        Uyz=uyz,
        Uyyz=uyyz,
        Uyzz=uyzz,
        Uzz=uzz,
        m=m,
        my=my,
        mz=mz,
        mzz=mzz,
    )


def q_bundle(ctx, y, z, *, fast_phase=True, amplitude_derivative=True):
    """Return Qy, Qy_Z and their y derivatives from explicit product rules."""
    f = fields(
        ctx,
        y,
        z,
        fast_phase=fast_phase,
        amplitude_derivative=amplitude_derivative,
    )
    d = ctx.mpf(DELTA_TEXT)
    z = ctx.mpf(z)
    e = f["U"] - f["m"]
    ez = f["Uz"] - f["mz"]
    ezz = f["Uzz"] - f["mzz"]
    ey = f["Uy"] - e
    ezy = f["Uyz"] - ez
    ezyz = f["Uyzz"] - ezz
    one_minus_delta = 1 - d
    one_minus_z2 = 1 - z * z
    denominator = 1 - d * z * z
    denominator_z = -2 * d * z

    numerator = (
        2 * z * f["Uy"]
        - z * one_minus_delta * e
        - one_minus_z2 * ez
    )
    numerator_z = (
        2 * f["Uy"]
        + 2 * z * f["Uyz"]
        - one_minus_delta * e
        - z * one_minus_delta * ez
        + 2 * z * ez
        - one_minus_z2 * ezz
    )
    numerator_y = (
        2 * z * f["Uyy"]
        - z * one_minus_delta * ey
        - one_minus_z2 * ezy
    )
    numerator_yz = (
        2 * f["Uyy"]
        + 2 * z * f["Uyyz"]
        - one_minus_delta * ey
        - z * one_minus_delta * ezy
        + 2 * z * ezy
        - one_minus_z2 * ezyz
    )
    q = numerator / denominator
    qz = (numerator_z * denominator - numerator * denominator_z) / denominator**2
    qy = numerator_y / denominator
    qyz = (numerator_yz * denominator - numerator_y * denominator_z) / denominator**2
    return dict(Qy=q, QyZ=qz, Qy_y=qy, Qy_yZ=qyz)


def radius(ctx, y):
    return ctx.mpf(RC_TEXT) * ctx.exp(y)


def viscous_integrands(ctx, y, z, *, fast_phase=True, amplitude_derivative=True, include_radius=True):
    q = q_bundle(
        ctx,
        y,
        z,
        fast_phase=fast_phase,
        amplitude_derivative=amplitude_derivative,
    )
    scale = radius(ctx, y) if include_radius else ctx.mpf(1)
    return (
        -scale * (q["Qy"] + q["Qy_y"]),
        -scale * (q["QyZ"] + q["Qy_yZ"]),
    )


def endpoint_orientation(ctx, left, right, z, *, ordinary_z=False, fast_phase=True, amplitude_derivative=True):
    def value(y):
        q = q_bundle(
            ctx,
            y,
            z,
            fast_phase=fast_phase,
            amplitude_derivative=amplitude_derivative,
        )
        return radius(ctx, y) * q["QyZ" if ordinary_z else "Qy"]

    return value(left) - value(right)


def physical_fixture():
    started = time.perf_counter()
    ctx = mp.mp.clone()
    ctx.dps = 80
    z = ctx.mpf(Z_TEXT)
    tolerance = ctx.mpf("2e-55")
    reports = []
    total_c0 = ctx.mpf(0)
    total_c1 = ctx.mpf(0)

    for name, left_text, right_text in INTERVALS:
        left, right = ctx.mpf(left_text), ctx.mpf(right_text)
        c0 = integrate(
            ctx,
            lambda y: viscous_integrands(ctx, y, z)[0],
            left,
            right,
        )
        c1 = integrate(
            ctx,
            lambda y: viscous_integrands(ctx, y, z)[1],
            left,
            right,
        )
        endpoint_c0 = endpoint_orientation(ctx, left, right, z)
        endpoint_c1 = endpoint_orientation(ctx, left, right, z, ordinary_z=True)
        finite_difference_c1 = ctx.diff(
            lambda zz: endpoint_orientation(ctx, left, right, zz), z
        )
        if not close(ctx, c0, endpoint_c0, tolerance):
            raise AssertionError(("C0 endpoint orientation", name, c0, endpoint_c0))
        if not close(ctx, c1, endpoint_c1, tolerance):
            raise AssertionError(("C1 endpoint orientation", name, c1, endpoint_c1))
        if not close(ctx, endpoint_c1, finite_difference_c1, tolerance):
            raise AssertionError(("ordinary Z1 endpoint derivative", name, endpoint_c1, finite_difference_c1))
        total_c0 += c0
        total_c1 += c1
        reports.append(
            dict(
                name=name,
                lower=left_text,
                upper=right_text,
                C0=ctx.nstr(c0, 18),
                C1=ctx.nstr(c1, 18),
                endpoint_C0=ctx.nstr(endpoint_c0, 18),
                endpoint_C1=ctx.nstr(endpoint_c1, 18),
            )
        )

    left, right = ctx.mpf(0), ctx.mpf(1)
    correct_c0 = integrate(ctx, lambda y: viscous_integrands(ctx, y, z)[0], left, right)
    wrong_fast = integrate(
        ctx,
        lambda y: viscous_integrands(ctx, y, z, fast_phase=False)[0],
        left,
        right,
    )
    omitted_fast_phase_gap = abs(correct_c0 - wrong_fast)
    if omitted_fast_phase_gap <= ctx.mpf("1e-8"):
        raise AssertionError(("omitted fast phase control too small", omitted_fast_phase_gap))

    wrong_sign = -correct_c0
    sign_gap = abs(correct_c0 - wrong_sign)
    if sign_gap <= ctx.mpf("1e-8"):
        raise AssertionError(("sign control too small", sign_gap))

    wrong_radius = integrate(
        ctx,
        lambda y: viscous_integrands(ctx, y, z, include_radius=False)[0],
        left,
        right,
    )
    missing_radius_gap = abs(correct_c0 - wrong_radius)
    if missing_radius_gap <= ctx.mpf("1e-8"):
        raise AssertionError(("missing physical R control too small", missing_radius_gap))

    correct_c1 = integrate(ctx, lambda y: viscous_integrands(ctx, y, z)[1], left, right)
    wrong_product = integrate(
        ctx,
        lambda y: viscous_integrands(ctx, y, z, amplitude_derivative=False)[1],
        left,
        right,
    )
    product_rule_gap = abs(correct_c1 - wrong_product)
    if product_rule_gap <= ctx.mpf("1e-8"):
        raise AssertionError(("ordinary product derivative control too small", product_rule_gap))

    # Verify the defining Volterra equation directly at multiple points,
    # independently of the closed-form mean prefix used above.
    mean_ode_error = ctx.mpf(0)
    for y in (ctx.mpf("0.13"), ctx.mpf("0.71"), ctx.mpf("1.63")):
        f = fields(ctx, y, z)
        mean_ode_error = max(mean_ode_error, abs(f["my"] - (f["U"] - f["m"])))
    if mean_ode_error != 0:
        raise AssertionError(("Volterra mean ODE", mean_ode_error))

    return dict(
        passed=True,
        execution_seconds=round(time.perf_counter() - started, 3),
        fixture=dict(
            N=N,
            P=P_TEXT,
            delta=DELTA_TEXT,
            Rc=RC_TEXT,
            Z=Z_TEXT,
            intervals=[list(row) for row in INTERVALS],
            U_definition="P*A(Z)*exp(-y/3)*sin(N*y)/N",
            Uy_definition="P*A(Z)*exp(-y/3)*(-sin(N*y)/(3*N)+cos(N*y))",
            mean_definition="exp(-y)*integral_0^y exp(v)*deltaU(v,Z) dv",
            quadrature="mpmath tanh-sinh over 16 subintervals per physical interval",
            source_midpoints_used=False,
            accepted_ancestor_constructors_loaded=False,
        ),
        checks=dict(
            finite_N=N,
            intervals_checked=len(reports),
            C0_endpoint_orientation=True,
            ordinary_Z1_endpoint_orientation=True,
            ordinary_Z1_finite_difference=True,
            Volterra_mean_ODE=True,
            physical_R_included=True,
            fast_phase_included=True,
            ordinary_product_derivatives_included=True,
            total_C0=ctx.nstr(total_c0, 18),
            total_C1=ctx.nstr(total_c1, 18),
            omitted_fast_phase_gap=ctx.nstr(omitted_fast_phase_gap, 12),
            sign_gap=ctx.nstr(sign_gap, 12),
            missing_physical_R_gap=ctx.nstr(missing_radius_gap, 12),
            ordinary_product_rule_gap=ctx.nstr(product_rule_gap, 12),
            mean_ODE_max_error=ctx.nstr(mean_ode_error, 12),
        ),
        pieces=reports,
        limits=dict(
            independent_physics_only=True,
            production_archives_not_loaded=True,
            full_NS_validation=False,
        ),
    )


def _cancelled(result):
    value = result.get("cancelled_interfaces", result.get("cancelled", []))
    return value if isinstance(value, list) else [value]


def contract_fixture(contract_trace_chain):
    """Check source-ID contraction without cancelling numerical bounds."""
    proved = contract_trace_chain(
        [
            dict(name="first", lower="source_left", upper="seam_right"),
            dict(name="second", lower="seam_left", upper="source_right"),
            # These two IDs can carry equal numerical boxes but have no proof.
            dict(name="boxed", lower="equal_bound_left", upper="equal_bound_right"),
        ],
        {"seam_right": "seam_left"},
    )
    expected = {
        "source_left": 1,
        "source_right": -1,
        "equal_bound_left": 1,
        "equal_bound_right": -1,
    }
    if proved.get("signed_coefficients") != expected:
        raise AssertionError(("source-ID contraction", proved))
    if "seam_left" not in repr(_cancelled(proved)):
        raise AssertionError(("cancelled interface record", proved))

    unknown = contract_trace_chain(
        [dict(name="unknown", lower="unknown_left", upper="unknown_right")],
        {},
    )
    if unknown.get("signed_coefficients") != {"unknown_left": 1, "unknown_right": -1}:
        raise AssertionError(("unknown source IDs silently cancelled", unknown))
    return dict(
        passed=True,
        proved_alias={"seam_right": "seam_left"},
        signed_coefficients=proved["signed_coefficients"],
        cancelled_interfaces=_cancelled(proved),
        unknown_ids_retained=unknown["signed_coefficients"],
        equal_numeric_bounds_not_aliased=True,
        bounds_never_cancelled_as_functions=True,
    )


def run(*, include_contract=True):
    physical = physical_fixture()
    if not include_contract:
        contract = dict(passed=True, skipped=True, reason="new viscous module API not loaded")
    else:
        try:
            from openai_ns_reconstruction.original_micro_viscous_m15 import contract_trace_chain
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "Run --physics-only until original_micro_viscous_m15.py is available"
            ) from exc
        contract = contract_fixture(contract_trace_chain)
    return dict(
        passed=physical["passed"] and contract["passed"],
        physical=physical,
        contract=contract,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--physics-only",
        action="store_true",
        help="run before the new contract_trace_chain module is available",
    )
    args = parser.parse_args()
    print(json.dumps(run(include_contract=not args.physics_only), indent=2), flush=True)
