"""Independent checks for the reserved m15 response algebra.

The fixture below deliberately bypasses ``OriginalReservedM15Response.__init__``.
It uses a finite synthetic ``mu``, the actual smooth compact bump, and an
angular matrix reconstructed by independent high precision quadrature.  No
accepted production archive or ancestor constructor is loaded.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "root_st073"))

from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from openai_ns_reconstruction.original_core_numeric import _bounds
from openai_ns_reconstruction.original_reserved_m15_response import (
    OriginalReservedM15Response,
    inverse3,
)


def smooth_bump(c, t):
    """The defining C-infinity bump, evaluated as a scalar on ``c``."""
    t = c.mpf(t)
    if abs(t) >= 1:
        return c.mpf(0)
    return c.exp(-1 / (1 - t * t))


def integrate(c, fn, left, right):
    if right <= left:
        return c.mpf(0)
    midpoint = (left + right) / 2
    return c.quad(fn, [left, midpoint, right])


def determinant3(matrix):
    a = matrix
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inverse3_scalar(matrix, determinant):
    """Independent scalar implementation, with a different loop layout."""
    out = [[None] * 3 for _ in range(3)]
    for row in range(3):
        for col in range(3):
            minor = [
                [matrix[r][c] for c in range(3) if c != row]
                for r in range(3)
                if r != col
            ]
            out[row][col] = (
                (-1) ** (row + col)
                * (minor[0][0] * minor[1][1] - minor[0][1] * minor[1][0])
                / determinant
            )
    return out


def close(c, left, right, tolerance):
    return abs(left - right) <= tolerance * max(c.mpf(1), abs(left), abs(right))


def contains(c, box, value, tolerance):
    lo, hi = _bounds(c, box)
    return lo - tolerance <= value <= hi + tolerance


def make_fixture():
    """Create a synthetic owner and independent scalar data."""
    c = mp.mp.clone()
    c.dps = 75
    iv = MPIntervalContext()
    iv.dps = 100

    # These are finite and well separated so the partial-support case has one
    # genuinely active bump and two disjoint supports.
    mu = c.mpf(".08")
    ell = c.mpf(".025")
    centers = [c.mpf(".18"), c.mpf(".36"), c.mpf(".54")]

    normal = integrate(
        c,
        lambda t: smooth_bump(c, t) * c.exp(c.mpf(".5") * ell * t),
        c.mpf(-1),
        c.mpf(1),
    )

    # M_j(k) = exp(k c_j) J(k + 1/2)/J(1/2).  The third row is precisely
    # the physical partial-row exponent k=-mu used by partial_row().
    row_exponents = [c.mpf(1), -(1 + mu), -mu]
    matrix_scalar = []
    for exponent in row_exponents:
        row = []
        for center in centers:
            integral = integrate(
                c,
                lambda t, exponent=exponent: smooth_bump(c, t)
                * c.exp((c.mpf(".5") + exponent) * ell * t),
                c.mpf(-1),
                c.mpf(1),
            )
            row.append(c.exp(exponent * center) * integral / normal)
        matrix_scalar.append(row)
    determinant = determinant3(matrix_scalar)
    if determinant <= 0:
        raise AssertionError("synthetic angular determinant must be positive")

    owner = object.__new__(OriginalReservedM15Response)
    owner.ctx = c
    owner.iv = iv
    owner.J = IntervalTaylor
    owner.beta_jets = beta_jets
    owner.cells = 128
    owner._pressure_integrals = {}
    owner.mu = iv.mpf(mu)
    owner.ell = iv.mpf(ell)
    # The production method intentionally addresses beta1, beta2, beta3 as
    # centers[1], centers[2], centers[3]; retain that one-based shape here.
    owner.centers = [iv.mpf(0)] + [iv.mpf(v) for v in centers]
    owner.normal = iv.mpf(normal)
    owner.matrix = [[iv.mpf(v) for v in row] for row in matrix_scalar]
    owner.inverse = inverse3(owner.matrix, iv.mpf(determinant))
    return owner, c, iv, centers, normal, matrix_scalar, determinant


def direct_prefix_and_suffix(c, mu, ell, centers, normal, x):
    """Independent scalar quadrature of the two physical support pieces."""
    s = c.log(x)
    rate = c.mpf(".5") - mu
    forward, suffix = [], []
    for center in centers:
        r = (s - center) / ell
        factor = c.exp(-mu * center) / normal

        def integrand(t):
            return smooth_bump(c, t) * c.exp(rate * ell * t)

        if r <= -1:
            f = c.mpf(0)
            b = factor * integrate(c, integrand, c.mpf(-1), c.mpf(1))
        elif r >= 1:
            f = factor * integrate(c, integrand, c.mpf(-1), c.mpf(1))
            b = c.mpf(0)
        else:
            f = factor * integrate(c, integrand, c.mpf(-1), r)
            b = factor * integrate(c, integrand, r, c.mpf(1))
        forward.append(f)
        suffix.append(b)
    return forward, suffix


def matrix_vector_row(row, matrix):
    return [sum(row[j] * matrix[j][k] for j in range(3)) for k in range(3)]


def direct_physical_response(c, mu, ell, centers, normal, inverse, x, z):
    """Direct -R theta0 Utheta1 d(log R), including C0 and ordinary C1."""
    # This is the original rational A_rc dependence used by the response
    # module, with nonconstant signed theta inputs so product-rule mistakes are
    # observable.
    base = c.exp(-c.mpf("1.7") - c.mpf("2.5") * mu)
    A = base / (1 + z * z)
    Az = -2 * z * base / (1 + z * z) ** 2
    theta = [
        c.mpf(".63") + c.mpf(".11") * z - c.mpf(".04") * z * z,
        -c.mpf(".28") + c.mpf(".17") * z + c.mpf(".06") * z * z,
        c.mpf(".41") - c.mpf(".09") * z + c.mpf(".025") * z * z * z,
    ]
    theta_z = [
        c.mpf(".11") - c.mpf(".08") * z,
        c.mpf(".17") + c.mpf(".12") * z,
        -c.mpf(".09") + c.mpf(".075") * z * z,
    ]
    beta = [sum(inverse[j][k] * theta[k] for k in range(3)) for j in range(3)]
    beta_z = [sum(inverse[j][k] * theta_z[k] for k in range(3)) for j in range(3)]
    s_end = c.log(x)
    alpha = c.mpf(".5") + mu
    rate = c.mpf(".5") - mu
    c0 = c.mpf(0)
    c1 = c.mpf(0)
    for j, center in enumerate(centers):
        upper = min(s_end, center + ell)
        lower = max(c.mpf(0), center - ell)

        def shape(s):
            return (
                c.exp(-center / 2)
                / (ell * normal)
                * smooth_bump(c, (s - center) / ell)
            )

        if upper <= lower:
            integral = c.mpf(0)
        else:
            integral = integrate(
                c,
                lambda s: c.exp(rate * s) * shape(s),
                lower,
                upper,
            )
        # theta0(R,Z) = A(Z) exp(-(1/2+mu) log(R/Rc)); R/Rc = exp(s),
        # and d(log R)=ds.  The extra R therefore produces exp(rate*s).
        c0 -= A * beta[j] * integral
        c1 -= (Az * beta[j] + A * beta_z[j]) * integral
    return c0, c1, (A, Az, theta, theta_z, beta, beta_z)


def run():
    started = time.perf_counter()
    owner, c, iv, centers, normal, matrix_scalar, determinant = make_fixture()
    tolerance = c.mpf("2e-48")

    # Independent inverse comparison catches row/column transpose and cofactor
    # indexing errors before any pressure calculation is attempted.
    scalar_inverse = inverse3_scalar(matrix_scalar, determinant)
    inverse_checks = 0
    max_inverse_error = c.mpf(0)
    for i in range(3):
        for j in range(3):
            lo, hi = _bounds(c, owner.inverse[i][j])
            expected = scalar_inverse[i][j]
            if not (lo - tolerance <= expected <= hi + tolerance):
                raise AssertionError(("inverse3", i, j, lo, hi, expected))
            max_inverse_error = max(max_inverse_error, abs(expected - (lo + hi) / 2))
            inverse_checks += 1

    # Independent M M^-1 and M^-1 M identity checks use the scalar quadrature
    # matrix and independent scalar inverse, rather than the interval matrix.
    identity_error = c.mpf(0)
    for left, right in ((matrix_scalar, scalar_inverse), (scalar_inverse, matrix_scalar)):
        for i in range(3):
            for j in range(3):
                value = sum(left[i][k] * right[k][j] for k in range(3))
                identity_error = max(identity_error, abs(value - (1 if i == j else 0)))
                if identity_error > c.mpf("1e-65"):
                    raise AssertionError(("inverse identity", i, j, value))

    # The full source row is independently reconstructed by quadrature, not
    # inferred from a midpoint or from the production matrix.
    full_forward, full_suffix = direct_prefix_and_suffix(
        c, c.mpf(".08"), c.mpf(".025"), centers, normal, c.mpf(2)
    )
    for j in range(3):
        if not close(c, full_forward[j], matrix_scalar[2][j], c.mpf("1e-60")):
            raise AssertionError(("full source row", j, full_forward[j], matrix_scalar[2][j]))
        if full_suffix[j] != 0:
            raise AssertionError(("full suffix", j, full_suffix[j]))

    cases = {"empty": "1", "partial": c.nstr(c.exp(centers[1]), 75), "full": "2"}
    case_reports = {}
    direct_response_by_case = {}
    for label, x_text in cases.items():
        x = c.mpf(x_text)
        response, diagnostic = owner.partial_row(x_text)
        forward, suffix = direct_prefix_and_suffix(
            c, c.mpf(".08"), c.mpf(".025"), centers, normal, x
        )
        expected = matrix_vector_row(forward, scalar_inverse)
        expected_suffix = matrix_vector_row(suffix, scalar_inverse)
        for k in range(3):
            if not contains(c, response[k], expected[k], tolerance):
                raise AssertionError(("partial row", label, k, response[k], expected[k]))

        if label == "empty":
            if not diagnostic["exact_empty"] or diagnostic["exact_full"]:
                raise AssertionError(("empty flags", diagnostic))
            if any(response[k]._mpi_ != iv.mpf(0)._mpi_ for k in range(3)):
                raise AssertionError(("empty response", response))
            if any(value != 0 for value in forward):
                raise AssertionError(("empty forward quadrature", forward))
            for j in range(3):
                if not close(c, suffix[j], full_forward[j], c.mpf("1e-60")):
                    raise AssertionError(("empty suffix/full source", j, suffix[j], full_forward[j]))
        elif label == "full":
            if diagnostic["exact_empty"] or not diagnostic["exact_full"]:
                raise AssertionError(("full flags", diagnostic))
            target = [c.mpf(0), c.mpf(0), c.mpf(1)]
            for k in range(3):
                if not contains(c, response[k], target[k], tolerance):
                    raise AssertionError(("full e2", k, response[k]))
            if any(abs(expected_suffix[k]) > c.mpf("1e-60") for k in range(3)):
                raise AssertionError(("full direct suffix", expected_suffix))
        else:
            if diagnostic["exact_empty"] or diagnostic["exact_full"]:
                raise AssertionError(("partial flags", diagnostic))
            if not any(abs(value) > c.mpf("1e-8") for value in expected):
                raise AssertionError(("partial response unexpectedly zero", expected))
            if max(abs(expected[k] - (1 if k == 2 else 0)) for k in range(3)) < c.mpf("1e-8"):
                raise AssertionError(("partial response falsely full", expected))
        direct_response_by_case[label] = expected
        case_reports[label] = dict(
            x=x_text,
            exact_empty=diagnostic["exact_empty"],
            exact_full=diagnostic["exact_full"],
            response=[
                [c.nstr(v, 16) for v in _bounds(c, row)] for row in response
            ],
            direct_response=[c.nstr(v, 18) for v in expected],
            direct_suffix_response=[c.nstr(v, 18) for v in expected_suffix],
        )

    # C0/C1 physical response at a nontrivial partial support.  The direct
    # integral below is a separate quadrature in physical log-radius s; it
    # includes -R*theta0*Utheta1*d(log R), with R=exp(s), and uses A(Z) and
    # A_Z(Z) from the original rational A_rc form.
    x_text = cases["partial"]
    z = c.mpf(".37")
    response, _ = owner.partial_row(x_text)
    c0_direct, c1_direct, source = direct_physical_response(
        c, c.mpf(".08"), c.mpf(".025"), centers, normal, scalar_inverse, c.mpf(x_text), z
    )
    A, Az, theta, theta_z, beta, beta_z = source
    input_c0 = [A * value for value in theta]
    input_c1 = [Az * theta[k] + A * theta_z[k] for k in range(3)]
    module_c0 = -sum(response[k] * iv.mpf(input_c0[k]) for k in range(3))
    module_c1 = -sum(response[k] * iv.mpf(input_c1[k]) for k in range(3))
    if not contains(c, module_c0, c0_direct, tolerance):
        raise AssertionError(("physical C0", module_c0, c0_direct))
    if not contains(c, module_c1, c1_direct, tolerance):
        raise AssertionError(("physical C1", module_c1, c1_direct))

    # Explicit negative controls: a transpose, a sign flip, omitted R, or an
    # omitted ordinary product-rule term must disagree by a visible amount.
    forward, _ = direct_prefix_and_suffix(
        c, c.mpf(".08"), c.mpf(".025"), centers, normal, c.mpf(x_text)
    )
    correct_response = matrix_vector_row(forward, scalar_inverse)
    transpose_response = [sum(forward[j] * scalar_inverse[k][j] for j in range(3)) for k in range(3)]
    transpose_gap = max(abs(correct_response[k] - transpose_response[k]) for k in range(3))
    if transpose_gap <= c.mpf("1e-8"):
        raise AssertionError(("transpose control too small", transpose_gap))

    wrong_sign = -c0_direct
    sign_gap = abs(c0_direct - wrong_sign)
    if sign_gap <= c.mpf("1e-8"):
        raise AssertionError(("sign control too small", sign_gap))

    # Omitting the extra R changes exp((1/2-mu)s) to exp(-(1/2+mu)s).
    alpha = c.mpf(".5") + c.mpf(".08")
    s_end = c.log(c.mpf(x_text))
    no_R = c.mpf(0)
    for j, center in enumerate(centers):
        lower = max(c.mpf(0), center - c.mpf(".025"))
        upper = min(s_end, center + c.mpf(".025"))
        if upper > lower:
            no_R += integrate(
                c,
                lambda s, center=center: c.exp(-alpha * s)
                * c.exp(-center / 2)
                / (c.mpf(".025") * normal)
                * smooth_bump(c, (s - center) / c.mpf(".025")),
                lower,
                upper,
            ) * (A * beta[j])
    no_R = -no_R
    missing_R_gap = abs(c0_direct - no_R)
    if missing_R_gap <= c.mpf("1e-8"):
        raise AssertionError(("missing-R control too small", missing_R_gap))

    wrong_c1_without_Az = -sum(correct_response[k] * A * theta_z[k] for k in range(3))
    wrong_c1_without_theta_z = -sum(correct_response[k] * Az * theta[k] for k in range(3))
    product_rule_gap = max(
        abs(c1_direct - wrong_c1_without_Az),
        abs(c1_direct - wrong_c1_without_theta_z),
    )
    if abs(wrong_c1_without_Az - c1_direct) <= c.mpf("1e-8") or abs(
        wrong_c1_without_theta_z - c1_direct
    ) <= c.mpf("1e-8"):
        raise AssertionError(("ordinary product-rule control too small", wrong_c1_without_Az, wrong_c1_without_theta_z, c1_direct))

    return dict(
        passed=True,
        execution_seconds=round(time.perf_counter() - started, 3),
        fixture=dict(
            mu="0.08",
            ell="0.025",
            centers=[c.nstr(v, 18) for v in centers],
            normal=c.nstr(normal, 18),
            angular_row_exponents=["1", "-1-mu", "-mu"],
            quadrature="mpmath tanh-sinh over each smooth-bump support",
            source_midpoints_used=False,
            accepted_ancestor_constructors_loaded=False,
            synthetic_owner="object.__new__(OriginalReservedM15Response)",
        ),
        checks=dict(
            inverse3_entries=inverse_checks,
            inverse3_max_interval_midpoint_error=c.nstr(max_inverse_error, 8),
            inverse_identity_max_error=c.nstr(identity_error, 8),
            independently_rebuilt_full_third_row=True,
            exact_empty_support_e0=True,
            exact_full_support_e2=True,
            nontrivial_partial_support=True,
            direct_physical_C0_contains=True,
            direct_physical_C1_contains=True,
            transpose_gap=c.nstr(transpose_gap, 12),
            sign_gap=c.nstr(sign_gap, 12),
            missing_extra_R_gap=c.nstr(missing_R_gap, 12),
            ordinary_product_rule_gap=c.nstr(product_rule_gap, 12),
            original_theta_inputs=True,
            rational_A_rc_and_A_Z=True,
        ),
        cases=case_reports,
        physical=dict(
            direct_C0=c.nstr(c0_direct, 18),
            direct_C1=c.nstr(c1_direct, 18),
            A=c.nstr(A, 18),
            A_Z=c.nstr(Az, 18),
            theta=[c.nstr(v, 18) for v in theta],
            theta_Z=[c.nstr(v, 18) for v in theta_z],
            no_extra_R_C0=c.nstr(no_R, 18),
        ),
        limits=dict(
            synthetic_math_check_only=True,
            production_archives_not_loaded=True,
            full_original_m15_closure=False,
            corrected_NS_validation=False,
        ),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2), flush=True)
