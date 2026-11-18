"""Fresh checks for the finite-mode Part II auxiliary operators."""
from __future__ import annotations

import cmath
import math

from scipy.special import roots_legendre

from openai_ns_reconstruction.lei_ren_auxiliary_torus import (
    inverse_directional_derivative,
    pullback,
    transported_radial_primitive,
)


def value(rows, y):
    return sum(c * cmath.exp(2j * math.pi * (n*y[0] + m*y[1])) for (n, m), c in rows.items())


def run():
    # Paper (2.6)-(2.9). These choices exercise adjacent reference bands
    # at an admissible delta, independently of the C1 and ST073 fields.
    delta = 1/256
    lr, lt = 4 - math.sqrt(2), 4 + math.sqrt(2)
    vr, vt = (1, 1 - math.sqrt(2)), (math.sqrt(2) - 1, 1)
    exponent = (2 + delta)*math.log(lr)/math.log(lt) - 1e-5*delta

    def band_level(k):
        return math.floor(((2 + delta)*k*math.log(2) - 2*math.log(k))/math.log(lt))

    bands = (64, 65)
    lower, upper = map(band_level, bands)
    assert upper == lower + 1
    scale = math.ldexp(1.0, -bands[0])
    radius = 1.43*scale
    y = (0.32, 0.21)
    ey = ((3*y[0]+y[1]) % 1, (y[0]+5*y[1]) % 1)
    rows = {(0, 0): 3, (2, -1): 0.5+0.25j, (-2, 1): 0.5-0.25j, (1, 0): -0.1j, (-1, 0): 0.1j}
    common = pullback(rows, level=upper, common_level=lower)
    assert set(common) == {(0, 0), (5, -3), (-5, 3), (3, 1), (-3, -1)}
    assert common[0, 0] == rows[0, 0]
    phase_error = abs(value(common, y) - value(rows, ey))
    assert phase_error < 3e-14
    oscillatory = {mode: -c for mode, c in rows.items() if mode != (0, 0)}
    common_osc = pullback(oscillatory, level=upper, common_level=lower)

    inverse_errors = {}
    for direction, lam, v in (("time", lt, vt), ("radial", lr, vr)):
        local = inverse_directional_derivative(oscillatory, level=upper, direction=direction)
        other = inverse_directional_derivative(common_osc, level=lower, direction=direction)
        pulled = pullback(local, level=upper, common_level=lower)
        relative = max(abs(other[mode] / c - 1) for mode, c in pulled.items())
        reconstructed = {mode: c*(2j*math.pi*lam**lower*(mode[0]*v[0]+mode[1]*v[1])) for mode, c in other.items()}
        residual = max(abs(reconstructed[mode] - common_osc[mode]) for mode in reconstructed)
        assert relative < 4e-14 and residual < 4e-14
        inverse_errors[direction] = {"cross_band_relative": relative, "source_recovery_abs": residual}

    try:
        inverse_directional_derivative(rows, level=upper)
    except ValueError as error:
        assert "zero mean" in str(error)
    else:
        raise AssertionError("the directional inverse accepted a mean")

    # This nonconstant Pell mode would otherwise lose the small divisor.
    a, b = 1, 0
    for _ in range(24):
        a, b = 3*a+4*b, 2*a+3*b
    assert a*a-2*b*b == 1
    # Evaluate the nearly cancelling time divisor of (-b,a-b).
    half = inverse_directional_derivative({(-b, a-b): 1}, level=0)[-b, a-b]
    reference = -1j * (a + math.sqrt(2)*b) / (2*math.pi)
    pell_error = abs(half / reference - 1)
    assert pell_error < 4e-15

    def cutoff(r):
        u = (r/scale - 1.15)/0.7
        if u <= 0:
            return 0, 0
        if u >= 1:
            return 1, 0
        # The fixture works inside this C² step. Production admits a smooth cutoff.
        return u**3*(10 + u*(-15+6*u)), 30*u**2*(1-u)**2/(0.7*scale)

    def profile(s):
        x = s/scale
        return (x-1)*(2-x)

    forcing = {mode: -c for mode, c in inverse_directional_derivative(oscillatory, level=upper).items()}
    # Keep amplitude O(1) for this binary64 quadrature fixture.
    forcing = {mode: c * lt**upper for mode, c in forcing.items()}
    common_forcing = pullback(forcing, level=upper, common_level=lower)

    def evaluate(r, a, i, modes):
        return transported_radial_primitive(
            lambda s, q: profile(s)*value(modes, q), radius=r, auxiliary=a, annulus=(scale, 2*scale),
            cutoff=cutoff, level=i, radial_exponent=exponent, atol=1e-12*scale,
        )

    first = evaluate(radius, ey, upper, forcing)
    second = evaluate(radius, y, lower, common_forcing)
    cross_band = abs(first.value - second.value)/scale
    cross_band_dr = abs(first.physical_radial_derivative - second.physical_radial_derivative)

    def physical_y(s):
        lift = lr**lower*(s**exponent-radius**exponent)
        return y[0] + lift*vr[0], y[1] + lift*vr[1]

    # Independent Gauss–Legendre reference along the realized radius path.
    nodes, weights = roots_legendre(48)

    def reference_integral(a, b):
        return (b-a)/2 * sum(
            w*profile(s)*value(common_forcing, physical_y(s))
            for x, w in zip(nodes, weights)
            for s in (a + (b-a)*(x+1)/2,)
        )

    w, dw = cutoff(radius)
    total = reference_integral(scale, 2*scale)
    reference_value = reference_integral(scale, radius) - w*total
    ref_error = abs(second.value-reference_value)/scale
    derivative_error = abs(second.physical_radial_derivative - (profile(radius)*value(common_forcing, y)-dw*total))
    step = 1e-5*scale
    p = evaluate(radius+step, physical_y(radius+step), lower, common_forcing).value
    m = evaluate(radius-step, physical_y(radius-step), lower, common_forcing).value
    fd_error = abs((p-m)/(2*step)-second.physical_radial_derivative)
    assert cross_band < 1e-13 and cross_band_dr < 1e-13
    assert ref_error < 1e-13 and derivative_error < 1e-13 and fd_error < 1e-9
    for r in (0, scale, 2*scale, 3*scale):
        exterior = evaluate(r, y, lower, common_forcing)
        assert exterior.value == exterior.physical_radial_derivative == 0

    return {
        "bands": bands, "covering_levels": [lower, upper], "delta": delta,
        "pullback_phase_abs_error": phase_error, "inverse": inverse_errors,
        "pell_relative_error": pell_error, "rejected_nonzero_mean": True,
        "scaled_primitive_cross_band_abs_error": cross_band,
        "radial_derivative_cross_band_abs_error": cross_band_dr,
        "scaled_independent_radial_quadrature_abs_error": ref_error,
        "independent_radial_derivative_abs_error": derivative_error,
        "physical_centered_radial_derivative_abs_error": fd_error,
        "exterior_value_and_physical_derivative_zero": True,
        "source": "manufactured Fourier data at paper-admissible neighboring bands",
    }
