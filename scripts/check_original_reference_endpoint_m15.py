"""Independent three-piece finite-N signed endpoint fixture.

The fixture uses one smooth physical source on three native charts.  Each
piece is represented by ``D(lower) - D(upper)``, where

``D = R * (Qy_modified - Qy_leading)``.

The leading trace is deliberately nonzero, the fast phase is global, and the
three native Jacobians differ.  No accepted source graph, numerical report,
or production constructor is imported.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import mpmath as mp


_N = mp.mpf(17)
_PHASE = mp.mpf("0.37")
_DELTA = mp.mpf("0.07")
_R0 = mp.mpf("2.4")
_R_RATE = mp.mpf(1)
_TOL = mp.mpf("1e-38")


def _amplitude(z):
    """Slow nonconstant amplitude and its ordinary-Z derivative."""
    return 1 + mp.mpf("0.2") * z + mp.mpf("0.1") * z * z


def _mean_amplitude(z):
    return _amplitude(z)


def _u(y, z):
    """Finite-N slow/fast velocity source, with a global unwrapped phase."""
    return (_amplitude(z) * mp.exp(-y / 3) *
            mp.sin(_N * y + _PHASE) / _N)


def _m(y, z):
    """Duhamel mean solving m_y = U - m, m(0)=0."""
    a = mp.mpf(2) / 3 + 1j * _N
    integral = mp.im(mp.exp(1j * _PHASE) * (mp.exp(a * y) - 1) / a)
    return _mean_amplitude(z) * mp.exp(-y) * integral / _N


def _u_y(y, z):
    phase = _N * y + _PHASE
    return (_amplitude(z) * mp.exp(-y / 3) *
            (mp.cos(phase) - mp.sin(phase) / (3 * _N)))


def _u_z(y, z):
    return (mp.mpf("0.2") + mp.mpf("0.2") * z) * mp.exp(-y / 3) * \
        mp.sin(_N * y + _PHASE) / _N


def _m_z(y, z):
    a = mp.mpf(2) / 3 + 1j * _N
    integral = mp.im(mp.exp(1j * _PHASE) * (mp.exp(a * y) - 1) / a)
    return (mp.mpf("0.2") + mp.mpf("0.2") * z) * mp.exp(-y) * integral / _N


def _qy_formula(y, z):
    """Physical Q_y formula using the independently specified U and m jets."""
    u = _u(y, z)
    m = _m(y, z)
    numerator = (2 * z * _u_y(y, z)
                 - z * (1 - _DELTA) * (u - m)
                 - (1 - z * z) * (_u_z(y, z) - _m_z(y, z)))
    return numerator / (1 - _DELTA * z * z)


def _qy_with_phase_offset(y, z, phase_offset):
    """Q_y with an explicitly supplied global phase offset.

    This is used only by the negative phase-reset control.  The production
    fixture always calls :func:`_qy_formula`, whose phase is global.
    """
    phase = _N * y + phase_offset
    amplitude = _amplitude(z)
    amplitude_z = mp.mpf("0.2") + mp.mpf("0.2") * z
    u = amplitude * mp.exp(-y / 3) * mp.sin(phase) / _N
    uy = amplitude * mp.exp(-y / 3) * (mp.cos(phase) - mp.sin(phase) / (3 * _N))
    uz = amplitude_z * mp.exp(-y / 3) * mp.sin(phase) / _N
    a = mp.mpf(2) / 3 + 1j * _N
    integral = mp.im(mp.exp(1j * phase_offset) * (mp.exp(a * y) - 1) / a)
    mean = amplitude * mp.exp(-y) * integral / _N
    mean_z = amplitude_z * mp.exp(-y) * integral / _N
    numerator = (2 * z * uy
                 - z * (1 - _DELTA) * (u - mean)
                 - (1 - z * z) * (uz - mean_z))
    return numerator / (1 - _DELTA * z * z)


def _qy_direct(y, z):
    """Recompute Q_y by differentiating the physical Q in y."""
    u = _u(y, z)
    m = _m(y, z)
    uz = mp.diff(lambda zz: _u(y, zz), z)
    mz = mp.diff(lambda zz: _m(y, zz), z)
    q = (2 * z * u - (1 - _DELTA) * z * m - (1 - z * z) * mz) / \
        (1 - _DELTA * z * z)

    # Rebuild Q at the perturbed physical coordinate so that the y derivative
    # also differentiates the independently computed m_Z term.
    def q_at(yy):
        uu = _u(yy, z)
        mm = _m(yy, z)
        mmz = mp.diff(lambda zz: _m(yy, zz), z)
        return (2 * z * uu - (1 - _DELTA) * z * mm -
                (1 - z * z) * mmz) / (1 - _DELTA * z * z)

    # Keep the local Q construction live in the direct path; the assignment
    # above is a sanity check that the unperturbed expression is well-defined.
    if not mp.isfinite(q):
        raise ArithmeticError("synthetic Q became non-finite")
    return mp.diff(q_at, y)


def _q_leading(y, z):
    """A nonzero leading trace retained inside every physical subtraction."""
    return (mp.mpf("0.6") + mp.mpf("0.08") * z + mp.mpf("0.03") * z * z
            + mp.mpf("0.12") * mp.exp(-mp.mpf("0.2") * y)
            * (1 + mp.mpf("0.1") * z)
            * mp.cos(mp.mpf("0.4") * y + mp.mpf("0.2") * z))


def _radius(y):
    return _R0 * mp.exp(_R_RATE * y)


def _D(y, z):
    return _radius(y) * (_qy_direct(y, z) - _q_leading(y, z))


def _D_formula(y, z):
    return _radius(y) * (_qy_formula(y, z) - _q_leading(y, z))


def _D_without_radius(y, z):
    return _qy_direct(y, z) - _q_leading(y, z)


def _D_without_leading(y, z):
    return _radius(y) * _qy_direct(y, z)


def _pieces():
    # Native x runs from 0 to 1 on every chart, while physical y advances by
    # a different positive Jacobian.  Seams are shared in physical y.
    return (
        {"name": "reference", "lower": mp.mpf("0"), "upper": mp.mpf("1.25"), "J": mp.mpf("1.25")},
        {"name": "restoration", "lower": mp.mpf("1.25"), "upper": mp.mpf("2.0"), "J": mp.mpf("0.75")},
        {"name": "postrestore", "lower": mp.mpf("2.0"), "upper": mp.mpf("3.5"), "J": mp.mpf("1.5")},
    )


def _chart_value(piece, x, z, function=_D):
    return function(piece["lower"] + piece["J"] * x, z)


def _sector(piece, z, function=_D):
    return _chart_value(piece, 0, z, function) - _chart_value(piece, 1, z, function)


def _chain(z, function=_D):
    return sum((_sector(piece, z, function) for piece in _pieces()), mp.mpf(0))


def _collapsed(z, function=_D):
    pieces = _pieces()
    return function(pieces[0]["lower"], z) - function(pieces[-1]["upper"], z)


def _max_abs(values):
    return max((abs(value) for value in values), default=mp.mpf(0))


def _run_impl():
    """Run the bounded source-independent C0/C1 endpoint check."""
    ctx = mp.mp
    # Keep all calculations in this context, including constants used above.
    # The functions read mpmath's active context and therefore remain free of
    # report/source imports.
    z0 = ctx.mpf("0.23")
    z1 = ctx.mpf("-0.31")
    pieces = _pieces()

    # Verify direct physical differentiation of Q_y before using it in the
    # signed chain.  This catches an accidentally copied stored derivative.
    qy_differences = []
    mean_ode_differences = []
    for y in (ctx.mpf("0.17"), pieces[0]["upper"], pieces[1]["upper"], pieces[2]["upper"]):
        for z in (z0, z1):
            qy_differences.append(_qy_direct(y, z) - _qy_formula(y, z))
            mean_ode_differences.append(
                mp.diff(lambda yy: _m(yy, z), y) - (_u(y, z) - _m(y, z))
            )

    # The direct composite evaluates each source piece in native coordinates;
    # the collapsed form uses only the two exterior physical endpoints.
    direct_c0 = _chain(z0)
    collapsed_c0 = _collapsed(z0)
    direct_c1 = mp.diff(_chain, z0)
    collapsed_c1 = mp.diff(_collapsed, z0)
    endpoint_c1 = sum((
        mp.diff(lambda zz, p=piece: _D(p["lower"], zz), z0)
        - mp.diff(lambda zz, p=piece: _D(p["upper"], zz), z0)
        for piece in pieces), ctx.mpf(0))

    # Native derivative / physical derivative conversion is checked at an
    # interior point on every chart, using three genuinely unequal J values.
    native_derivative_differences = []
    for piece in pieces:
        x = ctx.mpf("0.37")
        y = piece["lower"] + piece["J"] * x
        native = mp.diff(lambda xx, p=piece: _D_formula(
            p["lower"] + p["J"] * xx, z0), x)
        physical = mp.diff(lambda yy: _D_formula(yy, z0), y)
        native_derivative_differences.append(native / piece["J"] - physical)

    # Common seam values use the global fast phase N*y + phase.  A local phase
    # reset would change the second and third chart endpoint traces.
    seam_differences = []
    phase_reset_differences = []
    for left, right in zip(pieces, pieces[1:]):
        seam_differences.append(_D(left["upper"], z0) - _D(right["lower"], z0))
        reset_offset = _PHASE - _N * right["lower"]
        phase_reset_differences.append(
            _radius(right["lower"]) * (
                _qy_with_phase_offset(right["lower"], z0, reset_offset)
                - _q_leading(right["lower"], z0)
            ) - _D(left["upper"], z0)
        )

    # The exterior traces are nonzero at both tested slow coordinates, so
    # dropping the inlet cannot pass accidentally.
    exterior = {
        "reference_native0": _D(pieces[0]["lower"], z0),
        "postrestore_native1": _D(pieces[-1]["upper"], z0),
    }
    _require_nonzero = lambda value: abs(value) > ctx.mpf("1e-8")
    assert _require_nonzero(exterior["reference_native0"])
    assert _require_nonzero(exterior["postrestore_native1"])

    # Negative controls each alter exactly one physical rule.
    wrong_sign_c0 = -collapsed_c0
    wrong_sign_c1 = -collapsed_c1
    dropped_inlet_c0 = -_D(pieces[-1]["upper"], z0)
    dropped_inlet_c1 = -mp.diff(lambda zz: _D(pieces[-1]["upper"], zz), z0)
    no_radius_c0 = _collapsed(z0, _D_without_radius)
    no_radius_c1 = mp.diff(lambda zz: _collapsed(zz, _D_without_radius), z0)
    no_leading_c0 = _collapsed(z0, _D_without_leading)
    no_leading_c1 = mp.diff(lambda zz: _collapsed(zz, _D_without_leading), z0)

    # A positive enclosure is metadata, not a signed source value.  A toy
    # duplicate-center calculation must therefore disagree with the exact
    # source chain and is recorded as rejected rather than used as a bound.
    cap_ref = abs(exterior["reference_native0"]) + ctx.mpf("0.75")
    cap_rm = abs(exterior["postrestore_native1"]) + ctx.mpf("0.75")
    duplicate_cap_center = cap_ref + cap_ref - cap_rm

    differences = [
        *qy_differences,
        *mean_ode_differences,
        direct_c0 - collapsed_c0,
        direct_c1 - collapsed_c1,
        direct_c1 - endpoint_c1,
        *native_derivative_differences,
        *seam_differences,
    ]
    maximum = _max_abs(differences)
    controls = {
        "wrong_sign_rejected": abs(wrong_sign_c0 - collapsed_c0) > ctx.mpf("1e-8") and
        abs(wrong_sign_c1 - collapsed_c1) > ctx.mpf("1e-8"),
        "dropped_external_inlet_rejected": abs(dropped_inlet_c0 - collapsed_c0) > ctx.mpf("1e-8") and
        abs(dropped_inlet_c1 - collapsed_c1) > ctx.mpf("1e-8"),
        "missing_physical_R_rejected": abs(no_radius_c0 - collapsed_c0) > ctx.mpf("1e-8") and
        abs(no_radius_c1 - collapsed_c1) > ctx.mpf("1e-8"),
        "missing_leading_subtraction_rejected": abs(no_leading_c0 - collapsed_c0) > ctx.mpf("1e-8") and
        abs(no_leading_c1 - collapsed_c1) > ctx.mpf("1e-8"),
        "duplicate_cap_as_value_rejected": abs(duplicate_cap_center - collapsed_c0) > ctx.mpf("1e-8"),
        "caps_not_used_as_source_values": True,
        "global_fast_phase_reset_rejected": _max_abs(phase_reset_differences) > ctx.mpf("1e-8"),
    }
    if not all(controls.values()):
        raise AssertionError(("negative control", controls))
    if maximum > _TOL:
        raise AssertionError(("three-piece signed chain", mp.nstr(maximum, 30)))

    return {
        "passed": True,
        "synthetic_parameters_only": True,
        "original_source_values_evaluated": False,
        "finite_N": 17,
        "finite_N_slow_fast": True,
        "native_Jacobians": [str(piece["J"]) for piece in pieces],
        "unequal_J_naturality": True,
        "common_physical_seams": 2,
        "nonzero_leading_subtraction": True,
        "physical_R_once": True,
        "physical_y_and_ordinary_Z1_independent_differentiation": True,
        "both_external_endpoints_nonzero": True,
        "C0_chain_difference": mp.nstr(direct_c0 - collapsed_c0, 30),
        "Z1_chain_difference": mp.nstr(direct_c1 - collapsed_c1, 30),
        "maximum_absolute_difference": mp.nstr(maximum, 30),
        "controls": controls,
        "negative_controls": controls,
        "signed_trace": {
            "piece_orientation": "D(lower)-D(upper)",
            "delta_definition": "D=R*(Qy_modified-Qy_leading)",
            "external_remainder": "D(reference_native0)-D(postrestore_native1)",
            "internal_seams_cancelled": [
                "reference_native1=restoration_native0",
                "restoration_native1=postrestore_native0",
            ],
        },
    }


def run():
    """Run the fixture while restoring the caller's mpmath precision."""
    with mp.workdps(55):
        return _run_impl()


# The production runner uses the explicit fixture name; keep ``run`` as the
# standalone public entry point as well.
run_fixture = run


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
