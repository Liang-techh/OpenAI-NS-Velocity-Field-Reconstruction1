"""Independent finite-N and endpoint-chain fixture for the long reference seam.

This file deliberately contains only a synthetic source.  It does not import an
accepted producer, report, source archive, or reconstruction owner.  The
long-reference geometry is represented by two native coordinates with unequal
Jacobians, while the finite-N check uses an explicit periodic ``B`` source and
an independent ``p2_yZ`` term.  The four-piece signed chain then exercises the
external long-reshape/postrestore endpoints and three exact internal
cancellations.
"""

from __future__ import annotations

import json

import mpmath as mp


_TOL_TEXT = "1e-38"


def _c(value: str) -> mp.mpf:
    return mp.mpf(value)


def _smooth_step(q):
    """A C-infinity cutoff which is one on the long-reference endpoint."""
    u = (q - _c("0.25")) / _c("0.50")
    if u <= 0:
        return mp.mpf(0)
    if u >= 1:
        return mp.mpf(1)
    a = mp.exp(-1 / u)
    b = mp.exp(-1 / (1 - u))
    return a / (a + b)


def _smooth_step_prime(q):
    """Derivative of :func:`_smooth_step`, including its flat branches."""
    u = (q - _c("0.25")) / _c("0.50")
    if u <= 0 or u >= 1:
        return mp.mpf(0)
    a = mp.exp(-1 / u)
    b = mp.exp(-1 / (1 - u))
    ds_du = a * b * (1 / u**2 + 1 / (1 - u)**2) / (a + b) ** 2
    return ds_du / _c("0.50")


def _slow_B(z):
    """Nonzero slow source amplitude and its ordinary-Z derivative."""
    return _c("0.23") + _c("0.19") * z + _c("0.07") * z * z


def _slow_B_z(z):
    return _c("0.19") + _c("0.14") * z


def _logE_left(q, z):
    # Physical y is T*q.  The cutoff is flat at q=1, so the B term is
    # present in the interior but absent from the seam germ.
    return (_slow_B(z) * (1 - _smooth_step(q))
            - mp.log(1 + z * z)
            + _c("0.1") * _c("2.3") * q - _c("0.19"))


def _logE_right(t, z):
    # Physical y is T + Jref*t.
    return (-mp.log(1 + z * z)
            + _c("0.1") * (_c("2.3") + _c("5.7") * t)
            - _c("0.19"))


def _E_left(q, z):
    return mp.exp(_logE_left(q, z))


def _E_right(t, z):
    return mp.exp(_logE_right(t, z))


def _E_left_y(q, z):
    # d/dy = (1/T)d/dq.
    logarithmic_q_rate = -_slow_B(z) * _smooth_step_prime(q) + _c("0.1") * _c("2.3")
    return _E_left(q, z) * logarithmic_q_rate / _c("2.3")


def _E_right_y(t, z):
    # d/dy = (1/Jref)d/dt.
    return _E_right(t, z) * (_c("0.1") * _c("5.7")) / _c("5.7")


def _E_left_yy_at_seam(z):
    return _c("0.01") * _E_left(mp.mpf(1), z)


def _E_right_yy_at_seam(z):
    return _c("0.01") * _E_right(mp.mpf(0), z)


def _E_left_physical(y, z):
    return _E_left(y / _c("2.3"), z)


def _E_right_physical(y, z):
    return _E_right((y - _c("2.3")) / _c("5.7"), z)


def _A(z):
    return _c("0.31") + _c("0.17") * z + _c("0.06") * z * z


def _A_z(z):
    return _c("0.17") + _c("0.12") * z


def _H(phi, z):
    """Periodic finite-N profile with a genuine slow-Z dependence."""
    return (mp.sin(phi) + _c("0.18") * z * mp.cos(2 * phi)
            + _c("0.04") * z * z * mp.sin(3 * phi))


def _H_phi(phi, z):
    return (mp.cos(phi) - _c("0.36") * z * mp.sin(2 * phi)
            + _c("0.12") * z * z * mp.cos(3 * phi))


def _H_z(phi, z):
    return _c("0.18") * mp.cos(2 * phi) + _c("0.08") * z * mp.sin(3 * phi)


def _H_phi_z(phi, z):
    return -_c("0.36") * mp.sin(2 * phi) + _c("0.24") * z * mp.cos(3 * phi)


def _phase(y, phase=None):
    if phase is None:
        phase = _c("0.413")
    return _c(17) * y + phase


def _B(y, z, phase=None):
    """Synthetic finite-N velocity correction: the primitive is divided by N."""
    phi = _phase(y, phase)
    return _A(z) * mp.exp(-y / 6) * _H(phi, z) / _c(17)


def _B_y(y, z, phase=None):
    phi = _phase(y, phase)
    return (_A(z) * mp.exp(-y / 6)
            * (_H_phi(phi, z) - _H(phi, z) / (_c(6) * _c(17))))


def _B_y_without_fast_phase_term(y, z, phase=None):
    phi = _phase(y, phase)
    return -_A(z) * mp.exp(-y / 6) * _H(phi, z) / (_c(6) * _c(17))


def _B_y_with_extra_fast_division(y, z, phase=None):
    """Wrong control: divides the already-cancelled fast term by N again."""
    phi = _phase(y, phase)
    return (_A(z) * mp.exp(-y / 6)
            * (_H_phi(phi, z) / _c(17)
               - _H(phi, z) / (_c(6) * _c(17))))


def _B_y_z(y, z, phase=None):
    phi = _phase(y, phase)
    return (mp.exp(-y / 6)
            * (_A_z(z) * (_H_phi(phi, z) - _H(phi, z) / (_c(6) * _c(17)))
               + _A(z) * (_H_phi_z(phi, z) - _H_z(phi, z) / (_c(6) * _c(17)))))


def _P(z):
    return _c("0.21") + _c("0.13") * z - _c("0.05") * z * z


def _P_z(z):
    return _c("0.13") - _c("0.10") * z


def _p2(y, z, phase=None):
    phi = _phase(y, phase)
    return _P(z) * mp.exp(-y / 7) * mp.cos(phi) / _c(17)


def _p2_y(y, z, phase=None):
    phi = _phase(y, phase)
    return (_P(z) * mp.exp(-y / 7)
            * (-mp.sin(phi) - mp.cos(phi) / (_c(7) * _c(17))))


def _p2_y_z(y, z, phase=None):
    phi = _phase(y, phase)
    return (_P_z(z) * mp.exp(-y / 7)
            * (-mp.sin(phi) - mp.cos(phi) / (_c(7) * _c(17))))


def _Q_y(y, z, phase=None):
    return _B_y(y, z, phase) + _p2_y(y, z, phase)


def _Q_y_z(y, z, phase=None):
    return _B_y_z(y, z, phase) + _p2_y_z(y, z, phase)


def _leading(y, z):
    return ((_c("0.41") + _c("0.07") * z + _c("0.03") * z * z)
            * mp.exp(-y / 9) * mp.cos(_c("0.23") * y + _c("0.11") * z))


def _leading_z(y, z):
    base = _c("0.41") + _c("0.07") * z + _c("0.03") * z * z
    base_z = _c("0.07") + _c("0.06") * z
    angle = _c("0.23") * y + _c("0.11") * z
    return (mp.exp(-y / 9)
            * (base_z * mp.cos(angle) - _c("0.11") * base * mp.sin(angle)))


def _radius(y):
    return _c("2.1") * mp.exp(y)


def _D(y, z, phase=None):
    return _radius(y) * (_Q_y(y, z, phase) - _leading(y, z))


def _D_z(y, z, phase=None):
    return _radius(y) * (_Q_y_z(y, z, phase) - _leading_z(y, z))


def _D_without_radius(y, z, phase=None):
    return _Q_y(y, z, phase) - _leading(y, z)


def _D_without_radius_z(y, z, phase=None):
    return _Q_y_z(y, z, phase) - _leading_z(y, z)


def _D_without_leading(y, z, phase=None):
    return _radius(y) * _Q_y(y, z, phase)


def _pieces():
    endpoints = (_c("0"), _c("0.8"), _c("1.4"), _c("2.5"), _c("3.1"))
    names = ("long_reshape", "reference", "restoration", "postrestore")
    return tuple({"name": name, "lower": lo, "upper": hi, "J": hi - lo}
                 for name, lo, hi in zip(names, endpoints, endpoints[1:]))


def _max_abs(values):
    return max((abs(v) for v in values), default=mp.mpf(0))


def _run_impl():
    z0 = mp.mpf("0.27")
    z1 = mp.mpf("-0.34")
    seam_y = mp.mpf("2.3")
    tol = mp.mpf(_TOL_TEXT)

    # The two native endpoint descriptions are intentionally differentiated
    # in native coordinates and then divided by their unequal physical J.
    left_value = _E_left(mp.mpf(1), z0)
    right_value = _E_right(mp.mpf(0), z0)
    left_native_y = left_value * mp.mpf("0.1") * mp.mpf("2.3")
    right_native_y = right_value * mp.mpf("0.1") * mp.mpf("5.7")
    left_native_yy = left_value * (mp.mpf("0.1") * mp.mpf("2.3")) ** 2
    right_native_yy = right_value * (mp.mpf("0.1") * mp.mpf("5.7")) ** 2
    E_jets = [
        left_value - right_value,
        left_native_y / mp.mpf("2.3") - right_native_y / mp.mpf("5.7"),
        left_native_yy / mp.mpf("2.3") ** 2 - right_native_yy / mp.mpf("5.7") ** 2,
        mp.diff(lambda zz: _E_left(mp.mpf(1), zz), z0)
        - mp.diff(lambda zz: _E_right(mp.mpf(0), zz), z0),
        mp.diff(lambda zz: _E_left_y(mp.mpf(1), zz), z0)
        - mp.diff(lambda zz: _E_right_y(mp.mpf(0), zz), z0),
    ]
    physical_jets = [
        mp.diff(lambda yy: _E_left_physical(yy, z0), seam_y) - left_native_y / mp.mpf("2.3"),
        mp.diff(lambda yy: _E_right_physical(yy, z0), seam_y) - right_native_y / mp.mpf("5.7"),
        mp.diff(lambda yy: _E_left_physical(yy, z0), seam_y, 2)
        - left_native_yy / mp.mpf("2.3") ** 2,
        mp.diff(lambda yy: _E_right_physical(yy, z0), seam_y, 2)
        - right_native_yy / mp.mpf("5.7") ** 2,
    ]

    # Direct ordinary-Z derivatives of the two physical y jets agree at the
    # seam, while the cutoff makes the whole-long chart law false in its body.
    ordinary_z_rows = [
        mp.diff(lambda zz: _E_left(mp.mpf(1), zz), z0),
        mp.diff(lambda zz: _E_right(mp.mpf(0), zz), z0),
        mp.diff(lambda zz: _E_left_y(mp.mpf(1), zz), z0),
        mp.diff(lambda zz: _E_right_y(mp.mpf(0), zz), z0),
    ]
    inside_q = mp.mpf("0.50")
    whole_long_law_difference = _E_left_y(inside_q, z0) - _c("0.1") * _E_left(inside_q, z0)

    # Independent finite-N physical source checks.  B has a global fast phase
    # and p2_yZ is a separately specified, nonzero slow derivative.
    finite_differences = []
    p2_yz_seen = []
    for piece in _pieces():
        x = mp.mpf("0.37")
        y = piece["lower"] + piece["J"] * x
        J = piece["J"]
        native_B_y = mp.diff(lambda xx, p=piece: _B(p["lower"] + p["J"] * xx, z0), x) / J
        native_p2_y = mp.diff(lambda xx, p=piece: _p2(p["lower"] + p["J"] * xx, z0), x) / J
        finite_differences.extend((native_B_y - _B_y(y, z0), native_p2_y - _p2_y(y, z0)))
        finite_differences.append(
            mp.diff(lambda zz: _p2_y(y, zz), z0) - _p2_y_z(y, z0))
        finite_differences.append(
            mp.diff(lambda zz: _B_y(y, zz), z0) - _B_y_z(y, z0))
        p2_yz_seen.append(abs(_p2_y_z(y, z0)))

    # Four signed sectors use D(lower)-D(upper); all four native spans differ
    # from the long-reference T/J pair above and share physical seam traces.
    pieces = _pieces()
    direct_c0 = sum((_D(piece["lower"], z0) - _D(piece["upper"], z0)
                     for piece in pieces), mp.mpf(0))
    collapsed_c0 = _D(pieces[0]["lower"], z0) - _D(pieces[-1]["upper"], z0)
    direct_c1 = sum((_D_z(piece["lower"], z0) - _D_z(piece["upper"], z0)
                     for piece in pieces), mp.mpf(0))
    collapsed_c1 = _D_z(pieces[0]["lower"], z0) - _D_z(pieces[-1]["upper"], z0)
    direct_z_checks = []
    seam_checks = []
    for piece in pieces:
        for endpoint in (piece["lower"], piece["upper"]):
            direct_z_checks.append(
                _D_z(endpoint, z0) - mp.diff(lambda zz, yy=endpoint: _D(yy, zz), z0))
    for left, right in zip(pieces, pieces[1:]):
        seam_checks.append(_D(left["upper"], z0) - _D(right["lower"], z0))

    # Native derivative conversion on the complete physical D source catches a
    # forgotten J division without assuming equal chart lengths.
    native_D_checks = []
    for piece in pieces:
        x = mp.mpf("0.41")
        native = mp.diff(
            lambda xx, p=piece: _D(p["lower"] + p["J"] * xx, z0), x) / piece["J"]
        physical = mp.diff(lambda yy: _D(yy, z0), piece["lower"] + piece["J"] * x)
        native_D_checks.append(native - physical)

    # A local phase reset at each chart boundary is intentionally wrong.
    reset_chain = mp.mpf(0)
    for piece in pieces:
        reset_phase = mp.mpf("0.413") - mp.mpf(17) * piece["lower"]
        reset_chain += (
            _D(piece["lower"], z0, reset_phase)
            - _D(piece["upper"], z0, reset_phase))

    exterior = {
        "long_reshape_native0": _D(pieces[0]["lower"], z0),
        "postrestore_native1": _D(pieces[-1]["upper"], z0),
    }
    assert abs(exterior["long_reshape_native0"]) > mp.mpf("1e-8")
    assert abs(exterior["postrestore_native1"]) > mp.mpf("1e-8")

    # Negative controls must remain visibly different from the signed source
    # chain; none is interpreted as a source value or a certified cap.
    wrong_sign_c0 = -collapsed_c0
    wrong_sign_c1 = -collapsed_c1
    dropped_inlet_c0 = -exterior["postrestore_native1"]
    dropped_inlet_c1 = -_D_z(pieces[-1]["upper"], z0)
    no_radius_c0 = _D_without_radius(pieces[0]["lower"], z0) - _D_without_radius(pieces[-1]["upper"], z0)
    no_radius_c1 = _D_without_radius_z(pieces[0]["lower"], z0) - _D_without_radius_z(pieces[-1]["upper"], z0)
    no_leading_c0 = _D_without_leading(pieces[0]["lower"], z0) - _D_without_leading(pieces[-1]["upper"], z0)
    no_leading_c1 = mp.diff(lambda zz: _D_without_leading(pieces[0]["lower"], zz), z0) - mp.diff(
        lambda zz: _D_without_leading(pieces[-1]["upper"], zz), z0)
    omitted_p2_yz_c1 = _D_z(pieces[0]["lower"], z0) - _D_z(pieces[-1]["upper"], z0)
    omitted_p2_yz_c1 -= (_radius(pieces[0]["lower"]) * (_B_y_z(pieces[0]["lower"], z0) - _leading_z(pieces[0]["lower"], z0))
                         - _radius(pieces[-1]["upper"]) * (_B_y_z(pieces[-1]["upper"], z0) - _leading_z(pieces[-1]["upper"], z0)))
    no_fast_c0 = (_radius(pieces[0]["lower"]) * (_B_y_without_fast_phase_term(pieces[0]["lower"], z0) + _p2_y(pieces[0]["lower"], z0) - _leading(pieces[0]["lower"], z0))
                  - _radius(pieces[-1]["upper"]) * (_B_y_without_fast_phase_term(pieces[-1]["upper"], z0) + _p2_y(pieces[-1]["upper"], z0) - _leading(pieces[-1]["upper"], z0)))
    extra_fast_division_c0 = (
        _radius(pieces[0]["lower"]) * (
            _B_y_with_extra_fast_division(pieces[0]["lower"], z0)
            + _p2_y(pieces[0]["lower"], z0) - _leading(pieces[0]["lower"], z0))
        - _radius(pieces[-1]["upper"]) * (
            _B_y_with_extra_fast_division(pieces[-1]["upper"], z0)
            + _p2_y(pieces[-1]["upper"], z0) - _leading(pieces[-1]["upper"], z0)))
    wrong_J = right_native_y / mp.mpf("2.3") - right_native_y / mp.mpf("5.7")

    controls = {
        "wrong_sign_rejected": abs(wrong_sign_c0 - collapsed_c0) > mp.mpf("1e-8") and abs(wrong_sign_c1 - collapsed_c1) > mp.mpf("1e-8"),
        "dropped_external_inlet_rejected": abs(dropped_inlet_c0 - collapsed_c0) > mp.mpf("1e-8") and abs(dropped_inlet_c1 - collapsed_c1) > mp.mpf("1e-8"),
        "missing_physical_R_rejected": abs(no_radius_c0 - collapsed_c0) > mp.mpf("1e-8"),
        "missing_leading_subtraction_rejected": abs(no_leading_c0 - collapsed_c0) > mp.mpf("1e-8") and abs(no_leading_c1 - collapsed_c1) > mp.mpf("1e-8"),
        "missing_p2_yZ_rejected": abs(omitted_p2_yz_c1) > mp.mpf("1e-8"),
        "omitted_fast_phase_term_rejected": abs(no_fast_c0 - collapsed_c0) > mp.mpf("1e-8"),
        "extraneous_fast_Bphi_over_N_rejected": abs(extra_fast_division_c0 - collapsed_c0) > mp.mpf("1e-8"),
        "wrong_native_J_rejected": abs(wrong_J) > mp.mpf("1e-8"),
        "phase_endpoint_reset_rejected": abs(reset_chain - collapsed_c0) > mp.mpf("1e-8"),
        "whole_chart_Ey_equals_point_one_E_rejected": abs(whole_long_law_difference) > mp.mpf("1e-8"),
        "caps_not_used_as_source_values": True,
    }
    if not all(controls.values()):
        raise AssertionError(("negative control", controls))

    differences = [
        *E_jets,
        *physical_jets,
        ordinary_z_rows[0] - ordinary_z_rows[1],
        ordinary_z_rows[2] - ordinary_z_rows[3],
        *finite_differences,
        direct_c0 - collapsed_c0,
        direct_c1 - collapsed_c1,
        *direct_z_checks,
        *seam_checks,
        *native_D_checks,
    ]
    maximum = _max_abs(differences)
    if maximum > tol:
        raise AssertionError(("long reference synthetic chain", mp.nstr(maximum, 30)))

    return {
        "passed": True,
        "synthetic_parameters_only": True,
        "original_source_values_evaluated": False,
        "finite_N": 17,
        "long_reference_T": "2.3",
        "reference_native_J": "5.7",
        "arbitrary_global_phase": "0.413",
        "physical_y_seam": "T*q = T + Jref*t",
        "physical_E_two_jets_at_seam": True,
        "ordinary_Z1_E_and_Ey_at_seam": True,
        "whole_long_law_false_inside": True,
        "periodic_B_synthetic": True,
        "B_is_finite_N_velocity_correction": True,
        "B_primitive_divided_by_N": True,
        "genuine_p2_yZ": True,
        "unequal_native_J_divisions": True,
        "native_chain_checks": len(native_D_checks),
        "common_physical_seams": len(seam_checks),
        "both_external_endpoints_nonzero": True,
        "C0_chain_difference": mp.nstr(direct_c0 - collapsed_c0, 30),
        "Z1_chain_difference": mp.nstr(direct_c1 - collapsed_c1, 30),
        "maximum_absolute_difference": mp.nstr(maximum, 30),
        "maximum_p2_yZ_magnitude": mp.nstr(max(p2_yz_seen), 30),
        "controls": controls,
        "negative_controls": controls,
        "signed_trace": {
            "piece_orientation": "D(lower)-D(upper)",
            "delta_definition": "D=R*(Qy_modified-Qy_original_leading)",
            "external_remainder": "D(long_reshape_native0)-D(postrestore_native1=Rm)",
            "internal_seams_cancelled": [
                "long_reshape native1=reference native0",
                "reference native1=restoration native0",
                "restoration native1=postrestore native0",
            ],
            "terminal_scope": "postrestore native1=Rm; no Rm_patch identity claimed",
        },
    }


def run_fixture():
    """Run with temporary high precision and restore the caller context."""
    with mp.workdps(65):
        return _run_impl()


run = run_fixture


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
