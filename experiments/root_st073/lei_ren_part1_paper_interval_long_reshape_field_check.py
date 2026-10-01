"""Independent moderate-scale checks for the directed long-reshape kernels.

This is a fixture only.  It resolves the production kernel on a small
``T=10`` interval by high precision quadrature and checks the directed
negative-``B`` bounds, the first ``Z`` derivative, the ``y=T`` flat-kernel
identity, and the downstream physical packet on a synthetic field.  It does
not certify the very-large-``T`` source or any cone.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_interval_functional_defects import flat_kernel
from lei_ren_part1_paper_interval_long_reshape_field import (
    physical_packet,
    positive_decay_integral,
    reshape_kernel,
    sigma_value_derivative,
)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


HERE = Path(__file__).parent
PRECISION = 70
QUAD_PRECISION = 100
with mp.workdps(140):
    T_VALUE = mp.mpf("10")
    B_VALUE = mp.mpf("-2")
    BZ_VALUE = mp.mpf(".3")
    Z_VALUE = mp.mpf(".5")
    R0_VALUE = mp.mpf("1")
    U0_VALUE = mp.mpf("2")
    V0_VALUE = mp.mpf("1.1")
    VZ_VALUE = mp.mpf(".2")


def _record(value):
    lo, hi = endpoints(value)
    return [mp.nstr(lo, 90), mp.nstr(hi, 90)]


def _contains(box, value):
    lo, hi = endpoints(box)
    return lo <= value <= hi


def _sigma_prime(value):
    if value <= 0 or value >= 1:
        return mp.mpf("0")
    sig = _sigma_mp(value)
    return sig * (1 - sig) * (2 / value**3 + 2 / (1 - value) ** 3)


def _quad(function, points):
    return mp.fsum(mp.quad(function, [left, right]) for left, right in zip(points, points[1:]))


def _direct_kernel(k, m, B, BZ, T, y):
    """Return J and J_Z for the endpoint-normalized kernel."""

    endpoint_sigma = _sigma_mp(y / T)

    def integrand(ell):
        delta_sigma = endpoint_sigma - _sigma_mp((y - ell) / T)
        return mp.exp(-k * ell + m * B * delta_sigma)

    def derivative(ell):
        delta_sigma = endpoint_sigma - _sigma_mp((y - ell) / T)
        return m * BZ * delta_sigma * integrand(ell)

    cuts = [mp.mpf("0")]
    if 0 < y - T < y:
        cuts.append(y - T)
    cuts.append(y)
    return _quad(integrand, cuts), _quad(derivative, cuts)


def _direct_moments(y, T):
    """Resolve five physical moment increments and their Z derivatives."""

    B = B_VALUE
    BZ = BZ_VALUE
    V = V0_VALUE
    VZ = VZ_VALUE
    R = R0_VALUE * mp.exp(y)
    endpoint_sigma = _sigma_mp(y / T)
    endpoint_u = U0_VALUE * mp.exp(y / 10 - B * endpoint_sigma)

    def terms(ell):
        r = R * mp.exp(-ell)
        source_sigma = _sigma_mp((y - ell) / T)
        source_u = endpoint_u * mp.exp(-ell / 10 + B * (endpoint_sigma - source_sigma))
        source_uz = -BZ * source_sigma * source_u
        theta = mp.sqrt(2) * r ** mp.mpf("1.5") * source_u
        theta_z = theta * (-BZ * source_sigma)
        z = V * r
        z_z = VZ * r
        ztheta = V * V * r - source_u * source_u * r / 2
        ztheta_z = 2 * V * VZ * r - source_u * source_uz * r
        pressure = source_u * source_u / 2
        pressure_z = source_u * source_uz
        return (z, theta, V * theta, ztheta, pressure), (z_z, theta_z,
            VZ * theta + V * theta_z, ztheta_z, pressure_z)

    cuts = [mp.mpf("0")]
    if 0 < y - T < y:
        cuts.append(y - T)
    cuts.append(y)
    values = []
    derivatives = []
    for index in range(5):
        values.append(_quad(lambda ell, i=index: terms(ell)[0][i], cuts))
        derivatives.append(_quad(lambda ell, i=index: terms(ell)[1][i], cuts))
    return dict(zip(("z", "theta", "theta_z", "z_theta", "p"), values)), dict(
        zip(("z", "theta", "theta_z", "z_theta", "p"), derivatives)
    )


def _kernel_interval(c, k, m, B, T, y, flat):
    if y <= T:
        return reshape_kernel(c, c.mpf(k), m, B, T, y)
    correction = flat_kernel(c, c.mpf(k), m, B, T)["jet"]
    base = IntervalTaylor(c, [positive_decay_integral(c, c.mpf(k), y), c.mpf(0)])
    return correction * c.exp(-c.mpf(k) * (y - T)) + base


def _moment_intervals(c, B, z, y, T):
    R = c.exp(c.mpf(y))
    endpoint_sigma = _sigma_mp(mp.mpf(y) / T_VALUE)
    endpoint_prime = _sigma_prime(mp.mpf(y) / T_VALUE) / T_VALUE
    u_value = U0_VALUE * mp.exp(mp.mpf(y) / 10 - B_VALUE * endpoint_sigma)
    u_z = u_value * (-BZ_VALUE * endpoint_sigma)
    u_y = u_value * (mp.mpf(".1") - B_VALUE * endpoint_prime)
    u_y_z = u_z * (mp.mpf(".1") - B_VALUE * endpoint_prime) - u_value * BZ_VALUE * endpoint_prime
    u = IntervalTaylor(c, [c.mpf(u_value), c.mpf(u_z)])
    uy = IntervalTaylor(c, [c.mpf(u_y), c.mpf(u_y_z)])
    V = IntervalTaylor(c, [c.mpf(V0_VALUE), c.mpf(VZ_VALUE)])
    Bjet = IntervalTaylor(c, [c.mpf(B_VALUE), c.mpf(BZ_VALUE)])
    theta_kernel = _kernel_interval(c, "1.6", 1, Bjet, c.mpf(T_VALUE), c.mpf(y), flat=True)
    pressure_kernel = _kernel_interval(c, ".2", 2, Bjet, c.mpf(T_VALUE), c.mpf(y), flat=True)
    energy_kernel = _kernel_interval(c, "1.2", 2, Bjet, c.mpf(T_VALUE), c.mpf(y), flat=True)
    theta = u * (c.sqrt(2) * R ** c.mpf("1.5")) * theta_kernel
    pressure = u * u * pressure_kernel / 2
    swirl = u * u * R * energy_kernel / 2
    delta_r = R - c.mpf(R0_VALUE)
    moments = dict(
        z=V * delta_r,
        theta=theta,
        theta_z=V * theta,
        z_theta=V * V * delta_r - swirl,
        p=pressure,
    )
    return moments, u, uy, V, R


def run():
    c = MPIntervalContext()
    c.dps = PRECISION
    report = dict(
        passed=False,
        fixture_only=True,
        production_source_scope=False,
        whole_path_cone_claimed=False,
        T=str(T_VALUE),
        B=str(B_VALUE),
        B_Z=str(BZ_VALUE),
        kernel_cases=[],
        cutoff_cases=[],
        physical_cases=[],
    )

    with mp.workdps(QUAD_PRECISION):
        Bjet = IntervalTaylor(c, [c.mpf(B_VALUE), c.mpf(BZ_VALUE)])
        for s in (mp.mpf(".3"), mp.mpf(".5"), mp.mpf(".7")):
            sig_box, derivative_box = sigma_value_derivative(c, c.mpf(s))
            sig = _sigma_mp(s)
            derivative = _sigma_prime(s)
            if not (_contains(sig_box, sig) and _contains(derivative_box, derivative)):
                raise AssertionError("cutoff value or derivative escaped directed box")
            report["cutoff_cases"].append(dict(
                s=mp.nstr(s, 30), sigma=_record(sig_box), sigma_prime=_record(derivative_box),
                direct_sigma=mp.nstr(sig, 80), direct_sigma_prime=mp.nstr(derivative, 80),
            ))

        for k, m in (("1.6", 1), (".2", 2), ("1.2", 2)):
            for y in (mp.mpf("3"), T_VALUE, mp.mpf("15")):
                direct, direct_z = _direct_kernel(mp.mpf(k), mp.mpf(m), B_VALUE, BZ_VALUE, T_VALUE, y)
                if y <= T_VALUE:
                    box = reshape_kernel(c, c.mpf(k), m, Bjet, c.mpf(T_VALUE), c.mpf(y))
                else:
                    box = _kernel_interval(c, k, m, Bjet, c.mpf(T_VALUE), c.mpf(y), flat=True)
                if not (_contains(box[0], direct) and _contains(box[1], direct_z)):
                    raise AssertionError(f"kernel enclosure failed for k={k}, m={m}, y={y}")
                case = dict(k=k, m=m, y=mp.nstr(y, 30), interval_value=_record(box[0]),
                    interval_Z=_record(box[1]), direct_value=mp.nstr(direct, 80),
                    direct_Z=mp.nstr(direct_z, 80))
                if y == T_VALUE:
                    flat = flat_kernel(c, c.mpf(k), m, Bjet, c.mpf(T_VALUE))["jet"]
                    identity = flat + IntervalTaylor(
                        c, [positive_decay_integral(c, c.mpf(k), c.mpf(T_VALUE)), c.mpf(0)]
                    )
                    if not (_contains(identity[0], direct) and _contains(identity[1], direct_z)):
                        raise AssertionError("y=T flat-kernel identity escaped direct quadrature")
                    case["flat_identity_value"] = _record(identity[0])
                    case["flat_identity_Z"] = _record(identity[1])
                report["kernel_cases"].append(case)

        # A synthetic downstream packet checks that the same five moment
        # increments reach pressure and stress assembly without a midpoint
        # projection.  The direct quadrature is analytic in Z, so no nested
        # numerical differentiation is used.
        zjet = IntervalTaylor(c, [c.mpf(Z_VALUE), c.mpf(1)])
        delta = c.mpf("1e-4")
        # One resolved physical endpoint supplies five value and five Z
        # derivative moment checks while keeping this fixture bounded.
        for y in (mp.mpf("3"),):
            direct_values, direct_derivatives = _direct_moments(y, T_VALUE)
            moments, u, uy, V, R = _moment_intervals(c, Bjet, zjet, y, T_VALUE)
            for name in ("z", "theta", "theta_z", "z_theta", "p"):
                if not (_contains(moments[name][0], direct_values[name]) and
                        _contains(moments[name][1], direct_derivatives[name])):
                    raise AssertionError(f"physical moment enclosure failed for {name}, y={y}")
            packet = physical_packet(c, zjet, delta, R, u, uy, V, moments,
                                     IntervalTaylor(c, [c.mpf(0), c.mpf(0)]))
            if not (_contains(packet["P"][0], direct_values["p"]) and
                    _contains(packet["P"][1], direct_derivatives["p"])):
                raise AssertionError("pressure packet did not retain the physical pressure increment")
            shape_st = mp.mpf(".8") + 2 * B_VALUE * _sigma_prime(y / T_VALUE) / T_VALUE
            report["physical_cases"].append(dict(
                y=mp.nstr(y, 30),
                moment_intervals={name: [_record(moments[name][0]), _record(moments[name][1])]
                    for name in ("z", "theta", "theta_z", "z_theta", "p")},
                pressure_interval=[_record(packet["P"][0]), _record(packet["P"][1])],
                stress_keys=sorted(packet["stress"]),
                normalized_stress_keys=sorted(packet["normalized_stress"]),
                prescribed_shape_coefficient=mp.nstr(shape_st, 80),
                cone_status=packet["cone"]["status"],
            ))

    dependencies = (
        Path(__file__).name,
        "lei_ren_part1_paper_interval_long_reshape_field.py",
        "lei_ren_part1_paper_interval_functional_defects.py",
        "lei_ren_part1_paper_interval_taylor.py",
    )
    report["input_hashes"] = {
        name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in dependencies
    }
    report["passed"] = True
    out = HERE / "lei_ren_part1_paper_interval_long_reshape_field_check.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("Long-reshape kernel and physical packet fixture passed:",
          len(report["kernel_cases"]), "kernels;",
          len(report["physical_cases"]), "physical cases", flush=True)
    return report


if __name__ == "__main__":
    run()
