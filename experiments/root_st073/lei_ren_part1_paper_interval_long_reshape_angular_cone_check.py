"""Independent scalar check of the cancellation-safe long-shape angular stress.

The fixture uses a moderate radius and T so direct MP quadrature can resolve
the endpoint-normalized kernels.  It compares ``normalized_angular`` with
the original scalar stress evaluator and repeats the comparison after a
common angular-amplitude rescaling.  This is a local algebra check only;
the production cone receipt remains relaxed-branch and source-conditional.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_long_reshape_angular_cone import normalized_angular
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


HERE = Path(__file__).parent


def _contains(box, value):
    lo, hi = endpoints(box)
    return lo <= value <= hi


def _record(value):
    lo, hi = endpoints(value)
    return [mp.nstr(lo, 80), mp.nstr(hi, 80)]


def _sigma(value):
    if value <= 0:
        return mp.mpf(0)
    if value >= 1:
        return mp.mpf(1)
    phase = 1 / (value * value) - 1 / ((1 - value) * (1 - value))
    if phase >= 0:
        e = mp.exp(-phase)
        return e / (1 + e)
    e = mp.exp(phase)
    return 1 / (1 + e)


def _sigma_prime(value):
    if value <= 0 or value >= 1:
        return mp.mpf(0)
    sig = _sigma(value)
    return sig * (1 - sig) * (2 / value**3 + 2 / (1 - value) ** 3)


def _quad(function, cuts):
    return mp.fsum(mp.quad(function, [left, right]) for left, right in zip(cuts, cuts[1:]))


def _kernel(k, m, B, BZ, T, y):
    endpoint = _sigma(y / T)

    def value(ell):
        ds = endpoint - _sigma((y - ell) / T)
        return mp.exp(-k * ell + m * B * ds)

    def tangent(ell):
        ds = endpoint - _sigma((y - ell) / T)
        return m * BZ * ds * value(ell)

    cuts = [mp.mpf(0), y]
    if 0 < y - T < y:
        cuts.insert(1, y - T)
    return _quad(value, cuts), _quad(tangent, cuts)


def _stress_fixture(scale):
    """Build physical moments and scalar stress for one common amplitude."""

    R0 = mp.mpf("110")
    y = mp.mpf("3")
    R = R0 * mp.exp(y)
    Z = mp.mpf(".5")
    delta = mp.mpf(".01")
    T = mp.mpf("10")
    B = mp.mpf("-2")
    BZ = mp.mpf(".3")
    V = mp.mpf("1.1")
    VZ = mp.mpf(".2")
    mass0 = mp.mpf("1.23")
    mass0Z = mp.mpf(".17")
    u0 = mp.mpf("2") * scale
    u0Z = -mp.mpf(".3") * scale
    endpoint_sigma = _sigma(y / T)
    endpoint_prime = _sigma_prime(y / T)
    u = u0 * mp.exp(y / 10 - B * endpoint_sigma)
    zeta = u0Z / u0 - BZ * endpoint_sigma
    uZ = u * zeta
    uy = u * (mp.mpf(".1") - B * endpoint_prime / T)
    st = (2 * uy - u) / u
    r0 = mp.sqrt(2 * R0)
    a0 = r0 * u0
    ctheta0 = R0 * a0
    ntheta = mp.mpf(".35")
    nthetaZ = mp.mpf(".04")
    nmixed = mp.mpf(".12")
    nmixedZ = -mp.mpf(".02")
    source = dict(theta=ctheta0 * ntheta, thetaZ=ctheta0 * nthetaZ,
                  mixed=ctheta0 * nmixed, mixedZ=ctheta0 * nmixedZ)
    J, JZ = _kernel(mp.mpf("1.6"), 1, B, BZ, T, y)
    Jp, JpZ = _kernel(mp.mpf(".2"), 2, B, BZ, T, y)
    Je, JeZ = _kernel(mp.mpf("1.2"), 2, B, BZ, T, y)
    zeta_source = zeta
    ar = r0 * u0
    a = mp.sqrt(2 * R) * u
    dR = R - R0
    moments = {
        "z": R0 * mass0 + V * dR,
        "theta": source["theta"] + R * a * J,
        "theta_z": source["mixed"] + V * R * a * J,
        "z_theta": R0 * (V * V - u0 * u0 / 2)
            + V * V * dR - R * u * u * Je / 2,
        "p": -mp.mpf(".7") + u * u * Jp / 2,
    }
    moments_Z = {
        "z": R0 * mass0Z + VZ * dR,
        "theta": source["thetaZ"] + R * a * (zeta_source * J + JZ),
        "theta_z": source["mixedZ"] + R * a * (VZ * J + V * (zeta_source * J + JZ)),
        "z_theta": R0 * (2 * V * VZ - u0 * u0Z)
            + 2 * V * VZ * dR - R * u * u * (2 * zeta_source * Je + JeZ) / 2,
        "p": mp.mpf(".13") + u * u * (2 * zeta_source * Jp + JpZ) / 2,
    }
    pressure = moments["p"]
    pressure_Z = moments_Z["p"]
    scalar = evaluate_mp_stress(
        mp.log(R), Z, delta, Utheta=u, Uz=V, Utheta_y=uy,
        Utheta_Z=uZ, Uz_y=mp.mpf(0), Uz_Z=VZ, moments=moments,
        moments_Z=moments_Z, P=pressure, P_Z=pressure_Z,
        precision=130,
    )
    c = MPIntervalContext()
    c.dps = 95
    z = IntervalTaylor(c, [c.mpf(Z), c.mpf(0)])
    Vjet = IntervalTaylor(c, [c.mpf(V), c.mpf(VZ)])
    mass_source = IntervalTaylor(c, [c.mpf(mass0), c.mpf(mass0Z)])
    source_n = {"theta": c.mpf(ntheta), "theta_Z": c.mpf(nthetaZ),
                "mixed": c.mpf(nmixed), "mixed_Z": c.mpf(nmixedZ)}
    Jjet = IntervalTaylor(c, [c.mpf(J), c.mpf(JZ)])
    angular = normalized_angular(
        c, z, c.mpf(delta), c.mpf(R), Vjet, mass_source,
        c.mpf(mp.exp(-y)), source_n,
        c.mpf(mp.exp(-mp.mpf("1.6") * y + B * endpoint_sigma)),
        c.mpf(zeta), Jjet, c.mpf(st),
    )
    F = u / mp.sqrt(2 * R)
    direct_i = scalar["I_theta"] / F
    direct_t = scalar["T_theta"] / F
    if not (_contains(angular["I_theta_over_F"], direct_i) and
            _contains(angular["T_theta_over_F"], direct_t)):
        raise AssertionError("normalized angular formula escaped scalar stress")
    return dict(scale=mp.nstr(scale, 40), direct_I_theta_over_F=mp.nstr(direct_i, 80),
                direct_T_theta_over_F=mp.nstr(direct_t, 80),
                interval_I_theta_over_F=_record(angular["I_theta_over_F"]),
                interval_T_theta_over_F=_record(angular["T_theta_over_F"]),
                source_combo=_record(angular["source_angular_combo"]),
                increment_combo=_record(angular["increment_angular_combo"]),
                S_theta_over_F=mp.nstr(st, 80),
                transport_over_R=_record(angular["transport_over_R"]),
                cone_scope="relaxed angular branch only; no axial T_z supplied")


def run():
    with mp.workdps(160):
        rows = [_stress_fixture(mp.mpf("1")),
                _stress_fixture(mp.mpf("1e-25")),
                _stress_fixture(mp.mpf("1e25"))]
        base_i = mp.mpf(rows[0]["direct_I_theta_over_F"])
        base_t = mp.mpf(rows[0]["direct_T_theta_over_F"])
        for row in rows[1:]:
            if abs(mp.mpf(row["direct_I_theta_over_F"]) - base_i) > mp.mpf("1e-70"):
                raise AssertionError("I_theta/F changed under common angular amplitude scaling")
            if abs(mp.mpf(row["direct_T_theta_over_F"]) - base_t) > mp.mpf("1e-70"):
                raise AssertionError("T_theta/F changed under common angular amplitude scaling")
        parent_receipt = HERE / "lei_ren_part1_paper_interval_long_reshape_angular_cone.json"
        raw = json.loads(parent_receipt.read_text(encoding="utf-8"))
        cells = raw.get("cells", [])
        if len(cells) != 2:
            raise AssertionError("expected two whole-shape angular cone cells")
        scoped = []
        for cell in cells:
            cone = cell["cone"]
            if cone.get("branch") != "kappa<=2":
                raise AssertionError("whole-cell certificate unexpectedly claims another branch")
            if cone.get("admissible_cone_certified"):
                raise AssertionError("relaxed cell was mislabeled as strongly admissible")
            scoped.append(dict(status=cone.get("status"), branch=cone.get("branch"),
                               relaxed=cone.get("relaxed_cone_certified"),
                               strong=cone.get("admissible_cone_certified")))
        dependencies = (
            Path(__file__).name,
            "lei_ren_part1_paper_interval_long_reshape_angular_cone.py",
            "lei_ren_part1_paper_interval_long_reshape_field.py",
            "lei_ren_part1_paper_mp_stress.py",
            parent_receipt.name,
        )
        report = dict(passed=True, fixture_only=True, direct_quadrature_precision=130,
                      normalized_angular_scalar_containment=len(rows),
                      common_amplitude_scaling_invariant=True, scales=["1", "1e-25", "1e25"],
                      rows=rows, whole_cell_scope=scoped,
                      strong_admissibility_claimed=False,
                      axial_Tz_evaluated=False, pressure_source_replaced=False,
                      input_hashes={name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                    for name in dependencies})
        output = HERE / "lei_ren_part1_paper_interval_long_reshape_angular_cone_check.json"
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print("Normalized angular scalar fixture passed:", len(rows),
              "amplitude scales; whole-cell scope retained as relaxed-only", flush=True)
        return report


if __name__ == "__main__":
    run()
