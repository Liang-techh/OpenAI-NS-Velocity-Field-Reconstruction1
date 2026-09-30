"""Bounded multi-time geometry diagnostic for the regular corrected core.

This is deliberately a local chart experiment.  It builds the existing joined
field, samples one point in the regular core, maps that sample through the
actual physical field, and then inverts the physical callable at four extreme
q scales.  The receipt keeps values in signed-log form because delta is
1e-200.  It does not assert a global flow, a finite-energy construction, or a
recursive transition.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_exit_field import physical_chart
from lei_ren_part1_paper_joined_outer import build_joined_field


OUT_JSON = Path(__file__).with_suffix(".json")


def _text(value: mp.mpf, digits: int) -> str:
    """Serialize an MP value without converting it to binary float."""

    return mp.nstr(mp.mpf(value), digits)


def signed_log(value: mp.mpf, digits: int) -> dict[str, object]:
    """Return sign, log absolute value, and a readable MP value."""

    value = mp.mpf(value)
    if value == 0:
        return {
            "sign": 0,
            "log_abs": None,
            "arbitrary_exponent_value": "0",
        }
    sign = 1 if value > 0 else -1
    return {
        "sign": sign,
        "log_abs": _text(mp.log(abs(value)), digits),
        "arbitrary_exponent_value": _text(value, digits),
    }


def _error(value: mp.mpf, digits: int) -> dict[str, object]:
    """Signed-log representation for a roundtrip error."""

    return signed_log(abs(mp.mpf(value)), digits)


def run() -> dict[str, object]:
    # The joined object is intentional: it obtains the same corrected inner
    # candidate used by the outer prototype rather than rebuilding a toy core.
    print("building joined field", flush=True)
    joined = build_joined_field()
    inner = joined.inner

    precision = max(int(inner.precision), 260)
    digits = min(precision - 12, 120)
    with mp.workdps(precision):
        nu = mp.mpf(".01")
        Z = mp.mpf(".3")
        delta = mp.mpf(inner.core.delta)
        Lambda = mp.mpf(inner.core.Lambda)
        R = 1 / Lambda
        phi = mp.mpf(".37")

        # This point is inside the regular polynomial core: the core radius is
        # 4.1 / Lambda, so R = 1 / Lambda has a comfortable margin.
        sample = dict(inner.core.evaluate(R, Z))
        sample["R"] = R
        physical = inner.physical_field(nu=nu, T=0)

        raw_rows: list[dict[str, object]] = []
        for alpha_int in (0, 2, 4, 6):
            alpha = mp.mpf(alpha_int)
            logq = -alpha / delta
            # Supplying delta explicitly is essential.  physical_chart's
            # default is a diagnostic value (1e-32), while this candidate is
            # built with delta = 1e-200.
            chart = physical_chart(
                sample,
                Z,
                logq,
                delta=delta,
                nu=nu,
                phi=phi,
                precision=precision,
            )
            # This is the actual callable roundtrip.  It recovers the internal
            # point from the physical coordinates and time returned above.
            back = physical.evaluate(*chart["xyz"], -chart["tau"])
            q_back = back["tau"] / (1 - back["Z"] * back["Z"])

            x_phys, y_phys, z_phys = back["xyz"]
            rho = mp.sqrt(x_phys * x_phys + y_phys * y_phys)
            z_abs = abs(z_phys)
            ur, utheta, uz = back["cylindrical_velocity"]
            aspect = z_abs / rho
            ut_over_uz = utheta / uz if uz != 0 else mp.mpf(0)
            winding = (
                utheta / (2 * mp.pi * rho * uz)
                if uz != 0 and rho != 0
                else mp.mpf(0)
            )
            # The axial cylindrical-curl component can be obtained directly
            # from the actual core jets returned by the callable.  With
            # u_theta = sqrt(nu) q^(-(1+delta)/2) sqrt(2R) F and
            # r = sqrt(2 nu q R),
            # omega_z = r^(-1) d_r(r u_theta)
            #          = q^(-(2+delta)/2) (2 F + 2 R F_R).
            # This is one physical vorticity component; it is not a claim
            # about the full curl or a vorticity concentration estimate.
            q = q_back
            # LocalCoreExitField returns chart values and the recovered R, Z;
            # it intentionally does not copy the internal jets into its public
            # result. Re-evaluate the same actual recovered core point to use
            # those jets, rather than a formula-only value at the target.
            roundtrip_jets = inner.core.evaluate(back["R"], back["Z"])
            omega_z = q ** (-(2 + delta) / 2) * (
                2 * roundtrip_jets["F"] + 2 * back["R"] * roundtrip_jets["F_R"]
            )

            velocity_errors = [
                _error(actual - expected, digits)
                for actual, expected in zip(
                    back["cylindrical_velocity"], chart["cylindrical_velocity"]
                )
            ]
            velocity_relative_errors = [
                _error(
                    (actual - expected) / expected if expected != 0 else actual,
                    digits,
                )
                for actual, expected in zip(
                    back["cylindrical_velocity"], chart["cylindrical_velocity"]
                )
            ]
            raw_rows.append(
                {
                    "alpha": alpha_int,
                    "logq": _text(logq, digits),
                    "q": signed_log(mp.exp(logq), digits),
                    "q_roundtrip": signed_log(q_back, digits),
                    "tau": signed_log(chart["tau"], digits),
                    "region": back.get("region"),
                    "R_roundtrip": _text(back["R"], digits),
                    "Z_roundtrip": _text(back["Z"], digits),
                    "R_relative_error": _error((back["R"] - R) / R, digits),
                    "Z_absolute_error": _error(back["Z"] - Z, digits),
                    "velocity_roundtrip_absolute_error": velocity_errors,
                    "velocity_roundtrip_relative_error": velocity_relative_errors,
                    "rho": signed_log(rho, digits),
                    "axial_coordinate_abs": signed_log(z_abs, digits),
                    "aspect": signed_log(aspect, digits),
                    "cylindrical_velocity": {
                        "ur": signed_log(ur, digits),
                        "utheta": signed_log(utheta, digits),
                        "uz": signed_log(uz, digits),
                    },
                    "utheta_over_uz": signed_log(ut_over_uz, digits),
                    "local_winding_density": signed_log(winding, digits),
                    "omega_z": signed_log(omega_z, digits),
                    # Private MP values are removed before JSON serialization;
                    # they let us form ratio logs without binary float loss.
                    "_rho": rho,
                    "_z_abs": z_abs,
                    "_aspect": aspect,
                    "_ut_over_uz": ut_over_uz,
                    "_winding": winding,
                    "_omega_z": omega_z,
                    "_utheta": utheta,
                    "_uz": uz,
                }
            )

        baseline = raw_rows[0]
        base_rho = baseline["_rho"]
        base_z_abs = baseline["_z_abs"]
        base_aspect = baseline["_aspect"]
        base_ut_over_uz = baseline["_ut_over_uz"]
        base_winding = baseline["_winding"]
        base_omega_z = baseline["_omega_z"]
        base_utheta = baseline["_utheta"]
        base_uz = baseline["_uz"]

        rows: list[dict[str, object]] = []
        for raw in raw_rows:
            alpha = mp.mpf(raw["alpha"])
            # The ratios use actual values from the callable roundtrip, not
            # the chart's expected power laws.
            ratio_rho = raw["_rho"] / base_rho
            ratio_z = raw["_z_abs"] / base_z_abs
            ratio_aspect = raw["_aspect"] / base_aspect
            ratio_ut_over_uz = (
                raw["_ut_over_uz"] / base_ut_over_uz
                if base_ut_over_uz != 0
                else mp.mpf(0)
            )
            ratio_winding = (
                raw["_winding"] / base_winding if base_winding != 0 else mp.mpf(0)
            )
            ratio_omega_z = (
                raw["_omega_z"] / base_omega_z if base_omega_z != 0 else mp.mpf(0)
            )
            ratio_utheta = (
                raw["_utheta"] / base_utheta if base_utheta != 0 else mp.mpf(0)
            )
            ratio_uz = raw["_uz"] / base_uz if base_uz != 0 else mp.mpf(0)
            expected_log_aspect_gain = alpha / 2
            expected_log_winding_gain = alpha / (2 * delta)
            expected_log_omega_z_gain = alpha * (2 + delta) / (2 * delta)
            expected_log_velocity_gain = alpha * (1 + delta) / (2 * delta)
            row = {key: value for key, value in raw.items() if not key.startswith("_")}
            row.update(
                {
                    "radial_ratio_to_alpha0": signed_log(ratio_rho, digits),
                    "axial_ratio_to_alpha0": signed_log(ratio_z, digits),
                    "aspect_ratio_to_alpha0": signed_log(ratio_aspect, digits),
                    "utheta_over_uz_ratio_to_alpha0": signed_log(
                        ratio_ut_over_uz, digits
                    ),
                    "winding_ratio_to_alpha0": signed_log(ratio_winding, digits),
                    "omega_z_ratio_to_alpha0": signed_log(ratio_omega_z, digits),
                    "utheta_ratio_to_alpha0": signed_log(ratio_utheta, digits),
                    "uz_ratio_to_alpha0": signed_log(ratio_uz, digits),
                    "expected_log_aspect_gain": _text(expected_log_aspect_gain, digits),
                    "expected_log_winding_gain": _text(
                        expected_log_winding_gain, digits
                    ),
                    "expected_log_omega_z_gain": _text(
                        expected_log_omega_z_gain, digits
                    ),
                    "expected_log_velocity_gain": _text(
                        expected_log_velocity_gain, digits
                    ),
                    "measured_minus_expected_log_aspect_gain": _text(
                        mp.log(abs(ratio_aspect)) - expected_log_aspect_gain,
                        digits,
                    ),
                    "measured_minus_expected_log_winding_gain": _text(
                        mp.log(abs(ratio_winding)) - expected_log_winding_gain,
                        digits,
                    ),
                    "measured_minus_expected_log_omega_z_gain": _text(
                        mp.log(abs(ratio_omega_z)) - expected_log_omega_z_gain,
                        digits,
                    ),
                    "measured_minus_expected_log_utheta_gain": _text(
                        mp.log(abs(ratio_utheta)) - expected_log_velocity_gain,
                        digits,
                    ),
                    "measured_minus_expected_log_uz_gain": _text(
                        mp.log(abs(ratio_uz)) - expected_log_velocity_gain,
                        digits,
                    ),
                }
            )
            rows.append(row)

        receipt: dict[str, object] = {
            "kind": "bounded_local_core_scale_geometry",
            "status": "measured",
            "field": "build_joined_field().inner.physical_field()",
            "precision_digits": precision,
            "sample": {
                "Z": _text(Z, digits),
                "R": _text(R, digits),
                "R_over_core_radius": _text(R / (mp.mpf("4.1") / Lambda), digits),
                "Lambda": _text(Lambda, digits),
                "nu": _text(nu, digits),
                "delta": _text(delta, digits),
                "phi": _text(phi, digits),
            },
            "time_parameterization": {
                "definition": "logq = -alpha / delta; tau = q * (1 - Z^2)",
                "alphas": [0, 2, 4, 6],
                "interpretation": (
                    "These are arbitrary MP chart scales, not practical simulation times."
                ),
            },
            "rows": rows,
            "vorticity": {
                "available": True,
                "component": "omega_z",
                "formula": "q^(-(2+delta)/2) * (2*F + 2*R*F_R)",
                "full_vector_available": False,
                "reason": (
                    "The axial component is measured from the actual F and F_R jets. "
                    "The full cylindrical curl remains unavailable because the core "
                    "state does not expose Ur_R and Ur_Z physical-space jets; an "
                    "extreme-scale finite-difference or inferred-jet construction is "
                    "not used."
                ),
            },
            "scope": {
                "measured": [
                    "actual physical chart coordinates and cylindrical velocities",
                    "actual physical.evaluate roundtrip at each alpha",
                    "radial, axial, aspect, velocity-ratio, and local winding ratios",
                    "the axial vorticity component omega_z from returned F and F_R jets",
                ],
                "not_certified": [
                    "recursive transition across scales",
                    "finite energy",
                    "global flow or outer matching",
                    "continuous-in-time closure from these four samples",
                    "full vorticity vector (missing Ur_R and Ur_Z physical jets)",
                    "integrated streamline winding from the pointwise winding density",
                ],
            },
        }

    OUT_JSON.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_JSON}", flush=True)
    print("measured rows: 4; omega_z component: available; full vorticity vector: unavailable", flush=True)
    return receipt


if __name__ == "__main__":
    run()
