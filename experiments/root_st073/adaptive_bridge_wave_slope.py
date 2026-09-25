"""Fit a first local time-amplitude slope for the two exact-curl modes."""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from adaptive_bridge_recursive_defect import build_fields
from curl_wave_prototype import LocalizedCurlWave, WavePerturbedField
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    source = json.loads((ROOT / "compact_potential" /
                         "adaptive_bridge_wave_source.json").read_text())
    source["nu"] = base.nu
    tau = float(source["tau"])
    inner_radius = float(inner.from_similarity([inner.p.X_max], [-0.2], tau)[0, 0])
    halfwidth = 0.0047 * 15.0 * inner_radius
    time_halfwidth = 0.04 * tau
    common = dict(radial_halfwidth=halfwidth, axial_halfwidth=1.0e-4,
                  time_halfwidth=time_halfwidth, min_tau=0.5 * 2.0**-20)
    first = LocalizedCurlWave(source, **common)

    def field(rates):
        wave = LocalizedCurlWave(
            source, pulse_indices=first.pulse_indices,
            time_slope_rates=rates, **common)
        return WavePerturbedField(base, ScaledWave(wave, 0.1))

    def points(ys, angle_shift):
        angles = (np.arange(16) + angle_shift) * 2.0 * np.pi / 16.0
        blocks = []
        for y in ys:
            X = inner.p.X_max * (1.0 + 15.0 * y)**2
            location = inner.from_similarity([X], [-0.2], tau)[0]
            radius, z = float(location[0]), float(location[2])
            blocks.extend((radius * np.cos(angle), radius * np.sin(angle), z)
                          for angle in angles)
        return np.asarray(blocks)

    train = points((0.0475, 0.05, 0.0525), 0.0)
    holdout = points((0.04875, 0.05125), 0.5)
    h, ht = 0.0005 * np.sqrt(base.nu * tau), 0.0001 * tau

    def residual(current, pts):
        velocity, gradient, linear = kinematics(
            current, pts, tau, h, ht, time_min=0.5 * 2.0**-20)
        return linear + np.einsum("nij,nj->ni", gradient, velocity)

    zero_rates = np.zeros(2)
    before = residual(field(zero_rates), train)
    unit_rate = 1.0e6
    columns = []
    for j in range(2):
        rates = np.zeros(2)
        rates[j] = unit_rate
        columns.append(((residual(field(rates), train) - before)
                        / unit_rate).ravel())
    matrix = np.column_stack(columns)
    rates = np.linalg.lstsq(matrix, -before.ravel(), rcond=None)[0]
    after = residual(field(rates), train)
    old_holdout = residual(field(zero_rates), holdout)
    new_holdout = residual(field(rates), holdout)
    report = dict(
        source="First local time-slope fit for two compact exact-curl bridge modes",
        k=11.0, amplitude_scale=0.1,
        pulse_indices=first.pulse_indices,
        fitted_time_slope_rates=rates.tolist(),
        rate_times_halfwidth=(rates * time_halfwidth).tolist(),
        slope_matrix_singular_values=np.linalg.svd(matrix, compute_uv=False).tolist(),
        training=dict(before_max=float(np.max(np.linalg.norm(before, axis=1))),
                      after_max=float(np.max(np.linalg.norm(after, axis=1))),
                      before_rms=float(np.sqrt(np.mean(np.sum(before**2, axis=1)))),
                      after_rms=float(np.sqrt(np.mean(np.sum(after**2, axis=1))))),
        holdout=dict(before_max=float(np.max(np.linalg.norm(old_holdout, axis=1))),
                     after_max=float(np.max(np.linalg.norm(new_holdout, axis=1))),
                     before_rms=float(np.sqrt(np.mean(np.sum(old_holdout**2, axis=1)))),
                     after_rms=float(np.sqrt(np.mean(np.sum(new_holdout**2, axis=1))))),
        scope="Local first Taylor slope of each wave amplitude at one time, fitted on three radii and 16 angles, tested on two radial/phase holdouts. Does not solve the transported amplitude PDE or preserve smallness throughout the pulse window; no global momentum or scale-recursion acceptance.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_bridge_wave_slope.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rates=report["fitted_time_slope_rates"],
                          rate_times_halfwidth=report["rate_times_halfwidth"],
                          training=report["training"], holdout=report["holdout"]),
          indent=2))
    return report


if __name__ == "__main__":
    run()
