"""Two spatial wave/mean slope stages and actual interval midpoint checks."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, sample_residual)
from curl_wave_patch_trajectory import FirstIntervalField, HermiteIntervalField
from curl_wave_prototype import LocalizedCurlWave
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def complex_rows(rows):
    return [np.array([complex(re, im) for re, im in row]) for row in rows]


def encode(harmonic, mean, tau):
    return dict(
        tau=tau,
        harmonic_slope_coefficients=[[[float(c.real), float(c.imag)]
                                      for c in coefficients]
                                     for coefficients in harmonic],
        mean_slope_coefficients=mean.tolist())


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    source = json.loads((ROOT / "compact_potential" /
                         "adaptive_bridge_wave_source.json").read_text())
    source["nu"] = base.nu
    tau0 = float(source["tau"])
    ri = float(inner.from_similarity([inner.p.X_max], [-0.2], tau0)[0, 0])
    wave = LocalizedCurlWave(
        source, radial_halfwidth=0.0047 * 15.0 * ri,
        axial_halfwidth=1.0e-4, time_halfwidth=0.04 * tau0,
        min_tau=0.5 * 2.0**-20)
    previous = json.loads((ROOT / "adaptive_bridge_wave_spatial_slope.json").read_text())
    dt = float(previous["time_step"])
    amplitude = float(previous["wave_amplitude"])
    stage0 = previous["stage"]
    first_slopes = complex_rows(stage0["harmonic_slope_coefficients"])
    first_mean = np.asarray(stage0["mean_slope_coefficients"])
    initial_harmonic = [dt * row[:27] for row in first_slopes]
    initial_mean = dt * first_mean[:18]
    tau1 = tau0 + dt
    frozen1 = FrozenPotentialField(base, wave, amplitude,
                                   initial_harmonic, initial_mean)
    angles = np.arange(16) * 2.0 * np.pi / 16.0
    train_axis = (-0.6, -0.2, 0.2, 0.6)
    holdout_axis = (-0.45, 0.0, 0.45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    holdout_grid = [(x, y) for x in holdout_axis for y in holdout_axis]
    time_min = 0.5 * 2.0**-20
    train_rows = sample_residual(
        frozen1, wave, tau1, train_grid, angles,
        time_step=dt, time_min=time_min)
    holdout_rows = sample_residual(
        frozen1, wave, tau1, holdout_grid, angles,
        time_step=dt, time_min=time_min)
    harmonic1, mean1 = fit_slope(train_rows, wave, tau1, angles)
    stage1 = encode(harmonic1, mean1, tau1)
    projected_train1 = metrics(train_rows, wave, tau1, angles,
                               harmonic1, mean1)
    projected_holdout1 = metrics(holdout_rows, wave, tau1, angles,
                                 harmonic1, mean1)
    interval0_linear = FirstIntervalField(base, wave, amplitude, stage0)
    interval0_hermite = HermiteIntervalField(base, wave, amplitude,
                                             stage0, stage1, dt)
    interval1_linear = FirstIntervalField(base, wave, amplitude, stage1,
                                           initial_harmonic, initial_mean)
    frozen0 = FrozenPotentialField(
        base, wave, amplitude,
        [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
    shifted_angles = (np.arange(16) + 0.5) * 2.0 * np.pi / 16.0
    points = []
    for xi, eta in holdout_grid:
        r = wave.radius + wave.radial_halfwidth * xi
        z = wave.zcenter + wave.axial_halfwidth * eta
        points.extend((r * np.cos(angle), r * np.sin(angle), z)
                      for angle in shifted_angles)
    points = np.asarray(points)

    def direct(field, tau):
        hs = 0.0005 * np.sqrt(base.nu * tau)
        velocity, gradient, linear = kinematics(
            field, points, tau, hs, dt / 8.0, time_min=time_min)
        residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
        norms = np.linalg.norm(residual, axis=1)
        return dict(maximum=float(np.max(norms)),
                    rms=float(np.sqrt(np.mean(norms**2))),
                    fd_divergence_max=float(np.max(np.abs(np.trace(
                        gradient, axis1=1, axis2=2)))))

    first_mid = tau0 + 0.5 * dt
    second_mid = tau1 + 0.5 * dt
    midpoint_results = dict(
        first_frozen=direct(frozen0, first_mid),
        first_linear=direct(interval0_linear, first_mid),
        first_hermite=direct(interval0_hermite, first_mid),
        second_frozen_updated=direct(frozen1, second_mid),
        second_linear=direct(interval1_linear, second_mid),
    )
    report = dict(
        source="Two spatial wave/mean slope stages on moment-closed late bridge",
        tau0=tau0, time_step=dt,
        projected_stage0_train=previous["projected_train"],
        projected_stage0_holdout=previous["projected_holdout"],
        projected_stage1_train=projected_train1,
        projected_stage1_holdout=projected_holdout1,
        midpoint_direct=midpoint_results,
        stage1=stage1,
        scope="Two local coefficient-derivative stages and direct full Cartesian momentum at two interval midpoints, with shifted angular and off-grid radial/axial nodes. No pulse endpoint match, global support/energy proof, dyadic recursion, or full-domain 1e-3 gate.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_bridge_wave_two_stage.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          projected_stage1_train=projected_train1,
                          projected_stage1_holdout=projected_holdout1,
                          midpoint_direct=midpoint_results), indent=2))
    return report


if __name__ == "__main__":
    run()
