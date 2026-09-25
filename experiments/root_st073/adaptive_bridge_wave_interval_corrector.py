"""Damped trapezoid endpoint solve for the late-bridge curl-wave coefficients.

This attempts one interval beyond the explicit Taylor-step horizon. Its
fixed-point defect and direct interior residual are both reported; a small
endpoint defect alone is not a Navier--Stokes acceptance condition.
"""

import argparse
import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_wave_two_stage import encode
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, sample_residual)
from curl_wave_patch_trajectory import HermiteIntervalField
from curl_wave_prototype import LocalizedCurlWave
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def norm(harmonic, mean):
    return float(np.sqrt(sum(np.vdot(row, row).real for row in harmonic)
                         + np.dot(mean, mean)))


def run(dt=1e-8, damping=0.25, iterations=4):
    if dt <= 0 or not 0 < damping <= 1 or iterations < 1:
        raise ValueError("Need positive step, damping in (0,1], and iterations")
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
    if dt >= wave.time_halfwidth / 2:
        raise ValueError("Step must remain below half the wave time half-width")
    amplitude = 0.1
    angles = np.arange(16) * 2.0 * np.pi / 16.0
    train_axis = (-0.6, -0.2, 0.2, 0.6)
    holdout_axis = (-0.45, 0.0, 0.45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    holdout_grid = [(x, y) for x in holdout_axis for y in holdout_axis]
    time_min = 0.5 * 2.0**-20
    empty_harmonic = [np.zeros(27, complex) for _ in wave.waves]
    empty_mean = np.zeros(18)
    frozen0 = FrozenPotentialField(base, wave, amplitude,
                                   empty_harmonic, empty_mean)
    train0 = sample_residual(frozen0, wave, tau0, train_grid, angles,
                             time_step=dt, time_min=time_min)
    slope0_h, slope0_m = fit_slope(train0, wave, tau0, angles)
    stage0 = encode(slope0_h, slope0_m, tau0)
    xh = [dt * row[:27] for row in slope0_h]
    xm = dt * slope0_m[:18]
    tau1 = tau0 + dt
    history = []
    for iteration in range(iterations):
        frozen = FrozenPotentialField(base, wave, amplitude, xh, xm)
        train = sample_residual(frozen, wave, tau1, train_grid, angles,
                                time_step=dt, time_min=time_min)
        slope1_h, slope1_m = fit_slope(train, wave, tau1, angles)
        candidate_h = [0.5 * dt * (a[:27] + b[:27])
                       for a, b in zip(slope0_h, slope1_h)]
        candidate_m = 0.5 * dt * (slope0_m[:18] + slope1_m[:18])
        difference_h = [a - b for a, b in zip(candidate_h, xh)]
        difference_m = candidate_m - xm
        defect = norm(difference_h, difference_m)
        relative = defect / max(norm(xh, xm), 1e-30)
        history.append(dict(iteration=iteration + 1,
                            coefficient_norm=norm(xh, xm),
                            fixed_point_defect=defect,
                            relative_fixed_point_defect=relative,
                            projected_train=metrics(train, wave, tau1, angles,
                                                    slope1_h, slope1_m)))
        print(f"iteration {iteration + 1}: relative defect {relative:.6g}",
              flush=True)
        xh = [a + damping * d for a, d in zip(xh, difference_h)]
        xm = xm + damping * difference_m

    endpoint = FrozenPotentialField(base, wave, amplitude, xh, xm)
    endpoint_train = sample_residual(endpoint, wave, tau1, train_grid,
                                     angles, time_step=dt,
                                     time_min=time_min)
    endpoint_holdout = sample_residual(endpoint, wave, tau1, holdout_grid,
                                       angles, time_step=dt,
                                       time_min=time_min)
    slope1_h, slope1_m = fit_slope(endpoint_train, wave, tau1, angles)
    projected_holdout = metrics(endpoint_holdout, wave, tau1, angles,
                                slope1_h, slope1_m)
    candidate_h = [0.5 * dt * (a[:27] + b[:27])
                   for a, b in zip(slope0_h, slope1_h)]
    candidate_m = 0.5 * dt * (slope0_m[:18] + slope1_m[:18])
    final_defect = norm([a - b for a, b in zip(candidate_h, xh)],
                        candidate_m - xm) / max(norm(xh, xm), 1e-30)
    # Hermite endpoint state matches the iterated coefficients exactly.
    effective_h = [np.r_[2 * value / dt - old[:27], new[27:]]
                   for value, old, new in zip(xh, slope0_h, slope1_h)]
    effective_m = np.r_[2 * xm / dt - slope0_m[:18], slope1_m[18:]]
    stage1 = encode(effective_h, effective_m, tau1)
    trajectory = HermiteIntervalField(base, wave, amplitude,
                                       stage0, stage1, dt)
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
        values = np.linalg.norm(residual, axis=1)
        return dict(maximum=float(np.max(values)),
                    rms=float(np.sqrt(np.mean(values**2))),
                    fd_divergence_max=float(np.max(np.abs(np.trace(
                        gradient, axis1=1, axis2=2)))))

    direct_midpoint = direct(trajectory, tau0 + 0.5 * dt)
    direct_frozen_midpoint = direct(frozen0, tau0 + 0.5 * dt)
    report = dict(source="Damped trapezoid endpoint coefficient solve",
                  tau0=tau0, time_step=dt, damping=damping,
                  iterations=iterations, history=history,
                  final_relative_fixed_point_defect=float(final_defect),
                  endpoint_projected_holdout=projected_holdout,
                  direct_frozen_midpoint=direct_frozen_midpoint,
                  direct_corrected_midpoint=direct_midpoint,
                  scope="One supported local patch interval and sampled collocation. Nonzero fixed-point defect, no full pulse, dyadic recursion, continuum residual, or PDE acceptance.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "adaptive_bridge_wave_interval_corrector.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          final_relative_fixed_point_defect=final_defect,
                          direct_frozen_midpoint=direct_frozen_midpoint,
                          direct_corrected_midpoint=direct_midpoint,
                          endpoint_projected_holdout=projected_holdout),
                     indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--time-step", type=float, default=1e-8)
    parser.add_argument("--damping", type=float, default=0.25)
    parser.add_argument("--iterations", type=int, default=4)
    args = parser.parse_args()
    run(args.time_step, args.damping, args.iterations)
