"""Spatial wave/pressure/mean slope on the three-knot moment-closed bridge."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, sample_residual)
from curl_wave_patch_trajectory import FirstIntervalField
from curl_wave_prototype import LocalizedCurlWave
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run(dt=1.0e-10):
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
    amplitude = 0.1
    frozen = FrozenPotentialField(
        base, wave, amplitude,
        [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
    angles = np.arange(16) * 2.0 * np.pi / 16.0
    train_axis = (-0.6, -0.2, 0.2, 0.6)
    holdout_axis = (-0.45, 0.0, 0.45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    holdout_grid = [(x, y) for x in holdout_axis for y in holdout_axis]
    time_min = 0.5 * 2.0**-20
    training_rows = sample_residual(
        frozen, wave, tau0, train_grid, angles,
        time_step=dt, time_min=time_min)
    heldout_rows = sample_residual(
        frozen, wave, tau0, holdout_grid, angles,
        time_step=dt, time_min=time_min)
    harmonic_fit, mean_fit = fit_slope(
        training_rows, wave, tau0, angles)
    projected_train = metrics(training_rows, wave, tau0, angles,
                              harmonic_fit, mean_fit)
    projected_holdout = metrics(heldout_rows, wave, tau0, angles,
                                harmonic_fit, mean_fit)
    stage = dict(
        tau=tau0,
        harmonic_slope_coefficients=[[[float(c.real), float(c.imag)]
                                      for c in coefficients]
                                     for coefficients in harmonic_fit],
        mean_slope_coefficients=mean_fit.tolist(),
    )
    evolved = FirstIntervalField(base, wave, amplitude, stage)
    midpoint_tau = tau0 + dt / 2.0
    shifted_angles = (np.arange(16) + 0.5) * 2.0 * np.pi / 16.0
    points = []
    for xi, eta in holdout_grid:
        r = wave.radius + wave.radial_halfwidth * xi
        z = wave.zcenter + wave.axial_halfwidth * eta
        points.extend((r * np.cos(angle), r * np.sin(angle), z)
                      for angle in shifted_angles)
    points = np.asarray(points)
    hs, ht = 0.0005 * np.sqrt(base.nu * midpoint_tau), dt / 8.0

    def direct(field):
        velocity, gradient, linear = kinematics(
            field, points, midpoint_tau, hs, ht, time_min=time_min)
        residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
        norms = np.linalg.norm(residual, axis=1)
        return dict(maximum=float(np.max(norms)),
                    rms=float(np.sqrt(np.mean(norms**2))),
                    fd_divergence_max=float(np.max(np.abs(np.trace(
                        gradient, axis1=1, axis2=2)))),
                    speed_max=float(np.max(np.linalg.norm(velocity, axis=1))))

    direct_frozen = direct(frozen)
    direct_evolved = direct(evolved)
    report = dict(
        source="One spatially varying amplitude/pressure/mean Taylor step on late moment-closed bridge",
        tau0=tau0, time_step=dt,
        wave_amplitude=amplitude,
        train_nodes=len(train_grid), holdout_nodes=len(holdout_grid),
        angles_per_node=len(angles),
        projected_train=projected_train,
        projected_holdout=projected_holdout,
        harmonic_slope_norms=[float(np.linalg.norm(x[:27]))
                              for x in harmonic_fit],
        harmonic_pressure_norms=[float(np.linalg.norm(x[27:]))
                                 for x in harmonic_fit],
        mean_slope_norm=float(np.linalg.norm(mean_fit[:18])),
        mean_pressure_norm=float(np.linalg.norm(mean_fit[18:])),
        midpoint_direct_frozen=direct_frozen,
        midpoint_direct_evolved=direct_evolved,
        stage=stage,
        scope="One local projected slope using 16 radial-axial train nodes and 9 off-grid nodes, then direct full Cartesian momentum at a single midpoint with shifted angular nodes. No completed transported PDE, pulse endpoint match, global stress cone, scale recursion, or full-domain 1e-3 acceptance.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_bridge_wave_spatial_slope.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          projected_train=projected_train,
                          projected_holdout=projected_holdout,
                          midpoint_direct_frozen=direct_frozen,
                          midpoint_direct_evolved=direct_evolved,
                          slope_norms=report["harmonic_slope_norms"],
                          mean_slope_norm=report["mean_slope_norm"]),
          indent=2))
    return report


if __name__ == "__main__":
    run()
