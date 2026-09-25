"""Try a second pulse-time slope and direct quadratic potential replay."""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, polynomial_data,
                                       sample_residual)
from curl_wave_prototype import WavePerturbedField
from midplane_wider_cone_slope_projection import (LinearSlopeField,
                                                   build_case, row_stats)
from radial_continuation import ROOT


class QuadraticSlopeField:
    """Exact-curl quadratic potential, with linearly updated pressure."""

    def __init__(self, base, wave, initial_h, initial_mean,
                 delta_h, delta_mean, horizon):
        self.base, self.wave = base, wave
        self.initial_h, self.initial_mean = initial_h, initial_mean
        self.delta_h, self.delta_mean = delta_h, delta_mean
        self.horizon = horizon
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        times = np.broadcast_to(tau, (len(pts),))
        if not np.all(times == times[0]):
            raise ValueError("Evaluate one tau per spatial batch")
        dt = float(times[0] - self.wave.tau0)
        quadratic = .5 * dt**2 / self.horizon
        fraction = dt / self.horizon
        field = FrozenPotentialField(
            self.base, self.wave, .1,
            [dt * initial[:27] + quadratic * delta[:27]
             for initial, delta in zip(self.initial_h, self.delta_h)],
            dt * self.initial_mean[:18] + quadratic * self.delta_mean[:18])
        velocity, pressure = field.fields(pts, times)
        for i, (x, y, z) in enumerate(pts):
            r = np.hypot(x, y)
            if (r == 0 or abs(r - self.wave.radius) >= self.wave.radial_halfwidth
                    or abs(z - self.wave.zcenter) >= self.wave.axial_halfwidth):
                continue
            theta = np.arctan2(y, x)
            for mode, initial, delta in zip(
                    self.wave.waves, self.initial_h, self.delta_h):
                kr, _, kz = mode["normal"]
                phase = (mode["m"] * theta + kr * (r - self.wave.radius)
                         + kz * (z - self.wave.zcenter)
                         + mode["omega"] * dt)
                coefficient = initial[27:] + fraction * delta[27:]
                q = polynomial_data(self.wave, r, z, coefficient, 1)[0][0]
                pressure[i] += np.real(np.exp(1j * phase) * q)
            pressure[i] += polynomial_data(
                self.wave, r, z,
                self.initial_mean[18:] + fraction * self.delta_mean[18:],
                1)[0][0]
        return velocity, pressure


def run():
    angles = np.arange(16) * 2 * np.pi / 16
    train_axis = (-.6, -.3, 0., .3, .6)
    held_axis = (-.45, -.15, .15, .45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    held_grid = [(x, y) for x in held_axis for y in held_axis]
    prior = json.loads((ROOT / "midplane_wider_cone_slope_projection.json").read_text())
    scales = []
    for k in (11, 19):
        mean, wave, args, _ = build_case(k, angles)
        tau = wave.tau0
        horizon = .1 * args["time_halfwidth"]
        frozen = WavePerturbedField(mean, ScaledWave(wave, .1))
        rows0 = sample_residual(frozen, wave, tau, train_grid, angles,
                                time_step=args["time_halfwidth"],
                                time_min=.5 * 2.**-20)
        first_h, first_mean = fit_slope(rows0, wave, tau, angles,
                                        regularization=1e-3)
        linear = LinearSlopeField(mean, wave, first_h, first_mean)
        tau1 = tau + horizon
        rows1 = sample_residual(linear, wave, tau1, train_grid, angles,
                                time_step=args["time_halfwidth"],
                                time_min=.5 * 2.**-20)
        held1 = sample_residual(linear, wave, tau1, held_grid, angles,
                                time_step=args["time_halfwidth"],
                                time_min=.5 * 2.**-20)
        delta_h, delta_mean = fit_slope(rows1, wave, tau1, angles,
                                        regularization=1e-3)
        projected = metrics(held1, wave, tau1, angles, delta_h, delta_mean)
        quadratic = QuadraticSlopeField(mean, wave, first_h, first_mean,
                                        delta_h, delta_mean, horizon)
        direct = []
        for fraction in (.05, .1):
            rows = sample_residual(
                quadratic, wave, tau + fraction * args["time_halfwidth"],
                held_grid, angles, time_step=args["time_halfwidth"],
                time_min=.5 * 2.**-20)
            direct.append(dict(time_offset_fraction=fraction,
                               corrected=row_stats(rows)))
        baseline = next(s for s in prior["scales"] if s["k"] == k)
        prior_linear = baseline["direct_replay"][1]
        scales.append(dict(k=k, horizon=horizon,
                           linear_at_horizon=row_stats(held1),
                           second_slope_projection=projected,
                           quadratic_direct=direct,
                           prior_linear_at_horizon=prior_linear,
                           relative_delta_h_norms=[
                               float(np.linalg.norm(d[:27])
                                     / max(np.linalg.norm(f[:27]), 1e-30))
                               for d, f in zip(delta_h, first_h)],
                           relative_delta_mean_norm=float(
                               np.linalg.norm(delta_mean[:18])
                               / max(np.linalg.norm(first_mean[:18]), 1e-30))))
        print(f"k={k}: projected at +0.1 = "
              f"{projected['projected_derivative_and_pressure']['max']:.6g}; "
              f"quadratic direct = {[x['corrected']['max'] for x in direct]}",
              flush=True)
    report = dict(source="Second pulse-time slope on wider-cone exact-curl field",
                  train_axis=train_axis, held_axis=held_axis,
                  angles_per_node=len(angles), scales=scales,
                  scope="Two instantaneous full-momentum derivative/pressure "
                        "fits followed by quadratic-in-tau exact-curl direct "
                        "replay on held-out nodes. No nonlinear interval "
                        "solve, support endpoint control, or scale recursion.",
                  accepted=False, scale_recursion_established=False)
    (ROOT / "midplane_wider_cone_second_slope.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
