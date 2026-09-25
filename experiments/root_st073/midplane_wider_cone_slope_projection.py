"""Project full frozen-wave momentum onto a supported time-slope basis.

The projection includes exact-curl harmonic/mean velocity derivatives and
pressure gradients. It tests instantaneous representability, not whether
the derivative can be integrated over a pulse with endpoint control.
"""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, polynomial_data,
                                       sample_residual)
from curl_wave_prototype import (LocalizedCurlWave, WavePerturbedField,
                                 cylindrical_residual)
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


class LinearSlopeField:
    """Exact-curl velocity with a local linear-in-tau potential and pressure."""

    def __init__(self, base, wave, harmonic, mean):
        self.base, self.wave = base, wave
        self.harmonic, self.mean = harmonic, mean
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        times = np.broadcast_to(tau, (len(pts),))
        if not np.all(times == times[0]):
            raise ValueError("Evaluate one tau per spatial batch")
        dt = float(times[0] - self.wave.tau0)
        field = FrozenPotentialField(
            self.base, self.wave, .1,
            [dt * fit[:27] for fit in self.harmonic],
            dt * self.mean[:18])
        velocity, pressure = field.fields(pts, times)
        for i, (x, y, z) in enumerate(pts):
            r = np.hypot(x, y)
            if (r == 0 or abs(r - self.wave.radius) >= self.wave.radial_halfwidth
                    or abs(z - self.wave.zcenter) >= self.wave.axial_halfwidth):
                continue
            theta = np.arctan2(y, x)
            for mode, fit in zip(self.wave.waves, self.harmonic):
                kr, _, kz = mode["normal"]
                phase = (mode["m"] * theta + kr * (r - self.wave.radius)
                         + kz * (z - self.wave.zcenter)
                         + mode["omega"] * dt)
                q = polynomial_data(self.wave, r, z, fit[27:], 1)[0][0]
                pressure[i] += np.real(np.exp(1j * phase) * q)
            pressure[i] += polynomial_data(
                self.wave, r, z, self.mean[18:], 1)[0][0]
        return velocity, pressure


def row_stats(rows):
    values = np.concatenate([np.linalg.norm(row["residual"], axis=1)
                             for row in rows])
    return dict(max=float(np.max(values)),
                rms=float(np.sqrt(np.mean(values**2))))


def run():
    inner, fields = build_fields()
    repair = json.loads((ROOT / "midplane_axial_cone_all_knots_repair.json").read_text())
    mean = SeparatedMomentModes(
        fields["two_sided_cone"], repair["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    angles = np.arange(16) * 2 * np.pi / 16
    train_axis = (-.6, -.3, 0., .3, .6)
    held_axis = (-.45, -.15, .15, .45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    held_grid = [(x, y) for x in held_axis for y in held_axis]
    scales = []
    for k in (11, 19):
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{k}.json").read_text())
        source["nu"] = mean.nu
        tau = float(source["tau"])
        point = np.asarray(source["point"], float)
        source["velocity"] = mean.fields(point[None, :], tau)[0][0].tolist()
        target = stress_primitive(mean, float(point[0]), float(point[2]), tau)
        source["local_tangential_stress_primitive"] = target.tolist()
        selected = next(row["selected"]["indices"] for row in pairs["scales"]
                        if row["k"] == k)
        args = support(inner, source, mean.nu)
        args["axial_halfwidth"] *= 1.4
        args["time_halfwidth"] *= 1.4**2
        wave = LocalizedCurlWave(source, pulse_indices=selected, **args)
        center_points = np.array([
            (point[0] * np.cos(a), point[0] * np.sin(a), point[2])
            for a in angles])
        columns = []
        for j in range(2):
            wave.weights = np.eye(2)[j]
            values, _ = wave.fields(center_points, tau)
            cyl = cylindrical_residual(values, center_points)
            columns.append(np.mean(cyl[:, 0, None] * cyl[:, 1:], axis=0))
        weights = np.linalg.solve(np.column_stack(columns), target)
        if np.any(weights <= 0):
            raise ValueError(f"k={k}: covariance weights are not positive")
        wave.weights = weights
        frozen = WavePerturbedField(mean, ScaledWave(wave, .1))
        train = sample_residual(frozen, wave, tau, train_grid, angles,
                                time_step=args["time_halfwidth"],
                                time_min=.5 * 2.**-20)
        held = sample_residual(frozen, wave, tau, held_grid, angles,
                               time_step=args["time_halfwidth"],
                               time_min=.5 * 2.**-20)
        harmonic, mean_fit = fit_slope(train, wave, tau, angles,
                                       regularization=1e-3)
        train_stats = metrics(train, wave, tau, angles, harmonic, mean_fit)
        held_stats = metrics(held, wave, tau, angles, harmonic, mean_fit)
        realized = LinearSlopeField(mean, wave, harmonic, mean_fit)
        direct = []
        for offset in (0., .1):
            check_tau = tau + offset * args["time_halfwidth"]
            direct_rows = sample_residual(
                realized, wave, check_tau, held_grid, angles,
                time_step=args["time_halfwidth"], time_min=.5 * 2.**-20)
            frozen_rows = (held if offset == 0 else sample_residual(
                frozen, wave, check_tau, held_grid, angles,
                time_step=args["time_halfwidth"], time_min=.5 * 2.**-20))
            direct.append(dict(time_offset_fraction=offset,
                               frozen=row_stats(frozen_rows),
                               corrected=row_stats(direct_rows)))
        scales.append(dict(k=k, tau=tau, weights=weights.tolist(),
                           train=train_stats, heldout=held_stats,
                           direct_replay=direct,
                           harmonic_potential_slope_norms=[
                               float(np.linalg.norm(f[:27])) for f in harmonic],
                           harmonic_pressure_norms=[
                               float(np.linalg.norm(f[27:])) for f in harmonic],
                           mean_potential_slope_norm=float(
                               np.linalg.norm(mean_fit[:18])),
                           mean_pressure_norm=float(
                               np.linalg.norm(mean_fit[18:]))))
        print(f"k={k}: heldout max {held_stats['frozen']['max']:.6g} "
              f"-> {held_stats['projected_derivative_and_pressure']['max']:.6g}; "
              f"direct {[row['corrected']['max'] for row in direct]}",
              flush=True)
    report = dict(source="Wider-cone instantaneous full-momentum slope projection",
                  train_axis=train_axis, held_axis=held_axis,
                  angles_per_node=len(angles), scales=scales,
                  scope="Exact-curl linear-in-tau potential and pressure "
                        "directly evaluated at pulse center and one interior "
                        "time. No solved amplitude trajectory, endpoint "
                        "cutoff, uniform space-time residual, or recursive "
                        "acceptance.",
                  accepted=False, scale_recursion_established=False)
    (ROOT / "midplane_wider_cone_slope_projection.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
