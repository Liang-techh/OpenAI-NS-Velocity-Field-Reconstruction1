"""Direct interior-time screen of the transferred one-instant slope.

This deliberately tests the simplest constant-slope continuation before
attempting a pulse-coordinate inverse. It is not a complete correction.
"""

import json
import math

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, polynomial_data,
                                       sample_residual)
from curl_wave_prototype import LocalizedCurlWave
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def decode(row):
    harmonics = [np.array([complex(*pair) for pair in values])
                 for values in row["harmonic"]]
    return harmonics, np.asarray(row["mean"], float)


class LinearSlopeField:
    """Exact-curl potential linear in time, plus fitted harmonic/mean pressure."""

    def __init__(self, base, wave, amplitude, harmonic_slope, mean_slope,
                 temporal_profile="linear"):
        self.base, self.wave, self.amplitude = base, wave, amplitude
        self.harmonic_slope, self.mean_slope = harmonic_slope, mean_slope
        self.temporal_profile = temporal_profile
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        if not np.all(ts == ts[0]):
            raise ValueError("The trajectory screen expects one time per call")
        dt = float(ts[0] - self.wave.tau0)
        if self.temporal_profile == "linear":
            potential_time, pressure_time = dt, 1.0
        elif self.temporal_profile == "integrated_bump":
            f = dt / self.wave.time_halfwidth
            if abs(f) >= 1:
                raise ValueError("Integrated-bump diagnostic is pulse interior only")
            primitive = sum((-1)**k * math.comb(5, k) * f**(2*k+1)
                            / (2*k+1) for k in range(6))
            potential_time = self.wave.time_halfwidth * primitive
            pressure_time = (1 - f*f)**5
        else:
            raise ValueError(self.temporal_profile)
        field = FrozenPotentialField(
            self.base, self.wave, self.amplitude,
            [potential_time * row[:27] for row in self.harmonic_slope],
            potential_time * self.mean_slope[:18])
        velocity, pressure = field.fields(pts, ts)
        for i, ((x, y, z), t) in enumerate(zip(pts, ts)):
            r = np.hypot(x, y)
            if (r == 0 or abs(r - self.wave.radius) >= self.wave.radial_halfwidth
                    or abs(z - self.wave.zcenter) >= self.wave.axial_halfwidth):
                continue
            theta = np.arctan2(y, x)
            for mode, coeff in zip(self.wave.waves, self.harmonic_slope):
                q = polynomial_data(self.wave, r, z, coeff[27:], 1)[0][0]
                kr, _, kz = mode["normal"]
                phase = (mode["m"] * theta + kr * (r - self.wave.radius)
                         + kz * (z - self.wave.zcenter)
                         + mode["omega"] * (t - self.wave.tau0))
                pressure[i] += pressure_time * (np.exp(1j * phase) * q).real
            pressure[i] += pressure_time * polynomial_data(
                self.wave, r, z, self.mean_slope[18:], 1)[0][0]
        return velocity, pressure


def stats(rows):
    norm = np.concatenate([np.linalg.norm(row["residual"], axis=1)
                           for row in rows])
    return dict(max=float(norm.max()), rms=float(np.sqrt(np.mean(norm**2))))


def run(profile_only=False):
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    transfer = json.loads((ROOT / "midplane_wave_scale_transfer_slope.json").read_text())
    source = json.loads((ROOT / "compact_potential" /
                         "midplane_wave_source_k19.json").read_text())
    source["nu"] = base.nu
    selected = next(row["selected"] for row in pairs["scales"]
                    if row["k"] == 19)
    wave = LocalizedCurlWave(source, pulse_indices=selected["indices"],
                             **support(inner, source, base.nu))
    wave.weights = np.asarray(selected["positive_weights"], float)
    first_h, first_m = decode(transfer["scales"][0]["fitted_coefficients"])
    factor = transfer["tau_scaling_factor"]
    corrected = LinearSlopeField(base, wave, 0.1,
                                  [factor * row for row in first_h],
                                  factor * first_m,
                                  temporal_profile=("integrated_bump"
                                                    if profile_only else "linear"))
    frozen = FrozenPotentialField(
        base, wave, 0.1,
        [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
    axis = (-0.45, 0.0, 0.45)
    grid = [(x, y) for x in axis for y in axis]
    angles = np.arange(8) * 2 * np.pi / 8
    fractions = (-0.5, -0.2, 0.0, 0.2, 0.5)
    rows = []
    reference = (json.loads((ROOT / "midplane_wave_transfer_trajectory.json")
                            .read_text())["trajectory"] if profile_only else None)
    for fraction in fractions:
        tau = wave.tau0 + fraction * wave.time_halfwidth
        before = (reference[len(rows)]["frozen"] if profile_only else
                  stats(sample_residual(
                      frozen, wave, tau, grid, angles,
                      time_step=wave.time_halfwidth, time_min=0.5 * 2**-20)))
        after = stats(sample_residual(
            corrected, wave, tau, grid, angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20))
        rows.append(dict(fraction=fraction, tau=tau, frozen=before,
                         corrected=after,
                         max_ratio=after["max"] / before["max"],
                         rms_ratio=after["rms"] / before["rms"]))
        print(f"finished fraction={fraction}: max ratio={rows[-1]['max_ratio']:.3g}",
              flush=True)
    report = dict(source="Direct full nonlinear residual of k11-to-k19 slope continuation",
                  k=19, angle_count=len(angles), heldout_nodes=len(grid),
                  trajectory=rows,
                  temporal_profile=("integrated_bump" if profile_only else "linear"),
                  scope="Interior-time direct finite-difference check of a compact exact-curl correction with transferred coefficients and fixed spatial shape. No evolved amplitude inverse, time endpoint cutoff, moment closure, or uniform-in-scale estimate.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / ("midplane_wave_transfer_envelope.json" if profile_only
                     else "midplane_wave_transfer_trajectory.json")
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    import sys
    run(profile_only="--integrated-bump" in sys.argv)
