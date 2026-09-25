"""Implicit midpoint evolution in a transferred three-template potential space.

This is a low-rank numerical test of whether time integration, rather than
explicit Euler, controls the midplane wave defect. It is not the supported
pulse inverse of Proposition 7.2.
"""

import json

import numpy as np
from scipy.optimize import root

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_collocation import basis
from curl_wave_patch_evolution import FrozenPotentialField, sample_residual
from curl_wave_patch_mean import mean_basis
from curl_wave_patch_trajectory import FirstIntervalField, residual_metrics
from curl_wave_prototype import LocalizedCurlWave
from midplane_physical_covariance_pairs import support
from midplane_wave_short_interval import stage_data
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def decode(row):
    harmonics = [np.array([complex(*pair) for pair in values])
                 for values in row["harmonic"]]
    return harmonics, np.asarray(row["mean"], float)


def state_field(base, wave, slopes_h, slope_m, q):
    h = wave.time_halfwidth
    return FrozenPotentialField(
        base, wave, 0.1,
        [h * q[j] * slopes_h[j][:27] for j in range(len(wave.waves))],
        h * q[-1] * slope_m[:18])


def fit_restricted(rows, wave, tau, angles, slopes_h, slope_m):
    """Fit three pulse-coordinate derivatives and all compact pressures."""
    blocks, targets = [], []
    for row in rows:
        r, z = row["r"], row["z"]
        n = len(angles)
        matrix = np.zeros((n, 3, 3 + 18 * len(wave.waves) + 9))
        for j, (mode, slope) in enumerate(zip(wave.waves, slopes_h)):
            phase = np.exp(1j * (mode["omega"] * (tau - wave.tau0)
                                 + mode["m"] * angles))
            columns = basis(wave, mode, r, z)
            matrix[:, :, j] = (phase[:, None]
                               * (columns[:, :27] @ slope[:27])[None, :]).real
            p = phase[:, None, None] * columns[:, 27:][None, :, :]
            start = 3 + 18 * j
            matrix[:, :, start:start+9] = p.real
            matrix[:, :, start+9:start+18] = -p.imag
        mean_columns = mean_basis(wave, r, z)
        matrix[:, :, 2] = mean_columns[:, :18] @ slope_m[:18]
        matrix[:, :, -9:] = mean_columns[:, 18:][None, :, :]
        blocks.append(matrix.reshape(n * 3, -1))
        targets.append(-row["residual"].reshape(-1))
    matrix = np.vstack(blocks)
    target = np.concatenate(targets)
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    normalized = matrix / scale
    ridge = 1e-4 * np.linalg.norm(normalized, ord=2)
    augmented = np.vstack((normalized, ridge * np.eye(matrix.shape[1])))
    solution = np.linalg.lstsq(
        augmented, np.r_[target, np.zeros(matrix.shape[1])], rcond=None)[0]
    coefficients = solution / scale
    projected = target - matrix @ coefficients
    pressure_h = [coefficients[3+18*j:3+18*j+9]
                  + 1j * coefficients[3+18*j+9:3+18*j+18]
                  for j in range(len(wave.waves))]
    return (coefficients[:3], pressure_h, coefficients[-9:],
            float(np.linalg.norm(projected) / np.linalg.norm(target)))


def run(train_order=2):
    if train_order not in (2, 3):
        raise ValueError("train_order must be 2 or 3")
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    source = json.loads((ROOT / "compact_potential" /
                         "midplane_wave_source_k19.json").read_text())
    source["nu"] = base.nu
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    selected = next(row["selected"] for row in pairs["scales"]
                    if row["k"] == 19)
    wave = LocalizedCurlWave(source, pulse_indices=selected["indices"],
                             **support(inner, source, base.nu))
    wave.weights = np.asarray(selected["positive_weights"], float)
    transfer = json.loads((ROOT / "midplane_wave_scale_transfer_slope.json").read_text())
    slopes_h, slope_m = decode(transfer["scales"][0]["fitted_coefficients"])
    factor = transfer["tau_scaling_factor"]
    slopes_h = [factor * x for x in slopes_h]
    slope_m = factor * slope_m
    train_axis = (-0.35, 0.35) if train_order == 2 else (-0.35, 0., 0.35)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    angles = np.arange(8) * 2 * np.pi / 8
    held_axis = (-0.45, 0., 0.45)
    held_angles = (np.arange(8) + 0.5) * 2 * np.pi / 8
    held_points = np.vstack([
        np.column_stack((r * np.cos(held_angles), r * np.sin(held_angles),
                         np.full(len(held_angles), z)))
        for x in held_axis for y in held_axis
        for r, z in [(wave.radius + x * wave.radial_halfwidth,
                      wave.zcenter + y * wave.axial_halfwidth)]])
    fractions = (-0.1, 0.0, 0.1)
    q = np.zeros(3)
    frozen = state_field(base, wave, slopes_h, slope_m, q)
    intervals = []
    for index, (left, right) in enumerate(zip(fractions[:-1], fractions[1:])):
        mid = (left + right) / 2
        tau = wave.tau0 + mid * wave.time_halfwidth
        dt = (right - left) * wave.time_halfwidth
        cache = {}

        def fit(qmid):
            key = tuple(qmid)
            if key not in cache:
                rows = sample_residual(
                    state_field(base, wave, slopes_h, slope_m, qmid),
                    wave, tau, train_grid, angles, time_step=dt,
                    time_min=0.5 * 2**-20)
                cache[key] = fit_restricted(rows, wave, tau, angles,
                                            slopes_h, slope_m)
            return cache[key]

        def equation(qright):
            derivative = fit((q + qright) / 2)[0]
            return qright - q - (right - left) * derivative

        solved = root(equation, q, method="hybr", options={"xtol": 1e-6})
        qright = solved.x
        qmid = (q + qright) / 2
        derivative, pressure_h, pressure_m, projection_error = fit(qmid)
        actual_derivative = (qright - q) / (right - left)
        harmonic = [np.r_[actual_derivative[j] * slopes_h[j][:27],
                              pressure_h[j]] for j in range(len(wave.waves))]
        mean = np.r_[actual_derivative[-1] * slope_m[:18], pressure_m]
        interval = FirstIntervalField(
            base, wave, 0.1,
            stage_data(wave.tau0 + left * wave.time_halfwidth,
                       harmonic, mean),
            initial_harmonic=[wave.time_halfwidth * q[j] * slopes_h[j][:27]
                              for j in range(len(wave.waves))],
            initial_mean=wave.time_halfwidth * q[-1] * slope_m[:18])
        before = residual_metrics(frozen, held_points, tau, dt)
        after = residual_metrics(interval, held_points, tau, dt)
        intervals.append(dict(midpoint_fraction=mid,
                              q_left=q.tolist(), q_right=qright.tolist(),
                              implicit_success=bool(solved.success),
                              implicit_message=solved.message,
                              evaluations=len(cache),
                              fixed_point_defect=float(np.linalg.norm(equation(qright))),
                              train_projection_relative=projection_error,
                              frozen=before, corrected=after,
                              max_ratio=after["max_momentum"] /
                              before["max_momentum"]))
        q = qright
        print(f"finished interval {index}: success={solved.success}, ratio={intervals[-1]['max_ratio']:.3g}, evaluations={len(cache)}", flush=True)
        if not solved.success:
            break
    report = dict(source="Three-template implicit-midpoint midplane pulse trial",
                  k=19, pulse_fractions=fractions,
                  train_nodes=len(train_grid), heldout_nodes=9,
                  intervals=intervals,
                  scope="Low-rank transferred potential space, full frozen-state nonlinear residual for each implicit RHS, algebraic harmonic/mean pressure fit, direct midpoint full-momentum holdouts. No full amplitude inverse, support closure, moment restoration, or scale recursion.",
                  accepted=False, scale_recursion_established=False)
    output = (ROOT / ("midplane_wave_implicit_template.json"
                      if train_order == 2 else
                      "midplane_wave_implicit_template_9train.json"))
    output.write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-order", type=int, default=2)
    run(parser.parse_args().train_order)
