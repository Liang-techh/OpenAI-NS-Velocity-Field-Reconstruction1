"""Screen one damped interval correction on the widened-cone exact-curl field.

The first slope is fit at the pulse centre on a 5 x 5 grid. A second
full-momentum residual is fit after advancing the linear trajectory by 0.1
pulse half-width. A 2 x 2 grid selects a common damping alpha, and a
disjoint 4 x 4 grid verifies the selected correction.
"""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from curl_wave_patch_evolution import fit_slope, sample_residual
from curl_wave_prototype import WavePerturbedField
from midplane_wider_cone_second_slope import QuadraticSlopeField
from midplane_wider_cone_slope_projection import (LinearSlopeField,
                                                   build_case, row_stats)
from radial_continuation import ROOT


ALPHAS = (0.0, 0.1, 0.2)
TRAIN_AXIS = (-0.6, -0.3, 0.0, 0.3, 0.6)
# Every coordinate in this selection grid is disjoint from TRAIN_AXIS.
SWEEP_AXIS = (-0.5, 0.5)
# Independent verification coordinates are also disjoint from both grids.
VERIFY_AXIS = (-0.45, -0.15, 0.15, 0.45)
FRACTIONS = (0.05, 0.1)


def grid(axis):
    return [(x, y) for x in axis for y in axis]


def direct_stats(field, wave, tau, points, angles, time_halfwidth):
    rows = sample_residual(field, wave, tau, points, angles,
                           time_step=time_halfwidth,
                           time_min=.5 * 2.**-20)
    return row_stats(rows)


def run():
    angles = np.arange(16) * 2 * np.pi / 16
    train_grid = grid(TRAIN_AXIS)
    sweep_grid = grid(SWEEP_AXIS)
    verify_grid = grid(VERIFY_AXIS)
    train_sweeps = []
    scale_payloads = []

    for k in (11, 19):
        mean, wave, args, _ = build_case(k, angles)
        tau = wave.tau0
        horizon = .1 * args["time_halfwidth"]
        frozen = WavePerturbedField(mean, ScaledWave(wave, .1))

        # Centre fit on the prescribed 5 x 5 training grid.
        rows0 = sample_residual(
            frozen, wave, tau, train_grid, angles,
            time_step=args["time_halfwidth"], time_min=.5 * 2.**-20)
        first_h, first_mean = fit_slope(
            rows0, wave, tau, angles, regularization=1e-3)
        linear = LinearSlopeField(mean, wave, first_h, first_mean)

        # Fit the second correction at +0.1 pulse half-width from the
        # linearly advanced state, as in the prior second-slope experiment.
        tau1 = tau + horizon
        rows1 = sample_residual(
            linear, wave, tau1, train_grid, angles,
            time_step=args["time_halfwidth"], time_min=.5 * 2.**-20)
        delta_h, delta_mean = fit_slope(
            rows1, wave, tau1, angles, regularization=1e-3)

        # Keep every alpha evaluation independent of the verification grid.
        alpha_rows = {}
        for alpha in ALPHAS:
            scaled_h = [alpha * delta for delta in delta_h]
            scaled_mean = alpha * delta_mean
            trajectory = QuadraticSlopeField(
                mean, wave, first_h, first_mean, scaled_h, scaled_mean,
                horizon)
            alpha_rows[str(alpha)] = {}
            for fraction in FRACTIONS:
                check_tau = tau + fraction * args["time_halfwidth"]
                key = str(fraction)
                alpha_rows[str(alpha)][key] = direct_stats(
                    trajectory, wave, check_tau, sweep_grid, angles,
                    args["time_halfwidth"])

        train_sweeps.append(alpha_rows)
        objectives = {
            str(alpha): float(np.mean([
                alpha_rows[str(alpha)][str(fraction)]["rms"] /
                max(alpha_rows["0.0"][str(fraction)]["rms"], 1e-300)
                for fraction in FRACTIONS]))
            for alpha in ALPHAS
        }
        best_alpha = min(ALPHAS, key=lambda alpha: objectives[str(alpha)])
        scale_payloads.append(dict(
            k=k,
            tau=float(tau),
            horizon=float(horizon),
            pulse_time_halfwidth=float(args["time_halfwidth"]),
            second_fit_time_fraction=0.1,
            first_slope_potential_norms=[
                float(np.linalg.norm(f[:27])) for f in first_h],
            first_slope_pressure_norms=[
                float(np.linalg.norm(f[27:])) for f in first_h] + [
                    float(np.linalg.norm(first_mean[18:]))],
            second_slope_relative_potential_norms=[
                float(np.linalg.norm(d[:27]) /
                      max(np.linalg.norm(f[:27]), 1e-30))
                for d, f in zip(delta_h, first_h)],
            second_slope_relative_mean_potential_norm=float(
                np.linalg.norm(delta_mean[:18]) /
                max(np.linalg.norm(first_mean[:18]), 1e-30)),
            selection_alpha_sweep_2x2=alpha_rows,
            best_train_alpha=float(best_alpha),
            best_train_objective=objectives[str(best_alpha)],
            train_objectives=objectives,
            verification_4x4={},
        ))
        print(f"finished k={k}: best selection alpha={best_alpha}", flush=True)

    # One common alpha is selected from scale-normalized training RMS.  This
    # is the requested simultaneous two-scale choice; no holdout information
    # enters the selection.
    common_objectives = {}
    for alpha in ALPHAS:
        terms = []
        for payload, rows in zip(scale_payloads, train_sweeps):
            for fraction in FRACTIONS:
                key = str(fraction)
                terms.append(rows[str(alpha)][key]["rms"] /
                             max(rows["0.0"][key]["rms"], 1e-300))
        common_objectives[str(alpha)] = float(np.mean(terms))
    common_alpha = min(ALPHAS, key=lambda alpha:
                       common_objectives[str(alpha)])

    # Verify only the selected common alpha and the no-second-correction
    # baseline on the independent 4 x 4 grid.
    for payload in scale_payloads:
        k = payload["k"]
        mean, wave, args, _ = build_case(k, angles)
        tau = wave.tau0
        horizon = .1 * args["time_halfwidth"]
        frozen = WavePerturbedField(mean, ScaledWave(wave, .1))
        rows0 = sample_residual(
            frozen, wave, tau, train_grid, angles,
            time_step=args["time_halfwidth"], time_min=.5 * 2.**-20)
        first_h, first_mean = fit_slope(
            rows0, wave, tau, angles, regularization=1e-3)
        linear = LinearSlopeField(mean, wave, first_h, first_mean)
        tau1 = tau + horizon
        rows1 = sample_residual(
            linear, wave, tau1, train_grid, angles,
            time_step=args["time_halfwidth"], time_min=.5 * 2.**-20)
        delta_h, delta_mean = fit_slope(
            rows1, wave, tau1, angles, regularization=1e-3)
        verification = {}
        for alpha in sorted({0.0, float(common_alpha)}):
            trajectory = QuadraticSlopeField(
                mean, wave, first_h, first_mean,
                [alpha * delta for delta in delta_h], alpha * delta_mean,
                horizon)
            verification[str(alpha)] = {}
            for fraction in FRACTIONS:
                check_tau = tau + fraction * args["time_halfwidth"]
                verification[str(alpha)][str(fraction)] = direct_stats(
                    trajectory, wave, check_tau, verify_grid, angles,
                    args["time_halfwidth"])
        payload["verification_4x4"] = verification

    common_training_improves_both = all(
        all(payload["selection_alpha_sweep_2x2"][str(common_alpha)][str(fraction)]["rms"]
            < payload["selection_alpha_sweep_2x2"]["0.0"][str(fraction)]["rms"]
            for fraction in FRACTIONS)
        for payload in scale_payloads)
    common_verification_improves_both = all(
        all(payload["verification_4x4"][str(common_alpha)][str(fraction)]["rms"]
            < payload["verification_4x4"]["0.0"][str(fraction)]["rms"]
            for fraction in FRACTIONS)
        for payload in scale_payloads)

    for payload in scale_payloads:
        verify = payload["verification_4x4"]
        payload["verification_relative_rms_to_alpha0"] = {
            str(fraction): float(
                verify[str(common_alpha)][str(fraction)]["rms"] /
                max(verify["0.0"][str(fraction)]["rms"], 1e-300))
            for fraction in FRACTIONS
        }
        print(
            f"k={payload['k']}: best selection alpha={payload['best_train_alpha']}; "
            f"common alpha={common_alpha}; "
            f"4x4 RMS ratios={payload['verification_relative_rms_to_alpha0']}",
            flush=True)

    report = dict(
        source="One damped second-slope interval correction on widened-cone "
               "exact-curl fields",
        alphas=list(ALPHAS),
        fractions=list(FRACTIONS),
        train_axis=list(TRAIN_AXIS),
        selection_axis=list(SWEEP_AXIS),
        verification_axis=list(VERIFY_AXIS),
        angles_per_node=len(angles),
        alpha_selection=dict(
            rule="Minimize mean over k=11,19 and fractions +0.05,+0.1 of "
                 "2x2 selection RMS normalized by the alpha=0 trajectory",
            common_alpha=float(common_alpha),
            common_objectives=common_objectives,
            common_training_improves_both_scales=common_training_improves_both,
            common_verification_improves_both_scales=
                common_verification_improves_both,
        ),
        scales=scale_payloads,
        scope="At k=11 and k=19, fit the initial exact-curl potential/pressure "
              "slope on a 5 x 5 pulse-centre grid, fit one residual slope at "
              "+0.1 pulse half-width from the linear trajectory, and apply a "
              "single scalar alpha to that second slope and pressure in the "
              "quadratic replay. The 2 x 2 grid selects alpha; "
              "the 4 x 4 grid is independent verification at +0.05 and +0.1. "
              "All values are direct finite-difference full Cartesian momentum "
              "residual norms from sample_residual. This is an interval screen "
              "only: no endpoint control, uniform space-time bound, nonlinear "
              "interval solve, or scale recursion is established.",
        accepted=False,
        scale_recursion_established=False,
    )
    output = ROOT / "midplane_wider_cone_interval_corrector.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
