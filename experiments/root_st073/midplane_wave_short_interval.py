"""Short physical-time coefficient march for the k=19 midplane pulse.

The complete frozen-state residual determines a new compact exact-curl
potential slope and pressure at every stage. Direct midpoint checks decide
whether this provisional collocation ODE helps beyond its fit nodes.
"""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, sample_residual)
from curl_wave_patch_trajectory import FirstIntervalField, residual_metrics
from curl_wave_prototype import LocalizedCurlWave
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def stage_data(tau, harmonic, mean):
    return dict(tau=tau,
                harmonic_slope_coefficients=[
                    [[float(value.real), float(value.imag)] for value in row]
                    for row in harmonic],
                mean_slope_coefficients=mean.tolist())


def run(steps=2):
    if steps < 2:
        raise ValueError("At least two intervals are required")
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
    train_axis = (-0.35, 0.35)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    held_axis = (-0.45, 0.0, 0.45)
    held_grid = [(x, y) for x in held_axis for y in held_axis]
    train_angles = np.arange(8) * 2 * np.pi / 8
    held_angles = (np.arange(8) + 0.5) * 2 * np.pi / 8
    held_points = np.vstack([
        np.column_stack((r * np.cos(held_angles), r * np.sin(held_angles),
                         np.full(len(held_angles), z)))
        for x, y in held_grid
        for r, z in [(wave.radius + x * wave.radial_halfwidth,
                      wave.zcenter + y * wave.axial_halfwidth)]])
    fractions = np.linspace(-0.1, 0.1, steps + 1).tolist()
    times = [wave.tau0 + f * wave.time_halfwidth for f in fractions]
    dt = times[1] - times[0]
    harmonic_state = [np.zeros(27, complex) for _ in wave.waves]
    mean_state = np.zeros(18)
    frozen = FrozenPotentialField(base, wave, 0.1,
                                  [np.zeros(27, complex) for _ in wave.waves],
                                  np.zeros(18))
    stages = []
    for index, tau in enumerate(times[:-1]):
        current = FrozenPotentialField(base, wave, 0.1,
                                       harmonic_state, mean_state)
        samples = sample_residual(current, wave, tau, train_grid,
                                  train_angles, time_step=dt,
                                  time_min=0.5 * 2**-20)
        harmonic, mean = fit_slope(samples, wave, tau, train_angles)
        interval = FirstIntervalField(
            base, wave, 0.1, stage_data(tau, harmonic, mean),
            initial_harmonic=harmonic_state, initial_mean=mean_state)
        midpoint = tau + dt / 2
        before = residual_metrics(frozen, held_points, midpoint,
                                  dt)
        after = residual_metrics(interval, held_points, midpoint,
                                 dt)
        stages.append(dict(start_fraction=fractions[index],
                           midpoint_fraction=(fractions[index] + fractions[index+1])/2,
                           start_train=metrics(samples, wave, tau,
                                               train_angles, harmonic, mean),
                           heldout_midpoint_frozen=before,
                           heldout_midpoint_evolved=after,
                           midpoint_max_ratio=(after["max_momentum"] /
                                               before["max_momentum"]),
                           slope_norms=[float(np.linalg.norm(row[:27]))
                                        for row in harmonic],
                           mean_slope_norm=float(np.linalg.norm(mean[:18]))))
        harmonic_state = [state + dt * slope[:27]
                          for state, slope in zip(harmonic_state, harmonic)]
        mean_state = mean_state + dt * mean[:18]
        print(f"finished interval {index}: midpoint ratio={stages[-1]['midpoint_max_ratio']:.3g}",
              flush=True)
    report = dict(source="Explicit time-dependent compact potential and pressure fit",
                  k=19, pulse_fractions=fractions, time_step=dt,
                  steps=steps,
                  train_nodes=len(train_grid), heldout_nodes=len(held_grid),
                  stages=stages,
                  scope="A provisional physical-time collocation ODE. Pressure is recomputed at stage starts and held fixed within each interval. Full nonlinear momentum is evaluated directly at held-out interior midpoints. No stiff convergence, pulse endpoint support, coupled stress/moment correction, or uniform scale recursion.",
                  accepted=False, scale_recursion_established=False)
    output = (ROOT / ("midplane_wave_short_interval.json" if steps == 2
                      else f"midplane_wave_short_interval_{steps}step.json"))
    output.write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=2)
    run(parser.parse_args().steps)
