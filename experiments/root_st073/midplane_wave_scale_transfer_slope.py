"""Transfer one instantaneous curl/pressure/mean slope between scales.

The k=11 fit is scaled by tau_11/tau_19 and evaluated at k=19 without
refitting. This tests a finite-grid normalized correction template, not
an evolved field or a complete residual-improvement cycle.
"""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import (FrozenPotentialField, fit_slope,
                                       metrics, sample_residual)
from curl_wave_prototype import LocalizedCurlWave
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def encode(harmonic, mean):
    return dict(harmonic=[[[float(x.real), float(x.imag)] for x in row]
                          for row in harmonic], mean=mean.tolist())


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json")
                       .read_text())
    angles = np.arange(16) * 2 * np.pi / 16
    train_axis = (-0.6, -0.2, 0.2, 0.6)
    heldout_axis = (-0.45, 0.0, 0.45)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    heldout_grid = [(x, y) for x in heldout_axis for y in heldout_axis]
    rows = []
    for k in (11, 19):
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{k}.json").read_text())
        source["nu"] = base.nu
        chosen = next(row["selected"] for row in pairs["scales"]
                      if row["k"] == k)
        wave = LocalizedCurlWave(
            source, pulse_indices=chosen["indices"],
            **support(inner, source, base.nu))
        wave.weights = np.asarray(chosen["positive_weights"], float)
        tau = float(source["tau"])
        field = FrozenPotentialField(
            base, wave, 0.1,
            [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
        train = sample_residual(
            field, wave, tau, train_grid, angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20)
        heldout = sample_residual(
            field, wave, tau, heldout_grid, angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20)
        harmonic, mean = fit_slope(train, wave, tau, angles)
        rows.append(dict(k=k, tau=tau, wave=wave, train=train,
                         heldout=heldout, harmonic=harmonic, mean=mean,
                         own_train=metrics(train, wave, tau, angles,
                                           harmonic, mean),
                         own_holdout=metrics(heldout, wave, tau, angles,
                                             harmonic, mean)))
        print(f"fitted k={k}", flush=True)
    first, last = rows
    scale = first["tau"] / last["tau"]
    transferred_h = [scale * value for value in first["harmonic"]]
    transferred_m = scale * first["mean"]
    transfer_train = metrics(last["train"], last["wave"], last["tau"],
                             angles, transferred_h, transferred_m)
    transfer_holdout = metrics(last["heldout"], last["wave"], last["tau"],
                               angles, transferred_h, transferred_m)
    first_coeff = np.r_[np.concatenate([x.real for x in first["harmonic"]]),
                        np.concatenate([x.imag for x in first["harmonic"]]),
                        first["mean"]]
    last_coeff = np.r_[np.concatenate([x.real for x in last["harmonic"]]),
                       np.concatenate([x.imag for x in last["harmonic"]]),
                       last["mean"]]
    normalized_coefficient_drift = float(np.linalg.norm(
        last["tau"] * last_coeff - first["tau"] * first_coeff)
        / np.linalg.norm(first["tau"] * first_coeff))
    report = dict(source="Instantaneous spatial wave/pressure/mean slope transfer",
                  scales=[dict(k=row["k"], tau=row["tau"],
                               own_train=row["own_train"],
                               own_holdout=row["own_holdout"],
                               fitted_coefficients=encode(row["harmonic"],
                                                          row["mean"]))
                          for row in rows],
                  tau_scaling_factor=scale,
                  normalized_coefficient_drift=normalized_coefficient_drift,
                  transferred_k19_train=transfer_train,
                  transferred_k19_holdout=transfer_holdout,
                  scope="One instant and 16 train/9 held-out radial-axial nodes x16 angles per scale. Fits complete frozen-wave Fourier residual to a compact potential time derivative and pressure/mean basis. k=11 coefficients transferred with tau^-1 scale to k=19 without refitting. Projection only: no continuous-time exact-curl field, full nonlinear post-correction residual, support/moment closure, or PDE acceptance.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "midplane_wave_scale_transfer_slope.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          own_holdout=[row["own_holdout"] for row in rows],
                          transferred_k19_holdout=transfer_holdout,
                          normalized_coefficient_drift=(
                              normalized_coefficient_drift)), indent=2))
    return report


if __name__ == "__main__":
    run()
