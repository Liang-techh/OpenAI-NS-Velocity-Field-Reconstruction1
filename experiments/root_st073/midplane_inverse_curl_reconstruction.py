"""Reconstruct a compact exact-curl midpoint field from spatial pulse inverses.

Five node amplitudes and their principal time derivatives are interpolated
by a low-degree vector potential. The resulting full Cartesian momentum is
screened directly on independent spatial and angular nodes at two scales.
"""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_collocation import basis
from curl_wave_patch_mean import mean_basis
from curl_wave_patch_evolution import (FrozenPotentialField, polynomial_data,
                                       sample_residual)
from curl_wave_prototype import LocalizedCurlWave, bump
from joint_collar_fit import kinematics
from midplane_physical_covariance_pairs import support
from midplane_spatial_pulse_inverse import local_k_matrix
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


POLYNOMIAL_INDICES = (0, 3, 6, 1, 2)  # 1, xi, xi^2, eta, eta^2.
POTENTIAL_COLUMNS = tuple(9 * component + index
                          for component in range(3)
                          for index in POLYNOMIAL_INDICES)


def ridge(matrix, target, weight=1e-4, physical=False):
    if physical:
        penalty = weight * np.linalg.norm(matrix, ord=2)
        augmented = np.vstack((matrix, penalty * np.eye(matrix.shape[1])))
        solved = np.linalg.lstsq(
            augmented, np.r_[target, np.zeros(matrix.shape[1], complex)],
            rcond=None)[0]
    else:
        scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
        normalized = matrix / scale
        augmented = np.vstack((normalized,
                               weight * np.linalg.norm(normalized, ord=2)
                               * np.eye(matrix.shape[1])))
        solved = np.linalg.lstsq(
            augmented, np.r_[target, np.zeros(matrix.shape[1], complex)],
            rcond=None)[0] / scale
    return solved, float(np.linalg.norm(matrix @ solved - target)
                         / np.linalg.norm(target))


def fit_vector_potential(wave, mode, nodes, values, full_polynomial=False,
                         physical_ridge=False):
    columns = tuple(range(27)) if full_polynomial else POTENTIAL_COLUMNS
    matrix = np.vstack([-basis(wave, mode, node["r"], node["z"])
                        [:, columns] for node in nodes])
    target = np.concatenate(values)
    fitted, error = ridge(matrix, target, weight=(.05 if physical_ridge
                                                   else 1e-4),
                          physical=physical_ridge)
    full = np.zeros(27, complex)
    full[list(columns)] = fitted
    return full, error


def pressure_values(wave, mode, r, z):
    xi = (r - wave.radius) / wave.radial_halfwidth
    eta = (z - wave.zcenter) / wave.axial_halfwidth
    br = bump(r, wave.radius, wave.radial_halfwidth)[0]
    bz = bump(z, wave.zcenter, wave.axial_halfwidth)[0]
    kr, _, kz = mode["normal"]
    phase = np.exp(1j * (kr * (r - wave.radius) + kz * (z - wave.zcenter)))
    return phase * br * bz * np.array((1., xi, xi**2, eta, eta**2))


def fit_pressure(wave, mode, nodes, values):
    matrix = np.vstack([pressure_values(wave, mode, node["r"], node["z"])
                        for node in nodes])
    fitted, error = ridge(matrix, np.asarray(values))
    full = np.zeros(9, complex)
    full[list(POLYNOMIAL_INDICES)] = fitted
    return full, error


def fit_pressure_gradients(rows, wave, angles):
    harmonic = []
    for mode in wave.waves:
        matrix = np.vstack([basis(wave, mode, row["r"], row["z"])[:, 27:]
                            for row in rows])
        target = -np.concatenate([
            2 * np.mean(row["residual"]
                        * np.exp(-1j*mode["m"]*angles)[:, None], axis=0)
            for row in rows])
        fitted, _ = ridge(matrix, target, weight=.05, physical=True)
        harmonic.append(fitted)
    matrix = np.vstack([mean_basis(wave, row["r"], row["z"])[:, 18:]
                        for row in rows])
    target = -np.concatenate([row["residual"].mean(axis=0)
                              for row in rows])
    mean, _ = ridge(matrix, target, weight=.05, physical=True)
    return harmonic, mean.real


class ReconstructedField:
    """Compact spatial curl correction with local physical-time Taylor jet."""

    def __init__(self, base, wave, potentials, physical_slopes, pressures,
                 mean_pressure=None):
        self.base, self.wave = base, wave
        self.potentials, self.physical_slopes = potentials, physical_slopes
        self.pressures = pressures
        self.mean_pressure = (np.zeros(9) if mean_pressure is None
                              else mean_pressure)
        self.nu = base.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        if not np.all(ts == ts[0]):
            raise ValueError("Evaluate one physical time per spatial batch")
        dtau = float(ts[0] - self.wave.tau0)
        field = FrozenPotentialField(
            self.base, self.wave, 0.1,
            [a - dtau * adot for a, adot in zip(
                self.potentials, self.physical_slopes)], np.zeros(18))
        velocity, pressure = field.fields(pts, ts)
        for i, ((x, y, z), t) in enumerate(zip(pts, ts)):
            r = np.hypot(x, y)
            if (r == 0 or abs(r - self.wave.radius) >= self.wave.radial_halfwidth
                    or abs(z - self.wave.zcenter) >= self.wave.axial_halfwidth):
                continue
            theta = np.arctan2(y, x)
            for mode, coeff in zip(self.wave.waves, self.pressures):
                kr, _, kz = mode["normal"]
                phase = (mode["m"] * theta + kr * (r - self.wave.radius)
                         + kz * (z - self.wave.zcenter)
                         + mode["omega"] * (t - self.wave.tau0))
                q = polynomial_data(self.wave, r, z, coeff, 1)[0][0]
                pressure[i] += (np.exp(1j * phase) * q).real
            pressure[i] += polynomial_data(
                self.wave, r, z, self.mean_pressure, 1)[0][0]
        return velocity, pressure


def stats(rows):
    norms = np.concatenate([np.linalg.norm(row["residual"], axis=1)
                            for row in rows])
    return dict(max=float(norms.max()),
                rms=float(np.sqrt(np.mean(norms**2))))


def peak_term_budget(field, wave, grid, angles):
    points = np.vstack([
        np.column_stack((r*np.cos(angles), r*np.sin(angles),
                         np.full(len(angles), z)))
        for xi, eta in grid
        for r, z in [(wave.radius + xi*wave.radial_halfwidth,
                      wave.zcenter + eta*wave.axial_halfwidth)]])
    tau = wave.tau0
    hs = .0005*np.sqrt(field.nu*tau)
    ht = min(.0001*tau, wave.time_halfwidth/8)
    u, grad, part, terms = kinematics(
        field, points, tau, hs, ht, time_min=0.5*2**-20,
        return_terms=True)
    residual = part + terms["convection"]
    peak = int(np.argmax(np.linalg.norm(residual, axis=1)))
    return dict(peak_point=points[peak].tolist(),
                peak_residual=float(np.linalg.norm(residual[peak])),
                peak_speed=float(np.linalg.norm(u[peak])),
                components={name: float(np.linalg.norm(value[peak]))
                            for name, value in terms.items()})


def run(full_polynomial=False, physical_ridge=False,
        gradient_pressure=False):
    if gradient_pressure and not (full_polynomial and physical_ridge):
        raise ValueError("Gradient-pressure screen uses full regularized potential")
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    inverse = json.loads((ROOT / "midplane_spatial_pulse_inverse.json").read_text())
    train_angles = np.arange(16) * 2 * np.pi / 16
    held_angles = (np.arange(16) + .5) * 2 * np.pi / 16
    held_axis = (-.3, 0., .3)
    held_grid = [(x, y) for x in held_axis for y in held_axis
                 if (x, y) != (0., 0.)]
    scales = []
    for scale in inverse["scales"]:
        k = scale["k"]
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{k}.json").read_text())
        source["nu"] = base.nu
        selected = next(row["selected"] for row in pairs["scales"]
                        if row["k"] == k)
        wave = LocalizedCurlWave(source, pulse_indices=selected["indices"],
                                 **support(inner, source, base.nu))
        wave.weights = np.asarray(selected["positive_weights"], float)
        frozen = FrozenPotentialField(
            base, wave, 0.1,
            [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
        grid = [(node["xi"], node["eta"]) for node in scale["nodes"]]
        source_rows = sample_residual(
            frozen, wave, wave.tau0, grid, train_angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20)
        potentials, slopes, pressures, mode_fits = [], [], [], []
        for j, mode in enumerate(wave.waves):
            normal = np.asarray(mode["normal"])
            projection = np.eye(3) - np.outer(normal, normal) / np.dot(normal, normal)
            source_harmonics = [2 * np.mean(
                row["residual"]
                * np.exp(-1j * mode["m"] * train_angles)[:, None], axis=0)
                for row in source_rows]
            amplitudes, derivatives, pressure_targets = [], [], []
            for node, forcing in zip(scale["nodes"], source_harmonics):
                amplitude = np.array([complex(*pair)
                                      for pair in node["modes"][j]["midpoint_amplitude"]])
                k_matrix = np.asarray(node["local_k_matrix"])
                derivative = (-projection @ (k_matrix @ amplitude + forcing)
                              - base.nu * np.dot(normal, normal) * amplitude)
                pressure = 1j * np.dot(normal, k_matrix @ amplitude + forcing) \
                    / np.dot(normal, normal)
                amplitudes.append(amplitude)
                derivatives.append(derivative)
                pressure_targets.append(pressure)
            potential, amplitude_fit = fit_vector_potential(
                wave, mode, scale["nodes"], amplitudes,
                full_polynomial=full_polynomial,
                physical_ridge=physical_ridge)
            slope, derivative_fit = fit_vector_potential(
                wave, mode, scale["nodes"], derivatives,
                full_polynomial=full_polynomial,
                physical_ridge=physical_ridge)
            pressure, pressure_fit = fit_pressure(
                wave, mode, scale["nodes"], pressure_targets)
            potentials.append(potential)
            slopes.append(slope)
            pressures.append(pressure)
            mode_fits.append(dict(angular_mode=int(mode["m"]),
                                  amplitude_relative_error=amplitude_fit,
                                  derivative_relative_error=derivative_fit,
                                  pressure_relative_error=pressure_fit,
                                  potential_coefficient_norm=float(np.linalg.norm(potential)),
                                  physical_slope_coefficient_norm=float(np.linalg.norm(slope)),
                                  max_target_amplitude=float(max(
                                      np.linalg.norm(a) for a in amplitudes))))
        reconstructed = ReconstructedField(base, wave, potentials,
                                           slopes, pressures)
        if gradient_pressure:
            pressure_grid_axis = (-.6, -.2, .2, .6)
            pressure_grid = [(x, y) for x in pressure_grid_axis
                             for y in pressure_grid_axis]
            without_pressure = ReconstructedField(
                base, wave, potentials, slopes,
                [np.zeros(9, complex) for _ in wave.waves])
            pressure_rows = sample_residual(
                without_pressure, wave, wave.tau0, pressure_grid,
                train_angles, time_step=wave.time_halfwidth,
                time_min=0.5*2**-20)
            gradient_h, gradient_mean = fit_pressure_gradients(
                pressure_rows, wave, train_angles)
            reconstructed = ReconstructedField(
                base, wave, potentials, slopes, gradient_h,
                mean_pressure=gradient_mean)
        before_train = stats(source_rows)
        after_train = stats(sample_residual(
            reconstructed, wave, wave.tau0, grid, train_angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20))
        before_held = stats(sample_residual(
            frozen, wave, wave.tau0, held_grid, held_angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20))
        after_held = stats(sample_residual(
            reconstructed, wave, wave.tau0, held_grid, held_angles,
            time_step=wave.time_halfwidth, time_min=0.5 * 2**-20))
        scales.append(dict(k=k, mode_fits=mode_fits,
                           train_frozen=before_train,
                           train_reconstructed=after_train,
                           heldout_frozen=before_held,
                           heldout_reconstructed=after_held,
                           heldout_max_ratio=after_held["max"] / before_held["max"],
                           heldout_rms_ratio=after_held["rms"] / before_held["rms"],
                           peak_term_budget=(peak_term_budget(
                               reconstructed, wave, held_grid, held_angles)
                               if physical_ridge else None)))
        print(f"finished k={k}: heldout max ratio={scales[-1]['heldout_max_ratio']:.3g}",
              flush=True)
    report = dict(source="Exact-curl midpoint reconstruction from five-node principal inverse",
                  full_polynomial=full_polynomial,
                  physical_ridge=physical_ridge,
                  gradient_pressure=gradient_pressure,
                  pressure_method=("full-momentum gradient projection"
                                   if gradient_pressure else
                                   "five-node normal-pressure value interpolation"),
                  pressure_training_nodes=(16 if gradient_pressure else 5),
                  heldout_grid=held_grid, heldout_angles=len(held_angles),
                  scales=scales,
                  scope="Compact spatial exact-curl polynomial reconstruction of inverse midpoint amplitudes and their frozen-principal physical-time derivatives; pressure either five-node normal-value interpolation or separate full-momentum gradient projection. Direct full Cartesian momentum at center time only. No time endpoint cutoff, moving normal, full mean/stress/moment cycle, uniform spatial residual, or scale recursion.",
                  accepted=False, scale_recursion_established=False)
    output_name = ("midplane_inverse_curl_reconstruction_gradient_pressure.json"
                   if gradient_pressure else
                   ("midplane_inverse_curl_reconstruction_regularized.json"
                   if physical_ridge else
                   ("midplane_inverse_curl_reconstruction_full.json"
                    if full_polynomial else
                    "midplane_inverse_curl_reconstruction.json")))
    output = ROOT / output_name
    output.write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-polynomial", action="store_true")
    parser.add_argument("--physical-ridge", action="store_true")
    parser.add_argument("--gradient-pressure", action="store_true")
    args = parser.parse_args()
    run(args.full_polynomial, args.physical_ridge,
        args.gradient_pressure)
