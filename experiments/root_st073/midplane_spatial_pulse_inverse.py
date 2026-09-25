"""Spatially resolved frozen-principal pulse inverse on two dyadic scales.

Solves the transverse amplitude ODE independently at five r-z patch nodes
for both exact-curl wave harmonics. This is a spatial map of the principal
inverse, not yet a continuous supported amplitude or a corrected field.
"""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import FrozenPotentialField, sample_residual
from curl_wave_prototype import LocalizedCurlWave
from joint_collar_fit import kinematics
from midplane_physical_covariance_pairs import support
from pulse_amplitude_inverse import solve_mode
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def local_k_matrix(base, r, z, tau, time_halfwidth):
    point = np.array([[r, 0., z]])
    hs = .0005 * np.sqrt(base.nu * tau)
    ht = min(.0001 * tau, time_halfwidth / 8)
    velocity, gradient, _ = kinematics(
        base, point, tau, hs, ht, time_min=0.5 * 2**-20)
    v, g = velocity[0], gradient[0]
    swirl_rate = v[1] / r
    matrix = np.array([
        [0., -2. * swirl_rate, 0.],
        [swirl_rate + g[1, 0], 0., 0.],
        [g[2, 0], 0., 0.]])
    return matrix, float(np.linalg.norm(v))


def transfer_metrics(scales):
    first, last = scales
    ratios = [np.array([[mode["amplitude_to_background_speed"]
                         for mode in node["modes"]]
                        for node in scale["nodes"]]) for scale in scales]
    vectors = [np.array([complex(*pair)
                         for node in scale["nodes"]
                         for mode in node["modes"]
                         for pair in mode["midpoint_amplitude"]])
               for scale in scales]
    a, b = ratios
    x, y = vectors
    scalar = np.vdot(x, y) / np.vdot(x, x)
    expected_scale = np.sqrt(first["tau0"] / last["tau0"])
    return dict(ratio_map_relative_l2_drift=float(np.linalg.norm(b-a) /
                                                   np.linalg.norm(a)),
                ratio_map_max_pair_relative_drift=float(np.max(
                    np.abs(b-a) / np.maximum(a, 1e-12))),
                ratio_map_cosine=float(np.vdot(a, b) /
                                       (np.linalg.norm(a)*np.linalg.norm(b))),
                midpoint_vector_best_complex_scalar=[float(scalar.real),
                                                     float(scalar.imag)],
                midpoint_vector_shape_error_after_scalar=float(
                    np.linalg.norm(y-scalar*x)/np.linalg.norm(y)),
                expected_tau_minus_half_scale=float(expected_scale),
                midpoint_vector_tau_minus_half_relative_drift=float(
                    np.linalg.norm(y/expected_scale-x)/np.linalg.norm(x)))


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    grid = [(0., 0.), (-.45, 0.), (.45, 0.), (0., -.45), (0., .45)]
    angles = np.arange(8) * 2 * np.pi / 8
    scales = []
    for k in (11, 19):
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
        tau0 = wave.tau0
        times = np.linspace(0.5 - (tau0 + 0.5 * wave.time_halfwidth),
                            0.5 - (tau0 - 0.5 * wave.time_halfwidth), 9)
        source_values = [[[] for _ in wave.waves] for _ in grid]
        for t in times:
            tau = 0.5 - t
            rows = sample_residual(
                frozen, wave, tau, grid, angles,
                time_step=wave.time_halfwidth,
                time_min=0.5 * 2**-20)
            for node, row in enumerate(rows):
                for j, mode in enumerate(wave.waves):
                    harmonic = 2 * np.mean(
                        row["residual"]
                        * np.exp(-1j * mode["m"] * angles)[:, None], axis=0)
                    source_values[node][j].append(
                        harmonic * np.exp(-1j * mode["omega"] * (tau - tau0)))
        nodes = []
        for index, (xi, eta) in enumerate(grid):
            r = wave.radius + xi * wave.radial_halfwidth
            z = wave.zcenter + eta * wave.axial_halfwidth
            k_matrix, background_speed = local_k_matrix(
                base, r, z, tau0, wave.time_halfwidth)
            modes = []
            for mode, forcing in zip(wave.waves, source_values[index]):
                result = solve_mode(times, forcing, mode["normal"],
                                    k_matrix, base.nu)
                result["angular_mode"] = int(mode["m"])
                result["amplitude_to_background_speed"] = (
                    result["amplitude_max"] / background_speed)
                modes.append(result)
            nodes.append(dict(xi=xi, eta=eta, r=r, z=z,
                              background_speed=background_speed,
                              local_k_matrix=k_matrix.tolist(), modes=modes))
        max_ratio = max(mode["amplitude_to_background_speed"]
                        for node in nodes for mode in node["modes"])
        max_ode_error = max(mode["centered_equation_defect_relative"]
                            for node in nodes for mode in node["modes"])
        scales.append(dict(k=k, tau0=tau0,
                           radial_halfwidth=wave.radial_halfwidth,
                           axial_halfwidth=wave.axial_halfwidth,
                           time_halfwidth=wave.time_halfwidth,
                           nodes=nodes,
                           max_amplitude_to_background_speed=max_ratio,
                           max_centered_ode_error_relative=max_ode_error))
        checkpoint = dict(source="Two-scale five-node frozen-principal transverse pulse inverse",
                          spatial_grid=grid, source_time_samples=9,
                          scales=scales, completed_scales=len(scales),
                          accepted=False, scale_recursion_established=False)
        (ROOT / "midplane_spatial_pulse_inverse.json").write_bytes(
            (json.dumps(checkpoint, indent=2) + "\n").encode())
        print(f"finished k={k}: max amplitude/background={max_ratio:.3g}",
              flush=True)
    report = dict(source="Two-scale five-node frozen-principal transverse pulse inverse",
                  spatial_grid=grid, source_time_samples=9,
                  scales=scales, completed_scales=len(scales),
                  scale_transfer=transfer_metrics(scales),
                  scope="Independent five-node spatial paths, frozen local physical background and fixed carrier normals; the source comes from the full frozen-wave Cartesian residual. Equation (7.13)-inspired transverse ODE and normal pressure only. No spatial interpolation, moving phase, supported potential reconstruction, full post-correction momentum, mean/stress/moment cycle, or scale-recursion proof.",
                  accepted=False, scale_recursion_established=False)
    (ROOT / "midplane_spatial_pulse_inverse.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
