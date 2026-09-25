"""Frozen-background principal amplitude inverse for an actual ST073 harmonic.

This uses the transverse projection and normal pressure identity in the
paper's (7.13), but with a frozen physical background/phase normal. It is
an exploratory local inverse, not the paper's supported pulse theorem.
"""

import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_patch_evolution import FrozenPotentialField, sample_residual
from curl_wave_prototype import LocalizedCurlWave
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def complex_spline(times, values):
    return (CubicSpline(times, values.real, axis=0),
            CubicSpline(times, values.imag, axis=0))


def solve_mode(times, source, normal, k_matrix, viscosity):
    """Solve the physical-time transverse principal equation with zero data."""
    normal = np.asarray(normal, float)
    k_matrix = np.asarray(k_matrix, float)
    source = np.asarray(source, complex)
    n2 = float(np.dot(normal, normal))
    projection = np.eye(3) - np.outer(normal, normal) / n2
    real_spline, imag_spline = complex_spline(times, source)

    def forcing(t):
        return real_spline(t) + 1j * imag_spline(t)

    def rhs(t, amplitude):
        return -projection @ (k_matrix @ amplitude + forcing(t)) \
               - viscosity * n2 * amplitude

    solved = solve_ivp(rhs, (float(times[0]), float(times[-1])),
                       np.zeros(3, complex), method="BDF", dense_output=True,
                       rtol=1e-8, atol=1e-7)
    if not solved.success:
        raise RuntimeError(solved.message)
    sample_times = np.linspace(times[0], times[-1], 65)
    amplitudes = solved.sol(sample_times).T
    pressure = np.asarray([
        1j * np.dot(normal, k_matrix @ a + forcing(t)) / n2
        for t, a in zip(sample_times, amplitudes)])
    # Independent centered finite difference checks the integrated path.
    derivative = (amplitudes[2:] - amplitudes[:-2]) / (
        sample_times[2:] - sample_times[:-2])[:, None]
    residual = np.asarray([
        da + k_matrix @ a + viscosity * n2 * a
        + 1j * normal * p + forcing(t)
        for t, a, p, da in zip(sample_times[1:-1], amplitudes[1:-1],
                               pressure[1:-1], derivative)])
    scale = max(float(np.max(np.linalg.norm(source, axis=1))), 1.0)
    return dict(success=True, solver_steps=len(solved.t),
                source_max=scale,
                amplitude_max=float(np.max(np.linalg.norm(amplitudes, axis=1))),
                final_amplitude_norm=float(np.linalg.norm(amplitudes[-1])),
                pressure_max=float(np.max(np.abs(pressure))),
                transverse_constraint_max=float(np.max(np.abs(
                    amplitudes @ normal))),
                centered_equation_defect_max=float(np.max(np.linalg.norm(
                    residual, axis=1))),
                centered_equation_defect_relative=float(np.max(
                    np.linalg.norm(residual, axis=1)) / scale),
                endpoint_amplitude=[[float(x.real), float(x.imag)]
                                    for x in amplitudes[-1]],
                normal=normal.tolist())


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    source_data = json.loads((ROOT / "compact_potential" /
                              "adaptive_bridge_wave_source.json").read_text())
    source_data["nu"] = base.nu
    tau0 = float(source_data["tau"])
    ri = float(inner.from_similarity([inner.p.X_max], [-0.2], tau0)[0, 0])
    wave = LocalizedCurlWave(
        source_data, radial_halfwidth=0.0047 * 15.0 * ri,
        axial_halfwidth=1.0e-4, time_halfwidth=0.04 * tau0,
        min_tau=0.5 * 2.0**-20)
    frozen = FrozenPotentialField(
        base, wave, 0.1,
        [np.zeros(27, complex) for _ in wave.waves], np.zeros(18))
    radius = float(source_data["point"][0])
    velocity = np.asarray(source_data["velocity"], float)
    gradient = np.asarray(source_data["gradient"], float)
    swirl_rate = velocity[1] / radius
    k_matrix = np.array([
        [0.0, -2.0 * swirl_rate, 0.0],
        [2.0 * swirl_rate + gradient[1, 0] - swirl_rate, 0.0, 0.0],
        [gradient[2, 0], 0.0, 0.0]])
    half_window = 0.5 * wave.time_halfwidth
    times = np.linspace(0.5 - (tau0 + half_window),
                        0.5 - (tau0 - half_window), 17)
    angles = np.arange(16) * 2.0 * np.pi / 16.0
    sources = [[] for _ in wave.waves]
    for t in times:
        tau = 0.5 - t
        rows = sample_residual(frozen, wave, tau, [(0.0, 0.0)], angles,
                               time_step=1e-9, time_min=0.5 * 2.0**-20)
        residual = rows[0]["residual"]
        for j, mode in enumerate(wave.waves):
            harmonic = 2 * np.mean(
                residual * np.exp(-1j * mode["m"] * angles)[:, None],
                axis=0)
            sources[j].append(harmonic * np.exp(
                -1j * mode["omega"] * (tau - tau0)))
    background_speed = float(np.linalg.norm(source_data["velocity"]))
    results = []
    for mode, harmonic_source in zip(wave.waves, sources):
        result = solve_mode(times, harmonic_source, mode["normal"],
                            k_matrix, base.nu)
        result.update(angular_mode=int(mode["m"]),
                      carrier_omega=float(mode["omega"]),
                      carrier_norm_times_radial_halfwidth=float(
                          np.linalg.norm(mode["normal"])
                          * wave.radial_halfwidth),
                      carrier_norm_times_axial_halfwidth=float(
                          np.linalg.norm(mode["normal"])
                          * wave.axial_halfwidth),
                      radial_phase_change_halfwidth=float(
                          abs(mode["normal"][0])
                          * wave.radial_halfwidth),
                      axial_phase_change_halfwidth=float(
                          abs(mode["normal"][2])
                          * wave.axial_halfwidth),
                      viscous_damping_over_pulse_halfwidth=float(
                          base.nu * np.dot(mode["normal"], mode["normal"])
                          * wave.time_halfwidth),
                      amplitude_to_background_speed=float(
                          result["amplitude_max"] / background_speed))
        results.append(result)
    report = dict(source="Frozen-background transverse amplitude inverse for actual late-bridge Fourier residual",
                  tau0=tau0, physical_time_interval=[float(times[0]),
                                                       float(times[-1])],
                  source_time_samples=len(times),
                  radial_axial_point=[float(wave.radius),
                                      float(wave.zcenter)],
                  background_speed=background_speed,
                  radial_halfwidth=float(wave.radial_halfwidth),
                  axial_halfwidth=float(wave.axial_halfwidth),
                  pulse_time_halfwidth=float(wave.time_halfwidth),
                  radial_diffusion_time=float(wave.radial_halfwidth**2
                                              / base.nu),
                  pulse_to_radial_diffusion_time=float(
                      wave.time_halfwidth * base.nu
                      / wave.radial_halfwidth**2),
                  local_k_matrix=k_matrix.tolist(), modes=results,
                  scope="Single spatial point, two nonzero angular modes, frozen physical background and phase normals, cubic-interpolated source, physical-time principal ODE and pressure identity. Not the paper's moving phase/support inverse, exact curl reconstruction, full residual reduction, or dyadic recursion.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "pulse_amplitude_inverse.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
