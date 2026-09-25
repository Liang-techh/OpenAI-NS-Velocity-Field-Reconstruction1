"""Full-momentum trial of a compact exact-curl wave on the late bridge."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_prototype import (LocalizedCurlWave, WavePerturbedField,
                                 cylindrical_residual)
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


class ScaledWave:
    def __init__(self, wave, scale):
        self.wave = wave
        self.scale = float(scale)

    def fields(self, points, tau):
        velocity, pressure = self.wave.fields(points, tau)
        return self.scale * velocity, pressure


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    source = json.loads((ROOT / "compact_potential" /
                         "adaptive_bridge_wave_source.json").read_text())
    source["nu"] = base.nu
    tau = float(source["tau"])
    inner_radius = float(inner.from_similarity([inner.p.X_max], [-0.2], tau)[0, 0])
    bridge_width = 15.0 * inner_radius
    radial_halfwidth = 0.0047 * bridge_width
    axial_halfwidth = 1.0e-4
    time_halfwidth = 0.04 * tau
    wave = LocalizedCurlWave(
        source, radial_halfwidth=radial_halfwidth,
        axial_halfwidth=axial_halfwidth,
        time_halfwidth=time_halfwidth,
        min_tau=0.5 * 2.0**-20)
    angles = np.arange(16) * 2.0 * np.pi / 16.0
    points = []
    for y in (0.0475, 0.05, 0.0525):
        X = inner.p.X_max * (1.0 + 15.0 * y)**2
        location = inner.from_similarity([X], [-0.2], tau)[0]
        radius, z = float(location[0]), float(location[2])
        points.extend((radius * np.cos(angle), radius * np.sin(angle), z)
                      for angle in angles)
    points = np.asarray(points)
    h, ht = 0.0005 * np.sqrt(base.nu * tau), 0.0001 * tau
    cases = []
    for scale in (0.0, 0.1, 0.3, 1.0):
        current = WavePerturbedField(base, ScaledWave(wave, scale))
        velocity, gradient, linear, terms = kinematics(
            current, points, tau, h, ht, time_min=0.5 * 2.0**-20,
            return_terms=True)
        residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
        peak_index = int(np.argmax(np.linalg.norm(residual, axis=1)))
        cyl = cylindrical_residual(residual, points).reshape(3, 16, 3)
        wave_velocity, _ = wave.fields(points, tau)
        cases.append(dict(
            amplitude_scale=scale,
            peak_index=peak_index,
            peak_point=points[peak_index].tolist(),
            peak_momentum=residual[peak_index].tolist(),
            peak_terms={name: value[peak_index].tolist()
                        for name, value in terms.items()},
            full_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
            full_momentum_rms=float(np.sqrt(np.mean(np.sum(residual**2, axis=1)))),
            angle_mean_cylindrical_residual=cyl.mean(axis=1).tolist(),
            fd_divergence_max=float(np.max(np.abs(np.trace(
                gradient, axis1=1, axis2=2)))),
            wave_speed_max=float(scale * np.max(np.linalg.norm(wave_velocity, axis=1))),
        ))
    center = points[16:32]
    w, _ = wave.fields(center, tau)
    wc = cylindrical_residual(w, center)
    covariance = np.mean(wc[:, 0, None] * wc[:, 1:], axis=0)
    report = dict(
        source="Three-knot moment-closed bridge plus two-harmonic compact exact-curl wave",
        k=11.0, tau=tau, eta=-0.2,
        radial_halfwidth=radial_halfwidth,
        axial_halfwidth=axial_halfwidth,
        time_halfwidth=time_halfwidth,
        selected_pulse_indices=wave.pulse_indices,
        angular_modes=[row["m"] for row in wave.waves],
        target_stress=wave.target.tolist(),
        sampled_wave_covariance=covariance.tolist(),
        covariance_relative_error=float(np.linalg.norm(covariance-wave.target)
                                        / np.linalg.norm(wave.target)),
        carrier_radial_support_products=[float(abs(row["normal"][0])
                                                * radial_halfwidth)
                                         for row in wave.waves],
        carrier_axial_support_products=[float(abs(row["normal"][2])
                                               * axial_halfwidth)
                                        for row in wave.waves],
        carrier_viscous_exponents=[float(base.nu * np.dot(row["normal"],
                                                             row["normal"])
                                         * time_halfwidth)
                                   for row in wave.waves],
        cases=cases,
        scope="Three bridge radii and sixteen angles at one axial/time point; exact analytic curl perturbation, full Cartesian finite-difference nonlinear momentum. Frozen local Kelvin amplitudes, no transported wave equation, pressure correction, global compact stress, continuum bounds, or PDE acceptance.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_bridge_curl_wave_trial.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          pulse_indices=wave.pulse_indices,
                          angular_modes=report["angular_modes"],
                          covariance_relative_error=report["covariance_relative_error"],
                          radial_products=report["carrier_radial_support_products"],
                          axial_products=report["carrier_axial_support_products"],
                          viscous_exponents=report["carrier_viscous_exponents"],
                          cases=cases), indent=2))
    return report


if __name__ == "__main__":
    run()
