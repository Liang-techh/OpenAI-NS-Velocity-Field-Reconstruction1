"""Direct full-momentum trial of a support-adapted midplane curl wave."""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from adaptive_bridge_recursive_defect import build_fields
from curl_wave_prototype import (LocalizedCurlWave, WavePerturbedField,
                                 cylindrical_residual)
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    base = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    pairs_report = json.loads((ROOT / "midplane_physical_covariance_pairs.json")
                              .read_text())
    scales = []
    for k in (11.0, 19.0):
        chosen = next(row["selected"] for row in pairs_report["scales"]
                      if row["k"] == int(k))
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{int(k)}.json").read_text())
        source["nu"] = base.nu
        tau = float(source["tau"])
        X = inner.p.X_max * (1 + 15 * 0.325)**2
        ri = float(inner.from_similarity([inner.p.X_max], [-0.0125], tau)[0, 0])
        radial_halfwidth = 0.035 * 15 * ri
        below = inner.from_similarity([X], [-0.0375], tau)[0, 2]
        above = inner.from_similarity([X], [0.0125], tau)[0, 2]
        zcenter = float(source["point"][2])
        axial_halfwidth = 0.99 * min(zcenter - below, above - zcenter)
        time_halfwidth = 0.75 * axial_halfwidth**2 / base.nu
        wave = LocalizedCurlWave(
            source, radial_halfwidth=radial_halfwidth,
            axial_halfwidth=axial_halfwidth,
            time_halfwidth=time_halfwidth,
            pulse_indices=chosen["indices"],
            min_tau=0.5 * 2.0**-20)
        angles = np.arange(16) * 2 * np.pi / 16
        center_points = np.array([
            (wave.radius * np.cos(angle),
             wave.radius * np.sin(angle), wave.zcenter)
            for angle in angles])

        def sampled_covariance():
            values, _ = wave.fields(center_points, tau)
            cyl = cylindrical_residual(values, center_points)
            return np.mean(cyl[:, 0, None] * cyl[:, 1:], axis=0)

        target = wave.target
        ideal_weights = wave.weights.copy()
        ideal_covariance = sampled_covariance()
        physical_weights = np.asarray(chosen["positive_weights"], float)
        wave.weights = physical_weights
        covariance = sampled_covariance()
        spatial_nodes = ((0.0, 0.0), (-0.5, 0.0), (0.5, 0.0),
                         (0.0, -0.5), (0.0, 0.5))
        points = np.array([
            (r * np.cos(angle), r * np.sin(angle), z)
            for xi, ez in spatial_nodes
            for r, z in [(wave.radius + xi * radial_halfwidth,
                          wave.zcenter + ez * axial_halfwidth)]
            for angle in angles])
        hs = 0.0005 * np.sqrt(base.nu * tau)
        ht = min(0.0001 * tau, time_halfwidth / 8)
        cases = []
        for amplitude in (0.0, 0.1, 1.0):
            field = WavePerturbedField(base, ScaledWave(wave, amplitude))
            velocity, gradient, linear, terms = kinematics(
                field, points, tau, hs, ht, time_min=0.5 * 2.0**-20,
                return_terms=True)
            residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
            norms = np.linalg.norm(residual, axis=1)
            peak = int(np.argmax(norms))
            cases.append(dict(amplitude=amplitude,
                              peak_point=points[peak].tolist(),
                              peak_terms={name: value[peak].tolist()
                                          for name, value in terms.items()},
                              peak_advection=np.einsum(
                                  "ij,j->i", gradient[peak],
                                  velocity[peak]).tolist(),
                              momentum_max=float(np.max(norms)),
                              momentum_rms=float(np.sqrt(np.mean(norms**2))),
                              fd_divergence_max=float(np.max(np.abs(np.trace(
                                  gradient, axis1=1, axis2=2))))))
        wave_velocity, _ = wave.fields(center_points, tau)
        scales.append(dict(
            k=k, tau=tau, source_point=source["point"],
            source_cone=source["cone"],
            radial_halfwidth=float(radial_halfwidth),
            axial_halfwidth=float(axial_halfwidth),
            time_halfwidth=float(time_halfwidth),
            axial_pulse_to_diffusion=float(
                time_halfwidth * base.nu / axial_halfwidth**2),
            radial_pulse_to_diffusion=float(
                time_halfwidth * base.nu / radial_halfwidth**2),
            angular_modes=[mode["m"] for mode in wave.waves],
            selected_pulse_indices=chosen["indices"],
            ideal_weights=ideal_weights.tolist(),
            physical_covariance_weights=physical_weights.tolist(),
            ideal_sampled_covariance=ideal_covariance.tolist(),
            ideal_covariance_relative_error=float(np.linalg.norm(
                ideal_covariance - target) / np.linalg.norm(target)),
            sampled_covariance=covariance.tolist(),
            target_stress=target.tolist(),
            covariance_relative_error=float(np.linalg.norm(
                covariance - target) / np.linalg.norm(target)),
            wave_speed_center_max=float(np.max(np.linalg.norm(
                wave_velocity, axis=1))),
            cases=cases))
    report = dict(source="Two-scale support-adapted midplane exact-curl wave trial",
                  spatial_nodes=[list(pair) for pair in spatial_nodes],
                  angles_per_node=len(angles), scales=scales,
                  scope="Five radial/axial patch nodes and sixteen angles at the source time on k=11,19. Exact-curl wave plus full Cartesian FD momentum, with frozen Kelvin amplitudes and no pressure/mean wave correction. No interior-time or full-support residual bound, moment closure after wave, or scale-recursive acceptance.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "midplane_wave_two_scale_trial.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        {key: item[key] for key in (
            "k", "radial_halfwidth", "axial_halfwidth", "time_halfwidth",
            "covariance_relative_error", "wave_speed_center_max", "cases")}
        for item in scales]), indent=2))
    return report


if __name__ == "__main__":
    run()
