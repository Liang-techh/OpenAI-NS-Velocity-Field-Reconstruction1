"""Test whether a wider axial cutoff lowers direct frozen-wave momentum.

Uses the three-knot moment-repaired mean and recomputes the center physical
stress target and exact-curl covariance weights. This remains a midpoint
screen of a frozen wave, not a supported pulse inverse.
"""

import json

import numpy as np

from adaptive_bridge_curl_wave_trial import ScaledWave
from adaptive_bridge_recursive_defect import build_fields
from curl_wave_prototype import (LocalizedCurlWave, WavePerturbedField,
                                 cylindrical_residual)
from joint_collar_fit import kinematics
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run():
    inner, fields = build_fields()
    repair = json.loads((ROOT / "midplane_axial_cone_all_knots_repair.json").read_text())
    mean = SeparatedMomentModes(
        fields["two_sided_cone"], repair["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    pairs = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    angles = np.arange(16) * 2 * np.pi / 16
    scales = []
    for k in (11, 19):
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{k}.json").read_text())
        source["nu"] = mean.nu
        tau = float(source["tau"])
        point = np.asarray(source["point"], float)
        source["velocity"] = mean.fields(point[None, :], tau)[0][0].tolist()
        target = stress_primitive(mean, float(point[0]), float(point[2]), tau)
        source["local_tangential_stress_primitive"] = target.tolist()
        selected = next(row["selected"]["indices"] for row in pairs["scales"]
                        if row["k"] == k)
        old_support = support(inner, source, mean.nu)
        azimuthal_points = np.array([
            (point[0] * np.cos(a), point[0] * np.sin(a), point[2])
            for a in angles])
        axial_offsets = (-.7, -.5, -.5 / 1.4, 0., .5 / 1.4, .5, .7)
        spatial_nodes = [(x, z) for x in (-.5, 0., .5)
                         for z in axial_offsets]
        points = np.array([
            (r * np.cos(a), r * np.sin(a), z)
            for x, ez in spatial_nodes
            for r, z in [(point[0] + x * old_support["radial_halfwidth"],
                          point[2] + ez * 1.4 * old_support["axial_halfwidth"])]
            for a in angles])
        h = .0005 * np.sqrt(mean.nu * tau)
        cases = []
        for width_factor in (1., 1.4):
            args = dict(old_support)
            args["axial_halfwidth"] *= width_factor
            args["time_halfwidth"] *= width_factor**2
            wave = LocalizedCurlWave(source, pulse_indices=selected, **args)
            columns = []
            for j in range(2):
                wave.weights = np.eye(2)[j]
                values, _ = wave.fields(azimuthal_points, tau)
                cylindrical = cylindrical_residual(values, azimuthal_points)
                columns.append(np.mean(cylindrical[:, 0, None]
                                       * cylindrical[:, 1:], axis=0))
            matrix = np.column_stack(columns)
            weights = np.linalg.solve(matrix, target)
            if np.any(weights <= 0):
                cases.append(dict(width_factor=width_factor,
                                  positive_weights=False,
                                  weights=weights.tolist()))
                continue
            wave.weights = weights
            field = WavePerturbedField(mean, ScaledWave(wave, .1))
            time_step = min(.0001 * tau, args["time_halfwidth"] / 8)
            velocity, gradient, linear = kinematics(
                field, points, tau, h, time_step,
                time_min=.5 * 2.**-20)
            residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
            norms = np.linalg.norm(residual, axis=1)
            peak = int(np.argmax(norms))
            cases.append(dict(width_factor=width_factor,
                              positive_weights=True,
                              weights=weights.tolist(),
                              axial_halfwidth=args["axial_halfwidth"],
                              time_halfwidth=args["time_halfwidth"],
                              momentum_max=float(np.max(norms)),
                              momentum_rms=float(np.sqrt(np.mean(norms**2))),
                              peak_spatial_node=int(peak // len(angles)),
                              max_fd_divergence=float(np.max(np.abs(
                                  np.trace(gradient, axis1=1, axis2=2))))))
        scales.append(dict(k=k, target_stress=target.tolist(),
                           spatial_nodes=spatial_nodes, cases=cases))
        print(f"k={k}: {[c.get('momentum_max') for c in cases]}", flush=True)
    report = dict(source="Wider axial support frozen-wave screen",
                  mean="all-three-knot moment-repaired axial-cone mean",
                  scales=scales,
                  scope="Midpoint full Cartesian momentum on a common 21-node "
                        "radial-axial grid and 16 angles. Cone support between "
                        "nodes and across time is unproved; no pulse inverse, "
                        "pressure, mean correction, or recursive acceptance.",
                  accepted=False, scale_recursion_established=False)
    (ROOT / "midplane_wider_cone_wave_screen.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
