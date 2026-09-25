"""Estimate phase-normal transport neglected by the frozen Kelvin wave."""

import json

import numpy as np

from joint_collar_fit import kinematics
from midplane_wider_cone_slope_projection import build_case
from radial_continuation import ROOT


def run():
    angles = np.arange(16) * 2 * np.pi / 16
    scales = []
    for k in (11, 19):
        mean, wave, args, _ = build_case(k, angles)
        point = np.array([[wave.radius, 0., wave.zcenter]])
        tau = wave.tau0
        hs = .0005 * np.sqrt(mean.nu * tau)
        ht = min(.0001 * tau, args["time_halfwidth"] / 8)
        _, gradient, _ = kinematics(mean, point, tau, hs, ht,
                                    time_min=.5 * 2.**-20)
        J = gradient[0]
        interval = .1 * args["time_halfwidth"]
        modes = []
        for mode in wave.waves:
            normal = np.asarray(mode["normal"], float)
            # For physical time t, Kelvin transport is n_t=-J^T n;
            # tau is time-to-critical, so n_tau=+J^T n at this point.
            n_tau = J.T @ normal
            transverse = n_tau - normal * (np.dot(normal, n_tau)
                                            / np.dot(normal, normal))
            modes.append(dict(angular_mode=int(mode["m"]),
                              normal_norm=float(np.linalg.norm(normal)),
                              relative_first_order_normal_change=float(
                                  interval * np.linalg.norm(n_tau)
                                  / np.linalg.norm(normal)),
                              angular_first_order_change_radians=float(
                                  interval * np.linalg.norm(transverse)
                                  / np.linalg.norm(normal))))
        scales.append(dict(k=k, interval=interval, modes=modes))
    report = dict(source="Local frozen-normal transport-size diagnostic",
                  scales=scales,
                  scope="First-order Kelvin normal drift J^T n over one tenth "
                        "of the pulse half-width at the source point only. "
                        "It is not a transported phase solve or a full error "
                        "budget.")
    (ROOT / "midplane_wider_cone_normal_drift.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(report["scales"], indent=2))
    return report


if __name__ == "__main__":
    run()
