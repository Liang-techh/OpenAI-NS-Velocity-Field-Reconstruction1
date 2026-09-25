"""Refine the finite-difference stencil at the exact-curl wave hotspot."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_curl_wave_trial import ScaledWave
from curl_wave_prototype import LocalizedCurlWave, WavePerturbedField
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


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
    ri = float(inner.from_similarity([inner.p.X_max], [-0.2], tau)[0, 0])
    wave = LocalizedCurlWave(
        source, radial_halfwidth=0.0047 * 15.0 * ri,
        axial_halfwidth=1.0e-4, time_halfwidth=0.04 * tau,
        min_tau=0.5 * 2.0**-20)
    point = inner.from_similarity(
        [inner.p.X_max * (1.0 + 15.0 * 0.0525)**2], [-0.2], tau)
    field = WavePerturbedField(base, ScaledWave(wave, 0.1))
    h0, ht = 0.0005 * np.sqrt(base.nu * tau), 0.0001 * tau
    rows = []
    for factor in (1.0, 0.5, 0.25, 0.125):
        velocity, gradient, linear, terms = kinematics(
            field, point, tau, factor * h0, ht,
            time_min=0.5 * 2.0**-20, return_terms=True)
        residual = linear + np.einsum("nij,nj->ni", gradient, velocity)
        rows.append(dict(
            h_factor=factor,
            momentum_norm=float(np.linalg.norm(residual[0])),
            fd_divergence=float(np.trace(gradient[0])),
            viscosity_norm=float(np.linalg.norm(terms["viscosity"][0])),
        ))
    report = dict(
        source="Stencil refinement at one wave-perturbed bridge point",
        point=point[0].tolist(), rows=rows,
        scope="One point and one amplitude scale. The wave is an exact analytic curl, while finite-difference divergence and momentum depend on stencil resolution; no full-support bound or PDE acceptance.",
        accepted=False,
    )
    output = ROOT / "adaptive_bridge_wave_fd_refinement.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
