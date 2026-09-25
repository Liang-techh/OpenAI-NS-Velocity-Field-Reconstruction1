"""Axial cone extent near the broader midplane radial candidate."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from radial_peak_cone import operator, stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


ETA_NODES = np.arange(-0.20, 0.2001, 0.025)
CENTER_Y = 0.325


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    field = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    slices = []
    for k in (11.0, 19.0):
        tau = 0.5 * 2.0**-k
        X = np.full(len(ETA_NODES),
                    inner.p.X_max * (1.0 + 15.0 * CENTER_Y)**2)
        points = inner.from_similarity(X, ETA_NODES, tau)
        velocity, gradient, residual = operator(field, points, tau)
        rows = []
        for i, (r, _, z) in enumerate(points):
            F = velocity[i, 1] / r
            shear = np.array([gradient[i, 1, 0] - F,
                              gradient[i, 2, 0]])
            shear_norm = float(np.linalg.norm(shear))
            N = shear / shear_norm
            K = np.array([-N[1], N[0]])
            lam2 = float(-2 * F * N[0] * (2 * F * N[0] + shear_norm))
            row = dict(eta=float(ETA_NODES[i]), radius=float(r),
                       z=float(z),
                       momentum_norm=float(np.linalg.norm(residual[i])),
                       lambda_squared=lam2, cone_pass=False)
            if lam2 > 0 and abs(2 * F * N[0]) > 1e-14:
                target = stress_primitive(field, float(r), float(z), tau,
                                          order=12)
                dot_n, dot_k = float(target @ N), float(target @ K)
                row.update(target_dot_N=dot_n, target_dot_K=dot_k)
                if abs(dot_n) > 1e-14:
                    ratio = abs(np.sqrt(lam2) * dot_k
                                / (2 * F * N[0] * dot_n))
                    row.update(cone_ratio=float(ratio),
                               cone_pass=bool(dot_n < 0 and ratio < 1))
            rows.append(row)
        center = int(np.argmin(np.abs(ETA_NODES)))
        if rows[center]["cone_pass"]:
            left = right = center
            while left > 0 and rows[left - 1]["cone_pass"]:
                left -= 1
            while right + 1 < len(rows) and rows[right + 1]["cone_pass"]:
                right += 1
            sampled_halfwidth = min(rows[center]["z"] - rows[left]["z"],
                                    rows[right]["z"] - rows[center]["z"])
            failure_halfwidth = min(
                rows[center]["z"] - rows[left - 1]["z"] if left > 0
                else np.inf,
                rows[right + 1]["z"] - rows[center]["z"]
                if right + 1 < len(rows) else np.inf)
            passing_eta = [row["eta"] for row in rows[left:right + 1]]
        else:
            sampled_halfwidth = failure_halfwidth = 0.0
            passing_eta = []
        diffusion_equal = np.sqrt(field.nu * 0.04 * tau)
        slices.append(dict(k=k, tau=tau, rows=rows,
                           passing_eta=passing_eta,
                           sampled_symmetric_axial_halfwidth=float(
                               sampled_halfwidth),
                           nearest_failure_axial_halfwidth=float(
                               failure_halfwidth),
                           diffusion_equal_halfwidth=float(diffusion_equal),
                           required_to_sampled_width=(
                               float(diffusion_equal / sampled_halfwidth)
                               if sampled_halfwidth else None),
                           required_to_failure_width=(
                               float(diffusion_equal / failure_halfwidth)
                               if failure_halfwidth else None)))
    report = dict(source="Midplane axial cone extent for moment-closed bridge",
                  center_y=CENTER_Y, eta_nodes=ETA_NODES.tolist(),
                  slices=slices,
                  scope="17 axial similarity nodes at fixed radial similarity y=0.325 on two scales. Physical cone analogue; sampled passing nodes do not certify a 3D spacetime support or the paper's normalized cone.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "moment_bridge_midplane_axial_support.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        {key: row[key] for key in (
            "k", "passing_eta", "sampled_symmetric_axial_halfwidth",
            "nearest_failure_axial_halfwidth", "diffusion_equal_halfwidth",
            "required_to_sampled_width", "required_to_failure_width")}
        for row in slices]), indent=2))
    return report


if __name__ == "__main__":
    run()
