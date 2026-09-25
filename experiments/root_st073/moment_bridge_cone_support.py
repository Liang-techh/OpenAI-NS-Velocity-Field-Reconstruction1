"""Map usable radial cone support of the three-knot moment-closed bridge.

This tests a sampled physical cone analogue at two dyadic scales. A passing
node sequence is not a continuous cone or an admissible paper profile.
"""

import argparse
import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from radial_peak_cone import operator, stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


Y_NODES = np.arange(0.035, 0.0651, 0.0025)
CENTER_Y = 0.05


def run(y_nodes=Y_NODES, center_y=CENTER_Y, eta=-0.2,
        output_name="moment_bridge_cone_support.json"):
    y_nodes = np.asarray(y_nodes, float)
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    field = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    slices = []
    for k in (11.0, 19.0):
        tau = 0.5 * 2.0**-k
        X = inner.p.X_max * (1.0 + 15.0 * y_nodes)**2
        points = inner.from_similarity(X, np.full(len(X), eta), tau)
        velocity, gradient, residual = operator(field, points, tau)
        rows = []
        for i, (r, _, z) in enumerate(points):
            F = velocity[i, 1] / r
            shear = np.array([gradient[i, 1, 0] - F,
                              gradient[i, 2, 0]])
            magnitude = float(np.linalg.norm(shear))
            N = shear / magnitude
            K = np.array([-N[1], N[0]])
            lam2 = float(-2.0 * F * N[0] * (2.0 * F * N[0] + magnitude))
            row = dict(y=float(y_nodes[i]), radius=float(r),
                       momentum_norm=float(np.linalg.norm(residual[i])),
                       lambda_squared=lam2, cone_pass=False)
            if lam2 > 0:
                target = stress_primitive(field, float(r), float(z), tau,
                                          order=12)
                dot_n, dot_k = float(target @ N), float(target @ K)
                row.update(target_dot_N=dot_n, target_dot_K=dot_k)
                if abs(dot_n) > 1e-14 and abs(2.0 * F * N[0]) > 1e-14:
                    ratio = abs(np.sqrt(lam2) * dot_k
                                / (2.0 * F * N[0] * dot_n))
                    row.update(cone_ratio=float(ratio),
                               cone_pass=bool(dot_n < 0 and ratio < 1))
            rows.append(row)
        center = int(np.argmin(np.abs(y_nodes - center_y)))
        if not rows[center]["cone_pass"]:
            halfwidth = 0.0
            failure_bound = 0.0
            band_halfwidth = 0.0
            failure_band_bound = 0.0
            passing = []
        else:
            left = right = center
            while left > 0 and rows[left - 1]["cone_pass"]:
                left -= 1
            while right + 1 < len(rows) and rows[right + 1]["cone_pass"]:
                right += 1
            passing = rows[left:right + 1]
            band_halfwidth = 0.5 * (rows[right]["radius"]
                                    - rows[left]["radius"])
            failure_band_bound = (0.5 * (rows[right + 1]["radius"]
                                         - rows[left - 1]["radius"])
                                  if left > 0 and right + 1 < len(rows)
                                  else None)
            halfwidth = min(rows[center]["radius"] - rows[left]["radius"],
                            rows[right]["radius"] - rows[center]["radius"])
            failure_distances = []
            if left > 0:
                failure_distances.append(rows[center]["radius"]
                                         - rows[left - 1]["radius"])
            if right + 1 < len(rows):
                failure_distances.append(rows[right + 1]["radius"]
                                         - rows[center]["radius"])
            failure_bound = (min(failure_distances) if failure_distances
                             else None)
        ri = float(inner.from_similarity([inner.p.X_max], [eta], tau)[0, 0])
        current_halfwidth = 0.0047 * 15.0 * ri
        pulse_halfwidth = 0.04 * tau
        diffusion_equal_halfwidth = np.sqrt(field.nu * pulse_halfwidth)
        slices.append(dict(k=k, tau=tau, rows=rows,
                           passing_y=[row["y"] for row in passing],
                           sampled_symmetric_cone_halfwidth=float(halfwidth),
                           best_center_sampled_band_halfwidth=float(
                               band_halfwidth),
                           failure_bracket_band_halfwidth=(
                               float(failure_band_bound)
                               if failure_band_bound is not None else None),
                           nearest_sampled_failure_halfwidth=(
                               float(failure_bound)
                               if failure_bound is not None else None),
                           current_wave_radial_halfwidth=current_halfwidth,
                           diffusion_equal_halfwidth=float(
                               diffusion_equal_halfwidth),
                           pulse_to_diffusion_at_sampled_halfwidth=(
                               float(field.nu * pulse_halfwidth / halfwidth**2)
                               if halfwidth > 0 else None),
                           required_width_to_sampled_width=(
                               float(diffusion_equal_halfwidth / halfwidth)
                               if halfwidth > 0 else None),
                           required_width_to_failure_bound=(
                               float(diffusion_equal_halfwidth / failure_bound)
                               if failure_bound else None),
                           required_width_to_failure_band_bound=(
                               float(diffusion_equal_halfwidth
                                     / failure_band_bound)
                               if failure_band_bound else None)))
    report = dict(source="Three-knot moment-closed radial cone support at two scales",
                  y_nodes=y_nodes.tolist(), center_y=center_y,
                  eta=eta, slices=slices,
                  scope="Radial nodes at one axial similarity coordinate on each of two scales; physical full-residual cone analogue with Gauss12 stress primitive. Passing samples do not certify an interval, nor the paper's normalized cone, global wave support, momentum gate, or scale recursion.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / output_name
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        {key: row[key] for key in (
            "k", "passing_y", "sampled_symmetric_cone_halfwidth",
            "nearest_sampled_failure_halfwidth",
            "best_center_sampled_band_halfwidth",
            "failure_bracket_band_halfwidth",
            "current_wave_radial_halfwidth", "diffusion_equal_halfwidth",
            "pulse_to_diffusion_at_sampled_halfwidth",
            "required_width_to_sampled_width",
            "required_width_to_failure_bound",
            "required_width_to_failure_band_bound")}
        for row in slices]), indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--midplane", action="store_true")
    args = parser.parse_args()
    if args.midplane:
        run(y_nodes=np.arange(0.25, 0.4501, 0.01), center_y=0.35,
            eta=0.0, output_name="moment_bridge_midplane_cone_support.json")
    else:
        run()
