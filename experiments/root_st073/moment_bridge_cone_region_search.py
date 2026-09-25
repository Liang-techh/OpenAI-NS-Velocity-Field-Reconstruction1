"""Coarse radial-axial cone search for the moment-closed adaptive bridge.

The search looks for a replacement wave location at two separated scales.
Passing points only nominate regions for refinement; they do not certify
support or a normalized leading-profile stress cone.
"""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from radial_peak_cone import operator, stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


Y_NODES = (0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95)
ETA_NODES = (-0.30, -0.20, 0.0, 0.20, 0.30)


def run():
    inner, fields = build_fields()
    fitted = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    field = SeparatedMomentModes(
        fields["two_sided_cone"], fitted["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    scales = []
    for k in (11.0, 19.0):
        tau = 0.5 * 2.0**-k
        coordinates = [(y, eta) for eta in ETA_NODES for y in Y_NODES]
        y = np.array([pair[0] for pair in coordinates])
        eta = np.array([pair[1] for pair in coordinates])
        X = inner.p.X_max * (1.0 + 15.0 * y)**2
        points = inner.from_similarity(X, eta, tau)
        velocity, gradient, residual = operator(field, points, tau)
        rows = []
        for i, (r, _, z) in enumerate(points):
            F = velocity[i, 1] / r
            shear = np.array([gradient[i, 1, 0] - F,
                              gradient[i, 2, 0]])
            shear_norm = float(np.linalg.norm(shear))
            N = shear / shear_norm
            K = np.array([-N[1], N[0]])
            lam2 = float(-2.0 * F * N[0] * (2.0 * F * N[0]
                                                 + shear_norm))
            row = dict(y=float(y[i]), eta=float(eta[i]),
                       radius=float(r), z=float(z),
                       momentum_norm=float(np.linalg.norm(residual[i])),
                       lambda_squared=lam2, cone_pass=False)
            if lam2 > 0 and abs(2.0 * F * N[0]) > 1e-14:
                target = stress_primitive(field, float(r), float(z), tau,
                                          order=12)
                dot_n, dot_k = float(target @ N), float(target @ K)
                row.update(target_dot_N=dot_n, target_dot_K=dot_k)
                if abs(dot_n) > 1e-14:
                    ratio = abs(np.sqrt(lam2) * dot_k
                                / (2.0 * F * N[0] * dot_n))
                    row.update(cone_ratio=float(ratio),
                               cone_pass=bool(dot_n < 0 and ratio < 1))
            rows.append(row)
        passing = [row for row in rows if row["cone_pass"]]
        scales.append(dict(k=k, tau=tau, rows=rows,
                           lambda_positive_count=sum(
                               row["lambda_squared"] > 0 for row in rows),
                           cone_pass_count=len(passing),
                           passing_coordinates=[
                               [row["y"], row["eta"]] for row in passing]))
    overlap = sorted(set(map(tuple, scales[0]["passing_coordinates"]))
                     & set(map(tuple, scales[1]["passing_coordinates"])))
    report = dict(source="Coarse alternate-cone search on three-knot moment-closed bridge",
                  y_nodes=list(Y_NODES), eta_nodes=list(ETA_NODES),
                  scales=scales,
                  passing_coordinate_overlap=[list(pair) for pair in overlap],
                  scope="45 radial-axial nodes per scale at k=11,19. Full physical residual cone analogue and Gauss12 stress primitive where lambda squared is positive. This map nominates points only; no continuous cone, dyadic recursion, supported wave, or PDE acceptance.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "moment_bridge_cone_region_search.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          counts=[(x["k"], x["lambda_positive_count"],
                                   x["cone_pass_count"]) for x in scales],
                          passing_coordinate_overlap=report[
                              "passing_coordinate_overlap"]), indent=2))
    return report


if __name__ == "__main__":
    run()
