"""Physical local stress-cone overlap with the adaptive bridge defect."""

import json

import numpy as np

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge, sample_points
from joined_field import JoinedField
from radial_continuation import ROOT
from heat_exterior import physical
from radial_peak_cone import operator, stress_primitive
from wide_modes import WideJointModes


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    heat_amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                           / physical(point, tau0, c=1.0)["velocity"][0, 1])
    joined = JoinedField(inner=inner, join_X=1.0/64.0,
                         heat_amplitude=heat_amplitude, outer_ratio=16.0)
    base = CachedAdaptiveBridge(joined)
    coefficients = json.loads((ROOT/"adaptive_join_multiscale_fit.json").read_text())["amplitudes"]
    field = WideJointModes(base, coefficients)
    scales = []
    for k in (11.0, 19.0):
        points, tau = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
        velocity, gradient, residual = operator(field, points, tau)
        magnitudes = np.linalg.norm(residual, axis=1)
        rows = []
        for i, point in enumerate(points):
            r, _, z = point
            F = velocity[i, 1]/r
            J = gradient[i]
            shear = np.array([J[1, 0]-F, J[2, 0]])
            shear_norm = np.linalg.norm(shear)
            N = shear/shear_norm
            K = np.array([-N[1], N[0]])
            lam2 = -2.0*F*N[0]*(2.0*F*N[0]+shear_norm)
            row = dict(index=i, eta=(-0.2, 0.0, 0.2)[i//5],
                       radius=float(r), z=float(z),
                       momentum_norm=float(magnitudes[i]),
                       lambda_squared=float(lam2),
                       cone_pass=False)
            if lam2 > 0 and abs(2.0*F*N[0]) > 1e-14:
                target = stress_primitive(field, float(r), float(z), tau, order=12)
                dot_n, dot_k = float(target@N), float(target@K)
                row.update(target=target.tolist(), target_dot_N=dot_n,
                           target_dot_K=dot_k)
                if abs(dot_n) > 1e-14:
                    ratio = float(abs(np.sqrt(lam2)/(2.0*F*N[0])*dot_k/dot_n))
                    row["cone_ratio"] = ratio
                    row["cone_pass"] = bool(dot_n < 0 and ratio < 1.0)
            rows.append(row)
        sorted_rows = sorted(rows, key=lambda x:x["momentum_norm"], reverse=True)
        scales.append(dict(
            k=k, tau=tau, rows=rows,
            positive_lambda_count=sum(r["lambda_squared"]>0 for r in rows),
            cone_pass_count=sum(r["cone_pass"] for r in rows),
            top_five_positive_lambda_count=sum(r["lambda_squared"]>0 for r in sorted_rows[:5]),
            top_five_cone_pass_count=sum(r["cone_pass"] for r in sorted_rows[:5]),
            peak=sorted_rows[0],
        ))
    report = dict(
        source="Shared 24-mode adaptive bridge, physical local Section-7 cone analogue",
        scales=scales,
        scope="15 bridge nodes at each of two times. Lambda condition evaluated everywhere; radial stress primitive and cone ratio only where lambda squared is positive. This physical full-residual analogue is not the paper's normalized leading cone, and grid passes are not an open-patch certificate.",
        accepted=False, pde_validated=False,
    )
    output = ROOT/"adaptive_bridge_cone_overlap.json"
    output.write_bytes((json.dumps(report, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), counts=[
        {key: scale[key] for key in ("k", "positive_lambda_count", "cone_pass_count",
                                       "top_five_positive_lambda_count", "top_five_cone_pass_count")}
        for scale in scales]), indent=2))
    return report


if __name__ == "__main__":
    run()
