"""Refine the sampled cone-passing radial band near both bridge hotspots."""

import json

import numpy as np

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge
from joined_field import JoinedField
from radial_continuation import ROOT
from heat_exterior import physical
from radial_peak_cone import operator, stress_primitive
from wide_modes import WideJointModes


Y_NODES = np.array([0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
                    0.055, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20])


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    heat_amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                           / physical(point, tau0, c=1.0)["velocity"][0, 1])
    joined = JoinedField(inner=inner, join_X=1.0/64.0,
                         heat_amplitude=heat_amplitude, outer_ratio=16.0)
    report = json.loads((ROOT/"adaptive_bridge_cone_biside_fit.json").read_text())
    selected = next(row for row in report["candidates"]
                    if row["cone_weight"] == report["selected_cone_weight"])
    field = WideJointModes(CachedAdaptiveBridge(joined), selected["amplitudes"])
    slices = []
    for k in (11.0, 19.0):
        tau = 0.5*2.0**(-k)
        for eta in (-0.2, 0.2):
            X = inner.p.X_max*(1.0+15.0*Y_NODES)**2
            points = inner.from_similarity(X, np.full(len(X), eta), tau)
            velocity, gradient, residual = operator(field, points, tau)
            rows = []
            for i, point in enumerate(points):
                radius, _, z = point
                F = velocity[i, 1]/radius
                J = gradient[i]
                shear = np.array([J[1, 0]-F, J[2, 0]])
                norm = np.linalg.norm(shear)
                N = shear/norm
                K = np.array([-N[1], N[0]])
                lam2 = float(-2.0*F*N[0]*(2.0*F*N[0]+norm))
                row = dict(y=float(Y_NODES[i]), radius=float(radius),
                           momentum_norm=float(np.linalg.norm(residual[i])),
                           lambda_squared=lam2, cone_pass=False)
                if lam2 > 0:
                    target = stress_primitive(field, float(radius), float(z), tau, order=12)
                    dot_n, dot_k = float(target@N), float(target@K)
                    row.update(target_dot_N=dot_n, target_dot_K=dot_k)
                    if abs(dot_n)>1e-14 and abs(2.0*F*N[0])>1e-14:
                        ratio = float(abs(np.sqrt(lam2)/(2.0*F*N[0])*dot_k/dot_n))
                        row.update(cone_ratio=ratio,
                                   cone_pass=bool(dot_n<0 and ratio<1))
                rows.append(row)
            passing = [r for r in rows if r["cone_pass"]]
            sampled_span = (max(r["radius"] for r in passing)-min(r["radius"] for r in passing)
                            if len(passing) >= 2 else None)
            slices.append(dict(k=k, eta=eta, rows=rows,
                               passing_count=len(passing),
                               passing_y=[r["y"] for r in passing],
                               sampled_passing_radius_width=sampled_span,
                               sampled_span_diffusion_time=(sampled_span**2/inner.nu
                                                            if sampled_span is not None else None),
                               tenth_tau_to_diffusion_time=(0.1*tau/(sampled_span**2/inner.nu)
                                                            if sampled_span is not None else None),
                               sampled_passing_radius_span=(
                                   [min(r["radius"] for r in passing),
                                    max(r["radius"] for r in passing)]
                                   if passing else None)))
    result = dict(
        source="Selected two-sided grid-aware mean correction",
        y_nodes=Y_NODES.tolist(), slices=slices,
        scope="Fourteen radial nodes at each of eta=±.2 and k=11,19. A passing set of nodes is not a certified continuous support interval; no wave amplitude/phase solve or full momentum acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT/"adaptive_bridge_cone_patch.json"
    output.write_bytes((json.dumps(result, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        {key:s[key] for key in ("k", "eta", "passing_count", "passing_y")}
        for s in slices]), indent=2))
    return result


if __name__ == "__main__":
    run()
