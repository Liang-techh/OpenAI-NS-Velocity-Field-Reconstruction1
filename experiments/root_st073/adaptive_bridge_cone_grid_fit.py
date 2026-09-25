"""Trade full-grid momentum against cone geometry at bridge hotspots."""

import json

import numpy as np
from scipy.optimize import least_squares

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge, sample_points
from adaptive_join_scale_knots import assemble_slice
from affine_momentum import momentum, combine
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from heat_exterior import physical
from radial_peak_cone import operator, stress_primitive
from wide_modes import WideJointModes


def geometry(jet, radius):
    u, J, _ = jet
    F = u[0, 1]/radius
    shear = np.array([J[0, 1, 0]-F, J[0, 2, 0]])
    norm2 = max(float(np.dot(shear, shear)), 1e-30)
    g = float(2.0*F*shear[0]/norm2)
    return dict(F=float(F), shear=shear.tolist(), g=g,
                normalized_lambda_squared=float(-g*(g+1.0)))


def run(hotspot_indices=(0,), output_name="adaptive_bridge_cone_grid_fit.json"):
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    heat_amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                           / physical(point, tau0, c=1.0)["velocity"][0, 1])
    joined = JoinedField(inner=inner, join_X=1.0/64.0,
                         heat_amplitude=heat_amplitude, outer_ratio=16.0)
    cached = CachedAdaptiveBridge(joined)
    zero = WideJointModes(cached, np.zeros(24))
    slices = [assemble_slice(cached, zero, inner, 16.0, k) for k in (11.0, 19.0)]
    source = np.asarray(json.loads((ROOT/"adaptive_join_multiscale_fit.json").read_text())["amplitudes"])
    hotspot_sites = [(si, i, float(s["args"][0][i, 0]))
                     for si, s in enumerate(slices) for i in hotspot_indices]

    def features(a):
        evaluated = [combine(s["baseline"], s["modes"], a) for s in slices]
        return [geometry(tuple(part[i:i+1] for part in evaluated[si]), r)
                for si, i, r in hotspot_sites]

    candidates = []
    for cone_weight in (0.2, 1.0, 5.0):
        def objective(a):
            rows = features(a)
            residuals = [(momentum(combine(s["baseline"], s["modes"], a))/s["scale"]).ravel()
                         for s in slices]
            cone = np.array([cone_weight*(row["g"]+0.5) for row in rows])
            return np.concatenate([*residuals, cone, 0.002*(a-source)])

        fit = least_squares(objective, source, bounds=(-20.0, 20.0),
                            max_nfev=100, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        rows = features(fit.x)
        peaks = [float(np.max(np.linalg.norm(
            momentum(combine(s["baseline"], s["modes"], fit.x)), axis=1)))
            for s in slices]
        candidates.append(dict(
            cone_weight=cone_weight, amplitudes=fit.x.tolist(),
            fit_nfev=fit.nfev, fit_success=bool(fit.success),
            geometry=rows, predicted_grid_peaks=peaks,
            both_lambda_positive=all(row["normalized_lambda_squared"]>0 for row in rows),
            max_peak_ratio_to_shared=float(max(
                new/np.max(np.linalg.norm(momentum(combine(s["baseline"], s["modes"], source)), axis=1))
                for new, s in zip(peaks, slices))),
        ))
    feasible = [c for c in candidates if c["both_lambda_positive"]]
    selected = min(feasible, key=lambda c:c["max_peak_ratio_to_shared"]) if feasible else min(
        candidates, key=lambda c:sum(max(0.0, -r["normalized_lambda_squared"]) for r in c["geometry"]))
    field = WideJointModes(cached, selected["amplitudes"])
    training = []
    for s in slices:
        points, tau, h, ht = s["args"]
        direct, divergence = independent_fd(field, points, tau, h, ht)
        training.append(dict(k=s["k"],
                             shared_max=float(np.max(np.linalg.norm(
                                 momentum(combine(s["baseline"], s["modes"], source)), axis=1))),
                             selected_max=float(np.max(np.linalg.norm(direct, axis=1))),
                             divergence_max=float(np.max(np.abs(divergence)))))
    holdouts = []
    for k in (10.5, 11.5, 18.5, 19.5):
        points, tau = sample_points(inner, 16.0, k, (-0.1, 0.1), 6)
        h, ht = 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau
        old, _ = independent_fd(WideJointModes(cached, source), points, tau, h, ht)
        new, divergence = independent_fd(field, points, tau, h, ht)
        holdouts.append(dict(k=k,
                             shared_max=float(np.max(np.linalg.norm(old, axis=1))),
                             selected_max=float(np.max(np.linalg.norm(new, axis=1))),
                             divergence_max=float(np.max(np.abs(divergence)))))
    cone_direct = []
    for si, i, _ in hotspot_sites:
        s = slices[si]
        point = s["args"][0][i:i+1]
        radius, z = float(point[0, 0]), float(point[0, 2])
        tau = s["args"][1]
        u, J, residual = operator(field, point, tau)
        F = u[0, 1]/radius
        shear = np.array([J[0, 1, 0]-F, J[0, 2, 0]])
        norm = np.linalg.norm(shear)
        N = shear/norm
        K = np.array([-N[1], N[0]])
        lam2 = float(-2.0*F*N[0]*(2.0*F*N[0]+norm))
        row = dict(k=s["k"], hotspot_index=i, lambda_squared=lam2,
                   momentum_norm=float(np.linalg.norm(residual[0])), cone_pass=False)
        if lam2>0:
            target = stress_primitive(field, radius, z, tau, order=12)
            dot_n, dot_k = float(target@N), float(target@K)
            row.update(target_dot_N=dot_n, target_dot_K=dot_k)
            if abs(dot_n)>1e-14 and abs(2.0*F*N[0])>1e-14:
                ratio = float(abs(np.sqrt(lam2)/(2.0*F*N[0])*dot_k/dot_n))
                row.update(cone_ratio=ratio, cone_pass=bool(dot_n<0 and ratio<1))
        cone_direct.append(row)
    report = dict(
        source="Full-grid momentum with selected physical cone-geometry penalties",
        hotspot_indices=list(hotspot_indices),
        candidates=candidates, selected_cone_weight=selected["cone_weight"],
        training=training, holdouts=holdouts, cone_direct=cone_direct,
        scope="15 training bridge nodes per scale at k=11 and 19, four off-scale holdouts; local physical cone analogue only. No radial moment/stress target coupling in optimization, volume L2, finite energy, or PDE acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT/output_name
    output.write_bytes((json.dumps(report, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), candidates=[
        {key:c[key] for key in ("cone_weight", "both_lambda_positive", "max_peak_ratio_to_shared")}
        for c in candidates], training=training, holdouts=holdouts,
        cone_direct=cone_direct), indent=2))
    return report


if __name__ == "__main__":
    run()
