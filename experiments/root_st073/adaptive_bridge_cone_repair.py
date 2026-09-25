"""Bounded joint mean-mode search for cone geometry at bridge hotspots."""

import json

import numpy as np
from scipy.optimize import least_squares

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge
from affine_momentum import jets, momentum, combine
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from heat_exterior import physical
from radial_peak_cone import operator, stress_primitive
from wide_modes import WideJointModes


PEAKS = ((11.0, 0.00048027512363772447, -0.00033242031834369965),
         (19.0, 0.00003001719522735778, -0.000021360370346910557))


def cone_geometry(jet, radius):
    u, gradient, _ = jet
    F = u[0, 1]/radius
    shear = np.array([gradient[0, 1, 0]-F, gradient[0, 2, 0]])
    norm2 = float(np.dot(shear, shear))
    g = float(2.0*F*shear[0]/norm2)
    return dict(F=float(F), shear=shear.tolist(), g=g,
                normalized_lambda_squared=float(-g*(g+1.0)),
                momentum_norm=float(np.linalg.norm(momentum(jet)[0])))


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    heat_amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                           / physical(point, tau0, c=1.0)["velocity"][0, 1])
    joined = JoinedField(inner=inner, join_X=1.0/64.0,
                         heat_amplitude=heat_amplitude, outer_ratio=16.0)
    cached = CachedAdaptiveBridge(joined)
    source = np.asarray(json.loads((ROOT/"adaptive_join_multiscale_fit.json").read_text())["amplitudes"])
    zero = WideJointModes(cached, np.zeros(24))
    sites = []
    for k, radius, z in PEAKS:
        tau = 0.5*2.0**(-k)
        points = np.array([[radius, 0.0, z]])
        args = (points, tau, 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau)
        baseline = jets(zero, *args)
        unit_jets = []
        for j in range(24):
            amplitude = np.zeros(24)
            amplitude[j] = 1.0
            unit = jets(WideJointModes(cached, amplitude), *args)
            unit_jets.append(tuple(u-b for u, b in zip(unit, baseline)))
        modes = tuple(np.stack([row[i] for row in unit_jets]) for i in range(3))
        original = cone_geometry(combine(baseline, modes, source), radius)
        sites.append(dict(k=k, radius=radius, z=z, tau=tau,
                          baseline=baseline, modes=modes,
                          original=original, original_norm=original["momentum_norm"]))

    def features(a):
        return [cone_geometry(combine(site["baseline"], site["modes"], a),
                              site["radius"]) for site in sites]

    candidates = []
    for cone_weight in (0.2, 1.0, 5.0):
        def objective(a):
            rows = features(a)
            pieces = []
            for site, row in zip(sites, rows):
                R = momentum(combine(site["baseline"], site["modes"], a))[0]
                pieces.extend((R/site["original_norm"]).tolist())
                pieces.append(cone_weight*(row["g"]+0.5))
                pieces.append(cone_weight*max(0.0, -row["F"]/
                                               max(site["original"]["F"], 1e-30)))
            return np.r_[pieces, 0.005*(a-source)]

        fit = least_squares(objective, source, bounds=(-20.0, 20.0),
                            max_nfev=150, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        rows = features(fit.x)
        candidates.append(dict(cone_weight=cone_weight,
                               amplitudes=fit.x.tolist(), fit_nfev=fit.nfev,
                               fit_success=bool(fit.success),
                               coefficient_relative_change=float(np.linalg.norm(fit.x-source)/np.linalg.norm(source)),
                               geometry=rows,
                               both_lambda_positive=all(row["normalized_lambda_squared"]>0 for row in rows),
                               peak_momentum_inflation_max=float(max(
                                   row["momentum_norm"]/site["original_norm"]
                                   for row, site in zip(rows, sites)))))
    feasible = [c for c in candidates if c["both_lambda_positive"]]
    selected = min(feasible, key=lambda c:c["peak_momentum_inflation_max"]) if feasible else min(
        candidates, key=lambda c:sum(max(0.0, -r["normalized_lambda_squared"]) for r in c["geometry"]))
    field = WideJointModes(cached, selected["amplitudes"])
    direct = []
    for site in sites:
        point = np.array([[site["radius"], 0.0, site["z"]]])
        velocity, gradient, residual = operator(field, point, site["tau"])
        F = velocity[0, 1]/site["radius"]
        shear = np.array([gradient[0, 1, 0]-F, gradient[0, 2, 0]])
        norm = np.linalg.norm(shear)
        N = shear/norm
        K = np.array([-N[1], N[0]])
        lam2 = -2.0*F*N[0]*(2.0*F*N[0]+norm)
        row = dict(k=site["k"], momentum_norm=float(np.linalg.norm(residual[0])),
                   lambda_squared=float(lam2), cone_pass=False)
        if lam2 > 0:
            target = stress_primitive(field, site["radius"], site["z"], site["tau"], order=12)
            dot_n, dot_k = float(target@N), float(target@K)
            row.update(target=target.tolist(), target_dot_N=dot_n,
                       target_dot_K=dot_k)
            if abs(dot_n)>1e-14 and abs(2.0*F*N[0])>1e-14:
                ratio = float(abs(np.sqrt(lam2)/(2.0*F*N[0])*dot_k/dot_n))
                row["cone_ratio"] = ratio
                row["cone_pass"] = bool(dot_n<0 and ratio<1)
        direct.append(row)
    report = dict(
        source="Adaptive 24-mode bridge joint momentum/cone-geometry screen",
        sites=[{key:site[key] for key in ("k", "radius", "z", "original")}
               for site in sites],
        candidates=candidates,
        selected_cone_weight=selected["cone_weight"],
        selected_direct=direct,
        scope="Only two former maximum-residual locations, physical local stress-cone analogue. Opening lambda squared at isolated points does not establish target orientation, an open support patch, moment matching or full momentum improvement over the bridge.",
        accepted=False, pde_validated=False,
    )
    output = ROOT/"adaptive_bridge_cone_repair.json"
    output.write_bytes((json.dumps(report, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), candidates=[
        {key:c[key] for key in ("cone_weight", "both_lambda_positive", "peak_momentum_inflation_max")}
        for c in candidates], direct=direct), indent=2))
    return report


if __name__ == "__main__":
    run()
