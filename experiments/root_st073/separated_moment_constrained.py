"""Optimize momentum within the separated modes' closed-moment manifold."""

import json

import numpy as np
from scipy.optimize import minimize

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_multiscale_fit import sample_points
from affine_momentum import combine, momentum
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_fit import grid_slice
from separated_moment_modes import SeparatedMomentModes


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    source = json.loads((ROOT / "separated_moment_fit.json").read_text())
    start = np.asarray(source["exact_moment"]["amplitudes"])
    zero = SeparatedMomentModes(base, np.zeros(16))
    units = [SeparatedMomentModes(base, np.eye(16)[j]) for j in range(16)]
    moments_data = moment_slices(inner, base, zero, unit_fields=units)
    grids = [grid_slice(inner, zero, units, k) for k in (11.0, 15.0, 19.0)]
    initial_moments = np.array([outer_moments(s, np.zeros(16))
                                for s in moments_data])
    moment_scales = np.empty_like(initial_moments)
    for i in range(2):
        for component in (0, 1):
            moment_scales[i, component::2] = max(
                np.max(np.abs(initial_moments[i, component::2])), 1e-12)

    def normalized_moments(a):
        return (np.array([outer_moments(s, a) for s in moments_data])
                / moment_scales).ravel()

    def grid_norms(a):
        return [np.linalg.norm(momentum(combine(s["baseline"], s["modes"], a)),
                               axis=1) / s["scale"] for s in grids]

    def objective(a, power):
        values = np.concatenate(grid_norms(a))
        return float(np.mean(values**power) + 1e-7 * np.dot(a, a))

    candidates = []
    for power in (2, 4):
        result = minimize(lambda a: objective(a, power), start,
                          method="SLSQP", bounds=[(-100.0, 100.0)] * 16,
                          constraints=[dict(type="eq", fun=normalized_moments)],
                          options=dict(maxiter=350, ftol=1e-12))
        a = result.x
        norms = grid_norms(a)
        candidates.append(dict(
            power=power, amplitudes=a.tolist(),
            optimizer_success=bool(result.success),
            optimizer_message=str(result.message), iterations=int(result.nit),
            normalized_moment_max=float(np.max(np.abs(normalized_moments(a)))),
            normalized_grid_maxima=[float(np.max(r)) for r in norms],
            grid_maxima=[float(np.max(r) * s["scale"])
                         for r, s in zip(norms, grids)],
            objective=float(objective(a, power)),
        ))

    reference = dict(
        exact_closed=source["exact_moment"]["grid_maxima"],
        unmodified=source["initial"]["grid_maxima"],
    )
    holdouts = []
    for name, amplitudes in (("exact_closed", start),
                             ("constrained_l2", candidates[0]["amplitudes"]),
                             ("constrained_l4", candidates[1]["amplitudes"])):
        field = SeparatedMomentModes(base, amplitudes)
        rows = []
        for k in (13.0, 17.0):
            tau = 0.5 * 2.0**-k
            points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
            residual, divergence = independent_fd(
                field, points, tau,
                0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
            stress = []
            for eta in (-0.2, 0.2):
                end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                            [eta], tau)[0]
                stress.append(stress_primitive(
                    field, float(end[0]), float(end[2]), tau).tolist())
            rows.append(dict(k=k,
                             direct_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                             direct_divergence_max=float(np.max(np.abs(divergence))),
                             outer_stress=stress))
        holdouts.append(dict(name=name, rows=rows))

    report = dict(
        source="Constrained full-momentum optimization in 16 separated moment modes",
        reference=reference, candidates=candidates, holdouts=holdouts,
        scope="Exact sampled physical outer-moment equality at k=11 and 19; momentum trained at k=11,15,19, with k=13,17 holdouts. Neither the paper's five profile moments nor a global scale-recursive or PDE certificate.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_constrained.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          candidates=[{key: c[key] for key in (
                              "power", "optimizer_success", "optimizer_message",
                              "iterations", "normalized_moment_max",
                              "normalized_grid_maxima", "grid_maxima")}
                              for c in candidates], holdouts=holdouts), indent=2))
    return report


if __name__ == "__main__":
    run()
