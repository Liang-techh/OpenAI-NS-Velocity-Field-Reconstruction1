"""Test an additional independent radial window under exact sampled moments."""

import json

import numpy as np
from scipy.optimize import least_squares, minimize

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_multiscale_fit import sample_points
from affine_momentum import combine, momentum
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_fit import grid_slice
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    count = 2 * 2 * len(RADIAL_WINDOWS_THREE) * 2

    def field(a):
        return SeparatedMomentModes(base, a, windows=RADIAL_WINDOWS_THREE)

    zero = field(np.zeros(count))
    units = [field(np.eye(count)[j]) for j in range(count)]
    moment_data = moment_slices(inner, base, zero, unit_fields=units)
    grids = [grid_slice(inner, zero, units, k) for k in (11.0, 15.0, 19.0)]
    raw = np.array([outer_moments(s, np.zeros(count)) for s in moment_data])
    scales = np.empty_like(raw)
    for i in range(2):
        for component in (0, 1):
            scales[i, component::2] = max(np.max(np.abs(raw[i, component::2])), 1e-12)

    def normalized_moments(a):
        return (np.array([outer_moments(s, a) for s in moment_data])
                / scales).ravel()

    def grid_norms(a):
        return [np.linalg.norm(momentum(combine(s["baseline"], s["modes"], a)),
                               axis=1) / s["scale"] for s in grids]

    first = least_squares(
        lambda a: np.r_[normalized_moments(a), 0.001 * a],
        np.zeros(count), bounds=(-100.0, 100.0), max_nfev=200,
        ftol=1e-10, xtol=1e-10, gtol=1e-10)
    exact = least_squares(normalized_moments, first.x,
                          bounds=(-100.0, 100.0), max_nfev=100,
                          ftol=1e-11, xtol=1e-11, gtol=1e-11)

    def objective(a):
        values = np.concatenate(grid_norms(a))
        return float(np.mean(values**4) + 1e-7 * np.dot(a, a))

    constrained = minimize(
        objective, exact.x, method="SLSQP",
        bounds=[(-100.0, 100.0)] * count,
        constraints=[dict(type="eq", fun=normalized_moments)],
        options=dict(maxiter=400, ftol=1e-12))

    def summary(a):
        norms = grid_norms(a)
        return dict(amplitudes=a.tolist(),
                    normalized_moment_max=float(np.max(np.abs(normalized_moments(a)))),
                    grid_maxima=[float(np.max(r) * s["scale"])
                                 for r, s in zip(norms, grids)],
                    normalized_grid_maxima=[float(np.max(r)) for r in norms])

    holdouts = []
    for name, a in (("exact", exact.x), ("constrained_l4", constrained.x)):
        current = field(a)
        rows = []
        for k in (13.0, 17.0):
            tau = 0.5 * 2.0**-k
            points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
            residual, divergence = independent_fd(
                current, points, tau,
                0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
            rows.append(dict(k=k,
                             direct_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                             direct_divergence_max=float(np.max(np.abs(divergence)))))
        holdouts.append(dict(name=name, rows=rows))

    direct_outer = []
    for k in (11.0, 19.0):
        tau = 0.5 * 2.0**-k
        for eta in (-0.2, 0.2):
            end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                        [eta], tau)[0]
            target = stress_primitive(field(constrained.x),
                                      float(end[0]), float(end[2]), tau)
            direct_outer.append(dict(k=k, eta=eta,
                                     outer_stress=target.tolist()))

    report = dict(
        source="Third disjoint radial window plus exactly constrained physical outer moments",
        radial_windows=[list(x) for x in RADIAL_WINDOWS_THREE],
        mode_count=count,
        initial=summary(np.zeros(count)),
        exact={**summary(exact.x), "optimizer_success": bool(exact.success),
               "nfev": exact.nfev},
        constrained_l4={**summary(constrained.x),
                        "optimizer_success": bool(constrained.success),
                        "optimizer_message": str(constrained.message),
                        "iterations": int(constrained.nit)},
        holdouts=holdouts, direct_outer_stress=direct_outer,
        scope="Three smooth radial windows, two scale knots. Eight physical tangential moments constrained at k=11,19; momentum trained at k=11,15,19; k=13,17 holdouts. No paper five-moment or full-domain PDE acceptance.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_three_window.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          initial={k: report["initial"][k] for k in (
                              "normalized_moment_max", "grid_maxima")},
                          exact={k: report["exact"][k] for k in (
                              "normalized_moment_max", "grid_maxima",
                              "optimizer_success")},
                          constrained_l4={k: report["constrained_l4"][k] for k in (
                              "normalized_moment_max", "grid_maxima",
                              "optimizer_success", "iterations")},
                          holdouts=holdouts,
                          direct_outer_stress=direct_outer), indent=2))
    return report


if __name__ == "__main__":
    run()
