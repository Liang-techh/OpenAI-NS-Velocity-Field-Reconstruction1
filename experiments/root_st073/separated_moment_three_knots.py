"""Close sampled outer moments at three scales with smooth knot transfer."""

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


KNOTS = (11.0, 15.0, 19.0)


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    previous = json.loads((ROOT / "separated_moment_three_window.json").read_text())
    old = np.asarray(previous["constrained_l4"]["amplitudes"]).reshape(2, 2, 3, 2)
    start = np.stack((old[0], 0.5 * (old[0] + old[1]), old[1])).ravel()
    count = len(start)

    def field(a):
        return SeparatedMomentModes(base, a, windows=RADIAL_WINDOWS_THREE,
                                    knots=KNOTS)

    zero = field(np.zeros(count))
    units = [field(np.eye(count)[j]) for j in range(count)]
    moment_data = moment_slices(inner, base, zero, orders=KNOTS,
                                unit_fields=units)
    grids = [grid_slice(inner, zero, units, k) for k in KNOTS]
    raw = np.array([outer_moments(s, np.zeros(count)) for s in moment_data])
    scales = np.empty_like(raw)
    for i in range(len(KNOTS)):
        for component in (0, 1):
            scales[i, component::2] = max(np.max(np.abs(raw[i, component::2])), 1e-12)

    def normalized_moments(a):
        return (np.array([outer_moments(s, a) for s in moment_data])
                / scales).ravel()

    def grid_norms(a):
        return [np.linalg.norm(momentum(combine(s["baseline"], s["modes"], a)),
                               axis=1) / s["scale"] for s in grids]

    exact = least_squares(normalized_moments, start,
                          bounds=(-100.0, 100.0), max_nfev=120,
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

    transfer = []
    candidate = field(constrained.x)
    for k in (11.0, 13.0, 15.0, 17.0, 19.0):
        tau = 0.5 * 2.0**-k
        points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
        residual, divergence = independent_fd(
            candidate, points, tau,
            0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
        stress_rows = []
        for eta in (-0.2, 0.2):
            end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                        [eta], tau)[0]
            args = (float(end[0]), float(end[2]), tau)
            current = stress_primitive(candidate, *args)
            original = stress_primitive(base, *args)
            stress_rows.append(dict(
                eta=eta, outer_stress=current.tolist(),
                norm_ratio_to_unrepaired=float(np.linalg.norm(current)
                                               / np.linalg.norm(original))))
        transfer.append(dict(k=k,
                             direct_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                             direct_divergence_max=float(np.max(np.abs(divergence))),
                             outer_stress=stress_rows))

    report = dict(
        source="Three scale knots and three separated radial windows",
        knots=list(KNOTS), radial_windows=[list(x) for x in RADIAL_WINDOWS_THREE],
        mode_count=count,
        initial=summary(start),
        exact={**summary(exact.x), "optimizer_success": bool(exact.success),
               "nfev": exact.nfev},
        constrained_l4={**summary(constrained.x),
                        "optimizer_success": bool(constrained.success),
                        "optimizer_message": str(constrained.message),
                        "iterations": int(constrained.nit)},
        transfer=transfer,
        scope="Twelve sampled physical tangential moment equalities at k=11,15,19 and momentum training at those scales; k=13,17 holdouts. Smooth finite knot interpolation, not an infinite scale-recursive law or the paper's five profile moments.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_three_knots.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          initial={key: report["initial"][key] for key in (
                              "normalized_moment_max", "grid_maxima")},
                          exact={key: report["exact"][key] for key in (
                              "normalized_moment_max", "grid_maxima", "optimizer_success")},
                          constrained_l4={key: report["constrained_l4"][key] for key in (
                              "normalized_moment_max", "grid_maxima",
                              "optimizer_success", "iterations")},
                          transfer=[dict(k=row["k"],
                                         direct_momentum_max=row["direct_momentum_max"],
                                         outer_ratios=[s["norm_ratio_to_unrepaired"]
                                                       for s in row["outer_stress"]])
                                    for row in transfer]), indent=2))
    return report


if __name__ == "__main__":
    run()
