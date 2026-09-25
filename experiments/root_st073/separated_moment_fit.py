"""Try radially separated, smooth time-transferred moment repair directions."""

import json

import numpy as np
from scipy.optimize import least_squares

from adaptive_bridge_moment_fit import moment_slices, outer_moments, hotspot_geometry
from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_multiscale_fit import sample_points
from affine_momentum import combine, jets, momentum
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import SeparatedMomentModes


def grid_slice(inner, zero, unit_fields, k):
    points, tau = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
    args = (points, tau, 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
    baseline = jets(zero, *args)
    changes = []
    for field in unit_fields:
        sample = jets(field, *args)
        changes.append(tuple(a - b for a, b in zip(sample, baseline)))
    modes = tuple(np.stack([change[i] for change in changes])
                  for i in range(3))
    scale = float(np.max(np.linalg.norm(momentum(baseline), axis=1)))
    return dict(k=k, args=args, baseline=baseline, modes=modes, scale=scale)


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    zero = SeparatedMomentModes(base, np.zeros(16))
    units = [SeparatedMomentModes(base, np.eye(16)[j]) for j in range(16)]
    moment_data = moment_slices(inner, base, zero, unit_fields=units)
    grids = [grid_slice(inner, zero, units, k) for k in (11.0, 19.0)]
    initial = np.array([outer_moments(s, np.zeros(16)) for s in moment_data])
    scales = np.empty_like(initial)
    for i in range(2):
        for component in (0, 1):
            scales[i, component::2] = max(
                np.max(np.abs(initial[i, component::2])), 1e-12)

    def moments(a):
        return np.array([outer_moments(s, a) for s in moment_data])

    def describe(a):
        outer = moments(a)
        residuals = [momentum(combine(s["baseline"], s["modes"], a))
                     for s in grids]
        geometry = []
        for s in grids:
            jet = combine(s["baseline"], s["modes"], a)
            for index in (0, 10):
                point_jet = tuple(part[index:index + 1] for part in jet)
                geometry.append(hotspot_geometry(
                    point_jet, float(s["args"][0][index, 0])))
        return dict(
            amplitudes=a.tolist(), outer_moments=outer.tolist(),
            normalized_moment_max=float(np.max(np.abs(outer / scales))),
            grid_maxima=[float(np.max(np.linalg.norm(r, axis=1)))
                         for r in residuals],
            hotspot_g=geometry,
            hotspot_lambda_positive=bool(all(-1.0 < g < 0.0 for g in geometry)),
        )

    step = 1e-4
    eye = np.eye(16)
    jac = np.column_stack([
        ((moments(step * eye[j]) - moments(-step * eye[j]))
         / (2.0 * step) / scales).ravel()
        for j in range(16)])
    singular_values = np.linalg.svd(jac, compute_uv=False)
    initial_summary = describe(np.zeros(16))

    def moment_objective(a):
        return np.r_[(moments(a) / scales).ravel(), 0.001 * a]

    moment_fit = least_squares(moment_objective, np.zeros(16),
                               bounds=(-100.0, 100.0), max_nfev=200,
                               ftol=1e-9, xtol=1e-9, gtol=1e-9)
    moment_summary = describe(moment_fit.x)
    moment_summary.update(nfev=moment_fit.nfev,
                          optimizer_success=bool(moment_fit.success))

    exact_fit = least_squares(
        lambda a: (moments(a) / scales).ravel(), moment_fit.x,
        bounds=(-100.0, 100.0), max_nfev=100,
        ftol=1e-11, xtol=1e-11, gtol=1e-11)
    exact_summary = describe(exact_fit.x)
    exact_summary.update(nfev=exact_fit.nfev,
                         optimizer_success=bool(exact_fit.success))

    def joint_objective(a, weight):
        outer = (moments(a) / scales).ravel()
        grid = [momentum(combine(s["baseline"], s["modes"], a)).ravel()
                / s["scale"] for s in grids]
        return np.r_[*grid, weight * outer, 0.003 * a]

    joint_fit = least_squares(lambda a: joint_objective(a, 5.0), np.zeros(16),
                              bounds=(-100.0, 100.0), max_nfev=200,
                              ftol=1e-9, xtol=1e-9, gtol=1e-9)
    joint_summary = describe(joint_fit.x)
    joint_summary.update(nfev=joint_fit.nfev,
                         optimizer_success=bool(joint_fit.success))

    near_closed_fit = least_squares(
        lambda a: joint_objective(a, 100.0), moment_fit.x,
        bounds=(-100.0, 100.0), max_nfev=200,
        ftol=1e-9, xtol=1e-9, gtol=1e-9)
    near_closed_summary = describe(near_closed_fit.x)
    near_closed_summary.update(nfev=near_closed_fit.nfev,
                               optimizer_success=bool(near_closed_fit.success))

    middle = []
    k, tau = 15.0, 0.5 * 2.0**-15
    points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
    for name, amplitudes in (("initial", np.zeros(16)),
                             ("moment_only", moment_fit.x),
                             ("exact_moment", exact_fit.x),
                             ("joint", joint_fit.x),
                             ("near_closed", near_closed_fit.x)):
        field = SeparatedMomentModes(base, amplitudes)
        residual, divergence = independent_fd(
            field, points, tau, 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
        outer = []
        for eta in (-0.2, 0.2):
            end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                        [eta], tau)[0]
            outer.append(stress_primitive(field, float(end[0]),
                                          float(end[2]), tau).tolist())
        middle.append(dict(name=name, k=k,
                           momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                           divergence_max=float(np.max(np.abs(divergence))),
                           outer_moments=outer))

    report = dict(
        source="Two scale knots and disjoint radial bumps in divergence-free swirl/poloidal bridge corrections",
        initial=initial_summary,
        normalized_moment_jacobian_singular_values=singular_values.tolist(),
        moment_only=moment_summary, exact_moment=exact_summary,
        joint=joint_summary,
        near_closed=near_closed_summary,
        middle_scale=middle,
        scope="Sixteen smooth scale-interpolated modes on two disjoint radial windows. Training k=11,19 at eta=±.2 for physical outer moments, 15 Cartesian momentum nodes per scale; k=15 holdout. Not the paper's five profile moments, full stress cone, continuous scale recursion, global energy, or PDE acceptance.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_fit.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          singular_values=singular_values.tolist(),
                          initial={key: initial_summary[key] for key in (
                              "normalized_moment_max", "grid_maxima")},
                          moment_only={key: moment_summary[key] for key in (
                              "normalized_moment_max", "grid_maxima", "nfev",
                              "optimizer_success", "hotspot_lambda_positive")},
                          exact_moment={key: exact_summary[key] for key in (
                              "normalized_moment_max", "grid_maxima", "nfev",
                              "optimizer_success", "hotspot_lambda_positive")},
                          joint={key: joint_summary[key] for key in (
                              "normalized_moment_max", "grid_maxima", "nfev",
                              "optimizer_success", "hotspot_lambda_positive")},
                          near_closed={key: near_closed_summary[key] for key in (
                              "normalized_moment_max", "grid_maxima", "nfev",
                              "optimizer_success", "hotspot_lambda_positive")},
                          middle_scale=middle), indent=2))
    return report


if __name__ == "__main__":
    run()
