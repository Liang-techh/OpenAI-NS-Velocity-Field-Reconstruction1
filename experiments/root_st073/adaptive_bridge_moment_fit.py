"""Jointly fit bridge momentum, physical radial moments and cone geometry.

The physical moment constraints here are a diagnostic analogue of, not a
replacement for, the paper's normalized five-moment profile construction.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_scale_knots import assemble_slice
from affine_momentum import combine, jets, momentum
from radial_continuation import ROOT
from wide_modes import WideJointModes


def moment_slices(inner, base, zero, orders=(11.0, 19.0), n=12,
                  unit_fields=None):
    nodes, weights = leggauss(n)
    output = []
    for k in orders:
        tau = 0.5 * 2.0**-k
        points = []
        panels = []
        for eta in (-0.2, 0.2):
            ends = inner.from_similarity(
                inner.p.X_max * np.array([1.0, 16.0**2]),
                [eta, eta], tau)
            ri, ro = map(float, ends[:, 0])
            z = float(ends[0, 2])
            start = len(points)
            radial_weights = []
            for lo, hi in ((0.0, ri), (ri, ro)):
                rr = 0.5 * (lo + hi) + 0.5 * (hi - lo) * nodes
                ww = 0.5 * (hi - lo) * weights
                points.extend((float(r), 0.0, z) for r in rr)
                radial_weights.extend(ww.tolist())
            panels.append(dict(eta=eta, outer_radius=ro,
                               slice=slice(start, len(points)),
                               radii=np.asarray([p[0] for p in points[start:]]),
                               weights=np.asarray(radial_weights)))
        points = np.asarray(points)
        args = (points, tau, 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
        baseline = jets(zero, *args)
        changes = []
        fields = unit_fields if unit_fields is not None else [
            WideJointModes(base, np.eye(24)[j]) for j in range(24)]
        for field in fields:
            sample = jets(field, *args)
            changes.append(tuple(x - y for x, y in zip(sample, baseline)))
        modes = tuple(np.stack([change[i] for change in changes])
                      for i in range(3))
        output.append(dict(k=k, tau=tau, panels=panels,
                           baseline=baseline, modes=modes))
    return output


def outer_moments(slice_data, amplitudes):
    residual = momentum(combine(slice_data["baseline"],
                                slice_data["modes"], amplitudes))
    result = []
    for panel in slice_data["panels"]:
        rr = panel["radii"]
        ww = panel["weights"]
        local = residual[panel["slice"]]
        ro = panel["outer_radius"]
        result.extend((-np.dot(ww * rr**2, local[:, 1]) / ro**2,
                       -np.dot(ww * rr, local[:, 2]) / ro))
    return np.asarray(result)


def hotspot_geometry(jet, radius):
    u, grad, _ = jet
    F = u[0, 1] / radius
    shear = np.array([grad[0, 1, 0] - F, grad[0, 2, 0]])
    return float(2.0 * F * shear[0] / max(np.dot(shear, shear), 1e-30))


def run():
    inner, fields = build_fields()
    base = fields["unmodified"].base
    selected_report = json.loads((ROOT / "adaptive_bridge_cone_biside_fit.json").read_text())
    selected = next(row for row in selected_report["candidates"]
                    if row["cone_weight"] == selected_report["selected_cone_weight"])
    start = np.asarray(selected["amplitudes"])
    zero = fields["unmodified"]
    grids = [assemble_slice(base, zero, inner, 16.0, k)
             for k in (11.0, 19.0)]
    moments = moment_slices(inner, base, zero)
    raw = np.array([outer_moments(s, start) for s in moments])
    # Give theta and axial moments comparable weight at each scale.
    scales = np.maximum(np.max(np.abs(raw), axis=1)[:, None], 1e-12)
    scales = np.broadcast_to(scales, raw.shape).copy()
    for i in range(len(moments)):
        for component in (0, 1):
            scales[i, component::2] = max(
                np.max(np.abs(raw[i, component::2])), 1e-12)

    def evaluate(a):
        mom = np.array([outer_moments(s, a) for s in moments])
        grid = [momentum(combine(s["baseline"], s["modes"], a))
                for s in grids]
        geometry = []
        for s in grids:
            jet = combine(s["baseline"], s["modes"], a)
            for index in (0, 10):
                point_jet = tuple(part[index:index + 1] for part in jet)
                geometry.append(hotspot_geometry(
                    point_jet, float(s["args"][0][index, 0])))
        return mom, grid, np.asarray(geometry)

    def normalized_moment_vector(a):
        return (np.array([outer_moments(s, a) for s in moments])
                / scales).ravel()

    step = 1e-4
    moment_jacobian = np.column_stack([
        (normalized_moment_vector(start + step * np.eye(24)[j])
         - normalized_moment_vector(start - step * np.eye(24)[j]))
        / (2 * step) for j in range(24)])
    singular_values = np.linalg.svd(moment_jacobian, compute_uv=False)

    candidates = []
    for weight in (0.5, 2.0, 8.0):
        def objective(a):
            mom, grid, geometry = evaluate(a)
            normalized_grid = [residual.ravel() / s["scale"]
                               for residual, s in zip(grid, grids)]
            return np.concatenate([*normalized_grid,
                                   weight * (mom / scales).ravel(),
                                   geometry + 0.5,
                                   0.003 * (a - start)])

        fit = least_squares(objective, start, bounds=(-40.0, 40.0),
                            max_nfev=120, ftol=1e-9, xtol=1e-9, gtol=1e-9)
        mom, grid, geometry = evaluate(fit.x)
        candidates.append(dict(
            weight=weight, amplitudes=fit.x.tolist(), nfev=fit.nfev,
            optimizer_success=bool(fit.success),
            outer_moments=mom.tolist(),
            normalized_moment_max=float(np.max(np.abs(mom / scales))),
            grid_maxima=[float(np.max(np.linalg.norm(r, axis=1))) for r in grid],
            hotspot_g=geometry.tolist(),
            all_hotspots_lambda_positive=bool(np.all((geometry > -1) & (geometry < 0))),
        ))

    def moment_only_objective(a):
        return np.concatenate([normalized_moment_vector(a),
                               0.001 * (a - start)])

    moment_only = least_squares(moment_only_objective, start,
                                bounds=(-40.0, 40.0), max_nfev=300,
                                ftol=1e-9, xtol=1e-9, gtol=1e-9)
    moment_only_values, moment_only_grid, moment_only_geometry = evaluate(moment_only.x)
    moment_only_report = dict(
        amplitudes=moment_only.x.tolist(), nfev=moment_only.nfev,
        optimizer_success=bool(moment_only.success),
        outer_moments=moment_only_values.tolist(),
        normalized_moment_max=float(np.max(np.abs(moment_only_values / scales))),
        grid_maxima=[float(np.max(np.linalg.norm(r, axis=1)))
                     for r in moment_only_grid],
        hotspot_g=moment_only_geometry.tolist(),
        all_hotspots_lambda_positive=bool(np.all(
            (moment_only_geometry > -1) & (moment_only_geometry < 0))),
    )
    report = dict(
        source="24 fixed compact mean modes fitted to full momentum, four physical outer moment vectors and four hotspot cone geometries",
        initial_outer_moments=raw.tolist(), moment_scales=scales.tolist(),
        initial_normalized_moment_jacobian_singular_values=singular_values.tolist(),
        candidates=candidates, moment_only=moment_only_report,
        scope="Two scales, two axial positions, 12 Gauss nodes on each of core and bridge radial panels for moment surrogate; 15 momentum grid points per scale. Physical moment analogue only. No off-grid cone, off-scale holdout, direct post-fit radial integral, paper five-moment closure, or PDE acceptance yet.",
        accepted=False, pde_validated=False,
    )
    output = ROOT / "adaptive_bridge_moment_fit.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), initial_outer_moments=raw.tolist(),
                          candidates=[{key: candidate[key] for key in (
                              "weight", "nfev", "optimizer_success",
                              "normalized_moment_max", "grid_maxima",
                              "hotspot_g", "all_hotspots_lambda_positive")}
                              for candidate in candidates],
                          moment_only={key: moment_only_report[key] for key in (
                              "nfev", "optimizer_success", "normalized_moment_max",
                              "grid_maxima", "hotspot_g",
                              "all_hotspots_lambda_positive")},
                          singular_values=singular_values.tolist()), indent=2))
    return report


if __name__ == "__main__":
    run()
