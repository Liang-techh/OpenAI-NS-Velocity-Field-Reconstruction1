"""Joint sampled outer-moment and physical rectangular stress-cone repair."""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from affine_momentum import combine, jets, momentum
from joined_field import coordinates
from midplane_axial_cone_knob import cone_row
from midplane_physical_covariance_pairs import support as physical_support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def cone_cache(inner, field, make_delta, k, support, offsets, order=12):
    if k in (11, 19):
        source = json.loads((ROOT / 'compact_potential' /
                             f'midplane_wave_source_k{k}.json').read_text())
    else:
        tau = .5 * 2.**-k
        X = inner.p.X_max * (1 + 15 * .325)**2
        source = dict(tau=tau, point=inner.from_similarity([X], [-.0125], tau)[0])
        support = physical_support(inner, source, inner.nu)
    tau = float(source['tau'])
    radius, _, height = source['point']
    centers = [(radius + x * support['radial_halfwidth'], 0.,
                height + z * 1.4 * support['axial_halfwidth'])
               for x in offsets for z in offsets]
    points = list(centers)
    panels = []
    g, w = leggauss(order)
    for r, _, z in centers:
        q = float(coordinates(0., z / np.sqrt(inner.nu), tau, inner.h)['q'])
        ri = np.sqrt(2 * inner.nu * q * inner.p.X_max)
        cuts = sorted(set(np.clip([0., ri, r], 0., r)))
        start = len(points)
        radii, weights = [], []
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            rr = (hi + lo) / 2 + (hi - lo) / 2 * g
            ww = (hi - lo) / 2 * w
            points.extend((float(x), 0., float(z)) for x in rr)
            radii.extend(rr)
            weights.extend(ww)
        panels.append((slice(start, len(points)), np.array(radii),
                       np.array(weights), r))
    points = np.array(points)
    args = (points, tau, .0005 * np.sqrt(inner.nu * tau), .0001 * tau)
    baseline = jets(field, *args)
    offset = {11: 0, 15: 12, 19: 24}[k]
    indices = np.array([offset + j for j in (0, 1, 6, 7)])
    changes = []
    for index in indices:
        changes.append(jets(make_delta(index), *args))
    modes = tuple(np.stack([change[i] for change in changes]) for i in range(3))
    return dict(k=k, tau=tau, centers=centers, baseline=baseline,
                modes=modes, indices=indices, panels=panels)


def cached_cones(cache, delta):
    jet = combine(cache['baseline'], cache['modes'], delta[cache['indices']])
    residual = momentum(jet)
    rows = []
    for i, (panel, rr, ww, r) in enumerate(cache['panels']):
        u, J = jet[0][i], jet[1][i]
        F = u[1] / r
        shear = np.array([J[1, 0] - F, J[2, 0]])
        magnitude = np.linalg.norm(shear)
        N = shear / max(magnitude, 1e-30)
        K = np.array([-N[1], N[0]])
        lam2 = -2 * F * N[0] * (2 * F * N[0] + magnitude)
        local = residual[panel]
        target = np.array([-np.dot(ww * rr**2, local[:, 1]) / r**2,
                           -np.dot(ww * rr, local[:, 2]) / r])
        dn, dk = target @ N, target @ K
        ratio = np.sqrt(max(lam2, 0)) * abs(dk) / max(abs(2 * F * N[0] * dn), 1e-30)
        rows.append((lam2, dn, ratio, np.linalg.norm(residual[i])))
    return np.array(rows)


class ZeroBackground:
    """Metadata-compatible zero field for exact affine basis jets."""
    def __init__(self, base):
        for name in ('inner', 'nu', 'join_X', 'ratio'):
            setattr(self, name, getattr(base, name))

    def fields(self, points, tau):
        return np.zeros_like(np.asarray(points, float)), np.zeros(len(points))


def run(edge=False, resolved=False):
    inner, fields = build_fields()
    a = np.array(json.loads((ROOT / 'midplane_axial_cone_all_knots_repair.json')
                            .read_text())['amplitudes'])
    def make_field(amplitudes):
        return SeparatedMomentModes(fields['two_sided_cone'], amplitudes,
                                    windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    current = make_field(a)
    zero = ZeroBackground(fields['two_sided_cone'])
    def make_delta(index):
        return SeparatedMomentModes(zero, np.eye(len(a))[index],
                                    windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    units = [make_delta(i) for i in range(len(a))]
    breaks = sorted({v for lo, hi in RADIAL_WINDOWS_THREE for v in (lo, (lo + hi) / 2, hi)}) if resolved else None
    moments = moment_slices(inner, fields['two_sided_cone'], current,
                            orders=(11., 15., 19.), unit_fields=units,
                            unit_fields_are_deltas=True, radial_breaks=breaks,
                            n=24 if resolved else 12)
    def moment_values(delta):
        return np.concatenate([outer_moments(s, delta) for s in moments])
    m0 = moment_values(np.zeros_like(a))
    mscale = np.maximum(np.max(np.abs(np.array([
        moment_values(row) - m0 for row in np.eye(len(a))])), axis=0), 1e-8)
    supports = json.loads((ROOT / 'midplane_physical_covariance_pairs.json').read_text())
    offsets = (-.99, 0., .99) if edge else (-.8, 0., .8)
    caches = [cone_cache(inner, current, make_delta, k,
                         next((s['support'] for s in supports['scales'] if s['k'] == k), None), offsets,
                         order=48 if resolved else 12)
              for k in ((11, 15, 19) if edge else (11, 19))]
    base_cones = [cached_cones(c, np.zeros_like(a)) for c in caches]
    scale = np.maximum(np.abs(a), .02)
    def inequalities(z):
        delta = z * scale
        result = []
        for cache, baseline in zip(caches, base_cones):
            rows = cached_cones(cache, delta)
            result.extend((rows[:, 0] / np.maximum(abs(baseline[:, 0]), 1e-20) - .1,
                           -rows[:, 1] / np.maximum(abs(baseline[:, 1]), 1e-20) - .1,
                           .95 - rows[:, 2]))
        return np.concatenate(result)
    fit = minimize(lambda z: float(z @ z), np.zeros_like(a), method='SLSQP',
                   bounds=[(-2., 2.)] * len(a),
                   constraints=[dict(type='eq', fun=lambda z: moment_values(z * scale) / mscale),
                                dict(type='ineq', fun=inequalities)],
                   options=dict(maxiter=300, ftol=1e-10))
    delta = fit.x * scale
    repaired = make_field(a + delta)
    results = []
    for cache, baseline in zip(caches, base_cones):
        predicted = cached_cones(cache, delta)
        rows = [cone_row(repaired, float(r), float(z), cache['tau'], float(i),
                         stress_order=96 if resolved else 12)
                for i, (r, _, z) in enumerate(cache['centers'])]
        results.append(dict(k=cache['k'], baseline_cone_ratios=baseline[:, 2].tolist(),
                            predicted_cone_ratios=predicted[:, 2].tolist(), direct_rows=rows,
                            local_momentum_max_before=float(max(baseline[:, 3])),
                            local_momentum_max_after=float(max(predicted[:, 3]))))
        print(f"k={cache['k']}: {sum(r['cone_pass'] for r in rows)}/9 pass; "
              f"max ratio={max(r.get('cone_ratio', float('inf')) for r in rows):.6g}", flush=True)
    direct_slices = moment_slices(inner, fields['two_sided_cone'], repaired,
                                  orders=(11., 15., 19.), unit_fields=[repaired],
                                  radial_breaks=breaks, n=48 if resolved else 12)
    direct_moments = np.concatenate([outer_moments(s, np.zeros(1)) for s in direct_slices])
    report = dict(optimizer_success=bool(fit.success), message=str(fit.message),
                  training_offsets=offsets,
                  resolved_quadrature=resolved, radial_breaks=breaks,
                  fit_gauss_order=24 if resolved else 12,
                  direct_gauss_order=48 if resolved else 12,
                  iterations=int(fit.nit), objective=float(fit.fun),
                  normalized_moment_max=float(max(abs(moment_values(delta) / mscale))),
                  direct_normalized_moment_max=float(max(abs(direct_moments / mscale))),
                  direct_absolute_moment_max=float(max(abs(direct_moments))),
                  minimum_constraint=float(min(inequalities(fit.x))),
                  amplitudes=(a + delta).tolist(), scales=results,
                  scope='Joint constrained mean repair on 12 sampled outer moments and nine physical support nodes per reported scale. Direct cone replay; no continuous cone, wave residual gain, pulse dynamics, or scale recursion.',
                  accepted=False, scale_recursion_established=False)
    name = ('midplane_connected_cone_resolved_repair.json' if resolved else
            ('midplane_connected_cone_edge_repair.json' if edge else 'midplane_connected_cone_repair.json'))
    (ROOT / name).write_bytes((json.dumps(report, indent=2) + '\n').encode())
    print(json.dumps({k: report[k] for k in ('optimizer_success', 'normalized_moment_max', 'minimum_constraint')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--edge', action='store_true')
    parser.add_argument('--resolved', action='store_true')
    args = parser.parse_args()
    run(edge=args.edge or args.resolved, resolved=args.resolved)
