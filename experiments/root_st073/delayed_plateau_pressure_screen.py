"""Search smooth broad radial pressure plateaus on the ten-mode mean.

Each plateau rises inside the transition, stays nonzero over a wider
annulus, and falls smoothly before X=3. Cone admission is constrained at
the same twelve nodes; a linear program minimizes the largest full
Cartesian momentum component on training and nearby-time nodes.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import linprog
from scipy.special import beta, betainc

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, stats
from delayed_multimode_cone_fit import (BASE_NAME, evaluate_holdout,
                                         load_cache)
from delayed_multimode_pressure_admission import geometry
from delayed_outer_radial_extension import (CACHE_NAME, NEARBY_CACHE_NAME,
                                             signature)
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_similarity_pressure_screen import (AXIAL_INTERVALS,
                                                SimilarityPressurePatch)
from joined_field import coordinates
from paper_moment_bridge import bump
from radial_continuation import ROOT


RADIAL_SHAPES = (
    ((1., 1.012), (1.12, 1.2)),
    ((1., 1.02), (1.12, 1.2)),
    ((1.005, 1.04), (1.12, 1.2)),
    ((1., 1.012), (2.5, 3.)),
    ((1., 1.02), (2.5, 3.)),
    ((1.005, 1.04), (2.5, 3.)),
)


def smoothstep(value, interval):
    lo, hi = interval
    s = np.clip((value-lo)/(hi-lo), 0., 1.)
    shape = betainc(6, 6, s)
    derivative = np.where((value > lo) & (value < hi),
                          s**5*(1-s)**5/(beta(6, 6)*(hi-lo)), 0.)
    return shape, derivative


class PlateauPressurePatch:
    def __init__(self, base):
        self.base = base
        self.nu = base.nu
        self.beta = 2*base.A
        self.radial_shapes = RADIAL_SHAPES
        self.axial_intervals = AXIAL_INTERVALS
        self.count = len(RADIAL_SHAPES)*len(AXIAL_INTERVALS)

    def basis(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        radius = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(radius/sn, pts[:, 2]/sn, ts,
                         self.base.heat.h)
        q, X, eta = (np.asarray(co[key])
                     for key in ('q', 'X', 'eta'))
        factor = q**(-self.beta)
        Xr = np.asarray(co['X_r'])/sn
        Xz = np.asarray(co['X_z'])/sn
        etaz = np.asarray(co['eta_z'])/sn
        qz = np.asarray(co['q_z'])/sn
        ca = np.divide(pts[:, 0], radius, out=np.ones_like(radius),
                       where=radius > 0)
        sa = np.divide(pts[:, 1], radius, out=np.zeros_like(radius),
                       where=radius > 0)
        values, gradients = [], []
        for rise, fall in self.radial_shapes:
            up, up_d = smoothstep(X, rise)
            down, down_d = smoothstep(X, fall)
            radial = up*(1-down)
            radial_d = up_d*(1-down)-up*down_d
            for axial in self.axial_intervals:
                be, bed = bump(eta, *axial)
                value = factor*radial*be
                pr = factor*radial_d*Xr*be
                pz = factor*(radial_d*Xz*be+radial*bed*etaz
                             -self.beta*radial*be*qz/q)
                values.append(value)
                gradients.append(np.column_stack((pr*ca, pr*sa, pz)))
        return np.column_stack(values), np.stack(gradients, axis=2)


def plateau_stress_columns(patch, point, tau, order=64):
    radius, _, z = point
    q = float(coordinates(0., z/np.sqrt(patch.nu), tau,
                          patch.base.heat.h)['q'])
    cuts = [0., radius]
    radial_edges = [edge for rise, fall in patch.radial_shapes
                    for edge in (*rise, *fall)]
    cuts.extend(np.sqrt(2*patch.nu*q*np.asarray(radial_edges)))
    edges = sorted(set(np.clip(cuts, 0., radius)))
    g, w = leggauss(order)
    axial = np.zeros(patch.count)
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        rr = (lo+hi)/2+(hi-lo)/2*g
        weights = (hi-lo)/2*w
        points = np.column_stack((rr, np.zeros(order),
                                  np.full(order, z)))
        _, gradients = patch.basis(points, tau)
        axial -= (weights*rr)@gradients[:, 2, :]/radius
    return np.vstack((np.zeros_like(axial), axial))


def minmax_pressure(A, rhs, gradients, residual, active,
                    coefficient_limit=None):
    A = A[:, active]
    gradient = gradients[:, active]
    residual_scale = 1e8
    column_norm = np.maximum(np.max(np.abs(A), axis=0),
                             np.max(np.abs(gradient), axis=0)
                             /residual_scale)
    scales = 1/np.maximum(column_norm, 1e-14)
    scaled_A = A*scales
    scaled_G = gradient*scales/residual_scale
    n = len(active)
    fit = linprog(
        np.r_[np.zeros(n), 1.],
        A_ub=np.vstack((np.column_stack((scaled_A,
                                         np.zeros(len(rhs)))),
                        np.column_stack((scaled_G,
                                         -np.ones(len(residual)))),
                        np.column_stack((-scaled_G,
                                         -np.ones(len(residual)))))),
        b_ub=np.r_[rhs, -residual/residual_scale,
                   residual/residual_scale],
        bounds=([(-coefficient_limit/scale,
                  coefficient_limit/scale) for scale in scales]
                if coefficient_limit is not None else
                [(None, None)]*n)+[(0, None)], method='highs')
    report = dict(active_count=n, feasible=bool(fit.success),
                  coefficient_limit=coefficient_limit,
                  message=str(fit.message))
    if fit.success:
        coefficients = np.zeros(gradients.shape[1])
        coefficients[np.asarray(active)] = scales*fit.x[:n]
        actual = residual+gradients@coefficients
        report.update(coefficients=coefficients.tolist(),
                      optimized_component_max=float(
                          fit.x[-1]*residual_scale),
                      actual_component_max=float(np.max(np.abs(actual))),
                      min_cone_slack=float(
                          np.min(rhs-A@coefficients[np.asarray(active)])))
    return report


def run():
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    c = np.asarray(source['coefficients'])
    rows, _ = load_cache(ROOT/CACHE_NAME, signature())
    with np.load(ROOT/NEARBY_CACHE_NAME, allow_pickle=False) as data:
        nearby = {key: data[key] for key in ('tau', 'R0', 'L', 'Q')}
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    old = SimilarityPressurePatch(base, np.zeros(6))
    plateau = PlateauPressurePatch(base)
    tau = .5*2**(-5.5)
    train_points = base.compact.joined.inner.from_similarity(
        np.array([row['X'] for row in rows]),
        np.array([row['eta'] for row in rows]), tau)
    A, rhs = [], []
    for row, point in zip(rows, train_points):
        target, N, K, multiplier, lam2 = geometry(row, c)
        if lam2 <= 0 or multiplier is None:
            raise ValueError('Selected mean has invalid cone shear')
        columns = np.column_stack((
            pressure_stress_columns(old, point, tau),
            plateau_stress_columns(plateau, point, tau)))
        for sign in (-1., 1.):
            direction = .8*N+sign*multiplier*K
            A.append(direction@columns)
            rhs.append(-.1-direction@target)
    A, rhs = np.asarray(A), np.asarray(rhs)
    train_residual = np.stack([
        row['R0']+row['L']@c+
        np.einsum('aij,i,j->a', row['Q'], c, c)
        for row in rows])
    nearby_tau = float(nearby['tau'])
    nearby_points, _, _ = nodes(base, (1.01, 1.015, 1.025),
                                 (.22, .28, .32), nearby_tau)
    nearby_residual = evaluate_holdout(nearby, c)
    _, old_train_gradient = old.basis(train_points, tau)
    _, new_train_gradient = plateau.basis(train_points, tau)
    _, old_near_gradient = old.basis(nearby_points, nearby_tau)
    _, new_near_gradient = plateau.basis(nearby_points, nearby_tau)
    gradients = np.concatenate((
        np.concatenate((old_train_gradient, new_train_gradient),
                       axis=2).reshape(-1, 6+plateau.count),
        np.concatenate((old_near_gradient, new_near_gradient),
                       axis=2).reshape(-1, 6+plateau.count)))
    residual = np.r_[train_residual.ravel(), nearby_residual.ravel()]
    searches = {}
    for name, active in (('existing_six', list(range(6))),
                         ('plateau_twelve', list(range(6, 18))),
                         ('combined_eighteen', list(range(18))),
                         ('combined_cap_1e2', list(range(18))),
                         ('combined_cap_5e2', list(range(18))),
                         ('combined_cap_1e3', list(range(18))),
                         ('combined_cap_1e4', list(range(18))),
                         ('combined_cap_1e5', list(range(18)))):
        limit = (1e2 if name.endswith('1e2') else
                 5e2 if name.endswith('5e2') else
                 1e3 if name.endswith('1e3') else
                 1e4 if name.endswith('1e4') else
                 1e5 if name.endswith('1e5') else None)
        result = minmax_pressure(A, rhs, gradients, residual,
                                 active, coefficient_limit=limit)
        if result['feasible']:
            pressure_c = np.asarray(result['coefficients'])
            train_after = (train_residual+
                           np.concatenate((old_train_gradient,
                                           new_train_gradient),
                                          axis=2)@pressure_c)
            nearby_after = (nearby_residual+
                            np.concatenate((old_near_gradient,
                                            new_near_gradient),
                                           axis=2)@pressure_c)
            result.update(train_max=stats(train_after)['max'],
                          nearby_max=stats(nearby_after)['max'])
        searches[name] = result
        print(json.dumps(dict(name=name, feasible=result['feasible'],
                              component_max=result.get(
                                  'actual_component_max'),
                              train_max=result.get('train_max'),
                              nearby_max=result.get('nearby_max'))),
              flush=True)
    report = dict(source='delayed_outer_admission_search.json',
                  radial_shapes=RADIAL_SHAPES,
                  axial_intervals=AXIAL_INTERVALS,
                  shape='C5 beta(6,6) rise, plateau, C5 fall',
                  training_before=stats(train_residual),
                  nearby_before=stats(nearby_residual),
                  searches=searches,
                  scope='Twelve new smooth compact radial plateau '
                        'pressure modes plus six existing modes, '
                        'constrained sampled strict cone and full '
                        'Cartesian component minmax on training and '
                        'nearby-time nodes. No continuous cone, '
                        'pressure-Poisson compatibility, moment closure, '
                        'wave, volume L2, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_plateau_pressure_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
