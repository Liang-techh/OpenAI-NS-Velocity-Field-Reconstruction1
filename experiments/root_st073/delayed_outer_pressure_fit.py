"""Fit realizable compact pressure after ten-mode pointwise admission.

Cone inequalities are linear in the six pressure coefficients. Among fields
passing those sampled inequalities, minimize full Cartesian momentum on the
cone centers and nearby-time nodes. No continuous cone or wave is implied.
"""

import json

import numpy as np
from scipy.optimize import LinearConstraint, linprog, minimize

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, stats
from delayed_multimode_cone_fit import (BASE_NAME, cone,
                                         evaluate_holdout, load_cache)
from delayed_multimode_pressure_admission import geometry
from delayed_outer_radial_extension import (CACHE_NAME, NEARBY_CACHE_NAME,
                                             signature)
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    c = np.asarray(source['coefficients'])
    rows, _ = load_cache(ROOT/CACHE_NAME, signature())
    with np.load(ROOT/NEARBY_CACHE_NAME, allow_pickle=False) as data:
        nearby = {key: data[key] for key in ('tau', 'R0', 'L', 'Q')}
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    pressure = SimilarityPressurePatch(base, np.zeros(6))
    tau = .5*2**(-5.5)
    train_points = base.compact.joined.inner.from_similarity(
        np.array([row['X'] for row in rows]),
        np.array([row['eta'] for row in rows]), tau)
    A, b = [], []
    stress_columns = []
    for row, point in zip(rows, train_points):
        target, N, K, multiplier, lam2 = geometry(row, c)
        if lam2 <= 0 or multiplier is None:
            raise ValueError('Selected mean has invalid cone shear')
        columns = pressure_stress_columns(pressure, point, tau)
        stress_columns.append(columns)
        for sign in (-1., 1.):
            direction = .8*N+sign*multiplier*K
            A.append(direction@columns)
            b.append(-.1-direction@target)
    A, b = np.asarray(A), np.asarray(b)
    stress_columns = np.asarray(stress_columns)
    coefficient_scale = 1000.
    lp = linprog(np.zeros(6), A_ub=A*coefficient_scale,
                 b_ub=b, bounds=[(None, None)]*6, method='highs')
    if not lp.success:
        raise ValueError(f'Pressure cone LP infeasible: {lp.message}')

    train_residual = np.stack([
        row['R0']+row['L']@c+
        np.einsum('aij,i,j->a', row['Q'], c, c)
        for row in rows])
    _, train_gradient = pressure.basis(train_points, tau)
    nearby_tau = float(nearby['tau'])
    nearby_points, _, _ = nodes(base, (1.01, 1.015, 1.025),
                                 (.22, .28, .32), nearby_tau)
    nearby_residual = evaluate_holdout(nearby, c)
    _, nearby_gradient = pressure.basis(nearby_points, nearby_tau)
    train_scale = max(stats(train_residual)['max'], 1.)
    nearby_scale = max(stats(nearby_residual)['max'], 1.)
    matrix = np.concatenate((train_gradient.reshape(-1, 6)/train_scale,
                             nearby_gradient.reshape(-1, 6)/nearby_scale))
    vector = np.concatenate((train_residual.ravel()/train_scale,
                             nearby_residual.ravel()/nearby_scale))
    matrix *= coefficient_scale

    def objective(y):
        defect = vector+matrix@y
        return float(.5*defect@defect+1e-8*y@y)

    def gradient(y):
        defect = vector+matrix@y
        return matrix.T@defect+2e-8*y

    fit = minimize(objective, lp.x, jac=gradient, method='SLSQP',
                   constraints=LinearConstraint(A*coefficient_scale,
                                                -np.inf, b),
                   options=dict(maxiter=1000, ftol=1e-11))
    fitted_y = fit.x if (np.all(np.isfinite(fit.x)) and
                         np.max(A*coefficient_scale@fit.x-b) < 1e-6
                         and objective(fit.x) <= objective(lp.x)) else lp.x
    coefficients = coefficient_scale*fitted_y
    train_after = train_residual+train_gradient@coefficients
    nearby_after = nearby_residual+nearby_gradient@coefficients
    inequalities = b-A@coefficients
    residual_scale = 1e8
    gradient_all = np.concatenate((train_gradient.reshape(-1, 6),
                                   nearby_gradient.reshape(-1, 6)))
    residual_all = np.concatenate((train_residual.ravel(),
                                   nearby_residual.ravel()))
    normalized_gradient = gradient_all*coefficient_scale/residual_scale
    normalized_residual = residual_all/residual_scale
    max_lp_matrix = np.column_stack((normalized_gradient,
                                     -np.ones(len(residual_all))))
    minmax = linprog(
        np.r_[np.zeros(6), 1.],
        A_ub=np.vstack((np.column_stack((A*coefficient_scale,
                                         np.zeros(len(b)))),
                        max_lp_matrix,
                        np.column_stack((-normalized_gradient,
                                         -np.ones(len(residual_all)))))),
        b_ub=np.r_[b, -normalized_residual, normalized_residual],
        bounds=[(None, None)]*6+[(0, None)], method='highs')
    minmax_coefficients = (coefficient_scale*minmax.x[:6]
                           if minmax.success else None)
    minmax_actual_max = (float(np.max(np.abs(
        residual_all+gradient_all@minmax_coefficients)))
        if minmax.success else None)
    minmax_cone_slack = (float(np.min(b-A@minmax_coefficients))
                         if minmax.success else None)
    wider_support = []
    for name, active in (('second_and_third_radial', (2, 3, 4, 5)),
                         ('second_radial_only', (2, 3)),
                         ('third_radial_only', (4, 5))):
        restricted = linprog(
            np.zeros(len(active)),
            A_ub=A[:, active]*coefficient_scale, b_ub=b,
            bounds=[(None, None)]*len(active), method='highs')
        item = dict(active_pressure_indices=active, name=name,
                    feasible=bool(restricted.success))
        if restricted.success:
            trial = np.zeros(6)
            trial[list(active)] = coefficient_scale*restricted.x
            item.update(coefficients=trial.tolist(),
                        min_cone_slack=float(np.min(b-A@trial)),
                        train_max=stats(train_residual+
                                        train_gradient@trial)['max'],
                        nearby_max=stats(nearby_residual+
                                         nearby_gradient@trial)['max'])
        wider_support.append(item)
    report = dict(source='delayed_outer_admission_search.json',
                  pressure_coefficients=coefficients.tolist(),
                  lp_coefficients=(coefficient_scale*lp.x).tolist(),
                  lp_feasible=bool(lp.success),
                  optimizer_success=bool(fit.success),
                  optimizer_message=str(fit.message),
                  min_cone_inequality_slack=float(np.min(inequalities)),
                  cone_inequality_slacks=inequalities.tolist(),
                  minmax_lp_success=bool(minmax.success),
                  minmax_component_residual=(
                      float(minmax.x[-1]*residual_scale)
                      if minmax.success else None),
                  minmax_pressure_coefficients=(
                      minmax_coefficients.tolist()
                      if minmax.success else None),
                  minmax_actual_component_residual=minmax_actual_max,
                  minmax_cone_slack=minmax_cone_slack,
                  wider_pressure_support=wider_support,
                  training=dict(before=stats(train_residual),
                                after=stats(train_after)),
                  nearby=dict(before=stats(nearby_residual),
                              after=stats(nearby_after)),
                  scope='Six compact pressure modes satisfying twelve '
                        'sampled strict .8/.1 cone inequalities, with '
                        'full Cartesian momentum cost at cone centers '
                        'and nearby-time points. No continuous cone, '
                        'moment closure, wave, volume L2, or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_outer_pressure_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(lp_feasible=report['lp_feasible'],
                          optimizer_success=report['optimizer_success'],
                          min_cone_slack=(
                              report['min_cone_inequality_slack']),
                          coefficients=report['pressure_coefficients'],
                          minmax_component_residual=(
                              report['minmax_component_residual']),
                          wider_pressure_support=wider_support,
                          training=report['training'],
                          nearby=report['nearby'])), flush=True)


if __name__ == '__main__':
    run()
