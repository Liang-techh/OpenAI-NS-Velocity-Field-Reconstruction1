"""Screen a compact pressure patch for a physical stress-cone window.

The cone inequalities are linear in pressure coefficients because velocity
and shear remain fixed. A feasible window is only a precondition for the
paper's wave construction, not a momentum or pressure-Poisson solution.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import linprog, minimize

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_similarity_pressure_screen import (
    AXIAL_INTERVALS, RADIAL_INTERVALS, SimilarityPressurePatch,
)
from delayed_momentum_tangent_screen import residual, stats
from joined_field import coordinates
from radial_peak_cone import operator, stress_primitive
from radial_continuation import ROOT


BASE_NAME = 'delayed005_rise146_degree31.json'
TRAIN_X = (1.008, 1.012, 1.016, 1.02)
TRAIN_ETA = (.2, .3)
HOLDOUT_X = (1.01, 1.014, 1.018, 1.024)
HOLDOUT_ETA = (.22, .28)
CONE_RATIO_MARGIN = .8
STRESS_MARGIN = .1


def pressure_stress_columns(patch, point, tau, order=64):
    radius, _, z = point
    # q is recovered from the same physical similarity map as the field.
    q = float(coordinates(0., z/np.sqrt(patch.nu), tau,
                          patch.base.heat.h)['q'])
    cuts = [0., radius]
    for start, end in RADIAL_INTERVALS:
        cuts.extend(np.sqrt(2*patch.nu*q*np.array([start, end])))
    edges = sorted(set(np.clip(cuts, 0., radius)))
    g, w = leggauss(order)
    axial = np.zeros(len(RADIAL_INTERVALS)*len(AXIAL_INTERVALS))
    for lo, hi in zip(edges[:-1], edges[1:]):
        rr = (lo+hi)/2+(hi-lo)/2*g
        weights = (hi-lo)/2*w
        points = np.column_stack((rr, np.zeros(order),
                                  np.full(order, z)))
        _, gradients = patch.basis(points, tau)
        axial -= (weights*rr)@gradients[:, 2, :]/radius
    return np.vstack((np.zeros_like(axial), axial))


def cone_row(base, patch, X, eta, tau, order=64):
    point = base.compact.joined.inner.from_similarity(
        np.array([X]), np.array([eta]), tau)
    radius, _, z = point[0]
    u, grad, R = operator(base, point, tau)
    target = stress_primitive(base, radius, z, tau, order=order)
    columns = pressure_stress_columns(patch, point[0], tau, order)
    _, gradients = patch.basis(point, tau)
    F = u[0, 1]/radius
    shear = np.array([grad[0, 1, 0]-F, grad[0, 2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    c = (float(np.sqrt(lam2)/(2*F*N[0]))
         if lam2 > 0 and abs(2*F*N[0]) > 1e-14 else None)
    return dict(X=X, eta=eta, tau=tau, point=point[0],
                target=target, columns=columns, N=N, K=K,
                lambda_squared=float(lam2), c=c,
                residual=R[0], gradients=gradients[0])


def cone_metrics(row, coefficients):
    target = row['target']+row['columns']@coefficients
    dot_n = float(row['N']@target)
    dot_k = float(row['K']@target)
    ratio = (abs(row['c']*dot_k/dot_n)
             if row['c'] is not None and abs(dot_n) > 1e-14 else None)
    return dict(X=row['X'], eta=row['eta'], tau=row['tau'],
                lambda_squared=row['lambda_squared'],
                target_dot_N=dot_n, target_dot_K=dot_k,
                ratio=ratio,
                strict_pass=bool(row['lambda_squared'] > 0 and dot_n < 0
                                 and ratio is not None and ratio < 1),
                residual_norm=float(np.linalg.norm(
                    row['residual']+row['gradients']@coefficients)))


def free_axial_feasible(row):
    """Check whether arbitrary axial stress could enter the cone locally."""
    lower, upper = -np.inf, np.inf
    for sign in (-1., 1.):
        direction = CONE_RATIO_MARGIN*row['N']+sign*row['c']*row['K']
        rhs = -STRESS_MARGIN-direction[0]*row['target'][0]
        if direction[1] > 1e-14:
            upper = min(upper, rhs/direction[1])
        elif direction[1] < -1e-14:
            lower = max(lower, rhs/direction[1])
        elif rhs < 0:
            return False
    return bool(lower <= upper)


def run():
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    patch = SimilarityPressurePatch(base, np.zeros(6))
    tau = .5*2**(-5.5)
    train = [cone_row(base, patch, X, eta, tau)
             for eta in TRAIN_ETA for X in TRAIN_X]
    if any(row['c'] is None for row in train):
        raise RuntimeError('A training point has nonpositive lambda squared')
    A, b = [], []
    for row in train:
        for sign in (-1., 1.):
            direction = CONE_RATIO_MARGIN*row['N']+sign*row['c']*row['K']
            A.append(direction@row['columns'])
            b.append(-STRESS_MARGIN-direction@row['target'])
    A, b = np.asarray(A), np.asarray(b)
    scale = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
    feasible = linprog(np.zeros(6), A_ub=A*scale, b_ub=b,
                       bounds=[(None, None)]*6, method='highs')
    single_node = [bool(linprog(
        np.zeros(6), A_ub=A[2*i:2*i+2]*scale,
        b_ub=b[2*i:2*i+2], bounds=[(None, None)]*6,
        method='highs').success) for i in range(len(train))]
    prefix = []
    for count in range(1, len(TRAIN_X)+1):
        selected = [i for i, row in enumerate(train)
                    if row['X'] in TRAIN_X[:count]]
        inequalities = [2*i+j for i in selected for j in (0, 1)]
        fit_prefix = linprog(np.zeros(6),
                             A_ub=A[inequalities]*scale,
                             b_ub=b[inequalities],
                             bounds=[(None, None)]*6, method='highs')
        prefix.append(dict(max_X=TRAIN_X[count-1],
                           feasible=bool(fit_prefix.success)))
    relaxed_A, relaxed_b = [], []
    for row in train:
        for sign in (-1., 1.):
            direction = .999*row['N']+sign*row['c']*row['K']
            relaxed_A.append(direction@row['columns'])
            relaxed_b.append(-direction@row['target'])
    relaxed = linprog(np.zeros(6), A_ub=np.asarray(relaxed_A)*scale,
                      b_ub=np.asarray(relaxed_b),
                      bounds=[(None, None)]*6, method='highs')
    report = dict(base_slice=BASE_NAME, tau=tau,
                  train_x=TRAIN_X, train_eta=TRAIN_ETA,
                  radial_intervals=RADIAL_INTERVALS,
                  axial_intervals=AXIAL_INTERVALS,
                  cone_ratio_margin=CONE_RATIO_MARGIN,
                  stress_margin=STRESS_MARGIN,
                  linear_feasible=bool(feasible.success),
                  linear_message=feasible.message,
                  single_node_feasible=single_node,
                  free_axial_feasible=[free_axial_feasible(row)
                                       for row in train],
                  x_prefix_feasible=prefix,
                  relaxed_0999_feasible=bool(relaxed.success),
                  raw_geometry=[dict(X=row['X'], eta=row['eta'],
                                     target=row['target'].tolist(),
                                     pressure_columns=row['columns'].tolist(),
                                     N=row['N'].tolist(), K=row['K'].tolist(),
                                     c=row['c']) for row in train],
                  baseline=[cone_metrics(row, np.zeros(6)) for row in train],
                  scope='Compact pressure-only physical-cone feasibility. '
                        'Velocity and five moments unchanged. No supported '
                        'wave, pressure-Poisson matching, volume norm, or '
                        'PDE acceptance.', accepted=False)
    if feasible.success:
        R0 = np.stack([row['residual'] for row in train])
        G = np.stack([row['gradients'] for row in train])
        normalizer = np.linalg.norm(R0)

        def objective(v):
            c = scale*v
            R = R0+G@c
            return .5*(np.sum(R**2)/normalizer**2+1e-6*np.sum(v**2))

        def jac(v):
            c = scale*v
            R = R0+G@c
            return scale*np.einsum('nik,ni->k', G, R)/normalizer**2+1e-6*v

        fit = minimize(objective, feasible.x, jac=jac, method='SLSQP',
                       constraints=[dict(type='ineq',
                                         fun=lambda v: b-A@(scale*v),
                                         jac=lambda v: -A*scale)],
                       options=dict(maxiter=1000, ftol=1e-12))
        v = (fit.x if np.max(A@(scale*fit.x)-b) < 1e-6
             else feasible.x)
        coefficients = scale*v
        patch.coefficients = coefficients
        train_after = [cone_metrics(row, coefficients) for row in train]
        holdout_tau = .5*2**(-5.4)
        holdout = [cone_row(base, patch, X, eta, holdout_tau)
                   for eta in HOLDOUT_ETA for X in HOLDOUT_X]
        points = np.stack([row['point'] for row in holdout])
        R_fd, _ = residual(patch, points, holdout_tau)
        R_linear = np.stack([row['residual']+row['gradients']@coefficients
                             for row in holdout])
        report.update(coefficients=coefficients.tolist(),
                      optimizer_success=bool(fit.success),
                      optimizer_message=fit.message,
                      train_after=train_after,
                      train_momentum=dict(before=stats(R0),
                                          after=stats(R0+G@coefficients)),
                      holdout_before=[cone_metrics(row, np.zeros(6))
                                      for row in holdout],
                      holdout_after=[cone_metrics(row, coefficients)
                                     for row in holdout],
                      holdout_momentum=dict(before=stats(np.stack(
                          [row['residual'] for row in holdout])),
                                            after=stats(R_linear)),
                      analytic_vs_fd_max=float(np.max(np.abs(
                          R_fd-R_linear))))
    (ROOT/'delayed_pressure_cone_window.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report.get(key) for key in (
        'linear_feasible', 'linear_message', 'coefficients',
        'train_momentum', 'holdout_momentum', 'analytic_vs_fd_max')}),
        flush=True)


if __name__ == '__main__':
    run()
