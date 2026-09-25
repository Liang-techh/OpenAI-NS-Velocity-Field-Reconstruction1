"""Test nonlinear moment reachability before imposing the swirl floor.

The five-eta, forty-mode Jacobian is full rank, but this alone says
nothing about the nonlinear equations or positive swirl. This solve
removes the soft floor and coefficient regularizer to separate those
questions. The output is diagnostic, never an accepted NS candidate.
"""

import json

import numpy as np
from scipy.optimize import least_squares

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import (HOLDOUT_ETAS, MOMENT_SCALES,
                                          NEAR5_AXIAL_INTERVALS,
                                          NEAR_RADIAL_INTERVALS, TRAIN_ETAS,
                                          current_mean, defects, make_rows,
                                          quadrature, remote_modes)
from delayed_similarity_curl_screen import CurlPatchedLift
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


def moment_jacobian(row, coefficients):
    _, U, E = defects(row, coefficients)
    Bu, Be = row['Bu'], row['Be']
    X, w = row['X'], row['weights']
    H = np.sqrt(2*X)
    return np.stack(((w*H)@Be,
                     (w*H*E)@Bu+(w*H*U)@Be,
                     (2*w*U)@Bu-(w*E)@Be,
                     (w*E/X)@Be))/MOMENT_SCALES[:, None]


def run():
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    modes, kinds = remote_modes(base, NEAR_RADIAL_INTERVALS,
                                NEAR5_AXIAL_INTERVALS)
    X, weights = quadrature(radial_intervals=NEAR_RADIAL_INTERVALS)
    tau = .5*2**(-5.5)
    training = make_rows(mean, target, base, modes, X, weights,
                         TRAIN_ETAS, tau)
    prior = json.loads((ROOT/'delayed_remote_moment_repair_near5wide.json').read_text())
    start = np.asarray(prior['coefficients'])
    lower = np.array([-.5 if kind == 'swirl' else -.1 for kind in kinds])
    upper = -lower

    def objective(c):
        return np.array([defects(row, c)[0][1:]/MOMENT_SCALES
                         for row in training]).ravel()

    def jacobian(c):
        return np.vstack([moment_jacobian(row, c) for row in training])

    fit = least_squares(objective, start, jac=jacobian,
                        bounds=(lower, upper), max_nfev=1000,
                        ftol=1e-12, xtol=1e-12, gtol=1e-12)
    coefficients = fit.x
    holdout = make_rows(mean, target, base, modes, X, weights,
                        HOLDOUT_ETAS, tau)

    def summarize(rows):
        return [dict(eta=row['eta'],
                     defect=defects(row, coefficients)[0].tolist(),
                     min_relative_E=float(np.min(
                         defects(row, coefficients)[2]/row['E0'])))
                for row in rows]

    patched = CurlPatchedLift(mean, modes, coefficients)
    nearby_tau = .5*2**(-5.4)
    remote_points, _, _ = nodes(base, (1.2, 1.5, 2., 2.7),
                                (.225, .275), nearby_tau)
    before, _ = residual(mean, remote_points, nearby_tau)
    after, divergence = residual(patched, remote_points, nearby_tau)
    cone_points, _, _ = nodes(base, (1.008, 1.016, 1.02, 1.03),
                              (.2, .25, .3), tau)
    old_local, _ = mean.fields(cone_points, tau)
    new_local, _ = patched.fields(cone_points, tau)
    report = dict(source='delayed_remote_moment_repair_near5wide.json',
                  moment_order=['M', 'I', 'J', 'S', 'Cp'],
                  modes=len(modes), coefficients=coefficients.tolist(),
                  optimizer_success=bool(fit.success),
                  optimizer_message=str(fit.message),
                  function_evaluations=int(fit.nfev),
                  scaled_training_l2=float(np.linalg.norm(fit.fun)),
                  training=summarize(training), holdout=summarize(holdout),
                  max_abs_local_velocity_change=float(np.max(np.abs(
                      new_local-old_local))),
                  remote_momentum=dict(before=stats(before),
                                       after=stats(after),
                                       max_abs_fd_divergence=float(
                                           np.max(np.abs(divergence)))),
                  scope='Unregularized bounded forty-mode nonlinear '
                        'moment reachability at five eta slices. No '
                        'positivity constraint, continuous closure, '
                        'or PDE acceptance.', accepted=False)
    (ROOT/'delayed_remote_moment_feasibility.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(success=report['optimizer_success'],
                          message=report['optimizer_message'],
                          evaluations=report['function_evaluations'],
                          scaled_l2=report['scaled_training_l2'],
                          training_max=max(abs(x) for row in report['training']
                                           for x in row['defect']),
                          holdout_max=max(abs(x) for row in report['holdout']
                                          for x in row['defect']),
                          min_relative_E=min(row['min_relative_E']
                                             for row in report['training']
                                             +report['holdout']),
                          remote_momentum=report['remote_momentum'])),
          flush=True)


if __name__ == '__main__':
    run()
