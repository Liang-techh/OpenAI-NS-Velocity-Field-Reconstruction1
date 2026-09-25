"""Search fixed-moment U tangents using entrance, return, and two times.

This is a sampled momentum screen, not a spacetime or cone certificate.
"""

import json

import numpy as np
from scipy.optimize import least_squares

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import (
    SLICE, fit_geometry, nodes, repaired_coefficients, residual,
)
from moment_shear_slice_repair import moment_vector
from coupled_five_moment_slice import correction_modes
from radial_continuation import ROOT


def evaluate(field, node_sets):
    return [residual(field, points, tau)[0] for points, tau in node_sets]


def norms_summary(arrays):
    norms = [np.linalg.norm(array, axis=1) for array in arrays]
    return dict(max=float(max(np.max(n) for n in norms)),
                rms=float(np.sqrt(np.mean(np.concatenate([n*n for n in norms])))),
                per_time_max=[float(np.max(n)) for n in norms])


def run():
    tau0 = .5*2**(-5.5)
    tau1 = .5*2**(-5.25)
    field = CoupledMomentPhysicalLift(slice_filename=SLICE)
    Xq, weights, Q, factor, U, E, c0, wq, directions = fit_geometry(
        field, tau0)
    xs = (1.01, 1.015, 1.0175, 1.02, 1.025, 1.05, 1.08,
          1.12, 2.6, 2.975)
    node_sets = [(nodes(field, xs, (.3,), tau)[0], tau)
                 for tau in (tau0, tau1)]
    baseline = evaluate(field, node_sets)
    slopes = []
    eps = .035
    for direction in directions:
        paired = []
        for sign in (-1., 1.):
            repaired = repaired_coefficients(
                sign*eps, direction['vq'], wq, Q, factor, U, weights, c0)
            if repaired is None:
                raise ValueError('A tangent finite difference left S manifold')
            field.u_rows[1] = repaired[0]
            paired.append(evaluate(field, node_sets))
        slopes.append([(pos-neg)/(2*eps)
                       for neg, pos in zip(paired[0], paired[1])])
        field.u_rows[1] = c0
    # A normalised least-squares surrogate balances the entrance and return
    # samples. The far-return samples are weighted by their local baseline
    # scale so a small absolute residual cannot deteriorate unnoticed.
    scales = np.maximum(np.abs(np.concatenate(baseline, axis=0)), 1e4)
    scales = np.linalg.norm(scales, axis=1)
    jacobian = np.column_stack([
        np.concatenate(slope, axis=0).reshape(-1) for slope in slopes])
    target = np.concatenate(baseline, axis=0).reshape(-1)
    weights_flat = np.repeat(1/scales, 3)
    bounds = np.array([direction['bound'] for direction in directions])
    fit = least_squares(lambda a: np.r_[
        weights_flat*(target+jacobian@a), .1*a/bounds],
        np.zeros(len(directions)), bounds=(-bounds, bounds),
        xtol=1e-11, ftol=1e-11, gtol=1e-11)
    trials = []
    for fraction in (.025, .05, .1, .15, .2, .25, .5, .75, 1.):
        a = fraction*fit.x
        vq = sum((float(ai)*direction['vq']
                  for ai, direction in zip(a, directions)),
                 np.zeros_like(wq))
        repaired = repaired_coefficients(
            1., vq, wq, Q, factor, U, weights, c0)
        if repaired is None:
            continue
        field.u_rows[1] = repaired[0]
        arrays = evaluate(field, node_sets)
        summary = norms_summary(arrays)
        old_norms = np.linalg.norm(np.concatenate(baseline, axis=0), axis=1)
        new_norms = np.linalg.norm(np.concatenate(arrays, axis=0), axis=1)
        summary['far_return_ratio'] = float(
            max(new_norms[9]/old_norms[9],
                new_norms[19]/old_norms[19]))
        summary['outer_max'] = float(max(
            np.max(new_norms[5:10]), np.max(new_norms[15:20])))
        summary['a'] = a.tolist()
        summary['fraction_of_linear_fit'] = fraction
        summary['b'] = float(repaired[1])
        summary['coefficients'] = repaired[0].tolist()
        summary['node_norms'] = new_norms.tolist()
        trials.append(summary)
        print(json.dumps({key: summary[key] for key in
                          ('a', 'max', 'rms', 'far_return_ratio',
                           'outer_max')}), flush=True)
    field.u_rows[1] = c0
    before = norms_summary(baseline)
    before_norms = np.linalg.norm(np.concatenate(baseline, axis=0), axis=1)
    before['outer_max'] = float(max(np.max(before_norms[5:10]),
                                    np.max(before_norms[15:20])))
    before['node_norms'] = before_norms.tolist()
    admissible = [trial for trial in trials
                  if trial['max'] < before['max']
                  and trial['outer_max'] <= before['outer_max']
                  and trial['far_return_ratio'] <= 1.]
    selected = min(admissible, key=lambda trial: trial['max']) if admissible else None
    if selected is not None:
        U_trial = U+(np.asarray(selected['coefficients'])-c0)@correction_modes(
            Xq, field.width, field.degree, field.start_X)
        defect = moment_vector(U_trial, E, Xq, weights)-moment_vector(
            U, E, Xq, weights)
        selected['slice_moment_defect'] = defect.tolist()
        if np.max(np.abs(defect)) > 1e-7:
            raise ValueError('Selected profile left five-moment manifold')
    report = dict(slice_filename=SLICE, tau=[tau0, tau1], X=xs,
                  eta=.3, baseline=before, linear_fit=fit.x.tolist(),
                  trials=trials, selected=selected,
                  scope='Two times, ten similarity nodes per time, full '
                        'Cartesian residual; fixed eta=.3, fixed E. '
                        'The criterion requires entrance, outer max, and '
                        'far return improvement. No continuous or volume '
                        'claim.', accepted=False)
    (ROOT/'delayed_multitime_momentum_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'baseline': before['max'],
                      'selected': None if selected is None else selected['max']}))


if __name__ == '__main__':
    run()
