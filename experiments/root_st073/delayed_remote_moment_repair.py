"""Repair outgoing moments away from the fifteen-node stress cone.

Twenty-four compact exact-solenoidal velocity modes are supported at
X>=1.1. Five axial collocation slices constrain the J/S/I/Cp moments;
interlaced slices and remote momentum nodes check interpolation cost.
This is a mean-field experiment, not a completed Navier--Stokes field.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_axial_split_swirl_response import AXIAL_INTERVALS as SPLIT_AXIAL
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_outer_radial_extension import make_extended_modes
from delayed_shifted_swirl_cone_fit import AMPLITUDE, NEW_AXIAL, NEW_RADIAL
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


RADIAL_INTERVALS = ((1.1, 1.5), (1.3, 2.),
                    (1.7, 2.5), (2.2, 2.98))
NEAR_RADIAL_INTERVALS = ((1.04, 1.12), (1.04, 1.3),
                         (1.1, 1.5), (1.5, 2.5))
AXIAL_INTERVALS = ((.16, .24), (.21, .29), (.26, .34))
NEAR4_AXIAL_INTERVALS = AXIAL_INTERVALS+((.235, .315),)
NEAR5_AXIAL_INTERVALS = NEAR4_AXIAL_INTERVALS+((.185, .265),)
TRAIN_ETAS = (.2, .225, .25, .275, .3)
HOLDOUT_ETAS = (.2125, .2375, .2625, .2875)
MOMENT_SCALES = np.array((.001, .01, .01, .001))
RELATIVE_SWIRL_FLOOR = .11


def quadrature(order=32, radial_intervals=RADIAL_INTERVALS):
    edges = sorted(set((0., 3/32, 1., 1.005, 1.01, 1.02,
                        1.03, 1.06, 1.08, 1.12, 1.5, 1.75,
                        3., 3.5,
                        *(edge for interval in radial_intervals
                          for edge in interval))))
    g, w = leggauss(order)
    X = np.concatenate([(lo+hi)/2+(hi-lo)*g/2
                        for lo, hi in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(hi-lo)*w/2
                              for lo, hi in zip(edges[:-1], edges[1:])])
    return X, weights


def current_mean(base):
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    tuned = json.loads((ROOT/'delayed_axial_split_pressure_fit_tuned.json').read_text())
    broad = SimilaritySwirlMode(base, NEW_RADIAL,
                                eta_interval=NEW_AXIAL)
    splits = [SimilaritySwirlMode(base, NEW_RADIAL,
                                  eta_interval=axial)
              for axial in SPLIT_AXIAL]
    modes = make_extended_modes(base)+[broad]+splits
    coefficients = [*source['coefficients'], AMPLITUDE,
                    *tuned['split_amplitudes']]
    return CurlPatchedLift(base, modes, coefficients)


def remote_modes(base, radial_intervals=RADIAL_INTERVALS,
                 axial_intervals=AXIAL_INTERVALS):
    modes, kinds = [], []
    for axial in axial_intervals:
        for radial in radial_intervals:
            modes.append(SimilaritySwirlMode(base, radial,
                                              eta_interval=axial))
            kinds.append('swirl')
        for radial in radial_intervals:
            modes.append(SimilarityCurlMode(base, radial,
                                             eta_interval=axial))
            kinds.append('poloidal')
    return modes, kinds


def make_rows(mean, target, base, modes, X, weights, etas, tau):
    rows = []
    for eta in etas:
        U, E = profile(mean, X, eta, tau)
        _, E0 = profile(target, X, eta, tau)
        U0, _ = profile(target, X, eta, tau)
        basis = [mode_profile(mode, base, X, eta, tau)
                 for mode in modes]
        Bu = np.column_stack([item[0] for item in basis])
        Be = np.column_stack([item[1] for item in basis])
        rows.append(dict(eta=eta, X=X, weights=weights,
                         U=U, E=E, E0=E0, Bu=Bu, Be=Be,
                         target=moment_vector(U0, E0, X, weights)))
    return rows


def defects(row, coefficients):
    U = row['U']+row['Bu']@coefficients
    E = row['E']+row['Be']@coefficients
    return (moment_vector(U, E, row['X'], row['weights'])
            -row['target']), U, E


def run(near=False, near4=False, near5=False, near5wide=False):
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    radial_intervals = (NEAR_RADIAL_INTERVALS
                        if near or near4 or near5 or near5wide
                        else RADIAL_INTERVALS)
    axial_intervals = (NEAR5_AXIAL_INTERVALS if near5 or near5wide else
                       NEAR4_AXIAL_INTERVALS if near4 else
                       AXIAL_INTERVALS)
    modes, kinds = remote_modes(base, radial_intervals,
                                axial_intervals)
    tau = .5*2**(-5.5)
    X, weights = quadrature(radial_intervals=radial_intervals)
    floor_weight = (3000. if near5wide else
                    1000. if near or near4 or near5 else 100.)
    training = make_rows(mean, target, base, modes,
                         X, weights, TRAIN_ETAS, tau)
    scale = np.array([.05 if kind == 'swirl' else .01
                      for kind in kinds])
    lower = np.array([-.2 if kind == 'swirl' else -.05
                      for kind in kinds])
    if near5wide:
        lower = np.array([-.5 if kind == 'swirl' else -.1
                          for kind in kinds])
    upper = -lower

    def objective(c):
        moment_residuals, floor_residuals = [], []
        for row in training:
            delta, _, E = defects(row, c)
            moment_residuals.extend(delta[1:]/MOMENT_SCALES)
            floor_residuals.extend(floor_weight*np.maximum(
                0., RELATIVE_SWIRL_FLOOR-E/row['E0']))
        return np.r_[moment_residuals, floor_residuals,
                     .01*c/scale]

    def jacobian(c):
        blocks, floor_blocks = [], []
        for row in training:
            _, U, E = defects(row, c)
            Bu, Be = row['Bu'], row['Be']
            X, w = row['X'], row['weights']
            H = np.sqrt(2*X)
            derivative = np.stack((
                (w*H)@Be,
                (w*H*E)@Bu+(w*H*U)@Be,
                (2*w*U)@Bu-(w*E)@Be,
                (w*E/X)@Be))
            blocks.append(derivative/MOMENT_SCALES[:, None])
            active = E/row['E0'] < RELATIVE_SWIRL_FLOOR
            floor_blocks.append(np.where(active[:, None],
                                          -floor_weight*Be/row['E0'][:, None],
                                          0.))
        return np.vstack((*blocks, *floor_blocks,
                          np.diag(.01/scale)))

    start = np.zeros(len(modes))
    if near4 or near5 or near5wide:
        previous_name = ('delayed_remote_moment_repair_near5.json' if near5wide else
                         'delayed_remote_moment_repair_near4.json' if near5
                         else 'delayed_remote_moment_repair_near.json')
        previous = json.loads((ROOT/previous_name).read_text())
        start[:len(previous['coefficients'])] = previous['coefficients']
    fit = least_squares(objective, start,
                        jac=jacobian, bounds=(lower, upper),
                        max_nfev=(2000 if near5wide else
                                  1500 if near5 else
                                  1000 if near4 else 500), ftol=1e-11,
                        xtol=1e-11, gtol=1e-11)
    c = fit.x
    holdout = make_rows(mean, target, base, modes,
                        X, weights, HOLDOUT_ETAS, tau)

    def summarize(rows):
        return [dict(eta=row['eta'],
                     before=defects(row, np.zeros(len(modes)))[0].tolist(),
                     after=defects(row, c)[0].tolist(),
                     min_relative_E=float(np.min(
                         defects(row, c)[2]/row['E0'])))
                for row in rows]

    patched = CurlPatchedLift(mean, modes, c)
    nearby_tau = .5*2**(-5.4)
    remote_points, Xh, etah = nodes(
        base, (1.2, 1.5, 2., 2.7), (.225, .275), nearby_tau)
    before, _ = residual(mean, remote_points, nearby_tau)
    after, divergence = residual(patched, remote_points, nearby_tau)
    cone_points, _, _ = nodes(
        base, (1.008, 1.016, 1.02, 1.03),
        (.2, .25, .3), tau)
    old_local, _ = mean.fields(cone_points, tau)
    new_local, _ = patched.fields(cone_points, tau)
    report = dict(source='delayed_axial_split_pressure_fit_tuned.json',
                  tau=tau, radial_intervals=radial_intervals,
                  axial_intervals=axial_intervals,
                  floor_weight=floor_weight,
                  coefficient_bounds=dict(lower=lower.tolist(),
                                          upper=upper.tolist()),
                  kinds=kinds, coefficients=c.tolist(),
                  optimizer_success=bool(fit.success),
                  optimizer_message=str(fit.message),
                  objective_norm=float(np.linalg.norm(fit.fun)),
                  training=summarize(training),
                  holdout=summarize(holdout),
                  max_abs_local_velocity_change=float(np.max(np.abs(
                      new_local-old_local))),
                  remote_momentum=dict(tau=nearby_tau,
                                       X=Xh.tolist(), eta=etah.tolist(),
                                       before=stats(before),
                                       after=stats(after),
                                       max_abs_fd_divergence=float(
                                           np.max(np.abs(divergence)))),
                  scope='Five-eta compact remote exact-solenoidal mean '
                        'moment fit with interlaced eta and physical '
                        'remote momentum audit. It preserves the local '
                        'sampled velocity where X<=1.03, but has no '
                        'continuous moment proof, global full momentum '
                        'or volume L2 gate, wave, or PDE acceptance.',
                  accepted=False)
    output = ('delayed_remote_moment_repair_near5wide.json' if near5wide else
              'delayed_remote_moment_repair_near5.json' if near5 else
              'delayed_remote_moment_repair_near4.json' if near4 else
              'delayed_remote_moment_repair_near.json' if near else
              'delayed_remote_moment_repair.json')
    (ROOT/output).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(optimizer_success=report['optimizer_success'],
                          objective_norm=report['objective_norm'],
                          training_max=max(abs(x) for row in report['training']
                                           for x in row['after']),
                          holdout_max=max(abs(x) for row in report['holdout']
                                          for x in row['after']),
                          min_relative_E=min(row['min_relative_E']
                                             for row in report['training']
                                             +report['holdout']),
                          max_abs_local_velocity_change=(
                              report['max_abs_local_velocity_change']),
                          remote_momentum=report['remote_momentum'])),
          flush=True)


if __name__ == '__main__':
    import sys
    run(near='--near' in sys.argv[1:],
        near4='--near4' in sys.argv[1:],
        near5='--near5' in sys.argv[1:],
        near5wide='--near5wide' in sys.argv[1:])
