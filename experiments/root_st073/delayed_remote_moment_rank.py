"""Diagnose linearized reachability of five-slice remote moment repair."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import (MOMENT_SCALES,
                                           NEAR4_AXIAL_INTERVALS,
                                           NEAR5_AXIAL_INTERVALS,
                                           NEAR_RADIAL_INTERVALS,
                                           TRAIN_ETAS, current_mean,
                                           defects, make_rows,
                                           quadrature, remote_modes)
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


def run(near5=False):
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    axial_intervals = (NEAR5_AXIAL_INTERVALS if near5
                       else NEAR4_AXIAL_INTERVALS)
    modes, kinds = remote_modes(base, NEAR_RADIAL_INTERVALS,
                                axial_intervals)
    X, weights = quadrature(radial_intervals=NEAR_RADIAL_INTERVALS)
    tau = .5*2**(-5.5)
    rows = make_rows(mean, target, base, modes, X, weights,
                     TRAIN_ETAS, tau)
    matrix, baseline = [], []
    for row in rows:
        U, E = row['U'], row['E']
        Bu, Be = row['Bu'], row['Be']
        X, w = row['X'], row['weights']
        H = np.sqrt(2*X)
        derivative = np.stack((
            (w*H)@Be,
            (w*H*E)@Bu+(w*H*U)@Be,
            (2*w*U)@Bu-(w*E)@Be,
            (w*E/X)@Be))
        matrix.append(derivative/MOMENT_SCALES[:, None])
        baseline.extend(defects(row, np.zeros(len(modes)))[0][1:]
                        /MOMENT_SCALES)
    matrix = np.vstack(matrix)
    baseline = np.asarray(baseline)
    left_vectors, singular_values, _ = np.linalg.svd(
        matrix, full_matrices=False)
    missing_direction = left_vectors[:, -1]
    step = np.linalg.lstsq(matrix, -baseline, rcond=None)[0]
    predicted = baseline+matrix@step
    actual = np.array([defects(row, step)[0][1:]/MOMENT_SCALES
                       for row in rows]).ravel()
    relative_floor = min(float(np.min(
        defects(row, step)[2]/row['E0'])) for row in rows)
    report = dict(source=('proposed-near5' if near5 else
                          'delayed_remote_moment_repair_near4.json'),
                  axial_intervals=axial_intervals,
                  matrix_shape=matrix.shape,
                  rank=int(np.linalg.matrix_rank(matrix)),
                  singular_values=singular_values.tolist(),
                  missing_left_direction=(
                      missing_direction.reshape(len(TRAIN_ETAS), 4).tolist()),
                  baseline_missing_projection=float(
                      missing_direction@baseline),
                  baseline_scaled_l2=float(np.linalg.norm(baseline)),
                  linear_prediction_scaled_l2=float(np.linalg.norm(predicted)),
                  nonlinear_actual_scaled_l2=float(np.linalg.norm(actual)),
                  linear_solution_max_abs_coefficient=float(
                      np.max(np.abs(step))),
                  linear_solution_min_relative_E=relative_floor,
                  scope='Linearized rank and one unconstrained Newton '
                        'step for five-eta J/S/I/Cp moments. A full-rank '
                        'Jacobian is local reachability only; it does '
                        'not prove a bounded, positive-swirl nonlinear '
                        'solution or PDE acceptance.', accepted=False)
    output = ('delayed_remote_moment_rank_near5.json' if near5
              else 'delayed_remote_moment_rank.json')
    (ROOT/output).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(rank=report['rank'],
                          matrix_shape=report['matrix_shape'],
                          smallest_singular_values=(
                              singular_values[-5:].tolist()),
                          missing_left_direction=(
                              report['missing_left_direction']),
                          baseline_missing_projection=(
                              report['baseline_missing_projection']),
                          baseline_scaled_l2=(
                              report['baseline_scaled_l2']),
                          linear_prediction_scaled_l2=(
                              report['linear_prediction_scaled_l2']),
                          nonlinear_actual_scaled_l2=(
                              report['nonlinear_actual_scaled_l2']),
                          max_abs_coefficient=(
                              report['linear_solution_max_abs_coefficient']),
                          min_relative_E=(
                              report['linear_solution_min_relative_E']))),
          flush=True)


if __name__ == '__main__':
    import sys
    run(near5='--near5' in sys.argv[1:])
