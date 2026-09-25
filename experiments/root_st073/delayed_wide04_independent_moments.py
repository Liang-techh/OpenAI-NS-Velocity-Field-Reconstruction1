"""Independently integrate five moments of the actual lifted velocity."""
import argparse
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_taper_capacity_screen import grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


def run(slice_filename='delayed005_wide04_curvature_optimize.json',
        output_name='delayed_wide04_independent_moments.json'):
    tau = .5*2**(-5.5)
    X, weights = grid(order=48)
    target = make_field(16, 2.)
    lift = CoupledMomentPhysicalLift(
        slice_filename=slice_filename)
    rows = []
    for eta in (.2, .3):
        U0, E0 = profile(target, X, eta, tau)
        U, E = profile(lift, X, eta, tau)
        defect = moment_vector(U, E, X, weights)-moment_vector(
            U0, E0, X, weights)
        row = dict(eta=eta, five_moment_defect=defect.tolist(),
                   max_abs_defect=float(np.max(np.abs(defect))),
                   min_E=float(np.min(E)),
                   min_relative_E=float(np.min(E/np.maximum(E0, 1e-300))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, quadrature_per_piece=48,
                  target_field='make_field(16, 2.)',
                  physical_field='CoupledMomentPhysicalLift '+slice_filename,
                  rows=rows,
                  scope='Independent piecewise Gauss48 moments of actual '
                        'physical velocity at two fixed axial slices and '
                        'one time; no continuous eta/time identity.',
                  accepted=False)
    (ROOT/output_name).write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slice-filename',
                        default='delayed005_wide04_curvature_optimize.json')
    parser.add_argument('--output-name',
                        default='delayed_wide04_independent_moments.json')
    arguments = parser.parse_args()
    run(arguments.slice_filename, arguments.output_name)
