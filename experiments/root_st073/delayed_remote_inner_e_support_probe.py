"""Test whether E support approaching X=1.03 opens last-slice S capacity."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
import delayed_remote_positive_e as positive
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


def run():
    positive.U_SUPPORT = (1.03, 3.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    original = positive.RADIAL_INTERVALS
    rows = []
    for first_start in (1.04, 1.035, 1.031, 1.0301):
        positive.RADIAL_INTERVALS = ((first_start, original[0][1]),
                                     *original[1:])
        X, weights = positive.quadrature()
        rng = np.random.default_rng(73074)
        row = positive.solve_eta(mean, target, base, X, weights, .3, rng)
        rows.append(dict(first_E_start=first_start,
                         radial_intervals=positive.RADIAL_INTERVALS,
                         result=row))
    report = dict(source='delayed_remote_positive_e.json',
                  U_support=positive.U_SUPPORT, rows=rows,
                  scope='Fixed eta=.3 optimistic S capacity as first E '
                        'support approaches the protected cone edge. '
                        'No smooth U repair, full momentum, or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_remote_inner_e_support_probe.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(r['first_E_start'], r['result']['feasible'],
                       r['result']['optimistic_continuum_S_slack']['S_slack'],
                       r['result']['positivity']['min_relative_E'])
                      for r in rows]), flush=True)


if __name__ == '__main__':
    run()
