"""Screen positive E floors against eta=.3 optimistic J/S capacity."""

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
    X, weights = positive.quadrature()
    rows = []
    for floor in (.11, .10, .08, .05, .01):
        positive.RELATIVE_E_FLOOR = floor
        rng = np.random.default_rng(73074)
        row = positive.solve_eta(mean, target, base, X, weights, .3, rng)
        rows.append(dict(relative_E_floor=floor, result=row))
    report = dict(source='delayed_remote_positive_e.json',
                  U_support=positive.U_SUPPORT, rows=rows,
                  scope='Fixed eta=.3 I/Cp-positive E optimization '
                        'versus sampled relative floor. Positive '
                        'optimistic S slack is necessary, not '
                        'sufficient, for a smooth divergence-free '
                        'J/S correction or NS acceptance.',
                  accepted=False)
    (ROOT/'delayed_remote_swirl_floor_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(r['relative_E_floor'],r['result']['feasible'],
                       r['result']['optimistic_continuum_S_slack']['S_slack'],
                       r['result']['positivity']['min_relative_E'])
                      for r in rows]), flush=True)


if __name__ == '__main__':
    run()
