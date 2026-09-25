"""Re-optimize E capacity for a U repair beginning just outside X=1.03.

This is a fixed-slice necessary-condition search. The U correction is
arbitrary L2 on the proposed interval and is not constructed here.
"""

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
    rng = np.random.default_rng(73074)
    rows = [positive.solve_eta(mean, target, base, X, weights, eta, rng)
            for eta in (.275, .3)]
    report = dict(source='delayed_remote_positive_e.json',
                  U_support=positive.U_SUPPORT, rows=rows,
                  scope='Independent fixed-eta I/Cp-positive E fits '
                        'maximizing optimistic S slack for arbitrary '
                        'U on X=1.03..3.5. No smooth U mode, eta lift, '
                        'pressure, or full momentum acceptance.',
                  accepted=False)
    (ROOT/'delayed_remote_e_wide_support_probe.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(row['eta'], row['feasible'],
                       row['optimistic_continuum_S_slack']['S_slack'],
                       row['positivity']['min_relative_E'],
                       row['I_Cp_defect']) for row in rows]), flush=True)


if __name__ == '__main__':
    run()
