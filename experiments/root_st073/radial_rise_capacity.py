"""Screen the entrance streamfunction rise for cone and moment capacity."""

import argparse
import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS
from delayed_taper_capacity_screen import capacity, grid
from extended_physical_cone_map import physical_cone
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run(rise_ends=(1.1, 1.3, 1.5, 1.7, 1.9),
        output_name='radial_rise_capacity.json', cone_order=24):
    tau = .5*2**(-5.5)
    X, weights = grid(order=48)
    base = make_field(16, 2.)
    source = json.loads((ROOT/'delayed_e_capacity_optimize.json').read_text())
    e_coeff = {row['eta']: np.asarray(row['coefficients'])
               for row in source['rows']}
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    targets = {}
    for eta in (.2, .3):
        U0, E0 = profile(base, X, eta, tau)
        targets[eta] = dict(E=E0+e_coeff[eta]@eb,
                            target=moment_vector(U0, E0, X, weights))
    rows = []
    for rise_end in rise_ends:
        changed = RadialMomentStep(base, -.5, rise_end=rise_end)
        cone_rows = []
        moments = []
        for eta in (.2, .3):
            cone = physical_cone(changed, 1., eta, tau,
                                 order=cone_order)
            cone_rows.append(cone)
            U, _ = profile(changed, X, eta, tau)
            for degree in (19, 27, 31):
                result = capacity(X, weights, U, targets[eta]['E'],
                                  targets[eta]['target'], 1.005, .4,
                                  degree)
                moments.append(dict(eta=eta, **result))
        row = dict(rise_end=rise_end, cone_at_X1=cone_rows,
                   capacity=moments)
        rows.append(row)
        print(json.dumps(dict(rise_end=rise_end,
                              cone_pass=[c['strict_pass']
                                         for c in cone_rows],
                              cone_ratio=[c['ratio'] for c in cone_rows],
                              S_slack_eta03={m['degree']: m['S_slack']
                                             for m in moments
                                             if m['eta'] == .3})),
              flush=True)
    report = dict(tau=tau, amplitude=-.5,
                  restore_start=1.75, restore_end=3.,
                  E_source='delayed_e_capacity_optimize.json',
                  width=.4, start_X=1.005, quadrature_order=48,
                  cone_order=cone_order,
                  rows=rows, scope='One-point physical cone at X=1 and '
                        'fixed-slice finite-basis five-moment capacity '
                        'for alternate rise ends. No corrected physical '
                        'lift or continuous cone claim.', accepted=False)
    (ROOT/output_name).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rise-ends', default='1.1,1.3,1.5,1.7,1.9')
    parser.add_argument('--output-name', default='radial_rise_capacity.json')
    parser.add_argument('--cone-order', type=int, default=24)
    args = parser.parse_args()
    run(tuple(float(value) for value in args.rise_ends.split(',')),
        args.output_name, args.cone_order)
