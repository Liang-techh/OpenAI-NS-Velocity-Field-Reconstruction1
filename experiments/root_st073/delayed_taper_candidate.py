"""Build a delayed-start, curvature-optimized two-slice five-moment lift."""
import argparse
import json

from coupled_five_moment_slice import construct
from delayed_taper_capacity_screen import grid
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep
from wide_taper_curvature_optimize import optimize_slice


def run(start_X=1.01, width=.05, degree=31,
        output_name='delayed_taper_curvature_optimize.json'):
    tau = .5*2**(-5.5)
    X, weights = grid(order=64)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    source = json.loads((ROOT/'wide_taper_five_moment_slice.json').read_text())
    rows = []
    for prior in source['rows']:
        eta = prior['eta']
        raw = construct(base, changed, X, weights, eta, tau,
                        prior['e_coefficients'], width=width,
                        degree=degree, start_X=start_X)
        if not raw.get('five_moments_restored'):
            raise ValueError(f'Delayed five-moment repair failed at eta={eta}')
        optimized = optimize_slice(base, changed, X, weights,
                                   eta, tau, raw)
        if not optimized['five_moments_restored']:
            raise ValueError(f'Curvature search lost moments at eta={eta}')
        optimized['initial_S_slack'] = raw['S_slack']
        rows.append(optimized)
        print(json.dumps(dict(eta=eta, S_slack=raw['S_slack'],
                              five_moment_defect=optimized[
                                  'max_abs_five_moment_defect'],
                              curvature_ratio=optimized['curvature_ratio'],
                              optimizer_success=optimized[
                                  'optimizer_success'])), flush=True)
    report = dict(tau=tau, taper_width=width, u_degree=degree,
                  start_X=start_X, quadrature_per_piece=64,
                  rows=rows,
                  scope='Two fixed-eta exact five-moment slices with U '
                        'correction delayed until X=1.01 and radial '
                        'curvature proxy optimized. A physical lift and '
                        'continuous cone/momentum checks are separate.',
                  accepted=False)
    (ROOT/output_name).write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--start-X', type=float, default=1.01)
    parser.add_argument('--width', type=float, default=.05)
    parser.add_argument('--degree', type=int, default=31)
    parser.add_argument('--output-name',
                        default='delayed_taper_curvature_optimize.json')
    arguments = parser.parse_args()
    run(arguments.start_X, arguments.width, arguments.degree,
        arguments.output_name)
