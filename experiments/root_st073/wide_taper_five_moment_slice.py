"""Construct exact fixed-slice moments with a wider smooth U taper."""
import json

from coupled_five_moment_slice import construct
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep
from taper_width_capacity import grid


def run():
    tau = .5*2**(-5.5)
    width, degree = .05, 19
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid(width)
    prior = json.loads((ROOT/'taper_width_capacity.json').read_text())
    rows = []
    for item in prior['rows']:
        if item['width'] != width:
            continue
        row = dict(variant='maximum_slack', **construct(
            base, changed, X, weights, item['eta'], tau,
            item['coefficients'], width=width, degree=degree))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, quadrature_per_piece=96,
                  taper_width=width, u_degree=degree, rows=rows,
                  scope='Two fixed-eta normalized five-moment repairs with a wider 0.05 U taper and degree-19 basis. No smooth eta/time physical lift, pressure adjustment, continuous cone or complete momentum acceptance.',
                  accepted=False)
    (ROOT/'wide_taper_five_moment_slice.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
