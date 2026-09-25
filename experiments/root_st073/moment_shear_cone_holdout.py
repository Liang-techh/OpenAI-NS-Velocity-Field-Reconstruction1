"""Hold out nearby space/time points for the first strict physical cone node."""
import json

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    field = RadialMomentStep(make_field(16, 2.), -.5)
    rows = []
    for k, radii in ((5.5, (.95, 1., 1.05)),
                     (5.25, (1.,))):
        tau = .5*2**(-k)
        for eta in (.2, .3):
            for X in radii:
                row = dict(k=k,
                           **physical_cone(field, X, eta, tau, order=64))
                rows.append(row)
                print(json.dumps(row), flush=True)
    report = dict(coefficients=dict(shear=2., moment=-.5),
                  rows=rows,
                  pass_count=sum(row['strict_pass'] for row in rows),
                  scope='Eight physical full-residual cone nodes near one discovery point. A finite pass set is not an open-region or all-time certificate; no wave or full momentum acceptance.',
                  accepted=False)
    (ROOT/'moment_shear_cone_holdout.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
