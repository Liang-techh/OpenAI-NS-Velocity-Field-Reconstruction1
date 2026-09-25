"""Resolve the high-frequency physical stress primitive at X=1."""
import json

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    rows = []
    for label, amplitude, orders in (
            ('balanced', 0., (16, 32, 64)),
            ('shear_N16', 2., (16, 32, 64, 96))):
        field = make_field(16, amplitude)
        for order in orders:
            row = dict(field=label,
                       **physical_cone(field, 1., .2, tau, order=order))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Single-point Gauss-order convergence of the physical full-residual radial stress primitive. No spatially continuous stress cone or wave acceptance.',
                  accepted=False)
    (ROOT/'extended_cone_quadrature_audit.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
