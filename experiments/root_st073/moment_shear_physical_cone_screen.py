"""Couple radial M control and high-frequency shear at a physical point."""
import json

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    rows = []
    for shear in (-2., 0., 2., 4.):
        base = make_field(16, shear)
        for moment in (-1.5, -1., -.5, 0.):
            field = RadialMomentStep(base, moment)
            row = dict(shear=shear, moment=moment,
                       **physical_cone(field, 1., .2, tau, order=64))
            rows.append(row)
            print(json.dumps(row), flush=True)
    passing = [row for row in rows if row['strict_pass']]
    report = dict(tau=tau, rows=rows,
                  pass_count=len(passing),
                  scope='Sixteen-point parameter grid for one physical bridge point, with Gauss64 stress integration. Passing one point would not certify a continuous cone, moment match, or supported wave.',
                  accepted=False)
    (ROOT/'moment_shear_physical_cone_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(pass_count=len(passing),
                          best_by_projection=min(rows, key=lambda r: r['target_dot_N']))))


if __name__ == '__main__':
    run()
