"""Direct intermediate-scale support corners and sampled outer moments."""
import json
import numpy as np

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from midplane_axial_cone_knob import cone_row
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run():
    inner, fields = build_fields()
    data = json.loads((ROOT / 'midplane_connected_cone_edge_repair.json').read_text())
    field = SeparatedMomentModes(fields['two_sided_cone'], data['amplitudes'],
                                 windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    orders = (12, 13, 14, 16, 17, 18)
    slices = moment_slices(inner, fields['two_sided_cone'], field,
                            orders=orders, unit_fields=[field])
    scales = []
    for k, moment_slice in zip(orders, slices):
        tau = .5 * 2.**-k
        X = inner.p.X_max * (1 + 15 * .325)**2
        point = inner.from_similarity([X], [-.0125], tau)[0]
        args = support(inner, dict(tau=tau, point=point), field.nu)
        rows = []
        for x in (-.95, .95):
            for z in (-.95, .95):
                r = point[0] + x * args['radial_halfwidth']
                height = point[2] + z * 1.4 * args['axial_halfwidth']
                row = cone_row(field, float(r), float(height), tau, z)
                rows.append(dict(radial_offset=x, axial_offset=z, **row))
        moments = outer_moments(moment_slice, np.zeros(1))
        scales.append(dict(k=k, rows=rows, moments=moments.tolist(),
                           moment_max=float(max(abs(moments)))))
        print(f"k={k}: {sum(r['cone_pass'] for r in rows)}/4 pass; "
              f"cone max={max(r.get('cone_ratio', float('inf')) for r in rows):.6g}; "
              f"moment max={max(abs(moments)):.6g}", flush=True)
    report = dict(scales=scales,
                  scope='Four physical support corners and direct outer moments at six integer holdout scales. No continuous scale/time/space certificate or residual contraction.')
    (ROOT / 'midplane_connected_cone_interscale.json').write_bytes((json.dumps(report, indent=2) + '\n').encode())


if __name__ == '__main__':
    run()
