"""Direct physical-residual check of the three-time pressure refit."""
import json
import sys

import numpy as np

from curl_wave_cone_region import evaluate_height
from joined_field import ROOT
from radial_pressure_patch_screen import (RadialPressureCandidate,
                                          load_dense_candidate)


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'radial_pressure_time_fit.json').read_text())
    fit = report['constant_fit']
    if not fit['feasible']:
        raise ValueError('No constant pressure fit is registered')
    return RadialPressureCandidate(load_dense_candidate(),
                                   report['mode_degrees'], fit['amplitudes'])


def run():
    optimized = '--volume-constrained' in sys.argv
    poloidal = '--local-poloidal' in sys.argv
    poloidal_10pct = '--local-poloidal-10pct' in sys.argv
    if poloidal:
        from local_poloidal_basis_screen import load_robust_candidate
        field = load_robust_candidate()
    elif poloidal_10pct:
        from local_poloidal_basis_screen import load_robust_candidate
        field = load_robust_candidate(.1)
    elif optimized:
        from radial_pressure_volume_constrained import load_candidate as load_optimized
        field = load_optimized()
    else:
        field = load_candidate()
    times = [.00825, .008325, .0084, .008475, .00855]
    radii = [.0075, .00825, .009, .010, .011, .012, .013]
    heights = [.00355, .003625, .0037]
    rows = []
    for tau in times:
        _, r_support, _, z_support = field.support(tau)
        for z in heights:
            line = evaluate_height(field, np.asarray(radii), z, tau,
                                   r_support, z_support, quadrature_order=8)
            rows.extend([{'tau': tau, **row} for row in line])
            passing = sum(row['strict_local_pass'] and
                          row.get('cone_ratio', np.inf) < .8 for row in line)
            print(json.dumps({'tau': tau, 'z': z,
                              'passing': passing, 'total': len(line)}), flush=True)
    report = {'times': times, 'radii': radii, 'heights': heights,
              'pass_count': sum(row['strict_local_pass'] and
                                row.get('cone_ratio', np.inf) < .8 for row in rows),
              'total_count': len(rows),
              'max_cone_ratio': max(row.get('cone_ratio', np.inf) for row in rows),
              'rows': rows,
              'scope': 'Physical finite-difference residual and radial primitive on 5 sampled times, 7 radii and 3 heights. No continuum or full momentum certificate.',
              'accepted': False}
    path = ROOT/'compact_potential'/(
        'local_poloidal_cone_validate.json' if poloidal else
        'local_poloidal_10pct_cone_validate.json' if poloidal_10pct else
        'radial_pressure_volume_cone_validate.json' if optimized else
        'radial_pressure_time_validate.json')
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'pass_count': report['pass_count'],
                      'total_count': report['total_count'],
                      'max_cone_ratio': report['max_cone_ratio']}), flush=True)


if __name__ == '__main__':
    run()
