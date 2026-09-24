"""Check whether the fitted pressure cone patch persists with fixed amplitudes."""
import json

import numpy as np

from curl_wave_cone_region import evaluate_height
from joined_field import ROOT
from radial_pressure_patch_screen import load_dense_candidate


def run():
    candidate = load_dense_candidate()
    times = [.00825, .0084, .00855]
    radii = [.0075, .00825, .009, .010, .011, .012, .013]
    heights = [.00355, .003625, .0037]
    rows = []
    grid = []
    for tau in times:
        _, r_support, _, z_support = candidate.support(tau)
        slab = []
        for z in heights:
            line = evaluate_height(candidate, np.asarray(radii), z, tau,
                                   r_support, z_support, quadrature_order=8)
            rows.extend([{'tau': tau, **row} for row in line])
            slab.append([bool(row['strict_local_pass'] and
                              row.get('cone_ratio', np.inf) < .8)
                         for row in line])
            print(json.dumps({'tau': tau, 'z': z,
                              'passing': sum(slab[-1])}), flush=True)
        grid.append(slab)
    passing = np.asarray(grid)
    intervals = []
    for ti in range(len(times)-1):
        for tj in range(ti+1, len(times)):
            for zi in range(len(heights)-1):
                for zj in range(zi+1, len(heights)):
                    for ri in range(len(radii)-1):
                        for rj in range(ri+1, len(radii)):
                            if np.all(passing[ti:tj+1, zi:zj+1, ri:rj+1]):
                                intervals.append({
                                    'time_interval': [times[ti], times[tj]],
                                    'r_interval': [radii[ri], radii[rj]],
                                    'z_interval': [heights[zi], heights[zj]],
                                    'node_count': int((tj-ti+1)*(zj-zi+1)*(rj-ri+1)),
                                    'spacetime_box_volume': ((times[tj]-times[ti])
                                                             *(radii[rj]-radii[ri])
                                                             *(heights[zj]-heights[zi]))})
    intervals.sort(key=lambda row: row['spacetime_box_volume'], reverse=True)
    report = {'times': times, 'radii': radii, 'heights': heights,
              'passing_count': int(np.sum(passing)),
              'total_count': int(passing.size),
              'pass_grid': passing.tolist(),
              'largest_passing_sample_boxes': intervals[:20],
              'rows': rows,
              'scope': 'Fixed candidate at 3 times; finite sampled cone boxes only. No continuum or full slab certificate, no wave construction.',
              'accepted': False}
    path = ROOT/'compact_potential'/'radial_pressure_time_screen.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'passing_count': report['passing_count'],
                      'total_count': report['total_count'],
                      'largest': intervals[:2]}), flush=True)


if __name__ == '__main__':
    run()
