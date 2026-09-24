"""Map local cone patches for the coupled annular candidate at interior time."""
import json

import numpy as np

from annular_pressure_scale_screen import load_candidate
from curl_wave_cone_region import evaluate_height
from joined_field import ROOT


def rectangles(radii, heights, passing):
    results = []
    for r0 in range(len(radii)):
        for r1 in range(r0+1, len(radii)):
            for z0 in range(len(heights)):
                for z1 in range(z0+1, len(heights)):
                    if np.all(passing[z0:z1+1, r0:r1+1]):
                        dr = float(radii[r1]-radii[r0])
                        dz = float(heights[z1]-heights[z0])
                        results.append({'r_interval': [float(radii[r0]), float(radii[r1])],
                                        'z_interval': [float(heights[z0]), float(heights[z1])],
                                        'radial_halfwidth': dr/2,
                                        'axial_halfwidth': dz/2,
                                        'grid_rectangle_area': dr*dz,
                                        'node_count': int((r1-r0+1)*(z1-z0+1))})
    results.sort(key=lambda row: row['grid_rectangle_area'], reverse=True)
    return results


def run():
    field = load_candidate()
    tau = .0084
    radii = np.array([.0025, .0035, .0045, .0055, .0065, .0075,
                      .009, .011, .013])
    heights = np.array([.0015, .002, .0024, .0027, .003,
                        .0032, .0034, .00355, .0037])
    _, r_support, _, z_support = field.support(tau)
    rows = []
    for z in heights:
        rows.extend(evaluate_height(field, radii, float(z), tau,
                                    r_support, z_support, quadrature_order=8))
        print(json.dumps({'loaded_z': float(z)}), flush=True)
    passing = np.array([[row['strict_local_pass'] and
                         row.get('cone_ratio', float('inf')) < .8
                         for row in rows[j*len(radii):(j+1)*len(radii)]]
                        for j in range(len(heights))])
    patches = rectangles(radii, heights, passing)
    report = {'tau': tau, 'radii': radii.tolist(),
              'heights': heights.tolist(),
              'passing_count': int(np.sum(passing)),
              'total_count': len(rows),
              'pass_grid': passing.tolist(),
              'largest_passing_rectangles': patches[:20],
              'rows': rows,
              'scope': 'Coarse local physical-field cone analogue at one time. Rectangles mean all listed grid nodes pass ratio<0.8; no continuum or spacetime cone certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_cone_area_map.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'passing_count': report['passing_count'],
                      'total_count': report['total_count'],
                      'largest': patches[:5]}), flush=True)


if __name__ == '__main__':
    run()
