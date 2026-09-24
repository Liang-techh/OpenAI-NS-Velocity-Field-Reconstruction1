"""Find positive stress representations using distinct integer angular modes."""
import itertools
import json

import numpy as np

from joined_field import ROOT


def run():
    directory = ROOT/'compact_potential'
    sources = sorted(directory.glob('annular_pressure_scale_cone*.json'))
    rows = []
    for path in sources:
        source = json.loads(path.read_text())
        target = np.asarray(source['local_tangential_stress_primitive'])
        pulses = source['kelvin_pulses']
        pairs = []
        for i, j in itertools.combinations(range(len(pulses)), 2):
            a, b = pulses[i], pulses[j]
            if a['mode'] == b['mode']:
                continue
            matrix = np.array([a['mean_radial_tangential_covariance'],
                               b['mean_radial_tangential_covariance']]).T
            if np.linalg.cond(matrix) > 1e8:
                continue
            weights = np.linalg.solve(matrix, target)
            if np.min(weights) <= 0:
                continue
            normals = [a['peak_wavevector'], b['peak_wavevector']]
            pairs.append({'indices': [i, j], 'modes': [a['mode'], b['mode']],
                          'weights': weights.tolist(),
                          'weight_sum': float(sum(weights)),
                          'peak_normal_norms': [float(np.linalg.norm(n))
                                                for n in normals],
                          'matrix_condition': float(np.linalg.cond(matrix))})
        pairs.sort(key=lambda item: item['weight_sum'])
        row = {'source': path.name, 'z': source['point'][2],
               'cone_pass': source['cone']['strict_pass'],
               'cone_ratio': source['cone'].get('ratio'),
               'unconstrained_selected_modes': [pulses[item['index']]['mode']
                                                for item in source['covariance_nnls']['selected']],
               'positive_distinct_mode_pair_count': len(pairs),
               'lowest_weight_pairs': pairs[:5]}
        rows.append(row)
        print(json.dumps({'source': row['source'],
                          'best_pair': pairs[0] if pairs else None}), flush=True)
    report = {'rows': rows,
              'scope': 'Local frozen Kelvin covariance at one point per source. Distinct angular modes avoid same-m cross terms in the angular mean; weights are local and do not solve a supported wave amplitude equation.',
              'accepted': False}
    out = directory/'wave_pair_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
