"""Screen a solenoidal axial-collar mode against momentum hotspots."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import TransitionPoloidalMode, kinematics, nodes
from radial_pressure_volume_constrained import load_candidate


def compact_base(field):
    node = field
    while not isinstance(node, CompactPotentialField):
        node = node.base
    return node


def run():
    field = load_candidate()
    mode = TransitionPoloidalMode(compact_base(field), temporal_power=2.,
                                  reference_tau=.5/64)
    tau = .0084
    blocks = []
    for order in (6, 8, 10):
        points, weights = nodes(field, tau, order)
        hs = .0005*np.sqrt(field.nu*tau)
        ht = .0001*tau
        u, J, part = kinematics(field, points, tau, hs, ht)
        du, dJ, dpart = kinematics(mode, points, tau, hs, ht)
        base_residual = part+np.einsum('nij,nj->ni', J, u)
        linear = dpart+np.einsum('nij,nj->ni', dJ, u)+np.einsum('nij,nj->ni', J, du)
        quadratic = np.einsum('nij,nj->ni', dJ, du)
        blocks.append({'order': order, 'weights': weights,
                       'base_residual': base_residual,
                       'linear': linear, 'quadratic': quadratic,
                       'base_max_speed': float(np.max(np.linalg.norm(u, axis=1))),
                       'mode_max_speed': float(np.max(np.linalg.norm(du, axis=1)))})
        print(json.dumps({'loaded_order': order,
                          'base_max_speed': blocks[-1]['base_max_speed'],
                          'mode_max_speed': blocks[-1]['mode_max_speed']}), flush=True)
    speed_ratio = max(block['mode_max_speed']/block['base_max_speed']
                      for block in blocks)
    relative_strengths = [-1., -.5, -.2, -.1, -.05, -.02, 0.,
                          .02, .05, .1, .2, .5, 1.]
    rows = []
    for strength in relative_strengths:
        amplitude = strength/speed_ratio
        metrics = []
        for block in blocks:
            residual = (block['base_residual']+amplitude*block['linear']
                        +amplitude**2*block['quadratic'])
            norms = np.linalg.norm(residual, axis=1)
            metrics.append({'order': block['order'],
                            'max': float(np.max(norms)),
                            'physical_volume_l2': float(np.sqrt(block['weights']@norms**2))})
        rows.append({'relative_strength': strength,
                     'amplitude': amplitude, 'metrics': metrics})
    baseline = next(row for row in rows if row['relative_strength'] == 0.)
    for row in rows:
        row['worst_l2_ratio'] = max(
            metric['physical_volume_l2']/base['physical_volume_l2']
            for metric, base in zip(row['metrics'], baseline['metrics']))
        row['worst_max_ratio'] = max(
            metric['max']/base['max']
            for metric, base in zip(row['metrics'], baseline['metrics']))
    best = min(rows, key=lambda row: (row['worst_l2_ratio'], row['worst_max_ratio']))
    report = {'tau': tau, 'relative_strengths': relative_strengths,
              'baseline': baseline, 'best': best, 'rows': rows,
              'scope': 'One existing exact-solenoidal axial-collar mode, full quadratic momentum recombination at one time on Gauss6/8/10. No cone, time holdout, or continuum validation for nonzero amplitudes.',
              'accepted': False}
    path = ROOT/'compact_potential'/'poloidal_hotspot_screen.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'best': best, 'baseline': baseline}), flush=True)


if __name__ == '__main__':
    run()
