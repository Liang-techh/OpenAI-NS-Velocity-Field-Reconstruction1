"""Locate physical-volume momentum error before choosing new basis support."""
import json
import argparse

import numpy as np

from joined_field import ROOT, coordinates
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field


def run(order=8):
    field = current_field()
    tau = .5/64
    points, weights = nodes(field, tau, order)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, J, part = kinematics(field, points, tau, hs, ht)
    R = part+np.einsum('nij,nj->ni', J, u)
    norms = np.linalg.norm(R, axis=1)
    contributions = weights*norms**2
    ranking = np.argsort(contributions)[::-1]
    compact = field.base
    rflat, rsupp, zflat, zsupp = compact.support(tau)
    rows = []
    for i in ranking[:20]:
        r = float(np.hypot(points[i, 0], points[i, 1]))
        z = float(points[i, 2])
        co = coordinates(r/np.sqrt(field.nu), z/np.sqrt(field.nu),
                         tau, compact.base.inner.h)
        rows.append({'r': r, 'z': z, 'X': float(co['X']),
                     'eta': float(co['eta']),
                     'residual': R[i].tolist(), 'norm': float(norms[i]),
                     'volume_l2_squared_contribution': float(contributions[i]),
                     'fraction_of_total_l2_squared': float(contributions[i]/sum(contributions))})
    report = {'tau': tau, 'order': order, 'support': [rflat, rsupp, zflat, zsupp],
              'physical_volume_l2': float(np.sqrt(sum(contributions))),
              'sample_max': float(max(norms)), 'top_rows': rows,
              'scope': 'Gauss sample of full compact support at one time; hotspot ranking is quadrature-dependent, not a uniform bound.'}
    out = ROOT/'compact_potential'/f'transition_hotspot_map_gauss{order}.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'physical_volume_l2': report['physical_volume_l2'],
                      'sample_max': report['sample_max'],
                      'top': rows[:8]}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--order', type=int, default=8)
    run(parser.parse_args().order)
