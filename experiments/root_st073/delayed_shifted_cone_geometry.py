"""Locate the pointwise obstruction in the shifted-swirl 15-node cone."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_outer_radial_extension import make_extended_modes
from delayed_pressure_cone_geometry import line_interval
from delayed_shifted_swirl_cone_fit import (AMPLITUDE, NEW_AXIAL,
                                             NEW_RADIAL)
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_similarity_curl_screen import CurlPatchedLift
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_shifted_swirl_cone_fit.json').read_text())
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mode = SimilaritySwirlMode(base, NEW_RADIAL,
                               eta_interval=NEW_AXIAL)
    mean = CurlPatchedLift(base, make_extended_modes(base)+[mode],
                           [*mean_source['coefficients'], AMPLITUDE])
    rows = source['geometry']
    tau = source['tau']
    points = base.compact.joined.inner.from_similarity(
        np.array([row['X'] for row in rows]),
        np.array([row['eta'] for row in rows]), tau)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u, grad, _ = kinematics(mean, points, tau, hs, ht)
    diagnostics = []
    for i, row in enumerate(rows):
        target = np.asarray(row['target'])
        F = u[i, 1]/points[i, 0]
        shear = np.array([grad[i, 1, 0]-F, grad[i, 2, 0]])
        N = shear/np.linalg.norm(shear)
        K = np.array([-N[1], N[0]])
        lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        multiplier = (np.sqrt(lam2)/(2*F*N[0])
                      if lam2 > 0 and abs(2*F*N[0]) > 1e-14 else None)
        if multiplier is None:
            axial = None
            theta = None
        else:
            directions = [.8*N+sign*multiplier*K
                          for sign in (-1., 1.)]
            axial = line_interval(directions, target, .1, 1)
            theta = line_interval(directions, target, .1, 0)
        diagnostics.append(dict(X=row['X'], eta=row['eta'],
                                lambda_squared=float(lam2),
                                target=target.tolist(),
                                free_axial=axial is not None,
                                free_theta=theta is not None,
                                axial_interval=(list(axial)
                                                if axial is not None else None)))
    report = dict(source='delayed_shifted_swirl_cone_fit.json',
                  rows=diagnostics,
                  scope='Pointwise one-stress-component necessary cone '
                        'conditions on fifteen sampled nodes. No '
                        'realizable pressure, continuous cone, wave, '
                        'or PDE acceptance.', accepted=False)
    (ROOT/'delayed_shifted_cone_geometry.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(free_axial=sum(row['free_axial']
                                         for row in diagnostics),
                          free_theta=sum(row['free_theta']
                                         for row in diagnostics),
                          failures=[(row['X'], row['eta']) for row in diagnostics
                                    if not row['free_axial']])), flush=True)


if __name__ == '__main__':
    run()
