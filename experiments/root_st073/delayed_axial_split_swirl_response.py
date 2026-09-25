"""Quadratic physical stress response of two axial-split outer swirl modes.

The low and high eta supports are disjoint and both start at X=1.02,
so they target the two newly failing X=1.03 slices while leaving the
original X<=1.02 window and eta=.25 center untouched at this time.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_multimode_cone_fit import BASE_NAME, cone
from delayed_multimode_pressure_admission import geometry
from delayed_outer_radial_extension import make_extended_modes
from delayed_pressure_cone_geometry import line_interval
from delayed_shifted_swirl_cone_fit import (AMPLITUDE, NEW_AXIAL,
                                             NEW_RADIAL)
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_swirl_cone_response import blocks, stress_primitive
from joint_collar_fit import kinematics
from radial_continuation import ROOT


RADIAL = (1.02, 1.06)
AXIAL_INTERVALS = ((.16, .24), (.26, .34))
ETAS = (.2, .3)
X = 1.03
CACHE_NAME = 'delayed_axial_split_swirl_response.npz'


def signature():
    return json.dumps(dict(base=BASE_NAME, radial=RADIAL,
                           axial=AXIAL_INTERVALS, X=X, etas=ETAS,
                           mean_source='delayed_outer_admission_search.json',
                           shifted_amplitude=AMPLITUDE,
                           order=32, version=1), sort_keys=True)


def build_rows():
    source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    broad = SimilaritySwirlMode(base, NEW_RADIAL,
                                eta_interval=NEW_AXIAL)
    mean = CurlPatchedLift(base, make_extended_modes(base)+[broad],
                           [*source['coefficients'], AMPLITUDE])
    modes = [SimilaritySwirlMode(base, RADIAL, eta_interval=axial)
             for axial in AXIAL_INTERVALS]
    tau = .5*2**(-5.5)
    edges = [1.005, 1.01, 1.02, 1.04, 1.05, 1.06, 1.08, 1.12]
    points, segments = blocks(base, tau, order=32, xs=(X,),
                              etas=ETAS, extra_radial_edges=edges)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u0, g0, p0 = kinematics(mean, points, tau, hs, ht)
    R0 = p0+np.einsum('pab,pb->pa', g0, u0)
    mode_data = [kinematics(mode, points, tau, hs, ht)
                 for mode in modes]
    U = np.stack([item[0] for item in mode_data], axis=2)
    G = np.stack([item[1] for item in mode_data], axis=3)
    P = np.stack([item[2] for item in mode_data], axis=2)
    L = (P+np.einsum('pab,pbi->pai', g0, U)
         +np.einsum('pabi,pb->pai', G, u0))
    rows = []
    for block in segments:
        sl = slice(block['start'], block['stop'])
        rr, weights, radius = (block[key] for key in
                               ('rr', 'weights', 'radius'))
        quadratic = np.einsum('pabi,pbj->paij', G[sl], U[sl])
        qtheta = -np.einsum('p,pij->ij', weights*rr**2,
                            quadratic[:, 1])/radius**2
        qaxial = -np.einsum('p,pij->ij', weights*rr,
                            quadratic[:, 2])/radius
        Tlin = np.stack([stress_primitive(L[:, :, i], block)
                         for i in range(len(modes))], axis=1)
        i = block['center']
        rows.append(dict(X=block['X'], eta=block['eta'],
                         radius=radius, u0=u0[i], g0=g0[i],
                         R0=R0[i], U=U[i], G=G[i], L=L[i],
                         Q=np.einsum('abi,bj->aij', G[i], U[i]),
                         T0=stress_primitive(R0, block),
                         Tlin=Tlin,
                         Tquad=np.stack((qtheta, qaxial))))
    return rows


def save_rows(rows):
    keys = ('u0', 'g0', 'R0', 'U', 'G', 'L', 'Q',
            'T0', 'Tlin', 'Tquad')
    np.savez_compressed(
        ROOT/CACHE_NAME, signature=np.array(signature()),
        X=np.array([row['X'] for row in rows]),
        eta=np.array([row['eta'] for row in rows]),
        radius=np.array([row['radius'] for row in rows]),
        **{key: np.stack([row[key] for row in rows])
           for key in keys})


def load_rows():
    keys = ('u0', 'g0', 'R0', 'U', 'G', 'L', 'Q',
            'T0', 'Tlin', 'Tquad')
    with np.load(ROOT/CACHE_NAME, allow_pickle=False) as data:
        if str(data['signature']) != signature():
            raise ValueError('Split swirl cache signature mismatch')
        return [dict(X=float(data['X'][i]), eta=float(data['eta'][i]),
                     radius=float(data['radius'][i]),
                     **{key: data[key][i] for key in keys})
                for i in range(len(data['X']))]


def admission(row, c, margin=.1):
    target, N, K, multiplier, lam2 = geometry(row, c)
    if multiplier is None or lam2 <= 0:
        return False
    directions = [.8*N+sign*multiplier*K
                  for sign in (-1., 1.)]
    return line_interval(directions, target, margin, 1) is not None


def run():
    if (ROOT/CACHE_NAME).exists():
        rows = load_rows()
        print('Loaded split-swirl quadratic response', flush=True)
    else:
        rows = build_rows()
        save_rows(rows)
        print('Built split-swirl quadratic response', flush=True)
    results = []
    for i, row in enumerate(rows):
        options = []
        for amplitude in np.linspace(-.1, .1, 1601):
            c = np.zeros(2)
            c[i] = amplitude
            data = cone(row, c)
            if (data['lambda_squared'] > 0 and
                    admission(row, c, margin=.2)):
                options.append(dict(amplitude=float(amplitude),
                                    lambda_squared=data['lambda_squared'],
                                    residual_norm=data['residual_norm'],
                                    target=(geometry(row, c)[0]).tolist()))
        selected = (min(options, key=lambda item:
                        (abs(item['amplitude']), item['residual_norm']))
                    if options else None)
        results.append(dict(X=row['X'], eta=row['eta'],
                            baseline=cone(row, np.zeros(2)),
                            feasible_count=len(options),
                            selected=selected,
                            amplitude_range=(
                                [options[0]['amplitude'],
                                 options[-1]['amplitude']]
                                if options else None)))
        print(json.dumps(dict(eta=row['eta'],
                              feasible_count=len(options),
                              selected=selected)), flush=True)
    report = dict(radial=RADIAL, axial_intervals=AXIAL_INTERVALS,
                  X=X, etas=ETAS, tau=.5*2**(-5.5),
                  response_cache=CACHE_NAME,
                  rows=results,
                  scope='Exact quadratic full-physical stress and '
                        'momentum response at two formerly blocked '
                        'nodes. Pointwise free-axial stress admission '
                        'with .8 cone ratio and .2 stress margin only; '
                        'no fitted pressure across all nodes, moment '
                        'closure, continuous cone, wave or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_axial_split_swirl_response.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
