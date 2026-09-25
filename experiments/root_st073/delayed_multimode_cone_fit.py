"""Joint eight-mode mean-flow screen for a physical stress-cone window.

Four compact radial velocity modes are repeated on two overlapping axial
supports. The exact quadratic Cartesian response is precomputed, then a
bounded global search balances cone violations, two-slice moments, and
local momentum. Passing samples would still require independent holdouts,
continuous cone verification, and a nonaxisymmetric wave construction.
"""

import json

import numpy as np
from scipy.optimize import differential_evolution

from azimuthal_capacity_optimize import grid
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from delayed_swirl_cone_response import BASE_NAME, blocks, stress_primitive
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from joint_collar_fit import kinematics
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


AXIAL_INTERVALS = ((.14, .31), (.19, .36))
SWIRL_INTERVALS = ((1.005, 1.04), (1.005, 1.08))
POLOIDAL_INTERVALS = ((1.005, 1.08), (1.005, 1.12))
XS = (1.008, 1.012, 1.016, 1.02)
ETAS = (.2, .25, .3)
CACHE_NAME = 'delayed_multimode_cone_model.npz'
HOLDOUT_CACHE_NAME = 'delayed_multimode_holdout_model.npz'


def signature():
    return json.dumps(dict(base=BASE_NAME, axial=AXIAL_INTERVALS,
                           swirl=SWIRL_INTERVALS,
                           poloidal=POLOIDAL_INTERVALS,
                           xs=XS, etas=ETAS, order=48), sort_keys=True)


def save_cache(path, cone_rows, moments, cache_signature=None):
    cone_keys = ('u0', 'g0', 'R0', 'U', 'G', 'L', 'Q',
                 'T0', 'Tlin', 'Tquad')
    moment_keys = ('U', 'E', 'E0', 'Bu', 'Be', 'target')
    arrays = {f'cone_{key}': np.stack([row[key] for row in cone_rows])
              for key in cone_keys}
    arrays.update({f'moment_{key}': np.stack([row[key] for row in moments])
                   for key in moment_keys})
    arrays.update(cone_X=np.array([row['X'] for row in cone_rows]),
                  cone_eta=np.array([row['eta'] for row in cone_rows]),
                  cone_radius=np.array([row['radius'] for row in cone_rows]),
                  moment_eta=np.array([row['eta'] for row in moments]),
                  moment_X=moments[0]['X'],
                  moment_weights=moments[0]['weights'],
                  signature=np.array(cache_signature or signature()))
    np.savez_compressed(path, **arrays)


def load_cache(path, cache_signature=None):
    with np.load(path, allow_pickle=False) as data:
        if str(data['signature']) != (cache_signature or signature()):
            raise ValueError('Cached response belongs to different modes')
        cone_keys = ('u0', 'g0', 'R0', 'U', 'G', 'L', 'Q',
                     'T0', 'Tlin', 'Tquad')
        moment_keys = ('U', 'E', 'E0', 'Bu', 'Be', 'target')
        cone_rows = [dict(X=float(data['cone_X'][i]),
                          eta=float(data['cone_eta'][i]),
                          radius=float(data['cone_radius'][i]),
                          **{key: data[f'cone_{key}'][i]
                             for key in cone_keys})
                     for i in range(len(data['cone_X']))]
        moments = [dict(eta=float(data['moment_eta'][i]),
                        X=data['moment_X'],
                        weights=data['moment_weights'],
                        **{key: data[f'moment_{key}'][i]
                           for key in moment_keys})
                   for i in range(len(data['moment_eta']))]
    return cone_rows, moments


def holdout_response(base, modes):
    tau = .5*2**(-5.4)
    points, X, eta = nodes(base, (1.01, 1.015, 1.025),
                           (.22, .28, .32), tau)
    hs, ht = .001*np.sqrt(base.nu*tau), .00025*tau
    u0, g0, p0 = kinematics(base, points, tau, hs, ht)
    R0 = p0+np.einsum('pab,pb->pa', g0, u0)
    modes_data = [kinematics(mode, points, tau, hs, ht)
                  for mode in modes]
    U = np.stack([item[0] for item in modes_data], axis=2)
    G = np.stack([item[1] for item in modes_data], axis=3)
    P = np.stack([item[2] for item in modes_data], axis=2)
    L = (P+np.einsum('pab,pbi->pai', g0, U)
         +np.einsum('pabi,pb->pai', G, u0))
    Q = np.einsum('pabi,pbj->paij', G, U)
    return dict(tau=tau, X=X, eta=eta, R0=R0, L=L, Q=Q)


def evaluate_holdout(data, coefficients):
    return (data['R0']+data['L']@coefficients
            +np.einsum('paij,i,j->pa', data['Q'],
                       coefficients, coefficients))


def get_holdout_cache(base, modes):
    path = ROOT/HOLDOUT_CACHE_NAME
    cache_signature = signature()+'|holdout-v1'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            if str(data['signature']) != cache_signature:
                raise ValueError('Cached holdout belongs to different modes')
            return {key: data[key] for key in ('tau', 'X', 'eta',
                                                'R0', 'L', 'Q')}
    data = holdout_response(base, modes)
    np.savez_compressed(path, signature=np.array(cache_signature), **data)
    return data


def make_modes(base):
    modes, kinds = [], []
    for axial in AXIAL_INTERVALS:
        for radial in SWIRL_INTERVALS:
            modes.append(SimilaritySwirlMode(
                base, radial, eta_interval=axial))
            kinds.append('swirl')
        for radial in POLOIDAL_INTERVALS:
            modes.append(SimilarityCurlMode(base, radial,
                                            eta_interval=axial))
            kinds.append('poloidal')
    return modes, kinds


def precompute(base, modes, tau, extra_radial_edges=()):
    edges = [edge for interval in (*SWIRL_INTERVALS,
                                   *POLOIDAL_INTERVALS)
             for edge in interval]
    points, blocks_data = blocks(
        base, tau, order=48, xs=XS, etas=ETAS,
        extra_radial_edges=(*edges, *extra_radial_edges))
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u0, g0, p0 = kinematics(base, points, tau, hs, ht)
    R0 = p0+np.einsum('pab,pb->pa', g0, u0)
    mode_data = [kinematics(mode, points, tau, hs, ht)
                 for mode in modes]
    U = np.stack([item[0] for item in mode_data], axis=2)
    G = np.stack([item[1] for item in mode_data], axis=3)
    P = np.stack([item[2] for item in mode_data], axis=2)
    L = (P+np.einsum('pab,pbi->pai', g0, U)
         +np.einsum('pabi,pb->pai', G, u0))
    rows = []
    for block in blocks_data:
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
        center = block['center']
        rows.append(dict(X=block['X'], eta=block['eta'], radius=radius,
                         u0=u0[center], g0=g0[center], R0=R0[center],
                         U=U[center], G=G[center], L=L[center],
                         Q=np.einsum('abi,bj->aij', G[center], U[center]),
                         T0=stress_primitive(R0, block), Tlin=Tlin,
                         Tquad=np.stack((qtheta, qaxial))))
    return rows


def cone(row, c):
    u = row['u0']+row['U']@c
    g = row['g0']+np.einsum('abi,i->ab', row['G'], c)
    target = (row['T0']+row['Tlin']@c
              +np.einsum('aij,i,j->a', row['Tquad'], c, c))
    F = u[1]/row['radius']
    shear = np.array([g[1, 0]-F, g[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    dot_n, dot_k = float(target@N), float(target@K)
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14:
        multiplier = np.sqrt(lam2)/(2*F*N[0])
        side = abs(multiplier*dot_k)
        margin = -dot_n-side
        violation = max(0., -margin)/(1+abs(dot_n)+side)
        ratio = side/abs(dot_n) if abs(dot_n) > 1e-14 else None
    else:
        margin, ratio = -1e12, None
        violation = 1.+min(1., max(0., -lam2)/1e4)
    R = row['R0']+row['L']@c+np.einsum('aij,i,j->a', row['Q'], c, c)
    return dict(X=row['X'], eta=row['eta'],
                lambda_squared=float(lam2),
                target_dot_N=dot_n, target_dot_K=dot_k,
                margin=float(margin), ratio=ratio,
                cone_violation=float(violation),
                strict_pass=bool(margin > 0),
                residual_norm=float(np.linalg.norm(R)))


def moment_rows(base, target, modes, tau):
    X, weights = grid(order=48)
    rows = []
    for eta in ETAS:
        U0, E0 = profile(target, X, eta, tau)
        U, E = profile(base, X, eta, tau)
        basis = [mode_profile(mode, base, X, eta, tau)
                 for mode in modes]
        Bu = np.column_stack([item[0] for item in basis])
        Be = np.column_stack([item[1] for item in basis])
        rows.append(dict(eta=eta, X=X, weights=weights,
                         U=U, E=E, E0=E0, Bu=Bu, Be=Be,
                         target=moment_vector(U0, E0, X, weights)))
    return rows


def moment_summary(rows, c):
    result = []
    for row in rows:
        U, E = row['U']+row['Bu']@c, row['E']+row['Be']@c
        defect = (moment_vector(U, E, row['X'], row['weights'])
                  -row['target'])
        result.append(dict(eta=row['eta'], defect=defect.tolist(),
                           min_relative_E=float(np.min(E/row['E0']))))
    return result


def run(cache_only=False):
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes, kinds = make_modes(base)
    cache_path = ROOT/CACHE_NAME
    if cache_path.exists():
        cone_rows, moments = load_cache(cache_path)
        print('Loaded cached quadratic physical-cone response', flush=True)
    else:
        cone_rows = precompute(base, modes, tau)
        print('Precomputed full quadratic physical-cone response', flush=True)
        target = make_field(16, 2.)
        moments = moment_rows(base, target, modes, tau)
        save_cache(cache_path, cone_rows, moments)
    if cache_only:
        print(json.dumps(dict(response_cache=CACHE_NAME,
                              cone_rows=len(cone_rows),
                              moment_rows=len(moments))), flush=True)
        return
    holdout_model = get_holdout_cache(base, modes)
    baseline = [cone(row, np.zeros(len(modes))) for row in cone_rows]
    base_max = max(row['residual_norm'] for row in baseline)
    base_holdout_max = stats(holdout_model['R0'])['max']

    def objective(c):
        cone_data = [cone(row, c) for row in cone_rows]
        violations = np.array([row['cone_violation'] for row in cone_data])
        moment_data = moment_summary(moments, c)
        floor = min(row['min_relative_E'] for row in moment_data)
        debt = max(max(abs(x)/scale for x, scale in zip(
            row['defect'][1:], (.001, .01, .01, .001)))
            for row in moment_data)
        peak = max(row['residual_norm'] for row in cone_data)
        holdout_peak = stats(evaluate_holdout(holdout_model, c))['max']
        return (float(np.max(violations)+.25*np.mean(violations))
                +.02*debt+.03*peak/base_max
                +max(0., holdout_peak/base_holdout_max-1.)
                +100*max(0., .11005-floor))

    bounds = [(-.05, .3) if kind == 'swirl' else (-.05, .05)
              for kind in kinds]
    fit = differential_evolution(objective, bounds, seed=73074,
                                 maxiter=120, popsize=8, polish=True,
                                 tol=1e-5, updating='immediate')
    c = fit.x
    selected = [cone(row, c) for row in cone_rows]
    selected_moments = moment_summary(moments, c)
    patched = CurlPatchedLift(base, modes, c)
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(base, (1.01, 1.015, 1.025),
                             (.22, .28, .32), tau_holdout)
    before, _ = residual(base, points, tau_holdout)
    after, divergence = residual(patched, points, tau_holdout)
    report = dict(base_slice=BASE_NAME, tau=tau,
                  response_cache=CACHE_NAME,
                  holdout_response_cache=HOLDOUT_CACHE_NAME,
                  axial_intervals=AXIAL_INTERVALS,
                  swirl_intervals=SWIRL_INTERVALS,
                  poloidal_intervals=POLOIDAL_INTERVALS,
                  X=XS, eta=ETAS, kinds=kinds,
                  coefficients=c.tolist(),
                  optimizer_success=bool(fit.success),
                  optimizer_message=fit.message,
                  objective=float(fit.fun),
                  baseline=baseline, selected=selected,
                  moments=selected_moments,
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(), before=stats(before),
                               after=stats(after),
                               response_vs_fd_max=float(np.max(np.abs(
                                   evaluate_holdout(holdout_model, c)-after))),
                               max_abs_fd_divergence=float(np.max(np.abs(
                                   divergence)))),
                  scope='Eight-mode compact mean velocity screen. '
                        'Full nonlinear Cartesian stress response, '
                        'three-eta moment/floor penalty and separate '
                        'space/time momentum holdout. No continuous '
                        'cone, exact moment closure, supported wave, '
                        'volume norm, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_multimode_cone_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(optimizer_success=report['optimizer_success'],
                          objective=report['objective'],
                          coefficients=report['coefficients'],
                          pass_count=sum(row['strict_pass'] for row in selected),
                          min_margin=min(row['margin'] for row in selected),
                          max_train_residual=max(row['residual_norm']
                                                 for row in selected),
                          holdout=report['holdout'])), flush=True)


if __name__ == '__main__':
    import sys
    run(cache_only='--cache-only' in sys.argv[1:])
