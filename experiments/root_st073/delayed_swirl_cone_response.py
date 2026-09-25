"""Nonlinear physical-cone response to compact solenoidal swirl modes.

For each mode, full Cartesian momentum is quadratic in its amplitude.
The radial primitive, local shear, and cone are then screened cheaply on
a common two-eta radial window. This is a single-mode screen, not a wave.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from azimuthal_capacity_optimize import INTERVALS as E_INTERVALS, grid
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from joint_collar_fit import kinematics
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


BASE_NAME = 'delayed005_rise146_degree31.json'
RADIAL_INTERVALS = ((1.005, 1.04), (1.005, 1.08), (1.005, 1.12))
ETA_INTERVAL = (.14, .36)
XS = (1.008, 1.012, 1.016, 1.02)
ETAS = (.2, .3)


def blocks(base, tau, order=48):
    inner = base.compact.joined.inner
    joined = base.compact.joined
    g, w = leggauss(order)
    points, metadata = [], []
    radial_edges = [joined.join_X, 1., 1.005, base.start_X,
                    base.start_X+base.width, 3.-base.width, 3.]
    radial_edges += [edge for interval in E_INTERVALS for edge in interval]
    radial_edges += [edge for interval in RADIAL_INTERVALS
                     for edge in interval]
    for eta in ETAS:
        for X in XS:
            point = inner.from_similarity(np.array([X]),
                                          np.array([eta]), tau)[0]
            radius, _, z = point
            q = tau/(1-eta**2)
            cuts = [0., radius]
            cuts += list(np.sqrt(2*base.nu*q*np.array(radial_edges)))
            edges = sorted(set(np.clip(cuts, 0., radius)))
            rr, ww = [], []
            for lo, hi in zip(edges[:-1], edges[1:]):
                if hi <= lo:
                    continue
                rr.extend((lo+hi)/2+(hi-lo)/2*g)
                ww.extend((hi-lo)/2*w)
            rr, ww = np.asarray(rr), np.asarray(ww)
            start = len(points)
            points.extend(np.column_stack((rr, np.zeros(len(rr)),
                                           np.full(len(rr), z))))
            center = len(points)
            points.append(point)
            metadata.append(dict(X=X, eta=eta, radius=radius,
                                 start=start, stop=center, center=center,
                                 rr=rr, weights=ww))
    return np.asarray(points), metadata


def residual_terms(base, mode, points, tau, cached_base=None):
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    if cached_base is None:
        u0, grad0, part0 = kinematics(base, points, tau, hs, ht)
        R0 = part0+np.einsum('nij,nj->ni', grad0, u0)
    else:
        u0, grad0, R0 = cached_base
    u1, grad1, part1 = kinematics(mode, points, tau, hs, ht)
    L = (part1+np.einsum('nij,nj->ni', grad0, u1)
         +np.einsum('nij,nj->ni', grad1, u0))
    Q = np.einsum('nij,nj->ni', grad1, u1)
    return (u0, grad0, R0), (u1, grad1, L, Q)


def stress_primitive(R, row):
    sl = slice(row['start'], row['stop'])
    rr, weights, radius = row['rr'], row['weights'], row['radius']
    return np.array([-(weights*rr**2)@R[sl, 1]/radius**2,
                     -(weights*rr)@R[sl, 2]/radius])


def cone_metric(row, base_data, mode_data, amplitude):
    u0, grad0, R0 = base_data
    u1, grad1, L, Q = mode_data
    i, radius = row['center'], row['radius']
    u = u0[i]+amplitude*u1[i]
    grad = grad0[i]+amplitude*grad1[i]
    target = (stress_primitive(R0, row)
              +amplitude*stress_primitive(L, row)
              +amplitude**2*stress_primitive(Q, row))
    F = u[1]/radius
    shear = np.array([grad[1, 0]-F, grad[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    dot_n, dot_k = float(target@N), float(target@K)
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14:
        c = np.sqrt(lam2)/(2*F*N[0])
        margin = -dot_n-abs(c*dot_k)
        ratio = (abs(c*dot_k/dot_n) if abs(dot_n) > 1e-14
                 else None)
    else:
        margin, ratio = -1e12, None
    R = R0[i]+amplitude*L[i]+amplitude**2*Q[i]
    return dict(X=row['X'], eta=row['eta'],
                lambda_squared=float(lam2),
                target_dot_N=dot_n, target_dot_K=dot_k,
                margin=float(margin), ratio=ratio,
                strict_pass=bool(margin > 0),
                residual_norm=float(np.linalg.norm(R)))


def moment_profiles(base, target, X, weights, tau):
    rows = []
    for eta in ETAS:
        U0, E0 = profile(target, X, eta, tau)
        U, E = profile(base, X, eta, tau)
        rows.append(dict(eta=eta, U=U, E=E, E0=E0,
                         target=moment_vector(U0, E0, X, weights)))
    return rows


def run():
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    target = make_field(16, 2.)
    points, quadrature = blocks(base, tau)
    X_m, weights = grid(order=48)
    moment_rows = moment_profiles(base, target, X_m, weights, tau)
    amplitudes = np.unique(np.r_[np.linspace(-2., -.5, 31),
                                  np.linspace(-.5, .5, 101),
                                  np.linspace(.5, 2., 31)])
    mode_results = []
    cached_base = None
    for interval in RADIAL_INTERVALS:
        mode = SimilaritySwirlMode(base, interval,
                                   eta_interval=ETA_INTERVAL)
        base_data, mode_data = residual_terms(
            base, mode, points, tau, cached_base)
        cached_base = base_data
        Em = [mode_profile(mode, base, X_m, row['eta'], tau)[1]
              for row in moment_rows]
        baseline = [cone_metric(row, base_data, mode_data, 0.)
                    for row in quadrature]
        trials = []
        for amplitude in amplitudes:
            relative_floor = min(np.min((row['E']+amplitude*em)/row['E0'])
                                 for row, em in zip(moment_rows, Em))
            if relative_floor < .11:
                continue
            metrics = [cone_metric(row, base_data, mode_data, amplitude)
                       for row in quadrature]
            trials.append(dict(amplitude=float(amplitude),
                               pass_count=sum(item['strict_pass']
                                              for item in metrics),
                               min_margin=min(item['margin']
                                              for item in metrics),
                               max_residual=max(item['residual_norm']
                                                for item in metrics),
                               min_relative_E=float(relative_floor)))
        selected = max(trials, key=lambda item: (
            item['pass_count'], item['min_margin'],
            -item['max_residual']))
        amplitude = selected['amplitude']
        selected_metrics = [cone_metric(row, base_data, mode_data, amplitude)
                            for row in quadrature]
        moment_defects = []
        for row, em in zip(moment_rows, Em):
            achieved = moment_vector(row['U'], row['E']+amplitude*em,
                                     X_m, weights)
            moment_defects.append(dict(eta=row['eta'],
                                       defect=(achieved-row['target']).tolist()))
        mode_results.append(dict(radial_interval=interval,
                                 baseline=baseline,
                                 selected=selected,
                                 selected_metrics=selected_metrics,
                                 moment_defects=moment_defects,
                                 sampled_trials=len(trials)))
        print(json.dumps(dict(radial_interval=interval,
                              selected=selected)), flush=True)
    report = dict(base_slice=BASE_NAME, tau=tau,
                  eta_interval=ETA_INTERVAL,
                  radial_intervals=RADIAL_INTERVALS,
                  X=XS, eta=ETAS, mode_results=mode_results,
                  scope='Single compact toroidal velocity mode at a time. '
                        'Full quadratic Cartesian residual and radial '
                        'primitive on sampled nodes, with two-slice '
                        'positive swirl floor. No exact moment closure, '
                        'continuous cone, wave, volume norm, or PDE gate.',
                  accepted=False)
    (ROOT/'delayed_swirl_cone_response.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
