"""Screen coupled compact poloidal and swirl velocity cone response.

All quadratic and cross-advection terms in the physical residual are kept.
The two-dimensional amplitude screen is exploratory; outgoing moments,
continuous cone, and a nonaxisymmetric wave remain separate obligations.
"""

import json

import numpy as np

from azimuthal_capacity_optimize import grid
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_similarity_curl_screen import SimilarityCurlMode
from delayed_swirl_cone_response import (
    BASE_NAME, ETA_INTERVAL, ETAS, XS, blocks, residual_terms,
    stress_primitive,
)
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


SWIRL_INTERVAL = (1.005, 1.08)
POLOIDAL_INTERVAL = (1.005, 1.08)


def prepared_row(row, base_data, swirl_data, curl_data):
    u0, g0, R0 = base_data
    us, gs, Ls, Qss = swirl_data
    up, gp, Lp, Qpp = curl_data
    Qsp = (np.einsum('nij,nj->ni', gs, up)
           +np.einsum('nij,nj->ni', gp, us))
    i = row['center']
    stress = [stress_primitive(R, row)
              for R in (R0, Ls, Lp, Qss, Qpp, Qsp)]
    return dict(X=row['X'], eta=row['eta'], radius=row['radius'],
                u=(u0[i], us[i], up[i]),
                grad=(g0[i], gs[i], gp[i]),
                residual=(R0[i], Ls[i], Lp[i], Qss[i], Qpp[i], Qsp[i]),
                stress=stress)


def metric(row, swirl, poloidal):
    s, p = swirl, poloidal
    u = row['u'][0]+s*row['u'][1]+p*row['u'][2]
    g = row['grad'][0]+s*row['grad'][1]+p*row['grad'][2]
    t0, ts, tp, tss, tpp, tsp = row['stress']
    target = t0+s*ts+p*tp+s*s*tss+p*p*tpp+s*p*tsp
    F = u[1]/row['radius']
    shear = np.array([g[1, 0]-F, g[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    dot_n, dot_k = float(target@N), float(target@K)
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14:
        c = np.sqrt(lam2)/(2*F*N[0])
        margin = -dot_n-abs(c*dot_k)
        ratio = abs(c*dot_k/dot_n) if abs(dot_n) > 1e-14 else None
    else:
        margin, ratio = -1e12, None
    R0, Ls, Lp, Qss, Qpp, Qsp = row['residual']
    R = R0+s*Ls+p*Lp+s*s*Qss+p*p*Qpp+s*p*Qsp
    return dict(X=row['X'], eta=row['eta'],
                lambda_squared=float(lam2),
                target_dot_N=dot_n, target_dot_K=dot_k,
                margin=float(margin), ratio=ratio,
                strict_pass=bool(margin > 0),
                residual_norm=float(np.linalg.norm(R)))


def run():
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    target = make_field(16, 2.)
    points, quadrature = blocks(base, tau)
    swirl_mode = SimilaritySwirlMode(
        base, SWIRL_INTERVAL, eta_interval=ETA_INTERVAL)
    curl_mode = SimilarityCurlMode(base, POLOIDAL_INTERVAL,
                                   eta_interval=ETA_INTERVAL)
    base_data, swirl_data = residual_terms(base, swirl_mode, points, tau)
    _, curl_data = residual_terms(base, curl_mode, points, tau, base_data)
    rows = [prepared_row(row, base_data, swirl_data, curl_data)
            for row in quadrature]
    X_m, weights = grid(order=48)
    moments = []
    for eta in ETAS:
        U0, E0 = profile(target, X_m, eta, tau)
        U, E = profile(base, X_m, eta, tau)
        Up = mode_profile(curl_mode, base, X_m, eta, tau)[0]
        Es = mode_profile(swirl_mode, base, X_m, eta, tau)[1]
        moments.append(dict(eta=eta, U=U, E=E, E0=E0,
                            Up=Up, Es=Es,
                            target=moment_vector(U0, E0, X_m, weights)))
    baseline = [metric(row, 0., 0.) for row in rows]
    candidates = []
    positive_lambda = []
    per_eta_best = {eta: None for eta in ETAS}
    for s in np.linspace(-.1, .3, 81):
        floor = min(np.min((row['E']+s*row['Es'])/row['E0'])
                    for row in moments)
        if floor < .11:
            continue
        for p in np.linspace(-.1, .1, 81):
            metrics = [metric(row, s, p) for row in rows]
            summary = dict(swirl=float(s), poloidal=float(p),
                           pass_count=sum(item['strict_pass']
                                          for item in metrics),
                           min_margin=min(item['margin'] for item in metrics),
                           max_residual=max(item['residual_norm']
                                            for item in metrics),
                           min_relative_E=float(floor))
            candidates.append(summary)
            for eta in ETAS:
                local = [item for item in metrics if item['eta'] == eta]
                local_summary = dict(swirl=float(s), poloidal=float(p),
                                     pass_count=sum(item['strict_pass']
                                                    for item in local),
                                     min_margin=min(item['margin']
                                                    for item in local),
                                     max_residual=max(item['residual_norm']
                                                      for item in local))
                current = per_eta_best[eta]
                rank = lambda item: (item['pass_count'],
                                     item['min_margin'],
                                     -item['max_residual'])
                if current is None or rank(local_summary) > rank(current):
                    per_eta_best[eta] = local_summary
            if all(item['lambda_squared'] > 0 for item in metrics):
                positive_lambda.append(summary)
    selected = max(candidates, key=lambda item: (
        item['pass_count'], item['min_margin'], -item['max_residual']))
    best_positive_lambda = (max(positive_lambda,
                                key=lambda item: item['min_margin'])
                            if positive_lambda else None)
    selected_metrics = [metric(row, selected['swirl'], selected['poloidal'])
                        for row in rows]
    moment_defects = []
    for row in moments:
        achieved = moment_vector(row['U']+selected['poloidal']*row['Up'],
                                 row['E']+selected['swirl']*row['Es'],
                                 X_m, weights)
        moment_defects.append(dict(eta=row['eta'],
                                   defect=(achieved-row['target']).tolist()))
    report = dict(base_slice=BASE_NAME, tau=tau,
                  swirl_interval=SWIRL_INTERVAL,
                  poloidal_interval=POLOIDAL_INTERVAL,
                  eta_interval=ETA_INTERVAL, X=XS, eta=ETAS,
                  candidate_count=len(candidates),
                  positive_lambda_count=len(positive_lambda),
                  per_eta_best=list(per_eta_best.values()),
                  baseline=baseline, selected=selected,
                  best_positive_lambda=best_positive_lambda,
                  selected_metrics=selected_metrics,
                  moment_defects=moment_defects,
                  scope='Two compact solenoidal velocity modes with '
                        'quadratic cross-advection. Sampled physical '
                        'stress-cone response only; no exact moment '
                        'repair, supported wave, continuous cone, '
                        'volume norm, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_coupled_cone_response.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in (
        'candidate_count', 'positive_lambda_count', 'selected',
        'best_positive_lambda', 'per_eta_best')}), flush=True)


if __name__ == '__main__':
    run()
