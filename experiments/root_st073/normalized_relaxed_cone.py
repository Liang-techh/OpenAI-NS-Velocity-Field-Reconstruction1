"""Evaluate the paper's normalized relaxed cone on snapshot profiles.

The joined field is not the paper's exact leading profile. This diagnostic
reconstructs Pi from its axis datum and E via (4.7), then evaluates (4.15),
(4.16), (4.11), and (4.21) on finite-window snapshots.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from joined_field import JoinedField, ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import radial_boundaries


def profile(field, X, eta, tau):
    points = field.inner.from_similarity(np.asarray(X, float), eta, tau)
    velocity, _ = field.fields(points, tau)
    q = tau/(1-eta**2)
    scale = q**(.5+field.inner.h)/np.sqrt(field.nu)
    return scale*velocity[:, 2], scale*velocity[:, 1]


def moment_data(field, X, eta, tau, order=20):
    g, w = leggauss(order)
    edges = sorted(set(radial_boundaries(
        getattr(field, 'patch_start', 4.),
        getattr(field, 'patch_end', 4.))+[X]))
    nodes, weights = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo >= X:
            break
        hi = min(hi, X)
        nodes.append((lo+hi)/2+(hi-lo)/2*g)
        weights.append((hi-lo)/2*w)
    x = np.concatenate(nodes)
    weights = np.concatenate(weights)
    U, E = profile(field, x, eta, tau)
    H = np.sqrt(2*x)*E
    values = np.array([weights@U, weights@H, weights@(U*H),
                       weights@(U**2-E**2/2), weights@(E**2/(2*x))])
    axis = field.inner.from_similarity(np.array([0.]), eta, tau)
    _, pressure = field.fields(axis, tau)
    q = tau/(1-eta**2)
    A = .5+field.inner.h
    pi_axis = float(pressure[0]/(field.nu*q**(-2*A)))
    return values, pi_axis


def cone_point(field, X, eta, tau, order=20, eta_step=.001,
               radial_step_fraction=.001):
    h = field.inner.h
    D, A = .5-h, .5+h
    d, L = 1-eta**2, 1-2*h*eta**2
    he, hx = eta_step, radial_step_fraction*X
    m, pi_axis = moment_data(field, X, eta, tau, order)
    mp, pi_axis_p = moment_data(field, X, eta+he, tau, order)
    mm, pi_axis_m = moment_data(field, X, eta-he, tau, order)
    me = (mp-mm)/(2*he)
    pi = pi_axis+m[4]
    pi_eta = (pi_axis_p+mp[4]-pi_axis_m-mm[4])/(2*he)
    U, E = profile(field, np.array([X]), eta, tau)
    Up, Ep = profile(field, np.array([X+hx]), eta, tau)
    Um, Em = profile(field, np.array([X-hx]), eta, tau)
    U, E = float(U[0]), float(E[0])
    if E <= 0:
        return {'X': X, 'eta': eta, 'error': 'E nonpositive'}
    Ux = float((Up[0]-Um[0])/(2*hx))
    Ex = float((Ep[0]-Em[0])/(2*hx))
    M, I, J, S, Cp = m
    Me, Ie, Je, Se, _ = me
    H = np.sqrt(2*X)*E
    W = 1-2*D*eta*M/X-d*Me/X
    Qs = -W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*H)
    Ns = (-W*U+(D*(M-eta*Me)+4*h*eta*S-d*Se)/X
          +4*A*eta*pi-d*pi_eta)
    ps1, ps2 = X*Qs/L, X*Ns/(L*E)
    a, bs = 1-2*X*Ex/E, 2*X*Ux/E
    if a <= 0:
        return {'X': X, 'eta': eta, 'a': a, 'bs': bs,
                'ps': [ps1, ps2], 'relaxed_pass': False,
                'admissible_pass': False, 'failure': 'a nonpositive'}
    ts = -bs/a
    vs = a*(1+ts**2)
    Pc, Jc = ps1+ts*ps2, ps2-ts*ps1
    if Pc > 2:
        upper = (Pc+Jc**2/4
                 -abs(Jc)*np.sqrt((Pc-2)/2+Jc**2/16))
    else:
        upper = None
    relaxed = bool(Pc > 2 and vs < upper)
    return {'X': X, 'eta': eta, 'a': a, 'bs': bs,
            'ps': [ps1, ps2], 'pi': pi, 'W': W,
            'vs': vs, 'Pc': Pc, 'Jc': Jc, 'upper': upper,
            'relaxed_pass': relaxed,
            'admissible_pass': bool(relaxed and vs > 2)}


def run():
    base = JoinedField()
    fields = (('joined', base),
              ('two_moment', MomentMatchedJoinedField(base=base)),
              ('two_shape', MomentMatchedJoinedField(
                  base=base, meridional_null_amplitude=-3.,
                  meridional_even_amplitude=12.)))
    tau = .0084
    rows = []
    for name, field in fields:
        points = [cone_point(field, X, eta, tau)
                  for eta in (0., .2, .35)
                  for X in (6., 10., 14.)]
        row = {'field': name, 'points': points,
               'relaxed_pass_count': sum(p['relaxed_pass'] for p in points),
               'admissible_pass_count': sum(p['admissible_pass'] for p in points)}
        rows.append(row)
        print(json.dumps({'field': name,
                          'relaxed_pass_count': row['relaxed_pass_count'],
                          'admissible_pass_count': row['admissible_pass_count']}),
              flush=True)
    report = {'tau': tau, 'rows': rows,
              'scope': 'Paper-normalized (4.21) snapshot diagnostic reconstructed from U,E, consistent Pi radial integral, and finite-difference eta/X derivatives. The joined field is not proven to be an exact time-independent leading profile; sampled inequalities are not a continuum certificate. Physical pressure remains unchanged in these field objects.',
              'accepted': False}
    (ROOT/'normalized_relaxed_cone.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
