"""Paper-inspired relaxed-cone snapshot on the extended ST073 background.

This uses the Section 4 normalized formulas but the background is a
time-dependent approximate physical field, not the paper's exact leading
profile. The radial quadrature resolves its actual moving attachment edge.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

from extended_compact_join import load_extended_heated_candidate
from radial_continuation import ROOT


def profile(field, X, eta, tau):
    X = np.asarray(X, float)
    points = field.compact.joined.inner.from_similarity(
        X, np.full(len(X), eta), tau)
    velocity, _ = field.fields(points, tau)
    q = tau/(1-eta**2)
    scale = q**(.5+field.heat.h)/np.sqrt(field.nu)
    return scale*velocity[:, 2], scale*velocity[:, 1]


def moment_data(field, X, eta, tau, order):
    inner = field.compact.joined.inner
    radial = field.compact.joined
    g, w = leggauss(order)
    q = tau/(1-eta**2)
    pure_heat_X = float(field.attachment_radius(tau)**2/(2*field.nu*q))
    cuts = sorted(set([0., X]+[edge for edge in
                       (radial.join_X, radial.outer_ratio**2*radial.join_X,
                        pure_heat_X) if 0 < edge < X]))
    nodes, weights = [], []
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        nodes.append((lo+hi)/2+(hi-lo)/2*g)
        weights.append((hi-lo)/2*w)
    xs, ws = np.concatenate(nodes), np.concatenate(weights)
    U, E = profile(field, xs, eta, tau)
    H = np.sqrt(2*xs)*E
    moments = np.array([ws@U, ws@H, ws@(U*H),
                        ws@(U**2-E**2/2), ws@(E**2/(2*xs))])
    axis = inner.from_similarity(np.array([0.]), np.array([eta]), tau)
    _, pressure = field.fields(axis, tau)
    A = .5+field.heat.h
    pi_axis = float(pressure[0]/(field.nu*q**(-2*A)))
    return moments, pi_axis, pure_heat_X


def cone_point(field, X, eta, tau, order=8, eta_step=.0008):
    h = field.heat.h
    D, A = .5-h, .5+h
    d, L = 1-eta**2, 1-2*h*eta**2
    m, pi_axis, heat_edge = moment_data(field, X, eta, tau, order)
    mp, pi_p, _ = moment_data(field, X, eta+eta_step, tau, order)
    mm, pi_m, _ = moment_data(field, X, eta-eta_step, tau, order)
    me = (mp-mm)/(2*eta_step)
    pi = pi_axis+m[4]
    pi_eta = (pi_p+mp[4]-pi_m-mm[4])/(2*eta_step)
    hx = .001*X
    U, E = profile(field, np.array([X]), eta, tau)
    Up, Ep = profile(field, np.array([X+hx]), eta, tau)
    Um, Em = profile(field, np.array([X-hx]), eta, tau)
    U, E = float(U[0]), float(E[0])
    if E <= 0:
        return dict(X=X, eta=eta, failure='E nonpositive')
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
        return dict(X=X, eta=eta, a=a, failure='a nonpositive',
                    relaxed_pass=False)
    ts = -bs/a
    vs = a*(1+ts**2)
    Pc, Jc = ps1+ts*ps2, ps2-ts*ps1
    upper = (Pc+Jc**2/4-abs(Jc)*np.sqrt((Pc-2)/2+Jc**2/16)
             if Pc > 2 else None)
    relaxed = bool(Pc > 2 and vs < upper)
    return dict(X=X, eta=eta, pure_heat_X=heat_edge,
                moments=m.tolist(), pi_axis=pi_axis,
                a=a, vs=vs, Pc=Pc, Jc=Jc, upper=upper,
                relaxed_pass=relaxed,
                admissible_pass=bool(relaxed and vs > 2))


def run():
    field = load_extended_heated_candidate()
    rows = []
    for k in (3., 5.5, 6.):
        tau = .5*2**(-k)
        pairs = ([(e, X) for e in (.2, .3, .4)
                  for X in (1.5, 1.75, 2., 2.25, 2.5, 8.)]
                 if k == 5.5 else [(.2, 2.), (.3, 2.)])
        for eta, X in pairs:
            row = cone_point(field, X, eta, tau, order=12)
            row['k'] = k
            row['margin'] = (row['upper']-row['vs']
                             if row.get('upper') is not None else None)
            rows.append(row)
            print(json.dumps(dict(k=k, eta=eta, X=X,
                                  margin=row['margin'],
                                  relaxed_pass=row.get('relaxed_pass'))),
                  flush=True)
    report = dict(order=12, rows=rows,
                  scope='Sampled Section-4 relaxed-cone analogs with dynamic radial quadrature edges. Background is not an exact normalized leading profile. Positive nodes do not prove a continuous space-time stress cone or global PDE acceptance.',
                  accepted=False)
    (ROOT/'extended_relaxed_cone_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
