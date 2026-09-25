"""Physical exact-curl lift of the protected eta=.3 remote U slice.

The radial streamfunction returns to zero outside X=1.03..3.5.
The remaining tiny M defect is reported explicitly. This is a sampled
PDE screen, not a completed background matching construction.
"""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_u_basis_capacity import u_basis
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from joined_field import coordinates
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


START, END, WIDTH, DEGREE = 1.03, 3.5, .002, 63
AXIAL = (.26, .34)


class RemoteUPatchedField:
    def __init__(self, mean, base, e_intervals, e_coefficients,
                 u_coefficients, order=64):
        self.mean = mean
        self.base = base
        self.e_intervals = e_intervals
        self.e_coefficients = np.asarray(e_coefficients)
        self.u_coefficients = np.asarray(u_coefficients)
        self.nu = base.nu
        self.heat = base.heat
        self.A = base.A
        self.compact = base.compact
        self.g, self.w = leggauss(order)

    def primitive_basis(self, X):
        end = min(max(float(X), START), END)
        total = np.zeros(DEGREE+1)
        for lo, hi in ((START, START+WIDTH),
                       (START+WIDTH, END-WIDTH),
                       (END-WIDTH, END)):
            upper = min(end, hi)
            if upper <= lo:
                continue
            xx = (lo+upper)/2+(upper-lo)*self.g/2
            ww = (upper-lo)*self.w/2
            total += u_basis(xx, START, END, WIDTH, DEGREE)@ww
        return total

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.mean.fields(pts, ts)
        r = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(r/sn, pts[:, 2]/sn, ts, self.heat.h)
        for i, (radius, X, eta, q) in enumerate(zip(
                r, co['X'], co['eta'], co['q'])):
            if not (radius > 0 and START+1e-10 < X < END-1e-10
                    and AXIAL[0] < eta < AXIAL[1]):
                continue
            axial, axial_d = bump(np.array([eta]), *AXIAL)
            f, fd = float(axial[0]), float(axial_d[0])
            b = u_basis(np.array([X]), START, END, WIDTH, DEGREE)[:, 0]
            primitive = self.primitive_basis(X)
            G = f*float(self.u_coefficients@primitive)
            G_x = f*float(self.u_coefficients@b)
            G_eta = fd*float(self.u_coefficients@primitive)
            e_bumps = np.array([bump(np.array([X]), *interval)[0][0]
                                for interval in self.e_intervals])
            e_delta = f*float(self.e_coefficients@e_bumps)
            qz = float(co['q_z'][i])/sn
            Xz = float(co['X_z'][i])/sn
            etaz = float(co['eta_z'][i])/sn
            psi_z = self.nu**1.5*((1-self.A)*q**(-self.A)*qz*G
                      +q**(1-self.A)*(G_x*Xz+G_eta*etaz))
            ur = -psi_z/radius
            uz = sn*q**(-self.A)*G_x
            ut = sn*q**(-self.A)*e_delta
            ca, sa = pts[i, 0]/radius, pts[i, 1]/radius
            velocity[i] += np.array([ur*ca-ut*sa, ur*sa+ut*ca, uz])
        return velocity, pressure


def resolved_quadrature(order=96):
    # Rebuild the original piecewise rule with two extra taper edges.
    edges = sorted(set((0., 3/32, 1., START, START+WIDTH, 1.03, 1.04,
                        1.18, 1.2, 1.4, 1.45, 1.75, 1.8,
                        2.3, 2.35, 2.98, 3., END-WIDTH, END)))
    g, w = leggauss(order)
    X = np.concatenate([(lo+hi)/2+(hi-lo)*g/2
                        for lo, hi in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(hi-lo)*w/2
                              for lo, hi in zip(edges[:-1], edges[1:])])
    return X, weights


def solve_coefficients(mean, base, target_field, e_source, e_row):
    X, weights = resolved_quadrature()
    eta, tau = .3, .5*2**(-5.5)
    U, E = profile(mean, X, eta, tau)
    U0, E0 = profile(target_field, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    e_intervals = [tuple(interval)
                   for interval in e_source['radial_intervals']]
    source_factor = bump(np.array([eta]),
                         *e_source['mode_eta_interval'])[0][0]
    e_coefficients = source_factor*np.asarray(e_row['coefficients'])
    e_bumps = np.array([bump(X, *interval)[0]
                        for interval in e_intervals])
    E_new = E+e_coefficients@e_bumps
    B = u_basis(X, START, END, WIDTH, DEGREE)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    field_template = RemoteUPatchedField(mean, base, e_intervals,
                                         e_coefficients,
                                         np.zeros(DEGREE+1))
    mass_basis = field_template.primitive_basis(END)
    mass_q = np.linalg.solve(factor.T, mass_basis)
    H = np.sqrt(2*X)*E_new
    current = moment_vector(U, E_new, X, weights)
    linear = np.vstack((mass_q, Q@(weights*H)))
    required = np.array([0., target[2]-current[2]])
    unconstrained = -(Q@(weights*U))
    dual = linear@linear.T
    minimum = unconstrained+linear.T@np.linalg.solve(
        dual, required-linear@unconstrained)
    at_minimum = moment_vector(U+minimum@Q, E_new, X, weights)
    gap = float(target[3]-at_minimum[3])
    if gap < 0:
        raise RuntimeError(f'No finite-basis S capacity: {gap}')
    options = []
    for index in range(len(Q)):
        trial = np.eye(len(Q))[index]
        null = trial-linear.T@np.linalg.solve(dual, linear@trial)
        norm = float(null@null)
        if norm < 1e-12:
            continue
        for sign in (-1., 1.):
            cq = minimum+sign*np.sqrt(gap/norm)*null
            U_new = U+cq@Q
            options.append((float(np.max(np.abs(U_new))),cq,index,sign))
    _, cq, index, sign = min(options, key=lambda item:item[0])
    cb = np.linalg.solve(factor, cq)
    achieved = moment_vector(U+cb@B, E_new, X, weights)
    field = RemoteUPatchedField(mean, base, e_intervals,
                                e_coefficients, cb)
    report = dict(moment_defect=(achieved-target).tolist(),
                  S_slack_before_completion=gap,
                  selected_null_index=index, selected_null_sign=sign,
                  u_coefficients=cb.tolist(),
                  e_coefficients=e_coefficients.tolist(),
                  min_relative_E=float(np.min(E_new/E0)),
                  exterior_streamfunction_increment=float(mass_basis@cb),
                  quadrature_mass_increment=float(weights@(cb@B)),
                  max_abs_u_coefficient=float(np.max(np.abs(cb))))
    return field, report


def run(smooth=False):
    global START, WIDTH, DEGREE
    if smooth:
        START, WIDTH, DEGREE = 1.01, .02, 11
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    floor_source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    e_row = next(entry['result'] for entry in floor_source['rows']
                 if entry['relative_E_floor'] == .05)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target_field = make_field(16, 2.)
    field, slice_report = solve_coefficients(mean, base, target_field,
                                             e_source, e_row)
    tau = .5*2**(-5.5)
    screens = []
    center_xs = ((1.015, 1.02, 1.03, 1.04, 1.3, 2.0) if smooth
                 else (1.031, 1.04, 1.3, 2.0))
    for name, xs, etas, t in (
        ('center', center_xs, (.3,), tau),
        ('space_holdout', (1.034, 1.1, 2.4), (.28, .32), tau),
        ('time_holdout', (1.035, 1.5, 3.2), (.29, .31),
         .5*2**(-5.4)),
    ):
        points, Xh, etah = nodes(base, xs, etas, t)
        before, _ = residual(mean, points, t)
        after, divergence = residual(field, points, t)
        screens.append(dict(name=name, X=Xh.tolist(), eta=etah.tolist(),
                            before=stats(before), after=stats(after),
                            max_abs_fd_divergence=float(
                                np.max(np.abs(divergence)))))
    cone_points, _, _ = nodes(base, (1.008,1.016,1.02,1.03),
                              (.2,.25,.3), tau)
    old_local, _ = mean.fields(cone_points, tau)
    new_local, _ = field.fields(cone_points, tau)
    report = dict(source='delayed_remote_swirl_floor_capacity.json',
                  start=START,end=END,width=WIDTH,degree=DEGREE,
                  protected_inner_cone=not smooth,
                  axial_interval=AXIAL,slice=slice_report,
                  screens=screens,
                  max_abs_local_velocity_change=float(np.max(np.abs(
                      new_local-old_local))),
                  scope='Exact-curl U plus swirl E lift. '
                        'Fixed-slice J/S and exterior streamfunction '
                        'are solved; M defect is explicit. Sampled '
                        'full Cartesian momentum and FD divergence '
                        'do not prove a continuous or global PDE gate. '
                        'The smooth variant enters the old stress cone.',
                  accepted=False)
    output = ('delayed_remote_u_physical_smooth.json' if smooth
              else 'delayed_remote_u_physical.json')
    (ROOT/output).write_bytes(
        (json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(slice=slice_report,
                          screens=screens,
                          max_abs_local_velocity_change=(
                              report['max_abs_local_velocity_change']))),
          flush=True)


if __name__ == '__main__':
    import sys
    run(smooth='--smooth' in sys.argv[1:])
