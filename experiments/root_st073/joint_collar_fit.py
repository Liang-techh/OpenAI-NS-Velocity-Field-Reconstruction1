"""Six-mode solenoidal velocity/pressure screen for the compact axial collar.

This is a full physical NS residual fit, including the quadratic correction
transport. It is a finite-sample candidate screen, not a paper wave solution.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from compact_potential import CompactPotentialField
from joined_field import ROOT


class CollarMode:
    names = ('poloidal_even', 'poloidal_odd', 'swirl_even', 'swirl_odd',
             'pressure_even', 'pressure_odd')

    def __init__(self, base, mode, temporal_power=0., reference_tau=None):
        self.base = base
        self.mode = mode
        self.nu = base.nu
        self.temporal_power = temporal_power
        self.reference_tau = reference_tau

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        for i, (point, t) in enumerate(zip(pts, ts)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            _, radius, zflat, zsupp = self.base.support(t)
            if r >= radius or abs(z) <= zflat or abs(z) >= zsupp:
                continue
            y = r / radius
            s = (abs(z) - zflat) / (zsupp - zflat)
            axial = 1024 * s**5 * (1 - s)**5
            axial_derivative = (5120 * s**4 * (1 - s)**4 * (1 - 2 * s)
                                / (zsupp - zflat))
            if self.mode % 2:
                axial *= np.sign(z)
            else:
                axial_derivative *= np.sign(z)
            even_radial = (1 - y*y)**5
            if self.mode < 2:
                # psi=R^2 * 15*y^2*(1-y^2)^5 * axial.
                # u_r=-psi_z/r and u_z=psi_r/r, regular at the axis.
                ur = -15 * radius * y * even_radial * axial_derivative
                uz = 30 * (1 - y*y)**4 * (1 - 6*y*y) * axial
                velocity[i] = [ur*x/r if r else 0., ur*ycart/r if r else 0., uz]
            elif self.mode < 4:
                theta = 6 * y * even_radial * axial
                velocity[i] = [-theta*ycart/r if r else 0., theta*x/r if r else 0., 0.]
            else:
                pressure[i] = even_radial * axial
        if self.reference_tau is not None:
            factors = (self.reference_tau / ts)**self.temporal_power
            velocity *= factors[:, None]
            pressure *= factors
        return velocity, pressure


class TransitionPoloidalMode:
    """Compact solenoidal axial-shape bubble in the existing collar."""

    def __init__(self, base, temporal_power=2., reference_tau=None):
        self.base = base
        self.nu = base.nu
        self.temporal_power = temporal_power
        self.reference_tau = reference_tau

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        for i, (point, t) in enumerate(zip(pts, ts)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            _, radius, zflat, zsupp = self.base.support(t)
            if r >= radius or abs(z) <= zflat or abs(z) >= zsupp:
                continue
            y = r/radius
            s = (abs(z)-zflat)/(zsupp-zflat)
            bubble = 1024*s**5*(1-s)**5
            bubble_prime = 5120*s**4*(1-s)**4*(1-2*s)
            axial = bubble*(1+4*s)
            axial_z = ((bubble_prime*(1+4*s)+4*bubble)
                       * np.sign(z)/(zsupp-zflat))
            radial = (1-y*y)**5
            ur = -15*radius*y*radial*axial_z
            uz = 30*(1-y*y)**4*(1-6*y*y)*axial
            velocity[i] = [ur*x/r if r else 0.,
                           ur*ycart/r if r else 0., uz]
        if self.reference_tau is not None:
            velocity *= ((self.reference_tau/ts)**self.temporal_power)[:, None]
        return velocity, np.zeros(len(pts))


def kinematics(field, points, tau, hs, ht, time_min=.5/64,
               time_max=.5):
    u, p = field.fields(points, tau)
    grad = np.zeros((len(points), 3, 3))
    gp = np.zeros((len(points), 3))
    lap = np.zeros_like(u)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = hs
        um2, pm2 = field.fields(points - 2*step, tau)
        um1, pm1 = field.fields(points - step, tau)
        up1, pp1 = field.fields(points + step, tau)
        up2, pp2 = field.fields(points + 2*step, tau)
        grad[:, :, axis] = (um2 - 8*um1 + 8*up1 - up2) / (12*hs)
        gp[:, axis] = (pm2 - 8*pm1 + 8*pp1 - pp2) / (12*hs)
        lap += (-up2 + 16*up1 - 30*u + 16*um1 - um2) / (12*hs*hs)
    if tau - 2*ht < time_min:
        utau = (-25*u + 48*field.fields(points, tau+ht)[0]
                - 36*field.fields(points, tau+2*ht)[0]
                + 16*field.fields(points, tau+3*ht)[0]
                - 3*field.fields(points, tau+4*ht)[0]) / (12*ht)
    elif tau + 2*ht > time_max:
        utau = (25*u - 48*field.fields(points, tau-ht)[0]
                + 36*field.fields(points, tau-2*ht)[0]
                - 16*field.fields(points, tau-3*ht)[0]
                + 3*field.fields(points, tau-4*ht)[0]) / (12*ht)
    else:
        utau = (field.fields(points, tau-2*ht)[0]
                - 8*field.fields(points, tau-ht)[0]
                + 8*field.fields(points, tau+ht)[0]
                - field.fields(points, tau+2*ht)[0]) / (12*ht)
    return u, grad, -utau + gp - field.nu*lap


def nodes(base, tau, order):
    _, radius, _, height = base.support(tau)
    g, w = leggauss(order)
    rs = radius * (g + 1) / 2
    zs = height * g
    points = np.array([[r, 0., z] for z in zs for r in rs])
    weights = (2*np.pi*np.tile(rs, order)
               * np.tile(w*radius/2, order)
               * np.repeat(w*height, order))
    return points, weights


def screen(base, tau, order, temporal_power=0., reference_tau=None):
    points, weights = nodes(base, tau, order)
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    u, grad, linear = kinematics(base, points, tau, hs, ht)
    residual0 = linear + np.einsum('nij,nj->ni', grad, u)
    cols = []
    changes = []
    changes_grad = []
    for mode in range(6):
        mode_field = CollarMode(base, mode, temporal_power, reference_tau)
        du, dgrad, dlinear = kinematics(mode_field, points, tau, hs, ht)
        column = (dlinear + np.einsum('nij,nj->ni', dgrad, u)
                  + np.einsum('nij,nj->ni', grad, du))
        cols.append(column)
        changes.append(du)
        changes_grad.append(dgrad)
    return points, weights, residual0, np.stack(cols, axis=-1), np.stack(changes, axis=-1), np.stack(changes_grad, axis=-1)


def metrics(screen_data, amplitudes):
    _, weights, residual0, columns, changes, changes_grad = screen_data
    du = np.einsum('nik,k->ni', changes, amplitudes)
    dgrad = np.einsum('nijk,k->nij', changes_grad, amplitudes)
    residual = (residual0 + np.einsum('nik,k->ni', columns, amplitudes)
                + np.einsum('nij,nj->ni', dgrad, du))
    norm = np.linalg.norm(residual, axis=1)
    return {'max': float(np.max(norm)),
            'physical_volume_l2': float(np.sqrt(np.dot(weights, norm**2))),
            'max_radial': float(np.max(np.abs(residual[:, 0]))),
            'max_angular': float(np.max(np.abs(residual[:, 1]))),
            'max_axial': float(np.max(np.abs(residual[:, 2])))}


def run():
    base = CompactPotentialField()
    tau = .5*2**(-5.5)
    train = screen(base, tau, 5)
    _, weights, residual, columns, _, _ = train
    matrix = (np.sqrt(weights)[:, None, None]*columns).reshape(-1, 6)
    target = (-np.sqrt(weights)[:, None]*residual).reshape(-1)
    norms = np.linalg.norm(matrix, axis=0)
    scaled, *_ = np.linalg.lstsq(matrix/norms, target, rcond=1e-8)
    direction = scaled/norms
    holdout = screen(base, tau, 6)
    rows = []
    for strength in (0., .01, .03, .1, .3, 1.):
        amp = strength*direction
        rows.append({'strength': strength,
                     'train': metrics(train, amp),
                     'holdout': metrics(holdout, amp)})
        print(json.dumps(rows[-1]), flush=True)
    report = {'tau': tau, 'mode_names': CollarMode.names,
              'linear_fit_amplitudes': direction.tolist(),
              'rows': rows,
              'scope': 'Axisymmetric compact poloidal/swirl/pressure correction in axial cutoff collar, training Gauss5 and held-out Gauss6. Full nonlinear Cartesian momentum at each sampled point; no nonaxisymmetric stress or critical-time proof.',
              'accepted': False}
    out = ROOT/'compact_potential'/'joint_collar_fit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
