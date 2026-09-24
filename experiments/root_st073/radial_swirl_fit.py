"""Radial-cutoff swirl modes fitted against the full joint-collar momentum."""
import json

import numpy as np

from compact_control import cutoff
from compact_potential import CompactPotentialField
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from joint_collar_scale import ScaledJointCollarField


class RadialSwirlMode:
    def __init__(self, base, mode, reference_tau, temporal_power=1.5):
        self.base = base
        self.nu = base.nu
        self.mode = mode
        self.reference_tau = reference_tau
        self.temporal_power = temporal_power

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        for i, (point, t) in enumerate(zip(pts, ts)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            rflat, rsupp, zflat, zsupp = self.base.support(t)
            if r <= rflat or r >= rsupp or abs(z) >= zsupp:
                continue
            y = (r-rflat)/(rsupp-rflat)
            radial = 1024*y**5*(1-y)**5
            if self.mode:
                radial *= (2*y-1)
            axial = float(cutoff((abs(z)-zflat)/(zsupp-zflat))[0])
            theta = ((self.reference_tau/t)**self.temporal_power
                     * radial*axial)
            velocity[i] = [-theta*ycart/r, theta*x/r, 0.]
        return velocity, pressure


class RadialSwirlRepairedField:
    """Callable sum of joint collar field and two radial-collar swirl modes."""

    def __init__(self, joint, amplitudes, reference_tau):
        self.joint = joint
        self.base = joint.base
        self.nu = joint.nu
        self.amplitudes = np.asarray(amplitudes, float)
        self.modes = [RadialSwirlMode(self.base, i, reference_tau)
                      for i in range(2)]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        velocity, pressure = self.joint.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                du, _ = mode.fields(points, tau)
                velocity += amplitude*du
        return velocity, pressure


def linear_screen(field, modes, tau, order):
    points, weights = nodes(field.base, tau, order)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    baseline = part + np.einsum('nij,nj->ni', grad, u)
    columns, changes, changes_grad = [], [], []
    for mode in modes:
        du, dgrad, dpart = kinematics(mode, points, tau, hs, ht)
        columns.append(dpart + np.einsum('nij,nj->ni', dgrad, u)
                       + np.einsum('nij,nj->ni', grad, du))
        changes.append(du)
        changes_grad.append(dgrad)
    return (points, weights, baseline, np.stack(columns, axis=-1),
            np.stack(changes, axis=-1), np.stack(changes_grad, axis=-1))


def metrics(data, amplitudes):
    _, weights, baseline, columns, changes, changes_grad = data
    du = np.einsum('nik,k->ni', changes, amplitudes)
    dgrad = np.einsum('nijk,k->nij', changes_grad, amplitudes)
    residual = (baseline + np.einsum('nik,k->ni', columns, amplitudes)
                + np.einsum('nij,nj->ni', dgrad, du))
    norms = np.linalg.norm(residual, axis=1)
    return {'max': float(np.max(norms)),
            'physical_volume_l2': float(np.sqrt(np.dot(weights, norms**2))),
            'max_radial': float(np.max(np.abs(residual[:, 0]))),
            'max_angular': float(np.max(np.abs(residual[:, 1]))),
            'max_axial': float(np.max(np.abs(residual[:, 2])))}


def run():
    fit = json.loads((ROOT/'compact_potential'/'joint_collar_fit.json').read_text())
    base = CompactPotentialField()
    reference_tau = fit['tau']
    joint = ScaledJointCollarField(.1*np.array(fit['linear_fit_amplitudes']),
                                   reference_tau, base=base)
    modes = [RadialSwirlMode(base, i, reference_tau) for i in range(2)]
    train_tau = .5/64
    train = linear_screen(joint, modes, train_tau, 8)
    _, weights, residual, columns, _, _ = train
    matrix = (np.sqrt(weights)[:, None, None]*columns).reshape(-1, 2)
    target = (-np.sqrt(weights)[:, None]*residual).reshape(-1)
    norms = np.linalg.norm(matrix, axis=0)
    scaled, *_ = np.linalg.lstsq(matrix/norms, target, rcond=1e-8)
    direction = scaled/norms
    rows = []
    for k, order in ((6., 8), (6., 6), (5.5, 6), (4., 6), (.4, 6)):
        tau = .5*2**(-k)
        data = train if (k == 6. and order == 8) else linear_screen(joint, modes, tau, order)
        for strength in (0., .25, .5, 1.):
            row = {'k': k, 'order': order, 'strength': strength,
                   'metrics': metrics(data, strength*direction)}
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = {'reference_tau': reference_tau,
              'training_k': 6., 'linear_fit_amplitudes': direction.tolist(),
              'rows': rows,
              'scope': 'Two C4 radial-collar pure-swirl modes on top of the scale-weighted joint axial-collar candidate. Full nonlinear physical momentum, Gauss8 fit and Gauss6 time holdouts; no converged norm or pulse realization.',
              'accepted': False}
    out = ROOT/'compact_potential'/'radial_swirl_fit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
