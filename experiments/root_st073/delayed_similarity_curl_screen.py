"""Fit three compact similarity-scaled exact-curl velocity modes.

This probes whether a time-compatible solenoidal mean correction can absorb
the entrance momentum hotspot. Slice moments are audited before admission.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from joined_field import coordinates
from joint_collar_fit import kinematics
from paper_moment_bridge import bump
from radial_continuation import ROOT


INTERVALS = ((1.005, 1.06), (1.005, 1.12), (1.015, 1.18))
ETA_INTERVAL = (.15, .38)


class SimilarityCurlMode:
    def __init__(self, base, radial_interval, eta_interval=ETA_INTERVAL):
        self.base = base
        self.nu = base.nu
        self.radial_interval = radial_interval
        self.eta_interval = eta_interval

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        radius = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(radius/sn, pts[:, 2]/sn, ts,
                         self.base.heat.h)
        q, X, eta = (np.asarray(co[key]) for key in ('q', 'X', 'eta'))
        bx, bxd = bump(X, *self.radial_interval)
        be, bed = bump(eta, *self.eta_interval)
        phi = bx*be
        phi_x = bxd*be
        phi_eta = bx*bed
        qz = np.asarray(co['q_z'])/sn
        Xz = np.asarray(co['X_z'])/sn
        etaz = np.asarray(co['eta_z'])/sn
        A = self.base.A
        psi_z = self.nu**1.5*((1-A)*q**(-A)*qz*phi
                  +q**(1-A)*(phi_x*Xz+phi_eta*etaz))
        ur = np.divide(-psi_z, radius, out=np.zeros_like(radius),
                       where=radius > 0)
        uz = sn*q**(-A)*phi_x
        ca = np.divide(pts[:, 0], radius, out=np.ones_like(radius),
                       where=radius > 0)
        sa = np.divide(pts[:, 1], radius, out=np.zeros_like(radius),
                       where=radius > 0)
        velocity = np.column_stack((ur*ca, ur*sa, uz))
        return velocity, np.zeros(len(pts))


class CurlPatchedLift:
    def __init__(self, base, modes, coefficients):
        self.base = base
        self.modes = modes
        self.coefficients = np.asarray(coefficients, float)
        self.nu = base.nu
        self.heat = base.heat
        self.A = base.A
        self.compact = base.compact

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        for coefficient, mode in zip(self.coefficients, self.modes):
            if coefficient:
                velocity += coefficient*mode.fields(points, tau)[0]
        return velocity, pressure


def screen_data(field, modes, xs, etas, tau):
    points, X, eta = nodes(field, xs, etas, tau)
    hs = .001*np.sqrt(field.nu*tau)
    ht = .00025*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    base_residual = part+np.einsum('nij,nj->ni', grad, u)
    mode_data = [kinematics(mode, points, tau, hs, ht)
                 for mode in modes]
    columns = []
    for du, dgrad, dpart in mode_data:
        columns.append(dpart+np.einsum('nij,nj->ni', dgrad, u)
                       +np.einsum('nij,nj->ni', grad, du))
    return dict(points=points, X=X, eta=eta, tau=tau,
                u=u, grad=grad, base_residual=base_residual,
                mode_u=np.stack([item[0] for item in mode_data], axis=2),
                mode_grad=np.stack([item[1] for item in mode_data], axis=3),
                columns=np.stack(columns, axis=2))


def nonlinear_residual(data, coefficients):
    du = data['mode_u']@coefficients
    dgrad = data['mode_grad']@coefficients
    return (data['base_residual']+data['columns']@coefficients
            +np.einsum('nij,nj->ni', dgrad, du))


def run():
    base_name = 'delayed005_rise146_degree31.json'
    base = CoupledMomentPhysicalLift(slice_filename=base_name)
    modes = [SimilarityCurlMode(base, interval) for interval in INTERVALS]
    train = [screen_data(base, modes, (1.01, 1.015, 1.02),
                         (.2, .3), tau)
             for tau in (.5*2**(-5.5), .5*2**(-5.25))]
    R0 = np.concatenate([item['base_residual'] for item in train])
    C = np.concatenate([item['columns'] for item in train])
    matrix = C.reshape(-1, len(modes))
    norms = np.linalg.norm(matrix, axis=0)
    scaled = matrix/norms
    ridge = .03
    direction = np.linalg.solve(
        scaled.T@scaled+ridge**2*np.eye(len(modes)),
        -scaled.T@R0.reshape(-1))/norms
    trials = []
    for strength in (0., .01, .025, .05, .1, .2, .4, .7, 1.):
        coefficients = strength*direction
        R = np.concatenate([nonlinear_residual(item, coefficients)
                            for item in train])
        trials.append(dict(strength=strength,
                           coefficients=coefficients.tolist(),
                           metrics=stats(R)))
    selected = min(trials, key=lambda row: row['metrics']['max'])
    patched = CurlPatchedLift(base, modes, selected['coefficients'])
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(base, (1.0125, 1.0175, 1.025),
                             (.22, .28, .32), tau_holdout)
    before, _ = residual(base, points, tau_holdout)
    after, divergence = residual(patched, points, tau_holdout)
    report = dict(base_slice=base_name, radial_intervals=INTERVALS,
                  eta_interval=ETA_INTERVAL, ridge=ridge,
                  direction=direction.tolist(), trials=trials,
                  selected=selected,
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(),
                               before=stats(before), after=stats(after),
                               max_abs_divergence=float(np.max(np.abs(divergence))),
                               node_before=np.linalg.norm(before,
                                                          axis=1).tolist(),
                               node_after=np.linalg.norm(after,
                                                         axis=1).tolist()),
                  scope='Three compact similarity-scaled poloidal '
                        'streamfunction modes, exact curl. Two-time '
                        'train and disjoint space/time holdout, full '
                        'nonlinear Cartesian residual. Moment restoration '
                        'not yet imposed; no PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_similarity_curl_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(selected=selected,
                          holdout=report['holdout'])), flush=True)


if __name__ == '__main__':
    run()
