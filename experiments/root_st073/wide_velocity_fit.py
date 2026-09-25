"""Coupled solenoidal meridional/swirl trial on the wide ST073 transition."""
import json

import numpy as np
from scipy.optimize import least_squares

from joined_field import coordinates, independent_fd
from joint_collar_fit import kinematics
from radial_continuation import ROOT
from wide_pressure_fit import load_pressure_candidate, points_at


class BridgeVelocityMode:
    def __init__(self, joined, kind, radial_degree, axial_degree):
        self.joined = joined
        self.kind = kind
        self.radial_degree = radial_degree
        self.axial_degree = axial_degree
        self.nu = joined.nu

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        sn = np.sqrt(self.nu)
        r = np.hypot(pts[:, 0], pts[:, 1])
        co = coordinates(r/sn, pts[:, 2]/sn, tau, self.joined.inner.h)
        X, eta, q = (np.asarray(co[name]) for name in ('X', 'eta', 'q'))
        ratio = self.joined.outer_ratio
        y = (np.sqrt(X/self.joined.join_X)-1)/(ratio-1)
        active = (y > 0) & (y < 1) & (np.abs(eta) < .5)
        yy = np.clip(y, 0, 1)
        B = 1024*yy**5*(1-yy)**5 * active
        By = 5120*yy**4*(1-yy)**4*(1-2*yy) * active
        R = np.ones_like(y) if self.radial_degree == 0 else 2*y-1
        Ry = np.zeros_like(y) if self.radial_degree == 0 else 2*np.ones_like(y)
        s = np.clip(eta/.5, -1, 1)
        taper = (1-s*s)**4
        taper_e = -16*s*(1-s*s)**3
        if self.axial_degree == 0:
            T, Te = taper, taper_e
        else:
            T, Te = taper*s, taper_e*s + 2*taper
        shape = B*R*T
        sy = (By*R+B*Ry)*T
        se = B*R*Te
        safe_X = np.maximum(X, 1e-300)
        yr = np.asarray(co['X_r'])/(sn*2*(ratio-1)*np.sqrt(safe_X*self.joined.join_X))
        yz = np.asarray(co['X_z'])/(sn*2*(ratio-1)*np.sqrt(safe_X*self.joined.join_X))
        ez = np.asarray(co['eta_z'])/sn
        qz = np.asarray(co['q_z'])/sn
        A = .5+self.joined.inner.h
        ur = np.zeros_like(r)
        uz = np.zeros_like(r)
        uth = np.zeros_like(r)
        if self.kind == 'poloidal':
            scale = self.nu**1.5*q**(1-A)
            psi_r = scale*sy*yr
            psi_z = scale*(sy*yz+se*ez+(1-A)*qz/q*shape)
            ur = -psi_z/np.maximum(r, 1e-300)
            uz = psi_r/np.maximum(r, 1e-300)
        else:
            uth = sn*q**(-A)*shape
        cr = pts[:, 0]/np.maximum(r, 1e-300)
        sr = pts[:, 1]/np.maximum(r, 1e-300)
        velocity = np.column_stack((ur*cr-uth*sr,
                                    ur*sr+uth*cr, uz))
        return velocity, np.zeros(len(pts))


class VelocityCorrectedField:
    def __init__(self, base, modes, amplitudes):
        self.base = base
        self.modes = modes
        self.amplitudes = np.asarray(amplitudes, float)
        self.nu = base.nu

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            velocity += amplitude*mode.fields(points, tau)[0]
        return velocity, pressure


def make_modes(joined):
    return [BridgeVelocityMode(joined, kind, radial, axial)
            for kind in ('poloidal', 'swirl')
            for radial in (0, 1) for axial in (0, 1)]


def load_velocity_candidate():
    report = json.loads((ROOT/'wide_velocity_fit.json').read_text())
    base = load_pressure_candidate()
    return VelocityCorrectedField(base, make_modes(base.base),
                                  report['amplitudes'])


def run():
    base = load_pressure_candidate()
    joined = base.base
    modes = make_modes(joined)
    points, tau = points_at(joined, 5.5, [.2, .5, .8],
                            [-.3, -.15, 0, .15, .3])
    hs = .001*np.sqrt(joined.nu*tau)
    ht = .00025*tau
    u0, J0, part0 = kinematics(base, points, tau, hs, ht)
    pieces = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([item[0] for item in pieces], axis=-1)
    dJ = np.stack([item[1] for item in pieces], axis=-1)
    dpart = np.stack([item[2] for item in pieces], axis=-1)
    R0 = part0 + np.einsum('nij,nj->ni', J0, u0)
    scale = float(np.max(np.linalg.norm(R0, axis=1)))

    def trial(a):
        u = u0+np.einsum('nik,k->ni', du, a)
        J = J0+np.einsum('nijk,k->nij', dJ, a)
        part = part0+np.einsum('nik,k->ni', dpart, a)
        return part+np.einsum('nij,nj->ni', J, u)

    fit = least_squares(lambda a: np.r_[(trial(a)/scale).ravel(), .01*a],
                        np.zeros(len(modes)), bounds=(-2., 2.),
                        max_nfev=100, xtol=1e-9, ftol=1e-9, gtol=1e-9)
    candidate = VelocityCorrectedField(base, modes, fit.x)
    rows = []
    for name, k, ys, etas in (
            ('train', 5.5, [.2, .5, .8], [-.3, -.15, 0, .15, .3]),
            ('space_holdout', 5.5, [.3, .7], [-.25, .1, .25]),
            ('time_holdout', 5.25, [.3, .7], [-.25, .1, .25])):
        pts, t = points_at(joined, k, ys, etas)
        if name == 'train':
            before = R0
            after = trial(fit.x)
            div = None
        else:
            before = independent_fd(base, pts, t, .001*np.sqrt(joined.nu*t),
                                    .00025*t)[0]
            after, div = independent_fd(candidate, pts, t,
                                        .001*np.sqrt(joined.nu*t), .00025*t)
        row = dict(name=name, k=k, point_count=len(pts),
                   baseline_max=float(np.max(np.linalg.norm(before, axis=1))),
                   corrected_max=float(np.max(np.linalg.norm(after, axis=1))),
                   baseline_rms=float(np.sqrt(np.mean(before**2))),
                   corrected_rms=float(np.sqrt(np.mean(after**2))),
                   baseline_component_max=np.max(np.abs(before), axis=0).tolist(),
                   corrected_component_max=np.max(np.abs(after), axis=0).tolist(),
                   divergence_max=(None if div is None else
                                   float(np.max(np.abs(div)))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(join_X=joined.join_X, outer_ratio=joined.outer_ratio,
                  mode_ids=[dict(kind=m.kind, radial_degree=m.radial_degree,
                                 axial_degree=m.axial_degree) for m in modes],
                  amplitudes=fit.x.tolist(), optimizer_success=bool(fit.success),
                  optimizer_message=fit.message, rows=rows,
                  scope='Eight endpoint-preserving solenoidal velocity modes on '
                        'the pressure-corrected wide bridge. Fit at one late '
                        'time; independent spatial/time point checks. No axial '
                        'closure, exterior moment match, finite global energy '
                        'or complete-domain acceptance.', accepted=False)
    (ROOT/'wide_velocity_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
