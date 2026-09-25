"""Fit compact swirl modes to the pressure-invariant angular equation.

This is a velocity-correction diagnostic for the curvature-optimized lift.
It leaves meridional velocity and divergence unchanged but does not enforce
the five radial moments after the swirl correction.
"""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from curvature_pressure_screen import components, sample
from joined_field import coordinates, independent_fd
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import septic_step


INTERVALS = ((1., 1.12), (1., 1.5), (1., 3.), (2., 3.))


def swirl_basis(radius, z, tau, field):
    radius, z, tau = np.broadcast_arrays(radius, z, tau)
    sn = np.sqrt(field.nu)
    co = coordinates(radius/sn, z/sn, tau, field.heat.h)
    X, eta, q = (np.asarray(co[key]) for key in ('X', 'eta', 'q'))
    rise, _ = septic_step((eta-.1)/.1)
    fall, _ = septic_step((eta-.3)/.1)
    cut = rise*(1-fall)
    axial = (cut, cut*(eta-.25)/.05)
    values = [sn*q**(-field.A)*bump(X, lo, hi)[0]*ax
              for lo, hi in INTERVALS for ax in axial]
    return np.stack(values, axis=1)


class SwirlTransportField:
    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes, float)

    def fields(self, points, tau):
        points = np.asarray(points, float)
        velocity, pressure = self.base.fields(points, tau)
        radius = np.hypot(points[:, 0], points[:, 1])
        delta = swirl_basis(radius, points[:, 2], tau, self.base)@self.amplitudes
        safe = np.maximum(radius, 1e-300)
        velocity[:, 0] -= delta*points[:, 1]/safe
        velocity[:, 1] += delta*points[:, 0]/safe
        return velocity, pressure


def angular_operator(base, points, tau, h, ht):
    """Fourth-order cylindrical transport-diffusion for each swirl basis."""
    r = np.hypot(points[:, 0], points[:, 1])
    z = points[:, 2]
    u, _ = base.fields(points, tau)
    ur = (u[:, 0]*points[:, 0]+u[:, 1]*points[:, 1])/r
    uz = u[:, 2]

    def at(dr=0., dz=0., dt=0.):
        return swirl_basis(r+dr, z+dz, tau+dt, base)

    f = at()
    rm2, rm1, rp1, rp2 = at(dr=-2*h), at(dr=-h), at(dr=h), at(dr=2*h)
    zm2, zm1, zp1, zp2 = at(dz=-2*h), at(dz=-h), at(dz=h), at(dz=2*h)
    tm2, tm1, tp1, tp2 = at(dt=-2*ht), at(dt=-ht), at(dt=ht), at(dt=2*ht)
    fr = (rm2-8*rm1+8*rp1-rp2)/(12*h)
    fz = (zm2-8*zm1+8*zp1-zp2)/(12*h)
    frr = (-rp2+16*rp1-30*f+16*rm1-rm2)/(12*h*h)
    fzz = (-zp2+16*zp1-30*f+16*zm1-zm2)/(12*h*h)
    ft = -(tm2-8*tm1+8*tp1-tp2)/(12*ht)
    return (ft+ur[:, None]*fr+uz[:, None]*fz+ur[:, None]*f/r[:, None]
            -base.nu*(frr+fr/r[:, None]+fzz-f/r[:, None]**2))


def angular_residual(residual, points):
    return components(residual, points)[1]


def run():
    base = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    tau = .5*2**(-5.5)
    grid = (
        ('train', (1.01, 1.025, 1.05, 1.1, 1.25, 1.5, 2., 2.5,
                   2.975, 2.99), (.2, .25, .3), tau),
        ('space_holdout', (1.015, 1.075, 1.4, 2.75, 2.985),
         (.225, .275), tau),
        ('time_holdout', (1.015, 1.075, 1.4, 2.75, 2.985),
         (.225, .275), .5*2**(-5.25)))
    data = {}
    for label, xs, es, t in grid:
        points, _, _ = sample(base, xs, es, t)
        h, ht = .001*np.sqrt(base.nu*t), .00025*t
        residual, divergence = independent_fd(base, points, t, h, ht)
        operator = angular_operator(base, points, t, h, ht)
        data[label] = (points, t, residual, divergence, operator)

    train = data['train']
    matrix = train[4]
    target = -angular_residual(train[2], train[0])
    norms = np.maximum(np.linalg.norm(matrix, axis=0), 1e-300)
    trials = []
    for ridge in (.001, .01, .1, 1., 10.):
        system = np.vstack((matrix/norms, ridge*np.eye(matrix.shape[1])))
        rhs = np.r_[target, np.zeros(matrix.shape[1])]
        coefficients = np.linalg.lstsq(system, rhs, rcond=None)[0]/norms
        rows = []
        for label, (points, _, residual, _, op) in data.items():
            before = angular_residual(residual, points)
            after = before+op@coefficients
            rows.append(dict(name=label,
                             before_max=float(np.max(np.abs(before))),
                             after_max=float(np.max(np.abs(after))),
                             before_rms=float(np.sqrt(np.mean(before**2))),
                             after_rms=float(np.sqrt(np.mean(after**2)))))
        trials.append((ridge, coefficients, rows))
    # Select on independent space and nearby-time samples; retain all trials.
    chosen = min(trials, key=lambda item: max(item[2][1]['after_max'],
                                              item[2][2]['after_max']))
    field = SwirlTransportField(base, chosen[1])
    verification = []
    for label in ('train', 'space_holdout', 'time_holdout'):
        points, t, residual, divergence, operator = data[label]
        selected = points[:2]
        direct, direct_div = independent_fd(
            field, selected, t, .001*np.sqrt(base.nu*t), .00025*t)
        predicted = angular_residual(residual[:2], selected)+operator[:2]@chosen[1]
        verification.append(dict(name=label,
                                 direct_angular_max=float(np.max(np.abs(
                                     angular_residual(direct, selected)))),
                                 angular_model_fd_discrepancy=float(np.max(np.abs(
                                     angular_residual(direct, selected)-predicted))),
                                 direct_full_max=float(np.max(np.linalg.norm(
                                     direct, axis=1))),
                                 baseline_full_max=float(np.max(np.linalg.norm(
                                     residual[:2], axis=1))),
                                 divergence_max=float(np.max(np.abs(direct_div)))))
    report = dict(intervals=INTERVALS, axial_support=[.1, .4],
                  selected_ridge=chosen[0],
                  selected_amplitudes=chosen[1].tolist(),
                  trials=[dict(ridge=r, rows=rows) for r, _, rows in trials],
                  verification=verification,
                  scope='Sampled angular-equation fit only. Swirl correction '
                        'changes E/I/Cp moments and radial momentum; no '
                        'five-moment or global PDE admission.', accepted=False)
    (ROOT/'curvature_swirl_transport_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(selected_ridge=chosen[0],
                          selected_rows=chosen[2],
                          verification=verification)), flush=True)


if __name__ == '__main__':
    run()
