"""Callable outer-compensated reference and bounded annular momentum replay.

This is a candidate diagnostic, not an accepted global solution. Pressure is
unchanged and forcing is zero; the whole-reference scale generator is used.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import numpy as np
from grouped_joined_field import install_in_field
from scale_reference_candidate import load_reference
from scale_reference_outer_moment import OuterMomentNeutralizer
from scale_reference_actual_replay import _spatial_jets
from scale_transport_generator import solenoidal_scale_generator
from scale_generator_linearization import velocity_correction_response

ROOT = Path(__file__).resolve().parent


class OuterCompensatedReference:
    def __init__(self, reference, correction):
        self.reference, self.correction = reference, correction
        self.tau0, self.h, self.nu = reference.tau0, reference.h, reference.nu

    def velocity(self, points):
        return self.reference.velocity(points) + self.correction.velocity(points)

    def pressure(self, points):
        return self.reference.pressure(points)

    def fields(self, points, tau=None):
        u, p = self.reference.fields(points, tau)
        return u + self.correction.velocity(points), p


def load_compensated_reference():
    path = ROOT / 'scale_reference_outer_moment.json'
    report = json.loads(path.read_text())
    source = ROOT / report['sources']['reference_report']['path']
    if hashlib.sha256(source.read_bytes()).hexdigest() != report['sources']['reference_report']['sha256']:
        raise ValueError('Outer moment calibration belongs to another reference')
    reference = load_reference(step=source)
    return OuterCompensatedReference(reference, OuterMomentNeutralizer(report['parameters']['amplitude_C']))


def grid(correction, order):
    x, w = np.polynomial.legendre.leggauss(order)
    r = correction.radius_center + correction.radius_half_width*x
    z = correction.axial_center + correction.axial_half_width*x
    theta = np.arange(4)*np.pi/2
    rr, zz, tt = np.meshgrid(r, z, theta, indexing='ij')
    weights = (w[:, None, None]*w[None, :, None]*np.ones((1, 1, 4))*rr
               * correction.radius_half_width*correction.axial_half_width*np.pi/2)
    return np.column_stack(((rr*np.cos(tt)).ravel(), (rr*np.sin(tt)).ravel(), zz.ravel())), weights.ravel()


def metric(values, weights):
    return dict(volume_L2=float(np.sqrt(np.sum(weights[:, None]*values**2))),
                max_norm=float(np.linalg.norm(values, axis=1).max()))


def run():
    field = load_compensated_reference()
    base, correction = field.reference, field.correction
    install_in_field(base.base_field)
    rows = []
    for order in (12, 20):
        p, w = grid(correction, order)
        u, j, gp, lap = _spatial_jets(base, p, base.tau0, 1e-6)
        theta = np.arctan2(p[:, 1], p[:, 0])
        et = np.column_stack((-np.sin(theta), np.cos(theta), np.zeros(len(p))))
        swirl = np.sum(et*u, axis=1).reshape(-1, 4)
        mean = np.repeat(swirl.mean(axis=1), 4)
        residual = solenoidal_scale_generator(p, u, j, base.tau0, base.h, mean)
        residual += np.einsum('nij,nj->ni', j, u) - base.nu*lap + gp
        v, dv, lv = correction.velocity(p), correction.jacobian(p), correction.laplacian(p)
        linear, quadratic = velocity_correction_response(p, u, j, v, dv, lv,
            base.tau0, base.h, np.sum(et*v, axis=1), base.nu)
        delta = linear + quadratic
        cross = float(np.sum(w*np.sum(u*v, axis=1)))
        energy = float(np.sum(w*np.sum(v*v, axis=1))/2)
        row = dict(order=order, point_count=len(p), baseline=metric(residual, w),
                   compensated=metric(residual+delta, w), increment=metric(delta, w),
                   correction_self_energy=energy, energy_cross_term=cross,
                   total_energy_change=energy+cross,
                   swirl_ring_spread=float(np.max(np.ptp(swirl, axis=1))),
                   axial_moment_increment=float(np.sum(w*(p[:, 0]*v[:, 1]-p[:, 1]*v[:, 0]))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(status='completed', accepted=False, pde_validated=False,
                  scale_recursion_established=False, rows=rows,
                  sources={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                           for name in ('scale_reference_outer_moment.json', 'scale_reference_trust_nonlinear_fit.json')},
                  forcing='zero', pressure='unchanged from the frozen reference',
                  scope='Outer annulus only; base jets use spatial FD, correction jets are analytic. Four angular samples rely on disjointness from nonaxisymmetric wave support. This does not certify global residual, neutrality, or critical-time regularity.')
    (ROOT/'scale_reference_outer_replay.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    run()
