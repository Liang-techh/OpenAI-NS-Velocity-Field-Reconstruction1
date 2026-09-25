"""Test a restricted radial Reynolds-stress closure on the ST073 bridge.

If only Q_{r theta} and Q_{r z} carry the correction and both vanish at
the bridge endpoints, cylindrical divergence requires the weighted radial
integrals of the complete mean residual to vanish. This is a necessary
condition for that restricted ansatz, not for the paper's full wave system.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from high_frequency_shear_screen import make_field
from joined_field import independent_fd
from radial_continuation import ROOT


def flux_row(field, eta, tau, order=16):
    radial = field.compact.joined
    xi, xo = radial.join_X, radial.join_X*radial.outer_ratio**2
    q = tau/(1-eta**2)
    ri, ro = np.sqrt(2*field.nu*q*np.array([xi, xo]))
    g, w = leggauss(order)
    radii = (ri+ro)/2+(ro-ri)*g/2
    weights = (ro-ri)*w/2
    X = radii**2/(2*field.nu*q)
    points = radial.inner.from_similarity(
        X, np.full(order, eta), tau)
    R, divergence = independent_fd(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
    theta_weight = weights*radii**2
    axial_weight = weights*radii
    theta_integral = float(theta_weight@R[:, 1])
    axial_integral = float(axial_weight@R[:, 2])
    theta_total = float(theta_weight@np.abs(R[:, 1]))
    axial_total = float(axial_weight@np.abs(R[:, 2]))
    return dict(eta=eta, order=order, ri=float(ri), ro=float(ro),
                theta_weighted_integral=theta_integral,
                axial_weighted_integral=axial_integral,
                theta_cancellation_ratio=(abs(theta_integral)/theta_total
                                          if theta_total else 0.),
                axial_cancellation_ratio=(abs(axial_integral)/axial_total
                                          if axial_total else 0.),
                required_outer_Q_rtheta=-theta_integral/ro**2,
                required_outer_Q_rz=-axial_integral/ro,
                residual_max=float(np.max(np.linalg.norm(R, axis=1))),
                divergence_max=float(np.max(np.abs(divergence))))


def run():
    tau = .5*2**(-5.5)
    rows = []
    field = make_field(16, 0.)
    for order in (16, 32):
        for eta in (.2, .3):
            row = dict(field='balanced', **flux_row(field, eta, tau,
                                                   order=order))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  restricted_ansatz='Only radial-tangential Q_rtheta and Q_rz fluxes, zero at both bridge edges. Their divergence in cylindrical coordinates is ((d/dr+2/r)Q_rtheta,(d/dr+1/r)Q_rz).',
                  scope='Gauss16/Gauss32 weighted physical-residual diagnostic at two axial slices of the balanced mean. Nonzero closure integrals obstruct only this radial-flux-only supported stress ansatz. They do not rule out the paper’s other stress components, axial derivatives, mean corrections, or waves.',
                  accepted=False)
    (ROOT/'radial_flux_closure_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
