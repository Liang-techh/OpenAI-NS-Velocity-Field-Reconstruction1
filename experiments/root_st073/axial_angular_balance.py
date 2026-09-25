"""Resolve the pressure-independent swirl equation in the axial collar."""
import json

import numpy as np

from axial_compact_join import load_compact_candidate
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def run():
    field = load_compact_candidate()
    k = 5.5
    tau = .5*2**(-k)
    fractions = np.array([.25, .5, .75])
    X, s = np.meshgrid([.03, .5859375, 2.], fractions, indexing='ij')
    eta = field.eta_flat+(field.eta_outer-field.eta_flat)*s
    points = field.joined.inner.from_similarity(X.ravel(), eta.ravel(), tau)
    u, J, part, terms = kinematics(
        field, points, tau, .001*np.sqrt(field.nu*tau),
        .00025*tau, return_terms=True)
    residual = part+np.einsum('nij,nj->ni', J, u)
    rows = []
    for i in range(len(points)):
        row = dict(X=float(X.ravel()[i]), eta=float(eta.ravel()[i]),
                   angular_time=float(terms['time'][i, 1]),
                   angular_transport=float(terms['convection'][i, 1]),
                   angular_pressure=float(terms['pressure'][i, 1]),
                   angular_viscosity=float(terms['viscosity'][i, 1]),
                   angular_residual=float(residual[i, 1]),
                   complete_residual_norm=float(np.linalg.norm(residual[i])))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(k=k, tau=tau, rows=rows,
                  maximum_absolute_angular_residual=float(
                      np.max(np.abs(residual[:, 1]))),
                  maximum_absolute_by_term={name: float(np.max(np.abs(
                      term[:, 1]))) for name, term in terms.items()},
                  equation='For prescribed axisymmetric meridional flow, '
                           'u_theta obeys u_theta,t + u_r u_theta,r '
                           '+ u_z u_theta,z + (u_r/r)u_theta '
                           '= nu*(u_theta,rr+u_theta,r/r+u_theta,zz'
                           '-u_theta/r^2). Axisymmetric pressure has no '
                           'azimuthal gradient.',
                  scope='Fourth-order Cartesian/time term balance on nine '
                        'positive axial-collar points of the finite-energy '
                        'field. This defines a pressure-independent target '
                        'for a dynamically evolved swirl correction; no '
                        'full-domain or continuum bound.', accepted=False)
    (ROOT/'axial_angular_balance.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
