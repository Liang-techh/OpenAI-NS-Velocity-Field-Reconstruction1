"""Locate the remaining inner-collar force and local transport scale."""
import json
import numpy as np

from heated_interior_join import HeatedInteriorJoin
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def run():
    field = HeatedInteriorJoin()
    tau = .5*2**(-5.5)
    r1 = float(field.attachment_radius(tau))
    eta, fraction = np.meshgrid([.2725, .345, .4175],
                                [.2, .4, .6, .8], indexing='ij')
    eta, fraction = eta.ravel(), fraction.ravel()
    q = tau/(1-eta**2)
    z = np.sqrt(field.nu)*q**(.5-field.heat.h)*eta
    radius = fraction*r1
    points = np.column_stack((radius, np.zeros(len(radius)), z))
    u, J, part, terms, visc_axes = kinematics(
        field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau,
        return_terms=True, return_viscosity_axes=True)
    R = part+np.einsum('nij,nj->ni', J, u)
    rows = []
    for i in range(len(points)):
        transport_gradient = np.array([J[i, 1, 0]+u[i, 1]/radius[i],
                                       J[i, 1, 2]])
        den = np.dot(transport_gradient, transport_gradient)
        minimal_transport_delta = (-R[i, 1]*transport_gradient/den
                                   if den > 1e-20 else np.array([np.nan, np.nan]))
        rows.append(dict(eta=float(eta[i]), radius_fraction=float(fraction[i]),
                         radius=float(radius[i]),
                         velocity=u[i].tolist(), residual=R[i].tolist(),
                         residual_norm=float(np.linalg.norm(R[i])),
                         time=terms['time'][i].tolist(),
                         convection=terms['convection'][i].tolist(),
                         pressure=terms['pressure'][i].tolist(),
                         viscosity=terms['viscosity'][i].tolist(),
                         viscosity_z=visc_axes[2, i].tolist(),
                         swirl_transport_gradient=transport_gradient.tolist(),
                         minimal_meridional_velocity_delta_for_angular=\
                             minimal_transport_delta.tolist()))
    rows.sort(key=lambda row: row['residual_norm'], reverse=True)
    report = dict(tau=tau, rows=rows,
                  scope='Local linear transport target only; adding it pointwise does not solve coupled momentum, incompressibility, or boundary matching.',
                  accepted=False)
    (ROOT/'inner_collar_targets.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(rows[:4], indent=2))


if __name__ == '__main__':
    run()
