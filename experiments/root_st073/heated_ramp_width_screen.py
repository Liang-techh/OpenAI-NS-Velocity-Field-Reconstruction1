"""Diagnose whether radial ramp width resolves the collar swirl equation."""
import json

import numpy as np

from axial_compact_join import load_compact_candidate
from heated_global_join import HeatedGlobalJoin
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def run():
    compact = load_compact_candidate()
    k = 5.5
    tau = .5*2**(-k)
    eta = .345
    q = tau/(1-eta**2)
    z = np.sqrt(compact.nu)*q**(.5-compact.joined.inner.h)*eta
    fractions = np.array([.25, .5, .75])
    rows = []
    for outer_ratio in (1.5, 2., 4., 8.):
        field = HeatedGlobalJoin(compact=compact, outer_ratio=outer_ratio)
        r1, r2 = field.ramp_radii(tau)
        radius = r1+fractions*(r2-r1)
        points = np.column_stack((radius, np.zeros(len(radius)),
                                  np.full(len(radius), z)))
        u, grad, linear, terms = kinematics(
            field, points, tau, .001*np.sqrt(field.nu*tau),
            .00025*tau, return_terms=True)
        residual = linear+np.einsum('nij,nj->ni', grad, u)
        rows.append(dict(outer_ratio=outer_ratio, r1=float(r1),
                         r2=float(r2), fractions=fractions.tolist(),
                         radii=radius.tolist(),
                         complete_norms=np.linalg.norm(residual, axis=1).tolist(),
                         angular_residual=residual[:, 1].tolist(),
                         angular_time=terms['time'][:, 1].tolist(),
                         angular_transport=terms['convection'][:, 1].tolist(),
                         angular_viscosity=terms['viscosity'][:, 1].tolist()))
        print(json.dumps(rows[-1]), flush=True)
    report = dict(k=k, tau=tau, eta=eta, rows=rows,
                  scope='Twelve sampled points in the radial ramp at one time and axial similarity height; diagnostic only.',
                  accepted=False)
    (ROOT/'heated_ramp_width_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
