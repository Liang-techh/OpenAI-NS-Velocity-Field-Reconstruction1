"""Decompose the momentum hotspot produced by the coupled bridge bubble."""
import json
import numpy as np

from bridge_poloidal_mode import BridgePoloidalMode
from extended_compact_join import load_extended_heated_candidate
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def run():
    base = load_extended_heated_candidate()
    coupled = BridgePoloidalMode(
        load_extended_heated_candidate(outer_swirl_bubble_amplitude=1.1),
        -2.5, shape='curvature')
    tau = .5*2**(-5.5)
    point = base.compact.joined.inner.from_similarity(
        np.array([.5859375]), np.array([.2]), tau)
    rows = []
    for name, field in (('baseline', base), ('coupled', coupled)):
        u, J, part, terms, visc_axes = kinematics(
            field, point, tau, .001*np.sqrt(field.nu*tau), .00025*tau,
            return_terms=True, return_viscosity_axes=True)
        R = part+np.einsum('nij,nj->ni', J, u)
        rows.append(dict(name=name, velocity=u[0].tolist(),
                         residual=R[0].tolist(),
                         residual_norm=float(np.linalg.norm(R[0])),
                         terms={key: value[0].tolist()
                                for key, value in terms.items()},
                         viscous_axes_xyz=visc_axes[:, 0, :].tolist()))
    report = dict(tau=tau, X=.5859375, eta=.2, rows=rows,
                  scope='One finite-difference hotspot decomposition. The added poloidal mode raises axial radial-viscosity curvature; no momentum acceptance.',
                  accepted=False)
    (ROOT/'coupled_bridge_hotspot.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({row['name']: dict(residual_norm=row['residual_norm'],
                                          axial_viscosity=row['terms']['viscosity'][2],
                                          axial_z_viscosity=row['viscous_axes_xyz'][2][2])
                      for row in rows}))


if __name__ == '__main__':
    run()
