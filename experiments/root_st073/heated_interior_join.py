"""Move the heat-swirl attachment into the old radial join, experimentally.

The correction is purely azimuthal and therefore preserves divergence zero.
Its pressure correction is scalar. A flat radial ramp suppresses the singular
radial heat swirl at the axis and is exactly one by the old join radius.
"""
import json

import numpy as np

from axial_compact_join import load_compact_candidate
from axially_heated_exterior import AxiallyHeatedExterior
from heat_exterior import physical
from heated_global_join import smooth_ramp, HeatedGlobalJoin
from joined_field import coordinates, independent_fd
from radial_continuation import ROOT


class HeatedInteriorJoin:
    def __init__(self, compact=None, heat_age=32768.):
        self.compact = compact if compact is not None else load_compact_candidate()
        self.nu = self.compact.nu
        self.heat = AxiallyHeatedExterior(
            nu=self.nu, h=self.compact.joined.inner.h,
            amplitude=self.compact.joined.c, heat_age=heat_age)

    def attachment_radius(self, tau):
        qmax = np.asarray(tau)/(1-self.compact.eta_outer**2)
        return self.compact.joined.outer_ratio*np.sqrt(
            2*self.nu*qmax*self.compact.joined.join_X)

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.compact.fields(pts, ts)
        r = np.hypot(pts[:, 0], pts[:, 1])
        ramp = smooth_ramp(r/self.attachment_radius(ts))
        active = ramp > 0
        if not np.any(active):
            return velocity, pressure
        chosen = pts[active]
        radial = physical(chosen, ts[active], c=self.heat.amplitude,
                          h=self.heat.h, nu=self.nu)
        co = coordinates(r[active]/np.sqrt(self.nu),
                         chosen[:, 2]/np.sqrt(self.nu),
                         ts[active], self.heat.h)
        chi = self.compact.cutoff(co['eta'])[0]
        G, _ = self.heat.envelope(chosen[:, 2], ts[active])
        velocity[active] += (ramp[active]*(G-chi))[:, None]*radial['velocity']
        pressure[active] += ramp[active]*(G*G-chi)*radial['pressure']
        return velocity, pressure


def run():
    compact = load_compact_candidate()
    newer = HeatedInteriorJoin(compact=compact)
    older = HeatedGlobalJoin(compact=compact)
    k = 5.5
    tau = .5*2**(-k)
    r1 = float(newer.attachment_radius(tau))
    fractions = np.array([.25, .5, .75, 1., 1.25, 2.])
    rows = []
    for eta in (.2725, .345, .4175):
        q = tau/(1-eta**2)
        z = np.sqrt(newer.nu)*q**(.5-newer.heat.h)*eta
        points = np.column_stack((fractions*r1, np.zeros(len(fractions)),
                                  np.full(len(fractions), z)))
        hs = .001*np.sqrt(newer.nu*tau)
        ht = .00025*tau
        old_R, old_div = independent_fd(older, points, tau, hs, ht)
        new_R, div = independent_fd(newer, points, tau, hs, ht)
        row = dict(eta=eta, radius_fraction=fractions.tolist(),
                   baseline_norms=np.linalg.norm(old_R, axis=1).tolist(),
                   new_norms=np.linalg.norm(new_R, axis=1).tolist(),
                   baseline_angular=old_R[:, 1].tolist(),
                   new_angular=new_R[:, 1].tolist(),
                   baseline_divergence_max=float(np.max(np.abs(old_div))),
                   divergence_max=float(np.max(np.abs(div))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(k=k, tau=tau, attachment_radius=r1, rows=rows,
                  heat_age=newer.heat.heat_age,
                  scope='Point comparison only; pure heat exterior begins at old radial join radius. No whole-domain acceptance.',
                  accepted=False)
    tau_min = .5/64
    r_min = float(newer.attachment_radius(tau_min))
    p_edge = float(physical(np.array([[r_min, 0., 0.]]), tau_min,
                            c=newer.heat.amplitude, h=newer.heat.h,
                            nu=newer.nu)['pressure'][0])
    age = newer.heat.heat_age+.5-tau_min
    peak_point = np.array([[r_min, 0., np.sqrt(newer.nu*age)]])
    joined_u, joined_p = newer.fields(peak_point, tau_min)
    heat_u, heat_p = newer.heat.fields(peak_point, tau_min)
    velocity_gap = float(np.max(np.abs(joined_u-heat_u)))
    pressure_gap = float(np.max(np.abs(joined_p-heat_p)))
    if velocity_gap > 1e-12 or pressure_gap > 1e-12:
        raise AssertionError('The exterior edge does not equal the heat target')
    report['pure_exterior_spatial_time_bound'] = float(
        abs(p_edge)*newer.heat.heat_age/age*np.exp(-.5)
        /np.sqrt(newer.nu*age))
    report['bound_scope'] = ('Analytic spatial supremum for r>=attachment_radius(tau), '
                             'all z and tau in [1/128,1/2]; radial pressure '
                             'evaluated by quadrature, not interval-certified.')
    report['axial_energy_factor_at_tau_min'] = newer.heat.axial_energy_factor(tau_min)
    report['peak_edge_heat_velocity_gap'] = velocity_gap
    report['peak_edge_heat_pressure_gap'] = pressure_gap
    (ROOT/'heated_interior_join.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
