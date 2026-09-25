"""Graft an axially heated finite-energy swirl tail onto the ST073 field."""
import json

import numpy as np
from scipy.special import expit

from axial_compact_join import load_compact_candidate
from axially_heated_exterior import AxiallyHeatedExterior
from heat_exterior import physical
from joined_field import coordinates, independent_fd
from radial_continuation import ROOT


def smooth_ramp(s):
    s = np.asarray(s, float)
    value = np.zeros_like(s)
    value[s >= 1] = 1.
    middle = (s > 0) & (s < 1)
    ss = s[middle]
    value[middle] = expit(1/(1-ss)-1/ss)
    return value


class HeatedGlobalJoin:
    def __init__(self, compact=None, outer_ratio=2., heat_age=2048.):
        self.compact = compact if compact is not None else load_compact_candidate()
        self.nu = self.compact.nu
        self.outer_ratio = float(outer_ratio)
        if self.outer_ratio <= 1:
            raise ValueError('Outer heat ramp must have positive width')
        self.heat = AxiallyHeatedExterior(
            nu=self.nu, h=self.compact.joined.inner.h,
            amplitude=self.compact.joined.c, heat_age=heat_age)

    def ramp_radii(self, tau):
        qmax = np.asarray(tau)/(1-self.compact.eta_outer**2)
        r1 = self.compact.joined.outer_ratio*np.sqrt(
            2*self.nu*qmax*self.compact.joined.join_X)
        return r1, self.outer_ratio*r1

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.compact.fields(pts, ts)
        r = np.hypot(pts[:, 0], pts[:, 1])
        r1, r2 = self.ramp_radii(ts)
        ramp = smooth_ramp((r-r1)/(r2-r1))
        active = ramp > 0
        if not np.any(active):
            return velocity, pressure
        chosen = pts[active]
        radial = physical(chosen, ts[active], c=self.heat.amplitude,
                          h=self.heat.h, nu=self.nu)
        sn = np.sqrt(self.nu)
        co = coordinates(r[active]/sn, chosen[:, 2]/sn,
                         ts[active], self.heat.h)
        chi = self.compact.cutoff(co['eta'])[0]
        G, _ = self.heat.envelope(chosen[:, 2], ts[active])
        velocity[active] += (ramp[active]*(G-chi))[:, None]*radial['velocity']
        pressure[active] += ramp[active]*(G*G-chi)*radial['pressure']
        return velocity, pressure


def load_heated_global_candidate():
    return HeatedGlobalJoin()


def run():
    field = load_heated_global_candidate()
    rows = []
    for k in (3., 5.5):
        tau = .5*2**(-k)
        eta = np.array([.2725, .345, .4175])
        X, E = np.meshgrid([2., 8.], eta, indexing='ij')
        points = field.compact.joined.inner.from_similarity(
            X.ravel(), E.ravel(), tau)
        residual, div = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        r1, r2 = field.ramp_radii(tau)
        far_points = points[X.ravel() == 8.]
        far_u, far_p = field.fields(far_points, tau)
        heat_u, heat_p = field.heat.fields(far_points, tau)
        row = dict(k=k, tau=tau, r1=float(r1), r2=float(r2),
                   X=X.ravel().tolist(), eta=E.ravel().tolist(),
                   residual_norms=np.linalg.norm(residual, axis=1).tolist(),
                   angular_residuals=residual[:, 1].tolist(),
                   divergence_max=float(np.max(np.abs(div))),
                   far_heat_velocity_difference=float(np.max(
                       np.abs(far_u-heat_u))),
                   far_heat_pressure_difference=float(np.max(
                       np.abs(far_p-heat_p))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(rows=rows, heat_age=field.heat.heat_age,
                  scope='Full-space finite-energy solenoidal candidate. '
                        'The outer radial region equals an axial/radial heat '
                        'swirl with exact angular heat evolution; an annular '
                        'ramp preserves the compact inner candidate. '
                        'Radial ramp, axial collar and axial pressure remain '
                        'unmatched dynamically; no whole-domain acceptance.',
                  pde_validated=False, accepted=False)
    (ROOT/'heated_global_join.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
