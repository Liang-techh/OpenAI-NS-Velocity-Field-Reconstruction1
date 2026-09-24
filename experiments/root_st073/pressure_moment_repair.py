"""Trial compact pressure reconstruction for the axial weighted momentum defect.

This enforces int r (u_z^2+p) dr = 0 at each (z,tau), but must still be tested
against all three pointwise momentum components. It does not alter velocity.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from compact_balance_audit import radial_integrals
from compact_moment_audit import weighted_moments
from compact_potential import CompactPotentialField
from joined_field import ROOT, coordinates, independent_fd


class PressureMomentRepair:
    def __init__(self, strength=1., base=None, cache=None):
        self.base = base or CompactPotentialField()
        self.nu = self.base.nu
        self.strength = strength
        self.cache = {} if cache is None else cache

    def support(self, tau):
        return self.base.support(tau)

    def flux(self, tau, z):
        key = (float(tau), float(z))
        if key not in self.cache:
            values = radial_integrals(self.base, tau, z)
            self.cache[key] = float(values[3] + values[4])
        return self.cache[key]

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        if self.strength == 0:
            return velocity, pressure
        pressure = pressure.copy()
        for i, (point, t) in enumerate(zip(pts, ts)):
            r = np.hypot(point[0], point[1])
            z = point[2]
            _, rsupp, _, zsupp = self.base.support(t)
            if abs(z) >= zsupp or r >= rsupp:
                continue
            q = float(coordinates(0., z / np.sqrt(self.nu), t, self.base.base.inner.h)['q'])
            ri = np.sqrt(self.nu) * np.sqrt(2 * q * 3 / 64)
            if r <= ri:
                continue
            y = (r - ri) / (rsupp - ri)
            # Beta(6,6)=1/2772; int r*y^5*(1-y)^5 dr=(rsupp^2-ri^2)/5544.
            bubble = 5544 * y**5 * (1 - y)**5 / (rsupp**2 - ri**2)
            pressure[i] -= self.strength * self.flux(t, z) * bubble
        return velocity, pressure


def run():
    base = CompactPotentialField()
    repaired = PressureMomentRepair(base=base)
    tau = .5 * 2**(-5.5)
    _, rsupp, _, zsupp = base.support(tau)
    nodes, weights = leggauss(5)
    radii = rsupp * (nodes + 1) / 2
    zs = zsupp * nodes
    points = np.array([[r, 0., z] for z in zs for r in radii])
    volume_weights = (2 * np.pi * np.tile(radii, 5)
                      * np.tile(weights * rsupp / 2, 5)
                      * np.repeat(weights * zsupp, 5))
    hs = .0005 * np.sqrt(base.nu * tau)
    ht = .0001 * tau
    residual0, div0 = independent_fd(base, points, tau, hs, ht)
    residual1, div1 = independent_fd(repaired, points, tau, hs, ht)
    direction = residual1 - residual0
    denominator = float(np.sum(volume_weights[:, None] * direction**2))
    least_squares_strength = (-float(np.sum(volume_weights[:, None] * residual0 * direction))
                              / denominator if denominator else 0.)
    # Pressure enters linearly, so these represent exact affine combinations
    # up to the same finite-difference and quadrature errors.
    rows = []
    for strength in (0., .25, .5, 1., least_squares_strength):
        res = residual0 + strength * direction
        norms = np.linalg.norm(res, axis=1)
        rows.append({
            'strength': strength,
            'sampled_max': float(np.max(norms)),
            'sampled_physical_volume_l2': float(np.sqrt(np.dot(volume_weights, norms**2))),
            'sampled_max_radial': float(np.max(np.abs(res[:, 0]))),
            'sampled_max_angular': float(np.max(np.abs(res[:, 1]))),
            'sampled_max_axial': float(np.max(np.abs(res[:, 2]))),
        })
    _, _, zflat, zsupp = base.support(tau)
    zcollar = (zflat + zsupp) / 2
    flux_values = radial_integrals(repaired, tau, zcollar)
    direct_collar = weighted_moments(repaired, tau, zcollar, 12)
    report = {
        'tau': tau, 'rows': rows,
        'collar_z': float(zcollar),
        'repaired_axial_flux_integral': float(flux_values[3] + flux_values[4]),
        'repaired_direct_axial_moment': direct_collar['axial_moment'],
        'repaired_direct_angular_moment': direct_collar['angular_moment'],
        'max_divergence_baseline': float(np.max(np.abs(div0))),
        'max_divergence_repaired': float(np.max(np.abs(div1))),
        'scope': 'Pressure-only radial beta bubble, normalized to cancel the physical axial weighted flux at strength one. Five-by-five cylindrical Gauss grid over compact support; this is a sampled candidate screen, not a full-domain gate.',
        'accepted': False,
    }
    out = ROOT / 'compact_potential' / 'pressure_repair.json'
    out.write_bytes((json.dumps(report, indent=2) + '\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
