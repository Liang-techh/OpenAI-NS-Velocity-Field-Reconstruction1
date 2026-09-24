"""Physical weighted radial moments of the compact field's full NS residual.

These are diagnostics motivated by the weighted primitives in paper Section 8.
They are not the paper's normalized mean defects or a stress realization.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from compact_potential import CompactPotentialField
from joined_field import ROOT, independent_fd


def weighted_moments(field, tau, z, order):
    rflat, rsupp, _, _ = field.support(tau)
    # Split at moving inner/outer interfaces and the radial cutoff to avoid
    # hiding narrow interface and collar contributions in one quadrature cell.
    q = tau if z == 0 else None
    if q is None:
        from joined_field import coordinates
        source = field
        while not hasattr(source, 'inner'):
            source = source.base
        q = float(coordinates(0., z / np.sqrt(field.nu), tau, source.inner.h)['q'])
    ri = np.sqrt(field.nu) * np.sqrt(2 * q * 3 / 64)
    boundaries = [0., ri, 6 * ri, rflat, rsupp]
    boundaries = sorted(set(np.clip(boundaries, 0., rsupp)))
    nodes, weights = leggauss(order)
    parts = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        if hi <= lo:
            continue
        radius = (hi + lo) / 2 + (hi - lo) / 2 * nodes
        w = (hi - lo) / 2 * weights
        points = np.column_stack((radius, np.zeros(order), np.full(order, z)))
        residual, divergence = independent_fd(
            field, points, tau, .0005 * np.sqrt(field.nu * tau), .0001 * tau
        )
        parts.append({
            'interval': [float(lo), float(hi)],
            'angular_moment': float(np.dot(w * radius**2, residual[:, 1])),
            'axial_moment': float(np.dot(w * radius, residual[:, 2])),
            'max_momentum': float(np.max(np.linalg.norm(residual, axis=1))),
            'max_divergence': float(np.max(np.abs(divergence))),
        })
    return {
        'z': float(z), 'order_per_piece': order,
        'angular_moment': sum(p['angular_moment'] for p in parts),
        'axial_moment': sum(p['axial_moment'] for p in parts),
        'pieces': parts,
    }


def run():
    field = CompactPotentialField()
    tau = .5 * 2**(-5.5)
    _, _, zflat, zsupp = field.support(tau)
    rows = []
    for z in (0., (zflat + zsupp) / 2):
        for order in (6, 12):
            row = weighted_moments(field, tau, z, order)
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = {
        'tau': tau,
        'rows': rows,
        'scope': 'Physical r^2 R_theta and r R_z moments of full unforced Cartesian NS residual at fixed z; not the paper normalized mean defects. Six and twelve Gauss nodes per piece are a limited quadrature check, not a certificate.',
        'pde_validated': False,
    }
    out = ROOT / 'compact_potential' / 'moment_audit.json'
    out.write_bytes((json.dumps(report, indent=2) + '\n').encode())


if __name__ == '__main__':
    run()
