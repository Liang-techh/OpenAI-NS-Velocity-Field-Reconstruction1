"""Axisymmetric conservative balance behind the compact closure's radial defects.

For a compact divergence-free field at fixed z:
 int r^2 R_theta dr = d_t int r^2 u_theta dr
                     + d_z int r^2 u_z u_theta dr
                     - nu d_zz int r^2 u_theta dr;
 int r R_z dr = d_t int r u_z dr
                + d_z int r (u_z^2+p) dr - nu d_zz int r u_z dr.
Here d_t=-d_tau. These identities retain the full nonlinear flux and pressure.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from compact_potential import CompactPotentialField
from joined_field import ROOT, coordinates


def radial_integrals(field, tau, z, order=12):
    rflat, rsupp, _, _ = field.support(tau)
    q = float(coordinates(0., z / np.sqrt(field.nu), tau, field.base.inner.h)['q'])
    ri = np.sqrt(field.nu) * np.sqrt(2 * q * 3 / 64)
    boundaries = sorted(set(np.clip([0., ri, 6 * ri, rflat, rsupp], 0., rsupp)))
    nodes, weights = leggauss(order)
    result = np.zeros(5)
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        radius = (hi + lo) / 2 + (hi - lo) / 2 * nodes
        w = (hi - lo) / 2 * weights
        pts = np.column_stack((radius, np.zeros(order), np.full(order, z)))
        velocity, pressure = field.fields(pts, tau)
        theta, axial = velocity[:, 1], velocity[:, 2]
        result += [
            np.dot(w * radius**2, theta),
            np.dot(w * radius**2, axial * theta),
            np.dot(w * radius, axial),
            np.dot(w * radius, axial**2),
            np.dot(w * radius, pressure),
        ]
    return result


def decompose(field, tau, z, order=12):
    hs = .0005 * np.sqrt(field.nu * tau)
    ht = .0001 * tau
    center = radial_integrals(field, tau, z, order)
    tp2 = radial_integrals(field, tau + 2 * ht, z, order)
    tp1 = radial_integrals(field, tau + ht, z, order)
    tm1 = radial_integrals(field, tau - ht, z, order)
    tm2 = radial_integrals(field, tau - 2 * ht, z, order)
    zp2 = radial_integrals(field, tau, z + 2 * hs, order)
    zp1 = radial_integrals(field, tau, z + hs, order)
    zm1 = radial_integrals(field, tau, z - hs, order)
    zm2 = radial_integrals(field, tau, z - 2 * hs, order)
    dt = -(tm2 - 8 * tm1 + 8 * tp1 - tp2) / (12 * ht)
    dz = (zm2 - 8 * zm1 + 8 * zp1 - zp2) / (12 * hs)
    dzz = (-zp2 + 16 * zp1 - 30 * center + 16 * zm1 - zm2) / (12 * hs**2)
    angular_terms = [dt[0], dz[1], -field.nu * dzz[0]]
    axial_terms = [dt[2], dz[3], dz[4], -field.nu * dzz[2]]
    return {
        'z': float(z), 'order_per_piece': order,
        'angular_time': float(angular_terms[0]),
        'angular_axial_transport': float(angular_terms[1]),
        'angular_axial_viscosity': float(angular_terms[2]),
        'angular_sum': float(sum(angular_terms)),
        'axial_time': float(axial_terms[0]),
        'axial_kinetic_flux': float(axial_terms[1]),
        'axial_pressure_flux': float(axial_terms[2]),
        'axial_axial_viscosity': float(axial_terms[3]),
        'axial_sum': float(sum(axial_terms)),
        'axial_velocity_integral': float(center[2]),
    }


def run():
    field = CompactPotentialField()
    tau = .5 * 2**(-5.5)
    _, _, zflat, zsupp = field.support(tau)
    direct_rows = json.loads((ROOT / 'compact_potential' / 'moment_audit.json').read_text())['rows']
    rows = []
    for z in (0., (zflat + zsupp) / 2):
        row = decompose(field, tau, z)
        direct = next(r for r in direct_rows if r['z'] == z and r['order_per_piece'] == 12)
        row['angular_direct'] = direct['angular_moment']
        row['axial_direct'] = direct['axial_moment']
        row['angular_balance_difference'] = row['angular_sum'] - row['angular_direct']
        row['axial_balance_difference'] = row['axial_sum'] - row['axial_direct']
        print(json.dumps(row), flush=True)
        rows.append(row)
    report = {
        'tau': tau, 'rows': rows,
        'scope': 'Full physical conservative identities with 12 Gauss nodes per radial piece and fourth-order finite differences. Compare with separate residual moment audit; numerical quadrature and derivative errors remain.',
        'pde_validated': False,
    }
    out = ROOT / 'compact_potential' / 'balance_audit.json'
    out.write_bytes((json.dumps(report, indent=2) + '\n').encode())


if __name__ == '__main__':
    run()
