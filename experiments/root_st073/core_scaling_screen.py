"""Measure core swirl and relative axial extent for the current candidate."""
import json

import numpy as np
from scipy.optimize import minimize_scalar

from dynamic_poloidal_multitime_screen import load_dynamic_candidate
from joined_field import ROOT, coordinates
from local_poloidal_basis_screen import compact_base


def run():
    field = load_dynamic_candidate()
    compact = compact_base(field)
    rows = []
    for tau in (.128, .064, .032, .016, .0084):
        q = float(coordinates(0., 0., tau, compact.base.inner.h)['q'])
        inner_radius = np.sqrt(compact.nu)*np.sqrt(2*q*3/64)
        rflat, _, zflat, zsupport = field.support(tau)
        radii = np.linspace(.01*inner_radius, rflat, 120)
        points = np.column_stack((radii, np.zeros(len(radii)),
                                  np.zeros(len(radii))))
        velocity, _ = field.fields(points, tau)
        swirl = velocity[:, 1]
        peak = int(np.argmax(np.abs(swirl)))
        if peak in (0, len(radii)-1):
            raise ValueError('Swirl peak lies on the sampled radial boundary')
        def negative_abs_swirl(radius):
            point = np.array([[radius, 0., 0.]])
            return -abs(field.fields(point, tau)[0][0, 1])
        refined = minimize_scalar(negative_abs_swirl,
                                  bounds=(radii[peak-1], radii[peak+1]),
                                  method='bounded',
                                  options={'xatol': 1e-12})
        peak_radius = float(refined.x)
        peak_swirl = float(-refined.fun)
        row = {'tau': tau, 'inner_matching_radius': float(inner_radius),
               'midplane_swirl_peak_radius_grid': float(radii[peak]),
               'midplane_swirl_peak_abs': float(abs(swirl[peak])),
               'midplane_swirl_peak_radius_refined': peak_radius,
               'midplane_swirl_peak_abs_refined': peak_swirl,
               'midplane_peak_angular_speed_abs': peak_swirl/peak_radius,
               'axial_flat_half_length': float(zflat),
               'axial_support_half_length': float(zsupport),
               'axial_support_to_swirl_peak_radius': float(zsupport/peak_radius),
               'radial_flat_support': float(rflat),
               'peak_on_radial_grid_boundary': bool(peak in (0, len(radii)-1))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'rows': rows, 'radial_grid_count': 120,
              'scope': 'Current unaccepted finite-slab candidate. Midplane radial grid measures peak absolute azimuthal velocity, not vorticity-defined vortex-core radius or Lagrangian winding. Axial extent is mathematical support, not a measured coherent-vortex length.',
              'accepted': False}
    (ROOT/'compact_potential'/'core_scaling_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
