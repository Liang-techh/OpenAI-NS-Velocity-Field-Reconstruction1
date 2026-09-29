"""Small independent checks of imported Part I identities, not a PDE gate."""
import json
from pathlib import Path
import numpy as np
from lei_ren_part1 import (completed_stress_tensor, projected_divergence,
                          background_remainder, parameter_map,
                          linear_core_axis_slopes, sector_tau_bounds, SOURCE)


def run():
    points = np.array([[.7,.4,.3], [-.4,.8,-.2], [.3,-.5,.7]])
    def tensor(p):
        r, z = np.hypot(p[:,0], p[:,1]), p[:,2]
        return completed_stress_tensor(p, r*r*z, r*z*z, 2*r*z)
    # Manufactured stress jets; differentiate the Cartesian tensor independently.
    actual = np.zeros_like(points)
    step = 1e-5
    for axis in range(3):
        off = np.zeros(3); off[axis] = step
        derivative = (tensor(points-2*off)-8*tensor(points-off)
                      +8*tensor(points+off)-tensor(points+2*off))/(12*step)
        actual += derivative[:,:,axis]
    r, z = np.hypot(points[:,0],points[:,1]), points[:,2]
    cyl = projected_divergence(r, r*r*z, r*z*z, 2*r*z, z*z)
    expected = np.column_stack((-points[:,1]/r*cyl[:,1], points[:,0]/r*cyl[:,1], cyl[:,2]))
    error = float(np.max(np.abs(actual-expected)))
    assert error < 1e-8
    assert np.max(abs(tensor(points)-tensor(points).transpose(0,2,1))) < 1e-14
    residual = np.array([[3.,4.,5.]])
    assert background_remainder(residual,np.array([[0.,-4.,-5.]])).tolist() == [[3.,0.,0.]]
    slopes = linear_core_axis_slopes(np.array([0.]), .002, 1., 0., .01, 4., -1., 0.)
    assert np.allclose(slopes[0], (-3+.001)/4)
    assert parameter_map(.005)['in_part1_stated_delta_range'] is False
    assert parameter_map(.001)['in_part1_stated_delta_range'] is True
    try: sector_tau_bounds(1.,1.)
    except ValueError: pass
    else: raise AssertionError('Endpoint sector must be rejected')
    report = dict(status='completed', source=SOURCE, source_version='2609.35406v1',
        cartesian_tensor_divergence_max_error=error,
        radial_remainder_preserved=True, axis_slope_check=True,
        current_h_mapping=parameter_map(.005), interior_sector_check=True,
        pde_validated=False, scale_recursion_established=False,
        scope='Manufactured algebraic checks only; no existing candidate is certified as a Part I background.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__': run()
