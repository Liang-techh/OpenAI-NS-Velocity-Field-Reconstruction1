"""One synthetic analytic fixture for the new normalized pair quadrature."""
import importlib.util
import sys
from pathlib import Path
import mpmath as mp


def run():
    path = Path(__file__).resolve().parents[1] / 'src/openai_ns_reconstruction/pulse_scaled_pair_quadrature_v2.py'
    spec = importlib.util.spec_from_file_location('new_scaled_pair_fixture_module_v2', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    with mp.workdps(160):
        eps = mp.mpf('1e-100')
        z = mp.mpf('0.3')
        def provider(x):
            r = mp.exp(x)
            return module.ScaledPairSample((1+eps)*r, r, eps*r,
                1+z*x, x, eps*(2+z*x), eps*x,
                z+x, 1, eps*(z+2*x), eps, 1)
        opts = dict(log_radius_scale='1e25', difference_scale='1e-100', precision_dps=150,
            radius_and_endpoints_Z_independent=True, measure_Z_independent=True)
        forward = module.integrate_scaled_pair(provider, 0, 1, **opts)
        reverse = module.integrate_scaled_pair(provider, 1, 0, **opts)
        zero = module.integrate_scaled_pair(provider, 1, 1, **opts)
        i0 = (mp.exp(2)-1)/2
        i1 = (mp.exp(2)+1)/4
        c = (1+eps)**2*eps
        d = eps*(2+eps)
        exact0 = (c*(2*i0+z*i1)+d*(i0+z*i1))/2
        exact1 = (c+d)*i1/2
        def boundary(x):
            r=mp.exp(x)
            return (r*((1+eps)*eps*(z+2*x)+eps*(z+x)),
                    r*((1+eps)*eps+eps))
        lo, hi = tuple(x/eps for x in boundary(0)), tuple(x/eps for x in boundary(1))
        exact0, exact1 = exact0/eps, exact1/eps
        tol=mp.mpf('1e-140')
        def close(a,b): return abs(a-b) <= tol*abs(b) if b else a == 0
        checks = {
            'analytic_integral_Z0':close(forward.integral[0],exact0),
            'analytic_integral_Z1':close(forward.integral[1],exact1),
            'both_endpoints_Z0_Z1':all(close(a,b) for a,b in zip(forward.lower_boundary+forward.upper_boundary,lo+hi)),
            'reverse_integral_Z0_Z1':all(close(a,-b) for a,b in zip(reverse.integral,forward.integral)),
            'reverse_endpoints_swapped':all(close(a,b) for a,b in zip(reverse.lower_boundary+reverse.upper_boundary,hi+lo)),
            'zero_width_integral':zero.integral==(0,0),
            'scale_separation':abs(forward.integral_log_scale-(mp.mpf('2e25')+mp.log(eps)))<mp.mpf('1e-120') and abs(forward.boundary_log_scale-(mp.mpf('1e25')+mp.log(eps)))<mp.mpf('1e-120'),
            'uncertified':not forward.quadrature_error_certified and not forward.source_identity_proved,
        }
        try:
            module.integrate_scaled_pair(provider,0,1,log_radius_scale=0,difference_scale=1,precision_dps=50,
                radius_and_endpoints_Z_independent=False,measure_Z_independent=True)
            checks['reject_missing_Z_premise']=False
        except ValueError:
            checks['reject_missing_Z_premise']=True
        return {'checks': checks, 'all_passed': all(checks.values()),
            'precision_dps':150,'synthetic_difference':'1e-100',
            'relative_integral_discrepancy_Z0':mp.nstr(abs((forward.integral[0]-exact0)/exact0),12),
            'relative_integral_discrepancy_Z1':mp.nstr(abs((forward.integral[1]-exact1)/exact1),12),
            'quadrature_estimates_are_certified':False,
            'difference_scale':'1e-100','actual_source_values_computed':False,'source_identity_proved':False}
