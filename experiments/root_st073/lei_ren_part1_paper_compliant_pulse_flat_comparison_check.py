"""Independent anchored ODE differences and original-source flat-limit gates.

The reference keeps the exact, possibly nonzero endpoint histories. A
closed-form exponential forcing checks both orientations without using the
producer's recurrence to compute the target derivatives.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_flat_comparison import difference_transport
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import symmetric
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_pulse_flat_comparison.json'


def anchored_difference_fixture():
    with mp.workdps(90):
        c=MPIntervalContext(); c.dps=100
        mu,a,h=mp.mpf('.04'),mp.mpf('.17'),mp.mpf('.2')
        ap=lambda z:mp.sqrt(1+z*z)+mp.mpf('.02')*z
        checks=0
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            coefficients=mp.taylor(ap,Z,5)
            rows=[IntervalTaylor(c,[symmetric(c,c.mpf(abs(v)+mp.mpf('1e-70'))*
                c.exp(c.mpf(a*h))*c.mpf(a)**k) for v in coefficients]) for k in range(5)]
            actual=difference_transport(c,rows,c.mpf(mu),c.mpf(h))
            targets={
                'linear_m1':lambda y,z:ap(z)*(mp.exp(a*y)-mp.exp(-(.5-mu)*y))/(a+.5-mu),
                'linear_m2':lambda y,z:ap(z)*(mp.exp(a*y)-mp.exp(-(.5-2*mu)*y))/(a+.5-2*mu),
                'energy':lambda y,z:ap(z)**2*(mp.exp(2*a*y)-mp.exp(2*mu*y))/(2*a-2*mu)}
            for key,function in targets.items():
                for y in (-h,h):
                    for k in range(5):
                        for n in range(5-k):
                            value=mp.diff(function,(y,Z),(k,n))/math.factorial(n)
                            lo,hi=endpoints(actual[key][k][n])
                            if not lo<=value<=hi:raise ArithmeticError('Anchored difference failed: '+str((key,y,Z,k,n)))
                            checks+=1
            for flag in ('angular_difference_exact_zero','swirl_pressure_difference_exact_zero'):
                if not actual[flag]:raise ArithmeticError('Unchanged histories lost')
        return dict(independent_closed_form_derivative_checks=checks,
            both_endpoint_orientations_checked=True,
            constant_energy_source_cancels_only_in_difference=True,
            actual_Md40_source_admission=False,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Flat comparison source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    norm=lambda v:max(abs(t) for t in endpoints(read(v)))
    checks=zero_checks=monotonic_checks=0; groups={}
    for point in r['comparisons']:
        groups.setdefault(point['label'],[]).append(point)
        for key in ('actual_nonzero_histories_and_positive_terminal_energy_preserved',
                    'h_zero_comparison_is_a_difference_not_zero_actual_field',
                    'angular_swirl_pressure_differences_exact_zero'):
            if not point[key]:raise ArithmeticError('Same exact endpoint histories required')
        values=[]
        for jet in point['B_y_derivative_difference_bounds']:values+=jet['coefficients']
        for rows in point['normalized_primitive_difference_bounds'].values():
            for jet in rows:values+=jet['coefficients']
        for grid in point['profile_velocity_mixed_difference_bounds'].values():
            if len(grid)!=15:raise ValueError('Flat mixed-index coverage incomplete')
            values+=list(grid.values())
        for grid in point['physical_cylindrical_difference_bounds'].values():
            if len(grid)!=15:raise ValueError('Physical flat-index coverage incomplete')
            for bound in grid.values():
                values.append(bound['factored_difference_absolute_upper'])
                if bound['physical_difference_log_upper'] is not None:
                    if any(not mp.isfinite(v) for v in endpoints(read(bound['physical_difference_log_upper']))):
                        raise ArithmeticError('Physical flat log bound nonfinite')
        is_zero=endpoints(read(point['h']))[1]==0
        for value in values:
            ends=endpoints(read(value))
            if any(not mp.isfinite(v) for v in ends):raise ArithmeticError('Flat difference nonfinite')
            checks+=1
            if is_zero:
                if ends!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Endpoint difference not exactly zero')
                zero_checks+=1
        if not is_zero and norm(point['B_y_derivative_difference_bounds'][0]['coefficients'][0])<=0:
            raise ArithmeticError('Positive original source cap reset to zero')
    for label,points in groups.items():
        if len(points)!=5:raise ValueError('Endpoint distance coverage incomplete')
        rows=[]
        for p in points:
            row=[norm(j['coefficients'][0]) for j in p['B_y_derivative_difference_bounds']]
            row += [norm(j['coefficients'][0]) for js in p['normalized_primitive_difference_bounds'].values() for j in js]
            row += [norm(v) for grid in p['profile_velocity_mixed_difference_bounds'].values() for v in grid.values()]
            row += [norm(v['factored_difference_absolute_upper']) for grid in p['physical_cylindrical_difference_bounds'].values() for v in grid.values()]
            rows.append(row)
        for left,right in zip(rows,rows[1:]):
            for a,b in zip(left,right):
                if b>a:raise ArithmeticError('Flat upper bound does not decrease: '+label)
                monotonic_checks+=1
        if not any(b<a for a,b in zip(rows[0],rows[-2])):raise ArithmeticError('No strict pre-endpoint decay: '+label)
    flat=json.loads((HERE/'lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json').read_bytes())
    if not flat['all_passed']:raise ValueError('Original shape flat proof required')
    proof=r['underlying_flat_limit_proof']
    if not all(proof[key] for key in ('finite_source_scales_at_each_fixed_positive_tau',
        'each_required_difference_derivative_tends_to_zero','numerical_endpoint_caps_also_vanish')):
        raise ValueError('Exact flat-limit composition unavailable')
    for key in ('full_pulse_C4_installed','full_outer_C4_certified','two_sided_external_high_order_joins_certified',
        'physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
        if r[key]:raise ValueError('Flat comparison promoted beyond its scope: '+key)
    fixture=anchored_difference_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],all_passed=True,
        anchored_difference_fixture=fixture,finite_difference_bounds_checked=checks,
        exact_zero_difference_bounds_checked=zero_checks,monotone_flat_bounds_checked=monotonic_checks,
        quantitative_flat_velocity_moment_comparison_ledger_available=True,
        same_nonzero_interface_histories_preserved=True,
        flat_limit_scope='Fixed source family and fixed positive physical tau; not uniform tau->0 or a stress remainder',
        full_pulse_C4_installed=False,full_outer_C4_certified=False,
        two_sided_external_high_order_joins_certified=False,physical_energy_integral_certified=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Original pulse flat comparison: independent two-orientation ODE differences and retained-history endpoint bounds PASS',flush=True)
    return out


if __name__=='__main__':run()
