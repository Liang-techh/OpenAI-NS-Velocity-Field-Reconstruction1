"""Axis-correlation and first physical recovery check for transferred rows."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_finite_core import CompliantFiniteCore
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    f=CompliantFiniteCore();c=f.ctx
    with mp.workdps(220):
        zeros=0
        for field in ('A','Uz','P'):
            for m in range(len(f.rows[field][0])):
                if field=='P' and m==0:
                    if endpoints(f.cap(field,0,m))[0]<=0:
                        raise AssertionError('Pressure axis value shift lost')
                    continue
                if endpoints(f.cap(field,0,m))!=(mp.mpf(0),mp.mpf(0)):
                    raise AssertionError('Exact same-axis derivative correlation lost')
                if endpoints(f.coefficient(field,0,m))!=endpoints(f.rows[field][0][m]):
                    raise AssertionError('Unchanged axis row inflated')
                zeros+=1
        # Independently recover the change in the first physical axial slope:
        # Delta Uz_R=-(1+delta) Z Delta P0/L. The stored n=1 row divides by Lambda.
        dp=read_interval(c,f.major['pressure_perturbation']['physical_pressure_difference_abs_upper'])
        z=c.mpf(['.49','.51'])
        from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
        delta=CompliantPressureDatum().parameters.delta
        first_abs=f.eps*(1+delta)*z*dp/(1-delta*z**2)
        cap=f.cap('Uz',1,0)
        if endpoints(first_abs)[1]>endpoints(cap)[1]:
            raise AssertionError('First recovered axial sensitivity escapes transferred coefficient cap')
        samples=0
        for field in f.rows:
            for n in (1,72,144):
                for m in (0,len(f.rows[field][n])-1):
                    if endpoints(f.cap(field,n,m))[0]<=0:
                        raise AssertionError('Nonzero coefficient sensitivity cap lost')
                    value=f.coefficient(field,n,m);old=f.rows[field][n][m]
                    if not endpoints(value)[0]<=endpoints(old)[0]<=endpoints(old)[1]<=endpoints(value)[1]:
                        raise AssertionError('Transferred enclosure excludes old row')
                    samples+=1
        report=dict(exact_unchanged_axis_coefficients_checked=zeros,
            shifted_axis_pressure_value_remains_positive_cap=True,
            independent_first_axial_recovery_sensitivity_contained=True,
            sampled_early_middle_last_positive_caps_checked=samples,
            local_center_domain=['.49','.51'],radial_degree=144,
            row_units_verified=f.state['identity']['row_units'],
            new_finite_core_enclosure_view_checked=True,
            recomputed_new_point_coefficients=False,
            downstream_matching_or_temporal_recursion_claimed=False,all_passed=True,
            input_hashes={**f.hashes,
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
    print('Finite transfer:',zeros,'exact unchanged axis coefficients and independent first axial recovery PASS',flush=True)
    return report


if __name__=='__main__':run()
