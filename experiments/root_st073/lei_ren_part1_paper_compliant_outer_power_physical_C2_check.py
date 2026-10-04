"""Original preceding-power physical source join; reuse admitted general oracle."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_outer_power_physical_C2 import (
    CompliantOuterPowerStressC3,outer_power_physical_identities,
    outer_power_adapter_binding,angular_divergence_grids)
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import outer_power_shape
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_constant_K_divergence_fixture():
    """Differentiate full pressure and physical-factor source expressions.

    This bounded specialization tests nonzero P_Z with zero K_Z. The already
    accepted Cartesian oracle validates the unchanged full physical operator.
    """
    with mp.workdps(85):
        a,mu,q0,z0,KR,S=map(mp.mpf,('.1','.3','-.8','.4','2','.004'))
        delta=2*a;p=1+delta;d=a-mu;bh=mp.mpf('.5')+a;po=1+2*mu
        def K(q):return KR*mp.exp(d*q)
        def H(z):return mp.mpf('.1')+mp.mpf('.03')*z+mp.mpf('.02')*z*z
        def pressure(q,z):return K(q)**2/(2*po)+mp.exp(p*q)*H(z)
        def theta(q,z):
            return -d*K(q)/(1-delta*z*z)+2*S*mp.exp(-q)*mu*(1+mu)*K(q)
        def axial(q,z):
            return ((1-z*z)*mp.diff(lambda zz:pressure(q,zz),z)-2*p*z*pressure(q,z)+z*K(q)**2)/(1-delta*z*z)
        c=MPIntervalContext();c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),mu=c.mpf(mu),S=c.mpf(S),delta=c.mpf(delta),prate=c.mpf(p))
        z=IntervalTaylor.variable(c,c.mpf(z0),5);one=IntervalTaylor.constant(c,1,5)
        hj=one*c.mpf('.1')+z*c.mpf('.03')+z*z*c.mpf('.02')
        shape=outer_power_shape(heat,dict(Kright=c.mpf(KR)),c.mpf(q0))
        pd=shape['K_rows'][0]*shape['K_rows'][0]/(2*(1+2*heat.mu))+hj*c.exp(heat.prate*c.mpf(q0))-one/(2*heat.prate)
        grids=angular_divergence_grids(heat,dict(pressure_defect_rows=pd),shape,c.mpf(q0),c.mpf(z0))
        count=0;maxmiss=mp.mpf(0);tol=mp.mpf('1e-55')
        for label,rate,fn in (('theta',bh,theta),('axial',p,axial)):
            for j in range(3):
                for n in range(3-j):
                    value=mp.diff(lambda qq,zz:mp.exp(-rate*qq)*fn(qq,zz),(q0,z0),(j,n))*mp.exp(rate*q0)
                    lo,hi=endpoints(grids[label]['y'+str(j)+'_Z'+str(n)])
                    miss=max(lo-value,value-hi,mp.mpf(0));maxmiss=max(maxmiss,miss);count+=1
                    if miss>tol:raise ArithmeticError('Independent constant-K divergence specialization failed: '+str((label,j,n,miss)))
        if mp.diff(lambda zz:pressure(q0,zz),z0)==0:raise ArithmeticError('Fixture pressure axial dependence missing')
        return dict(all_passed=True,mixed2_comparisons=count,maximum_positive_enclosure_miss=str(maxmiss),
            comparison_tolerance=str(tol),K_Z_exact_zero=True,P_Z_nonzero=True,fixture_is_global_NS_validation=False)


def accepted_general_Cartesian_oracle(hashes):
    name=PREFIX+'angular_physical_C2_check.json';record=json.loads((HERE/name).read_bytes())
    if not record['all_passed']:raise ValueError('Accepted general physical oracle missing')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('General physical oracle source changed: '+path)
    fixture=record['independent_native_Cartesian_fixture']
    if not fixture['passed'] or fixture['nonzero_axial_viscosity_remainders_exercised']!=2 or fixture['reduced_divergence_mixed2_checks']!=24 or fixture['nonzero_remainder_mixed2_checks']!=12:
        raise ValueError('Required general-K Cartesian coverage missing')
    hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    return dict(accepted_receipt=name,passed=True,unchanged_general_operator_oracle_consumed=True,
        fixture_not_rerun=True,actual_outer_power_source_admitted_by_new_adapter_and_stress_receipt=True,
        nonzero_axial_viscosity_cases=2,divergence_mixed2_comparisons=24,remainder_mixed2_comparisons=12,
        original_fixture_tolerance=fixture['numeric_tolerance'],viscosities=fixture['viscosities'])


def run():
    with mp.workdps(300):
        name=PREFIX+'outer_power_physical_C2.json';record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Outer power physical source changed: '+path)
        provider=CompliantOuterPowerStressC3()
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(provider.family,provider.source):
            raise ValueError('Outer power physical family changed')
        if record['outer_power_physical_identities']!=outer_power_physical_identities():raise ValueError('Outer power physical identities changed')
        if record['outer_power_adapter_binding']!=outer_power_adapter_binding(provider):raise ValueError('Outer power adapter changed')
        c=MPIntervalContext();c.dps=300;finite=zeros=0
        def rowcheck(row):
            nonlocal finite,zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Nonfinite outer power physical row')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid power physical zero row')
                zeros+=1
            else:
                if any(not mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):
                    raise ArithmeticError('Nonfinite outer power physical logarithmic bound')
                finite+=1
        for point in record['samples']+[record['whole_outer_power']]:
            for flag in ('actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained',
                'actual_regional_physical_outer_power_stress_remainder_identity_verified',
                'full_axial_moment_history_and_pressure_baselines_cancelled_before_enclosure',
                'outer_power_angular_physical_stress_mixed3_and_remainder_mixed2_join_verified','original_phase_to_ordinary_q_factor_retained'):
                if not point[flag]:raise ValueError('Outer power physical gate missing: '+flag)
            for group,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6)):
                for grid in point[group].values():
                    if len(grid)!=count:raise ValueError('Outer power physical mixed coverage incomplete')
                    for row in grid.values():rowcheck(row)
            for group in ('completed_theta_theta_stress_mixed2','physical_angular_axial_viscosity_remainder_mixed2'):
                if len(point[group])!=6:raise ValueError('Outer power physical operator mixed2 coverage incomplete')
                for row in point[group].values():rowcheck(row)
            for label in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Outer power structural zero lost')
                zeros+=1
            for flag in ('outer_power_cone_certified','outer_power_regional_remainder_exact_zero','whole_outer_cone_certified',
                'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
                'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion','left_flatten_physical_stress_join_verified'):
                if point[flag]:raise ValueError('Outer power physical scope overclaimed: '+flag)
        oracle=accepted_general_Cartesian_oracle(hashes);fixture=independent_constant_K_divergence_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
            implicit_source_sha256=record['implicit_source_sha256'],
            actual_regional_physical_outer_power_stress_remainder_identity_verified=True,
            outer_power_angular_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,original_phase_to_ordinary_q_factor_retained=True,
            actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
            accepted_general_Cartesian_oracle=oracle,independent_constant_K_divergence_fixture=fixture,
            outer_power_regional_remainder_exact_zero=False,outer_power_cone_certified=False,
            left_flatten_physical_stress_join_verified=False,global_admissible_stress_lift_constructed=False,
            independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original preceding-power completed physical stress, full axial history, nonzero remainder and angular right join; cone/flatten pending',flush=True)
        return result


if __name__=='__main__':run()
