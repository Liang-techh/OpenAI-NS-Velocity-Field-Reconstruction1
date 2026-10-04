"""Focused original-flatten full-moment/source-join check and one bounded fixture."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_flatten_stress_C3 import (
    CompliantFlattenStressC3,flatten_remaining_kernels,flatten_shape,flatten_defect_rows)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_original_sigma_fixture():
    """Direct scalar quadrature/derivatives of original sigma and f powers.

    Independent binomial coefficients recover six axial integral coefficients.
    This checks the new variable-F kernels/shape/full zeroth moments. The
    unchanged general-K stress/pressure operators consume their prior oracle.
    """
    with mp.workdps(85):
        a,mu,KR,t0,z0=map(mp.mpf,('.01','.025','1.3','50','.3'))
        d=a-mu;delta=2*a;p=1+delta;y=t0-100
        def sigma(v):
            x=v/100
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        def f(z):return (1+z*z)/2
        def G(v,z):return mp.exp((sigma(v)-1)*mp.log(f(z)))
        def density_coefficient(v,n):
            exponent=2*(sigma(v)-1);base=f(z0);lin=z0/base;quad=mp.mpf('.5')/base
            return base**exponent*sum((mp.binomial(exponent,j)*mp.binomial(j,n-j)
                *lin**(2*j-n)*quad**(n-j) for j in range((n+1)//2,n+1)),mp.mpf(0))
        direct={}
        for label,rate,scale in (('energy',2*mu,mp.mpf(1)),('pressure',1+2*mu,mp.mpf('.5'))):
            direct[label]=[scale*mp.quad(lambda v:mp.exp(-rate*(v-t0))*density_coefficient(v,n),
                [t0,75,100]) for n in range(6)]
        c=MPIntervalContext();c.dps=110
        heat=SimpleNamespace(ctx=c,a=c.mpf(a),mu=c.mpf(mu),delta=c.mpf(delta),k=c.mpf(1-a),prate=c.mpf(p))
        z=IntervalTaylor.variable(c,c.mpf(z0),5);one=IntervalTaylor.constant(c,1,5)
        ER=one*3+z*c.mpf('.2')+z*z*c.mpf('.08')
        PR=one*2-z*c.mpf('.1')+z*z*c.mpf('.03')
        terminal=dict(Kright=c.mpf(KR),energy_defect_rows=[ER-one/heat.delta],
            pressure_defect_rows=[PR-one/(2*heat.prate)])
        sig=sigma_jets(c,c.mpf(t0)/100);sj=[sig[j]*math.factorial(j)/100**j for j in range(5)]
        packet=dict(Z=c.mpf(z0),F_Taylor=one,actual_original_flatten_right_endpoint=False,sigma_y_derivatives=sj)
        # Direct F axial coefficients, independent of production exp/log jets.
        packet['F_Taylor']=IntervalTaylor(c,[c.mpf(mp.diff(lambda zz:mp.exp(sigma(t0)*mp.log(f(zz))),z0,n)/math.factorial(n)) for n in range(6)])
        shape=flatten_shape(heat,terminal,c.mpf(y),packet)
        kernels=flatten_remaining_kernels(c,c.mpf(mu),c.mpf(z0),c.mpf(t0),128)
        tolerance=mp.mpf('1e-55');miss=mp.mpf(0);counts=dict(kernel_axial5=0,K_mixed4=0,full_E_P_axial5=0)
        def compare(value,box,label):
            nonlocal miss
            lo,hi=endpoints(box);gap=max(lo-value,value-hi,mp.mpf(0));miss=max(miss,gap)
            if gap>tolerance:raise ArithmeticError('Independent original flatten mismatch: '+label)
        for label in direct:
            for n in range(6):
                compare(direct[label][n],kernels[label][n],label+' axial'+str(n));counts['kernel_axial5']+=1
        for j in range(5):
            for n in range(5-j):
                value=mp.diff(lambda tt,zz:KR*mp.exp(d*(tt-100))*G(tt,zz),(t0,z0),(j,n))
                compare(value,shape['K_rows'][j][n]*math.factorial(n),'K'+str((j,n)));counts['K_mixed4']+=1
        defects=flatten_defect_rows(heat,terminal,c.mpf(y),shape,one*c.mpf('.8'),kernels)
        # Full moments derived from independent integral coefficients and
        # arbitrary axial right data, rather than a baseline reset.
        er=[mp.mpf(3)+mp.mpf('.2')*z0+mp.mpf('.08')*z0*z0,mp.mpf('.2')+mp.mpf('.16')*z0,mp.mpf('.08'),0,0,0]
        pr=[mp.mpf(2)-mp.mpf('.1')*z0+mp.mpf('.03')*z0*z0,-mp.mpf('.1')+mp.mpf('.06')*z0,mp.mpf('.03'),0,0,0]
        for label,right,rate,baseline in (('energy',er,delta,1/delta),('pressure',pr,p,1/(2*p))):
            full=defects[label+'_defect_rows'][0]+one*c.mpf(baseline)
            for n in range(6):
                value=mp.exp(rate*y)*right[n]+KR**2*mp.exp(2*d*y)*direct[label][n]
                compare(value,full[n],label+' full axial'+str(n));counts['full_E_P_axial5']+=1
        if mp.diff(lambda zz:KR*mp.exp(d*y)*G(t0,zz),z0)==0:raise ArithmeticError('Variable K_Z fixture lost')
        return dict(all_passed=True,fixture='Original100-unit sigma, direct scalar quadrature and binomial axial coefficients',
            comparison_counts=counts,maximum_positive_enclosure_miss=str(miss),comparison_tolerance=str(tolerance),
            actual_nonzero_K_Z_exercised=True,moderate_parameters_only=True,
            actual_source_history_admission_is_separate=True,fixture_is_global_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'flatten_stress_C3.json';record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten source changed: '+source)
        fresh=CompliantFlattenStressC3().report()
        if encode(pack(fresh))!=record:raise ValueError('Current original flatten full source/moments/join differ')
        c=MPIntervalContext();c.dps=300;signed=0;pressure=0;zeros=0
        for point in record['samples']+[record['whole_flatten']]:
            for grid in point['flatten_similarity_stress_mixed3_factored'].values():
                for value in grid.values():
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite flatten stress')
                    signed+=1
            for j,jet in enumerate(point['pressure_y_derivative_axial5_Taylor']):
                for value in jet['coefficients'][:5-j]:
                    if not all(mp.isfinite(v) for v in endpoints(read_interval(c,value))):raise ArithmeticError('Nonfinite flatten pressure')
                    pressure+=1
            for jet in point['flatten_meridional_moments_and_velocities'].values():
                for value in jet['coefficients']:
                    if endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Flatten meridional zero lost')
                    zeros+=1
            if not point['ordinary_q_derivatives_not_rescaled']:raise ValueError('Original flatten ordinary derivatives lost')
        for flag in ('actual_original_flatten_similarity_stress_recovered','actual_flatten_absolute_pressure_same_source_mixed4_available',
            'flatten_power_stress_mixed3_and_pressure_mixed4_join_verified','actual_variable_K_Z_and_full_moment_axial_dependence_retained',
            'original_selected_outer_angular_full_future_retained','ordinary_q_derivatives_not_rescaled'):
            if not record[flag]:raise ValueError('Flatten admission missing: '+flag)
        for flag in ('flatten_physical_decomposition_constructed','flatten_cone_certified','left_pulse_stress_companion_constructed',
            'global_admissible_stress_lift_constructed','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Flatten similarity scope overclaimed: '+flag)
        if not record['source_bridge']['native_energy_function_join_proved_not_interval_overlap']:raise ValueError('Actual energy functional bridge missing')
        prior=json.loads((HERE/(PREFIX+'angular_stress_C3_check.json')).read_bytes())
        if not prior['all_passed'] or not prior['independent_original_beta_fixture']['all_passed']:raise ValueError('Accepted unchanged general-K operators required')
        for source,digest in prior['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('General-K oracle source changed: '+source)
        print('Current original flatten source/full moments/right join admitted; bounded original-sigma fixture running',flush=True)
        fixture=independent_original_sigma_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=fresh['actual_five_defect_family_sha256'],
            implicit_source_sha256=fresh['implicit_source_sha256'],finite_signed_stress_mixed_rows=signed,
            finite_pressure_mixed_rows=pressure,exact_meridional_zero_coefficients=zeros,
            actual_original_flatten_similarity_stress_recovered=True,actual_flatten_absolute_pressure_same_source_mixed4_available=True,
            flatten_power_stress_mixed3_and_pressure_mixed4_join_verified=True,actual_variable_K_Z_and_full_moment_axial_dependence_retained=True,
            original_selected_outer_angular_full_future_retained=True,actual_full_native_history_and_source_AST_bridge_verified=True,
            actual_native_flatten_power_energy_Z5_function_join_verified=True,
            actual_flatten_power_formula_join_identities=len(fresh['source_bridge']['actual_flatten_power_formula_join']['identities']),
            accepted_unchanged_general_K_stress_pressure_fixture_consumed=True,
            independent_original_sigma_fixture=fixture,flatten_physical_decomposition_constructed=False,
            flatten_cone_certified=False,left_pulse_stress_companion_constructed=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original100-unit flatten full moments/stress/absolute pressure and right power join; physical/cone/left pulse pending',flush=True)
        return result


if __name__=='__main__':run()
