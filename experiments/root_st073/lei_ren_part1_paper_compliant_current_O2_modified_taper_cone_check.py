"""Focused O2 source/pressure seam and open taper cone checks."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O2_modified_taper_cone as source
from lei_ren_part1_paper_compliant_current_O3_modified_transition_cone_check import (
    independent_signed_weighted_fixtures,independent_source_scale_identities)
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


def independent_negative_buffer_fixtures():
    """Complete signed original program with nonzero functional memory."""
    c=MPIntervalContext();c.dps=90;count=0
    def overlap(actual,expected,name):
        lo,hi=source.endpoints(actual);u,v=source.endpoints(expected)
        if hi<u or v<lo:raise ArithmeticError('Original signed O2 source fixture: '+name)
    with mp.workdps(120):
        for dv in ('.001','.05'):
            delta_scalar=c.mpf(dv)
            for tv in ('-2','-.731','0'):
                t=c.mpf(tv);Ua=c.mpf('.4');Ma=c.mpf('.17');da=c.mpf('.03')
                us=Ua*c.exp(-t/2);ms=Ma*c.exp(-t);ds=da*c.exp(-t)
                az=c.mpf('-.21');aq=c.mpf('-.37')-t/2;H0=c.mpf('.87');hf=1+(H0-1)*c.exp(t)
                for zv in ('-1','-.373','0','.51','1'):
                    constant=lambda v:IntervalTaylor(c,[v]+[c.mpf(0)]*5)
                    delta,U,M,D,AZ,AQ,Hf=map(constant,(delta_scalar,us,ms,ds,az,aq,hf))
                    z=IntervalTaylor(c,[c.mpf(zv),c.mpf(1)]+[c.mpf(0)]*4);zero=z*0
                    C=(1+z*z).reciprocal();L=1-delta*z*z
                    A=(1-delta/2)*C+(1-delta)*z*z*C*C
                    B=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
                    memory=(1+z*z)*c.mpf('.017')
                    u=[U*C*(-c.mpf('.5'))**j for j in range(5)]
                    hist=dict(m=[M*z*(-1)**j for j in range(5)],
                        h=[U*C*(-D*(-c.mpf('1.5'))**j+(-c.mpf('.5'))**j) for j in range(5)],
                        k=[U*(M-4*D)*z*C*(-c.mpf('1.5'))**j for j in range(5)],
                        e=[U**2*(-1)**j*(z*z*AZ+C*C*(AQ+c.mpf(j)/2)) for j in range(5)])
                    pressure=[-U**2*C*C*Hf/2+memory]+[-U**2*C*C*(-1)**j/2 for j in range(1,5)]
                    raw=raw_pre_stress_rows(c,delta,z,u,[zero]*5,hist,pressure)
                    theta=sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U
                    expected_theta=(A-C)/L+(-A/L-4*B)*D+(C+B)*M
                    overlap(theta[0],expected_theta[0],'theta')
                    overlap((raw['theta']['variable_radial_shear']['shape'][0]/U)[0],(-2*C)[0],'radial shear')
                    CE=z*z*AZ+C*C*AQ;Pbase=-C*C*Hf/2
                    expected_axial=(2*delta*z*CE-(1-z*z)*axial_derivative(CE))/L
                    expected_axial+=(2*(1+delta)*z*Pbase-(1-z*z)*axial_derivative(Pbase))/L
                    expected_axial+=(2*(1+delta)*z*memory-(1-z*z)*axial_derivative(memory))/(L*U**2)
                    actual=(raw['axial']['retained_full_energy']['shape'][0]+raw['axial']['actual_absolute_pressure']['shape'][0])/U**2
                    overlap(actual[0],expected_axial[0],'full energy/absolute pressure')
                    for name in ('local_axial_transport','nonlinear_meridional_transport','retained_linear_axial_moment','axial_radial_shear'):
                        for row in raw['axial'][name]['shape']:
                            lo,hi=source.endpoints(row[0])
                            floor=c.mpf(2)**(-c.prec+20)*(1+ms)
                            if not lo<=0<=hi or hi-lo>source.endpoints(floor)[1]:raise ArithmeticError('Source baseline axial zero fixture lost')
                    count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed current O2 source: '+name)
    field=source.CurrentModifiedO2TaperCone(require_checked=False);c=field.ctx
    for key,value in (('exact_original_O2_baseline_and_weighted_source_theorem',field.theorem),
        ('original_whole_closed_O2_taper_baseline_bounds',field.baseline),
        ('whole_modified_O2_taper_cone',field.proof),('complete_leading_stress_error_log_ledger',field.error_rows)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('O2 source replay changed: '+key)
    for name,value in {**field.baseline['positive_margins'],**field.proof['positive_margins']}.items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('O2 source inequality unresolved: '+name)
    for key,minimum in (('weighted_H_over_theta_lower','1.9'),('original_signed_D_over_theta_lower','.9'),
        ('original_signed_Q_over_theta_squared_lower','1.9'),('reduced_weighted_quadratic_over_theta_squared_lower','3.9')):
        if source.endpoints(field.proof[key])[0]<=source.endpoints(c.mpf(minimum))[0]:raise ArithmeticError('O2 signed cone reserve: '+key)
    specs=(('closed_direction_domain',(-2,0),22),('exact_degenerate_left_edge','-2',22),
        ('very_small_left_taper','-1.999',22),('whole_positive_taper',('-1.99','0'),22),
        ('old_seam','0',22),('current_frequency',('-1.999','0'),10**12))
    for name,t,N in specs:
        q=field.query(t,N)
        if source.parameters.encoded(q)!=data['examples'][name]:raise ValueError('O2 current query changed: '+name)
        strict=source.endpoints(c.mpf(t))[0]>-2
        if q['complete_modified_signed_cone_certified_for_entire_query_box']!=strict:
            raise ArithmeticError('O2 degenerate edge or open domain relabeled')
        if strict and source.endpoints(q['actual_vs_minus2_source_lower'])[0]<=0:
            raise ArithmeticError('Positive cutoff primitive floor lost')
    q,caps=field.errors.inputs('O2','-2',22);zeros=0
    for label,parts in field.errors.model['stress'].items():
        for name,p in parts.items():
            cap=source.o3.errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),field.errors.delta,0,p['mode'][0])[0,0]
            if source.endpoints(cap)!=(0,0):raise ArithmeticError('Flat source edge has nonzero stress modification: '+name)
            zeros+=1
    native_count=0
    for N in (37,10**12):
        for label,parts in field.leading_error_rows(N).items():
            for name,row in parts.items():
                if source.endpoints(row['native_coefficient_absolute_cap'])[1]>source.endpoints(field.error_rows[label][name]['native_coefficient_absolute_cap'])[1]:
                    raise ArithmeticError('O2 N increase enlarged native0 error')
                native_count+=1
    invalid=((-3,22),(1,22),(('-2.001','0'),22),(0,21),(0,True),(0,1.5),(mp.inf,22))
    for t,N in invalid:
        try:field.query(t,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid O2 source query admitted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Open O2 taper promoted closed/global gate')
    fixtures,admitted,rejected=independent_signed_weighted_fixtures();scales=independent_source_scale_identities(field)
    original=independent_negative_buffer_fixtures()
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_original_O2_source_pressure_seam_and_weighted_identities=len(field.theorem['identities']),
        strict_O2_baseline_and_transfer_inequalities=len(field.baseline['positive_margins'])+len(field.proof['positive_margins']),
        independent_original_negative_buffer_signed_fixtures=original,independent_original_signed_cone_fixtures=fixtures,
        signed_admitted_fixtures=admitted,signed_rejected_fixtures=rejected,
        independent_current_Pstar_radius_Ad_weighted_scale_identities=scales,exact_flat_edge_stress_error_zeros=zeros,
        larger_N_native_error_comparisons=native_count,current_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),
        all_N_source_monotonicity_not_sampled_frequency_argument=True,
        left_flat_endpoint_kept_degenerate_and_not_strict=True,
        current_modified_O2_open_taper_signed_two_vector_cone_certified=True,
        **{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name),
            Path(independent_signed_weighted_fixtures.__code__.co_filename).name:source.sha(Path(independent_signed_weighted_fixtures.__code__.co_filename).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified O2 open taper cone PASS: all N>=22;',original,'original negative-buffer fixtures; flat edge degenerate',flush=True)
    return result


if __name__=='__main__':run()
