"""Current reference/restoration histories in raw current-radius paper units."""
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import (
    SourceAST,copy_jet,IntervalTaylor,endpoints)
from lei_ren_part1_paper_compliant_reference_restore_mixed_C4 import binomial_rate
from lei_ren_part1_paper_compliant_long_reshape_profiles import relative_amplitude_jet
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import positive_exp
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z


def raw_restore_rows(c,packet,invP2):
    """Only enclose actual exp(logu); never use a cap as the source function."""
    parent=packet['actual_inherited_axial5_packet'];z=IntervalTaylor.variable(c,c.mpf(packet['Z']),5)
    logu=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['log_Utheta_over_Pstar_axial5_coefficients']])
    amp=relative_amplitude_jet(logu);amp2=relative_amplitude_jet(logu,2)
    u0=positive_exp(c,logu[0]);u20=positive_exp(c,2*logu[0])
    u=amp*u0;u2=amp2*u20
    source={key:[copy_jet(c,row) for row in rows] for key,rows in packet['actual_centered_moment_y_derivative_axial5'].items()}
    zero=z*0
    m=[source['mean_error'][j]+(z*4 if j==0 else zero) for j in range(5)]
    H=[source['angular_error'][j]+(c.mpf('.625') if j==0 else 0) for j in range(5)]
    K=[z*H[j]*4+source['mixed_error'][j] for j in range(5)]
    A=[z*source['mean_error'][j]*8+source['axial_square'][j]+(square(z)*16 if j==0 else zero) for j in range(5)]
    b=[source['swirl_error'][j]+(c.mpf(5)/6 if j==0 else 0) for j in range(5)]
    p=[source['pressure_error'][j]+(5 if j==0 else 0) for j in range(5)]
    histories=dict(m=m,h=[u*binomial_rate(H,c.mpf('.1'),j) for j in range(5)],
        k=[u*binomial_rate(K,c.mpf('.1'),j) for j in range(5)],
        e=[A[j]*invP2-u2*binomial_rate(b,c.mpf('.2'),j)/2 for j in range(5)],
        p=[u2*binomial_rate(p,c.mpf('.2'),j)/2 for j in range(5)])
    E=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['actual_E_V110_minus_4Z_axial5_coefficients']])
    alpha=[c.mpf(endpoints(value)) for value in packet['actual_alpha_ordinary_y_derivatives']]
    V=[z*4+E*alpha[0]]+[E*value for value in alpha[1:]]
    p0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['pressure_axis_axial5_coefficients']])
    P=[p0+histories['p'][0]]+histories['p'][1:]
    Q=[copy_jet(c,row) for row in packet['actual_Q_y_derivative_axial4']]
    return dict(histories=histories,absolute_pressure=P,
        velocity=dict(theta=[u*c.mpf('.1')**j for j in range(5)],axial=V,radial=shifted_rows(Q,c.mpf('.5'),4)),
        actual_positive_swirl_base_enclosure=u0,actual_positive_squared_swirl_base_enclosure=u20,
        exact_actual_positive_source_log=logu[0],positive_caps_enclose_original_exponential_not_define_it=True)


def actual_restore_normalization_theorem():
    """Replay both row programs in frozen-basepoint positive-unit jet algebra."""
    z=Z;delta=s.Symbol('delta',real=True);Ps,B=s.symbols('Pstar frozen_positive_swirl_unit',positive=True)
    c=SimpleNamespace(mpf=lambda v:(v[0] if isinstance(v,tuple) else s.Rational(str(v)) if isinstance(v,(str,int,float)) else v))
    J=lambda expr,order=5:FunctionJet.function(c,expr,order)
    ell=s.Function('actual_source_logu')(z);F=s.Function('actual_normalized_swirl_jet')(z)
    E=J(s.Function('actual_current_E')(z));p0=J(s.Function('same_actual_P0')(z))
    centered={key:J(s.Function('actual_'+key)(z)) for key in
        ('mean_error','angular_error','mixed_error','axial_square','swirl_error','pressure_error')}
    alpha=list(s.symbols('original_alpha0:5',real=True))
    def positive(ctx,value):
        power=s.cancel(value/ell)
        if power not in (1,2):raise ValueError('Only the same original amplitude and its square may be mapped')
        return B**power
    env=dict(math=math,IntervalTaylor=FunctionJet,square=lambda value:value*value,
        derivative=lambda value:J(s.diff(value.expr,z),value.order-1),copy_jet=lambda ctx,value:value,
        endpoints=lambda value:(value,value),shifted_rows=shifted_rows,
        relative_amplitude_jet=lambda logjet,m=1:J(F**m),positive_exp=positive)
    asts=SourceAST();asts.replay('reference_restore_mixed_C4','binomial_rate',env)
    original=asts.replay('reference_restore_mixed_C4','source_mixed',env)
    adapter=asts.replay('current_restore_stress_operator','raw_restore_rows',env)
    packet=original(c,z,delta,E,alpha,centered,J(ell),p0,1/Ps**2)
    packet.update(Z=z,actual_inherited_axial5_packet=dict(log_Utheta_over_Pstar_axial5_coefficients=J(ell).coefficients,
        actual_E_V110_minus_4Z_axial5_coefficients=E.coefficients,pressure_axis_axial5_coefficients=p0.coefficients))
    result=adapter(c,packet,1/Ps**2);checks={};raw=result['histories'];velocity=result['velocity'];P=result['absolute_pressure']
    mapping=(('Utheta_over_current_Utheta',velocity['theta'],B),('Uz',velocity['axial'],1),
        ('Ur_over_current_sqrt_R_over_2',velocity['radial'],1),('P_over_Pstar2',P,1))
    for label,rows,factor in mapping:
        for j in range(5):
            for n in range(5-j):
                left=packet['physical_velocity_pressure_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)]*factor
                right=rows[j][n]*math.factorial(n)
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual restore velocity/P0 source units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    mapping=(('Mz_over_current_R','m',1,1),('Mtheta_over_current_sqrt2_R_1p5_Utheta','h',s.Rational(3,2),B),
        ('Mtheta_z_over_current_sqrt2_R_1p5_Utheta','k',s.Rational(3,2),B),('Mztheta_over_current_R_Pstar2','e',1,1),('Mp_over_Pstar2','p',0,1))
    for label,key,rate,factor in mapping:
        rows=shifted_rows(raw[key],rate,4)
        for j in range(5):
            for n in range(5-j):
                left=packet['physical_five_primitive_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)]*factor
                right=rows[j][n]*math.factorial(n)
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual restore five-moment source units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    u,v=velocity['theta'][0].expr,velocity['axial'][0].expr
    rhs=dict(m=v-raw['m'][0].expr,h=u-s.Rational(3,2)*raw['h'][0].expr,
        k=u*v-s.Rational(3,2)*raw['k'][0].expr,e=v*v/Ps**2-u*u/2-raw['e'][0].expr,p=u*u/2)
    for key,value in rhs.items():
        if s.cancel(raw[key][1].expr-value)!=0:raise ArithmeticError('Current restoration raw five ODE differs: '+key)
        checks['raw_five_history_ODE_'+key]=True
    if s.cancel(P[1].expr-u*u/2)!=0:raise ArithmeticError('Actual restore absolute pressure ODE differs')
    checks['absolute_pressure_y_ODE']=True
    asts.method('long_reshape_profiles','relative_amplitude_jet')
    asts.expression('reference_restore_mixed_C4','source_mixed','ratio0',wanted='positive_exp(c,2*logu[0])')
    return dict(identities=checks,original_source_mixed_and_actual_raw_unit_adapter_AST_replayed=True,
        arbitrary_positive_frozen_basepoint_unit_and_normalized_axial_source_jet_schema=True,
        amplitude_unit_constant_not_differentiated_as_Z_function=True,
        actual_original_relative_exponential_jet_and_positive_enclosure_programs_bound=True,
        actual_generic_physical_fixture_consumed_separately_not_relabelled_current_source=True,
        full_axial_baseline_and_energy_cross_terms_retained=True,
        actual_P0_added_once_before_physical_differentiation=True,
        positive_cap_not_selected_as_source_amplitude=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes=asts.hashes,passed=True)
