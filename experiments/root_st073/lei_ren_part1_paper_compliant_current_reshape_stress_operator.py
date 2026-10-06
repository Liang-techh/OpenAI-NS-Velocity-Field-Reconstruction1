"""Actual variable long-reshape source in raw current-radius stress units."""
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import (
    SourceAST,copy_jet,IntervalTaylor,endpoints)
from lei_ren_part1_paper_compliant_long_reshape_profiles import relative_amplitude_jet
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import (
    exponential_derivatives,scaled_positive_source)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z


def raw_reshape_rows(c,packet,invP2):
    """Cancel raw primitive prefactors before original source exp bounds."""
    parent=packet['actual_inherited_axial5_packet']
    logu=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['log_Utheta_over_Pstar_axial5_coefficients']])
    shapes={key:IntervalTaylor(c,[c.mpf(endpoints(value)) for value in row])
        for key,row in parent['actual_normalized_moment_shape_axial5_coefficients'].items()}
    V=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['Uz_actual_axial5_coefficients']])
    log_y=[copy_jet(c,row) for row in packet['actual_log_Utheta_y_derivative_axial5']]
    bell=exponential_derivatives(log_y);bell2=exponential_derivatives([row*2 for row in log_y])
    amp=relative_amplitude_jet(logu);amp2=relative_amplitude_jet(logu,2)
    proofs=[]
    scaled=lambda power,row:scaled_positive_source(c,logu[0]*power,(amp if power==1 else amp2)*row,proofs)
    u=[scaled(1,row) for row in bell]
    u2=[scaled(2,row) for row in bell2]
    h=[shapes['theta']];k=[shapes['theta_z']];b=[shapes['swirl']]
    for j in range(1,5):
        h.append(bell[j-1]-h[j-1]*c.mpf('1.5'))
        k.append(bell[j-1]*V-k[j-1]*c.mpf('1.5'))
        b.append(bell2[j-1]-b[j-1])
    m=[shapes['mean']]+[(V-shapes['mean'])*((-1)**(j-1)) for j in range(1,5)]
    A=[shapes['axial']]+[(square(V)-shapes['axial'])*((-1)**(j-1)) for j in range(1,5)]
    histories=dict(m=m,h=[scaled(1,row) for row in h],k=[scaled(1,row) for row in k],
        e=[A[j]*invP2-scaled(2,b[j])/2 for j in range(5)],
        p=[scaled(2,shapes['pressure'])/2]+[u2[j-1]/2 for j in range(1,5)])
    p0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['pressure_axis_axial5_coefficients']])
    Q=[copy_jet(c,row) for row in packet['actual_Q_y_derivative_axial4']]
    return dict(histories=histories,absolute_pressure=[p0+histories['p'][0]]+histories['p'][1:],
        velocity=dict(theta=u,axial=[V]+[V*0]*4,radial=shifted_rows(Q,c.mpf('.5'),4)),
        original_variable_log_amplitude_ordinary_y_axial5=log_y,
        original_variable_exponential_Bell_rows=bell,
        original_variable_squared_exponential_Bell_rows=bell2,
        exact_actual_positive_source_log=logu[0],
        original_scaled_positive_source_cap_proofs=proofs,
        source_derivative_factors_combined_before_positive_caps=True,
        positive_caps_enclose_original_exponential_not_define_it=True)


def actual_reshape_normalization_theorem():
    """Replay original mixed rows and raw adapter for arbitrary source jets."""
    z=Z;delta=s.Symbol('delta',real=True);Ps,B=s.symbols('Pstar frozen_positive_swirl_unit',positive=True)
    c=SimpleNamespace(mpf=lambda value:(value[0] if isinstance(value,tuple) else s.Rational(str(value)) if isinstance(value,(str,int,float)) else value))
    J=lambda expr,order=5:FunctionJet.function(c,expr,order)
    ell=s.Function('actual_source_logu')(z);F=s.Function('actual_normalized_swirl_jet')(z)
    L=[J(s.Function('original_log_y'+str(j))(z)) for j in range(1,5)]
    V=J(s.Function('same_actual_V110')(z));p0=J(s.Function('same_actual_P0')(z))
    shapes={key:J(s.Function('actual_shape_'+key)(z)) for key in ('theta','theta_z','mean','axial','swirl','pressure')}
    def scaled(ctx,logbase,factor,proofs):
        power=s.cancel(logbase/ell)
        if power not in (1,2):raise ValueError('Only same original amplitude and square allowed')
        return factor*B**power
    env=dict(math=math,IntervalTaylor=FunctionJet,square=lambda value:value*value,
        derivative=lambda value:J(s.diff(value.expr,z),value.order-1),copy_jet=lambda ctx,value:value,
        endpoints=lambda value:(value,value),shifted_rows=shifted_rows,
        relative_amplitude_jet=lambda logjet,m=1:J(F**m),scaled_positive_source=scaled)
    asts=SourceAST()
    asts.replay('reference_restore_mixed_C4','binomial_rate',env)
    asts.replay('long_reshape_mixed_C4','exponential_derivatives',env)
    original=asts.replay('long_reshape_mixed_C4','reshape_mixed',env)
    adapter=asts.replay('current_reshape_stress_operator','raw_reshape_rows',env)
    packet=original(c,z,delta,J(ell),L,V,shapes,p0,1/Ps**2)
    packet.update(actual_inherited_axial5_packet=dict(
        log_Utheta_over_Pstar_axial5_coefficients=J(ell).coefficients,
        actual_normalized_moment_shape_axial5_coefficients={key:row.coefficients for key,row in shapes.items()},
        Uz_actual_axial5_coefficients=V.coefficients,pressure_axis_axial5_coefficients=p0.coefficients))
    result=adapter(c,packet,1/Ps**2);checks={};raw=result['histories'];velocity=result['velocity'];P=result['absolute_pressure']
    mapping=(('Utheta_over_current_Utheta',velocity['theta'],B),('Uz',velocity['axial'],1),
        ('Ur_over_current_sqrt_R_over_2',velocity['radial'],1),('P_over_Pstar2',P,1))
    for label,rows,factor in mapping:
        for j in range(5):
            for n in range(5-j):
                left=packet['physical_velocity_pressure_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)]*factor
                right=rows[j][n]*math.factorial(n)
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual reshape velocity/P0 source units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    mapping=(('Mz_over_current_R','m',1,1),('Mtheta_over_current_sqrt2_R_1p5_Utheta','h',s.Rational(3,2),B),
        ('Mtheta_z_over_current_sqrt2_R_1p5_Utheta','k',s.Rational(3,2),B),('Mztheta_over_current_R_Pstar2','e',1,1),('Mp_over_Pstar2','p',0,1))
    for label,key,rate,factor in mapping:
        rows=shifted_rows(raw[key],rate,4)
        for j in range(5):
            for n in range(5-j):
                left=packet['physical_five_primitive_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)]*factor
                right=rows[j][n]*math.factorial(n)
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual reshape five-moment source units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    u,v=velocity['theta'][0].expr,velocity['axial'][0].expr
    rhs=dict(m=v-raw['m'][0].expr,h=u-s.Rational(3,2)*raw['h'][0].expr,
        k=u*v-s.Rational(3,2)*raw['k'][0].expr,e=v*v/Ps**2-u*u/2-raw['e'][0].expr,p=u*u/2)
    for key,value in rhs.items():
        if s.cancel(raw[key][1].expr-value)!=0:raise ArithmeticError('Current reshape raw five ODE differs: '+key)
        checks['raw_five_history_ODE_'+key]=True
    if s.cancel(P[1].expr-u*u/2)!=0:raise ArithmeticError('Actual reshape absolute pressure ODE differs')
    checks['absolute_pressure_y_ODE']=True
    asts.method('long_reshape_profiles','relative_amplitude_jet')
    asts.method('long_reshape_mixed_C4','scaled_positive_source')
    asts.expression('long_reshape_mixed_C4','evaluate','log_y',
        wanted="[B*(-cutoff[k]*math.factorial(k)/T**k)+(c.mpf('.1') if k==1 else 0) for k in range(1,5)]")
    return dict(identities=checks,original_reshape_mixed_and_actual_raw_unit_adapter_AST_replayed=True,
        arbitrary_variable_ordinary_log_amplitude_y_axial5_rows_retained=True,
        arbitrary_positive_frozen_basepoint_unit_and_normalized_axial_source_jet_schema=True,
        original_Bell_and_inverse_radial_prefactor_binomial_identities_verified=True,
        full_current_V110_energy_baseline_and_analytic_P0_retained=True,
        original_source_derivative_factors_combined_before_positive_caps=True,
        constant_rate_used_only_at_flat_endpoints_not_interior=True,
        positive_caps_are_source_enclosures_not_point_selection=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes=asts.hashes,passed=True)
